"""Train the MLPhoton mass regressor (m/E) on the packed diphoton sample.

Loss
----
The target spans four orders of magnitude (m/E from ~1e-4 at mA=0.1 GeV up to
~0.9 at mA=10), so a plain MSE would be driven entirely by the high-mass end and
would leave the sub-GeV region -- the part this analysis needs -- untrained.

A RELATIVE loss, L = mean(((pred-target)/target)^2), looks like the fix and is
not: predicting zero costs exactly 1.0 per sample while overshooting costs far
more, so the optimiser parks the output at zero. Measured on this dataset, from
scratch it converges to bias = -0.001 with val loss 1.00 -- worse than the
untrained deployed model (0.63). Do not use it as the objective.

What is used instead: MSE on a standardised log target, with the network
emitting exp(raw*sigma + mu) so the prediction is positive by construction (see
MLPhotonRegressorLog in models.py). Every mass point then contributes on the
same scale and there is no zero-collapse basin.

Only the head changes; the convolutional trunk is untouched, and the exported
graph still returns a single m/E scalar, which is all MLPhotonProducer reads.
Pass --head deployed to evaluate the shipped model under its original head.

Input : <packed>/regressor_{train,val}.npz   (from build_dataset.py)
Output: <out>/regressor_best.pt , <out>/regressor_last.pt , <out>/history.csv

Usage:
    python train_regressor.py --packed <dir> --out <dir> --epochs 200
    # add --warm-start <regressor.onnx> to initialise from the current model
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import (MLPhotonRegressor, MLPhotonRegressorLinear,  # noqa: E402
                    MLPhotonRegressorLog, load_onnx_weights)


def load_split(packed, split, keep_zero=False, neg_band=0.0, seed=0):
    """Returns (img, eta, moe_true, label) tensors plus numpy (energy, m_gen).

    label is the TRAINING target: equal to moe_true except that, with
    neg_band > 0, a fixed random half of the rows with 0 < m/E < neg_band get
    -m/E. Every metric is computed against moe_true, never against label.
    keep_zero keeps the m/E = 0 single-photon anchors (linear head only; the
    log head cannot represent them).
    """
    path = os.path.join(packed, f"regressor_{split}.npz")
    with np.load(path) as z:
        img = z["image"].astype(np.float32)
        moe = z["moe_true"].astype(np.float32)
        eta = z["eta"].astype(np.float32)
        energy = z["energy"].astype(np.float32)
        mgen = (z["m_gen"].astype(np.float32) if "m_gen" in z.files
                else np.full(len(moe), np.nan, dtype=np.float32))
    # regressor input is the sum-normalised image (see cluster_py)
    img = img.reshape(-1, 1, 30, 30)
    sums = img.sum(axis=(1, 2, 3), keepdims=True)
    img = np.divide(img, sums, out=np.zeros_like(img), where=sums > 0)
    good = np.isfinite(moe) & ((moe >= 0) if keep_zero else (moe > 0))
    print(f"[data] {split}: {int(good.sum())} rows "
          f"(dropped {int((~good).sum())} with {'negative' if keep_zero else 'non-positive'} "
          f"target; {int(np.sum(moe[good] == 0))} m/E=0 anchors)")
    img, eta, moe, energy, mgen = img[good], eta[good], moe[good], energy[good], mgen[good]
    label = moe.copy()
    if neg_band > 0:
        rng = np.random.default_rng(seed)
        flip = (moe > 0) & (moe < neg_band) & (rng.random(len(moe)) < 0.5)
        label[flip] = -moe[flip]
        print(f"[data] {split}: {int(flip.sum())} rows in (0, {neg_band}) "
              f"carry a negative label")
    return ((torch.from_numpy(img), torch.from_numpy(eta).reshape(-1, 1),
             torch.from_numpy(moe), torch.from_numpy(label)),
            energy, mgen)


def log_mse(model, img, eta, target):
    """MSE on the standardised log target -- see MLPhotonRegressorLog."""
    z = (torch.log(target) - model.mu) / model.sigma
    return F.mse_loss(model.trunk(img, eta), z)


def linear_mse(model, img, eta, label):
    """MSE on the standardised LINEAR target -- see MLPhotonRegressorLinear."""
    z = (label - model.mu) / model.sigma
    return F.mse_loss(model.trunk(img, eta), z)


def relative_mse(pred, target):
    """Kept only for reporting: it is comparable across runs, but must NOT be
    used as the training objective (it has a zero-collapse minimum -- see
    MLPhotonRegressorLog)."""
    return torch.mean(((pred - target) / target) ** 2)


@torch.no_grad()
def evaluate(model, loader, device, head="log"):
    """Model selection uses the SAME objective the model is trained on.

    The first run selected on relative MSE while training on log MSE. Those
    disagree badly: relative MSE is dominated by the smallest targets, so it
    climbed from 19.5 to 37.5 while the actual mass scale (bias) stayed between
    1.0 and 1.3. Selecting on a metric the optimiser is not minimising picks
    essentially arbitrary epochs.
    """
    model.eval()
    preds, targs, labs = [], [], []
    for img, eta, moe, lab in loader:
        p = model(img.to(device), eta.to(device)).cpu()
        preds.append(p)
        targs.append(moe)
        labs.append(lab)
    pred_all = torch.cat(preds).numpy().astype(np.float64)
    targ_all = torch.cat(targs).numpy().astype(np.float64)
    lab_all = torch.cat(labs).numpy().astype(np.float64)
    # the training objective for the linear head, on the same signed labels
    lin_loss = float(np.mean(((pred_all - lab_all) / float(model.sigma)) ** 2)) \
        if head == "linear" else float("nan")
    anchor = targ_all == 0
    anchor_med = float(np.median(pred_all[anchor])) if anchor.any() else float("nan")
    # everything below describes the real diphoton rows only (m/E > 0)
    pred, targ = pred_all[~anchor], targ_all[~anchor]
    ratio = pred / targ
    finite = np.isfinite(ratio) & (targ > 0)
    ratio = ratio[finite]
    rel = float(np.mean(((pred[finite] - targ[finite]) / targ[finite]) ** 2))
    # log-space MSE: the training objective, and a sane scale-free error measure
    pos = finite & (pred > 0)
    log_loss = (float(np.mean((np.log(pred[pos]) - np.log(targ[pos])) ** 2))
                if pos.any() else float("inf"))
    # a prediction <= 0 cannot be scored in log space; penalise rather than drop
    frac_bad = float(np.mean(~pos[finite])) if finite.any() else 1.0
    log_loss += 10.0 * frac_bad
    loss = {"log": log_loss, "linear": lin_loss}.get(head, rel)
    q16, q50, q84 = np.quantile(ratio, [0.16, 0.5, 0.84])
    return {
        "loss": loss,
        "log_loss": log_loss,
        "rel_loss": rel,
        "bias": float(q50),                       # median pred/true; 1.0 is perfect
        "resolution": float((q84 - q16) / 2.0),   # half the central 68% spread
        "frac_neg": float(np.mean(pred[finite] <= 0)),
        "anchor_median": anchor_med,
        "pred_all": pred_all, "targ_all": targ_all,
        "pred": pred[finite], "targ": targ[finite],
    }


def report_by_moe(res, edges=(0, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05, 1.0)):
    """Bias/resolution per m/E band -- the old model's failure was band-specific,
    so a single average would hide whether that is fixed."""
    pred, targ = res["pred"], res["targ"]
    print(f"    {'m/E band':>16s} {'N':>8s} {'bias':>7s} {'res':>7s}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (targ >= lo) & (targ < hi)
        if m.sum() < 20:
            continue
        r = pred[m] / targ[m]
        q16, q50, q84 = np.quantile(r, [0.16, 0.5, 0.84])
        print(f"    [{lo:6.3f},{hi:6.3f}) {int(m.sum()):>8d} {q50:>7.3f} "
              f"{(q84-q16)/2:>7.3f}")


def report_by_mgen(res, energy, mgen):
    """Reconstructed MASS per generated mass point -- what the per-mA ROI sees.
    The val loader is not shuffled, so its order matches energy/m_gen."""
    m_rec = res["pred_all"] * energy
    print(f"    {'m_gen':>6s} {'N':>7s} {'median':>8s} {'q16':>8s} {'q84':>8s} "
          f"{'width':>7s} {'med/m':>6s} {'frac<0':>6s}")
    for mg in sorted(set(np.round(mgen[np.isfinite(mgen)], 3))):
        if mg > 1.0:
            continue
        sel = np.isclose(mgen, mg)
        if sel.sum() < 20:
            continue
        q16, q50, q84 = np.quantile(m_rec[sel], [0.16, 0.5, 0.84])
        name = "mono" if mg == 0 else f"{mg:.1f}"
        print(f"    {name:>6s} {int(sel.sum()):>7d} {q50:>8.4f} {q16:>8.4f} {q84:>8.4f} "
              f"{q84-q16:>7.4f} {(q50/mg if mg > 0 else float('nan')):>6.2f} "
              f"{np.mean(m_rec[sel] < 0):>6.3f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--packed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--lr", type=float, default=1e-3)     # AN used torch defaults
    ap.add_argument("--warm-start", default=None, help="regressor.onnx to start from")
    ap.add_argument("--head", choices=("log", "linear", "deployed"), default="log",
                    help="'log' = standardised-log head (positive only); "
                         "'linear' = standardised linear head that may go "
                         "negative and uses the m/E=0 anchors; 'deployed' = the "
                         "original LeakyReLU head, for evaluating the shipped "
                         "model only")
    ap.add_argument("--neg-band", type=float, default=0.0,
                    help="linear head only: a fixed random half of the rows with "
                         "0 < m/E < NEG_BAND are trained on -m/E, so the target "
                         "distribution is symmetric about zero below the "
                         "resolution instead of being cut off there. Keep it "
                         "well under M0p1's m/E (~0.0027), or M0p1 itself gets "
                         "pulled to zero.")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--patience", type=int, default=30,
                    help="stop if val loss has not improved for this many epochs")
    ap.add_argument("--eval-only", action="store_true",
                    help="evaluate the warm-start (or --resume) weights and exit; "
                         "this is how the pre-training baseline is measured")
    ap.add_argument("--resume", default=None, help="a .pt to evaluate or continue from")
    ap.add_argument("--max-train", type=int, default=0,
                    help="cap training rows (0 = all); for quick timing tests")
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    # Line-buffer stdout: these runs are long and are watched through a log
    # file, where the default block buffering shows nothing for hours.
    sys.stdout.reconfigure(line_buffering=True)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    os.makedirs(args.out, exist_ok=True)
    device = torch.device(args.device)
    print(f"[train] device={device}")

    if args.neg_band > 0 and args.head != "linear":
        raise SystemExit("--neg-band needs --head linear (the log head cannot "
                         "represent a negative target)")
    linear = args.head == "linear"
    va, va_energy, va_mgen = load_split(args.packed, "val", keep_zero=linear,
                                        neg_band=args.neg_band, seed=args.seed + 1)
    val_loader = DataLoader(TensorDataset(*va), batch_size=1024)

    if linear:
        with np.load(os.path.join(args.packed, "regressor_train.npz")) as _z:
            _t = _z["moe_true"].astype(np.float64)
        _t = _t[np.isfinite(_t) & (_t >= 0)]
        model = MLPhotonRegressorLinear()
        model.set_target_scale(_t)
        print(f"[train] linear-target scale: mu={float(model.mu):.5f} "
              f"sigma={float(model.sigma):.5f}")
        model = model.to(device)
    elif args.head == "log":
        # scale set from the TRAINING targets; needed before any evaluation
        with np.load(os.path.join(args.packed, "regressor_train.npz")) as _z:
            _t = _z["moe_true"].astype(np.float64)
        _t = _t[np.isfinite(_t) & (_t > 0)]
        model = MLPhotonRegressorLog()
        model.set_target_scale(np.log(_t))
        print(f"[train] log-target scale: mu={float(model.mu):.4f} "
              f"sigma={float(model.sigma):.4f}")
        model = model.to(device)
    else:
        model = MLPhotonRegressor().to(device)
    if args.warm_start:
        load_onnx_weights(model, args.warm_start)
        print(f"[train] warm-started from {args.warm_start}")
    if args.resume:
        model.load_state_dict(torch.load(args.resume, map_location=device))
        print(f"[train] resumed from {args.resume}")
    if args.warm_start or args.resume:
        base = evaluate(model, val_loader, device, head=args.head)
        print(f"[train] starting point: loss {base['loss']:.4f} "
              f"bias {base['bias']:.3f} res {base['resolution']:.3f}")
        report_by_moe(base)
        print("[train] reconstructed mass per generated mass point (val):")
        report_by_mgen(base, va_energy, va_mgen)

    if args.eval_only:
        print("[train] --eval-only: done")
        return 0

    tr, _, _ = load_split(args.packed, "train", keep_zero=linear,
                          neg_band=args.neg_band, seed=args.seed)
    if args.max_train and args.max_train < len(tr[0]):
        tr = tuple(t[:args.max_train] for t in tr)
        print(f"[train] capped training rows to {args.max_train}")
    train_loader = DataLoader(TensorDataset(*tr), batch_size=args.batch_size,
                              shuffle=True, drop_last=True)

    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=10)

    hist_path = os.path.join(args.out, "history.csv")
    hist = open(hist_path, "w", newline="")
    writer = csv.writer(hist)
    writer.writerow(["epoch", "train_loss", "val_loss", "log_loss", "rel_loss",
                     "bias", "resolution", "lr", "sec"])

    best = float("inf")
    since_best = 0
    for ep in range(1, args.epochs + 1):
        t0 = time.time()
        model.train()
        tot, nb = 0.0, 0
        for img, eta, moe, lab in train_loader:
            img, eta, moe, lab = img.to(device), eta.to(device), moe.to(device), lab.to(device)
            opt.zero_grad()
            if args.head == "log":
                loss = log_mse(model, img, eta, moe)
            elif linear:
                loss = linear_mse(model, img, eta, lab)
            else:
                loss = relative_mse(model(img, eta), moe)
            loss.backward()
            opt.step()
            tot += float(loss.item())
            nb += 1
        train_loss = tot / max(nb, 1)

        res = evaluate(model, val_loader, device, head=args.head)
        sched.step(res["loss"])
        lr_now = opt.param_groups[0]["lr"]
        dt = time.time() - t0
        writer.writerow([ep, f"{train_loss:.6f}", f"{res['loss']:.6f}",
                         f"{res['log_loss']:.6f}", f"{res['rel_loss']:.4f}",
                         f"{res['bias']:.4f}", f"{res['resolution']:.4f}",
                         f"{lr_now:.2e}", f"{dt:.1f}"])
        hist.flush()

        flag = ""
        if res["loss"] < best:
            best = res["loss"]
            since_best = 0
            torch.save(model.state_dict(), os.path.join(args.out, "regressor_best.pt"))
            flag = "  *best"
        else:
            since_best += 1

        if ep % 5 == 0 or flag or ep == 1:
            print(f"[ep {ep:4d}] train {train_loss:.4f}  val {res['loss']:.4f}  "
                  f"bias {res['bias']:.3f}  res {res['resolution']:.3f}  "
                  f"lr {lr_now:.1e}  {dt:.0f}s{flag}")

        if since_best >= args.patience:
            print(f"[train] no improvement for {args.patience} epochs; stopping at {ep}")
            break

    torch.save(model.state_dict(), os.path.join(args.out, "regressor_last.pt"))
    hist.close()

    model.load_state_dict(torch.load(os.path.join(args.out, "regressor_best.pt")))
    final = evaluate(model, val_loader, device, head=args.head)
    print(f"\n[train] best model: val loss {final['loss']:.4f}  "
          f"bias {final['bias']:.3f}  resolution {final['resolution']:.3f}  "
          f"frac_neg {final['frac_neg']:.3f}  m/E=0 anchors -> median {final['anchor_median']:.5f}")
    report_by_moe(final)
    print("\n[train] reconstructed mass per generated mass point (val):")
    report_by_mgen(final, va_energy, va_mgen)
    print(f"[train] history -> {hist_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
