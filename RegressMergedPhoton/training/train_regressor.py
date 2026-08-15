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
from models import (MLPhotonRegressor, MLPhotonRegressorLog,  # noqa: E402
                    load_onnx_weights)


def load_split(packed, split):
    path = os.path.join(packed, f"regressor_{split}.npz")
    with np.load(path) as z:
        img = z["image"].astype(np.float32)
        moe = z["moe_true"].astype(np.float32)
        eta = z["eta"].astype(np.float32)
        energy = z["energy"].astype(np.float32)
    # regressor input is the sum-normalised image (see cluster_py)
    img = img.reshape(-1, 1, 30, 30)
    sums = img.sum(axis=(1, 2, 3), keepdims=True)
    img = np.divide(img, sums, out=np.zeros_like(img), where=sums > 0)
    good = np.isfinite(moe) & (moe > 0)
    print(f"[data] {split}: {int(good.sum())} rows "
          f"(dropped {int((~good).sum())} with non-positive target)")
    return (torch.from_numpy(img[good]), torch.from_numpy(eta[good]).reshape(-1, 1),
            torch.from_numpy(moe[good]), torch.from_numpy(energy[good]))


def log_mse(model, img, eta, target):
    """MSE on the standardised log target -- see MLPhotonRegressorLog."""
    z = (torch.log(target) - model.mu) / model.sigma
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
    preds, targs = [], []
    for img, eta, moe, _ in loader:
        p = model(img.to(device), eta.to(device)).cpu()
        preds.append(p)
        targs.append(moe)
    pred = torch.cat(preds).numpy().astype(np.float64)
    targ = torch.cat(targs).numpy().astype(np.float64)
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
    loss = log_loss if head == "log" else rel
    q16, q50, q84 = np.quantile(ratio, [0.16, 0.5, 0.84])
    return {
        "loss": loss,
        "log_loss": log_loss,
        "rel_loss": rel,
        "bias": float(q50),                       # median pred/true; 1.0 is perfect
        "resolution": float((q84 - q16) / 2.0),   # half the central 68% spread
        "frac_neg": float(np.mean(pred[finite] <= 0)),
        "pred": pred[finite], "targ": targ[finite],
    }


def report_by_moe(res, edges=(0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 1.0)):
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--packed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--lr", type=float, default=1e-3)     # AN used torch defaults
    ap.add_argument("--warm-start", default=None, help="regressor.onnx to start from")
    ap.add_argument("--head", choices=("log", "deployed"), default="log",
                    help="'log' = standardised-log head (trainable); 'deployed' "
                         "= the original LeakyReLU head, for evaluating the "
                         "shipped model only")
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

    va = load_split(args.packed, "val")
    val_loader = DataLoader(TensorDataset(*va), batch_size=1024)

    if args.head == "log":
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

    if args.eval_only:
        print("[train] --eval-only: done")
        return 0

    tr = load_split(args.packed, "train")
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
        for img, eta, moe, _ in train_loader:
            img, eta, moe = img.to(device), eta.to(device), moe.to(device)
            opt.zero_grad()
            loss = (log_mse(model, img, eta, moe) if args.head == "log"
                    else relative_mse(model(img, eta), moe))
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
          f"bias {final['bias']:.3f}  resolution {final['resolution']:.3f}")
    report_by_moe(final)
    print(f"[train] history -> {hist_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
