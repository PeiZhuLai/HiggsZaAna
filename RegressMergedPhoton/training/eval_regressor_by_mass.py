"""Per-mass-point comparison of regressors on one validation pack, in the terms the
mass plots use (Plot/scripts/plot_mlphoton_mass_versions.py): median of
m_rec / m_gen, sigma_eff of m_rec, and f20 = fraction within +-20% of m_gen.

Every model sees the same clusters, so this isolates the regressor (no classifier,
no tagger). Caveat: a pack is resampled flat in log(m/E) inside its training window,
so the energy mix within a mass point is not the analysis one, and the window cuts
the tails of the edge points (M0p1 below m/E 0.002, M1 above the upper cut). And a
model whose training files overlap the pack's val files is favoured -- use the
pack of the model you are testing, so the reference model carries that advantage.

    python eval_regressor_by_mass.py --packed <dir> name=log:<run>/regressor_best.pt \
        name=deployed:<regressor.onnx> [--out table.txt]
"""
import argparse
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import (MLPhotonRegressor, MLPhotonRegressorLog, build_from_onnx,  # noqa: E402
                    load_onnx_weights)
from train_regressor import load_split  # noqa: E402


@torch.no_grad()
def predict(model, img, eta, bs=4096):
    model.eval()
    return torch.cat([model(img[i:i + bs], eta[i:i + bs])
                      for i in range(0, len(img), bs)]).numpy().astype(np.float64)


def sigma_eff(v, frac=0.683):
    v = np.sort(v)
    k = int(np.ceil(frac * len(v)))
    return float((v[k - 1:] - v[:len(v) - k + 1]).min()) / 2.0


def load(spec):
    name, rest = spec.split("=", 1)
    kind, path = rest.split(":", 1)
    if kind == "log":
        model = MLPhotonRegressorLog()
        model.load_state_dict(torch.load(path, map_location="cpu"))
    elif kind == "deployed":
        model = MLPhotonRegressor()
        load_onnx_weights(model, path)
    else:
        raise SystemExit(f"unknown kind {kind} (log | deployed)")
    return name, model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packed", required=True)
    ap.add_argument("--max-mass", type=float, default=1.5)
    ap.add_argument("--out", default=None)
    ap.add_argument("--classifier", default=None,
                    help="classifier.onnx: keep only clusters with diphotonScore > --score-cut, "
                         "as the analysis does (the packs carry every truth-matched cluster, "
                         "and the ones the classifier rejects are badly reconstructed)")
    ap.add_argument("--score-cut", type=float, default=0.5)
    ap.add_argument("models", nargs="+", help="name=log:<.pt> | name=deployed:<.onnx>")
    args = ap.parse_args()

    (img, eta, moe, _), energy, mgen = load_split(args.packed, "val")
    keep = np.ones(len(mgen), dtype=bool)
    if args.classifier:
        # the classifier takes the RAW energy image; load_split hands out the
        # sum-normalised one, so re-read the raw rows (same filter: all val rows
        # here have m/E > 0, checked)
        with np.load(os.path.join(args.packed, "regressor_val.npz")) as z:
            raw = z["image"].astype(np.float32).reshape(-1, 1, 30, 30)
            moe_raw = z["moe_true"]
        assert len(raw) == len(mgen) and np.all(np.isfinite(moe_raw) & (moe_raw > 0))
        cls = build_from_onnx("classifier", args.classifier)
        with torch.no_grad():
            score = torch.cat([cls(torch.from_numpy(raw[i:i + 4096])).exp()[:, 1]
                               for i in range(0, len(raw), 4096)]).numpy()
        keep = score > args.score_cut
        print(f"[cls] {os.path.basename(args.classifier)}: {keep.sum()}/{len(keep)} "
              f"clusters with diphotonScore > {args.score_cut}")
    mg_round = np.where(keep, np.round(mgen, 3), -1.0)
    masses = [m for m in sorted(set(mg_round)) if 0 < m <= args.max_mass]
    models = [load(s) for s in args.models]
    preds = {name: predict(model, img, eta) * energy for name, model in models}

    names = [n for n, _ in models]
    head = f"{'m_gen':>5s} {'N':>6s} | " + " | ".join(
        f"{n:>8s} med/m  s_eff   f20" for n in names)
    sel_txt = (f"; diphotonScore > {args.score_cut} ({os.path.basename(args.classifier)})"
               if args.classifier else "; no classifier selection")
    lines = [f"# {args.packed} (val){sel_txt}; med/m = median(m_rec)/m_gen, s_eff [GeV], "
             f"f20 = |m_rec - m_gen| < 0.2 m_gen", head]
    for m in masses:
        sel = mg_round == m
        cells = []
        for n in names:
            v = preds[n][sel]
            cells.append(f"{'':>8s} {np.median(v) / m:5.3f} {sigma_eff(v):6.3f} "
                         f"{np.mean(np.abs(v - m) < 0.2 * m):5.3f}")
        lines.append(f"{m:5.2f} {int(sel.sum()):6d} | " + " | ".join(cells))
    print("\n".join(lines))
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print("wrote", args.out)


if __name__ == "__main__":
    main()
