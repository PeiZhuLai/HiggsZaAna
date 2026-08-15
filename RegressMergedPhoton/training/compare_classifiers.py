"""Compare classifiers the way the analysis actually uses them: ROC, not a fixed cut.

Reading `di_recall(>0.9)` off two different models is misleading, because the
threshold means different things to them. The deployed model pushes almost
everything to the extremes (it labels 76% of clusters hadronic), so a 0.9 cut
keeps a different slice of its output than of a better-calibrated model's. The
honest comparison is: at the SAME hadronic fake rate, which model keeps more
real diphotons?

Prints, for each model, the diphoton efficiency at a set of fixed hadronic fake
rates, plus the ROC AUC over the diphoton-vs-rest problem.

Input : <packed>/classifier_val.npz + one or more models (.onnx or .pt)
Output: stdout table

Usage (CERN, env `hzgdna`):
    python compare_classifiers.py \\
        --onnx ../RecoEgamma/EgammaMLPhotonProducers/data/classifier.onnx \\
        --pt   /eos/.../train_runs/classifier_v2/classifier_best.pt
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import MLPhotonClassifier, load_onnx_weights  # noqa: E402

DEFAULT_PACKED = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/train_packed"
FAKE_RATES = [0.001, 0.005, 0.01, 0.02, 0.0448, 0.10]   # 0.0448 = deployed model at >0.9


def scores(model, img, batch=2048):
    out = []
    model.eval()
    with torch.no_grad():
        for i in range(0, len(img), batch):
            out.append(torch.exp(model(torch.from_numpy(img[i:i + batch]))).numpy())
    return np.concatenate(out)


def eff_at_fake(dip, lab, fake_target):
    """Diphoton efficiency when the threshold is set to give this hadronic fake rate."""
    had = dip[lab == 2]
    di = dip[lab == 1]
    if had.size == 0 or di.size == 0:
        return float("nan"), float("nan")
    thr = np.quantile(had, 1.0 - fake_target)
    return float(np.mean(di > thr)), float(thr)


def auc_di_vs_rest(dip, lab):
    y = (lab == 1).astype(np.int8)
    order = np.argsort(dip)
    y = y[order]
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # rank-based AUC (ties broken by order; fine for a comparison at this precision)
    ranks = np.arange(1, len(y) + 1)
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--packed", default=DEFAULT_PACKED)
    ap.add_argument("--onnx", action="append", default=[], help="repeatable")
    ap.add_argument("--pt", action="append", default=[], help="repeatable")
    args = ap.parse_args()

    with np.load(os.path.join(args.packed, "classifier_val.npz")) as z:
        img = z["image"].astype(np.float32).reshape(-1, 1, 30, 30)
        lab = z["label"]
    print(f"val: {len(lab)} rows  (mono {np.sum(lab==0)}, di {np.sum(lab==1)}, "
          f"had {np.sum(lab==2)})\n")

    models = []
    for p in args.onnx:
        m = MLPhotonClassifier()
        load_onnx_weights(m, p)
        models.append((f"onnx:{os.path.basename(os.path.dirname(p))}/{os.path.basename(p)}", m))
    for p in args.pt:
        m = MLPhotonClassifier()
        m.load_state_dict(torch.load(p, map_location="cpu"))
        models.append((f"pt:{os.path.basename(os.path.dirname(p))}", m))

    if not models:
        print("[FATAL] give at least one --onnx or --pt")
        return 2

    results = {}
    for name, m in models:
        p = scores(m, img)          # one pass; both numbers come out of it
        dip = p[:, 1]
        results[name] = dip
        print(f"=== {name}")
        print(f"    AUC (di vs rest): {auc_di_vs_rest(dip, lab):.4f}")
        print(f"    argmax di recall: "
              f"{np.mean(p.argmax(1)[lab == 1] == 1):.4f}")

    print()
    hdr = f"{'hadronic fake rate':>20s} " + " ".join(f"{n[:22]:>24s}" for n, _ in models)
    print(hdr)
    print("-" * len(hdr))
    for fr in FAKE_RATES:
        cells = []
        for name, _ in models:
            eff, thr = eff_at_fake(results[name], lab, fr)
            cells.append(f"{eff:.4f} (thr {thr:.3f})".rjust(24))
        print(f"{fr:>20.4f} " + " ".join(cells))
    print()
    print("Read this as: at the same hadronic contamination, how many real")
    print("diphotons does each model keep. Higher is better.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
