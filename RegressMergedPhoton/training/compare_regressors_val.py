"""Compare trained regressors on one validation pack, in the terms the per-mA
ROI analysis uses: for each generated mass, the signal window [q16, q84] of the
reconstructed mass, and how many m/E = 0 single-photon clusters (a stand-in
for the single-photon background) land inside it.

    python compare_regressors_val.py --packed <dir> name=head:path.pt [...]
"""
import argparse
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import MLPhotonRegressorLinear, MLPhotonRegressorLog  # noqa: E402
from train_regressor import load_split  # noqa: E402


@torch.no_grad()
def predict(model, img, eta, bs=4096):
    model.eval()
    return torch.cat([model(img[i:i + bs], eta[i:i + bs])
                      for i in range(0, len(img), bs)]).numpy().astype(np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packed", required=True)
    ap.add_argument("models", nargs="+", help="name=head:path.pt (head = log|linear)")
    args = ap.parse_args()

    (img, eta, moe, _), energy, mgen = load_split(args.packed, "val", keep_zero=True)
    mono = moe.numpy() == 0
    masses = [m for m in sorted(set(np.round(mgen[~mono], 3))) if 0 < m < 1.0]

    for spec in args.models:
        name, rest = spec.split("=", 1)
        head, path = rest.split(":", 1)
        model = MLPhotonRegressorLinear() if head == "linear" else MLPhotonRegressorLog()
        model.load_state_dict(torch.load(path, map_location="cpu"))
        m_rec = predict(model, img, eta) * energy
        print(f"\n=== {name} ({head})   single-photon m_rec median "
              f"{np.median(m_rec[mono]):.3f}")
        print(f"  {'m_gen':>5s} {'median':>7s} {'q16':>7s} {'q84':>7s} {'width':>6s} "
              f"{'mono_in':>8s} {'0.68/sqrt':>9s} {'overlap_next':>12s}")
        win = {}
        for mg in masses:
            sel = np.isclose(mgen, mg) & ~mono
            q16, q50, q84 = np.quantile(m_rec[sel], [0.16, 0.5, 0.84])
            win[mg] = (q16, q84, sel)
        for i, mg in enumerate(masses):
            q16, q84, sel = win[mg]
            f_mono = float(np.mean((m_rec[mono] > q16) & (m_rec[mono] < q84)))
            fom = 0.68 / np.sqrt(f_mono) if f_mono > 0 else float("inf")
            ov = ""
            if i + 1 < len(masses):
                nxt = win[masses[i + 1]][2]
                ov = f"{np.mean((m_rec[nxt] > q16) & (m_rec[nxt] < q84)):.3f}"
            print(f"  {mg:>5.1f} {np.median(m_rec[sel]):>7.3f} {q16:>7.3f} {q84:>7.3f} "
                  f"{q84 - q16:>6.3f} {f_mono:>8.3f} {fom:>9.2f} {ov:>12s}")


if __name__ == "__main__":
    main()
