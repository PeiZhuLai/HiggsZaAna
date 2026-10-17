#!/usr/bin/env python3
"""Diagnostic for the closure pseudo-data: the DY MC selected at a BDT cut, with the closure selection (select_lib,
CLOSURE_DROP_TRAIN honored, CLOSURE_SCORED_BASE picks the model). Prints, per cut: sum w and effective number of MC
events (Neff = (sum w)^2 / sum w^2) in the signal window 120-130 GeV and in the sidebands, the window excess over an
exponential fit to the sidebands (115-135 blinded) in units of the MC statistical error, and the five largest
weights inside the window as a fraction of the window sum.
Env: higgs-alp-ana. Usage: CLOSURE_SCORED_BASE=<dir> CLOSURE_DROP_TRAIN=1 diag_window_mc.py <mA> <cut> [<cut> ...]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "closure_scripts"))
from select_lib import load_combined, select  # noqa: E402


def main():
    mA, cuts = int(sys.argv[1]), [float(c) for c in sys.argv[2:]]
    m, w, s = load_combined(mA)
    for cut in cuts:
        k = select(mA, m, w, s, cut) & (m > 95)
        mm, ww = m[k], w[k]
        win = (mm > 120) & (mm < 130)
        sb = (mm < 115) | (mm > 135)
        edges = np.arange(95, 181, 1.0)
        h, _ = np.histogram(mm[sb], bins=edges, weights=ww[sb])
        c = 0.5 * (edges[:-1] + edges[1:])
        use = ((c < 115) | (c > 135)) & (h > 0)
        p = np.polyfit(c[use], np.log(h[use]), 1, w=np.sqrt(h[use]))
        x = np.arange(120.5, 130, 1.0)
        exp_win = float(np.exp(np.polyval(p, x)).sum())
        sw = float(ww[win].sum()); sw2 = float((ww[win] ** 2).sum())
        neff = sw * sw / sw2 if sw2 > 0 else 0.0
        top = np.sort(ww[win])[::-1][:5] / sw if sw > 0 else []
        print("mA%d cut %.5f | window sumW %.2f Neff %.0f (n=%d) | sideband sumW %.2f Neff %.0f | exp-fit window %.2f "
              "-> excess %.2f = %.1f sigma_MC (R=%.2f) | top5 w/window %s"
              % (mA, cut, sw, neff, win.sum(), ww[sb].sum(), ww[sb].sum() ** 2 / (ww[sb] ** 2).sum(), exp_win,
                 sw - exp_win, (sw - exp_win) / np.sqrt(sw2), sw / exp_win, " ".join("%.3f" % t for t in top)), flush=True)


if __name__ == "__main__":
    main()
