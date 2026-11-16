#!/usr/bin/env python3
"""Paired toy-by-toy comparison of mA = 1 (or 10) bias fits on identical toys (side study, different from production).

The x0 (no signal injected) production-style fits and the side-study fits (selffit: pdfindex frozen to the
truth function; strat1: minimizer strategy 1, tolerance 0.001) use the same seeds, so chunk c / split s holds
the same toys in each. For every truth function this prints, on the toys converged in both, the mean pull of
each configuration (raw mean, not the Gaussian-fit mean of RunBiasStudy.py), and the mean paired difference
with its uncertainty, on the pairs with |pull| < CLIP (default 4) in both fits.
Pull = (r_hat - r_inj) / sigma, sigma on the side facing r_inj (as in RunBiasStudy.py).
Usage: paired_compare.py <x0_root> <side_root> [chunks] [mA]   (roots are BIAS_HS_ROOT directories; mA default 1)
"""
import sys, glob, os, math
import ROOT as r

x0_root, side_root = sys.argv[1], sys.argv[2]
chunks = [int(c) for c in sys.argv[3].split(",")] if len(sys.argv) > 3 else [0, 1]
MA = int(sys.argv[4]) if len(sys.argv) > 4 else 1
R_INJ = 0.0
# Pairs with |pull| >= CLIP in either fit are dropped (RunBiasStudy.py histograms the pull in [-4, 4] before its
# Gaussian fit); a few toys with a tiny interval otherwise dominate the raw mean.
CLIP = float(os.environ.get("CLIP", "4"))


def pulls(path):
    """Return {iToy: (pull, rhat)} for converged toys of one split file."""
    out = {}
    f = r.TFile(path)
    t = f.Get("limit")
    rows = [(e.iToy, e.quantileExpected, e.r) for e in t]
    f.Close()
    by = {}
    for it, q, rv in rows:
        by.setdefault(it, []).append((q, rv))
    for it, v in by.items():
        d = {round(q, 2): rv for q, rv in v}
        if len(v) != 3 or -1.0 not in d or -0.32 not in d or 0.32 not in d:
            continue
        bf, lo, hi = d[-1.0], d[-0.32], d[0.32]
        unc = (hi - bf) if bf < R_INJ else (bf - lo)
        if unc > 0:
            out[it] = ((bf - R_INJ) / unc, bf)
    return out


funcs = {1: ["Bern3", "Exp1", "Lau1", "Pow1"], 10: ["Bern5", "Exp3", "Lau1", "Pow1"]}[MA]
print("%-6s %6s %9s %9s %9s %7s" % ("truth", "nPair", "x0", "side", "side-x0", "err"))
for fn in funcs:
    a_all, b_all = [], []
    for c in chunks:
        for p in sorted(glob.glob("%s/mA_%d/chunk%d/BiasFits/biasStudy_%s_split*_fits.root" % (side_root, MA, c, fn))):
            q = p.replace(side_root, x0_root)
            if not os.path.exists(q):
                continue
            pa, pb = pulls(q), pulls(p)
            for it in sorted(set(pa) & set(pb)):
                if abs(pa[it][0]) >= CLIP or abs(pb[it][0]) >= CLIP:
                    continue
                a_all.append(pa[it][0])
                b_all.append(pb[it][0])
    n = len(a_all)
    if n == 0:
        print("%-6s %6d" % (fn, 0))
        continue
    ma, mb = sum(a_all) / n, sum(b_all) / n
    dd = [b - a for a, b in zip(a_all, b_all)]
    md = sum(dd) / n
    sd = math.sqrt(sum((x - md) ** 2 for x in dd) / (n - 1)) / math.sqrt(n)
    print("%-6s %6d %+9.3f %+9.3f %+9.3f %7.3f" % (fn, n, ma, mb, md, sd))
