#!/usr/bin/env python3
"""Recompute the bias-study pulls from the merged high-stat fit trees with several definitions.

RunBiasStudy.py used a symmetric uncertainty until 2026-09-26, pull = (r_fit - r_true) / (0.5*(hi-lo)),
and quotes the mean of a Gaussian fit. Near the expected limit the r uncertainty is asymmetric,
so this checks how much of the common negative offset seen at every mass point comes from the
pull definition rather than from the envelope.

  sym       : (bf - mu) / (0.5*(hi - lo))                          [flashggFinalFit upstream]
  toward    : error on the side facing the truth: bf < mu -> hi - bf, bf > mu -> bf - lo  [Combine tutorial plot_bias_pull.py; used since 2026-09-26]
  away      : the opposite side (what the commented-out block in RunBiasStudy.py would do)

For each: Gaussian-fit mean in [-4,4] (80 bins, as RunBiasStudy), plain mean, and median.
Input : Combine/Checks/Bias_nominal/bias_outputs_highstat/mA_<m>/merged/BiasFits/biasStudy_<fn>_fits.root
        injected mu from .../merged/BiasJson/<m>_gaussfit.json ("exp")
Output: stdout table (and logs_fsrfix/bias_pull_definitions.txt when run via the wrapper)
"""
import glob, json, os, sys
import numpy as np
import uproot

BN = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal"
masses = [int(x) for x in sys.argv[1:]] or list(range(1, 31))


def triplets(path):
    t = uproot.open(path)["limit"]
    r = t["r"].array(library="np")
    q = t["quantileExpected"].array(library="np")
    out = []
    i = 0
    while i + 2 < len(r):  # same sequential walk as RunBiasStudy, but resync on a broken triplet
        if q[i] == -1 and abs(q[i + 1] + 0.32) < 1e-3 and abs(q[i + 2] - 0.32) < 1e-3:
            out.append((r[i], r[i + 1], r[i + 2]))
            i += 3
        else:
            i += 1
    return np.array(out)


def gaus_mean(p):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    h = ROOT.TH1F("h%d" % np.random.randint(1e9), "", 80, -4, 4)
    for v in p:
        h.Fill(v)
    h.Fit("gaus", "Q0")
    return h.GetFunction("gaus").GetParameter(1)


print("%3s %-6s %6s | %8s %8s %8s | %8s %8s | %8s | %s" % (
    "mA", "truth", "ntoy", "sym_gaus", "sym_mean", "sym_med", "tow_gaus", "tow_med", "away_gaus", "asym(hi-bf)/(bf-lo)"))
for m in masses:
    d = "%s/bias_outputs_highstat/mA_%d/merged" % (BN, m)
    j = "%s/BiasJson/%d_gaussfit.json" % (d, m)
    if not os.path.exists(j):
        continue
    mu = json.load(open(j))["exp"]
    for f in sorted(glob.glob("%s/BiasFits/biasStudy_*_fits.root" % d)):
        if "_split" in f:
            continue
        fn = os.path.basename(f)[len("biasStudy_"):-len("_fits.root")]
        a = triplets(f)
        bf, lo, hi = a[:, 0], a[:, 1], a[:, 2]
        ok = (hi > bf) & (bf > lo)
        bf, lo, hi = bf[ok], lo[ok], hi[ok]
        diff = bf - mu
        sym = diff / (0.5 * (hi - lo))
        toward = diff / np.where(bf < mu, hi - bf, bf - lo)
        away = diff / np.where(bf < mu, bf - lo, hi - bf)
        asym = np.median((hi - bf) / (bf - lo))
        print("%3d %-6s %6d | %+8.3f %+8.3f %+8.3f | %+8.3f %+8.3f | %+8.3f | %.2f" % (
            m, fn, len(bf), gaus_mean(sym), sym.mean(), np.median(sym),
            gaus_mean(toward), np.median(toward), gaus_mean(away), asym))
