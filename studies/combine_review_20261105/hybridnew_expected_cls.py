#!/usr/bin/env python3
"""Expected CLs limits from the HybridNew grid, read directly from the stored HypoTestResults.

For every grid point r and quantile q of the limit: the test statistic of the b-only toys at
quantile 1-q is taken as the 'observed' value (a larger q_LHC is less signal-like), and CLs(r) and its toy uncertainty are computed. The limit is where
CLs = 0.05, from a weighted straight-line fit of ln CLs versus r over the five grid points closest
to the crossing (uncertainty from the fit covariance). Compared with the asymptotic expected limits.
Usage: python3 hybridnew_expected_cls.py [subdir]   (default hybridnew_v2); run inside cmsenv.
"""
import os, re, sys
import numpy as np
import ROOT

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
FF = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
SUB = sys.argv[1] if len(sys.argv) > 1 else "hybridnew_v2"
HN = "%s/Combine/Checks/review_20261105/%s" % (FF, SUB)
W = os.path.dirname(os.path.abspath(__file__))
QS = (0.025, 0.16, 0.5, 0.84, 0.975)


def grid(m):
    f = ROOT.TFile.Open("%s/grid_mA%d.root" % (HN, m))
    res = {}
    for k in f.Get("toys").GetListOfKeys():
        mt = re.search(r"_r([0-9.]+)_", k.GetName() + "_")
        if not mt:
            continue
        r = float(mt.group(1)); o = k.ReadObj()
        if r in res:
            res[r].Append(o)
        else:
            res[r] = o.Clone()
    return f, res


def limit(points):
    rs = np.array([p[0] for p in points]); cl = np.array([p[1] for p in points]); er = np.array([p[2] for p in points])
    ok = (cl > 0) & (er > 0)
    rs, cl, er = rs[ok], cl[ok], er[ok]
    i = int(np.argmin(np.abs(np.log(cl) - np.log(0.05))))
    sel = slice(max(0, i - 2), min(len(rs), i + 3))
    x, y, w = rs[sel], np.log(cl[sel]), (cl[sel] / er[sel]) ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    a, b = cov @ (A.T @ (w * y))
    r0 = (np.log(0.05) - a) / b
    g = np.array([-1 / b, -(np.log(0.05) - a) / b ** 2])
    return r0, float(np.sqrt(g @ cov @ g))


def main():
    out = ["mA  quantile  toys            asymptotic  ratio"]
    for m in (5, 20, 30):
        f, res = grid(m)
        fa = ROOT.TFile.Open("%s/Combine/output_combine_results_xsec/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (FF, m))
        asy = [e.limit for e in fa.Get("limit")]
        for iq, q in enumerate(QS):
            pts = []
            for r in sorted(res):
                h = res[r]
                b = np.array(list(h.GetNullDistribution().GetSamplingDistribution()))
                # a larger q_LHC is less signal-like: the q quantile of the limit is the (1-q) quantile of q_LHC
                h.SetTestStatisticData(float(np.quantile(b, 1 - q)))
                pts.append((r, h.CLs(), h.CLsError()))
            r0, e0 = limit(pts)
            out.append("%-3d %-8s  %.5f +- %.5f  %.5f     %.3f +- %.3f" % (m, q, r0, e0, asy[iq], r0 / asy[iq], e0 / asy[iq]))
        n = [len(list(res[r].GetNullDistribution().GetSamplingDistribution())) for r in res]
        out.append("    (mA%d: %d grid points, %d-%d b-only toys per point)" % (m, len(res), min(n), max(n)))
    open("%s/%s_expected_cls.txt" % (W, SUB), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
