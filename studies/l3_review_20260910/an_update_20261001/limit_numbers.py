"""Expected limits (fb) per mA from the current combine outputs, +-1/2 sigma, and the
Run3-vs-Run2 ratios for xs and Wilson, reusing makeLimitsPlot.py helpers.
Input: Combine/output_combine_results/higgsCombine<m>.AsymptoticLimits.mH125.38.root, Plots/run2/*.json
Output: stdout"""
import sys, os
P = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Plots"
sys.path.insert(0, P); os.chdir(P)
import makeLimitsPlot as M
masses = list(range(1, 31))
R = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results"
jx = M.read_hepdata_limits_json(M.run2XsLimitsJSON); jw = M.read_hepdata_limits_json(M.run2WilsonLimitsJSON)
def tomap(g): return {int(round(g.GetPointX(i))): g.GetPointY(i) for i in range(g.GetN())}
r2x = tomap(M.build_graph_from_hepdata_values(jx, masses, "exp")); r2w = tomap(M.build_graph_from_hepdata_values(jw, masses, "exp"))
g3w, _ = M.build_graph_from_combine(masses, mode="wilson"); r3w = tomap(g3w)
print("mA  exp[fb]  -2s  -1s  +1s  +2s | run2xs  ratio-1 | wilson3 wilson2 ratio-1")
for m in masses:
    ok, q025, q16, q50, q84, q975, obs, hasObs = M.read_limits_from_file(f"{R}/higgsCombine{m}.AsymptoticLimits.mH125.38.root")
    x = q50 * 100
    rx = (x / r2x[m] - 1) if m in r2x else float('nan')
    rw = (r3w[m] / r2w[m] - 1) if m in r2w and m in r3w else float('nan')
    print("%2d %7.2f %6.2f %6.2f %6.2f %6.2f | %6.2f %+6.3f | %.4g %.4g %+6.3f" % (m, x, q025*100, q16*100, q84*100, q975*100, r2x.get(m, float('nan')), rx, r3w.get(m, float('nan')), r2w.get(m, float('nan')), rw))
