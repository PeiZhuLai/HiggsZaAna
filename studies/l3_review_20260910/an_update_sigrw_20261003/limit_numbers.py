"""Expected limits (fb) per mA, +-1/2 sigma, and Run3-vs-Run2 ratios (xs, Wilson), reusing makeLimitsPlot.py helpers.
usage: limit_numbers.py <combine_results_dir>
Input: <dir>/higgsCombine<m>.AsymptoticLimits.mH125.38.root, Plots/run2/*.json. Output: stdout"""
import sys, os
P = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Plots"
sys.path.insert(0, P); os.chdir(P)
import makeLimitsPlot as M
R = sys.argv[1]
masses = list(range(1, 31))
jx = M.read_hepdata_limits_json(M.run2XsLimitsJSON); jw = M.read_hepdata_limits_json(M.run2WilsonLimitsJSON)
def tomap(g): return {int(round(g.GetPointX(i))): g.GetPointY(i) for i in range(g.GetN())}
r2x = tomap(M.build_graph_from_hepdata_values(jx, masses, "exp")); r2w = tomap(M.build_graph_from_hepdata_values(jw, masses, "exp"))
print("mA  exp[fb]  -2s  -1s  +1s  +2s | run2xs  ratio-1 | wilson3 wilson2 ratio-1")
for m in masses:
    ok, q025, q16, q50, q84, q975, obs, hasObs = M.read_limits_from_file(f"{R}/higgsCombine{m}.AsymptoticLimits.mH125.38.root")
    x = q50 * 100
    br = x / (59703.914 * float(M.ZToll_br)); w3 = float(M.calc_Wilson_coupling(br, m))
    rx = (x / r2x[m] - 1) if m in r2x else float('nan')
    rw = (w3 / r2w[m] - 1) if m in r2w else float('nan')
    print("%2d %7.2f %6.2f %6.2f %6.2f %6.2f | %6.2f %+6.3f | %.4g %.4g %+6.3f" % (m, x, q025*100, q16*100, q84*100, q975*100, r2x.get(m, float('nan')), rx, w3, r2w.get(m, float('nan')), rw))
