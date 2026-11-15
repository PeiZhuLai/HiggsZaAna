#!/usr/bin/env python3
"""Numbers quoted in AN Sec. 10 from the expected limits after the deployment (blind, median expected only):
range of the cross-section limit, its minimum, and the Run 3 vs Run 2 comparison of the cross-section and C_ZH
limits, built with the functions of Plots/makeLimitsPlot.py so that the numbers match the plotted ones.
Input : Combine/output_combine_results_{xsec,wilson}/ (AsymptoticLimits, blind), Plots/run2/{xs,wilson}_limits.json
Output: logs/an_numbers_limits.txt
Usage : an_numbers_limits.py  (inside cmsenv)"""
import os, sys
FF = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit'
sys.path.insert(0, os.path.join(FF, 'Plots'))
import makeLimitsPlot as L

H = os.path.dirname(os.path.abspath(__file__))
M = list(range(1, 31))
xs = L.build_graph_from_combine(M, mode='xs', assume_xs=100.0, higgs_xs=59703.914,
                                results_dir=os.path.join(FF, 'Combine', 'output_combine_results_xsec'))
wi = L.build_graph_from_combine(M, mode='wilson', assume_xs=100.0, higgs_xs=59703.914,
                                results_dir=os.path.join(FF, 'Combine', 'output_combine_results_wilson'))
g = lambda gr: {int(round(gr.GetX()[i])): gr.GetY()[i] for i in range(gr.GetN())}
xs, wi = g(xs[0] if isinstance(xs, tuple) else xs), g(wi[0] if isinstance(wi, tuple) else wi)
r2x = g(L.build_graph_from_hepdata_values(L.read_hepdata_limits_json(L.run2XsLimitsJSON), M, which='exp'))
r2w = g(L.build_graph_from_hepdata_values(L.read_hepdata_limits_json(L.run2WilsonLimitsJSON), M, which='exp'))
out = ['mA | xs [fb] Run3  Run2  (Run3-Run2)/Run2 | C_ZH Run3  Run2  (Run3-Run2)/Run2']
for m in M:
    rx = (xs[m] - r2x[m]) / r2x[m] if m in r2x else float('nan')
    rw = (wi[m] - r2w[m]) / r2w[m] if m in r2w else float('nan')
    out.append('%2d | %6.2f %6.2f %+6.1f%% | %.4g %.4g %+6.1f%%' % (m, xs[m], r2x.get(m, float('nan')), 100 * rx,
                                                              wi[m], r2w.get(m, float('nan')), 100 * rw))
v = [xs[m] for m in M if m >= 2]
best = min(M, key=lambda m: xs[m])
weaker = [m for m in M if m >= 4 and m in r2x and xs[m] > r2x[m]]
out += ['', 'xs range mA>=2: %.2f-%.2f fb; mA1 %.2f fb; best %.2f fb at mA%d' % (min(v), max(v), xs[1], xs[best], best),
        'xs best improvement vs Run2: %.1f%%; Run3 weaker for mA>=4 at %s' % (
            -100 * min((xs[m] - r2x[m]) / r2x[m] for m in M if m in r2x),
            ', '.join('mA%d %+.1f%%' % (m, 100 * (xs[m] - r2x[m]) / r2x[m]) for m in weaker)),
        'C_ZH best improvement vs Run2: %.1f%%' % (-100 * min((wi[m] - r2w[m]) / r2w[m] for m in M if m in r2w))]
txt = '\n'.join(out) + '\n'; open(os.path.join(H, 'logs', 'an_numbers_limits.txt'), 'w').write(txt); print(txt)
