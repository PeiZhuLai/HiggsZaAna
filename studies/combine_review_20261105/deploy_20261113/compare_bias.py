#!/usr/bin/env python3
"""High-statistics bias (Gaussian mean of the pull per truth function) of mA9, 11, 12, 13, 16, 17, 20 after the
deployment of the wide parameter ranges, compared with
  narrow : per-mass windows, bias fit with the production ranges (run of 11-12, moved to the variants directory)
  wide   : per-mass windows in the truth, bias fit with the wide ranges (run of 11-12)
The truth functions differ between the runs because the envelope was rebuilt. Flags |mean| above 0.2 by more than one
standard deviation (FAIL) or by less (AT).
Input : Combine/Checks/Bias_nominal/bias_outputs_highstat/mA_<m>/merged/BiasJson/<m>_gaussfit.json
        Combine/Checks/Bias_nominal/bias_outputs_highstat_variants/mA_<m>_{narrow,widefit}_20261112/merged/BiasJson/
Output: logs/bias_compare.txt
Usage : compare_bias.py"""
import json, os
B = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs', 'bias_compare.txt')
MASSES = (9, 11, 12, 13, 16, 17, 20)


def load(p):
    return json.load(open(p))['fit_results'] if os.path.exists(p) else None


def flag(x):
    ex = abs(x['mean']) - 0.2
    return 'FAIL' if ex > x['mean_err'] else ('AT' if ex > 0 else '')


def fmt(fr):
    if fr is None: return 'missing'
    return ', '.join('%s %+.3f' % (k, v['mean']) for k, v in sorted(fr.items()))


lines, n_fail, n_at = [], 0, 0
for m in MASSES:
    a = load('%s/bias_outputs_highstat/mA_%d/merged/BiasJson/%d_gaussfit.json' % (B, m, m))
    lines.append('== mA%d' % m)
    if a is None: lines.append('  after deployment: MISSING'); n_fail += 1; continue
    for fn, x in sorted(a.items()):
        f = flag(x); n_fail += f == 'FAIL'; n_at += f == 'AT'
        lines.append('  after  %-6s %+.3f +- %.3f  sigma %.3f +- %.3f  %s' % (fn, x['mean'], x['mean_err'],
                                                                             x.get('sigma', float('nan')),
                                                                             x.get('sigma_err', float('nan')), f))
    for tag in ('narrow', 'widefit'):
        lines.append('  %-6s %s' % (tag, fmt(load('%s/bias_outputs_highstat_variants/mA_%d_%s_20261112/merged/BiasJson/'
                                                  '%d_gaussfit.json' % (B, m, tag, m)))))
lines += ['', 'FAIL %d, AT %d' % (n_fail, n_at)]
open(OUT, 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
