#!/usr/bin/env python3
"""High-statistics bias (Gaussian mean of the pull per truth function) with the bias fit profiling the production
parameter ranges of each mass point (bias_fit_config.sh default), compared with the previous runs, in which every mass
point was profiled with the wide ranges. Flags truth functions whose mean exceeds 0.2 in absolute value by more than
one standard deviation (FAIL) or by less (AT).
Input : Combine/Checks/Bias_nominal/bias_outputs_highstat/mA_<m>/merged/BiasJson/<m>_gaussfit.json
        Combine/Checks/Bias_nominal/bias_outputs_highstat_variants/mA_<m>_widefit_20261112/merged/BiasJson/<m>_gaussfit.json
Output: logs/bias_compare_consistent.txt
Usage : compare_bias_consistent.py"""
import json, os
B = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs', 'bias_compare_consistent.txt')
MASSES = list(range(1, 21)) + [22, 23, 24, 25, 26, 29, 30]


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def flag(x):
    if not x: return ''
    ex = abs(x['mean']) - 0.2
    return 'FAIL' if ex > x['mean_err'] else ('AT' if ex > 0 else '')


lines = ['mA | truth | production ranges: mean +- err | previous (wide everywhere) | flag']
for m in MASSES:
    a = load('%s/bias_outputs_highstat/mA_%d/merged/BiasJson/%d_gaussfit.json' % (B, m, m))
    b = load('%s/bias_outputs_highstat_variants/mA_%d_widefit_20261112/merged/BiasJson/%d_gaussfit.json' % (B, m, m))
    if a is None: lines.append('%2d | MISSING' % m); continue
    fa, fb = a['fit_results'], (b or {}).get('fit_results', {})
    for fn in sorted(set(fa) | set(fb)):
        x, y = fa.get(fn), fb.get(fn)
        lines.append('%2d | %-6s | %s | %s | %s' % (
            m, fn, ('%+.3f +- %.3f' % (x['mean'], x['mean_err'])) if x else 'not in envelope',
            ('%+.3f +- %.3f' % (y['mean'], y['mean_err'])) if y else 'missing', flag(x)))
open(OUT, 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
