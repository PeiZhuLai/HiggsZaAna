#!/usr/bin/env python3
"""Compare the high-statistics bias (Gaussian mean of the pull per truth function) after the deployment with the
pre-deployment results, for the 11 redeployed mass points. Flags |mean| >= 0.2 (criterion |pull| < 0.2) and
truth functions that appear or disappear.
Input : Combine/Checks/Bias_nominal/bias_outputs_highstat/mA_<m>/merged/BiasJson/<m>_gaussfit.json
        Combine/Checks/Bias_nominal/bias_outputs_highstat_variants/mA_<m>_before_deploy_20261111/merged/BiasJson/<m>_gaussfit.json
Output: logs/bias_compare.txt
Usage : compare_bias.py"""
import json, os
B = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs', 'bias_compare.txt')
MASSES = [10, 15, 21, 22, 23, 24, 25, 26, 27, 28, 29]


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


lines = ['mA | r_inj after (before) | truth: mean +- err (sigma) after | before | flag']
flags = []
for m in MASSES:
    a = load('%s/bias_outputs_highstat/mA_%d/merged/BiasJson/%d_gaussfit.json' % (B, m, m))
    b = load('%s/bias_outputs_highstat_variants/mA_%d_before_deploy_20261111/merged/BiasJson/%d_gaussfit.json' % (B, m, m))
    if a is None:
        flags.append('mA%d: no merged json after deployment' % m); lines.append('%2d | MISSING' % m); continue
    fa, fb = a['fit_results'], (b or {}).get('fit_results', {})
    for i, fn in enumerate(sorted(set(fa) | set(fb))):
        x, y = fa.get(fn), fb.get(fn)
        sa = '%+.3f +- %.3f (%.3f)' % (x['mean'], x['mean_err'], x['sigma']) if x else 'not in envelope'
        sb = '%+.3f +- %.3f' % (y['mean'], y['mean_err']) if y else 'not in envelope'
        fl = ''
        if x and abs(x['mean']) >= 0.2: fl = 'BIAS>=0.2'; flags.append('mA%d %s: mean %+.3f' % (m, fn, x['mean']))
        head = '%2d | %.4f (%s)' % (m, a.get('exp', float('nan')), '%.4f' % b['exp'] if b else '-') if i == 0 else '   |' + ' ' * 17
        lines.append('%s | %-6s %s | %s | %s' % (head, fn, sa, sb, fl))
txt = '\n'.join(lines + ['', 'flags (|mean| >= 0.2 or missing): %d' % len(flags)] + ['  ' + f for f in flags]) + '\n'
open(OUT, 'w').write(txt); print(txt)
