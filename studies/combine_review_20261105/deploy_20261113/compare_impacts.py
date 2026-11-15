#!/usr/bin/env python3
"""Compare the expected impacts (q50, fixidx, prod ranges) after the deployment of mA9, 11, 12, 13, 16, 17, 20 (wide ranges) with the pre-deployment copies.
For every mass point and card: POI fit (r, -1/+1), number of parameters, number of parameters with a failed or
degenerate scan (impact == 0 on both sides, non-finite values, fit at the range edge), and the top-5 ranking.
Input : Combine/output_impacts{,_xsec}_q50_fixidx_prodrange/<m>_impacts.json and the *_bak_deploy_20261113 copies
Output: logs/impacts_compare.txt
Usage : compare_impacts.py"""
import json, math, os
FF = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs', 'impacts_compare.txt')
CARDS = [('wilson', 'output_impacts_q50_fixidx_prodrange'), ('xsec', 'output_impacts_xsec_q50_fixidx_prodrange')]


def load(path):
    if not os.path.exists(path): return None
    d = json.load(open(path))
    poi = d['POIs'][0]['fit']
    pars = {}
    for p in d['params']:
        imp = p.get('impact_r', p.get('impact_' + d['POIs'][0]['name']))
        r = p.get('r', p.get(d['POIs'][0]['name']))
        bad = (r is None or any(not math.isfinite(x) for x in r) or (r[0] == r[1] == r[2]))
        pars[p['name']] = (imp if imp is not None else float('nan'), bad, p.get('fit'))
    return poi, pars


lines, issues = [], []
for card, d in CARDS:
    lines.append('== %s (%s)' % (card, d))
    lines.append('mA | r_hi-r_lo after (before) | npar after/before | failed scans after/before | top5 after | top5 changed')
    for m in (9, 11, 12, 13, 16, 17, 20):
        a = load(os.path.join(FF, d, '%d_impacts.json' % m))
        b = load(os.path.join(FF, d + '_bak_deploy_20261113', '%d_impacts.json' % m))
        if a is None:
            issues.append('%s mA%d: no json after deployment' % (card, m)); lines.append('%2d | MISSING' % m); continue
        (pa, ga), (pb, gb) = a, (b if b else ((float('nan'),) * 3, {}))
        top = lambda g: [n for n, _ in sorted(g.items(), key=lambda kv: -abs(kv[1][0]) if math.isfinite(kv[1][0]) else 0)[:5]]
        fa = sorted(n for n, v in ga.items() if v[1]); fb = sorted(n for n, v in gb.items() if v[1])
        ta, tb = top(ga), top(gb)
        lines.append('%2d | %.4f (%.4f) | %d/%d | %d/%d %s | %s | %s' % (
            m, pa[2] - pa[0], pb[2] - pb[0], len(ga), len(gb), len(fa), len(fb), ','.join(fa)[:80], ' '.join(ta),
            'same' if ta == tb else 'CHANGED (before: %s)' % ' '.join(tb)))
        if fa: issues.append('%s mA%d: %d failed/degenerate scans: %s' % (card, m, len(fa), ', '.join(fa)))
txt = '\n'.join(lines + ['', 'issues: %d' % len(issues)] + ['  ' + i for i in issues]) + '\n'
open(OUT, 'w').write(txt); print(txt)
