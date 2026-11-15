#!/usr/bin/env python3
"""Summary of the Pow1 truth diagnostic (run.sh): per mass, pull mean/median of the envelope fit and of the
Pow1-only fit, and for the envelope fit the fraction of toys per selected function with the pull mean in each group.
Pull as in RunBiasStudy.py: (r_hat - r_inj) / (hi - r_hat) if r_hat < r_inj else / (r_hat - lo).
Input : out/mA<m>/{env,pow1}.root
Output: logs/summary.txt
Usage : summarize.py"""
import os, statistics as st
import ROOT
D = os.path.dirname(os.path.abspath(__file__))
CASES = {9: (0.0168, ['Pow1', 'Exp3', 'Lau1', 'Bern6']), 11: (0.0223, ['Lau1', 'Exp1', 'Pow1', 'Bern5']),
         12: (0.0238, ['Lau1', 'Exp1', 'Bern5', 'Pow1']), 17: (0.0285, ['Exp1', 'Pow1', 'Lau1', 'Bern5'])}
IDX = 'pdfindex_cat0_13p6TeV'


def pulls(path, r_inj, with_idx):
    f = ROOT.TFile.Open(path); t = f.Get('limit'); toys = {}
    for e in t:
        d = toys.setdefault(e.iToy, {})
        q = round(e.quantileExpected, 3)
        d[{-1.0: 'bf', -0.32: 'lo', 0.32: 'hi'}.get(q, q)] = e.r
        if with_idx and q == -1.0: d['idx'] = int(getattr(e, IDX))
    out = []
    for d in toys.values():
        if not all(k in d for k in ('bf', 'lo', 'hi')): continue
        bf = d['bf']; unc = (d['hi'] - bf) if bf < r_inj else (bf - d['lo'])
        if unc <= 0: continue
        out.append(((bf - r_inj) / unc, d.get('idx')))
    return out


lines = []
for m, (r, names) in CASES.items():
    o = '%s/out/mA%d' % (D, m)
    env, pw = pulls(o + '/env.root', r, True), pulls(o + '/pow1.root', r, False)
    pe, pp = [x for x, _ in env], [x for x, _ in pw]
    lines.append('mA%d  env: n %d mean %+.3f +- %.3f median %+.3f | Pow1 only: n %d mean %+.3f +- %.3f median %+.3f' % (
        m, len(pe), st.mean(pe), st.stdev(pe) / len(pe) ** .5, st.median(pe),
        len(pp), st.mean(pp), st.stdev(pp) / len(pp) ** .5, st.median(pp)))
    for i, n in enumerate(names):
        g = [x for x, j in env if j == i]
        if g: lines.append('      env picks %-5s %5.1f%%  pull mean %+.3f' % (n, 100. * len(g) / len(env), st.mean(g)))
open(D + '/logs/summary.txt', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
