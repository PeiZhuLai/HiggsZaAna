#!/usr/bin/env python3
"""Per-mass summary of fit_quality.py next to the limit comparison.
Input : results/mA<m>.json, ../impact_shape_test_mA10/checks/all_masses/all_masses.txt
Output: summary.txt, one line per function, then one line per mass with its largest dS"""
import glob, json
lim = {}
for l in open('../impact_shape_test_mA10/checks/all_masses/all_masses.txt'):
    if l[0].isdigit(): lim[int(l.split()[0])] = float(l.split('|')[2])
recs = [r for f in glob.glob('results/mA*.json') for r in json.load(open(f))]; recs.sort(key=lambda r: (r['mA'], r['S_prod']))
step = lambda r, k: [v for n, v in r[k].items() if 'turnon' in n or n == 'Bern_step'][0]
soft = lambda r, k: [v for n, v in r[k].items() if 'width' in n or n == 'Bern_stepWidth'][0]
out = ['mA  pdf    npar nbins  S_prod  S_wide     dS  P/ndf_prod P/ndf_wide  turn-on prod->wide  width prod->wide  bound prod->wide st sr_ratio']
for r in recs:
    ndf = r['nbins'] - r['npar'] - 1
    out.append('%-3d %-6s %3d %5d %7.2f %7.2f %6.2f %9.2f %10.2f %9.1f %7.1f %7.1f %7.1f %8d %3d   %d/%d %6.3f' % (r['mA'], r['pdf'], r['npar'], r['nbins'], r['S_prod'], r['S_wide'], r['dS'], r['P_prod']/ndf, r['P_wide']/ndf, step(r, 'values_prod'), step(r, 'values_wide'), soft(r, 'values_prod'), soft(r, 'values_wide'), len(r['bound_prod']), len(r['bound_wide']), r['status_prod'], r['status_wide'], r['sr_ratio']))
out += ['', 'mA  limit_ratio  n_pdf  best_prod  (S_prod S_wide dS sr_ratio)   largest dS  (pdf)   not_nested  fit_status!=0']
for m in sorted(lim):
    rs = [r for r in recs if r['mA'] == m]
    if not rs: continue
    b = min(rs, key=lambda r: r['S_prod'] + r['npar']); w = max(rs, key=lambda r: r['dS'])
    out.append('%-3d %6.3f %6d   %-6s %7.2f %7.2f %6.2f %6.3f   %6.2f %-6s %3d %3d' % (m, lim[m], len(rs), b['pdf'], b['S_prod'], b['S_wide'], b['dS'], b['sr_ratio'], w['dS'], w['pdf'], sum(not r['nested'] for r in rs), sum(r['status_prod'] != 0 or r['status_wide'] != 0 for r in rs)))
open('summary.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
