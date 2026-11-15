#!/usr/bin/env python3
"""Table of refit_limit.py results.
Input : results/mA<m>.json, limits/mA<m>.txt, ../impact_shape_test_mA10/checks/all_masses/all_masses.txt (profile-only ratio)
Output: summary.txt
Columns: q50 ratios to prod for refitP (closure), wideP (production pdfindex), wideB (pdfindex re-chosen), and the earlier
profile-only widening; function of each pdfindex; turn-on/width of the re-chosen function after the wide refit; number
of parameters at a bound before/after; sideband S (no 115-135 GeV bin) of the re-chosen function, stored -> wide refit;
fits with a nonzero status (wide refits)."""
import json, os
H = os.path.dirname(os.path.abspath(__file__))
prof = {}
for l in open(os.path.join(H, '..', 'impact_shape_test_mA10', 'checks', 'all_masses', 'all_masses.txt')):
    if l[0].isdigit(): prof[int(l.split()[0])] = float(l.split()[-1])
out = ['mA  q50_prod  refitP  wideP  wideB  profile | idx_prod idx_wide | turnon width (idx_wide) | bound st->wide | S_sb(idx_wide) st->wide | st!=0']
for m in range(1, 31):
    fj, fl = os.path.join(H, 'results', 'mA%d.json' % m), os.path.join(H, 'limits', 'mA%d.txt' % m)
    if not (os.path.exists(fj) and os.path.exists(fl)): out.append('%-3d MISSING' % m); continue
    d = json.load(open(fj)); L = {}
    for l in open(fl).read().splitlines()[1:]:
        p = l.split(); L[p[0]] = [float(x) for x in p[1:6]]
    P = {r['i']: r for r in d['pdfs']}; c = d['choice']; b = P[c['best_wide']]
    to = [v for k, v in b['values_wide'].items() if 'turnon' in k or k == 'Bern_step']
    wi = [v for k, v in b['values_wide'].items() if 'width' in k or k == 'Bern_stepWidth']
    r = lambda k: L[k][2]/L['prod'][2]
    out.append('%-3d %.5f  %.3f  %.3f  %.3f  %.3f | %-6s %-6s | %6.1f %5.2f | %3d -> %3d | %6.2f -> %6.2f | %d' % (
        m, L['prod'][2], r('refitP'), r('wideP'), r('wideB'), prof.get(m, float('nan')), P[c['prod_index']]['pdf'], b['pdf'],
        to[0] if to else float('nan'), wi[0] if wi else float('nan'),
        sum(len(x['bound_stored']) for x in d['pdfs']), sum(len(x['bound_wide']) for x in d['pdfs']),
        b['S_sb_stored'], b['S_sb_wide'], sum(x['status_wide'] != 0 for x in d['pdfs'])))
    if c['best_stored'] != c['prod_index']: out[-1] += '  [best_stored %s != prod]' % P[c['best_stored']]['pdf']
open(os.path.join(H, 'summary.txt'), 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
