#!/usr/bin/env python3
"""Truth models for the per-mass range choice: wide-range truth at the 13 narrow-range points whose wide refit is physical
(turn-on below 115 GeV and off the 125 GeV edge). Same construction as bias_wide_truth/make_truth.py.

For each mass point the truth is the function that the production rule min(2NLL + npar) picks after the full-range
refit with the wider ranges (range_full_refit/results/mA<m>.json, entry choice.best_wide), at its refitted values.
The injected signal is r = 0 and r = q50 of the cross-section card (the q50 used by the 11-05 bias study).

Input : ../range_full_refit/results/mA<m>.json, ../q50_at_bias_submission.txt
Output: truth.json  {mA: {pdf, index, set, ranges, q50}}
"""
import json, os, sys
H = os.path.dirname(os.path.abspath(__file__))
MASSES = [3, 4, 5, 8, 9, 11, 12, 13, 14, 16, 17, 18, 20]
q50 = {int(l.split()[0]): float(l.split()[1]) for l in open(os.path.join(H, '..', 'q50_at_bias_submission.txt')) if l.strip()}
out = {}
for m in MASSES:
    d = json.load(open(os.path.join(H, '..', 'range_full_refit', 'results', 'mA%d.json' % m)))
    k = d['choice']['best_wide']; r = [x for x in d['pdfs'] if x['i'] == k][0]
    out[m] = dict(pdf=r['pdf'], index=k, q50=q50[m],
                  set=','.join('%s=%.10g' % kv for kv in r['values_wide'].items()),
                  ranges=':'.join('%s=%g,%g' % (n, a, b) for n, (a, b) in r['ranges_wide'].items()))
    print(m, r['pdf'], k, q50[m], out[m]['set'])
json.dump(out, open(os.path.join(H, 'truth.json'), 'w'), indent=1)
