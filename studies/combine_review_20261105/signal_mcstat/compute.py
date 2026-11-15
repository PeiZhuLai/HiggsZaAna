#!/usr/bin/env python3
"""Relative signal MC statistical uncertainty per mass point, era and lepton channel: sqrt(sumw2) / nominal_yield
of makeYields.py (test tree after the BDT working point; the 1/TEST_TREE_FRAC scaling cancels in the ratio).
Interpolated mass points (makeYields.interploate_ma_list) take their yield from an anchor mass; for them the table
gives the anchor used (anchor_mass column of the pickle) and the larger value of the two bracketing anchors.
Input : flashggFinalFit/Datacard/yields_leptons/<m>_pruned_Datacard_leptons.pkl (production, read only)
Output: mcstat.txt (per mass x process), mcstat.json {m: {proc: value used}}"""
import json, math, os, pickle
H = os.path.dirname(os.path.abspath(__file__))
Y = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/yields_leptons'
INTERP = [11, 12, 13, 14, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29]
ANCH = [m for m in range(1, 31) if m not in INTERP]
def load(m):
    d = pickle.load(open(os.path.join(Y, '%d_pruned_Datacard_leptons.pkl' % m), 'rb'))
    d = d[d['type'] == 'sig']
    return {r['proc']: (float(r['nominal_yield']), float(r['sumw2']), int(r['anchor_mass'])) for _, r in d.iterrows()}
rel = lambda y, w2: math.sqrt(w2) / y if y > 0 else float('nan')
procs = None; out = []; used = {}; anchors = {m: load(m) for m in ANCH}
for m in range(1, 31):
    rows = load(m); procs = procs or list(rows)
    lo = max(a for a in ANCH if a <= m); hi = min(a for a in ANCH if a >= m)
    used[m] = {}
    for p in procs:
        y, w2, am = rows[p]; r = rel(y, w2)
        if m in INTERP:
            ra = {a: rel(anchors[a][p][0], anchors[a][p][1]) for a in (lo, hi)}
            v = max(ra.values())
            out.append('%-3d %-22s own/anchor(mA%d) %.4f | anchors mA%d %.4f, mA%d %.4f -> max %.4f' % (m, p, am, r, lo, ra[lo], hi, ra[hi], v))
        else:
            v = r; out.append('%-3d %-22s %.4f' % (m, p, r))
        used[m][p] = round(v, 5)
# compact table: value used, percent
LAB = {'2022preEE': '2022', '2022postEE': '2022EE', '2023preBPix': '2023', '2023postBPix': '2023BPix', '2024': '2024'}
lab = lambda p: '%s_%s' % (LAB[p.split('_')[1]], p.split('_')[2])
hdr = 'mA  ' + ' '.join('%9s' % lab(p) for p in procs) + '   max'
tab = [hdr] + ['%-3d ' % m + ' '.join('%8.2f%%' % (100 * used[m][p]) for p in procs) + '  %5.2f%%' % (100 * max(used[m].values())) + ('  (interp)' if m in INTERP else '') for m in range(1, 31)]
open(os.path.join(H, 'mcstat.txt'), 'w').write('\n'.join(tab) + '\n\n' + '\n'.join(out) + '\n')
json.dump(used, open(os.path.join(H, 'mcstat.json'), 'w'), indent=1)
print('\n'.join(tab))
