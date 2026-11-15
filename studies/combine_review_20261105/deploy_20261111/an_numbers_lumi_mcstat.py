#!/usr/bin/env python3
"""Numbers for AN Table tab:uncertainty_summary after the deployment: luminosity (LUM POG 3-nuisance scheme) and
signal MC statistics, as signal-yield-weighted averages over the ten era x channel processes (the convention of the
table), plus the effect on the total signal yield (nuisances added coherently over processes, then in quadrature).
Input : flashggFinalFit/Datacard/yields_leptons/<m>_pruned_Datacard_leptons.pkl (nominal_yield per process)
        studies/combine_review_20261105/signal_mcstat/mcstat.json (mcstat value per mass and process)
Output: logs/an_numbers_lumi_mcstat.txt
Usage : an_numbers_lumi_mcstat.py"""
import json, math, os, pickle
H = os.path.dirname(os.path.abspath(__file__))
Y = '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/yields_leptons'
MC = json.load(open(os.path.join(H, '..', 'signal_mcstat', 'mcstat.json')))
YEAR = {'2022preEE': '2022', '2022postEE': '2022', '2023preBPix': '2023', '2023postBPix': '2023', '2024': '2024'}
# LUM POG 3-nuisance scheme (lnN - 1, per year), as written in the cards
LUMI = {'lumi_13p6TeV_222324': {'2022': .0138, '2023': .0017, '2024': .0020},
        'lumi_13p6TeV_2324': {'2023': .0127, '2024': .0068},
        'lumi_13p6TeV_2024': {'2024': .0144}}
rows, out = [], []
for m in range(1, 31):
    d = pickle.load(open(os.path.join(Y, '%d_pruned_Datacard_leptons.pkl' % m), 'rb'))
    d = d[d['type'] == 'sig']
    y = {r['proc']: float(r['nominal_yield']) for _, r in d.iterrows()}
    tot = sum(y.values()); w = {p: v / tot for p, v in y.items()}
    yr = {p: YEAR[p.split('_')[1]] for p in y}
    per_proc_lumi = {p: math.sqrt(sum(L.get(yr[p], 0) ** 2 for L in LUMI.values())) for p in y}
    lumi_avg = sum(w[p] * per_proc_lumi[p] for p in y)
    lumi_tot = math.sqrt(sum(sum(w[p] * L.get(yr[p], 0) for p in y) ** 2 for L in LUMI.values()))
    mc_avg = sum(w[p] * MC[str(m)][p] for p in y)
    mc_tot = math.sqrt(sum((w[p] * MC[str(m)][p]) ** 2 for p in y))
    rows.append((m, lumi_avg, lumi_tot, mc_avg, mc_tot))
    out.append('%2d | lumi avg %.2f%% total %.2f%% | mcstat avg %.2f%% total %.2f%%' % (m, 100 * lumi_avg, 100 * lumi_tot, 100 * mc_avg, 100 * mc_tot))
rng = lambda i: '%.2f-%.2f%% (mean %.2f%%)' % (100 * min(r[i] for r in rows), 100 * max(r[i] for r in rows), 100 * sum(r[i] for r in rows) / len(rows))
out += ['', 'lumi   yield-weighted average of per-process totals: ' + rng(1), 'lumi   effect on the total yield: ' + rng(2),
        'mcstat yield-weighted average of per-process values: ' + rng(3), 'mcstat effect on the total yield: ' + rng(4)]
txt = '\n'.join(out) + '\n'; open(os.path.join(H, 'logs', 'an_numbers_lumi_mcstat.txt'), 'w').write(txt); print(txt)
