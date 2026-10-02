#!/usr/bin/env python3
"""Run 3 signal preselection efficiency after the FSR-fix re-production (input to AN Eqn. 33).

Three estimates, all per (mA anchor, era), then luminosity-weighted over the 5 eras:
 (1) unw_fsrfix : sum n_events_selected['nominal'] / sum n_events, from the FSR-fix staging summaries
                  /eos/project/.../parquet_DNA_tmp_fsrfix/Sig_MC/mA_M<m>_<era>/job_*/*_summary_job*.json
                  (same definition as the unweighted cutflow 'all cuts'/'all').
 (2) cf_prefsr  : pre-FSR-fix cutflow JSONs used for the current AN Fig. (cutflowVmA),
                  HiggsDNA/cutflow/cutflow_list/cutflow_Sig_MC_mA_M<m>_<era>.json, zgammas_w and zgammas.
 (3) sf_fsrfix  : sum(weight_central_no_lumi)/sigma[fb] over the 'inclusive' tree of the FSR-fix scored
                  ROOT (includes all object scale factors); sigma = 0.1 pb = 100 fb. Pre-FSR-fix
                  (run3_bdt_scored_nominal) is computed the same way for comparison.
"""
import json, glob, os, re, statistics as st
import uproot, numpy as np

ERAS = ['2022preEE', '2022postEE', '2023preBPix', '2023postBPix', '2024']
LUMI = {'2022preEE': 7.99, '2022postEE': 26.68, '2023preBPix': 17.96, '2023postBPix': 9.68, '2024': 109.82}
ANCH = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
STG = '/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC'
CF = '/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/cutflow/cutflow_list'
SC = '/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_{}'

res = {}
for m in ANCH:
    for e in ERAS:
        r = {}
        ne = ns = 0.0; nj = 0
        for f in glob.glob(f'{STG}/mA_M{m}_{e}/job_*/*_summary_job*.json'):
            j = json.load(open(f)); ne += j['n_events']; ns += j['n_events_selected'].get('nominal', 0); nj += 1
        r['unw_fsrfix'] = 100 * ns / ne if ne else float('nan'); r['njobs'] = nj; r['n_events'] = ne
        c = json.load(open(f'{CF}/cutflow_Sig_MC_mA_M{m}_{e}.json'))['cutflows']
        r['cfw_prefsr'] = 100 * c['zgammas_w']['all cuts'] / c['zgammas_w']['all']
        r['cfu_prefsr'] = 100 * c['zgammas']['all cuts'] / c['zgammas']['all']
        for tag in ['fsrfix', 'nominal']:
            a = uproot.open(f'{SC.format(tag)}/mA_M{m}/{e}.root')['inclusive'].arrays(['weight_central_no_lumi'], library='np')
            r[f'sf_{tag}'] = float(a['weight_central_no_lumi'].sum()) / 100.0 * 100  # [%]
        res[f'{m}_{e}'] = r
        print(m, e, ' '.join(f'{k}={v:.3f}' if isinstance(v, float) else f'{k}={v}' for k, v in r.items()), flush=True)

L = sum(LUMI.values())
print('\n# lumi-weighted over eras, per mA  [%]')
keys = ['unw_fsrfix', 'cfu_prefsr', 'cfw_prefsr', 'sf_nominal', 'sf_fsrfix']
summ = {k: [] for k in keys}
for m in ANCH:
    row = {k: sum(LUMI[e] * res[f'{m}_{e}'][k] for e in ERAS) / L for k in keys}
    for k in keys:
        summ[k].append(row[k])
    print(m, ' '.join(f'{k}={row[k]:.2f}' for k in keys))
print('\n# mean over 14 anchors (and range)')
for k in keys:
    print(f'{k:12s} mean={st.mean(summ[k]):.2f}  range={min(summ[k]):.2f}--{max(summ[k]):.2f}  mean(mA<=20)={st.mean(summ[k][:12]):.2f}')
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'presel_eff_fsrfix.json'), 'w'), indent=1)
