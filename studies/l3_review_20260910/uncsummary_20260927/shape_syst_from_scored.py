#!/usr/bin/env python3
"""Physical effect of the shape systematics (FNUF, Material, e/mu/gamma scale & smear) on the
signal yield, peak position and resolution, from the FSR-fix BDT-scored signal ROOT files.

Input : /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M<m>{,_<Syst>_<up|down>}/<era>.root
        tree 'inclusive' (all events; the syst samples are the same events with shifted energies,
        so up/nominal ratios are not limited by the 30% test split used in the signal model).
Selection: MVA_Score_mA_M<m> >= MVAcut (Plot/output/MVAcut_points_run3.json), 95 < H_mass < 180.
Per (mA, era, channel, syst):
  rate   = mean(|N_up/N_nom - 1|, |N_dn/N_nom - 1|)                (also without MVA cut: rate_presel)
  dpeak  = mean(|mean_up - mean_nom|, |mean_dn - mean_nom|)  [GeV]   (weighted mean of H_mass in 100-180)
  dres   = mean(|seff_up - seff_nom|, |seff_dn - seff_nom|)  [GeV]   (sigma_eff: smallest 68.3% interval)
usage: python3 shape_syst_from_scored.py <out.json> <mA,mA,...>
"""
import sys, json, numpy as np, uproot

BASE = '/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix'
ERAS = ['2022preEE', '2022postEE', '2023preBPix', '2023postBPix', '2024']
SYSTS = ['FNUF', 'Material', 'Electron_scale', 'Electron_smear', 'Muon_scale', 'Muon_smear',
         'Photon_scale', 'Photon_smear']
CUTS = {r['mA']: r['MVAcut'] for r in json.load(open(
    '/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json'))['results']}


def seff(x, w):
    if len(x) < 5:
        return float('nan')
    o = np.argsort(x); x = x[o]; w = np.clip(w[o], 0, None)
    c = np.cumsum(w); tot = c[-1]
    if tot <= 0:
        return float('nan')
    c = c / tot
    best = np.inf
    j = 0
    for i in range(len(x)):
        while j < len(x) and c[j] - (c[i - 1] if i > 0 else 0.0) < 0.683:
            j += 1
        if j >= len(x):
            break
        best = min(best, x[j] - x[i])
    return 0.5 * best


def load(sample, era, m):
    t = uproot.open(f'{BASE}/{sample}/{era}.root')['inclusive']
    a = t.arrays(['H_mass', 'weight', 'z_ee', 'z_mumu', f'MVA_Score_mA_M{m}'], library='np')
    return a


def stats(a, m, ch):
    chm = (a['z_ee'] == 1) if ch == 'ele' else (a['z_mumu'] == 1)
    win = (a['H_mass'] > 95) & (a['H_mass'] < 180)
    mva = a[f'MVA_Score_mA_M{m}'] >= CUTS[m]
    s = chm & win & mva
    s0 = chm & win
    w = a['weight']
    x = a['H_mass'][s]; ww = w[s]
    pk = (x > 100) & (x < 180)
    mean = float(np.sum(x[pk] * ww[pk]) / np.sum(ww[pk])) if np.sum(ww[pk]) > 0 else float('nan')
    return dict(N=float(ww.sum()), Npre=float(w[s0].sum()), mean=mean, seff=seff(x[pk], ww[pk]))


out = {}
masses = [int(v) for v in sys.argv[2].split(',')]
for m in masses:
    for era in ERAS:
        nom = load(f'mA_M{m}', era, m)
        for ch in ['ele', 'mu']:
            sn = stats(nom, m, ch)
            rec = {'nominal': sn}
            for sy in SYSTS:
                su = stats(load(f'mA_M{m}_{sy}_up', era, m), m, ch)
                sd = stats(load(f'mA_M{m}_{sy}_down', era, m), m, ch)
                rec[sy] = dict(
                    rate=0.5 * (abs(su['N'] / sn['N'] - 1) + abs(sd['N'] / sn['N'] - 1)),
                    rate_presel=0.5 * (abs(su['Npre'] / sn['Npre'] - 1) + abs(sd['Npre'] / sn['Npre'] - 1)),
                    dpeak=0.5 * (abs(su['mean'] - sn['mean']) + abs(sd['mean'] - sn['mean'])),
                    dres=0.5 * (abs(su['seff'] - sn['seff']) + abs(sd['seff'] - sn['seff'])),
                    up=su, down=sd)
            out[f'{m}_{era}_{ch}'] = rec
            print(m, era, ch, f"N={sn['N']:.3f} mean={sn['mean']:.2f} seff={sn['seff']:.2f}",
                  ' '.join(f"{sy}:{rec[sy]['rate']*100:.2f}%/{rec[sy]['dpeak']:.3f}/{rec[sy]['dres']:.3f}" for sy in SYSTS),
                  flush=True)
        del nom
json.dump(out, open(sys.argv[1], 'w'), indent=1)
