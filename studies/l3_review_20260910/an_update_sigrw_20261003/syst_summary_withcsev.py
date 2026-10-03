#!/usr/bin/env python3
"""Recompute tab:uncertainty_summary (AN-25-172 Sec 9) from HZa lepton datacards.

Signal-yield impact per nuisance group at each mass point, from the lnN lines of
<mA>_pruned_datacard_leptons.txt (10 signal processes = 5 eras x {ele,mu}).
Per-process weights = nominal signal yield from yields_leptons/<mA>_pruned_Datacard_leptons.pkl.

Definitions computed (the AN table never had a committed script; these variants are
compared against the pre-FSR-fix datacards to identify which reproduces the old table):
  wavg : yield-weighted average over all 10 processes of |delta| (symmetric avg of up/down);
         processes without the nuisance contribute 0.
  tot  : shift of the total signal yield when the (correlated) nuisance moves by +-1 sigma,
         |sum_i y_i (k_i - 1)| / sum_i y_i, averaged over up/down.
  max  : max |delta| over processes.
Groups combine member nuisances in quadrature per process before averaging.

usage: python3 syst_summary.py <datacard_dir> [<pkl_dir>]
"""
import sys, os, re, math, glob
import pandas as pd

GROUPS = {
    'QCD scale': ['QCDscale_ggH'],
    'alpha_s': ['alphaS_ggH'],
    'Pileup': ['CMS_hza_pileup'],
    'Luminosity': ['lumi_13p6TeV_Uncorrelated_2022preEE', 'lumi_13p6TeV_Uncorrelated_2022postEE',
                   'lumi_13p6TeV_Uncorrelated_2023preBPix', 'lumi_13p6TeV_Uncorrelated_2023postBPix',
                   'lumi_13p6TeV_Uncorrelated_2024', 'lumi_13p6TeV_Correlated_2223',
                   'lumi_13p6TeV_Correlated_2324', 'lumi_13p6TeV_Correlated_222324'],
    'PDF': ['pdf_Higgs_ggH'],
    'Shape reweighting': ['CMS_hza_mva_reweight'],
    'Electron sel. eff.': ['CMS_hza_electron_reco', 'CMS_hza_electron_id', 'CMS_hza_electron_id_nomatch'],
    'Photon sel. eff.': ['CMS_hza_photon_id'],
    'Muon sel. eff.': ['CMS_hza_muon_reco', 'CMS_hza_muon_id', 'CMS_hza_muon_id_nomatch'],
    'Trigger eff.': ['CMS_hza_trigger'],
    'CSEV': ['CMS_hza_photon_csev_2022preEE','CMS_hza_photon_csev_2022postEE','CMS_hza_photon_csev_2023preBPix','CMS_hza_photon_csev_2023postBPix','CMS_hza_photon_csev_2024'],
}
# also report individual electron pieces
EXTRA = ['CMS_hza_electron_reco', 'CMS_hza_electron_id', 'CMS_hza_electron_id_nomatch',
         'CMS_hza_muon_reco', 'CMS_hza_muon_id', 'CMS_hza_muon_id_nomatch']


def parse_kappa(tok):
    """return (k_up, k_down) as multiplicative factors; None if '-'"""
    if tok == '-':
        return None
    if '/' in tok:
        dn, up = tok.split('/')
        return float(up), float(dn)
    k = float(tok)
    return k, 1.0 / k


def parse_card(path):
    procs = None
    lnN = {}
    for line in open(path):
        w = line.split()
        if not w:
            continue
        if w[0] == 'process' and procs is None and not re.match(r'^-?\d', w[1]):
            procs = w[1:]
        if len(w) > 2 and w[1] == 'lnN':
            lnN[w[0]] = [parse_kappa(t) for t in w[2:]]
    return procs, lnN


def main():
    dcdir = sys.argv[1]
    pkldir = sys.argv[2] if len(sys.argv) > 2 else \
        '/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/yields_leptons'
    masses = list(range(1, 31))
    rows = {}
    for m in masses:
        cands = [os.path.join(dcdir, f'{m}.txt'), os.path.join(dcdir, f'{m}_pruned_datacard_leptons.txt')]
        card = [c for c in cands if os.path.exists(c)][0]
        procs, lnN = parse_card(card)
        df = pd.read_pickle(os.path.join(pkldir, f'{m}_pruned_Datacard_leptons.pkl'))
        ymap = dict(zip(df['proc'], df['nominal_yield']))
        sig = [i for i, p in enumerate(procs) if p.startswith('ggH')]
        y = [float(ymap[procs[i]]) for i in sig]
        Y = sum(y)
        res = {}
        for gname, members in list(GROUPS.items()) + [(e, [e]) for e in EXTRA]:
            dsym = []   # per-process quadrature-combined symmetric |delta|
            dup = []    # per-process signed up shift (quadrature not meaningful for sign; use sum)
            ddn = []
            for i in sig:
                s2 = 0.0; su = 0.0; sd = 0.0
                for nm in members:
                    if nm not in lnN:
                        continue
                    k = lnN[nm][i]
                    if k is None:
                        continue
                    ku, kd = k
                    d = 0.5 * (abs(ku - 1) + abs(kd - 1))
                    s2 += d * d; su += ku - 1; sd += kd - 1
                dsym.append(math.sqrt(s2)); dup.append(su); ddn.append(sd)
            wavg = sum(yi * di for yi, di in zip(y, dsym)) / Y
            if len(members) == 1:
                tot = 0.5 * (abs(sum(yi * u for yi, u in zip(y, dup))) + abs(sum(yi * d for yi, d in zip(y, ddn)))) / Y
            else:
                # correlated within a nuisance, uncorrelated between nuisances: quadrature of each nuisance's total shift
                t2 = 0.0
                for nm in members:
                    if nm not in lnN:
                        continue
                    su = sum(yi * ((lnN[nm][i][0] - 1) if lnN[nm][i] else 0) for yi, i in zip(y, sig))
                    sd = sum(yi * ((lnN[nm][i][1] - 1) if lnN[nm][i] else 0) for yi, i in zip(y, sig))
                    t2 += (0.5 * (abs(su) + abs(sd)) / Y) ** 2
                tot = math.sqrt(t2)
            mx = max(dsym)
            res[gname] = (100 * wavg, 100 * tot, 100 * mx)
        rows[m] = res
    names = list(GROUPS) + EXTRA
    import json
    json.dump({str(m): {n: rows[m][n] for n in names} for m in masses}, open(os.environ.get('SUMMARY_JSON', '/dev/null'), 'w'), indent=1)
    print(f'# datacards: {dcdir}\n# weights: {pkldir}')
    print('# per-mass wavg / tot / max  [%]')
    hdr = 'mA  ' + ' '.join(f'{n[:14]:>20s}' for n in names)
    print(hdr)
    for m in masses:
        print(f'{m:<3d} ' + ' '.join(f'{rows[m][n][0]:.2f}/{rows[m][n][1]:.2f}/{rows[m][n][2]:.2f}'.rjust(20) for n in names))
    print('\n# summary over 30 mass points: min--max (mean), per definition')
    for n in names:
        out = []
        for j, lab in enumerate(['wavg', 'tot', 'max']):
            v = [rows[m][n][j] for m in masses]
            out.append(f'{lab}: {min(v):5.2f}--{max(v):5.2f} (mean {sum(v)/len(v):5.2f})')
        print(f'{n:22s} ' + ' | '.join(out))


if __name__ == '__main__':
    main()
