#!/usr/bin/env python3
"""Signal-yield impact of the photon conversion-safe electron veto (CSEV) SF uncertainty.

CMS_hza_photon_csev is NOT a nuisance in the HZa lepton datacards, so the AN's "0--2%" cannot be
read from them. This recomputes it from the post-MVA-cut signal trees
(/eos/home-p/pelai/HZa/root_MVAcut/sig/mA_M<m>/output_<era>.root, the input of Tree2WS/makeYields),
with the same up/central ratio used for the other factory nuisances. photon_id is computed the same
way as a closure check against the datacard CMS_hza_photon_id values.

delta_X(proc) = sum(w * X_Up / X_central) / sum(w) - 1   (and Down), w = weight
"""
import uproot, numpy as np, json, sys

BASE = '/eos/home-p/pelai/HZa/root_MVAcut/sig'
ERAS = ['2022preEE', '2022postEE', '2023preBPix', '2023postBPix', '2024']
ANCHORS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
SYSTS = ['photon_csev_sf_SelectedPhoton', 'photon_id_sf_SelectedPhoton']

out = {}
for m in ANCHORS:
    for era in ERAS:
        f = uproot.open(f'{BASE}/mA_M{m}/output_{era}.root')
        for ch in ['ele', 'mu']:
            t = f[f'DiphotonTree/ggh_125_Za_{ch}_13p6TeV_cat0']
            br = ['weight'] + [f'weight_{s}_{v}' for s in SYSTS for v in ['central', 'Up', 'Down']]
            a = t.arrays(br, library='np')
            w = a['weight']; W = w.sum()
            rec = {'sumw': float(W)}
            for s in SYSTS:
                c = a[f'weight_{s}_central']
                with np.errstate(divide='ignore', invalid='ignore'):
                    ru = np.where(c != 0, a[f'weight_{s}_Up'] / c, 1.0)
                    rd = np.where(c != 0, a[f'weight_{s}_Down'] / c, 1.0)
                rec[s] = [float((w * ru).sum() / W), float((w * rd).sum() / W)]
            out[f'{m}_{era}_{ch}'] = rec
            print(m, era, ch, f"sumw={W:.4f}", ' '.join(f"{s.split('_sf')[0]}: {rec[s][0]:.4f}/{rec[s][1]:.4f}" for s in SYSTS), flush=True)

json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else 'csev_from_trees.json', 'w'), indent=1)

print('\n# yield-weighted |delta| over 10 processes per anchor (same definition as syst_summary.py wavg)')
for s in SYSTS:
    vals = []
    for m in ANCHORS:
        num = 0.0; den = 0.0; mx = 0.0
        for era in ERAS:
            for ch in ['ele', 'mu']:
                r = out[f'{m}_{era}_{ch}']
                d = 0.5 * (abs(r[s][0] - 1) + abs(r[s][1] - 1))
                num += r['sumw'] * d; den += r['sumw']; mx = max(mx, d)
        vals.append(100 * num / den)
        print(s, m, f'wavg={100*num/den:.3f}%  max_proc={100*mx:.3f}%')
    print(s, f'range {min(vals):.2f}--{max(vals):.2f}  mean(anchors) {np.mean(vals):.2f}')
