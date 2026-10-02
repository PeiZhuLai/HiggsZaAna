#!/usr/bin/env python3
"""Aggregate shape_json/*.json (from shape_syst_from_scored.py) into AN-level numbers.
yield: yield-weighted average over the 10 processes (5 eras x ele/mu) of the per-process rate, per mA;
       reported as range over the 14 anchors.  peak/resolution: per-process GeV shifts, median and max,
       separately for mA=1 and mA>=2.  Groups: FNUF, Material, lepton/photon scale+smear (max of members
       per process for GeV numbers; quadrature for yield)."""
import json, glob, os, statistics as st, math
W = os.path.dirname(os.path.abspath(__file__))
d = {}
for f in glob.glob(f'{W}/shape_json/shape_g*.json') + [f'{W}/shape_json/shape_5.json']:
    d.update(json.load(open(f)))
ERAS = ['2022preEE', '2022postEE', '2023preBPix', '2023postBPix', '2024']
ANCH = sorted({int(k.split('_')[0]) for k in d})
print('anchors', ANCH, 'entries', len(d))
GROUPS = {'FNUF': ['FNUF'], 'Material': ['Material'],
          'Electron scale/smear': ['Electron_scale', 'Electron_smear'],
          'Muon scale/smear': ['Muon_scale', 'Muon_smear'],
          'Photon scale/smear': ['Photon_scale', 'Photon_smear'],
          'All e/mu/gamma scale+smear': ['Electron_scale', 'Electron_smear', 'Muon_scale', 'Muon_smear', 'Photon_scale', 'Photon_smear']}
for g, mem in GROUPS.items():
    yv = []; ymax_proc = 0; ypre = []
    pk1 = []; pk2 = []; rs1 = []; rs2 = []
    for m in ANCH:
        num = den = numpre = 0.0
        for e in ERAS:
            for ch in ['ele', 'mu']:
                r = d[f'{m}_{e}_{ch}']; N = r['nominal']['N']
                rate = math.sqrt(sum(r[s]['rate'] ** 2 for s in mem))
                ratepre = math.sqrt(sum(r[s]['rate_presel'] ** 2 for s in mem))
                num += N * rate; den += N; numpre += N * ratepre; ymax_proc = max(ymax_proc, rate)
                pk = max(r[s]['dpeak'] for s in mem); rs = max(r[s]['dres'] for s in mem)
                (pk1 if m == 1 else pk2).append(pk); (rs1 if m == 1 else rs2).append(rs)
        yv.append(100 * num / den); ypre.append(100 * numpre / den)
    print(f"\n[{g}]")
    print(f"  yield (after MVA cut, yield-weighted over procs): {min(yv):.3f}--{max(yv):.3f}%  mean {st.mean(yv):.3f}%   max single process {100*ymax_proc:.2f}%")
    print(f"  yield (preselection only)                        : {min(ypre):.3f}--{max(ypre):.3f}%  mean {st.mean(ypre):.3f}%")
    print(f"  peak shift  mA>=2: median {st.median(pk2):.3f}  max {max(pk2):.3f} GeV | mA=1: median {st.median(pk1):.3f} max {max(pk1):.3f} GeV")
    print(f"  resolution  mA>=2: median {st.median(rs2):.3f}  max {max(rs2):.3f} GeV | mA=1: median {st.median(rs1):.3f} max {max(rs1):.3f} GeV")
