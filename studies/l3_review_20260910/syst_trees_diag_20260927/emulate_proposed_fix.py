#!/usr/bin/env python3
"""READ-ONLY. Emulate the proposed apply_bdt_sig.py fix for syst trees:
syst INCLUSIVE tree -> keep (run,lumi,event) in NOMINAL test-tree keys -> syst's own MVA cut (no join to nominal pass set).
Compare with nominal MVAcut yield (unchanged by the fix). Also the rate calcPhotonSyst would compute."""
import sys, json, os, uproot, numpy as np
SC="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
cuts={int(r["mA"]):float(r["MVAcut"]) for r in json.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"))["results"]}
ID=["run","luminosityBlock","event"]
def k(d): return list(zip(*[d[c].to_numpy().astype(np.int64) for c in ID]))
for cell in sys.argv[1:]:
  ma,era=cell.split(":"); mva="MVA_Score_"+ma; cut=cuts[int(ma.split("_M")[1])]; cols=ID+[mva,"weight","n_electrons","n_muons"]
  nt=uproot.open(f"{SC}/{ma}/{era}.root")["test"].arrays(cols,library="pd"); ntk=set(k(nt))
  n=nt[nt[mva]>cut]; y0={"ele":n[n.n_electrons==2].weight.sum(),"mu":n[n.n_muons==2].weight.sum()}
  print("="*80); print(f"{ma} {era} nominal ele {y0['ele']:.4f} mu {y0['mu']:.4f}")
  for s in ["FNUF","Material","Electron_scale","Electron_smear","Muon_scale","Muon_smear","Photon_scale","Photon_smear"]:
    y={}
    for ud in ("up","down"):
      d=uproot.open(f"{SC}/{ma}_{s}_{ud}/{era}.root")["inclusive"].arrays(cols,library="pd")
      d=d[np.array([x in ntk for x in k(d)])]; d=d[d[mva]>cut]
      y[ud]={"ele":d[d.n_electrons==2].weight.sum(),"mu":d[d.n_muons==2].weight.sum()}
    line=f"{s:15s}"
    for lep in ("ele","mu"):
      u,dn=y["up"][lep],y["down"][lep]; rate=abs(u-dn)/(u+dn) if u+dn else 0
      line+=f" | {lep}: up {100*(u/y0[lep]-1):+6.2f}% down {100*(dn/y0[lep]-1):+6.2f}% rate {100*rate:5.2f}%"
    print(line)
