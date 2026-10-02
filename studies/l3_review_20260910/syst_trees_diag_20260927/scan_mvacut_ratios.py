#!/usr/bin/env python3
"""READ-ONLY. For every mA x era, syst/nominal sumw ratio of the MVAcut signal trees (= Trees2WS input)."""
import uproot, numpy as np, os, collections
B="/eos/home-p/pelai/HZa/root_MVAcut/sig"
MAS=[1,2,3,4,5,6,7,8,9,10,15,20,25,30]; ERAS=["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
S=["FNUF","Material","ElectronScale","ElectronSmear","MuonScale","MuonSmear","PhotonScale","PhotonSmear"]
rows=[]
for m in MAS:
  for e in ERAS:
    f=uproot.open(f"{B}/mA_M{m}/output_{e}.root")
    for lep in ("ele","mu"):
      t0=f"DiphotonTree/ggh_125_Za_{lep}_13p6TeV_cat0"; w0=f[t0]["weight"].array(library="np").sum()
      for s in S:
        for d in ("Up","Down"):
          w=f[f"{t0}_{s}{d}01sigma"]["weight"].array(library="np").sum()
          rows.append((m,e,lep,s,d,w/w0))
with open("mvacut_syst_ratio_all.tsv","w") as o:
  o.write("mA\tera\tlep\tsyst\tdir\tratio\n")
  for r in rows: o.write("%d\t%s\t%s\t%s\t%s\t%.5f\n"%r)
r=np.array([x[5] for x in rows])
print("N=%d  |ratio-1|<1%%: %d  <2%%: %d  ratio<0.9: %d  ratio<0.5: %d  median %.3f  min %.3f  max %.3f"%(len(r),(abs(r-1)<0.01).sum(),(abs(r-1)<0.02).sum(),(r<0.9).sum(),(r<0.5).sum(),np.median(r),r.min(),r.max()))
by=collections.defaultdict(list)
for x in rows: by[x[3]].append(x[5])
for s in S:
  a=np.array(by[s]); print("%-14s median %.3f  frac(|r-1|<1%%) %.2f  min %.3f max %.3f"%(s,np.median(a),(abs(a-1)<0.01).mean(),a.min(),a.max()))
