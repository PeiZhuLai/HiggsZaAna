# Raw (unsmoothed) MC bkg yield, same inputs/selection as table_interpolate_bkgYield_1.py:
# inclusive tree, weight branch 'weight', 115<H_m<135, MVA_Score_mA_M<m> > cut.
import uproot, numpy as np, os
P="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
S={"2022preEE":["DYGto2LG_10to100","DYJetsToLL"],"2022postEE":["DYGto2LG_10to100","DYJetsToLL"],
   "2023preBPix":["DYGto2LG_10to100","DYJetsToLL"],"2023postBPix":["DYGto2LG_10to100","DYJetsToLL"],
   "2024":["DYGto2LG_10to100","DYJetsTo2E","DYJetsTo2Mu","DYJetsTo2Tau"]}
cache={}
def arr(s,y,m):
    k=(s,y)
    if k not in cache:
        t=uproot.open(f"{P}/{s}/{y}.root")["inclusive"]
        cache[k]=t.arrays(["H_m","weight"]+[f"MVA_Score_mA_M{i}" for i in (10,11,20,21,25,29,30)],library="np")
    return cache[k]
for m,cut in ((10,0.99),(11,0.99),(20,0.985),(21,0.985),(25,0.98),(29,0.975),(30,0.975)):
    tot=0; n=0
    for y,ss in S.items():
        for s in ss:
            a=arr(s,y,m); sel=(a["H_m"]>115)&(a["H_m"]<135)&(a[f"MVA_Score_mA_M{m}"]>cut)
            tot+=a["weight"][sel].sum(); n+=sel.sum()
    print(f"mA{m} cut {cut}: raw weighted yield {tot:.2f}  (unweighted {n})")
