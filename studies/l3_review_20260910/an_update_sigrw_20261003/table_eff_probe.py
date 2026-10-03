"""Probe which definition reproduces the tab:boundary signal efficiency (Plot/output/latexTable_MVAcut_eff_bkg.txt)."""
import json, numpy as np, uproot, math
B="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"; ERAS=["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
tab={1:0.094,2:0.412,3:0.603,4:0.668,5:0.672,6:0.775,7:0.752,8:0.724,9:0.703,10:0.689,15:0.578,20:0.582,25:0.573,30:0.589}
cuts={e["mA"]:e["MVAcut"] for e in json.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"))["results"]}
for m in tab:
    c=math.floor(cuts[m]/0.005+1e-6)*0.005; res=[]
    for tr,wn in (("test","weight"),("inclusive","weight"),("inclusive","factor")):
        num=den=0
        for e in ERAS:
            a=uproot.open(f"{B}/mA_M{m}/{e}.root")[tr].arrays([wn,"H_mass",f"MVA_Score_mA_M{m}"],library="np")
            sr=(a["H_mass"]>115)&(a["H_mass"]<135); w=a[wn]
            den+=w[sr].sum(); num+=w[sr&(a[f"MVA_Score_mA_M{m}"]>=c)].sum()
        res.append(num/den)
    print(m, cuts[m], c, "table %.3f"%tab[m], " ".join("%.3f"%r for r in res))
