#!/usr/bin/env python3
"""What-if test: keep the trained model fixed, but at *apply* time compute the
parametrized feature with a FIXED denominator (m_H nominal = 125.38) instead of
the per-event H_m. Compare the 125-GeV sculpting enhancement R for both.
Read-only; analysis repo untouched."""
import numpy as np, uproot, pickle, json

FEATURES = ["pho1Pt","pho1R9","pho1IetaIeta55","pho1PIso_noCorr","pho2Pt","pho2R9",
            "pho2IetaIeta55","pho2PIso_noCorr","ALP_calculatedPhotonIso","var_dR_Za",
            "var_dR_g1g2","var_dR_g1Z","var_PtaOverMh","H_pt","pho_pt_asym","param"]
BKG = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal/All_Bkg/run3.root"
MODEL = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_BDT_run3.pkl"
CUTS = {r["mA"]: r["MVAcut"] for r in json.load(open(
    "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"))["results"]}
MH_NOM = 125.38

base = [f for f in FEATURES if f not in ("param","pho_pt_asym")]
arr = uproot.open(BKG)["inclusive"].arrays(base+["H_m","ALP_m","factor"], library="np")
Hm = arr["H_m"].astype(float); ALPm = arr["ALP_m"].astype(float); w = arr["factor"].astype(float)
sel = (Hm>95)&(Hm<180)
for k in arr: arr[k] = arr[k][sel]
Hm, ALPm, w = Hm[sel], ALPm[sel], w[sel]
p1, p2 = arr["pho1Pt"].astype(float), arr["pho2Pt"].astype(float)
pho_pt_asym = (p1-p2)/(p1+p2+1e-6)
model = pickle.load(open(MODEL,"rb"))

def build_X(mhyp, denom):
    cols = []
    for f in FEATURES:
        if f=="param": cols.append((ALPm - mhyp)/denom)
        elif f=="pho_pt_asym": cols.append(pho_pt_asym)
        else: cols.append(arr[f].astype(float))
    return np.column_stack(cols)

def wcorr(x,y,wt):
    m=np.isfinite(x)&np.isfinite(y); x,y,wt=x[m],y[m],wt[m]
    mx=np.average(x,weights=wt); my=np.average(y,weights=wt)
    cov=np.average((x-mx)*(y-my),weights=wt)
    sx=np.sqrt(np.average((x-mx)**2,weights=wt)); sy=np.sqrt(np.average((y-my)**2,weights=wt))
    return cov/(sx*sy) if sx>0 and sy>0 else np.nan

peak=(Hm>120)&(Hm<130); frac_all=w[peak].sum()/w.sum()
print(f"background events: {len(Hm)}  | 120-130 baseline fraction: {frac_all:.3f}\n")
print(f"{'mhyp':>4} {'cut':>6} | {'denom':>8} {'corr(score,Hm)':>15} {'pass_frac':>10} {'R(125 enh)':>11}")
print("-"*70)
for m in (1,3,20):
    cut=CUTS.get(m,0.99)
    for tag,denom in (("H_m",Hm),("125.38",MH_NOM)):
        sc=model.predict_proba(build_X(m,denom))[:,1]
        pas=sc>cut; pw=w[pas]
        R=(w[pas&peak].sum()/pw.sum())/frac_all if pw.sum()>0 else np.nan
        c=wcorr(sc,Hm,w)
        print(f"{m:>4} {cut:>6.3f} | {tag:>8} {c:>+15.3f} {pw.sum()/w.sum():>10.4f} {R:>11.2f}")
    print("-"*70)
print("R>1 = background sculpted toward 125 GeV. Compare H_m vs 125.38 rows per mhyp.")
