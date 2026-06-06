#!/usr/bin/env python3
"""Quantitative sculpting diagnostic for the HZa low-ma BDT.
Read-only: loads background DY (All_Bkg inputs) + trained model, computes
- feature importance (gain)
- weighted Pearson corr(feature, H_m) for background
- per-mass-hypothesis corr(score, H_m) and the 125-GeV enhancement after the BDT cut
No code in the analysis repo is modified.
"""
import numpy as np, uproot, pickle, json

FEATURES = ["pho1Pt_oHm","pho1R9","pho1IetaIeta55","pho1PIso_noCorr","pho2Pt_oHm","pho2R9",
            "pho2IetaIeta55","pho2PIso_noCorr","ALP_calculatedPhotonIso","var_dR_Za",
            "var_dR_g1g2","var_dR_g1Z","var_PtaOverMh","H_pt_oHm","pho_pt_asym","param"]
BKG = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal/All_Bkg/run3.root"
MODEL = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_BDT_run3.pkl"
CUTS = {r["mA"]: r["MVAcut"] for r in json.load(open(
    "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"))["results"]}

def wpearson(x, y, w):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
    x, y, w = x[m], y[m], w[m]
    mx = np.average(x, weights=w); my = np.average(y, weights=w)
    cov = np.average((x-mx)*(y-my), weights=w)
    sx = np.sqrt(np.average((x-mx)**2, weights=w)); sy = np.sqrt(np.average((y-my)**2, weights=w))
    return cov/(sx*sy) if sx>0 and sy>0 else np.nan

base = [f for f in FEATURES if f not in ("param","pho_pt_asym")]
arr = uproot.open(BKG)["inclusive"].arrays(base+["H_m","ALP_m","factor"], library="np")
Hm = arr["H_m"].astype(float); ALPm = arr["ALP_m"].astype(float); w = arr["factor"].astype(float)
sel = (Hm>95)&(Hm<180)
for k in arr: arr[k] = arr[k][sel]
Hm, ALPm, w = Hm[sel], ALPm[sel], w[sel]
p1, p2 = arr["pho1Pt_oHm"].astype(float), arr["pho2Pt_oHm"].astype(float)
pho_pt_asym = (p1-p2)/(p1+p2+1e-6)
N = len(Hm); print(f"background events (95<H_m<180): {N}, sum_w={w.sum():.1f}\n")

# A. feature importance (gain)
model = pickle.load(open(MODEL,"rb"))
imp = model.feature_importances_
order = np.argsort(imp)[::-1]
print("=== A. BDT feature importance (gain), sorted ===")
for i in order: print(f"  {FEATURES[i]:<26} {imp[i]:.4f}")

# B. corr(feature, H_m) for background
print("\n=== B. weighted |corr(feature, H_m)| for background (sculpting risk) ===")
corrs = {f: wpearson(arr[f].astype(float), Hm, w) for f in base}
corrs["pho_pt_asym"] = wpearson(pho_pt_asym, Hm, w)
for f,_ in sorted(corrs.items(), key=lambda kv: -abs(kv[1])):
    print(f"  {f:<26} corr={corrs[f]:+.3f}")

# C. per-hypothesis score vs H_m + 125-peak enhancement after cut
def build_X(mhyp):
    cols = []
    for f in FEATURES:
        if f=="param": cols.append((ALPm - mhyp)/Hm)
        elif f=="pho_pt_asym": cols.append(pho_pt_asym)
        else: cols.append(arr[f].astype(float))
    return np.column_stack(cols)

print("\n=== C. per mass-hypothesis: corr(score,H_m) & 125-peak enhancement after cut ===")
peak = (Hm>120)&(Hm<130)
frac_all = w[peak].sum()/w.sum()
print(f"  (background fraction in 120-130 before any cut: {frac_all:.3f})")
print(f"  {'mhyp':>4} {'cut':>6} {'corr(score,Hm)':>15} {'pass_frac':>10} {'peakfrac_pass':>13} {'enhance(R)':>11}")
for m in (1,3,20):
    X = build_X(m)
    score = model.predict_proba(X)[:,1]
    c = wpearson(score, Hm, w)
    cut = CUTS.get(m, 0.99)
    pas = score>cut
    pw = w[pas]
    pass_frac = pw.sum()/w.sum()
    peakfrac_pass = w[pas & peak].sum()/pw.sum() if pw.sum()>0 else np.nan
    R = peakfrac_pass/frac_all if frac_all>0 else np.nan
    print(f"  {m:>4} {cut:>6.3f} {c:>+15.3f} {pass_frac:>10.4f} {peakfrac_pass:>13.3f} {R:>11.2f}")
print("\n  R>1 => BDT cut concentrates background into the 125 GeV Higgs window (sculpting).")
