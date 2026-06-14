#!/usr/bin/env python3
"""Verify the two-model add_mva_scores: route ma1-3 -> low model, ma>=4 -> high model, _oHm
features. Cross-check against directly applying the deployed models. Run in higgs-alp-ana."""
import numpy as np, pandas as pd, uproot, pickle
import Parque2Root_BDT as P

ROOT = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal/mA_M5/run3.root"
BASE = ["pho1Pt_oHm","pho1R9","pho1IetaIeta55","pho1PIso_noCorr","pho2Pt_oHm","pho2R9",
        "pho2IetaIeta55","pho2PIso_noCorr","ALP_calculatedPhotonIso","var_dR_Za","var_dR_g1g2",
        "var_dR_g1Z","var_PtaOverMh","H_pt_oHm","H_m","ALP_m"]
a = uproot.open(ROOT)["inclusive"].arrays(BASE, library="np")
df = pd.DataFrame({k: a[k] for k in a}); df["H_mass"] = df["H_m"]; df["ALP_mass"] = df["ALP_m"]

out = P.add_mva_scores(df.copy(), masses=[1, 5])           # 1 -> low model, 5 -> high model
asym = (a["pho1Pt_oHm"]-a["pho2Pt_oHm"])/(a["pho1Pt_oHm"]+a["pho2Pt_oHm"]+1e-6)

# direct low model @ ma=1
low = pickle.load(open(P.LOW_MODEL_FILE, "rb"))
param1 = (a["ALP_m"]-1)/a["H_m"]
Xl = np.column_stack([a["ALP_calculatedPhotonIso"], a["var_dR_g1g2"], a["pho1R9"], a["pho1Pt_oHm"],
                      a["pho1PIso_noCorr"], param1])
s_low = low.predict_proba(Xl)[:, 1]
# direct high model @ ma=5
high = pickle.load(open(P.HIGH_MODEL_FILE, "rb"))
param5 = (a["ALP_m"]-5)/a["H_m"]
Xh = np.column_stack([a[c] for c in ["pho1Pt_oHm","pho1R9","pho1IetaIeta55","pho1PIso_noCorr",
        "pho2Pt_oHm","pho2R9","pho2IetaIeta55","pho2PIso_noCorr","ALP_calculatedPhotonIso",
        "var_dR_Za","var_dR_g1g2","var_dR_g1Z","var_PtaOverMh","H_pt_oHm"]] + [asym, param5])
s_high = high.predict_proba(Xh)[:, 1]

m = np.isfinite(out["MVA_Score_mA_M1"].to_numpy())
d1 = np.nanmax(np.abs(out["MVA_Score_mA_M1"].to_numpy()[m] - s_low[m]))
d5 = np.nanmax(np.abs(out["MVA_Score_mA_M5"].to_numpy()[m] - s_high[m]))
print(f"rows={m.sum()}")
print(f"MVA_Score_mA_M1 (low)  mean={np.nanmean(out['MVA_Score_mA_M1']):.3f}  max|diff vs direct|={d1:.2e}")
print(f"MVA_Score_mA_M5 (high) mean={np.nanmean(out['MVA_Score_mA_M5']):.3f}  max|diff vs direct|={d5:.2e}")
print("VERDICT:", "OK (routing+features consistent)" if max(d1, d5) < 1e-5 else "MISMATCH!")
