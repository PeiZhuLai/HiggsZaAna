import pickle, uproot, numpy as np
m = pickle.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_lowmass_run3.pkl","rb"))
feats = ["ALP_calculatedPhotonIso","var_dR_g1g2","pho1R9","pho1Pt_oHm"]
for era in ["2024","2022preEE"]:
    t = uproot.open("/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M1/%s.root"%era)["inclusive"]
    a = t.arrays(feats+["ALP_mass","H_mass","MVA_Score_mA_M1"], library="np", entry_stop=3000)
    X = np.column_stack([a[f] for f in feats]+[(a["ALP_mass"]-1)/a["H_mass"]])
    s = m.predict_proba(X)[:,1]
    ok = np.isfinite(a["MVA_Score_mA_M1"])
    print(era, ok.sum(), "max|diff|", np.max(np.abs(s[ok]-a["MVA_Score_mA_M1"][ok])))
for bak in ["bak_preDYveto_20260929"]:
    mo = pickle.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_lowmass_run3.pkl."+bak,"rb"))
    so = mo.predict_proba(X)[:,1]; print("old model max|diff|", np.max(np.abs(so[ok]-a["MVA_Score_mA_M1"][ok])))
