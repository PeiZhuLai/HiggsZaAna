#!/usr/bin/env python3
"""Controlled A/B: retrain a light parametrized XGBoost twice with identical
settings, differing ONLY in whether pho1Pt/pho2Pt/H_pt are normalized by H_m.
Compare 125-GeV sculpting R at a FIXED background efficiency (fair comparison).
Read-only w.r.t. the analysis repo (self-contained training in /tmp)."""
import numpy as np, uproot, json
from xgboost import XGBClassifier

ROOT = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal"
MASSES = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
BASE = ["pho1Pt","pho1R9","pho1IetaIeta55","pho1PIso_noCorr","pho2Pt","pho2R9",
        "pho2IetaIeta55","pho2PIso_noCorr","ALP_calculatedPhotonIso","var_dR_Za",
        "var_dR_g1g2","var_dR_g1Z","var_PtaOverMh","H_pt","pho_pt_asym"]  # + param appended
NORM_BY_HM = ["pho1Pt","pho2Pt","H_pt"]   # variables to divide by H_m in the "normalized" variant
READ = [v for v in BASE if v != "pho_pt_asym"] + ["H_m","ALP_m","factor"]
rng = np.random.RandomState(7)
HP = dict(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
          colsample_bytree=0.8, min_child_weight=20, eval_metric="logloss", n_jobs=8)

def load(path, tree):
    a = uproot.open(path)[tree].arrays(READ, library="np")
    m = (a["H_m"]>95)&(a["H_m"]<180)
    return {k:a[k][m].astype(float) for k in a}

def feat(a, mhyp, normalize):
    Hm=a["H_m"]; cols=[]
    for v in BASE:
        if v=="pho_pt_asym":
            x=(a["pho1Pt"]-a["pho2Pt"])/(a["pho1Pt"]+a["pho2Pt"]+1e-6)
        else:
            x=a[v]/Hm if (normalize and v in NORM_BY_HM) else a[v]
        cols.append(x)
    cols.append((a["ALP_m"]-mhyp)/Hm)   # param (denominator H_m in both variants)
    return np.column_stack(cols)

# ---- build training set (signal train trees + background train tree) ----
sig=[load(f"{ROOT}/mA_M{m}/run3.root","train") for m in MASSES]
sig_hyp=[np.full(len(s["H_m"]),m,float) for s,m in zip(sig,MASSES)]
bkg=load(f"{ROOT}/All_Bkg/run3.root","train")
bkg_hyp=rng.choice(MASSES,size=len(bkg["H_m"])).astype(float)

def stack(normalize):
    Xs=[feat(s,h[0],normalize) for s,h in zip(sig,sig_hyp)]
    Xsig=np.vstack(Xs); Xbkg=feat(bkg,0,normalize); Xbkg[:,-1]=(bkg["ALP_m"]-bkg_hyp)/bkg["H_m"]
    X=np.vstack([Xsig,Xbkg])
    y=np.concatenate([np.ones(len(Xsig)),np.zeros(len(Xbkg))])
    wsig=np.concatenate([s["factor"] for s in sig]); wbkg=bkg["factor"]
    wsig=wsig*(wbkg.sum()/wsig.sum())   # balance to equal yield
    w=np.concatenate([wsig,wbkg])
    pos=w>0                              # drop negative-weight (NLO) events for training
    return X[pos],y[pos],w[pos]

# ---- background inclusive for sculpting eval ----
bkinc=load(f"{ROOT}/All_Bkg/run3.root","inclusive")
Hm=bkinc["H_m"]; wB=np.clip(bkinc["factor"],0,None); peak=(Hm>120)&(Hm<130); frac_all=wB[peak].sum()/wB.sum()

def R_at_eff(score, eff):
    thr=weighted_quantile(score, wB, 1-eff)
    pas=score>thr; pw=wB[pas].sum()
    return (wB[pas&peak].sum()/pw)/frac_all if pw>0 else np.nan

def weighted_quantile(v,w,q):
    o=np.argsort(v); v,w=v[o],w[o]; c=np.cumsum(w)/w.sum()
    return np.interp(q,c,v)

print(f"signal train evts={sum(len(s['H_m']) for s in sig)}  bkg train={len(bkg['H_m'])}  bkg inclusive={len(Hm)}")
print(f"120-130 baseline fraction={frac_all:.3f}\n")
results={}
for tag,normalize in (("baseline",False),("pT/H_m norm",True)):
    X,y,w=stack(normalize)
    model=XGBClassifier(**HP); model.fit(X,y,sample_weight=w)
    results[tag]={}
    for mh in (1,3,20):
        sc=model.predict_proba(feat(bkinc,mh,normalize))[:,1]
        results[tag][mh]={eff:R_at_eff(sc,eff) for eff in (0.002,0.005)}

print(f"{'mhyp':>4} {'bkg_eff':>8} | {'baseline R':>11} {'pT/H_m R':>10} {'Δ':>7}")
print("-"*48)
for mh in (1,3,20):
    for eff in (0.002,0.005):
        b=results['baseline'][mh][eff]; n=results['pT/H_m norm'][mh][eff]
        print(f"{mh:>4} {eff:>8.3f} | {b:>11.2f} {n:>10.2f} {n-b:>+7.2f}")
    print("-"*48)
print("R>1 = sculpted toward 125. Compare baseline vs pT/H_m at SAME bkg efficiency.")
