#!/usr/bin/env python3
"""Torch-free shared data layer for the HZa classifiers (NN and BDT).

Holds the feature definitions + ROOT loading so that BOTH environments can import the SAME
code without pulling in torch:
  - hza_NN (torch)        : train_nn.py / apply_nn.py
  - higgs-alp-ana (xgboost): compare_nn_vs_bdt.py (BDT scoring)
The 16-feature order here is identical to run3_Za_BDT.py `variables + ['param']`, so one
_feat() matrix scores either model (NN standardizes by mu/sd; BDT feeds it raw).
"""
from __future__ import annotations
import numpy as np, uproot

ROOT_DIR = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal"
MASSES   = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
# BDT input variables (mirror run3_Za_BDT.py `variables` + param; photon/Higgs pT are H_m-normalized)
BASE_VARS = ["pho1Pt_oHm","pho1R9","pho1IetaIeta55","pho1PIso_noCorr","pho2Pt_oHm","pho2R9",
             "pho2IetaIeta55","pho2PIso_noCorr","ALP_calculatedPhotonIso","var_dR_Za",
             "var_dR_g1g2","var_dR_g1Z","var_PtaOverMh","H_pt_oHm"]   # + pho_pt_asym + param
FEATURES  = BASE_VARS + ["pho_pt_asym", "param"]
NFEAT     = len(FEATURES)
READ      = BASE_VARS + ["H_m","ALP_m","factor"]
# Extra branches the sideband reweighter needs on the background frame
# (Z_m is a reweight step variable; event drives the per-event mass hypothesis).
REWEIGHT_EXTRA = ["Z_m", "event"]

def _load(path, tree, nmax=None, seed=0, extra=()):
    a = uproot.open(path)[tree].arrays(list(READ)+list(extra), library="np")
    m = (a["H_m"]>95)&(a["H_m"]<180)&np.isfinite(a["factor"])
    a = {k:a[k][m].astype(np.float64) for k in a}
    if nmax and len(a["H_m"])>nmax:
        rng = np.random.default_rng(seed); idx = rng.choice(len(a["H_m"]), nmax, replace=False)
        a = {k:v[idx] for k,v in a.items()}
    return a

def _feat(a, mhyp):
    asym  = (a["pho1Pt_oHm"]-a["pho2Pt_oHm"])/(a["pho1Pt_oHm"]+a["pho2Pt_oHm"]+1e-6)
    param = (a["ALP_m"]-mhyp)/a["H_m"]
    return np.column_stack([a[b] for b in BASE_VARS]+[asym, param])

def _sideband_reweight(bkg):
    """Per-event sideband reweight factor for the background frame (1.0 if JSON/helper missing).
    Mirrors run3_Za_BDT.apply_sideband_reweight_to_bkg; needs Z_m/event + the _oHm pT branches."""
    n = len(bkg["H_m"])
    try:
        import pandas as pd
        from sideband_reweight import load_sideband_reweighter
    except Exception as e:
        print("[reweight] helper import failed, using weight=1:", e); return np.ones(n)
    rw = load_sideband_reweighter()
    if rw is None:
        print("[reweight] no sideband JSON found, using weight=1"); return np.ones(n)
    frame = pd.DataFrame({k: bkg[k] for k in bkg})     # has _oHm pT, Z_m, H_m, ALP_m, event
    rwgt = rw.weights_for_dataframe(frame)
    print(f"[reweight] applied sideband reweight: mean={np.nanmean(rwgt):.3f}")
    return np.where(np.isfinite(rwgt), rwgt, 1.0)

def build_training_set(n_sig_per=8000, n_bkg=150000, tree="train"):
    """Build a balanced (sig vs bkg) feature/label/weight set from a given tree.
    tree="train" for training; tree="validation"/"test" gives an independent set for the
    overtraining check (the p2root files carry inclusive/train/validation/test splits)."""
    sig = [_load(f"{ROOT_DIR}/mA_M{m}/run3.root",tree,n_sig_per) for m in MASSES]
    Xs  = np.vstack([_feat(s,m) for s,m in zip(sig,MASSES)])
    ws  = np.concatenate([s["factor"] for s in sig])
    bkg = _load(f"{ROOT_DIR}/All_Bkg/run3.root",tree,n_bkg, extra=REWEIGHT_EXTRA)
    bhyp= np.random.default_rng(1).choice(MASSES, size=len(bkg["H_m"]))
    Xb  = _feat(bkg,0); Xb[:,-1] = (bkg["ALP_m"]-bhyp)/bkg["H_m"]
    wb  = bkg["factor"] * _sideband_reweight(bkg)   # apply sideband reweight to the background weight
    X   = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    wsig= ws*(np.clip(wb,0,None).sum()/np.clip(ws,0,None).sum())          # balance to equal yield
    w   = np.clip(np.concatenate([wsig, wb]), 0, None)                     # NN needs w>=0
    mass= np.concatenate([np.zeros(len(Xs)), bkg["H_m"]])                  # bkg m_llgg for DisCo
    isb = np.concatenate([np.zeros(len(Xs)), np.ones(len(Xb))]).astype(bool)
    # ALP_m / H_m per event -> needed to recompute param at a FIXED mass hypothesis for per-mass DisCo
    alp = np.concatenate([np.concatenate([s["ALP_m"] for s in sig]), bkg["ALP_m"]])
    hm  = np.concatenate([np.concatenate([s["H_m"]   for s in sig]), bkg["H_m"]])
    return X, y, w, mass, isb, alp, hm
