#!/usr/bin/env python3
"""DisCo proxy step 1/2 (higgs-alp-ana, has uproot): read _oHm root, dump npz for the
torch step. No analysis-repo writes."""
import numpy as np, uproot
np.random.seed(0)
ROOT = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal"
MASSES = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
BASE = ["pho1Pt_oHm","pho1R9","pho1IetaIeta55","pho1PIso_noCorr","pho2Pt_oHm","pho2R9",
        "pho2IetaIeta55","pho2PIso_noCorr","ALP_calculatedPhotonIso","var_dR_Za",
        "var_dR_g1g2","var_dR_g1Z","var_PtaOverMh","H_pt_oHm"]
READ = BASE + ["H_m","ALP_m","factor"]
N_SIG_PER, N_BKG = 8000, 150000

def load(path, tree, nmax=None):
    a = uproot.open(path)[tree].arrays(READ, library="np")
    m = (a["H_m"]>95)&(a["H_m"]<180)&np.isfinite(a["factor"])
    a = {k:a[k][m].astype(np.float64) for k in a}
    if nmax and len(a["H_m"])>nmax:
        idx = np.random.choice(len(a["H_m"]), nmax, replace=False); a={k:v[idx] for k,v in a.items()}
    return a

def feat(a, mhyp):
    asym=(a["pho1Pt_oHm"]-a["pho2Pt_oHm"])/(a["pho1Pt_oHm"]+a["pho2Pt_oHm"]+1e-6)
    param=(a["ALP_m"]-mhyp)/a["H_m"]
    return np.column_stack([a[b] for b in BASE]+[asym, param])

sig=[load(f"{ROOT}/mA_M{m}/run3.root","train",N_SIG_PER) for m in MASSES]
Xs=np.vstack([feat(s,m) for s,m in zip(sig,MASSES)]); ws=np.concatenate([s["factor"] for s in sig])
bkg=load(f"{ROOT}/All_Bkg/run3.root","train",N_BKG)
bhyp=np.random.choice(MASSES,size=len(bkg["H_m"]))
Xb=feat(bkg,0); Xb[:,-1]=(bkg["ALP_m"]-bhyp)/bkg["H_m"]
X=np.vstack([Xs,Xb]); y=np.concatenate([np.ones(len(Xs)),np.zeros(len(Xb))])
w=np.concatenate([ws*(np.clip(bkg["factor"],0,None).sum()/np.clip(ws,0,None).sum()), bkg["factor"]])
w=np.clip(w,0,None)
mass=np.concatenate([np.zeros(len(Xs)), bkg["H_m"]])
isb=np.concatenate([np.zeros(len(Xs)), np.ones(len(Xb))]).astype(bool)

bk=load(f"{ROOT}/All_Bkg/run3.root","inclusive")
bk_base=np.column_stack([bk[b] for b in BASE])
np.savez("/tmp/disco_data.npz", X=X, y=y, w=w, mass=mass, isb=isb,
         bk_base=bk_base, bk_Hm=bk["H_m"], bk_ALPm=bk["ALP_m"],
         bk_w=np.clip(bk["factor"],0,None), BASE=np.array(BASE))
print("saved /tmp/disco_data.npz  X=",X.shape," bkinc=",bk_base.shape)
