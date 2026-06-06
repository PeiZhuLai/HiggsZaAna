#!/usr/bin/env python3
"""Skeleton: parametrized MLP classifier for HZa (XGBoost -> NN migration).

loss = BCE + lambda_disco * dCorr^2(score_bkg, m_llgammagamma_bkg)
  - lambda_disco=0 -> plain NN (DisCo off)
  - DisCo penalty decorrelates the score from the Higgs mass to suppress low-ma sculpting.

Runs in the hza_NN env (torch + uproot). Uses the H_m-normalized branches (_oHm).
GPU is used automatically if available (submit via condor_nn/train_nn_gpu.sub on lxplus901).

This is a SKELETON: the training loop / model save / AUC+R eval are functional, but
(1) sideband reweight integration and (2) apply-to-ROOT (writing the score branch for the
downstream MVA-score p2root) are left as clearly marked TODOs.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, uproot, torch, torch.nn as nn
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")                       # headless (GPU/condor node, no display)
import matplotlib.pyplot as plt
from scipy import stats
# Torch-free shared data layer (feature defs + ROOT loading), reused by the BDT comparison
# in the higgs-alp-ana env so both classifiers build features identically.
from hza_features import (ROOT_DIR, MASSES, BASE_VARS, FEATURES, NFEAT, READ, REWEIGHT_EXTRA,
                          _load, _feat, _sideband_reweight, build_training_set)

# --------------------------------------------------------------------------- DisCo loss
def distance_corr(a_var, b_var, nw):
    """Weighted distance correlation (Kasieczka-Shih). a_var,b_var,nw: 1-D tensors; nw normalized to mean 1."""
    a = (a_var.view(-1,1)-a_var.view(1,-1)).abs(); b = (b_var.view(-1,1)-b_var.view(1,-1)).abs()
    am=(a*nw).mean(1); A=a-am.view(-1,1)-am.view(1,-1)+(am*nw).mean()
    bm=(b*nw).mean(1); B=b-bm.view(-1,1)-bm.view(1,-1)+(bm*nw).mean()
    return ((A*B*nw).mean(1)*nw).mean()/torch.sqrt(((A*A*nw).mean(1)*nw).mean()*((B*B*nw).mean(1)*nw).mean()+1e-12)

# ------------------------------------------------------------------------------- model
class ParamMLP(nn.Module):
    def __init__(self, n_in, hidden=(128,128,64), p_drop=0.1):
        super().__init__()
        layers=[]; d=n_in
        for h in hidden:
            layers += [nn.Linear(d,h), nn.ReLU(), nn.Dropout(p_drop)]; d=h
        layers += [nn.Linear(d,1)]
        self.net=nn.Sequential(*layers)
    def forward(self,x): return self.net(x).squeeze(-1)

# ----------------------------------------------------------------------------- training
def train(X, y, w, mass, isb, alp, hm, *, lambda_disco=0.0, disco_masses=(1, 2, 3),
          epochs=30, batch=4000, lr=1e-3, device="cpu", seed=1, val=None):
    """Per-mass DisCo: at each step, for the background, recompute the score at FIXED low-mass
    hypotheses (disco_masses) and penalize dCorr(score@m_hyp, m_llgammagamma). This matches what
    R(mA=1) measures (background scored under the mA=1 hypothesis), unlike the old random-hypo DisCo.

    val=(Xv,yv,wv): optional independent set (e.g. the validation tree); if given, per-epoch
    validation BCE is recorded so the loss curve shows train-vs-val divergence (overtraining).
    Returns (net, mu, sd, history) where history holds per-epoch bce/disco/total/val_bce."""
    torch.manual_seed(seed)
    mu, sd = X.mean(0), X.std(0)+1e-9; Xn=(X-mu)/sd
    Xt=torch.tensor(Xn,dtype=torch.float32,device=device); yt=torch.tensor(y,dtype=torch.float32,device=device)
    wt=torch.tensor(w,dtype=torch.float32,device=device)
    mt=torch.tensor((mass-125)/30,dtype=torch.float32,device=device); ibt=torch.tensor(isb,device=device)
    alpt=torch.tensor(alp,dtype=torch.float32,device=device); hmt=torch.tensor(hm,dtype=torch.float32,device=device)
    if val is not None:
        Xv,yv,wv = val
        Xvt=torch.tensor((Xv-mu)/sd,dtype=torch.float32,device=device)
        yvt=torch.tensor(yv,dtype=torch.float32,device=device); wvt=torch.tensor(wv,dtype=torch.float32,device=device)
    pidx = NFEAT-1                                   # param is the last feature column
    mu_p = float(mu[pidx]); sd_p = float(sd[pidx])
    net=ParamMLP(NFEAT).to(device); opt=torch.optim.Adam(net.parameters(),lr)
    bce=nn.BCEWithLogitsLoss(reduction='none'); n=len(Xt)
    history={"epoch":[], "bce":[], "disco":[], "total":[], "val_bce":[]}
    for ep in range(epochs):
        net.train(); perm=torch.randperm(n,device=device)
        ep_bce=ep_dis=ep_tot=ep_wsum=0.0           # weight-averaged loss components over the epoch
        for i in range(0,n,batch):
            idx=perm[i:i+batch]; opt.zero_grad()
            out=net(Xt[idx]); wl=wt[idx]
            bce_term=(bce(out,yt[idx])*wl).sum()/(wl.sum()+1e-9)
            loss=bce_term; dterm=torch.tensor(0.0,device=device)
            if lambda_disco>0:
                bm=ibt[idx]
                if int(bm.sum())>200:
                    Xb=Xt[idx][bm]; msb=mt[idx][bm]; ab=alpt[idx][bm]; hb=hmt[idx][bm]
                    nw=torch.ones(int(bm.sum()),device=device)
                    dloss=0.0
                    for mh in disco_masses:          # score the SAME bkg events under each fixed hypothesis
                        pm=(((ab-mh)/hb)-mu_p)/sd_p   # standardized param at hypothesis mh
                        Xbm=Xb.clone(); Xbm[:,pidx]=pm
                        scm=torch.sigmoid(net(Xbm))
                        dloss=dloss+distance_corr(scm,msb,nw)**2
                    dterm=lambda_disco*dloss/len(disco_masses)
                    loss=loss+dterm
            loss.backward(); opt.step()
            bw=float(wl.sum())                       # accumulate weighted by batch weight sum
            ep_bce+=float(bce_term)*bw; ep_dis+=float(dterm)*bw; ep_tot+=float(loss)*bw; ep_wsum+=bw
        history["epoch"].append(ep+1)
        history["bce"].append(ep_bce/ep_wsum); history["disco"].append(ep_dis/ep_wsum)
        history["total"].append(ep_tot/ep_wsum)
        if val is not None:
            net.eval()
            with torch.no_grad():
                vb=(bce(net(Xvt),yvt)*wvt).sum()/(wvt.sum()+1e-9)
            history["val_bce"].append(float(vb))
        # TODO(early stopping): could stop on val_bce plateau; currently fixed epoch count.
    net.eval()
    return net, mu, sd, history

# --------------------------------------------------------------------------- evaluation
def evaluate(net, mu, sd, device="cpu"):
    """AUC (sig vs bkg) + 125-peak sculpting R per mass hypothesis on background inclusive."""
    with torch.no_grad():
        # AUC on the training-style sample (rebuild a balanced eval set)
        Xev,yev,wev,_,_,_,_ = build_training_set(n_sig_per=4000, n_bkg=60000)
        s=torch.sigmoid(net(torch.tensor((Xev-mu)/sd,dtype=torch.float32,device=device))).cpu().numpy()
        auc=roc_auc_score(yev, s, sample_weight=wev)
        bk=_load(f"{ROOT_DIR}/All_Bkg/run3.root","inclusive")
        Hm=bk["H_m"]; wB=np.clip(bk["factor"],0,None); peak=(Hm>120)&(Hm<130); frac=wB[peak].sum()/wB.sum()
        def wq(v,wt,q): o=np.argsort(v); c=np.cumsum(wt[o])/wt.sum(); return np.interp(q,c,v[o])
        R={}
        for mh in (1,3,20):
            Xe=_feat(bk,mh); sc=torch.sigmoid(net(torch.tensor((Xe-mu)/sd,dtype=torch.float32,device=device))).cpu().numpy()
            thr=wq(sc,wB,1-0.002); p=sc>thr
            R[mh]=(wB[p&peak].sum()/wB[p].sum())/frac if wB[p].sum()>0 else float('nan')
    return auc, R

# ------------------------------------------------------------------------------- saving
def save_model(net, mu, sd, path, lambda_disco):
    torch.save({"state_dict":net.state_dict(), "mu":mu, "sd":sd,
                "features":FEATURES, "base_vars":BASE_VARS, "nfeat":NFEAT,
                "lambda_disco":lambda_disco, "arch":"ParamMLP(128,128,64)"}, path)
    print("[save]", path)

# --------------------------------------------------------------- diagnostics / plots
def _score(net, X, mu, sd, device="cpu"):
    with torch.no_grad():
        return torch.sigmoid(net(torch.tensor((X-mu)/sd,dtype=torch.float32,device=device))).cpu().numpy()

def weighted_ks_2samp(x1, w1, x2, w2):
    """Weighted two-sample KS (same approximation as run3_Za_BDT.weighted_ks_2samp)."""
    x1=np.asarray(x1,float); x2=np.asarray(x2,float); w1=np.abs(np.asarray(w1,float)); w2=np.abs(np.asarray(w2,float))
    m1=np.isfinite(x1)&np.isfinite(w1)&(w1>0); m2=np.isfinite(x2)&np.isfinite(w2)&(w2>0)
    x1,w1=x1[m1],w1[m1]; x2,w2=x2[m2],w2[m2]
    if len(x1)==0 or len(x2)==0: return np.nan, np.nan
    o1=np.argsort(x1); o2=np.argsort(x2); x1,w1=x1[o1],w1[o1]; x2,w2=x2[o2],w2[o2]
    cdf1=np.cumsum(w1)/w1.sum(); cdf2=np.cumsum(w2)/w2.sum()
    xa=np.sort(np.concatenate([x1,x2]))
    i1=np.searchsorted(x1,xa,side="right")-1; i2=np.searchsorted(x2,xa,side="right")-1
    c1=np.where(i1>=0,cdf1[np.maximum(i1,0)],0.0); c2=np.where(i2>=0,cdf2[np.maximum(i2,0)],0.0)
    stat=np.max(np.abs(c1-c2))
    n1=w1.sum()**2/(w1**2).sum(); n2=w2.sum()**2/(w2**2).sum(); en=np.sqrt(n1*n2/(n1+n2))
    return stat, stats.kstwobign.sf((en+0.12+0.11/en)*stat)

def ks_capped(x_tr, w_tr, x_va, w_va, cap=4000, reps=25, seed=0):
    """KS on fixed-size random subsamples, median over `reps` draws.
    The raw weighted KS p-value collapses to ~0 with O(1e6) events even when the train/valid
    shapes overlap (p ~ sf(sqrt(n_eff)*D), and n_eff is huge). Capping each side at `cap`
    events restores a sane statistical power so p reflects shape agreement, not sample size."""
    n_tr, n_va = len(x_tr), len(x_va)
    if n_tr == 0 or n_va == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed); ds = []; ps = []
    for _ in range(reps):
        itr = rng.choice(n_tr, min(cap, n_tr), replace=False)
        iva = rng.choice(n_va, min(cap, n_va), replace=False)
        d, p = weighted_ks_2samp(x_tr[itr], w_tr[itr], x_va[iva], w_va[iva])
        if np.isfinite(p): ds.append(d); ps.append(p)
    return (np.median(ds) if ds else np.nan), (np.median(ps) if ps else np.nan)

def plot_loss(history, out):
    """Loss vs epoch: total / BCE / DisCo (train) + validation BCE -> divergence = overtraining."""
    ep=history["epoch"]; fig,ax=plt.subplots(figsize=(8,6))
    ax.plot(ep, history["total"], "-",  color="k",         lw=2, label="train total")
    ax.plot(ep, history["bce"],   "--", color="tab:blue",  lw=2, label="train BCE")
    if any(d>0 for d in history["disco"]):
        ax.plot(ep, history["disco"], "-.", color="tab:red", lw=2, label=r"$\lambda\cdot$DisCo")
    if history["val_bce"]:
        ax.plot(ep, history["val_bce"], "-", color="tab:green", lw=2, label="val BCE")
    ax.set_xlabel("Epoch",fontsize=16); ax.set_ylabel("Loss",fontsize=16)
    ax.tick_params(labelsize=13,direction="in",top=True,right=True); ax.legend(fontsize=13,frameon=False)
    fig.tight_layout(); fig.savefig(out,dpi=200); plt.close(fig); print("[plot]", out)

def plot_overtraining(net, mu, sd, train_set, val_set, out, device="cpu", bins=40, ks_cap=4000):
    """NN analogue of run3_Za_BDT.compare_train_test: train (filled) vs validation (points)
    NN-score distributions for sig/bkg, with KS p-values computed on capped subsamples
    (ks_cap events/side) so the p-value is not crushed by the ~1e6 sample size."""
    (Xtr,ytr,wtr)=train_set; (Xva,yva,wva)=val_set
    str_=_score(net,Xtr,mu,sd,device); sva=_score(net,Xva,mu,sd,device)
    sig_tr,bk_tr = str_[ytr>0.5], str_[ytr<0.5]; sig_va,bk_va = sva[yva>0.5], sva[yva<0.5]
    wsig_tr,wbk_tr = wtr[ytr>0.5], wtr[ytr<0.5]; wsig_va,wbk_va = wva[yva>0.5], wva[yva<0.5]
    _,p_sig=ks_capped(sig_tr,wsig_tr,sig_va,wsig_va,cap=ks_cap,seed=1)
    _,p_bk =ks_capped(bk_tr,wbk_tr,bk_va,wbk_va,cap=ks_cap,seed=2)
    fig=plt.figure(figsize=(8,6)); rng=(0.0,1.0)
    plt.hist(0,color="w",range=rng,bins=bins,density=True,label=f"Sig K-S p = {p_sig:.3g} (n={ks_cap})")
    plt.hist(0,color="w",range=rng,bins=bins,density=True,label=f"Bkg K-S p = {p_bk:.3g} (n={ks_cap})")
    plt.hist(sig_tr,color="r",alpha=0.5,range=rng,bins=bins,density=True,weights=wsig_tr,histtype="stepfilled",label="Sig (Train)")
    plt.hist(bk_tr ,color="b",alpha=0.5,range=rng,bins=bins,density=True,weights=wbk_tr ,histtype="stepfilled",label="Bkg (Train)")
    for vals,wv,c,lab in ((sig_va,wsig_va,"r","Sig (Valid)"),(bk_va,wbk_va,"b","Bkg (Valid)")):
        h,edges=np.histogram(vals,bins=bins,range=rng,density=True,weights=wv)
        ctr=(edges[:-1]+edges[1:])/2; scale=len(vals)/max(h.sum(),1e-9); err=np.sqrt(np.clip(h*scale,0,None))/scale
        plt.errorbar(ctr,h,yerr=err,fmt=".",c=c,label=lab,markersize=8)
    plt.xlabel("NN Score",fontsize=18); plt.ylabel("Arbitrary Units",fontsize=18)
    plt.tick_params(labelsize=14,direction="in",top=True,right=True)
    plt.legend(loc="upper center",fontsize=13,frameon=False)
    plt.tight_layout(); fig.savefig(out,dpi=200); plt.close(fig)
    print(f"[plot] {out}  (KS p: sig={p_sig:.3g} bkg={p_bk:.3g})")
    return p_sig, p_bk

# Companion: apply_nn.py loads this .pt, scores each bdt_scored p2root file (per-event param
# from the mass hypothesis), and writes NN_Score_mA_M{ma} branches -- mirroring
# Parque2Root_BDT.py -- so flashgg (CMSSW, no torch) only needs to read the branches.

# --------------------------------------------------------------------------------- main
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--lambda-disco", type=float, default=0.0)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch", type=int, default=4000)
    ap.add_argument("--n-sig-per", type=int, default=8000)
    ap.add_argument("--n-bkg", type=int, default=150000)
    ap.add_argument("--disco-masses", default="1,2,3",
                    help="comma-separated low-mass hypotheses for per-mass DisCo")
    ap.add_argument("--out", default="model_Za_NN_run3.pt")
    ap.add_argument("--plot-dir", default="plots_nn", help="dir for loss + overtraining PDFs")
    ap.add_argument("--seed", type=int, default=1, help="training seed (DisCo R(mA=1) has run-to-run variance)")
    ap.add_argument("--device", default=("cuda" if torch.cuda.is_available() else "cpu"))
    a=ap.parse_args()
    disco_masses=tuple(int(x) for x in str(a.disco_masses).split(","))
    print(f"[cfg] device={a.device} lambda_disco={a.lambda_disco} disco_masses={disco_masses} epochs={a.epochs} feat={NFEAT}")
    X,y,w,mass,isb,alp,hm = build_training_set(a.n_sig_per, a.n_bkg)
    print(f"[data] sig+bkg={len(X)}  bkg={int(isb.sum())}")
    # Independent validation set (validation tree) for the loss curve + overtraining check.
    Xv,yv,wv,_,_,_,_ = build_training_set(max(a.n_sig_per//2,1000), a.n_bkg//2, tree="validation")
    net,mu,sd,history = train(X,y,w,mass,isb,alp,hm, lambda_disco=a.lambda_disco, disco_masses=disco_masses,
                              epochs=a.epochs, batch=a.batch, device=a.device, seed=a.seed, val=(Xv,yv,wv))
    auc,R = evaluate(net,mu,sd,device=a.device)
    print(f"[result] AUC={auc:.4f}  R(mA=1)={R[1]:.2f}  R(mA=3)={R[3]:.2f}  R(mA=20)={R[20]:.2f}")
    save_model(net,mu,sd,a.out,a.lambda_disco)
    pdir=Path(a.plot_dir); pdir.mkdir(parents=True, exist_ok=True)
    tag=Path(a.out).stem
    plot_loss(history, str(pdir/f"loss_{tag}.pdf"))
    plot_overtraining(net,mu,sd,(X,y,w),(Xv,yv,wv), str(pdir/f"overtrain_{tag}.pdf"), device=a.device)

if __name__=="__main__":
    main()
