#!/usr/bin/env python3
"""DisCo proxy step 2/2 (hza_NFlow, has torch): train parametrized MLP with
loss = BCE + lambda*dCorr^2(score_bkg, m_llgammagamma_bkg). Report 125-peak R (fixed bkg
eff) and signal-vs-bkg AUC per lambda. NN = DisCo ceiling, NOT the XGBoost analysis."""
import numpy as np, torch, torch.nn as nn
torch.manual_seed(1); np.random.seed(0)

d=np.load("/tmp/disco_data.npz", allow_pickle=True)
X=d["X"]; y=d["y"]; w=d["w"]; mass=d["mass"]; isb=d["isb"].astype(bool)
NFEAT=X.shape[1]
mu,sd=X.mean(0),X.std(0)+1e-9; Xn=(X-mu)/sd
bk_base=d["bk_base"]; bk_Hm=d["bk_Hm"]; bk_ALPm=d["bk_ALPm"]; wB=d["bk_w"]
peak=(bk_Hm>120)&(bk_Hm<130); frac=wB[peak].sum()/wB.sum()

def eval_feat(mh):                      # rebuild 16 feat for bkinc (BASE14 + asym + param)
    asym=(bk_base[:,0]-bk_base[:,4])/(bk_base[:,0]+bk_base[:,4]+1e-6)
    param=(bk_ALPm-mh)/bk_Hm
    return (np.column_stack([bk_base, asym, param])-mu)/sd

def wq(v,wt,q): o=np.argsort(v); c=np.cumsum(wt[o])/wt.sum(); return np.interp(q,c,v[o])
def R_at_eff(score,eff):
    thr=wq(score,wB,1-eff); p=score>thr; pw=wB[p].sum()
    return (wB[p&peak].sum()/pw)/frac if pw>0 else np.nan
def auc_np(score,label,n=40000):
    i=np.random.choice(len(score),min(n,len(score)),replace=False); s,l=score[i],label[i]
    o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(1,len(s)+1)
    npos=l.sum(); nneg=len(l)-npos
    return (r[l==1].sum()-npos*(npos+1)/2)/(npos*nneg)

def dcorr(a_var,b_var,nw):
    a=(a_var.view(-1,1)-a_var.view(1,-1)).abs(); b=(b_var.view(-1,1)-b_var.view(1,-1)).abs()
    am=(a*nw).mean(1); A=a-am.view(-1,1)-am.view(1,-1)+(am*nw).mean()
    bm=(b*nw).mean(1); B=b-bm.view(-1,1)-bm.view(1,-1)+(bm*nw).mean()
    return ((A*B*nw).mean(1)*nw).mean()/torch.sqrt(((A*A*nw).mean(1)*nw).mean()*((B*B*nw).mean(1)*nw).mean()+1e-12)

class Net(nn.Module):
    def __init__(s,n): super().__init__(); s.f=nn.Sequential(nn.Linear(n,64),nn.ReLU(),nn.Linear(64,64),nn.ReLU(),nn.Linear(64,1))
    def forward(s,x): return s.f(x).squeeze(-1)

Xt=torch.tensor(Xn,dtype=torch.float32); yt=torch.tensor(y,dtype=torch.float32)
wt=torch.tensor(w,dtype=torch.float32); mt=torch.tensor((mass-125)/30,dtype=torch.float32)
ibt=torch.tensor(isb); BS=4000; EP=8; bce=nn.BCEWithLogitsLoss(reduction='none'); n=len(Xt)

def run(lam):
    torch.manual_seed(1); net=Net(NFEAT); opt=torch.optim.Adam(net.parameters(),1e-3)
    for ep in range(EP):
        perm=torch.randperm(n)
        for i in range(0,n,BS):
            idx=perm[i:i+BS]; opt.zero_grad(); out=net(Xt[idx]); wl=wt[idx]
            loss=(bce(out,yt[idx])*wl).sum()/(wl.sum()+1e-9)
            if lam>0:
                bm=ibt[idx]
                if int(bm.sum())>200:
                    sc=torch.sigmoid(out[bm]); ms=mt[idx][bm]; nw=torch.ones_like(sc)
                    loss=loss+lam*dcorr(sc,ms,nw)**2
            loss.backward(); opt.step()
    net.eval()
    with torch.no_grad():
        auc=auc_np(torch.sigmoid(net(Xt)).numpy(), y)
        R={mh:R_at_eff(torch.sigmoid(net(torch.tensor(eval_feat(mh),dtype=torch.float32))).numpy(),0.002) for mh in (1,3,20)}
    return auc,R

print(f"{'lambda':>7} {'AUC':>7} | {'R(mA=1)':>9} {'R(mA=3)':>9} {'R(mA=20)':>9}")
print("-"*52)
for lam in (0.0,10.0,50.0,150.0):
    auc,R=run(lam); print(f"{lam:>7.0f} {auc:>7.4f} | {R[1]:>9.2f} {R[3]:>9.2f} {R[20]:>9.2f}")
print("-"*52)
print("R>1=sculpted. lambda=0 = NN baseline (no DisCo).")
print("XGBoost pho-pT/H_m model: R(mA=1)=1.60, AUC~0.966 (different model class).")
