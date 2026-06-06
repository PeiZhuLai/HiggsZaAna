#!/usr/bin/env python3
"""Per-variable low-mA sculpting diagnostic (no model, no GPU; runs in higgs-alp-ana or hza_NN).
For each input feature, measure:
  - |corr(feature, m_llgammagamma)| in BACKGROUND  -> how much it can sculpt the 125 peak
  - single-feature separation AUC (signal mA=1 vs background) -> how much it discriminates
Sculpting drivers = high mass-correlation; safe/useful = high separation, low mass-correlation.
Also a weighted distance-correlation (nonlinear) with m_llgg for the top suspects.
"""
import numpy as np
from hza_features import _load, _feat, BASE_VARS, ROOT_DIR

NAMES = BASE_VARS + ["pho_pt_asym"]                 # 15 hypothesis-independent features (drop param)

def weighted_auc(vals, y, w):
    """Single-feature weighted AUC via rank formula (no sklearn; EOS .so imports are flaky)."""
    o = np.argsort(vals, kind="mergesort"); vals=vals[o]; y=y[o]; w=w[o]
    is_b = y < 0.5; cum_b = np.cumsum(w*is_b)                 # bkg weight at-or-below each point
    below_b = cum_b - w*is_b                                  # strictly-below bkg weight
    Ws = w[~is_b].sum(); Wb = w[is_b].sum()
    auc = (w[~is_b]*below_b[~is_b]).sum()/(Ws*Wb+1e-12)
    return auc

def wpearson(x, y, w):
    w = np.clip(w, 0, None); mx = np.average(x, weights=w); my = np.average(y, weights=w)
    cov = np.average((x-mx)*(y-my), weights=w)
    sx = np.sqrt(np.average((x-mx)**2, weights=w)); sy = np.sqrt(np.average((y-my)**2, weights=w))
    return cov/(sx*sy+1e-12)

def dcorr(a, b, w, cap=3000, seed=0):              # weighted distance correlation on a subsample
    rng = np.random.default_rng(seed); n = len(a)
    idx = rng.choice(n, min(cap, n), replace=False); a = a[idx]; b = b[idx]; w = np.clip(w[idx],0,None)
    nw = w/ w.mean()
    A = np.abs(a[:,None]-a[None,:]); B = np.abs(b[:,None]-b[None,:])
    Am = (A*nw).mean(1); A = A-Am[:,None]-Am[None,:]+(Am*nw).mean()
    Bm = (B*nw).mean(1); B = B-Bm[:,None]-Bm[None,:]+(Bm*nw).mean()
    num=((A*B*nw).mean(1)*nw).mean(); da=((A*A*nw).mean(1)*nw).mean(); db=((B*B*nw).mean(1)*nw).mean()
    return np.sqrt(max(num,0)/ (np.sqrt(da*db)+1e-12))

bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")
sg = _load(f"{ROOT_DIR}/mA_M1/run3.root", "inclusive")
Hmb = bk["H_m"]; wb = np.clip(bk["factor"], 0, None); ws = np.clip(sg["factor"], 0, None)
Xb = _feat(bk, 1)[:, :15]; Xs = _feat(sg, 1)[:, :15]

rows = []
for i, name in enumerate(NAMES):
    r = abs(wpearson(Xb[:, i], Hmb, wb))
    vals = np.concatenate([Xs[:, i], Xb[:, i]]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws*(wb.sum()/ws.sum()), wb])
    auc = weighted_auc(vals, y, w); sep = max(auc, 1-auc)
    rows.append((name, r, sep))

rows.sort(key=lambda t: -t[1])                      # sort by mass-correlation (sculpting driver) desc
print(f"{'feature':<26}{'|corr(H_m)|':>12}{'sep AUC(mA1)':>14}{'sep-corr':>10}")
for name, r, sep in rows:
    flag = "  <== sculpt driver" if (r > 0.10 and sep < 0.62) else ("  (safe/useful)" if (r < 0.06 and sep > 0.6) else "")
    print(f"{name:<26}{r:>12.3f}{sep:>14.3f}{sep-r:>10.3f}{flag}")

print("\nNonlinear check (weighted dCorr with m_llgg, top-5 by |corr|):")
for name, r, sep in rows[:5]:
    i = NAMES.index(name)
    print(f"  {name:<26} dCorr={dcorr(Xb[:, i], Hmb, wb):.3f}")
