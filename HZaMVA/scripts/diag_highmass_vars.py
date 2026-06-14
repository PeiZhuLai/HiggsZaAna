#!/usr/bin/env python3
"""Per-variable sculpting diagnostic for the HIGH-mass model, evaluated at a chosen high m_a
(default 30, the point with R=1.21). For each of the 15 hypothesis-independent inputs:
  - |corr(feature, m_llgammagamma)| in BACKGROUND  -> sculpting potential (mass-independent)
  - single-feature separation AUC (signal m_a vs background) -> discrimination AT THAT MASS
A variable worth dropping from the high-mass model = HIGH mass-correlation AND LOW discrimination
at high m_a (a "pure burden" that only sculpts). A high-corr but high-sep variable is kept: at high
m_a the two photons are resolved, so sub-leading-photon / kinematic variables recover discrimination.
Runs in higgs-alp-ana. NO model, NO retraining.
"""
import sys
import numpy as np
from hza_features import _load, _feat, BASE_VARS, ROOT_DIR

NAMES = BASE_VARS + ["pho_pt_asym"]   # 15 hypothesis-independent features (param excluded)
MASSES = [int(x) for x in (sys.argv[1:] or [30, 25])]

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); vals=vals[o]; y=y[o]; w=w[o]
    is_b = y < 0.5; cum_b = np.cumsum(w*is_b); below_b = cum_b - w*is_b
    Ws = w[~is_b].sum(); Wb = w[is_b].sum()
    return (w[~is_b]*below_b[~is_b]).sum()/(Ws*Wb+1e-12)

def wpearson(x, y, w):
    w = np.clip(w, 0, None); mx = np.average(x, weights=w); my = np.average(y, weights=w)
    cov = np.average((x-mx)*(y-my), weights=w)
    sx = np.sqrt(np.average((x-mx)**2, weights=w)); sy = np.sqrt(np.average((y-my)**2, weights=w))
    return cov/(sx*sy+1e-12)

# background is mass-independent -> |corr(H_m)| computed once
bk = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")
Hmb = bk["H_m"]; wb = np.clip(bk["factor"], 0, None)
Xb = _feat(bk, 0)[:, :15]
corr = {name: abs(wpearson(Xb[:, i], Hmb, wb)) for i, name in enumerate(NAMES)}

for mA in MASSES:
    sg = _load(f"{ROOT_DIR}/mA_M{mA}/run3.root", "inclusive")
    ws = np.clip(sg["factor"], 0, None); Xs = _feat(sg, mA)[:, :15]
    rows = []
    for i, name in enumerate(NAMES):
        vals = np.concatenate([Xs[:, i], Xb[:, i]]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
        w = np.concatenate([ws*(wb.sum()/ws.sum()), wb])
        sep = max(weighted_auc(vals, y, w), 1-weighted_auc(vals, y, w))
        rows.append((name, corr[name], sep))
    rows.sort(key=lambda t: -t[1])
    print(f"\n================  m_a = {mA} GeV  ================")
    print(f"{'feature':<26}{'|corr(H_m)|':>12}{'sep AUC':>10}{'verdict':>22}")
    for name, r, sep in rows:
        if r > 0.10 and sep < 0.58:
            v = "<== DROP candidate"          # high corr, weak discrimination -> pure sculpt burden
        elif r > 0.10 and sep >= 0.62:
            v = "high-corr but discriminating" # keep: earns its sculpting
        elif r < 0.06 and sep > 0.60:
            v = "safe + useful"
        else:
            v = ""
        print(f"{name:<26}{r:>12.3f}{sep:>10.3f}{v:>22}")
