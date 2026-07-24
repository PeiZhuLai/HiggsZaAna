#!/usr/bin/env python3
"""Combinatorial search over low-correlation variables for the dedicated mA=1,2,3 BDT:
find the subset giving R(mA=1,2,3) closest to 1 (flat background, no sculpting) at the
highest AUC. CPU only (higgs-alp-ana); data loaded ONCE, ~119 quick XGBoost fits.
Pool = the 7 diagnosed low-|corr(m_llgg)| / useful variables (param always added).
"""
import itertools, numpy as np, xgboost as xgb
from hza_features import _load, _feat, BASE_VARS, _sideband_reweight, signal_weight, ROOT_DIR

LOWM = [1, 2, 3]
FULL16 = BASE_VARS + ["pho_pt_asym", "param"]
POOL = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "var_PtaOverMh", "H_pt_oHm", "pho2R9", "pho1R9", "pho1Pt_oHm"]

def weighted_auc(vals, y, w):
    o = np.argsort(vals, kind="mergesort"); y = y[o]; w = w[o]; is_b = y < 0.5
    below_b = np.cumsum(w * is_b) - w * is_b
    return (w[~is_b] * below_b[~is_b]).sum() / (w[~is_b].sum() * w[is_b].sum() + 1e-12)

def wq(v, wt, q): o = np.argsort(v); c = np.cumsum(wt[o]) / wt.sum(); return np.interp(q, c, v[o])

print("[load] signal mA=1,2,3 + All_Bkg ...")
sig_tr = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "train", 40000, extra=["Z_m", "event"]) for m in LOWM}
sig_in = {m: _load(f"{ROOT_DIR}/mA_M{m}/run3.root", "inclusive", extra=["Z_m", "event"]) for m in LOWM}
bk_tr = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "train", 350000, extra=["Z_m", "event"])
bk_in = _load(f"{ROOT_DIR}/All_Bkg/run3.root", "inclusive")

Xs16 = np.vstack([_feat(sig_tr[m], m) for m in LOWM])
ws = np.concatenate([np.clip(signal_weight(sig_tr[m], m), 0, None) for m in LOWM])  # signal: factor x sideband reweight
bhyp = np.random.default_rng(1).choice(LOWM, size=len(bk_tr["H_m"]))
Xb16 = _feat(bk_tr, 0); Xb16[:, FULL16.index("param")] = (bk_tr["ALP_m"] - bhyp) / bk_tr["H_m"]
wb = np.clip(bk_tr["factor"] * _sideband_reweight(bk_tr), 0, None)
y = np.concatenate([np.ones(len(Xs16)), np.zeros(len(Xb16))]); w = np.concatenate([ws * (wb.sum()/ws.sum()), wb])
X16 = np.vstack([Xs16, Xb16])
Hmb = bk_in["H_m"]; wBe = np.clip(bk_in["factor"], 0, None)
peak = (Hmb > 120) & (Hmb < 130); frac = wBe[peak].sum() / wBe.sum()
Xb_in16 = {m: _feat(bk_in, m) for m in LOWM}; Xs_in16 = {m: _feat(sig_in[m], m) for m in LOWM}

def evaluate(idx):
    clf = xgb.XGBClassifier(max_depth=4, n_estimators=200, learning_rate=0.1, subsample=0.8,
                            min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
    clf.fit(X16[:, idx], y, sample_weight=w)
    R = []
    for m in LOWM:
        s = clf.predict_proba(Xb_in16[m][:, idx])[:, 1]; thr = wq(s, wBe, 1 - 0.002); p = s > thr
        R.append((wBe[p & peak].sum() / wBe[p].sum()) / frac if wBe[p].sum() > 0 else np.nan)
    sv, yv, wv = [], [], []
    for m in LOWM:
        ss = clf.predict_proba(Xs_in16[m][:, idx])[:, 1]; sb = clf.predict_proba(Xb_in16[m][:, idx])[:, 1]
        wss = np.clip(signal_weight(sig_in[m], m), 0, None)
        sv += [ss, sb]; yv += [np.ones(len(ss)), np.zeros(len(sb))]; wv += [wss * (wBe.sum()/wss.sum()), wBe]
    return weighted_auc(np.concatenate(sv), np.concatenate(yv), np.concatenate(wv)), R

pidx = FULL16.index("param")
rows = []
combos = [c for k in range(2, 7) for c in itertools.combinations(POOL, k)]
print(f"[search] {len(combos)} subsets (sizes 2-6) ...")
for keep in combos:
    idx = [FULL16.index(n) for n in keep] + [pidx]
    auc, R = evaluate(idx)
    maxdev = max(abs(r - 1) for r in R)
    rows.append((keep, len(keep), auc, R[0], R[1], R[2], maxdev))

# Rank: among subsets where all R in [1-tol, 1+tol], highest AUC. Print several tolerances.
def show(tol, n=10):
    sel = sorted([r for r in rows if r[6] <= tol], key=lambda r: -r[2])[:n]
    print(f"\n=== all |R-1| <= {tol}  (best AUC first) ===")
    print(f"{'AUC':>8}{'R1':>7}{'R2':>7}{'R3':>7}{'nvar':>5}  variables")
    for keep, nv, auc, r1, r2, r3, md in sel:
        print(f"{auc:>8.4f}{r1:>7.2f}{r2:>7.2f}{r3:>7.2f}{nv+1:>5}  {'+'.join(keep)}")

for tol in (0.20, 0.30):
    show(tol)
best = min(rows, key=lambda r: (r[6], -r[2]))   # most-flat, tie-break AUC
print(f"\n>>> most-flat overall: AUC={best[2]:.4f} R={best[3]:.2f}/{best[4]:.2f}/{best[5]:.2f} "
      f"maxdev={best[6]:.2f} :: {'+'.join(best[0])}")
