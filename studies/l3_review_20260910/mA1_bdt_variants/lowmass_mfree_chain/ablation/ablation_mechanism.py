#!/usr/bin/env python3
"""Which input lets the win_mfree low-mass BDT suppress the m_a = 1 overlapping-shower population?

Population: mA_M1 signal, test tree (30%), 5 Run 3 eras, 95 < m_llgg < 180, weight factor/0.3.
"Overlap" = gen dR(gamma,gamma) < 0.05 (the second m_llgg population near 140 GeV).
Both models (production 5-variable, win_mfree) are evaluated at their working points
(production 0.963, win_mfree 0.9645 = sandbox WP; the two give the same data sideband count).

  A. swap      for overlap events, replace ONE input by the value of a randomly paired non-overlap
               in-window (115-133) signal event; pass rate of the overlap events at the WP, per input.
               A large rise means that input is what keeps (production) or rejects (win_mfree) them.
  B. SHAP      xgboost pred_contribs: mean contribution (log-odds) per input, overlap vs non-overlap.
  C. drop-one  retrain win_mfree without one input (same recipe, train_variants.build), cut at the
               same data sideband count as production; overlap share and window fraction.
  D. inputs    weighted medians of each input, overlap vs non-overlap in-window.
Env: higgs-alp-ana. Output: ablation_mechanism.{txt,json}, models_dropone/
"""
import json
import os
import pickle
import sys

import numpy as np
import uproot
import xgboost as xgb

HERE = os.path.dirname(os.path.abspath(__file__))
VAR = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, VAR)
from train_variants import build  # noqa: E402

SB = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
PROD = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_lowmass_run3.pkl"
NEW = os.path.join(VAR, "model_win_mfree.pkl")
F_PROD = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]
F_NEW = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt", "param_rel"]
WP = {"prod": 0.963, "new": 0.9645}
RAW = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "ALP_mass", "H_mass", "factor"]


def load(sample, tree, extra=()):
    out = {}
    for e in ERAS:
        a = uproot.open(f"{SB}/{sample}/{e}.root")[tree].arrays(RAW + list(extra), library="np")
        for k, v in a.items():
            out.setdefault(k, []).append(v)
    a = {k: np.concatenate(v) for k, v in out.items()}
    ok = (a["H_mass"] > 95) & (a["H_mass"] < 180) & np.isfinite(a["factor"].astype(float))
    return {k: v[ok] for k, v in a.items()}


def X(a, feats, m=1.0):
    d = {"param": (a["ALP_mass"] - m) / a["H_mass"], "param_rel": (a["ALP_mass"] - m) / m,
         "pho1Pt": a["pho1Pt_oHm"] * a["H_mass"]}
    return np.column_stack([d[f] if f in d else a[f].astype(float) for f in feats])


def wmed(x, w):
    o = np.argsort(x); c = np.cumsum(w[o]); return float(x[o][np.searchsorted(c, 0.5 * c[-1])])


def main():
    sig = load("mA_M1", "test", extra=["GenALP_dR_gg"])
    w = sig["factor"].astype(float) / 0.3
    gdr = np.array([x[0] if len(x) else np.nan for x in sig["GenALP_dR_gg"]], float)
    ov = gdr < 0.05
    m = sig["H_mass"]
    inwin = (m > 115) & (m < 133)
    ref = (~ov) & inwin
    data = load("Data", "inclusive")
    dm = data["H_mass"]; dsb = ((dm > 95) & (dm < 115)) | ((dm > 135) & (dm < 180))
    rng = np.random.default_rng(11)
    lines, res = [], {}

    lines.append("D. input medians (weighted), overlap vs non-overlap in-window signal")
    for f in sorted(set(F_PROD + F_NEW)):
        xo, xr = X(sig, [f])[ov, 0], X(sig, [f])[ref, 0]
        lines.append(f"   {f:24s} overlap {wmed(xo, w[ov]):10.4f}   non-overlap in-window {wmed(xr, w[ref]):10.4f}")
    lines.append(f"   overlap: m_llgg median {wmed(m[ov], w[ov]):.1f} GeV, share of signal {w[ov].sum() / w.sum():.3f}")

    for tag, path, feats in (("prod", PROD, F_PROD), ("new", NEW, F_NEW)):
        clf = pickle.load(open(path, "rb"))
        Xs = X(sig, feats)
        s0 = clf.predict_proba(Xs)[:, 1]
        cut = WP[tag]
        p0_ov = w[ov & (s0 > cut)].sum() / w[ov].sum()
        p0_ref = w[ref & (s0 > cut)].sum() / w[ref].sum()
        lines.append(f"\n[{tag}] WP {cut}: pass rate overlap {p0_ov:.3f}, non-overlap in-window {p0_ref:.3f}; "
                     f"overlap share after cut {w[ov & (s0 > cut)].sum() / w[s0 > cut].sum():.3f}")
        # A. swap one input for overlap events
        lines.append(f"  A. swap one input of the overlap events with a random non-overlap in-window event")
        idx_ref = np.flatnonzero(ref); idx_ov = np.flatnonzero(ov)
        res[tag] = {"pass_overlap": p0_ov, "pass_ref": p0_ref, "swap": {}, "shap": {}}
        for j, f in enumerate(feats):
            Xc = Xs[idx_ov].copy()
            Xc[:, j] = Xs[rng.choice(idx_ref, len(idx_ov)), j]
            s = clf.predict_proba(Xc)[:, 1]
            p = w[idx_ov][s > cut].sum() / w[idx_ov].sum()
            res[tag]["swap"][f] = p
            lines.append(f"     {f:24s} overlap pass rate {p0_ov:.3f} -> {p:.3f}  ({p / p0_ov if p0_ov > 0 else np.nan:.2f}x)")
        # all inputs swapped except one (which single input alone carries the decision)
        # B. SHAP
        booster = clf.get_booster()
        contrib = booster.predict(xgb.DMatrix(Xs), pred_contribs=True)
        lines.append("  B. mean SHAP contribution (log-odds): overlap | non-overlap in-window | difference")
        for j, f in enumerate(feats):
            a_ov = np.average(contrib[ov, j], weights=w[ov]); a_ref = np.average(contrib[ref, j], weights=w[ref])
            res[tag]["shap"][f] = (float(a_ov), float(a_ref))
            lines.append(f"     {f:24s} {a_ov:+8.3f} | {a_ref:+8.3f} | {a_ov - a_ref:+8.3f}")

    # C. drop-one retraining of win_mfree
    lines.append("\nC. win_mfree retrained without one input (cut at the production data sideband count)")
    prod = pickle.load(open(PROD, "rb"))
    n_target = int((prod.predict_proba(X(data, F_PROD))[:, 1][dsb] > WP["prod"]).sum())
    lines.append(f"   data sideband count target (production at {WP['prod']}): {n_target}")
    os.makedirs(os.path.join(HERE, "models_dropone"), exist_ok=True)
    res["dropone"] = {}
    for drop in [None] + F_NEW:
        feats = [f for f in F_NEW if f != drop]
        tag = "none" if drop is None else drop
        mp = os.path.join(HERE, "models_dropone", f"win_mfree_minus_{tag}.pkl")
        if not os.path.exists(mp):
            Xt, yt, wt = build("train", feats, True)
            c = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8,
                                  min_child_weight=5, eval_metric="logloss", n_jobs=4, verbosity=0)
            c.fit(Xt, yt, sample_weight=wt); pickle.dump(c, open(mp, "wb"))
        c = pickle.load(open(mp, "rb"))
        sd = c.predict_proba(X(data, feats))[:, 1][dsb]
        cut = np.sort(sd)[::-1][n_target - 1] if len(sd) >= n_target else 0.0
        s = c.predict_proba(X(sig, feats))[:, 1]
        p = s >= cut
        r = {"cut": float(cut), "overlap_share": float(w[p & ov].sum() / w[p].sum()),
             "win_frac": float(w[p & inwin].sum() / w[p].sum()), "S_win": float(w[p & inwin].sum()),
             "pass_overlap": float(w[p & ov].sum() / w[ov].sum())}
        res["dropone"][tag] = r
        lines.append(f"   minus {tag:24s} cut {cut:.4f}: overlap share {r['overlap_share']:.3f}  window fraction "
                     f"{r['win_frac']:.3f}  S(115-133) {r['S_win']:.1f}  overlap pass rate {r['pass_overlap']:.3f}")
        print(lines[-1], flush=True)
    open(os.path.join(HERE, "ablation_mechanism.txt"), "w").write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(HERE, "ablation_mechanism.json"), "w"), indent=1, default=float)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
