#!/usr/bin/env python3
"""One candidate of the HZa BDT input-variable search: train (production recipe) and evaluate.

Config (one JSON object per candidate): {"id": str, "model": "low"|"high", "feats": [...],
  "win_mA1": bool (m_a = 1 signal trained on 115-133 GeV only, weight rescaled), "plane": bool (background
  training weights flattened in m_llgg, 1 GeV bins), "stored_prod": bool (no training: evaluate the stored
  production scores; reference for the metric implementation), "ntrees": int (default 300)}.
The mass-hypothesis variable param_rel = (m_gg - m_a)/m_a is always added (low and high models share the
definition). Derived features: pho_pt_asym (GeV pT), pho{1,2}_relEerr = energyErr / (pT cosh eta).
Production-only features (baseline candidates): pho1Pt_oHm, pho2Pt_oHm, var_PtaOverMh, H_pt_oHm, param
(= (m_gg - m_a)/m_llgg), pho_pt_asym_oHm.

Per anchor mass of the model, at the cut giving the production data-sideband count:
  R, R_lo, R_hi   DY f(120-130), f(105-115), f(135-145) after/before the cut, 100 Poisson bootstrap sd
  dy_chi2/ndf     DY passing vs inclusive normalized shape, 5 GeV bins 95-180 (MC stat of the passing set)
  data_chi2/ndf   data sideband passing vs preselected shape, bands 95-105,105-115,135-145,145-160,160-180
  S_win, win_frac, overlap, sec      signal test tree (weight factor/0.3): yield in 115-133, its share of
                  95-180, share with gen dR(gg) < 0.05, S(135-160)/S(115-133); presel values for reference
  Z               S_win / sqrt(B_win), B_win = data sideband count x DY(115-133)/DY(sidebands)
Env: any python with numpy + xgboost (IHEP hzgdna, CERN higgs-alp-ana).
Usage: search_job.py --samples <dir> --config <file.jsonl> --index <i> --out <dir>
"""
import argparse
import json
import math
import os
import time

import numpy as np

LOW = [1, 2, 3]
HIGH = [4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
WIN = (115.0, 133.0)
BANDS = [(95, 105), (105, 115), (135, 145), (145, 160), (160, 180)]
EDGES = np.arange(95, 181, 5.0)
NBOOT = 100


class Sample(object):
    def __init__(self, path):
        self.z = np.load(path)
        self.c = {}

    def __getitem__(self, k):
        if k not in self.c:
            self.c[k] = self.z[k].astype(np.float64)
        return self.c[k]

    def __len__(self):
        return len(self["H_m"])


def column(s, f, mhyp):
    if f == "param_rel":
        return (s["ALP_m"] - mhyp) / mhyp
    if f == "param":
        return (s["ALP_m"] - mhyp) / s["H_m"]
    if f == "pho_pt_asym":
        return (s["pho1Pt"] - s["pho2Pt"]) / (s["pho1Pt"] + s["pho2Pt"] + 1e-6)
    if f == "pho_pt_asym_oHm":
        return (s["pho1Pt_oHm"] - s["pho2Pt_oHm"]) / (s["pho1Pt_oHm"] + s["pho2Pt_oHm"] + 1e-6)
    if f in ("pho1_relEerr", "pho2_relEerr"):
        w = "lead" if f.startswith("pho1") else "sublead"
        pt = s["pho1Pt"] if w == "lead" else s["pho2Pt"]
        return s["ALP_%s_photon_energyErr" % w] / (pt * np.cosh(s["ALP_%s_photon_eta" % w]))
    return s[f]


def matrix(s, feats, mhyp):
    return np.column_stack([column(s, f, mhyp) for f in feats])


def inr(x, lo, hi):
    return (x > lo) & (x < hi)


def frac(w, m, lo, hi, sel):
    t = w[sel].sum()
    return w[sel & inr(m, lo, hi)].sum() / t if t > 0 else float("nan")


def shape_chi2(m, w, passed, edges):
    h0, _ = np.histogram(m, edges, weights=w)
    h1, _ = np.histogram(m[passed], edges, weights=w[passed])
    e1, _ = np.histogram(m[passed], edges, weights=w[passed] ** 2)
    if h1.sum() <= 0 or h0.sum() <= 0:
        return float("nan"), 0
    n0, n1, var = h0 / h0.sum(), h1 / h1.sum(), e1 / h1.sum() ** 2
    ok = (var > 0) & (h0 > 0)
    return float(np.sum((n1[ok] - n0[ok]) ** 2 / var[ok])), int(ok.sum() - 1)


def chi2_pvalue(c2, ndf):
    if not (ndf > 0 and np.isfinite(c2)):
        return float("nan")
    # regularized upper incomplete gamma Q(ndf/2, c2/2) by series/continued fraction (no scipy needed)
    a, x = ndf / 2.0, c2 / 2.0
    if x <= 0:
        return 1.0
    if x < a + 1:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(500):
            ap += 1; d *= x / ap; s += d
            if abs(d) < abs(s) * 1e-12: break
        return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - math.lgamma(a)))
    b, c, d = x + 1 - a, 1e300, 1.0 / (x + 1 - a); h = d
    for i in range(1, 500):
        an = -i * (i - a); b += 2; d = an * d + b; d = 1.0 / d if abs(d) > 1e-300 else 1e300
        c = b + an / c; h *= d * c
        if abs(d * c - 1) < 1e-12: break
    return min(1.0, math.exp(-x + a * math.log(x) - math.lgamma(a)) * h)


SB_DY_LO = np.concatenate([np.arange(95, 115, 2.5), np.arange(135, 180, 2.5)])   # 2.5 GeV sideband bins
SB_DATA_LO = np.array([95, 100, 105, 110] + list(range(135, 180, 5)), float)     # 5 GeV sideband bins
PULL_WINDOWS = [(115, 125), (120, 130), (125, 135)]
MAXORDER = 3


def design(x, order):
    t = (np.asarray(x, float) - 137.5) / 42.5
    return np.column_stack([t ** k for k in range(order + 1)])


def glm_fit(X, y, var_scale=None):
    """Poisson GLM with log link (IRLS): mu = exp(X b); y may be non-integer (scaled MC).
    var_scale = Var(y)/mu per bin (overdispersion from MC weights), used in the fit weights and covariance.
    Returns b, cov(b), Poisson deviance (saturated); (None, None, nan) if singular."""
    vs = np.ones_like(y) if var_scale is None else var_scale
    b = np.zeros(X.shape[1]); b[0] = math.log(max(y.mean(), 1e-3))
    try:
        for _ in range(200):
            mu = np.exp(np.clip(X @ b, -50, 50))
            W = mu / vs
            z = X @ b + (y - mu) / np.maximum(mu, 1e-12)
            nb = np.linalg.solve(X.T @ (W[:, None] * X), X.T @ (W * z))
            done = np.max(np.abs(nb - b)) < 1e-9
            b = nb
            if done:
                break
        mu = np.exp(np.clip(X @ b, -50, 50))
        cov = np.linalg.inv(X.T @ ((mu / vs)[:, None] * X))
    except np.linalg.LinAlgError:
        return None, None, float("nan")
    dev = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0) - (y - mu))
    return b, cov, float(dev)


def ftest_order(xc, y, var_scale=None):
    """Lowest order 1..MAXORDER for which one more order lowers the deviance by < 3.84 (95% chi2_1)."""
    fits = [glm_fit(design(xc, o), y, var_scale) for o in range(1, MAXORDER + 1)]
    o = 1
    while o < MAXORDER and fits[o][0] is not None and \
            (fits[o - 1][0] is None or fits[o - 1][2] - fits[o][2] > 3.84):
        o += 1
    return o, fits[o - 1]


def closure_proxy(bm, bw, p, n_data_sb):
    """Gate A, cheap MC closure: passing DY scaled to the data-sideband count (as the closure pseudo-data),
    log-polynomial (F-test order) fitted to the sidebands 95-115, 135-180 in 2.5 GeV bins and extrapolated into
    the blinded 115-135. pull = (DY - fit) / sqrt(fit + MC var + fit var) per 10 GeV window: a fake peak, a dip or a
    turn-on at the window edge that a smooth function cannot follow, measured at data-level statistics."""
    m, w = bm[p], bw[p]
    nan = {"A_pulls": [float("nan")] * len(PULL_WINDOWS), "A_order": 0, "A_fit_p": float("nan")}
    sbw = w[inr(m, 95, 115) | inr(m, 135, 180)].sum()
    if sbw <= 0:
        return nan
    k = n_data_sb / sbw
    def hist(lo, ww):
        return np.array([ww[inr(m, a, a + 2.5)].sum() for a in lo])
    y, v_mc = hist(SB_DY_LO, w) * k, hist(SB_DY_LO, w ** 2) * k * k
    xc = SB_DY_LO + 1.25
    vs = 1.0 + v_mc / np.maximum(y, 1e-6)     # Var = mu (data-level Poisson) + MC variance
    order, (b, cov, dev) = ftest_order(xc, y, vs)
    if b is None:
        return nan
    mu = np.exp(design(xc, order) @ b)
    out = {"A_order": order, "A_k": float(k), "A_fit_chi2": float(np.sum((y - mu) ** 2 / (mu + v_mc))),
           "A_fit_ndf": len(y) - order - 1}
    out["A_fit_p"] = chi2_pvalue(out["A_fit_chi2"], out["A_fit_ndf"])
    pulls, info = [], []
    for a, c in PULL_WINDOWS:
        X = design(np.arange(a, c, 2.5) + 1.25, order)
        mw = np.exp(X @ b)
        pred, g = float(mw.sum()), (mw[:, None] * X).sum(0)
        sel = inr(m, a, c)
        obs, vmc, vfit = float(w[sel].sum() * k), float((w[sel] ** 2).sum() * k * k), float(g @ cov @ g)
        pulls.append((obs - pred) / math.sqrt(pred + vmc + vfit) if pred > 0 else float("nan"))
        info.append([obs, pred, vmc, vfit])
    out["A_pulls"], out["A_windows"] = pulls, info
    return out


def data_gof(dm, dp):
    """Gate B: passing data sidebands (5 GeV bins, 95-115 and 135-180), Poisson log-polynomial fit
    (F-test order), saturated-deviance GOF p-value."""
    y = np.array([(dp & inr(dm, a, a + 5)).sum() for a in SB_DATA_LO], float)
    xc = SB_DATA_LO + 2.5
    order, (b, cov, dev) = ftest_order(xc, y)
    ndf = len(y) - order - 1
    if b is None:
        return {"B_order": order, "B_p": float("nan"), "B_obs": y.tolist()}
    return {"B_order": order, "B_dev": dev, "B_ndf": ndf, "B_p": chi2_pvalue(dev, ndf), "B_obs": y.tolist(),
            "B_fit": np.exp(design(xc, order) @ b).tolist()}


def train(S, cfg, masses):
    import xgboost as xgb
    feats = cfg["feats_full"]
    Xs, ws = [], []
    for m in masses:
        s = S["train_sig_m%d" % m]
        x = matrix(s, feats, float(m)); w = np.clip(s["w_train"], 0, None)
        if (cfg.get("win_mA1") and m == 1) or m in cfg.get("win_masses", []):   # minscan3 copy: window training for any listed m_a
            keep = inr(s["H_m"], *WIN)
            x, w = x[keep], w[keep] * (w.sum() / w[keep].sum())
        Xs.append(x); ws.append(w)
    Xs = np.vstack(Xs); ws = np.concatenate(ws)
    b = S["train_bkg"]
    bhyp = np.random.default_rng(1).choice(masses, size=len(b)).astype(np.float64)
    Xb = matrix(b, feats, bhyp)
    wb = np.clip(b["w_train"], 0, None)
    if cfg.get("plane"):
        h, e = np.histogram(b["H_m"], np.arange(95, 181, 1.0), weights=wb)
        k = np.clip(np.digitize(b["H_m"], e) - 1, 0, len(h) - 1)
        wb = wb * np.where(h[k] > 0, h.mean() / h[k], 0.0)
    X = np.vstack([Xs, Xb]); y = np.concatenate([np.ones(len(Xs)), np.zeros(len(Xb))])
    w = np.concatenate([ws * (wb.sum() / ws.sum()), wb])
    clf = xgb.XGBClassifier(max_depth=4, n_estimators=int(cfg.get("ntrees", 300)), learning_rate=0.1, subsample=0.8,
                            min_child_weight=5, eval_metric="logloss", n_jobs=int(os.environ.get("NTHREADS", "1")),
                            verbosity=0)
    clf.fit(X, y, sample_weight=w)
    return clf


def evaluate(S, cfg, clf, m, n_target, rng):
    feats = cfg["feats_full"]
    d, dy, sg = S["data_sb"], S["dy"], S["test_sig_m%d" % m]
    if cfg.get("stored_prod"):
        sd, sb, ss = d["MVA_Score_mA_M%d" % m], dy["MVA_Score_mA_M%d" % m], sg["prod_score"]
    else:
        sd = clf.predict_proba(matrix(d, feats, float(m)))[:, 1]
        sb = clf.predict_proba(matrix(dy, feats, float(m)))[:, 1]
        ss = clf.predict_proba(matrix(sg, feats, float(m)))[:, 1]
    order = np.sort(sd[np.isfinite(sd)])[::-1]
    cut = float(0.5 * (order[n_target - 1] + order[n_target])) if len(order) > n_target else float(order[-1] - 1e-9)
    r = {"cut": cut, "n_data_sb": int((sd > cut).sum())}
    bm, bw = dy["H_m"], dy["factor"]
    ok = np.isfinite(bw) & np.isfinite(sb)
    bm, bw, sb = bm[ok], bw[ok], sb[ok]
    p = sb > cut; allm = np.ones_like(p)
    def ratios(ww):
        return [frac(ww, bm, lo, hi, p) / frac(ww, bm, lo, hi, allm) for lo, hi in ((120, 130), (105, 115), (135, 145))]
    r["R"], r["R_lo"], r["R_hi"] = ratios(bw)
    bt = np.array([ratios(bw * rng.poisson(1.0, len(bw))) for _ in range(NBOOT)])
    r["R_err"], r["R_lo_err"], r["R_hi_err"] = [float(np.nanstd(bt[:, i])) for i in range(3)]
    r["n_eff_120_130"] = float(bw[p & inr(bm, 120, 130)].sum() ** 2 / max((bw[p & inr(bm, 120, 130)] ** 2).sum(), 1e-30))
    r["dy_chi2"], r["dy_ndf"] = shape_chi2(bm, bw, p, EDGES)
    r["dy_p"] = chi2_pvalue(r["dy_chi2"], r["dy_ndf"])
    dm = d["H_m"]; dp = sd > cut
    n0 = np.array([inr(dm, lo, hi).sum() for lo, hi in BANDS], float)
    n1 = np.array([(dp & inr(dm, lo, hi)).sum() for lo, hi in BANDS], float)
    exp_ = n1.sum() * n0 / n0.sum()
    r["data_band_obs"] = n1.tolist(); r["data_band_exp"] = exp_.tolist()
    r["data_chi2"] = float(np.sum((n1 - exp_) ** 2 / np.where(exp_ > 0, exp_, 1))); r["data_ndf"] = len(BANDS) - 1
    r["data_p"] = chi2_pvalue(r["data_chi2"], r["data_ndf"])
    r.update(closure_proxy(bm, bw, p, r["n_data_sb"]))
    r.update(data_gof(dm, dp))
    sm, sw, g = sg["H_m"], sg["factor"] / 0.3, sg["gen_dR_gg"]
    ok = np.isfinite(ss) & np.isfinite(sw); sm, sw, g, ss = sm[ok], sw[ok], g[ok], ss[ok]
    ps = ss > cut
    S_all, S_win = sw[ps].sum(), sw[ps & inr(sm, *WIN)].sum()
    r.update(S_all=float(S_all), S_win=float(S_win), win_frac=float(S_win / S_all) if S_all > 0 else float("nan"),
             overlap=float(sw[ps & (g < 0.05)].sum() / S_all) if S_all > 0 else float("nan"),
             sec=float(sw[ps & inr(sm, 135, 160)].sum() / S_win) if S_win > 0 else float("nan"),
             presel_win_frac=float(sw[inr(sm, *WIN)].sum() / sw.sum()),
             presel_overlap=float(sw[g < 0.05].sum() / sw.sum()),
             presel_sec=float(sw[inr(sm, 135, 160)].sum() / sw[inr(sm, *WIN)].sum()))
    dyw, dys = bw[p & inr(bm, *WIN)].sum(), bw[p & (inr(bm, 95, 115) | inr(bm, 135, 180))].sum()
    B = r["n_data_sb"] * dyw / dys if dys > 0 else float("nan")
    r["B_win"] = float(B); r["Z"] = float(S_win / math.sqrt(B)) if B > 0 else float("nan")
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cfg = [json.loads(l) for l in open(a.config) if l.strip()][a.index]
    out = os.path.join(a.out, "%s.json" % cfg["id"])
    if os.path.exists(out):
        print("exists", out); return
    t0 = time.time()
    meta = json.load(open(os.path.join(a.samples, "meta.json")))
    masses = LOW if cfg["model"] == "low" else HIGH
    feats = list(cfg["feats"])
    if not cfg.get("stored_prod") and not any(f in feats for f in ("param", "param_rel")):
        feats = feats + ["param_rel"]
    cfg["feats_full"] = feats
    names = ["data_sb", "dy"] + ["test_sig_m%d" % m for m in masses]
    if not cfg.get("stored_prod"):
        names += ["train_bkg"] + ["train_sig_m%d" % m for m in masses]
    S = {n: Sample(os.path.join(a.samples, n + ".npz")) for n in names}
    clf = None if cfg.get("stored_prod") else train(S, cfg, masses)
    t1 = time.time()
    rng = np.random.default_rng(20261013)
    res = {"config": cfg, "per_mass": {}}
    for m in masses:
        res["per_mass"][str(m)] = evaluate(S, cfg, clf, m, int(meta["prod_data_sb_count"][str(m)]), rng)
    res["time_train_s"], res["time_eval_s"] = t1 - t0, time.time() - t1
    if clf is not None:
        imp = clf.get_booster().get_score(importance_type="gain")
        res["gain"] = {feats[int(k[1:])]: v for k, v in imp.items()}
    tmp = out + ".partial"
    json.dump(res, open(tmp, "w"), indent=1)
    os.replace(tmp, out)
    print("DONE", cfg["id"], "train %.0fs eval %.0fs" % (res["time_train_s"], res["time_eval_s"]))


if __name__ == "__main__":
    main()
