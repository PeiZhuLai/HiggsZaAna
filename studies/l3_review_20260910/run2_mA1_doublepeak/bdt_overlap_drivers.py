#!/usr/bin/env python3
"""Why does the low-mass BDT working point enrich m_a = 1 GeV overlapping-shower events?

Signal: run3_bdt_scored_fsrfix/mA_M1/<era>.root (inclusive = preselection), weight = factor,
split by generator dR(gamma, gamma) < 0.05 ("overlap") or >= 0.05 ("normal").
Background: run3_bdt_inputs_fsrfix/All_Bkg/run3.root (inclusive), param at m_a,hyp = 1 GeV.

1. The deployed low-mass model is re-evaluated and checked against the stored MVA_Score_mA_M1.
2. Pass rate of the m_a = 1 working point (score > 0.963) for overlap / normal signal events.
3. Weighted quantiles of the five inputs for overlap, normal and background.
4. Ablation: in the overlap events, one input at a time is replaced by values drawn from the
   normal events (same weights), and the pass rate is recomputed. The input whose replacement
   lowers the pass rate most is the one through which the BDT prefers the overlap events.

Env: higgs-alp-ana (xgboost 2.1.1, the env the model was trained in).
"""

import glob
import pickle

import numpy as np
import uproot

SIG = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M1"
BKG = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix/All_Bkg/run3.root"
MODEL = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_lowmass_run3.pkl"
FEATS = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]
RAW = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "ALP_m", "H_m", "factor"]
WP, MHYP = 0.963, 1.0
rng = np.random.default_rng(1)


def load(paths, extra=()):
    out = {}
    for p in paths:
        a = uproot.open(p)["inclusive"].arrays(RAW + list(extra), library="np")
        for k, v in a.items():
            out.setdefault(k, []).append(v)
    a = {k: np.concatenate(v) for k, v in out.items()}
    a["param"] = (a["ALP_m"] - MHYP) / a["H_m"]
    sel = (a["H_m"] > 95) & (a["H_m"] < 180)
    return {k: v[sel] for k, v in a.items()}


def wquant(x, w, qs=(0.16, 0.5, 0.84)):
    o = np.argsort(x); c = np.cumsum(w[o]); c /= c[-1]
    return [float(np.interp(q, c, x[o])) for q in qs]


def main():
    clf = pickle.load(open(MODEL, "rb"))
    s = load(sorted(glob.glob(SIG + "/*.root")), extra=("GenALP_dR_gg", "MVA_Score_mA_M1"))
    b = load([BKG])
    X = np.column_stack([s[f] for f in FEATS])
    score = clf.predict_proba(X)[:, 1]
    print("model check: max |score - stored| = %.2e over %d events"
          % (np.abs(score - s["MVA_Score_mA_M1"]).max(), len(score)))

    w = s["factor"]
    ov = s["GenALP_dR_gg"] < 0.05
    no = ~ov
    print("\nshare of overlap events: preselection %.3f, after WP %.3f" % (
        w[ov].sum() / w.sum(), w[ov & (score > WP)].sum() / w[score > WP].sum()))
    eff_ov = w[ov & (score > WP)].sum() / w[ov].sum()
    eff_no = w[no & (score > WP)].sum() / w[no].sum()
    print("WP pass rate: overlap %.3f   normal %.3f   ratio %.2f" % (eff_ov, eff_no, eff_ov / eff_no))

    wb = np.clip(b["factor"], 0, None)
    print("\nweighted 16/50/84%% quantiles          overlap                    normal                     background")
    for f in FEATS:
        q = [wquant(s[f][ov], w[ov]), wquant(s[f][no], w[no]), wquant(b[f], wb)]
        print("%-24s" % f + "".join("  [%7.3f %7.3f %7.3f]" % tuple(x) for x in q))

    # ablation: overlap events with one input replaced by draws from the normal events
    print("\nablation (overlap events, input replaced by values of normal events): WP pass rate")
    pn = w[no] / w[no].sum()
    Xov = X[ov]
    print("  %-26s %.3f" % ("none", eff_ov))
    for i, f in enumerate(FEATS):
        r = []
        for rep in range(5):
            Xr = Xov.copy()
            Xr[:, i] = X[no][rng.choice(no.sum(), size=ov.sum(), p=pn), i]
            sr = clf.predict_proba(Xr)[:, 1]
            r.append(w[ov][sr > WP].sum() / w[ov].sum())
        print("  %-26s %.3f +- %.3f" % (f, np.mean(r), np.std(r)))
    # and the reverse: normal events with one input taken from the overlap events
    print("\nreverse (normal events, input replaced by values of overlap events): WP pass rate")
    po = w[ov] / w[ov].sum()
    Xno = X[no]
    print("  %-26s %.3f" % ("none", eff_no))
    for i, f in enumerate(FEATS):
        Xr = Xno.copy()
        Xr[:, i] = X[ov][rng.choice(ov.sum(), size=no.sum(), p=po), i]
        sr = clf.predict_proba(Xr)[:, 1]
        print("  %-26s %.3f" % (f, w[no][sr > WP].sum() / w[no].sum()))


if __name__ == "__main__":
    main()
