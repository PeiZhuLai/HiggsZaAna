#!/usr/bin/env python3
"""Step 1 of the validation for a minimal-change low-mass candidate (user request 2026-10-14; candidates from
minscan2): float64 retrain with train_selected.build (same events, weights, recipe; m_a = 1 signal on 115-133 GeV),
then the search metrics and gates on float64 evaluation sets with the DY of the search ("incl") and the DY without
the BDT training events ("test"), for the candidate and for the stored production low-mass scores, exactly as
eval_selected.py does. The evaluation sets are rebuilt privately per candidate (the shared cache lacks some of the
production low-mass columns), only for m_a = 1, 2, 3.
Runs on condor (packed higgs-alp-ana env) or locally. Production code and files are only read.
Output (EOS): <out>/<name>/{model.pkl, model.meta.json, eval_<new|prod>_<incl|test>.json, summary.txt, cache/}
Usage: validate_lowcand.py <name> <out> <feature> [<feature> ...]
"""
import json
import os
import pickle
import sys
import time

import numpy as np
import xgboost as xgb

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import train_selected as ts  # noqa: E402
import eval_selected as es  # noqa: E402
sj = ts.sj
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "bdt_varsearch"))
import driver_condor as dc  # noqa: E402


def main():
    name, out, feats = sys.argv[1], sys.argv[2], sys.argv[3:]
    od = os.path.join(out, name)
    os.makedirs(od, exist_ok=True)
    cfg = {"feats": feats, "masses": [1, 2, 3], "win_mA1": True, "source": "minscan2 (%s)" % name}
    ts.SEL = {name: cfg}
    es.SEL = {name: cfg}
    es.ANCHORS = [1, 2, 3]
    es.SCORES = ["MVA_Score_mA_M%d" % m for m in es.ANCHORS]

    def raw_columns():   # eval_selected.raw_columns would read "param" as a branch
        cols = {"H_m", "ALP_m", "factor", "event", "pho1Pt", "pho2Pt"}
        for f in feats:
            if f in ("param", "param_rel", "pho_pt_asym", "pho_pt_asym_oHm"):
                continue
            if f in ("pho1_relEerr", "pho2_relEerr"):
                wh = "lead" if f.startswith("pho1") else "sublead"
                cols |= {"ALP_%s_photon_energyErr" % wh, "ALP_%s_photon_eta" % wh}
            else:
                cols.add(f)
        return sorted(cols)
    es.raw_columns = raw_columns
    t0 = time.time()
    print("===== %s: %s" % (name, feats), flush=True)
    X, y, w = ts.build(name)
    clf = xgb.XGBClassifier(max_depth=4, n_estimators=300, learning_rate=0.1, subsample=0.8, min_child_weight=5,
                            eval_metric="logloss", n_jobs=1, verbosity=0)
    clf.fit(X, y, sample_weight=w)
    pickle.dump(clf, open(os.path.join(od, "model.pkl.partial"), "wb"))
    os.replace(os.path.join(od, "model.pkl.partial"), os.path.join(od, "model.pkl"))
    json.dump({"model": name, "features": feats, "masses": [1, 2, 3], "mA1_window": sj.WIN, "inputs": ts.INP,
               "dtype": "float64", "xgboost": xgb.__version__,
               "recipe": "XGB depth4/300/lr0.1/subsample0.8/mcw5 (search_job.train = train_*_final.py)"},
              open(os.path.join(od, "model.meta.json"), "w"), indent=1)
    print("  [train] %d events, %.0f s" % (len(y), time.time() - t0), flush=True)
    cdir = os.path.join(od, "cache")
    es.build_cache(cdir)
    S = es.load_cache(cdir)
    meta = json.load(open(os.path.join(es.SLIM, "meta.json")))
    n_target = meta["prod_data_sb_count"]
    for m in (1, 2, 3):
        n = int((S["data_sb"]["MVA_Score_mA_M%d" % m] > meta["production_wp"][str(m)]).sum())
        assert n == int(n_target[str(m)]), (m, n, n_target[str(m)])
    res = {}
    for who, c in (("prod", None), ("new", clf)):
        for v in ("incl", "test"):
            res[(who, v)] = es.run(S, name, c, v, n_target)
            json.dump(res[(who, v)], open(os.path.join(od, "eval_%s_%s.json" % (who, v)), "w"), indent=1)
    zref = {m: dc.z_fit(r) for m, r in res[("prod", "incl")]["per_mass"].items()}
    L = ["%s: %s (m_a = 1 signal trained on %s GeV)" % (name, ", ".join(feats), sj.WIN)]
    for who in ("new", "prod"):
        for v in ("incl", "test"):
            r = res[(who, v)]
            L.append("%-4s %-4s violation %6.2f  objective %.3f" % (who, v, dc.gate_violation(r), dc.objective(r, zref)))
    L.append(" mA | Z/Zprod incl | R test new   prod        | R_lo test new | A pulls test (new)  | B p    | mA1 win/ovl test (prod, presel)")
    for m in ("1", "2", "3"):
        t, p, i = res[("new", "test")]["per_mass"][m], res[("prod", "test")]["per_mass"][m], res[("new", "incl")]["per_mass"][m]
        L.append("%3s | %6.2f | %5.2f+-%.2f %5.2f+-%.2f | %5.2f+-%.2f | %5.1f %5.1f %5.1f | %.4f | %.2f/%.2f (%.2f/%.2f, %.2f/%.2f)"
                 % (m, dc.z_fit(i) / zref[m], t["R"], t["R_err"], p["R"], p["R_err"], t["R_lo"], t["R_lo_err"], *t["A_pulls"], t["B_p"],
                    t["win_frac"], t["overlap"], p["win_frac"], p["overlap"], t["presel_win_frac"], t["presel_overlap"]))
    txt = "\n".join(L)
    open(os.path.join(od, "summary.txt"), "w").write(txt + "\n")
    print(txt, flush=True)
    print("DONE %s %.0f s" % (name, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
