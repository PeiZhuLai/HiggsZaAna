#!/usr/bin/env python3
"""Working points of the variable-search BDTs for the step-4 full chain, from the step-3 closure scans
(closure_scan/new_drop_mA<m>/summary.txt, 5 cuts = data-sideband count x 2, 1.5, 1, 0.7, 0.5 of production).
Rule, fixed before the mA4-7, 9, 10, 20, 25 scans were looked at:
  * mA1: the x0.5 cut (user decision 2026-10-14, option 1 after the nodR comparison);
  * every other anchor: among the cuts with closure Z < 1 and every data-envelope GOF > 0.01, the one with the
    largest relative sensitivity Zrel = S/sqrt(N) (S = test-signal weight above the cut in 120-130 GeV, N = data
    count of the closure data workspace; same definition as closure_table.py). No passing cut -> exit 1.
Interpolated masses take the nearest anchor's cut downstream (apply_bdt_data.find_nearest_cut), as in production.
Env: higgs-alp-ana. Output: chain4/plot_output/MVAcut_points_run3_chain4.json, eval/wp_chain4.txt
"""
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import closure_table as ct  # noqa: E402
es, sj = ct.es, ct.sj
ANCHORS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30)
OUT = os.path.join(HERE, "chain4", "plot_output", "MVAcut_points_run3_chain4.json")


def main():
    models = {k: pickle.load(open(os.path.join(HERE, "models", "model_%s.pkl" % k), "rb")) for k in ("low", "high")}
    res, lines, bad = [], [], []
    for m in ANCHORS:
        k = "low" if m <= 3 else "high"
        a = dict(np.load(os.path.join(ct.CACHE, "test_sig_m%d.npz" % m)))
        win = (a["H_m"] > 120) & (a["H_m"] < 130)
        s_new = np.full(len(a["H_m"]), np.nan)
        X = sj.matrix(a, es.SEL[k]["feats"], float(m)); ok = np.isfinite(X).all(axis=1)
        s_new[ok] = models[k].predict_proba(X[ok])[:, 1]
        rows = ct.parse(os.path.join(ct.B, "new_drop_mA%d" % m, "summary.txt"))
        if len(rows) != 5:
            bad.append("mA%d: %d scan rows" % (m, len(rows))); continue
        for r in rows:
            r["Zrel"] = float(a["factor"][win & (s_new > r["cut"])].sum()) / np.sqrt(r["N"])
            r["pass"] = r["Z"] < 1.0 and r["gof_min"] > 0.01
        if m == 1:
            pick, why = rows[4], "x0.5 (user decision 2026-10-14)"
        else:
            cand = [r for r in rows if r["pass"]]
            if not cand:
                bad.append("mA%d: no cut passes closure Z < 1 and data GOF > 0.01" % m); continue
            pick = max(cand, key=lambda r: r["Zrel"])
            why = "max Zrel among %d passing cuts (x%s)" % (len(cand), ("2", "1.5", "1", "0.7", "0.5")[rows.index(pick)])
        res.append({"mA": m, "MVAcut": pick["cut"], "note": "varsearch step 4: " + why})
        lines.append("mA%-2d cut %.5f  N %4d  closureZ %.3f  GOFmin %.3f  %s" % (m, pick["cut"], pick["N"], pick["Z"], pick["gof_min"], why))
    print("\n".join(lines))
    if bad:
        print("\n".join(bad)); sys.exit(1)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump({"results": res}, open(OUT, "w"), indent=1)
    open(os.path.join(HERE, "eval", "wp_chain4.txt"), "w").write("\n".join(lines) + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
