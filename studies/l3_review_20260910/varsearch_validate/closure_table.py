#!/usr/bin/env python3
"""Table of the step-3 closure scans (closure_scan/{new,prod}_drop_mA<m>/summary.txt) with a relative sensitivity per
cut: S = signal 'test' weight (factor) above the cut in 120 < m_llgg < 130 (eval_samples64 cache, new scores from
models/model_<low|high>.pkl, production scores stored), N = data count of the closure data workspace (dataN),
Zrel = S / sqrt(N) normalized to the production working point. Gates printed per cut: closure Z < 1, every data
envelope member GOF > 0.01.
Env: higgs-alp-ana. Output: eval/closure_table.txt
"""
import glob
import os
import pickle
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import eval_selected as es  # noqa: E402
sj = es.sj
B = "/eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_scan"
CACHE = "/eos/home-p/pelai/HZa/sandbox_varsearch_v1/eval_samples64"


def parse(path):
    rows = []
    for l in open(path):
        if "dataN=" not in l:
            continue
        g = lambda k: re.search(k + r"=(\S+)", l).group(1)
        dg = [float(x.split(":")[1]) for x in l.split("data GOF:")[1].split("|")[0].split()]
        rows.append({"cut": float(l.split()[0]), "N": int(g("dataN")), "Z": float(g("closureZ")), "rhat": l.split("rhat=")[1].split("data_fallback")[0].strip(),
                     "gof_min": min(dg), "fallback": int(g("data_fallback"))})
    return rows


def main():
    models = {k: pickle.load(open(os.path.join(HERE, "models", "model_%s.pkl" % k), "rb")) for k in ("low", "high")}
    lines = []
    for m in (1, 2, 3, 8, 15, 30):
        k = "low" if m <= 3 else "high"
        a = dict(np.load(os.path.join(CACHE, "test_sig_m%d.npz" % m)))
        win = (a["H_m"] > 120) & (a["H_m"] < 130)
        s_new = np.full(len(a["H_m"]), np.nan)
        X = sj.matrix(a, es.SEL[k]["feats"], float(m)); ok = np.isfinite(X).all(axis=1)
        s_new[ok] = models[k].predict_proba(X[ok])[:, 1]
        prod = parse(os.path.join(B, "prod_drop_mA%d" % m, "summary.txt"))[0]
        sp = float(a["factor"][win & (a["prod_score"] > prod["cut"])].sum())
        zp = sp / np.sqrt(prod["N"])
        lines.append("mA%-2d %-9s %7s %5s %7s %-24s %8s %5s  %s" % (m, "model", "cut", "N", "Z", "rhat", "GOFmin", "Zrel", "gates"))
        for who, r in [("prod", prod)] + [("new", r) for r in parse(os.path.join(B, "new_drop_mA%d" % m, "summary.txt"))]:
            sc = a["prod_score"] if who == "prod" else s_new
            s = float(a["factor"][win & (sc > r["cut"])].sum())
            ok = r["Z"] < 1.0 and r["gof_min"] > 0.01
            lines.append("     %-9s %7.5f %5d %7.3f %-24s %8.3f %5.2f  %s%s" % (who, r["cut"], r["N"], r["Z"], r["rhat"], r["gof_min"],
                         s / np.sqrt(r["N"]) / zp, "PASS" if ok else "fail", " (data fallback)" if r["fallback"] else ""))
    txt = "\n".join(lines)
    open(os.path.join(HERE, "eval", "closure_table.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
