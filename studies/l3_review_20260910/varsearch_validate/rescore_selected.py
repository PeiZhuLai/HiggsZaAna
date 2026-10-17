#!/usr/bin/env python3
"""Sandbox re-scoring with the two BDTs of the variable search (models/model_{low,high}.pkl, train_selected.py):
replace MVA_Score_mA_M1..30 in a copy of run3_bdt_scored_fsrfix (m_a 1-3 low-mass model, 4-30 high-mass model, the
same routing as Parque2Root_BDT.add_mva_scores). Every other branch is copied unchanged.

Gate per file, before anything is written: the deployed models (HZaMVA/using/model_Za_BDT_{low,high}mass_run3.pkl)
are re-evaluated with the production features (HIGH: Parque2Root_BDT HIGH_FEATURE_COLUMNS, param = (ALP_mass - m)/H_mass;
LOW: the deployed meta.json features through search_job.matrix)
and must reproduce all 30 stored scores of every tree (same NaN pattern, max |diff| <= 1e-6), which proves the
columns read here are the ones the production scoring used.
New scores: features from search_job.column (the definitions of the search and of the training), NaN where a
feature is not finite or H_m is 0 (production rule). Output written as <out>.partial, renamed after a full read-back
(every branch of every tree: scores equal to the computed ones, all other branches identical to the source).
--src-base DIR reads <DIR>/<sample>/<era>.root instead of the production copy (same file staged locally).
--low-model NAME takes the m_a 1-3 model from models/model_NAME.pkl (features from model_NAME.meta.json, e.g. low_nodR).
Env: higgs-alp-ana. Usage: rescore_selected.py --sample <dir> --era <era> --out-base <dir> [--low-model NAME] [--verify-only]
"""
import argparse
import json
import os
import pickle
import sys

import awkward as ak
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "bdt_varsearch"))
import search_job as sj  # noqa: E402

SRC = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
USING = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using"
SEL = json.load(open(os.path.join(HERE, "selected.json")))
MASSES = range(1, 31)
LOW = {1, 2, 3}
SCORES = ["MVA_Score_mA_M%d" % m for m in MASSES]
TOL = 1e-6
PROD_HIGH = ["pho1Pt_oHm", "pho1R9", "pho1IetaIeta55", "pho1PIso_noCorr", "pho2Pt_oHm", "pho2R9", "pho2IetaIeta55",
             "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "var_dR_Za", "var_dR_g1g2", "var_dR_g1Z", "var_PtaOverMh",
             "H_pt_oHm", "pho_pt_asym", "param"]
# deployed low-mass features (low_swapPt_relE2 since the deployment; pho2_relEerr has no plain branch, so the low-mass
# gate builds its matrix with search_job.matrix, the definition the deployed m_a 1-3 scores were verified against)
PROD_LOW = json.load(open(os.path.join(USING, "model_Za_BDT_lowmass_run3.meta.json")))["features"]


def prod_matrix(a, cols, m):
    out = []
    for c in cols:
        if c == "param":
            out.append((a["ALP_mass"] - m) / a["H_mass"])
        elif c == "pho_pt_asym":
            out.append((a["pho1Pt_oHm"] - a["pho2Pt_oHm"]) / (a["pho1Pt_oHm"] + a["pho2Pt_oHm"] + 1e-6))
        else:
            out.append(a[c])
    return np.column_stack(out)


def score(clf, X, base):
    s = np.full(len(X), np.nan)
    ok = base & np.isfinite(X).all(axis=1)
    if ok.any():
        s[ok] = clf.predict_proba(X[ok])[:, 1]
    return s


def _eq(x, y):
    if isinstance(x, np.ndarray) and isinstance(y, np.ndarray):
        if x.dtype.kind == "f":
            return x.shape == y.shape and np.array_equal(x, y, equal_nan=True)
        return np.array_equal(x, y)
    fx, fy = ak.flatten(x, axis=None), ak.flatten(y, axis=None)
    return (ak.to_list(ak.num(x)) == ak.to_list(ak.num(y))
            and np.array_equal(np.asarray(fx), np.asarray(fy), equal_nan=True))


def full_compare(src, dst, expect=None):
    bad = []
    fs, fd = uproot.open(src), uproot.open(dst)
    ts = [k for k in fs.keys(cycle=False) if fs[k].classname == "TTree"]
    td = [k for k in fd.keys(cycle=False) if fd[k].classname == "TTree"]
    if sorted(ts) != sorted(td):
        return ["trees %s vs %s" % (ts, td)]
    for tn in ts:
        a, b = fs[tn], fd[tn]
        if a.num_entries != b.num_entries:
            bad.append("%s entries %d vs %d" % (tn, a.num_entries, b.num_entries)); continue
        if set(a.keys()) != set(b.keys()):
            bad.append("%s branch names differ" % tn); continue
        for br in a.keys():
            jag = "Jagged" in str(a[br].interpretation)
            try:
                y = b[br].array(library="ak" if jag else "np")
            except Exception as e:  # noqa: BLE001
                bad.append("%s:%s unreadable (%s)" % (tn, br, type(e).__name__)); continue
            if br in SCORES:
                if expect is not None and not np.array_equal(y, expect[(tn, br)], equal_nan=True):
                    bad.append("%s:%s differs from the computed score" % (tn, br))
                continue
            if not _eq(a[br].array(library="ak" if jag else "np"), y):
                bad.append("%s:%s differs from source" % (tn, br))
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--era", required=True)
    ap.add_argument("--out-base", required=True)
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--src-base", default=SRC, help="input base (default the production scored copy; condor jobs pass a local staged copy)")
    ap.add_argument("--high-model", default="high", help="models/model_<NAME>.pkl for m_a 4-30; features from its .meta.json")
    ap.add_argument("--low-model", default="low", help="models/model_<NAME>.pkl for m_a 1-3; features from its .meta.json")
    args = ap.parse_args()
    feats = {"low": json.load(open(os.path.join(HERE, "models", "model_%s.meta.json" % args.low_model)))["features"],
             "high": SEL["high"]["feats"] if args.high_model == "high" else
                     json.load(open(os.path.join(HERE, "models", "model_%s.meta.json" % args.high_model)))["features"]}
    src = os.path.join(args.src_base, args.sample, "%s.root" % args.era)
    dst = os.path.join(args.out_base, args.sample, "%s.root" % args.era)
    if args.verify_only:
        bad = full_compare(src, dst)
        print(("VERIFY BAD %s: " % dst + "; ".join(bad[:5])) if bad else "VERIFY OK %s" % dst, flush=True)
        sys.exit(1 if bad else 0)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    prod = {"low": pickle.load(open(os.path.join(USING, "model_Za_BDT_lowmass_run3.pkl"), "rb")),
            "high": pickle.load(open(os.path.join(USING, "model_Za_BDT_highmass_run3.pkl"), "rb"))}
    new = {k: pickle.load(open(os.path.join(HERE, "models", "model_%s.pkl" % n), "rb")) for k, n in (("low", args.low_model), ("high", args.high_model))}
    need = set(PROD_HIGH + ["ALP_mass", "H_mass", "H_m", "ALP_m", "pho1Pt", "pho2Pt"]) - {"param", "pho_pt_asym"}
    for k in ("low", "high", "prodlow"):
        for f in (PROD_LOW if k == "prodlow" else feats[k]):
            if f in ("pho1_relEerr", "pho2_relEerr"):
                w = "lead" if f.startswith("pho1") else "sublead"
                need |= {"ALP_%s_photon_energyErr" % w, "ALP_%s_photon_eta" % w}
            elif f not in ("param_rel", "param", "pho_pt_asym", "pho_pt_asym_oHm"):
                need.add(f)

    fin = uproot.open(src)
    trees = [k for k in fin.keys(cycle=False) if fin[k].classname == "TTree"]
    out_trees, expect = {}, {}
    for tn in trees:
        t = fin[tn]
        jag = [b.name for b in t.branches if "Jagged" in str(b.interpretation)]
        counters = {"n" + j for j in jag}
        names = [b.name for b in t.branches if b.name not in counters]
        cols = {n: t[n].array(library="ak" if n in jag else "np") for n in names}
        a = {k: np.asarray(cols[k], dtype=np.float64) for k in need}
        if not (np.array_equal(a["H_m"], a["H_mass"], equal_nan=True) and np.array_equal(a["ALP_m"], a["ALP_mass"], equal_nan=True)):
            sys.exit("GATE FAIL %s:%s H_m/ALP_m differ from H_mass/ALP_mass" % (src, tn))
        base = np.isfinite(a["H_mass"]) & (a["H_mass"] != 0)
        nnan_new = 0
        for m in MASSES:
            br = "MVA_Score_mA_M%d" % m
            k = "low" if m in LOW else "high"
            stored = np.asarray(cols[br], dtype=np.float64)
            rep = score(prod[k], sj.matrix(a, PROD_LOW, float(m)) if k == "low" else prod_matrix(a, PROD_HIGH, m), base)
            fs, fr = np.isfinite(stored), np.isfinite(rep)
            if not np.array_equal(fs, fr):
                sys.exit("GATE FAIL %s:%s %s: NaN pattern differs (%d stored vs %d recomputed)" % (src, tn, br, (~fs).sum(), (~fr).sum()))
            d = float(np.max(np.abs(stored[fs] - rep[fs]))) if fs.any() else 0.0
            if d > TOL:
                sys.exit("GATE FAIL %s:%s %s: production model does not reproduce stored scores (max diff %.3g)" % (src, tn, br, d))
            s_new = score(new[k], sj.matrix(a, feats[k], float(m)), base)
            nnan_new = max(nnan_new, int((~np.isfinite(s_new) & fs).sum()))
            cols[br] = s_new.astype(cols[br].dtype)
            expect[(tn, br)] = cols[br]
        out_trees[tn] = cols
        print("  %s/%s:%s entries=%d branches=%d repro OK (30 scores); new-score NaN where production finite: %d"
              % (args.sample, args.era, tn, t.num_entries, len(t.branches), nnan_new), flush=True)

    tmp = dst + ".partial"
    with uproot.recreate(tmp) as fo:
        for tn, cols in out_trees.items():
            fo[tn] = cols
    bad = full_compare(src, tmp, expect)
    if bad:
        sys.exit("READBACK FAIL %s: %s" % (tmp, "; ".join(bad[:5])))
    os.replace(tmp, dst)
    print("DONE %s" % dst, flush=True)


if __name__ == "__main__":
    main()
