#!/usr/bin/env python3
"""Replace MVA_Score_mA_M1..3 of one production scored file (run3_bdt_scored_fsrfix) with the low_swapPt_relE2
low-mass BDT; every other branch, including MVA_Score_mA_M4..30, is copied unchanged.

Recipe = the updated Parque2Root_BDT.add_mva_scores for m_a = 1, 2, 3:
  features LOW_FEATURE_COLUMNS = ALP_calculatedPhotonIso, var_dR_g1g2, pho1R9, pho2_relEerr, param
  pho2_relEerr = ALP_sublead_photon_energyErr / (ALP_sublead_photon_pt cosh ALP_sublead_photon_eta) in float64
  param = (ALP_mass - m_a) / H_mass; NaN where H_mass is 0 / not finite or any feature is not finite.
Gates, all before the output is accepted:
  1. the previous low-mass model with the previous features reproduces the stored M1..3 scores (same NaN pattern,
     max |diff| <= 1e-6): the columns read here are the ones the production scoring used;
  2. with --ref, the new M1..3 equal the scores of the sandbox chain copy (run3_bdt_scored_low_swapPt_relE2) that
     produced the validated chain results (same NaN pattern, max |diff| <= 1e-6);
  3. full read-back of the output: M1..3 equal the computed scores, every other branch identical to the source.
Env: higgs-alp-ana. Usage: patch_scores.py --src <in.root> --out <out.root> --old-low <pkl> --new-low <pkl> [--ref <root>]
"""
import argparse
import os
import pickle
import sys

import awkward as ak
import numpy as np
import uproot

LOW = (1, 2, 3)
OLD_FEATS = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]
NEW_FEATS = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho2_relEerr", "param"]
NEED = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "ALP_mass", "H_mass",
        "ALP_sublead_photon_energyErr", "ALP_sublead_photon_pt", "ALP_sublead_photon_eta"]
TOL = 1e-6


def matrix(a, feats, m):
    cols = []
    for f in feats:
        if f == "param":
            cols.append((a["ALP_mass"] - m) / a["H_mass"])
        elif f == "pho2_relEerr":
            cols.append(a["ALP_sublead_photon_energyErr"] / (a["ALP_sublead_photon_pt"] * np.cosh(a["ALP_sublead_photon_eta"])))
        else:
            cols.append(a[f])
    return np.column_stack(cols)


def score(clf, X, base):
    s = np.full(len(X), np.nan)
    ok = base & np.isfinite(X).all(axis=1)
    if ok.any():
        s[ok] = clf.predict_proba(X[ok])[:, 1]
    return s


def same(x, y):
    fx, fy = np.isfinite(x), np.isfinite(y)
    if not np.array_equal(fx, fy):
        return False, "NaN pattern %d vs %d" % ((~fx).sum(), (~fy).sum())
    d = float(np.max(np.abs(x[fx] - y[fy]))) if fx.any() else 0.0
    return d <= TOL, "max diff %.3g" % d


def eq(x, y):
    if isinstance(x, np.ndarray) and isinstance(y, np.ndarray):
        if x.dtype.kind == "f":
            return x.dtype == y.dtype and x.shape == y.shape and np.array_equal(x, y, equal_nan=True)
        return x.dtype == y.dtype and np.array_equal(x, y)
    fx, fy = ak.flatten(x, axis=None), ak.flatten(y, axis=None)
    return ak.to_list(ak.num(x)) == ak.to_list(ak.num(y)) and np.array_equal(np.asarray(fx), np.asarray(fy), equal_nan=True)


def trees_of(f):
    return [k for k in f.keys(cycle=False) if f[k].classname == "TTree"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--old-low", required=True)
    ap.add_argument("--new-low", required=True)
    ap.add_argument("--ref", default=None)
    a = ap.parse_args()
    old = pickle.load(open(a.old_low, "rb"))
    new = pickle.load(open(a.new_low, "rb"))
    assert old.n_features_in_ == len(OLD_FEATS) and new.n_features_in_ == len(NEW_FEATS)
    fin = uproot.open(a.src)
    ref = uproot.open(a.ref) if a.ref else None
    out_trees, expect = {}, {}
    for tn in trees_of(fin):
        t = fin[tn]
        jag = [b.name for b in t.branches if "Jagged" in str(b.interpretation)]
        counters = {"n" + j for j in jag}
        names = [b.name for b in t.branches if b.name not in counters]
        cols = {n: t[n].array(library="ak" if n in jag else "np") for n in names}
        v = {k: np.asarray(cols[k], dtype=np.float64) for k in NEED}
        base = np.isfinite(v["H_mass"]) & (v["H_mass"] != 0)
        msg = []
        for m in LOW:
            br = "MVA_Score_mA_M%d" % m
            stored = np.asarray(cols[br], dtype=np.float64)
            ok, why = same(stored, score(old, matrix(v, OLD_FEATS, m), base))
            if not ok:
                sys.exit("GATE1 FAIL %s:%s %s: previous model does not reproduce the stored score (%s)" % (a.src, tn, br, why))
            s_new = score(new, matrix(v, NEW_FEATS, m), base)
            if ref is not None:
                r = np.asarray(ref[tn][br].array(library="np"), dtype=np.float64)
                ok, why = same(s_new, r)
                if not ok:
                    sys.exit("GATE2 FAIL %s:%s %s: differs from the sandbox chain score (%s)" % (a.src, tn, br, why))
            cols[br] = s_new.astype(cols[br].dtype)
            expect[(tn, br)] = cols[br]
            msg.append("M%d nan %d->%d" % (m, int((~np.isfinite(stored)).sum()), int((~np.isfinite(s_new)).sum())))
        out_trees[tn] = cols
        print("  %s:%s entries=%d branches=%d gate1 OK%s; %s" % (a.src, tn, t.num_entries, len(t.branches),
              ", gate2 OK" if ref is not None else "", ", ".join(msg)), flush=True)
    tmp = a.out + ".partial"
    with uproot.recreate(tmp) as fo:
        for tn, cols in out_trees.items():
            fo[tn] = cols
    fs, fd = uproot.open(a.src), uproot.open(tmp)
    bad = []
    if sorted(trees_of(fs)) != sorted(trees_of(fd)):
        bad.append("trees differ")
    for tn in trees_of(fs):
        x, y = fs[tn], fd[tn]
        if x.num_entries != y.num_entries or set(x.keys()) != set(y.keys()):
            bad.append("%s entries/branch names differ" % tn); continue
        for br in x.keys():
            jag = "Jagged" in str(x[br].interpretation)
            yv = y[br].array(library="ak" if jag else "np")
            if (tn, br) in expect:
                if not np.array_equal(yv, expect[(tn, br)], equal_nan=True):
                    bad.append("%s:%s differs from the computed score" % (tn, br))
            elif not eq(x[br].array(library="ak" if jag else "np"), yv):
                bad.append("%s:%s differs from the source" % (tn, br))
    if bad:
        sys.exit("READBACK FAIL %s: %s" % (tmp, "; ".join(bad[:5])))
    os.replace(tmp, a.out)
    print("DONE %s" % a.out, flush=True)


if __name__ == "__main__":
    main()
