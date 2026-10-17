#!/usr/bin/env python3
"""Sandbox re-scoring: replace the low-mass BDT scores MVA_Score_mA_M{1,2,3} with the win_mfree
model (../model_win_mfree.pkl) in a copy of run3_bdt_scored_fsrfix. Every other branch is copied
unchanged; the high-mass scores (mA4-30) are untouched.

Gate per file, before anything is written: the deployed low-mass model
(HZaMVA/using/model_Za_BDT_lowmass_run3.pkl) is re-evaluated with the production feature recipe and
must reproduce the stored MVA_Score_mA_M{1,2,3} of every tree (same NaN pattern, max |diff| <= 1e-6).
This proves that the feature construction used for the new model reads the same inputs as the
production scoring.

Features
  production : ALP_calculatedPhotonIso, var_dR_g1g2, pho1R9, pho1Pt_oHm, (ALP_mass - m_a)/H_mass
  win_mfree  : ALP_calculatedPhotonIso, var_dR_g1g2, pho1R9, pho1Pt_oHm*H_mass, (ALP_mass - m_a)/m_a

Output files are written as <out>.partial and renamed only after a read-back check (entries per tree,
branch set, new score branches equal to the computed ones).
Env: higgs-alp-ana. Usage: rescore_lowmass.py --sample <dir> --era <era> --out-base <dir>
"""
import argparse
import os
import pickle
import sys

import awkward as ak
import numpy as np
import uproot

SRC = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
PROD_MODEL = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_lowmass_run3.pkl"
NEW_MODEL = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/model_win_mfree.pkl"
LOWM = (1, 2, 3)
TOL = 1e-6


def feats_prod(a, m):
    return np.column_stack([a["ALP_calculatedPhotonIso"], a["var_dR_g1g2"], a["pho1R9"], a["pho1Pt_oHm"],
                            (a["ALP_mass"] - m) / a["H_mass"]])


def feats_new(a, m):
    return np.column_stack([a["ALP_calculatedPhotonIso"], a["var_dR_g1g2"], a["pho1R9"],
                            a["pho1Pt_oHm"] * a["H_mass"], (a["ALP_mass"] - m) / float(m)])


def predict(clf, X):
    if len(X) == 0:
        return np.zeros(0)
    return clf.predict_proba(X)[:, 1].astype(np.float64)


def _eq(x, y):
    if isinstance(x, np.ndarray) and isinstance(y, np.ndarray):
        if x.dtype.kind == "f":
            return x.shape == y.shape and np.array_equal(x, y, equal_nan=True)
        return np.array_equal(x, y)
    fx, fy = ak.flatten(x, axis=None), ak.flatten(y, axis=None)
    return (ak.to_list(ak.num(x)) == ak.to_list(ak.num(y))
            and np.array_equal(np.asarray(fx), np.asarray(fy), equal_nan=True))


def full_compare(src, dst, expect=None):
    """List of problems between src and dst (empty = identical apart from the low-mass scores)."""
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
            if br in ("MVA_Score_mA_M1", "MVA_Score_mA_M2", "MVA_Score_mA_M3"):
                if expect is not None and not np.array_equal(y, expect[(tn, br)], equal_nan=True):
                    bad.append("%s:%s differs from the computed score" % (tn, br))
                continue
            x = a[br].array(library="ak" if jag else "np")
            if not _eq(x, y):
                bad.append("%s:%s differs from source" % (tn, br))
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--era", required=True)
    ap.add_argument("--out-base", required=True)
    ap.add_argument("--verify-only", action="store_true",
                    help="only run the full comparison of an existing output against its source")
    args = ap.parse_args()
    src = os.path.join(SRC, args.sample, f"{args.era}.root")
    dst = os.path.join(args.out_base, args.sample, f"{args.era}.root")
    if args.verify_only:
        bad = full_compare(src, dst)
        print(("VERIFY BAD %s: " % dst + "; ".join(bad[:5])) if bad else "VERIFY OK %s" % dst, flush=True)
        sys.exit(1 if bad else 0)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    prod = pickle.load(open(PROD_MODEL, "rb"))
    new = pickle.load(open(NEW_MODEL, "rb"))

    fin = uproot.open(src)
    trees = [k for k in fin.keys(cycle=False) if fin[k].classname == "TTree"]
    out_trees, expect = {}, {}
    for tn in trees:
        t = fin[tn]
        jag = [b.name for b in t.branches if "Jagged" in str(b.interpretation)]
        counters = {"n" + j for j in jag}
        names = [b.name for b in t.branches if b.name not in counters]
        cols = {}
        for n in names:
            cols[n] = t[n].array(library="ak" if n in jag else "np")
        a = {k: np.asarray(cols[k], dtype=np.float64) for k in
             ("ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "ALP_mass", "H_mass")}
        for m in LOWM:
            br = f"MVA_Score_mA_M{m}"
            stored = np.asarray(cols[br], dtype=np.float64)
            rep = predict(prod, feats_prod(a, m))
            fs, fr = np.isfinite(stored), np.isfinite(rep)
            if not np.array_equal(fs, fr):
                sys.exit(f"GATE FAIL {src}:{tn} {br}: NaN pattern differs ({(~fs).sum()} stored vs {(~fr).sum()} recomputed)")
            d = float(np.max(np.abs(stored[fs] - rep[fs]))) if fs.any() else 0.0
            if d > TOL:
                sys.exit(f"GATE FAIL {src}:{tn} {br}: production model does not reproduce stored scores (max diff {d:.3g})")
            s_new = predict(new, feats_new(a, m))
            s_new[~fs] = np.nan
            cols[br] = s_new.astype(cols[br].dtype)
            expect[(tn, br)] = cols[br]
        out_trees[tn] = (cols, len(t.branches), t.num_entries)
        print(f"  {args.sample}/{args.era}:{tn} entries={t.num_entries} branches={len(t.branches)} "
              f"repro OK (mA1-3)", flush=True)

    tmp = dst + ".partial"
    with uproot.recreate(tmp) as fo:
        for tn, (cols, nb, ne) in out_trees.items():
            fo[tn] = cols
    # read-back check: every branch of every tree is read in full; the low-mass scores must equal the
    # computed ones and every other branch must be bit-identical to the source. (2026-10-13: a partial
    # check -- scores + two reference branches -- passed a file whose Z_noFSR_lead_lepton_eta basket
    # had been written corrupt to EOS.)
    bad = full_compare(src, tmp, expect)
    if bad:
        sys.exit("READBACK FAIL %s: %s" % (tmp, "; ".join(bad[:5])))
    os.replace(tmp, dst)
    print(f"DONE {dst}", flush=True)


if __name__ == "__main__":
    main()
