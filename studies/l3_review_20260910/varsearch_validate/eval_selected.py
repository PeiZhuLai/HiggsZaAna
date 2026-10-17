#!/usr/bin/env python3
"""Evaluate the float64 retrained models (train_selected.py) with the metrics and gates of the variable search, on
float64 evaluation sets, with two DY versions:
  incl  the DY set of the search (run3_bdt_scored_fsrfix 'inclusive'): reproduces the search numbers
  test  the same events minus those of the run3_bdt_inputs_fsrfix 'train' trees. All_Bkg 'train' (the background
        training pool) is exactly DYJetsToLL + DYGto2LG_10to100 'train', so after this removal no DY event of the
        evaluation was available to the training. In the inputs, DYJetsToLL/2024 is the merge of 2024
        DYJetsTo2E/2Mu/2Tau (same 188313 events), so those three are tagged against DYJetsToLL/2024 'train'.
        Each sample x era is rescaled by sum w(inclusive) / sum w(kept), so the sample mix is unchanged.
The production model is evaluated from the stored scores on the same sets (reference Z and gates).
Signal: 'test' trees (30 %); data: sidebands only (95-115, 135-180), never the window.

Env: higgs-alp-ana. Output: eval/<model>_<new|prod>_<incl|test>.json, eval/summary.txt; the float64 evaluation
sets are cached in <cache>/*.npz.
Usage: eval_selected.py [--cache DIR]
"""
import argparse
import json
import os
import pickle
import sys

import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
SEARCH = os.path.join(os.path.dirname(HERE), "bdt_varsearch")
sys.path.insert(0, SEARCH)
import search_job as sj  # noqa: E402
import driver_condor as dc  # noqa: E402

INP = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
SCO = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
SLIM = "/eos/home-p/pelai/HZa/bdt_varsearch_20261013/samples"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
DY = {e: ["DYGto2LG_10to100", "DYJetsToLL"] for e in ERAS[:4]}
DY["2024"] = ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
TRAIN_SRC = {"DYGto2LG_10to100": "DYGto2LG_10to100", "DYJetsToLL": "DYJetsToLL",
             "DYJetsTo2E": "DYJetsToLL", "DYJetsTo2Mu": "DYJetsToLL", "DYJetsTo2Tau": "DYJetsToLL"}
ANCHORS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
SEL = json.load(open(os.path.join(HERE, "selected.json")))
SCORES = ["MVA_Score_mA_M%d" % m for m in ANCHORS]


def raw_columns():
    cols = {"H_m", "ALP_m", "factor", "event", "pho1Pt", "pho2Pt"}
    for cfg in SEL.values():
        for f in cfg["feats"]:
            if f in ("param_rel",):
                continue
            if f in ("pho1_relEerr", "pho2_relEerr"):
                w = "lead" if f.startswith("pho1") else "sublead"
                cols |= {"ALP_%s_photon_energyErr" % w, "ALP_%s_photon_eta" % w}
            else:
                cols.add(f)
    return sorted(cols)


def mask(a):
    return (a["H_m"] > 95) & (a["H_m"] < 180) & np.isfinite(a["factor"])


def key(a):
    """Per-event key for matching the same event between the input and scored copies."""
    return a["event"].astype(np.uint64) ^ (a["H_m"].astype(np.float32).view(np.uint32).astype(np.uint64) << np.uint64(32))


def build_cache(cdir):
    os.makedirs(cdir, exist_ok=True)
    cols = raw_columns()
    def save(name, d):
        np.savez(os.path.join(cdir, name + ".partial.npz"), **d)
        os.replace(os.path.join(cdir, name + ".partial.npz"), os.path.join(cdir, name + ".npz"))
        print("  [cache] %-12s %d events" % (name, len(d["H_m"])), flush=True)
    for m in ANCHORS:
        if os.path.exists(os.path.join(cdir, "test_sig_m%d.npz" % m)):
            continue
        acc = {}
        for e in ERAS:
            t = uproot.open("%s/mA_M%d/%s.root" % (SCO, m, e))["test"]
            a = t.arrays(cols + ["GenALP_dR_gg", "MVA_Score_mA_M%d" % m], library="np")
            g = np.array([x[0] if len(x) else np.nan for x in a.pop("GenALP_dR_gg")], float)
            k = mask(a)
            a = {c: v[k].astype(np.float64) for c, v in a.items()}
            a["gen_dR_gg"] = g[k]; a["prod_score"] = a.pop("MVA_Score_mA_M%d" % m)
            for c, v in a.items():
                acc.setdefault(c, []).append(v)
        save("test_sig_m%d" % m, {c: np.concatenate(v) for c, v in acc.items()})
    if not os.path.exists(os.path.join(cdir, "data_sb.npz")):
        acc = {}
        for e in ERAS:
            a = uproot.open("%s/Data/%s.root" % (SCO, e))["inclusive"].arrays(cols + SCORES, library="np")
            h = a["H_m"]; k = mask(a) & (((h > 95) & (h < 115)) | ((h > 135) & (h < 180)))
            for c, v in a.items():
                acc.setdefault(c, []).append(v[k].astype(np.float64))
        save("data_sb", {c: np.concatenate(v) for c, v in acc.items()})
    if not os.path.exists(os.path.join(cdir, "dy.npz")):
        acc = {}
        for e in ERAS:
            tagged, ktrs = {}, {}
            for s in DY[e]:
                a = uproot.open("%s/%s/%s.root" % (SCO, s, e))["inclusive"].arrays(cols + SCORES, library="np")
                k = mask(a)
                a = {c: v[k].astype(np.float64) for c, v in a.items()}
                src = TRAIN_SRC[s]
                if src not in ktrs:
                    tr = uproot.open("%s/%s/%s.root" % (INP, src, e))["train"].arrays(["event", "H_m", "factor"], library="np")
                    ktrs[src] = key({c: v[mask(tr)] for c, v in tr.items()})
                intrain = np.isin(key(a), ktrs[src])
                tagged[src] = tagged.get(src, 0) + int(intrain.sum())
                print("  %-18s %-13s inclusive %7d, tagged as train %7d (against %s/%s train)"
                      % (s, e, len(intrain), intrain.sum(), src, e), flush=True)
                w, kept = a["factor"], ~intrain
                a["w_test"] = np.where(kept, w * (w.sum() / w[kept].sum()), 0.0)
                a["in_train"] = intrain.astype(np.float64)
                for c, v in a.items():
                    acc.setdefault(c, []).append(v)
            for src, ktr in ktrs.items():
                print("  %s %s: tagged %d of %d train events" % (src, e, tagged[src], len(ktr)), flush=True)
                if tagged[src] != len(ktr):
                    sys.exit("train tagging incomplete for %s %s" % (src, e))
        save("dy", {c: np.concatenate(v) for c, v in acc.items()})


def load_cache(cdir):
    return {f[:-4]: dict(np.load(os.path.join(cdir, f))) for f in os.listdir(cdir) if f.endswith(".npz") and "partial" not in f}


def run(S, model, clf, variant, n_target):
    feats = SEL[model]["feats"]
    cfg = {"feats_full": feats, "stored_prod": clf is None}
    dy = dict(S["dy"])
    if variant == "test":
        keep = dy["in_train"] == 0
        dy = {c: v[keep] for c, v in dy.items()}
        dy["factor"] = dy["w_test"]
    SS = dict(S); SS["dy"] = dy
    rng = np.random.default_rng(20261013)
    res = {"config": {"model": model, "feats": feats, "variant": variant, "stored_prod": clf is None}, "per_mass": {}}
    for m in SEL[model]["masses"]:
        res["per_mass"][str(m)] = sj.evaluate(SS, cfg, clf, m, n_target[str(m)], rng)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default="/eos/home-p/pelai/HZa/sandbox_varsearch_v1/eval_samples64")
    a = ap.parse_args()
    build_cache(a.cache)
    S = load_cache(a.cache)
    n_target = json.load(open(os.path.join(SLIM, "meta.json")))["prod_data_sb_count"]
    for m in ANCHORS:   # the float64 data sideband reproduces the production count used by the search
        n = int((S["data_sb"]["MVA_Score_mA_M%d" % m] > json.load(open(os.path.join(SLIM, "meta.json")))["production_wp"][str(m)]).sum())
        assert n == int(n_target[str(m)]), (m, n, n_target[str(m)])
    lines = []
    for model in ("low", "high"):
        clf = pickle.load(open(os.path.join(HERE, "models", "model_%s.pkl" % model), "rb"))
        search = json.load(open(SEL[model]["source"]))
        out = {}
        for who, c in (("prod", None), ("new", clf)):
            for variant in ("incl", "test"):
                r = run(S, model, c, variant, n_target)
                json.dump(r, open(os.path.join(HERE, "eval", "%s_%s_%s.json" % (model, who, variant)), "w"), indent=1)
                out[(who, variant)] = r
        zref = {m: dc.z_fit(r) for m, r in out[("prod", "incl")]["per_mass"].items()}
        lines.append("=" * 20 + " %s-mass model: %s" % (model, ", ".join(SEL[model]["feats"])))
        for tag, r in (("search (float32)", search), ("new incl", out[("new", "incl")]), ("new test", out[("new", "test")]),
                       ("prod incl", out[("prod", "incl")]), ("prod test", out[("prod", "test")])):
            lines.append("%-17s violation %6.2f  objective %.3f" % (tag, dc.gate_violation(r), dc.objective(r, zref)))
        lines.append(" mA | Z/Zprod srch  f64 |  R srch   incl            test            | R_lo test       | A pulls (test)        | B p    | win/ovl (presel)")
        for m in out[("new", "incl")]["per_mass"]:
            s, i, t = search["per_mass"][m], out[("new", "incl")]["per_mass"][m], out[("new", "test")]["per_mass"][m]
            pt = out[("prod", "test")]["per_mass"][m]
            lines.append("%3s | %8.2f %5.2f | %6.2f  %5.2f+-%.2f  %5.2f+-%.2f (prod %4.2f) | %5.2f+-%.2f | %5.1f %5.1f %5.1f (prod %4.1f %4.1f %4.1f) | %.4f | %.2f/%.2f (%.2f/%.2f)"
                         % (m, dc.z_fit(s) / zref[m], dc.z_fit(i) / zref[m], s["R"], i["R"], i["R_err"], t["R"], t["R_err"], pt["R"],
                            t["R_lo"], t["R_lo_err"], *t["A_pulls"], *pt["A_pulls"], t["B_p"], t["win_frac"], t["overlap"],
                            t["presel_win_frac"], t["presel_overlap"]))
    txt = "\n".join(lines)
    open(os.path.join(HERE, "eval", "summary.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
