#!/usr/bin/env python3
"""Slim, self-contained samples for the HZa BDT input-variable search on IHEP batch (2026-10-13).

The IHEP jobs must not depend on HZaMVA code or CERN EOS, so everything that the production training
derives from CERN-only helpers is precomputed here:
  * training sets are selected EXACTLY as hza_features._load does (95 < m_llgg < 180, finite factor;
    background: 600k events drawn with default_rng(0) from All_Bkg 'train'), so a search job that uses the
    production inputs reproduces the production training set event by event;
  * training weights: signal factor x sideband reweight with the true-mass param (hza_features.signal_weight),
    background factor x sideband reweight (hza_features._sideband_reweight).
Evaluation sets (5 Run 3 eras, run3_bdt_scored_fsrfix):
  * signal 'test' tree (30%) of each anchor m_a, with gen dR(gamma,gamma) (first entry) for the overlap tag;
  * DY simulation, inclusive tree, sculpt-R combination (raw factor);
  * data, inclusive tree, SIDEBANDS ONLY (95-115, 135-180): the 115-135 window is not copied;
  * production scores MVA_Score_mA_M<m> and the production working points, for the equal-data-count
    reference (data sideband count of the production model at its working point).
Columns: the variable pool (all free of m_llgg) + the production *_oHm/param inputs (baseline only) +
bookkeeping. float32, one .npz per set. Env: higgs-alp-ana.
Output: <out>/{train_sig_m<m>, train_bkg, test_sig_m<m>, dy, data_sb}.npz, meta.json
"""
import argparse
import json
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
import hza_features as hf  # noqa: E402

INP = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
SCO = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
WPJ = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
DY = {"2022preEE": ["DYGto2LG_10to100", "DYJetsToLL"], "2022postEE": ["DYGto2LG_10to100", "DYJetsToLL"],
      "2023preBPix": ["DYGto2LG_10to100", "DYJetsToLL"], "2023postBPix": ["DYGto2LG_10to100", "DYJetsToLL"],
      "2024": ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]}
ANCHORS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
PH = ["mvaID", "eta", "energyErr", "energyRaw", "hoe", "s4", "sieip", "etaWidth", "phiWidth", "esEffSigmaRR",
      "chiso", "ecalPFClusterIso", "hcalPFClusterIso", "pfChargedIsoWorstVtx", "trkSumPtHollowConeDR03",
      "lep_near_dR"]
POOL_RAW = (["pho1Pt", "pho2Pt", "pho1R9", "pho2R9", "pho1IetaIeta55", "pho2IetaIeta55", "pho1PIso_noCorr",
             "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "ALP_PhotonIso", "ALP_pt", "ALP_eta", "var_PtaOverMa",
             "var_dR_g1g2", "var_dR_Za", "var_dR_g1Z", "var_dEta_g1Z", "var_dPhi_g1Z", "min_lg_deltaR",
             "max_lg_deltaR", "H_pt"]
            + [f"ALP_{w}_photon_{v}" for w in ("lead", "sublead") for v in PH])
BASELINE = ["pho1Pt_oHm", "pho2Pt_oHm", "var_PtaOverMh", "H_pt_oHm"]
BOOK = ["H_m", "ALP_m", "factor"]
COLS = BOOK + POOL_RAW + BASELINE


def sel_mask(a):
    return (a["H_m"] > 95) & (a["H_m"] < 180) & np.isfinite(a["factor"])


def save(out, name, d):
    p = os.path.join(out, name + ".npz")
    np.savez_compressed(p + ".partial.npz", **{k: np.asarray(v, dtype=np.float32) for k, v in d.items()})
    os.replace(p + ".partial.npz", p)
    n = len(next(iter(d.values())))
    print(f"  {name}: {n} events, {len(d)} columns", flush=True)
    return n


def train_set(path, tree, nmax, extra=()):
    """Same events as hza_features._load (mask, then default_rng(0) subsample), all needed columns."""
    t = uproot.open(path)[tree]
    a = t.arrays(sorted(set(COLS + list(hf.READ) + list(hf.REWEIGHT_EXTRA) + list(extra))), library="np")
    m = sel_mask(a)
    a = {k: v[m].astype(np.float64) for k, v in a.items()}
    if nmax and len(a["H_m"]) > nmax:
        idx = np.random.default_rng(0).choice(len(a["H_m"]), nmax, replace=False)
        a = {k: v[idx] for k, v in a.items()}
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    meta = {"pool_raw": POOL_RAW, "baseline_cols": BASELINE, "anchors": ANCHORS, "eras": ERAS, "counts": {}}
    wp = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(WPJ))["results"]}
    meta["production_wp"] = wp

    print("[train] signal anchors", flush=True)
    for m in ANCHORS:
        a = train_set(f"{INP}/mA_M{m}/run3.root", "train", None)
        ref = hf._load(f"{INP}/mA_M{m}/run3.root", "train", None, extra=hf.REWEIGHT_EXTRA)
        assert np.array_equal(ref["H_m"], a["H_m"]) and np.array_equal(ref["event"], a["event"]), f"mA{m} selection differs"
        a["w_train"] = hf.signal_weight(ref, m)
        meta["counts"][f"train_sig_m{m}"] = save(args.out, f"train_sig_m{m}", {k: a[k] for k in COLS + ["w_train"]})
    print("[train] background", flush=True)
    a = train_set(f"{INP}/All_Bkg/run3.root", "train", 600000)
    ref = hf._load(f"{INP}/All_Bkg/run3.root", "train", 600000, extra=hf.REWEIGHT_EXTRA)
    assert np.array_equal(ref["H_m"], a["H_m"]) and np.array_equal(ref["event"], a["event"]), "background selection differs"
    a["w_train"] = ref["factor"] * hf._sideband_reweight(ref)
    meta["counts"]["train_bkg"] = save(args.out, "train_bkg", {k: a[k] for k in COLS + ["w_train"]})

    print("[eval] signal test trees", flush=True)
    for m in ANCHORS:
        acc = {}
        for e in ERAS:
            t = uproot.open(f"{SCO}/mA_M{m}/{e}.root")["test"]
            a = t.arrays(COLS + ["GenALP_dR_gg", f"MVA_Score_mA_M{m}"], library="np")
            g = np.array([x[0] if len(x) else np.nan for x in a.pop("GenALP_dR_gg")], float)
            msk = sel_mask(a)
            a = {k: v[msk] for k, v in a.items()}; a["gen_dR_gg"] = g[msk]; a["era"] = np.full(msk.sum(), ERAS.index(e))
            a["prod_score"] = a.pop(f"MVA_Score_mA_M{m}")
            for k, v in a.items():
                acc.setdefault(k, []).append(v)
        meta["counts"][f"test_sig_m{m}"] = save(args.out, f"test_sig_m{m}", {k: np.concatenate(v) for k, v in acc.items()})

    print("[eval] DY and data sidebands", flush=True)
    scores = [f"MVA_Score_mA_M{m}" for m in ANCHORS]
    for name, pairs, sb_only in (("dy", [(s, e) for e in ERAS for s in DY[e]], False), ("data_sb", [("Data", e) for e in ERAS], True)):
        acc = {}
        for s, e in pairs:
            a = uproot.open(f"{SCO}/{s}/{e}.root")["inclusive"].arrays(COLS + scores, library="np")
            msk = sel_mask(a)
            if sb_only:
                h = a["H_m"]; msk &= ((h > 95) & (h < 115)) | ((h > 135) & (h < 180))
            a = {k: v[msk] for k, v in a.items()}; a["era"] = np.full(msk.sum(), ERAS.index(e))
            for k, v in a.items():
                acc.setdefault(k, []).append(v)
        meta["counts"][name] = save(args.out, name, {k: np.concatenate(v) for k, v in acc.items()})

    # equal-data-count reference: production data sideband count at the production working point
    d = np.load(os.path.join(args.out, "data_sb.npz"))
    meta["prod_data_sb_count"] = {m: int((d[f"MVA_Score_mA_M{m}"] > wp[m]).sum()) for m in ANCHORS}
    print("  production data sideband counts:", meta["prod_data_sb_count"], flush=True)
    json.dump(meta, open(os.path.join(args.out, "meta.json"), "w"), indent=1)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
