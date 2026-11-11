#!/usr/bin/env python3
"""Reconstructed merged-photon mass per MLPhoton version, for plot_mlphoton_mass_versions.py
and plot_mlphoton_score_scan.py.

Per m_a, takes the leading MLPhoton (highest diphotonScore) of the events passing
pass_allcuts_merged_ML from the 2024 signal parquet of each version and writes one npz per m_a.
Runs in hza_ana (pyarrow); the plots run in LCG_104, which has no pyarrow.

Versions (key = parquet production, --versions key=dir overrides or adds):
  old       parquet_merged_DNA_lepclean_old  EXO-22-022 model, current tagger (lepton cleaning,
                                             strict loader): what the analysis chain uses
  old_prelc parquet_merged_DNA_tmp           EXO-22-022 model, July tagger (no lepton cleaning).
                                             Its efficiency is 3-24x lower than `old`; with the
                                             same selection it reproduces the M0p1 npz of
                                             2026-08-23 (620 / 25919 / 25920 events)
  v3, v4    parquet_merged_DNA_v3 / _v4      retrained models, August tagger (no lepton cleaning)

Input : /eos/cms/store/group/phys_susy/pelai/HZa_merged/<dir>/Sig_MC_MLNANO_{all,M1}/
        mA_MLNANO_<tag>_2024/merged_nominal.parquet
Output: <out>/<tag>_mass.npz  keys <ver> (mass), <ver>_score (diphotonScore), <ver>_eta,
        <ver>_phi (of the same MLPhoton), <ver>_key (lumi:event packed into int64, to match
        events across versions)

    /eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python extract_mlphoton_mass_versions.py
"""
import argparse
import os

import numpy as np
import pyarrow.parquet as pq

BASE = "/eos/cms/store/group/phys_susy/pelai/HZa_merged"
VERSIONS = {"old": "parquet_merged_DNA_lepclean_old", "old_prelc": "parquet_merged_DNA_tmp",
            "v3": "parquet_merged_DNA_v3", "v4": "parquet_merged_DNA_v4"}
TAGS = ["M0p1", "M0p2", "M0p3", "M0p4", "M0p5", "M0p6", "M0p7", "M0p8", "M0p9", "M1"]
# M1 was produced on its own (Sig_MC_MLNANO_M1) in the July production
SUBDIRS = ["Sig_MC_MLNANO_all", "Sig_MC_MLNANO_M1"]
OUT = "/eos/home-p/pelai/HZa_merged/output_plots/MLPhoton_mass_versions"
FIELDS = {"": "mass", "_score": "diphotonScore", "_eta": "eta", "_phi": "phi"}


def find(sub, tag, era):
    for s in SUBDIRS:
        p = os.path.join(BASE, sub, s, f"mA_MLNANO_{tag}_{era}", "merged_nominal.parquet")
        if os.path.exists(p):
            return p
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", nargs="+", default=TAGS)
    ap.add_argument("--era", default="2024")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--versions", nargs="*", default=[], metavar="KEY=DIR")
    args = ap.parse_args()
    versions = dict(VERSIONS)
    versions.update(v.split("=", 1) for v in args.versions)
    os.makedirs(args.out, exist_ok=True)
    for tag in args.tags:
        arrays = {}
        for key, sub in versions.items():
            path = find(sub, tag, args.era)
            if path is None:
                print(f"{tag} {key}: no parquet under {sub}")
                continue
            cols = ["run", "luminosityBlock", "event", "pass_allcuts_merged_ML"]
            cols += [f"MLPhoton_lead_{f}" for f in FIELDS.values()]
            t = pq.read_table(path, columns=cols)
            sel = t["pass_allcuts_merged_ML"].to_numpy()
            for suffix, f in FIELDS.items():
                arrays[key + suffix] = t[f"MLPhoton_lead_{f}"].to_numpy()[sel].astype(np.float32)
            # MC is all run 1, so lumi:event is the event id; packed as lumi * 2^40 + event,
            # which fits an int64 while lumi < 2^23 and event < 2^40 (checked)
            run, lumi, evt = (t[c].to_numpy().astype(np.int64)[sel]
                              for c in ("run", "luminosityBlock", "event"))
            assert (run == run[0]).all() if len(run) else True, f"{tag} {key}: more than one run"
            assert lumi.max(initial=0) < (1 << 23) and evt.max(initial=0) < (1 << 40)
            arrays[f"{key}_key"] = (lumi << 40) + evt
            assert len(np.unique(arrays[f"{key}_key"])) == sel.sum(), f"{tag} {key}: duplicate keys"
        if not arrays:
            continue
        dst = os.path.join(args.out, f"{tag}_mass.npz")
        np.savez(dst, **arrays)
        print(f"{tag} -> {dst}  " + "  ".join(
            f"{k} N={len(arrays[k])} med={np.median(arrays[k]):.4f}" for k in versions if k in arrays))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
