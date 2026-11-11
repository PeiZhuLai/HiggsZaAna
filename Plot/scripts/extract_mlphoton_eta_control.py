#!/usr/bin/env python3
"""Leading MLPhoton of Z + X events outside the merged signal selection, for the
eta -> gamma gamma control-sample study (plot_mlphoton_eta_control.py).

Idea (AN2020-159 validated its regressor on the eta-meson peak): photon-like
objects in Z + jets include merged eta -> gamma gamma (m = 0.548 GeV). If the
regressor sees a peak there in data, its position and width measured in data and
in MC give the data/MC mass scale and resolution of the regressor.

Blinding: events passing pass_allcuts_merged_ML are dropped (that selection
contains the Higgs mass window). Everything else with a valid MLPhoton is kept;
the plotting script applies the diphotonScore / pT cuts.

Old MLPhoton model only: these are the friends the analysis chain uses
(current tagger, lepton cleaning). v4/v6 would need MLNanoAOD for data.

Input : data  /eos/.../HZa_merged/parquet_friend_lepclean/Data_*/**/merged_nominal.parquet
        MC    /eos/.../HZa_merged/parquet_friend_ML_lepclean_covered/Bkg_*_2024/*/merged_nominal.parquet
Output: <out>/eta_control_{data,mc}.npz  (mass, pt, eta, diphotonScore, hadronScore,
        pfIsolation, weight, era index for data / sample index for MC)

    /eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python extract_mlphoton_eta_control.py
"""
import argparse
import glob
import os

import numpy as np
import pyarrow.parquet as pq

BASE = "/eos/cms/store/group/phys_susy/pelai/HZa_merged"
DATA_GLOB = f"{BASE}/parquet_friend_lepclean/Data_*/**/merged_nominal.parquet"
MC_GLOB = f"{BASE}/parquet_friend_ML_lepclean_covered/Bkg_*_2024/*/merged_nominal.parquet"
OUT = "/eos/home-p/pelai/HZa_merged/output_plots/MLPhoton_eta_control"
FIELDS = {"mass": "mass", "pt": "pt", "eta": "eta", "ds": "diphotonScore",
          "hs": "hadronScore", "iso": "pfIsolation"}


def read(files, label_of):
    out = {k: [] for k in list(FIELDS) + ["weight", "src"]}
    labels = []
    for f in files:
        cols = ["pass_allcuts_merged_ML", "weight_central"] + [f"MLPhoton_lead_{v}" for v in FIELDS.values()]
        t = pq.read_table(f, columns=cols)
        keep = (t["MLPhoton_lead_pt"].to_numpy() > 0) & ~t["pass_allcuts_merged_ML"].to_numpy()
        for k, v in FIELDS.items():
            out[k].append(t[f"MLPhoton_lead_{v}"].to_numpy()[keep].astype(np.float32))
        out["weight"].append(t["weight_central"].to_numpy()[keep].astype(np.float32))
        lab = label_of(f)
        if lab not in labels:
            labels.append(lab)
        out["src"].append(np.full(keep.sum(), labels.index(lab), dtype=np.int16))
        print(f"  {lab:28s} {keep.sum():8d} / {t.num_rows}")
    arrays = {k: np.concatenate(v) for k, v in out.items()}
    arrays["src_names"] = np.array(labels)
    return arrays


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    data = sorted(glob.glob(DATA_GLOB, recursive=True))
    # 2022/23 tags were merged at IHEP into <tag>/merged/, 2024 ones into <tag>/<tag>_2024/
    era = lambda f: f.split("/parquet_friend_lepclean/")[1].split("_")[1]
    print(f"data: {len(data)} files")
    assert len(data) == 68, f"expected 68 data friends (32 + 36), found {len(data)}"
    np.savez(os.path.join(args.out, "eta_control_data.npz"), **read(data, era))

    mc = sorted(glob.glob(MC_GLOB))
    print(f"MC: {len(mc)} files")
    sample = lambda f: f.split("/")[-2]
    np.savez(os.path.join(args.out, "eta_control_mc.npz"), **read(mc, sample))
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
