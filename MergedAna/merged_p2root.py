#!/usr/bin/env python3
"""Convert merged low-mA friend parquet -> flashggFinalFit ROOT trees.

Produces the exact `DiphotonTree` format the resolved analysis writes
(`root_MVAcut/sig/mA_M<m>/output_<year>.root`, data `.../data/<tag>/run3.root`),
so flashgg Trees2WS reads them unchanged:

  signal : DiphotonTree/ggh_125_Za_<ele|mu>_13p6TeV_cat0[_<Syst><Up|Down>01sigma]
  data   : DiphotonTree/Data_13p6TeV

Branches kept (what Trees2WS config reads): CMS_hza_mass(=H_mass), weight(=weight_central), dZ(=0).
Single category cat0 (merged has no MVA split). Lepton split via |Z_lead_lepton_id| (11=ele, 13=mu).
The 16 systematic-shifted CMS_hza_mass come from the merged_<syst>.parquet files (HiggsDNA already
re-applied each shift). Run in env `higgs-alp-ana` (uproot + pyarrow).
"""
import argparse
import glob
import os

import numpy as np
import pyarrow.parquet as pq
import uproot

EOS = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_friend"
SQRTS = "13p6TeV"
SEL = "pass_allcuts_merged_AN2020"

# merged_<key>.parquet  ->  flashgg tree-name suffix
SYST_MAP = {
    "FNUF_up": "FNUFUp01sigma", "FNUF_down": "FNUFDown01sigma",
    "Material_up": "MaterialUp01sigma", "Material_down": "MaterialDown01sigma",
    "Electron_scale_up": "ElectronScaleUp01sigma", "Electron_scale_down": "ElectronScaleDown01sigma",
    "Electron_smear_up": "ElectronSmearUp01sigma", "Electron_smear_down": "ElectronSmearDown01sigma",
    "Muon_scale_up": "MuonScaleUp01sigma", "Muon_scale_down": "MuonScaleDown01sigma",
    "Muon_smear_up": "MuonSmearUp01sigma", "Muon_smear_down": "MuonSmearDown01sigma",
    "Photon_scale_up": "PhotonScaleUp01sigma", "Photon_scale_down": "PhotonScaleDown01sigma",
    "Photon_smear_up": "PhotonSmearUp01sigma", "Photon_smear_down": "PhotonSmearDown01sigma",
}


def load(path, cols, mlo, mhi, use_weight=True):
    """Read one parquet, apply merged selection + H_mass window, return dict of arrays."""
    fs = glob.glob(path)
    if not fs:
        return None
    have = pq.ParquetFile(fs[0]).schema.names
    want = [c for c in cols if c in have]
    d = pq.read_table(fs[0], columns=want).to_pandas()
    if SEL in d.columns:
        d = d[d[SEL] == 1]
    h = d["H_mass"].to_numpy()
    m = np.isfinite(h) & (h > mlo) & (h < mhi)
    d = d[m]
    out = {"CMS_hza_mass": d["H_mass"].to_numpy().astype("float64")}
    if use_weight and "weight_central" in d.columns:
        out["weight"] = d["weight_central"].to_numpy().astype("float64")
    else:
        out["weight"] = np.ones(len(d), dtype="float64")
    out["dZ"] = np.zeros(len(d), dtype="float64")
    out["_lepid"] = (d["Z_lead_lepton_id"].to_numpy() if "Z_lead_lepton_id" in d.columns
                     else np.full(len(d), 11))
    return out


def split_lep(arrs, lep):
    """Select ele (|id|==11) or mu (|id|==13) and drop the helper column."""
    aid = np.abs(arrs["_lepid"])
    m = (aid == 11) if lep == "ele" else (aid == 13)
    return {k: v[m] for k, v in arrs.items() if k != "_lepid"}


def build_signal(tag, indir, mlo, mhi):
    base = glob.glob(f"{indir}/mA_{tag}/*/")
    if not base:
        print(f"[{tag}] no signal dir under {indir}/mA_{tag}/")
        return None
    base = base[0]
    cols = ["H_mass", "weight_central", SEL, "Z_lead_lepton_id"]
    trees = {}
    # nominal
    nom = load(f"{base}/merged_nominal.parquet", cols, mlo, mhi)
    if nom is None:
        print(f"[{tag}] missing merged_nominal.parquet")
        return None
    for lep in ("ele", "mu"):
        trees[f"DiphotonTree/ggh_125_Za_{lep}_{SQRTS}_cat0"] = split_lep(nom, lep)
    # systematics
    for key, suf in SYST_MAP.items():
        s = load(f"{base}/merged_{key}.parquet", cols, mlo, mhi)
        if s is None:
            print(f"[{tag}]   [warn] missing merged_{key}.parquet -> tree skipped")
            continue
        for lep in ("ele", "mu"):
            trees[f"DiphotonTree/ggh_125_Za_{lep}_{SQRTS}_cat0_{suf}"] = split_lep(s, lep)
    n_e = len(trees[f"DiphotonTree/ggh_125_Za_ele_{SQRTS}_cat0"]["CMS_hza_mass"])
    n_m = len(trees[f"DiphotonTree/ggh_125_Za_mu_{SQRTS}_cat0"]["CMS_hza_mass"])
    print(f"[{tag}] nominal ele={n_e} mu={n_m}; {len([k for k in trees if 'sigma' in k])} syst trees")
    return trees


def build_data(data_glob, mlo, mhi):
    """Single combined Data_13p6TeV tree (ele+mu, weight=1)."""
    fs = sorted(set(sum([glob.glob(p) for p in data_glob], [])))
    if not fs:
        print("[data] no DATA merged_nominal parquet found")
        return None
    cols = ["H_mass", SEL]
    parts = []
    for f in fs:
        a = load(f, cols, mlo, mhi, use_weight=False)
        if a:
            parts.append(a)
    if not parts:
        return None
    cat = {k: np.concatenate([p[k] for p in parts]) for k in ("CMS_hza_mass", "weight", "dZ")}
    print(f"[data] Data_{SQRTS} events in [{mlo},{mhi}] = {len(cat['CMS_hza_mass'])}")
    return {f"DiphotonTree/Data_{SQRTS}": cat}


def write_root(trees, out):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with uproot.recreate(out) as f:
        for name, arr in trees.items():
            f[name] = {k: v for k, v in arr.items()}
    print(f"  wrote {out}  ({len(trees)} trees)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default=EOS)
    ap.add_argument("--outdir", default="/eos/home-p/pelai/HZa/root_MVAcut")
    ap.add_argument("--mass-lo", type=float, default=100.0)
    ap.add_argument("--mass-hi", type=float, default=180.0)
    ap.add_argument("--signals", default="M0p1,M0p2,M0p3,M0p4,M0p5,M0p6,M0p7,M0p8,M0p9")
    ap.add_argument("--data-glob",
                    default=f"{EOS}/Data_2024/*/merged_nominal.parquet")
    ap.add_argument("--year", default="2024")
    ap.add_argument("--do-signal", action="store_true")
    ap.add_argument("--do-data", action="store_true")
    args = ap.parse_args()
    if not (args.do_signal or args.do_data):
        args.do_signal = args.do_data = True

    if args.do_signal:
        for tag in [t for t in args.signals.split(",") if t]:
            trees = build_signal(tag, args.indir, args.mass_lo, args.mass_hi)
            if trees:
                write_root(trees, f"{args.outdir}/sig/mA_{tag}/output_{args.year}.root")

    if args.do_data:
        trees = build_data([args.data_glob], args.mass_lo, args.mass_hi)
        if trees:
            # one combined data file (mirrors resolved data/<tag>/run3.root naming)
            write_root(trees, f"{args.outdir}/data/Data_merged/run3.root")


if __name__ == "__main__":
    main()
