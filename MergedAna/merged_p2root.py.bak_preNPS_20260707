#!/usr/bin/env python3
"""Convert merged low-mA ML friend parquet -> flashggFinalFit ROOT trees,
with a PER-mA ROI on the ML-regressed merged-photon mass.

This is the CORRECT per-mA merged analysis (supersedes the mA-inclusive
`pass_allcuts_merged_AN2020` + single-multipdf path used in progress-doc s.11):

  * Selection : pass_allcuts_merged_ML == 1  (reco ML merged candidate).
  * ROI       : per-mA window on MLPhoton_lead_mass (reco lead ML-photon mass,
                available in BOTH signal MC and data -- NOTE MergedML_mass is
                gen-truth-matched and MC-only, so it MUST NOT be used for data).
                Each m_a selects a different sub-sample -> per-mA background,
                exactly like the resolved analysis bins by m_gamma-gamma ~ m_a.
  * Observable: CMS_hza_mass = H_merged_ML_mass = (Z + MLPhoton_lead).mass,
                recomputed offline from the stored Z_* and MLPhoton_lead_*
                4-vectors (it is NOT saved in the parquet; the parquet's own
                H_mass is the resolved candidate and peaks ~135, wrong here).

Output = the exact flashgg DiphotonTree format Trees2WS reads unchanged:
  signal : sig/mA_<tag>/output_<year>.root
             DiphotonTree/ggh_125_Za_<ele|mu>_13p6TeV_cat0[_<Syst><Up|Down>01sigma]
  data   : data/Data_merged_<tag>/run3.root          (PER-mA now)
             DiphotonTree/Data_13p6TeV

Branches kept (what Trees2WS config reads): CMS_hza_mass, weight(=weight_central,
data=1), dZ(=0). Single category cat0 (merged has no MVA split). Lepton split via
|Z_lead_lepton_id| (11=ele, 13=mu). Run in env `higgs-alp-ana` (uproot + pyarrow).
"""
import argparse
import glob
import os

import numpy as np
import pyarrow.parquet as pq
import uproot

SQRTS = "13p6TeV"
SEL = "pass_allcuts_merged_ML"
# opt-1 verification toggle: when SEL_MERGED_GE1=1, recompute the merged
# selection as (n_MLPhoton_diphoton >= 1) directly from the parquet instead of
# reading the stored pass_allcuts_merged_ML flag. For SIGNAL this is exactly
# equivalent to re-running HiggsDNA with the >=1 tagger change (verified:
# pass_allcuts_merged_ML == (n_dip==1) for every mass, and all n_dip events pass
# the upstream trigger/Z/H cuts). Default OFF so DATA (which has genuine upstream
# failures) keeps using the full stored flag. Once HiggsDNA is re-run with the
# >=1 tagger, the stored flag is already correct and this toggle is unnecessary.
SEL_GE1 = os.environ.get("SEL_MERGED_GE1", "0") == "1"
ROI_VAR = "MLPhoton_lead_mass"

# Signal ML parquet (all 9 sub-GeV mass points; nominal + 16 syst per tag).
SIG_ML = "/eos/home-p/pelai/HZa/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all"
# Data ML friend parquet (per-dataset dirs from run_merged_data_ml_friend.sh).
DATA_ML_GLOB = ("/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_friend/"
                "Data_2024_*/*/merged_nominal.parquet")

# Per-mA ROI window on MLPhoton_lead_mass = [q16, q84] of the signal reco
# distribution (contains the central 68% of signal; tracks m_a). Derived
# 2026-07-01 from Sig_MC_MLNANO_all (pass_allcuts_merged_ML).
ROI_WINDOWS = {
    "M0p1": (0.154, 0.947), "M0p2": (0.176, 0.887), "M0p3": (0.226, 0.917),
    "M0p4": (0.326, 0.796), "M0p5": (0.413, 0.821), "M0p6": (0.489, 0.878),
    "M0p7": (0.583, 0.925), "M0p8": (0.666, 0.996), "M0p9": (0.747, 1.079),
}

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

_Z = ["Z_pt", "Z_eta", "Z_phi", "Z_mass"]
_ML = ["MLPhoton_lead_pt", "MLPhoton_lead_eta", "MLPhoton_lead_phi", "MLPhoton_lead_mass"]
_COLS = _Z + _ML + [SEL, "n_MLPhoton_diphoton", "weight_central", "Z_lead_lepton_id"]


def _p4(pt, eta, phi, m):
    px = pt * np.cos(phi); py = pt * np.sin(phi); pz = pt * np.sinh(eta)
    E = np.sqrt(px * px + py * py + pz * pz + m * m)
    return E, px, py, pz


def _h_merged_ml_mass(d):
    """CMS_hza_mass = (Z + MLPhoton_lead).mass, recomputed per event."""
    Ez, zx, zy, zz = _p4(d["Z_pt"].to_numpy(), d["Z_eta"].to_numpy(),
                         d["Z_phi"].to_numpy(), d["Z_mass"].to_numpy())
    Eg, gx, gy, gz = _p4(d["MLPhoton_lead_pt"].to_numpy(), d["MLPhoton_lead_eta"].to_numpy(),
                         d["MLPhoton_lead_phi"].to_numpy(), d["MLPhoton_lead_mass"].to_numpy())
    E = Ez + Eg; px = zx + gx; py = zy + gy; pz = zz + gz
    return np.sqrt(np.clip(E * E - px * px - py * py - pz * pz, 0.0, None))


def load(path, tag, mlo, mhi, use_weight=True):
    """Read one parquet -> pass_ML + per-mA ROI + recomputed observable."""
    fs = glob.glob(path)
    if not fs:
        return None
    have = pq.ParquetFile(fs[0]).schema.names
    want = [c for c in _COLS if c in have]
    d = pq.read_table(fs[0], columns=want).to_pandas()
    if SEL_GE1 and "n_MLPhoton_diphoton" in d.columns:
        d = d[d["n_MLPhoton_diphoton"] >= 1]
    else:
        if SEL not in d.columns:
            return None
        d = d[d[SEL] == 1]
    if ROI_VAR not in d.columns:
        return None
    lo, hi = ROI_WINDOWS[tag]
    roi = d[ROI_VAR].to_numpy()
    d = d[(roi > lo) & (roi < hi)]
    if len(d) == 0:
        return {"CMS_hza_mass": np.array([], "float64"), "weight": np.array([], "float64"),
                "dZ": np.array([], "float64"), "_lepid": np.array([], "float64")}
    h = _h_merged_ml_mass(d)
    m = np.isfinite(h) & (h > mlo) & (h < mhi)
    out = {"CMS_hza_mass": h[m].astype("float64")}
    if use_weight and "weight_central" in d.columns:
        out["weight"] = d["weight_central"].to_numpy()[m].astype("float64")
    else:
        out["weight"] = np.ones(int(m.sum()), dtype="float64")
    out["dZ"] = np.zeros(int(m.sum()), dtype="float64")
    out["_lepid"] = (d["Z_lead_lepton_id"].to_numpy()[m] if "Z_lead_lepton_id" in d.columns
                     else np.full(int(m.sum()), 11))
    return out


def split_lep(arrs, lep):
    aid = np.abs(arrs["_lepid"])
    m = (aid == 11) if lep == "ele" else (aid == 13)
    return {k: v[m] for k, v in arrs.items() if k != "_lepid"}


def build_signal(tag, indir, mlo, mhi):
    base = f"{indir}/mA_MLNANO_{tag}_2024"
    if not os.path.isdir(base):
        print(f"[{tag}] no signal ML dir {base}")
        return None
    trees = {}
    nom = load(f"{base}/merged_nominal.parquet", tag, mlo, mhi)
    if nom is None:
        print(f"[{tag}] missing/invalid merged_nominal.parquet")
        return None
    for lep in ("ele", "mu"):
        trees[f"DiphotonTree/ggh_125_Za_{lep}_{SQRTS}_cat0"] = split_lep(nom, lep)
    for key, suf in SYST_MAP.items():
        s = load(f"{base}/merged_{key}.parquet", tag, mlo, mhi)
        if s is None:
            print(f"[{tag}]   [warn] missing merged_{key}.parquet -> tree skipped")
            continue
        for lep in ("ele", "mu"):
            trees[f"DiphotonTree/ggh_125_Za_{lep}_{SQRTS}_cat0_{suf}"] = split_lep(s, lep)
    n_e = len(trees[f"DiphotonTree/ggh_125_Za_ele_{SQRTS}_cat0"]["CMS_hza_mass"])
    n_m = len(trees[f"DiphotonTree/ggh_125_Za_mu_{SQRTS}_cat0"]["CMS_hza_mass"])
    print(f"[{tag}] ROI={ROI_WINDOWS[tag]} nominal ele={n_e} mu={n_m}; "
          f"{len([k for k in trees if 'sigma' in k])} syst trees")
    return trees


def build_data(tag, data_glob, mlo, mhi):
    """Single combined Data_13p6TeV tree for this m_a (ele+mu, weight=1)."""
    fs = sorted(set(glob.glob(data_glob)))
    if not fs:
        print(f"[data {tag}] no DATA ML parquet found ({data_glob})")
        return None
    parts = []
    for f in fs:
        a = load(f, tag, mlo, mhi, use_weight=False)
        if a and len(a["CMS_hza_mass"]):
            parts.append(a)
    if not parts:
        print(f"[data {tag}] 0 events in ROI")
        cat = {k: np.array([], "float64") for k in ("CMS_hza_mass", "weight", "dZ")}
        return {f"DiphotonTree/Data_{SQRTS}": cat}
    cat = {k: np.concatenate([p[k] for p in parts]) for k in ("CMS_hza_mass", "weight", "dZ")}
    print(f"[data {tag}] ROI={ROI_WINDOWS[tag]} Data_{SQRTS} events in [{mlo},{mhi}] = {len(cat['CMS_hza_mass'])}")
    return {f"DiphotonTree/Data_{SQRTS}": cat}


def write_root(trees, out):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with uproot.recreate(out) as f:
        for name, arr in trees.items():
            f[name] = {k: v for k, v in arr.items()}
    print(f"  wrote {out}  ({len(trees)} trees)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default=SIG_ML, help="signal ML parquet base (Sig_MC_MLNANO_all)")
    ap.add_argument("--outdir", default="/eos/home-p/pelai/HZa/root_MVAcut")
    ap.add_argument("--mass-lo", type=float, default=100.0)
    ap.add_argument("--mass-hi", type=float, default=180.0)
    ap.add_argument("--signals", default="M0p1,M0p2,M0p3,M0p4,M0p5,M0p6,M0p7,M0p8,M0p9")
    ap.add_argument("--data-glob", default=DATA_ML_GLOB)
    ap.add_argument("--year", default="2024")
    ap.add_argument("--do-signal", action="store_true")
    ap.add_argument("--do-data", action="store_true")
    args = ap.parse_args()
    if not (args.do_signal or args.do_data):
        args.do_signal = args.do_data = True

    tags = [t for t in args.signals.split(",") if t]
    for tag in tags:
        if tag not in ROI_WINDOWS:
            print(f"[{tag}] no ROI window defined -> skip")
            continue
        if args.do_signal:
            trees = build_signal(tag, args.indir, args.mass_lo, args.mass_hi)
            if trees:
                write_root(trees, f"{args.outdir}/sig/mA_{tag}/output_{args.year}.root")
        if args.do_data:
            trees = build_data(tag, args.data_glob, args.mass_lo, args.mass_hi)
            if trees:
                # PER-mA data file (ROI selects a different sub-sample per m_a).
                # Path mirrors signal (sig/mA_<tag>) so makeYields' _inputWSFile
                # (data/mA_M<mass>/ws/run3.root) resolves to a real file.
                write_root(trees, f"{args.outdir}/data/mA_{tag}/run3.root")


if __name__ == "__main__":
    main()
