#!/usr/bin/env python
import os
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import uproot


EVENT_COLUMNS = [
    "run",
    "event",
    "luminosityBlock",
    "weight_central",
    "fixedGridRhoAll",
    "rho",
    "dZ",
]

SELECTION_COLUMNS = [
    "n_electrons",
    "n_muons",
    "n_photons",
    "n_leptons",
    "passing_dimuon_trigger",
    "has_dimuon_candidate",
    "has_photon_candidate",
    "has_zmmg_candidate",
    "pass_tnp_presel",
    "z_mumu",
    "z_ee",
]

CANDIDATE_COLUMNS = [
    "dimuon_pt",
    "dimuon_eta",
    "dimuon_phi",
    "dimuon_mass",
    "zmmg_pt",
    "zmmg_eta",
    "zmmg_phi",
    "zmmg_mass",
    "zmmg_minus_dimuon_mass",
    "abs_zmmg_mass_minus_z",
    "minMuonGammaDR",
    "farMuonGammaDR",
    "tag_muon_1_is_near",
    "tag_muon_2_is_near",
    "tag_muon_1_probe_photon_dr",
    "tag_muon_2_probe_photon_dr",
]

TAG_MUON_FIELDS = [
    "pt",
    "eta",
    "phi",
    "mass",
    "charge",
    "ptErr",
    "ptE_error",
    "pfRelIso03_chg",
    "pfRelIso03_all",
    "miniPFRelIso_all",
    "dxy",
    "dz",
    "sip3d",
    "tightId",
    "mediumId",
    "looseId",
    "highPtId",
    "isGlobal",
    "isTracker",
    "nTrackerLayers",
]

PROBE_PHOTON_FIELDS = [
    "pt",
    "eta",
    "phi",
    "mass",
    "mvaID",
    "mvaID_WP80",
    "mvaID_WP90",
    "energyRaw",
    "energyErr",
    "r9",
    "sieie",
    "hoe",
    "hoe_PUcorr",
    "hcalPFClusterIso",
    "ecalPFClusterIso",
    "pfRelIso03_chg",
    "pfRelIso03_all",
    "pfRelIso03_chg_quadratic",
    "pfRelIso03_all_quadratic",
    "sieip",
    "etaWidth",
    "phiWidth",
    "s4",
    "trkSumPtHollowConeDR03",
    "trkSumPtSolidConeDR04",
    "pfChargedIso",
    "pfChargedIsoWorstVtx",
    "esEffSigmaRR",
    "esEnergyOverRawE",
    "chiso",
    "alliso",
    "e_veto",
    "pass_csev",
    "pass_pixel_veto",
    "lep_near_dR",
    "lep_far_dR",
    "electronVeto",
    "pixelSeed",
    "isScEtaEB",
    "isScEtaEE",
    "genPartFlav",
    "pass_ph_kinematic",
    "pass_phid_official_loose",
    "pass_phid_official_medium",
    "pass_phid_official_tight",
]

OPTIONAL_COLUMNS = {
    "dZ",
    "probe_photon_genPartFlav",
}


def get_args():
    parser = ArgumentParser(
        description="Convert flat tnp_zmmg parquet files into ROOT trees."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Path to the input parquet file."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to the output ROOT file."
    )
    parser.add_argument(
        "-s",
        "--split",
        action="store_true",
        help="Also write train/test trees split by event index parity.",
    )
    parser.add_argument(
        "--keep-all",
        action="store_true",
        help="Keep all flat parquet columns instead of the default slim TnP set.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=500000000,
        help="Kept for backward compatibility; parquet is read in one shot.",
    )
    return parser.parse_args()


def requested_columns():
    columns = EVENT_COLUMNS + SELECTION_COLUMNS + CANDIDATE_COLUMNS
    for prefix in ("tag_muon_1", "tag_muon_2"):
        columns.extend(f"{prefix}_{field}" for field in TAG_MUON_FIELDS)
    columns.extend(f"probe_photon_{field}" for field in PROBE_PHOTON_FIELDS)
    return columns


def is_non_scalar(value):
    if isinstance(value, (str, bytes)):
        return False
    return isinstance(value, (list, tuple, dict, set, np.ndarray))


def select_columns(data, keep_all):
    requested = requested_columns()
    if keep_all:
        target_columns = list(data.columns)
    else:
        target_columns = [column for column in requested if column in data.columns]

    missing_columns = []
    if not keep_all:
        missing_columns = sorted(set(requested) - set(target_columns) - OPTIONAL_COLUMNS)

    return data[target_columns].copy(), missing_columns


def sanitize_dataframe(data):
    clean_columns = {}
    dropped_columns = []

    for column in data.columns:
        series = data[column]
        sample = series.dropna().head(50)

        if not sample.empty and sample.map(is_non_scalar).any():
            dropped_columns.append(column)
            continue

        if pd.api.types.is_bool_dtype(series):
            clean_columns[column] = series.fillna(False).astype(np.int8)
            continue

        if pd.api.types.is_integer_dtype(series):
            clean_columns[column] = series.fillna(0).astype(np.int64)
            continue

        if pd.api.types.is_float_dtype(series):
            clean_columns[column] = pd.to_numeric(series, errors="coerce").astype(np.float64)
            continue

        if sample.empty:
            clean_columns[column] = pd.to_numeric(series, errors="coerce").astype(np.float64)
            continue

        if sample.map(lambda value: isinstance(value, (bool, np.bool_))).all():
            clean_columns[column] = series.fillna(False).astype(np.int8)
            continue

        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().sum() == series.notna().sum():
            finite_values = numeric.dropna().to_numpy()
            if finite_values.size > 0 and np.all(np.isclose(finite_values, np.round(finite_values))):
                clean_columns[column] = numeric.fillna(0).astype(np.int64)
            else:
                clean_columns[column] = numeric.astype(np.float64)
            continue

        dropped_columns.append(column)

    return pd.DataFrame(clean_columns), dropped_columns


def dataframe_to_tree_payload(data):
    return {column: data[column].to_numpy() for column in data.columns}


def write_output(output_path, data, split):
    if os.path.isfile(output_path):
        os.remove(output_path)

    payload = dataframe_to_tree_payload(data)

    with uproot.recreate(output_path) as root_file:
        root_file["inclusive"] = payload
        if split:
            indices = np.arange(len(data))
            train_mask = (indices % 2) == 0
            test_mask = ~train_mask
            root_file["train"] = dataframe_to_tree_payload(data.loc[train_mask])
            root_file["test"] = dataframe_to_tree_payload(data.loc[test_mask])


def main():
    args = get_args()

    print(f"Processing {args.input}")
    data = pd.read_parquet(args.input)
    input_events = len(data)

    data, missing_columns = select_columns(data, keep_all=args.keep_all)
    data, dropped_columns = sanitize_dataframe(data)

    if data.shape[1] == 0:
        raise RuntimeError("No writable flat branches found in the input parquet file.")

    write_output(args.output, data, split=args.split)

    print(f"Input events: {input_events}")
    print(f"Written events: {len(data)}")
    print(f"Written branches: {len(data.columns)}")

    if missing_columns:
        print(
            "Missing requested branches (not found in this parquet): "
            + ", ".join(missing_columns)
        )

    if dropped_columns:
        print(
            "Dropped non-flat or non-numeric branches: "
            + ", ".join(dropped_columns)
        )

    print(f"Finished writing {args.output}")


if __name__ == "__main__":
    main()
