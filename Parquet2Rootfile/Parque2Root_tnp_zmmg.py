#!/usr/bin/env python
import os
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import uproot


TREE_PATH = "tnpPhoIDs/fitter_tree"

# histUtils.cpp directly binds pair_mass to a C++ float buffer via SetBranchAddress.
# Keep this branch as float32 in the ROOT output to match that expectation.
ROOT_DTYPE_OVERRIDES = {
    "pair_mass": np.float32,
}

# Other branches are accessed through TTreeFormula, so they only need to be present
# and numeric rather than matching one exact C++ primitive type.
HISTUTILS_FORMULA_COLUMNS = (
    "totWeight",
    "ph_passElectronVeto",
    "ph_r9",
    "ph_et",
    "event_nPV",
    "event_met_pfmet",
    "event_met_pfphi",
    "tag_Ele_pt",
    "tag_Ele_phi",
    "mcTrue",
)

EVENT_SOURCE_COLUMNS = [
    "run",
    "event",
    "luminosityBlock",
    "weight_central",
    "fixedGridRhoAll",
    "rho",
    "PV",
    "PV_z",
    "event_nPV",
    "dZ",
    "pass_tnp_presel",
    "passing_dimuon_trigger",
    "has_zmmg_candidate",
    "zmmg_minus_dimuon_mass",
    "abs_zmmg_mass_minus_z",
    "minMuonGammaDR",
    "farMuonGammaDR",
]

DIMUON_SOURCE_COLUMNS = [
    "dimuon_pt",
    "dimuon_eta",
    "dimuon_phi",
    "dimuon_mass",
    "zmmg_pt",
    "zmmg_eta",
    "zmmg_phi",
    "zmmg_mass",
]

TAG_MUON_SOURCE_COLUMNS = [
    "tag_muon_1_pt",
    "tag_muon_1_eta",
    "tag_muon_1_phi",
    "tag_muon_1_mass",
    "tag_muon_1_charge",
    "tag_muon_1_probe_photon_dr",
    "tag_muon_2_pt",
    "tag_muon_2_eta",
    "tag_muon_2_phi",
    "tag_muon_2_mass",
    "tag_muon_2_charge",
    "tag_muon_2_probe_photon_dr",
]

PROBE_PHOTON_SOURCE_COLUMNS = [
    "probe_photon_pt",
    "probe_photon_eta",
    "probe_photon_phi",
    "probe_photon_mass",
    "probe_photon_mvaID",
    "probe_photon_r9",
    "probe_photon_sieie",
    "probe_photon_sieip",
    "probe_photon_etaWidth",
    "probe_photon_phiWidth",
    "probe_photon_s4",
    "probe_photon_hoe",
    "probe_photon_hoe_PUcorr",
    "probe_photon_hcalPFClusterIso",
    "probe_photon_ecalPFClusterIso",
    "probe_photon_chiso",
    "probe_photon_alliso",
    "probe_photon_pfChargedIso",
    "probe_photon_pfChargedIsoWorstVtx",
    "probe_photon_trkSumPtHollowConeDR03",
    "probe_photon_trkSumPtSolidConeDR04",
    "probe_photon_esEffSigmaRR",
    "probe_photon_esEnergyOverRawE",
    "probe_photon_electronVeto",
    "probe_photon_pixelSeed",
    "probe_photon_pass_csev",
    "probe_photon_pass_pixel_veto",
    "probe_photon_pass_mva_min",
    "probe_photon_isScEtaEB",
    "probe_photon_isScEtaEE",
    "probe_photon_lep_near_dR",
    "probe_photon_lep_far_dR",
    "probe_photon_genPartFlav",
]

REQUESTED_SOURCE_COLUMNS = list(
    dict.fromkeys(
        EVENT_SOURCE_COLUMNS
        + DIMUON_SOURCE_COLUMNS
        + TAG_MUON_SOURCE_COLUMNS
        + PROBE_PHOTON_SOURCE_COLUMNS
    )
)

OPTIONAL_SOURCE_COLUMNS = {
    "PV",
    "PV_z",
    "dZ",
    "probe_photon_genPartFlav",
}


def get_args():
    parser = ArgumentParser(
        description=(
            "Convert flat tnp_zmmg parquet files into egm_tnp_analysis-style ROOT trees."
        )
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
        help="Also write event-parity split trees under tnpPhoIDs_train/test.",
    )
    parser.add_argument(
        "--keep-all",
        action="store_true",
        help="Keep the original flat parquet columns in addition to the EGM aliases.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=500000000,
        help="Kept for backward compatibility; parquet is read in one shot.",
    )
    return parser.parse_args()


def is_non_scalar(value):
    if isinstance(value, (str, bytes)):
        return False
    return isinstance(value, (list, tuple, dict, set, np.ndarray))


def select_columns(data, keep_all):
    if keep_all:
        target_columns = list(data.columns)
    else:
        target_columns = [column for column in REQUESTED_SOURCE_COLUMNS if column in data.columns]

    missing_columns = sorted(
        set(REQUESTED_SOURCE_COLUMNS) - set(target_columns) - OPTIONAL_SOURCE_COLUMNS
    )
    return data[target_columns].copy(), missing_columns


def source_series(data, column, default):
    if column in data.columns:
        return data[column].copy()
    return pd.Series(default, index=data.index)


def first_existing_series(data, columns, default):
    for column in columns:
        if column in data.columns:
            return data[column].copy()
    return pd.Series(default, index=data.index)


def alias_series(output, name, values):
    output[name] = values


def build_egm_compatible_dataframe(data, keep_all):
    output = data.copy() if keep_all else pd.DataFrame(index=data.index)

    rho = first_existing_series(data, ["rho", "fixedGridRhoAll"], 0.0)
    pv = first_existing_series(data, ["PV", "PV_z"], np.nan)
    weight = source_series(data, "weight_central", 1.0)
    tag_pt = source_series(data, "dimuon_pt", np.nan)
    tag_eta = source_series(data, "dimuon_eta", np.nan)
    tag_phi = source_series(data, "dimuon_phi", np.nan)
    tag_mass = source_series(data, "dimuon_mass", np.nan)
    probe_eta = source_series(data, "probe_photon_eta", np.nan)
    probe_phi = source_series(data, "probe_photon_phi", np.nan)
    probe_csev = source_series(data, "probe_photon_pass_csev", False)
    probe_pixel_veto = source_series(data, "probe_photon_pass_pixel_veto", False)
    probe_mva_min = source_series(data, "probe_photon_pass_mva_min", False)
    probe_gen_part_flav = source_series(data, "probe_photon_genPartFlav", 1)

    alias_series(output, "run", source_series(data, "run", 0))
    alias_series(output, "event", source_series(data, "event", 0))
    alias_series(output, "luminosityBlock", source_series(data, "luminosityBlock", 0))
    alias_series(output, "weight_central", weight)
    alias_series(output, "weight", weight)
    alias_series(output, "totWeight", weight)
    alias_series(output, "fixedGridRhoAll", source_series(data, "fixedGridRhoAll", rho))
    alias_series(output, "rho", rho)
    alias_series(output, "event_rho", rho)
    alias_series(output, "PV", pv)
    alias_series(output, "PV_z", source_series(data, "PV_z", pv))
    alias_series(output, "event_nPV", source_series(data, "event_nPV", np.nan))
    alias_series(output, "event_met_pfmet", 0.0)
    alias_series(output, "event_met_pfphi", 0.0)
    alias_series(output, "dZ", source_series(data, "dZ", 0.0))
    alias_series(output, "pass_tnp_presel", source_series(data, "pass_tnp_presel", True))
    alias_series(
        output,
        "passing_dimuon_trigger",
        source_series(data, "passing_dimuon_trigger", True),
    )
    alias_series(
        output,
        "has_zmmg_candidate",
        source_series(data, "has_zmmg_candidate", True),
    )

    alias_series(output, "pair_mass", source_series(data, "zmmg_mass", np.nan))
    alias_series(output, "z_mass", source_series(data, "dimuon_mass", np.nan))
    alias_series(
        output,
        "zmmg_minus_dimuon_mass",
        source_series(data, "zmmg_minus_dimuon_mass", np.nan),
    )
    alias_series(
        output,
        "abs_zmmg_mass_minus_z",
        source_series(data, "abs_zmmg_mass_minus_z", np.nan),
    )
    alias_series(output, "tag_mumu_pt", tag_pt)
    alias_series(output, "tag_mumu_eta", tag_eta)
    alias_series(output, "tag_mumu_phi", tag_phi)
    alias_series(output, "tag_mumu_mass", tag_mass)
    alias_series(output, "tag_pt", tag_pt)
    alias_series(output, "tag_eta", tag_eta)
    alias_series(output, "tag_phi", tag_phi)
    alias_series(output, "tag_mass", tag_mass)
    alias_series(output, "tag_sc_eta", tag_eta)
    alias_series(output, "tag_sc_abseta", np.abs(tag_eta))
    alias_series(output, "tag_sc_phi", tag_phi)
    alias_series(output, "tag_Ele_pt", tag_pt)
    alias_series(output, "tag_Ele_eta", tag_eta)
    alias_series(output, "tag_Ele_phi", tag_phi)
    alias_series(output, "tag_Ele_mass", tag_mass)

    alias_series(output, "tag_muon_1_pt", source_series(data, "tag_muon_1_pt", np.nan))
    alias_series(output, "tag_muon_1_eta", source_series(data, "tag_muon_1_eta", np.nan))
    alias_series(output, "tag_muon_1_phi", source_series(data, "tag_muon_1_phi", np.nan))
    alias_series(output, "tag_muon_1_mass", source_series(data, "tag_muon_1_mass", np.nan))
    alias_series(
        output, "tag_muon_1_charge", source_series(data, "tag_muon_1_charge", 0)
    )
    alias_series(
        output,
        "tag_muon_1_probe_photon_dr",
        source_series(data, "tag_muon_1_probe_photon_dr", np.nan),
    )
    alias_series(output, "tag_muon_2_pt", source_series(data, "tag_muon_2_pt", np.nan))
    alias_series(output, "tag_muon_2_eta", source_series(data, "tag_muon_2_eta", np.nan))
    alias_series(output, "tag_muon_2_phi", source_series(data, "tag_muon_2_phi", np.nan))
    alias_series(output, "tag_muon_2_mass", source_series(data, "tag_muon_2_mass", np.nan))
    alias_series(
        output, "tag_muon_2_charge", source_series(data, "tag_muon_2_charge", 0)
    )
    alias_series(
        output,
        "tag_muon_2_probe_photon_dr",
        source_series(data, "tag_muon_2_probe_photon_dr", np.nan),
    )
    alias_series(
        output, "minMuonGammaDR", source_series(data, "minMuonGammaDR", np.nan)
    )
    alias_series(
        output, "farMuonGammaDR", source_series(data, "farMuonGammaDR", np.nan)
    )

    alias_series(output, "ph_et", source_series(data, "probe_photon_pt", np.nan))
    alias_series(output, "ph_eta", probe_eta)
    alias_series(output, "ph_phi", probe_phi)
    alias_series(output, "ph_sc_eta", probe_eta)
    alias_series(output, "ph_sc_abseta", np.abs(probe_eta))
    alias_series(output, "ph_sc_phi", probe_phi)
    alias_series(output, "ph_mass", source_series(data, "probe_photon_mass", 0.0))
    alias_series(output, "ph_mvaID", source_series(data, "probe_photon_mvaID", np.nan))
    alias_series(
        output, "ph_mva122XV1", source_series(data, "probe_photon_mvaID", np.nan)
    )
    alias_series(output, "ph_r9", source_series(data, "probe_photon_r9", np.nan))
    alias_series(output, "ph_sieie", source_series(data, "probe_photon_sieie", np.nan))
    alias_series(output, "ph_sieip", source_series(data, "probe_photon_sieip", np.nan))
    alias_series(
        output, "ph_etaWidth", source_series(data, "probe_photon_etaWidth", np.nan)
    )
    alias_series(
        output, "ph_phiWidth", source_series(data, "probe_photon_phiWidth", np.nan)
    )
    alias_series(output, "ph_s4", source_series(data, "probe_photon_s4", np.nan))
    alias_series(output, "ph_hoe", source_series(data, "probe_photon_hoe", np.nan))
    alias_series(
        output,
        "ph_hoe_PUcorr",
        source_series(data, "probe_photon_hoe_PUcorr", np.nan),
    )
    alias_series(
        output,
        "ph_chIso",
        source_series(data, "probe_photon_chiso", np.nan),
    )
    alias_series(
        output,
        "ph_neuIso",
        source_series(data, "probe_photon_hcalPFClusterIso", np.nan),
    )
    alias_series(
        output,
        "ph_ecalIso",
        source_series(data, "probe_photon_ecalPFClusterIso", np.nan),
    )
    alias_series(
        output,
        "ph_allIso",
        source_series(data, "probe_photon_alliso", np.nan),
    )
    alias_series(
        output,
        "ph_pfChargedIso",
        source_series(data, "probe_photon_pfChargedIso", np.nan),
    )
    alias_series(
        output,
        "ph_pfChargedIsoWorstVtx",
        source_series(data, "probe_photon_pfChargedIsoWorstVtx", np.nan),
    )
    alias_series(
        output,
        "ph_trkSumPtHollowConeDR03",
        source_series(data, "probe_photon_trkSumPtHollowConeDR03", np.nan),
    )
    alias_series(
        output,
        "ph_trkSumPtSolidConeDR04",
        source_series(data, "probe_photon_trkSumPtSolidConeDR04", np.nan),
    )
    alias_series(
        output,
        "ph_esEffSigmaRR",
        source_series(data, "probe_photon_esEffSigmaRR", np.nan),
    )
    alias_series(
        output,
        "ph_esEnergyOverRawE",
        source_series(data, "probe_photon_esEnergyOverRawE", np.nan),
    )
    alias_series(output, "ph_passElectronVeto", probe_csev)
    alias_series(output, "ph_passCSEV", probe_csev)
    alias_series(output, "ph_passPixelVeto", probe_pixel_veto)
    alias_series(output, "probe_photon_pass_mva_min", probe_mva_min)
    alias_series(output, "ph_passMVAMin", probe_mva_min)
    alias_series(
        output,
        "ph_electronVeto",
        source_series(data, "probe_photon_electronVeto", False),
    )
    alias_series(
        output, "ph_pixelSeed", source_series(data, "probe_photon_pixelSeed", False)
    )
    alias_series(
        output, "ph_isEB", source_series(data, "probe_photon_isScEtaEB", False)
    )
    alias_series(
        output, "ph_isEE", source_series(data, "probe_photon_isScEtaEE", False)
    )
    alias_series(output, "ph_lepNearDR", source_series(data, "probe_photon_lep_near_dR", np.nan))
    alias_series(output, "ph_lepFarDR", source_series(data, "probe_photon_lep_far_dR", np.nan))
    alias_series(output, "passingProbeCSEV", probe_csev)
    alias_series(output, "passingProbePixelVeto", probe_pixel_veto)
    alias_series(output, "ph_genPartFlav", probe_gen_part_flav)
    alias_series(output, "mcTrue", probe_gen_part_flav.eq(1))

    return output


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
    payload = {}
    for column in data.columns:
        values = data[column].to_numpy()
        target_dtype = ROOT_DTYPE_OVERRIDES.get(column)
        if target_dtype is not None:
            values = values.astype(target_dtype, copy=False)
        payload[column] = values
    return payload


def validate_histutils_payload(payload):
    issues = []

    pair_mass = payload.get("pair_mass")
    if pair_mass is None:
        issues.append("missing required branch 'pair_mass'")
    elif pair_mass.dtype != np.dtype(np.float32):
        issues.append(
            "branch 'pair_mass' must be float32 for libPython/histUtils.cpp, "
            f"got {pair_mass.dtype}"
        )

    for column in HISTUTILS_FORMULA_COLUMNS:
        values = payload.get(column)
        if values is None:
            issues.append(f"missing formula branch '{column}'")
            continue
        if values.dtype.kind not in "biuf":
            issues.append(
                f"formula branch '{column}' must be numeric, got {values.dtype}"
            )

    return issues


def write_output(output_path, data, split):
    if os.path.isfile(output_path):
        os.remove(output_path)

    payload = dataframe_to_tree_payload(data)
    compatibility_issues = validate_histutils_payload(payload)
    if compatibility_issues:
        raise RuntimeError(
            "histUtils compatibility check failed: "
            + "; ".join(compatibility_issues)
        )

    with uproot.recreate(output_path) as root_file:
        root_file[TREE_PATH] = payload

        if split:
            indices = np.arange(len(data))
            train_mask = (indices % 2) == 0
            test_mask = ~train_mask
            root_file["tnpPhoIDs_train/fitter_tree"] = dataframe_to_tree_payload(
                data.loc[train_mask]
            )
            root_file["tnpPhoIDs_test/fitter_tree"] = dataframe_to_tree_payload(
                data.loc[test_mask]
            )


def main():
    args = get_args()

    print(f"Processing {args.input}")
    data = pd.read_parquet(args.input)
    input_events = len(data)

    data, missing_columns = select_columns(data, keep_all=args.keep_all)
    data = build_egm_compatible_dataframe(data, keep_all=args.keep_all)
    data, dropped_columns = sanitize_dataframe(data)

    if data.shape[1] == 0:
        raise RuntimeError("No writable flat branches found in the input parquet file.")

    write_output(args.output, data, split=args.split)

    print(f"Input events: {input_events}")
    print(f"Written events: {len(data)}")
    print(f"Written branches: {len(data.columns)}")
    print(f"Primary tree: {TREE_PATH}")
    print(
        "histUtils compatibility checks passed for pair_mass and formula branches."
    )

    if missing_columns:
        print(
            "Missing requested parquet branches (filled with defaults when possible): "
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
