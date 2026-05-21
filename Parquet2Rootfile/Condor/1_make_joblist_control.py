#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

# ===== Base paths =====
INPUT_BASE = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_control"
OUTPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_control"

# Control parquet merges are stored in dedicated *_control folders.
SIGNAL_INPUT_DIR = "Sig_MC_control"
BKG_INPUT_DIR = "Bkg_MC_control"
DATA_INPUT_DIR = "Data_control"

# ===== Switches =====
DO_SIGNAL_NOMINAL = True
DO_SIGNAL_SYST = True
DO_BKG_NOMINAL = True
DO_DATA_NOMINAL = True

DO_BKG_SYST = False

# Systematic variations
systs = [
    "FNUF",
    "Material",
    "Electron_scale",
    "Electron_smear",
    "Muon_scale",
    "Muon_smear",
    "Photon_scale",
    "Photon_smear",
]
updown = ["up", "down"]

# Year of Signal
year_sig_2022 = ["2022preEE", "2022postEE"]
year_sig_2023 = ["2023preBPix", "2023postBPix"]
year_sig_2024 = ["2024"]
# Year of Bkg
year_DYG_2022 = ["2022preEE", "2022postEE"]
year_DYG_2023 = ["2023preBPix", "2023postBPix"]
year_DYG_2024 = ["2024"]
years_DYJet_2022 = ["2022preEE", "2022postEE"]
year_DYJet_2023 = ["2023preBPix", "2023postBPix"]
years_DYJet_2024 = ["2024"]
# Year of Data
years_Data_2022 = ["2022preEE", "2022postEE"]
year_Data_2023 = ["2023preBPix", "2023postBPix"]
years_Data_2024 = ["2024"]

# Name of Signal Sample
name_sig_2022 = ["mA_M1", "mA_M2", "mA_M3", "mA_M4", "mA_M5", "mA_M6", "mA_M7", "mA_M8", "mA_M9", "mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
name_sig_2023 = ["mA_M1", "mA_M2", "mA_M3", "mA_M4", "mA_M5", "mA_M6", "mA_M7", "mA_M8", "mA_M9", "mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
name_sig_2024 = ["mA_M1", "mA_M2", "mA_M3", "mA_M4", "mA_M5", "mA_M6", "mA_M7", "mA_M8", "mA_M9", "mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
# Name of Bkg Sample
name_DYG_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
name_DYG_2023 = ["DYGto2LG_10to100"]
name_DYG_2024 = ["DYGto2LG_10to100"]
name_DYJet_2022 = ["DYJetsToLL"]
name_DYJet_2023 = ["DYJetsToLL"]
name_DYJet_2024 = ["DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
# Name of Data Sample
name_Data_2022 = ["Data"]
name_Data_2023 = ["Data"]
name_Data_2024 = ["Data"]

SIGNAL_BY_YEAR = {
    **{y: name_sig_2022 for y in year_sig_2022},
    **{y: name_sig_2023 for y in year_sig_2023},
    **{y: name_sig_2024 for y in year_sig_2024},
}

DYG_BY_YEAR = {
    **{y: name_DYG_2022 for y in year_DYG_2022},
    **{y: name_DYG_2023 for y in year_DYG_2023},
    **{y: name_DYG_2024 for y in year_DYG_2024},
}

DYJET_BY_YEAR = {
    **{y: name_DYJet_2022 for y in years_DYJet_2022},
    **{y: name_DYJet_2023 for y in year_DYJet_2023},
    **{y: name_DYJet_2024 for y in years_DYJet_2024},
}

DATA_BY_YEAR = {
    **{y: name_Data_2022 for y in years_Data_2022},
    **{y: name_Data_2023 for y in year_Data_2023},
    **{y: name_Data_2024 for y in years_Data_2024},
}


def iter_year_sample(mapping):
    for year, samples in mapping.items():
        for sample in samples:
            yield year, sample


def line(input_pq, output_root, corr, split_flag):
    return f"{input_pq}\t{output_root}\t{corr}\t{split_flag}\n"


os.makedirs("logs_control", exist_ok=True)

with open("joblist_control.tsv", "w") as f:
    # ---------- SIGNAL: nominal ----------
    if DO_SIGNAL_NOMINAL:
        for y, s in iter_year_sample(SIGNAL_BY_YEAR):
            inp = os.path.join(INPUT_BASE, SIGNAL_INPUT_DIR, f"{s}_{y}", "merged_nominal.parquet")
            out_dir = os.path.join(OUTPUT_BASE, s)
            os.makedirs(out_dir, exist_ok=True)
            out = os.path.join(out_dir, f"{y}.root")
            f.write(line(inp, out, "nominal", 1))

    # ---------- SIGNAL: systematics ----------
    if DO_SIGNAL_SYST:
        for y, s in iter_year_sample(SIGNAL_BY_YEAR):
            for ud in updown:
                for syst in systs:
                    corr = f"{syst}_{ud}"
                    inp = os.path.join(INPUT_BASE, SIGNAL_INPUT_DIR, f"{s}_{y}", f"merged_{corr}.parquet")
                    out_dir = os.path.join(OUTPUT_BASE, f"{s}_{syst}_{ud}")
                    os.makedirs(out_dir, exist_ok=True)
                    out = os.path.join(out_dir, f"{y}.root")
                    f.write(line(inp, out, corr, 1))

    # ---------- BACKGROUND: nominal ----------
    if DO_BKG_NOMINAL:
        for y, s in iter_year_sample(DYG_BY_YEAR):
            inp = os.path.join(INPUT_BASE, BKG_INPUT_DIR, f"{s}_{y}", "merged_nominal.parquet")
            out_dir = os.path.join(OUTPUT_BASE, s)
            os.makedirs(out_dir, exist_ok=True)
            out = os.path.join(out_dir, f"{y}.root")
            f.write(line(inp, out, "nominal", 0))

        for y, s in iter_year_sample(DYJET_BY_YEAR):
            inp = os.path.join(INPUT_BASE, BKG_INPUT_DIR, f"{s}_{y}", "merged_nominal.parquet")
            out_dir = os.path.join(OUTPUT_BASE, s)
            os.makedirs(out_dir, exist_ok=True)
            out = os.path.join(out_dir, f"{y}.root")
            f.write(line(inp, out, "nominal", 0))

    # ---------- BACKGROUND: systematics ----------
    if DO_BKG_SYST:
        for y, s in iter_year_sample(DYG_BY_YEAR):
            for ud in updown:
                for syst in systs:
                    corr = f"{syst}_{ud}"
                    inp = os.path.join(INPUT_BASE, BKG_INPUT_DIR, f"{s}_{y}", f"merged_{corr}.parquet")
                    out_dir = os.path.join(OUTPUT_BASE, f"{s}_{syst}_{ud}")
                    os.makedirs(out_dir, exist_ok=True)
                    out = os.path.join(out_dir, f"{y}.root")
                    f.write(line(inp, out, corr, 0))

        for y, s in iter_year_sample(DYJET_BY_YEAR):
            for ud in updown:
                for syst in systs:
                    corr = f"{syst}_{ud}"
                    inp = os.path.join(INPUT_BASE, BKG_INPUT_DIR, f"{s}_{y}", f"merged_{corr}.parquet")
                    out_dir = os.path.join(OUTPUT_BASE, f"{s}_{syst}_{ud}")
                    os.makedirs(out_dir, exist_ok=True)
                    out = os.path.join(out_dir, f"{y}.root")
                    f.write(line(inp, out, corr, 0))

    # ---------- DATA ----------
    if DO_DATA_NOMINAL:
        for y, s in iter_year_sample(DATA_BY_YEAR):
            inp = os.path.join(INPUT_BASE, DATA_INPUT_DIR, f"{s}_{y}", "merged_nominal.parquet")
            out_dir = os.path.join(OUTPUT_BASE, s)
            os.makedirs(out_dir, exist_ok=True)
            out = os.path.join(out_dir, f"{y}.root")
            f.write(line(inp, out, "nominal", 0))

print("Wrote joblist_control.tsv with tasks.")
