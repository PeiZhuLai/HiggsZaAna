#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate job list TSV for HTCondor queue-from.
每行一条作业，字段：
INPUT \t OUTPUT \t MA \t CORR \t SPLIT_FLAG

SPLIT_FLAG:
  - 信号样本: 1（使用 --split，并且 --ma=样本名中的质量 M）
  - 背景样本/数据: 0（背景仍然会扫 ma；如不需要可把 ma_list 设为只含一个占位 NA）

CORR:
  - "nominal" 或 "{syst}_{up|down}"

按你的实际需求增删样本/年份/系统误差即可。
"""

import os

# ===== Base paths =====
INPUT_BASE  = "/eos/home-p/pelai/HZa/parquet_DNA"
OUTPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_BDT"

# ===== Switches =====
DO_SIGNAL_NOMINAL = False
DO_SIGNAL_SYST    = True       # 先 False，稳定后再 True
DO_BKG_NOMINAL    = False
DO_BKG_SYST       = False       # 如后续需要也可 True
DO_DATA_NOMINAL   = False       # 你目前注释掉 Data，可按需改 True


# Year of Signal
year_sig_2022 = ["2022preEE", "2022postEE"]
year_sig_2023 = ["2023preBPix", "2023postBPix"]
year_sig_2024 = ["2024"]
# Year of Bkg
year_DYG_2022 = ["2022preEE", "2022postEE"]
year_DYG_2023 = ["2023preBPix", "2023postBPix"]
year_DYG_2024 = ["2024"]
years_DYJet_2022 = ["2022preEE","2022postEE"]
year_DYJet_2023  = ["2023preBPix", "2023postBPix"]
years_DYJet_2024 = ["2024"]
# Year of Data
years_Data_2022 = ["2022preEE","2022postEE"]
year_Data_2023  = ["2023preBPix", "2023postBPix"]
years_Data_2024 = ["2024"]

# Name of Signal Sample
name_sig_2022 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
name_sig_2023 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
name_sig_2024 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
# Name of Bkg Sample
name_DYG_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
name_DYG_2023 = ["DYGto2LG_10to100"]
name_DYG_2024 = ["DYGto2LG_10to100"]
name_DYJet_2022 = ["DYJetsToLL"]
name_DYJet_2023 = ["DYJetsToLL"]
name_DYJet_2024 = ["DYJetsTo2E","DYJetsTo2Mu","DYJetsTo2Tau"]
# Name of Data Sample
name_Data_2022 = ["Data"]
name_Data_2023 = ["Data"]
name_Data_2024 = ["Data"]

ma_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30] # Bkg MC / Data

# 系统误差类型
# systs = ["FNUF","Material","Electron_scale","Electron_smear","Muon_scale","Muon_smear","Photon_scale","Photon_smear"]
systs = ["FNUF","Material","Electron_scale","Electron_smear"]
# systs = ["Muon_scale","Muon_smear","Photon_scale","Photon_smear"]
updown = ["up","down"]

# 以你提供的清單建立「year -> samples」mapping（統一驅動 loop）
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
    # mapping: dict[str, list[str]]  -> yields (year, sample)
    for y, samples in mapping.items():
        for s in samples:
            yield y, s

# ===== Helpers =====
def mass_from_sig(sample: str) -> str:
    # ALP_M5 -> "5"
    if "_M" in sample:
        return sample.split("_M", 1)[1]
    return "NA"

def line(input_pq, output_root, ma, corr, split_flag):
    return f"{input_pq}\t{output_root}\t{ma}\t{corr}\t{split_flag}\n"

# ===== Build job list =====
os.makedirs("logs", exist_ok=True)

with open("joblist.tsv", "w") as f:
    # ---------- SIGNAL: nominal ----------
    if DO_SIGNAL_NOMINAL:
        for y, s in iter_year_sample(SIGNAL_BY_YEAR):
            s_ma = mass_from_sig(s)
            inp = os.path.join(INPUT_BASE, "Sig_MC", f"{s}_{y}", "merged_nominal.parquet")
            out_dir = os.path.join(OUTPUT_BASE, s)
            os.makedirs(out_dir, exist_ok=True)
            out = os.path.join(out_dir, f"{y}.root")
            f.write(line(inp, out, s_ma, "nominal", 1))

    # ---------- SIGNAL: systematics ----------
    if DO_SIGNAL_SYST:
        for y, s in iter_year_sample(SIGNAL_BY_YEAR):
            s_ma = mass_from_sig(s)
            for ud in updown:
                for syst in systs:
                    corr = f"{syst}_{ud}"
                    inp = os.path.join(INPUT_BASE, "Sig_MC", f"{s}_{y}", f"merged_{corr}.parquet")
                    out_dir = os.path.join(OUTPUT_BASE, f"{s}_{syst}_{ud}")
                    os.makedirs(out_dir, exist_ok=True)
                    out = os.path.join(out_dir, f"{y}.root")
                    f.write(line(inp, out, s_ma, corr, 1))

    # ---------- BACKGROUND: nominal ----------
    if DO_BKG_NOMINAL:
        # DYG
        for y, s in iter_year_sample(DYG_BY_YEAR):
            for ma in ma_list:
                inp = os.path.join(INPUT_BASE, "Bkg_MC", f"{s}_{y}", "merged_nominal.parquet")
                out_dir = os.path.join(OUTPUT_BASE, s, f"mA_M{ma}")
                os.makedirs(out_dir, exist_ok=True)
                out = os.path.join(out_dir, f"{y}.root")
                f.write(line(inp, out, str(ma), "nominal", 0))

        # DYJets
        for y, s in iter_year_sample(DYJET_BY_YEAR):
            for ma in ma_list:
                inp = os.path.join(INPUT_BASE, "Bkg_MC", f"{s}_{y}", "merged_nominal.parquet")
                out_dir = os.path.join(OUTPUT_BASE, s, f"mA_M{ma}")
                os.makedirs(out_dir, exist_ok=True)
                out = os.path.join(out_dir, f"{y}.root")
                f.write(line(inp, out, str(ma), "nominal", 0))

    # ---------- BACKGROUND: systematics（如需要） ----------
    if DO_BKG_SYST:
        # DYG
        for y, s in iter_year_sample(DYG_BY_YEAR):
            for ud in updown:
                for syst in systs:
                    corr = f"{syst}_{ud}"
                    for ma in ma_list:
                        inp = os.path.join(INPUT_BASE, "Bkg_MC", f"{s}_{y}", f"merged_{corr}.parquet")
                        out_dir = os.path.join(OUTPUT_BASE, f"{s}_{syst}_{ud}", f"mA_M{ma}")
                        os.makedirs(out_dir, exist_ok=True)
                        out = os.path.join(out_dir, f"{y}.root")
                        f.write(line(inp, out, str(ma), corr, 0))

        # DYJets
        for y, s in iter_year_sample(DYJET_BY_YEAR):
            for ud in updown:
                for syst in systs:
                    corr = f"{syst}_{ud}"
                    for ma in ma_list:
                        inp = os.path.join(INPUT_BASE, "Bkg_MC", f"{s}_{y}", f"merged_{corr}.parquet")
                        out_dir = os.path.join(OUTPUT_BASE, f"{s}_{syst}_{ud}", f"mA_M{ma}")
                        os.makedirs(out_dir, exist_ok=True)
                        out = os.path.join(out_dir, f"{y}.root")
                        f.write(line(inp, out, str(ma), corr, 0))

    # ---------- DATA（如需要） ----------
    if DO_DATA_NOMINAL:
        # 你目前資料樣本名其實都叫 Data，但仍保留 mapping 以便未來擴充
        for y, s in iter_year_sample(DATA_BY_YEAR):
            for ma in ma_list:
                inp = os.path.join(INPUT_BASE, "Data", f"{s}_{y}", "merged_nominal.parquet")
                out_dir = os.path.join(OUTPUT_BASE, s, f"mA_M{ma}")
                os.makedirs(out_dir, exist_ok=True)
                out = os.path.join(out_dir, f"{y}.root")
                f.write(line(inp, out, str(ma), "nominal", 0))

print("Wrote joblist.tsv with tasks.")
