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
INPUT_BASE  = "/eos/home-p/pelai/HZa/parquet_DNA/NanoV12/run3/"
OUTPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/"

# ===== Switches =====
DO_SIGNAL_NOMINAL = True
DO_SIGNAL_SYST    = True       # 先 False，稳定后再 True
DO_BKG_NOMINAL    = True
DO_BKG_SYST       = False       # 如后续需要也可 True
DO_DATA_NOMINAL   = True       # 你目前注释掉 Data，可按需改 True

# ===== Common config =====
years_all = {
    "Run2":      ["2016preVFP","2016postVFP","2017","2018"],
    "Run3_2022": ["2022preEE","2022postEE"],
    "Run3_2023": ["2023preBPix","2023postBPix"],
}

# 你之前实际使用的组合：
years_sig  = ["2022preEE","2022postEE","2023preBPix","2023postBPix", "2024"]
years_22   = ["2022preEE", "2022postEE"]  # 背景（2022）
years_23   = ["2023preBPix","2023postBPix", "2024"]  # 背景（2023）
years_dyll = ["2022preEE","2022postEE","2023preBPix","2023postBPix", "2024"]

# 系统误差类型
systs = ["FNUF","Material","Electron_scale","Electron_smear",
         "Muon_scale","Muon_smear","Photon_scale","Photon_smear"]
updown = ["up","down"]

# 扫描的 ma 列表（背景也按你的原脚本扫）
# ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30] # MC
ma_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30] # Bkg MC / Data

# ===== Samples =====
sig_samples = [f"mA_M{m}" for m in [1,2,3,4,5,6,7,8,9,10,15,20,25,30]]

# 背景
bkg_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
bkg_2023 = ["DYGto2LG_10to100"]
bkg_dyll = ["DYJetsToLL"]
# 如需其它：["ZGToLLG","WGToLNuG","ZG2JToG2L2J","EWKZ2J","TT","TTGJets","TGJets","ttWJets","ttZJets","WW","WZ","ZZ"]

# Data（如需）
data_samples = ["Data"]  # 输出路径合并为 /Data/<year>.root

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
        for s in sig_samples:
            s_ma = mass_from_sig(s)
            for y in years_sig:
                inp = os.path.join(INPUT_BASE, "Sig_MC", f"{s}_{y}", "merged_nominal.parquet")
                out = os.path.join(OUTPUT_BASE, s, f"{y}.root")
                f.write(line(inp, out, s_ma, "nominal", 1))

    # ---------- SIGNAL: systematics ----------
    if DO_SIGNAL_SYST:
        for s in sig_samples:
            s_ma = mass_from_sig(s)
            for y in years_sig:
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
        # 2022
        for s in bkg_2022:
            for y in years_22:
                for ma in ma_list:
                    inp = os.path.join(INPUT_BASE, "Bkg_MC", f"{s}_{y}", "merged_nominal.parquet")
                    out_dir = os.path.join(OUTPUT_BASE, s, f"ALP_M{ma}")
                    os.makedirs(out_dir, exist_ok=True)
                    out = os.path.join(out_dir, f"{y}.root")
                    f.write(line(inp, out, str(ma), "nominal", 0))
        # 2023
        for s in bkg_2023:
            for y in years_23:
                for ma in ma_list:
                    inp = os.path.join(INPUT_BASE, "Bkg_MC", f"{s}_{y}", "merged_nominal.parquet")
                    out_dir = os.path.join(OUTPUT_BASE, s, f"ALP_M{ma}")
                    os.makedirs(out_dir, exist_ok=True)
                    out = os.path.join(out_dir, f"{y}.root")
                    f.write(line(inp, out, str(ma), "nominal", 0))
        # DYJetsToLL
        for s in bkg_dyll:
            for y in years_dyll:
                for ma in ma_list:
                    inp = os.path.join(INPUT_BASE, "Bkg_MC", f"{s}_{y}", "merged_nominal.parquet")
                    out_dir = os.path.join(OUTPUT_BASE, s, f"ALP_M{ma}")
                    os.makedirs(out_dir, exist_ok=True)
                    out = os.path.join(out_dir, f"{y}.root")
                    f.write(line(inp, out, str(ma), "nominal", 0))

    # ---------- BACKGROUND: systematics（如需要） ----------
    if DO_BKG_SYST:
        for s in (bkg_2022 + bkg_2023 + bkg_dyll):
            # 用哪些年？按样本归属做简单划分：
            if s in bkg_2022:
                yrs = years_22
            elif s in bkg_2023:
                yrs = years_23
            else:
                yrs = years_dyll
            for y in yrs:
                for ud in updown:
                    for syst in systs:
                        corr = f"{syst}_{ud}"
                        for ma in ma_list:
                            inp = os.path.join(INPUT_BASE, "Bkg_MC", f"{s}_{y}", f"merged_{corr}.parquet")
                            out_dir = os.path.join(OUTPUT_BASE, f"{s}_{syst}_{ud}", f"ALP_M{ma}")
                            os.makedirs(out_dir, exist_ok=True)
                            out = os.path.join(out_dir, f"{y}.root")
                            f.write(line(inp, out, str(ma), corr, 0))

    # ---------- DATA（如需要） ----------
    if DO_DATA_NOMINAL:
        for y in years_dyll:  # 或按你的实际年份集合
            for ma in ma_list:
                inp = os.path.join(INPUT_BASE, "Data", f"Data_{y}", "merged_nominal.parquet")
                out_dir = os.path.join(OUTPUT_BASE, "Data", f"ALP_M{ma}")
                os.makedirs(out_dir, exist_ok=True)
                out = os.path.join(out_dir, f"{y}.root")
                # 注意：这里第三列传 str(ma)，SPLIT_FLAG=0
                f.write(line(inp, out, str(ma), "nominal", 0))


print("Wrote joblist.tsv with tasks.")
