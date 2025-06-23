#!/usr/bin/env python3
"""
Cut‑flow aggregator for Run‑3 background MC samples.

Features (compared with original script):
─────────────────────────────────────────
1.  **year_to_datasets mapping**
    – One place to configure which datasets should be summed for every year.
2.  **Outer loop on `year`, inner loop on `dataset`**
    – Guarantees that yields are accumulated across all datasets *within* the same year.
3.  **`yields_dict` reset only once per year**
    – Prevents accidental overwriting; keeps a clean per‑year summary.
4.  **Robust parquet loading & weight handling**
    – Graceful fallback to weight = 1 when file or columns are missing.
5.  **Cleaner log parsing**
    – Regex compiled once; string clean‑up consolidated.
6.  **Structured file output**
    – Appends yearly sections to `output_file`.
7.  **Timing & memory hygiene**
    – `del data` after each parquet to free RAM; prints total runtime.

Author: ChatGPT (modified for pelai)
Date  : 2025‑06‑18
"""

import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ────────────────────────────────
# User‑Configurable Paths & Flags
# ────────────────────────────────
EOS_PATH   = "/eos/home-p/pelai/HZa/Parquet/NanoV12/run3/"
LOG_PATH   = "/afs/cern.ch/work/p/pelai/HZa/Output/cutflow_outfile/run3/"
DATASET_TYP = "Bkg_MC"
OUTPUT_FILE = "cutflow_list/cut_yields_output_DYG_Bkg_MC_run3.txt"

# Which datasets to aggregate for each year
YEAR_TO_DATASETS = {
    "2022preEE"  : ["DYGto2LG_10to50", "DYGto2LG_50to100"],
    "2022postEE" : ["DYGto2LG_10to50", "DYGto2LG_50to100"],
    "2023preBPix": ["DYGto2LG_10to100"],
    "2023postBPix":["DYGto2LG_10to100"],
}

# ────────────────────────────────
# Cut‑flow Meta‑Information
# ────────────────────────────────
CUTFLOW_TYPE = [
    "zgammas", "zgammas_ele", "zgammas_mu",
    "zgammas_w", "zgammas_ele_w", "zgammas_mu_w",
]

CUT_TYPE = [
    "all", "N_lep_sel", "trig_cut", "lep_pt_cut", "has_z_cand",
    "has_2g_cand", "sel_h_1", "sel_h_2", "event",
]

CUT_NAME = {
    "zgammas": [
        "Initial Events", r"$N_{l} \geq 2$",
        r"e, ee trigger || $\mu\mu$, $\mu$ triggers ",
        r"lepton trigger pT cut",
        r"$N_{\gamma} \geq 2$",
        r"$80 \\text{ GeV} < m_{ll} < 100 \\text{ GeV}$",
        r"$m_{ll} + m_{ll\gamma\gamma} > 185 \\text{ GeV}$",
        r"$95 \\text{ GeV} < m_{ll\gamma\gamma} < 180 \\text{ GeV}$",
        "Event Filtering", "Total Baseline Events",
    ],
    "zgammas_ele": [
        "Initial Events", r"$N_{\\textrm{e}} \geq 2$", "ee, e triggers ",
        r"ee, e trigger pT cut", r"$N_{\\gamma} \geq 2$",
        r"$80 \\text{ GeV} < m_{\\textrm{ee}} < 100 \\text{ GeV}$",
        r"$m_{\\textrm{ee}} + m_{\\textrm{ee}\\gamma\\gamma} > 185 \\text{ GeV}$",
        r"$95 \\text{ GeV} < m_{\\textrm{ee}\\gamma\\gamma} < 180 \\text{ GeV}$",
        "Event Filtering", "Total Baseline Events",
    ],
    "zgammas_mu": [
        "Initial Events", r"$N_{\mu} \geq 2$", r"$\mu\mu$, $\mu$ triggers ",
        r"$\mu\mu$, $\mu$ trigger pT cut", r"$N_{\gamma} \geq 2$",
        r"$80 \\text{ GeV} < m_{\mu\mu} < 100 \\text{ GeV}$",
        r"$m_{\mu\mu} + m_{\mu\mu\gamma\gamma} > 185 \\text{ GeV}$",
        r"$95 \\text{ GeV} < m_{\mu\mu\gamma\gamma} < 180 \\text{ GeV}$",
        "Event Filtering", "Total Baseline Events",
    ],
    "zgammas_w": [
        "Initial Events", r"$N_{l} \geq 2$", r"e, ee trigger || $\mu\mu$, $\mu$ triggers ",
        r"lepton trigger pT cut", r"$N_{\gamma} \geq 2$",
        r"$80 \\text{ GeV} < m_{ll} < 100 \\text{ GeV}$",
        r"$m_{ll} + m_{ll\gamma\gamma} > 185 \\text{ GeV}$",
        r"$95 \\text{ GeV} < m_{ll\gamma\gamma} < 180 \\text{ GeV}$",
        "Event Filtering", "Total Baseline Events",
    ],
    "zgammas_ele_w": [
        "Initial Events", r"$N_{\\textrm{e}} \geq 2$", "ee, e triggers ",
        r"ee, e trigger pT cut", r"$N_{\\gamma} \geq 2$",
        r"$80 \\text{ GeV} < m_{\\textrm{ee}} < 100 \\text{ GeV}$",
        r"$m_{\\textrm{ee}} + m_{\\textrm{ee}\\gamma\\gamma} > 185 \\text{ GeV}$",
        r"$95 \\text{ GeV} < m_{\\textrm{ee}\\gamma\\gamma} < 180 \\text{ GeV}$",
        "Event Filtering", "Total Baseline Events",
    ],
    "zgammas_mu_w": [
        "Initial Events", r"$N_{\mu} \geq 2$", r"$\mu\mu$, $\mu$ triggers ",
        r"$\mu\mu$, $\mu$ trigger pT cut", r"$N_{\gamma} \geq 2$",
        r"$80 \\text{ GeV} < m_{\mu\mu} < 100 \\text{ GeV}$",
        r"$m_{\mu\mu} + m_{\mu\mu\gamma\gamma} > 185 \\text{ GeV}$",
        r"$95 \\text{ GeV} < m_{\mu\mu\gamma\gamma} < 180 \\text{ GeV}$",
        "Event Filtering", "Total Baseline Events",
    ],
}

CUT_NUM = len(CUT_TYPE)

# Template for zeroed arrays (if someday you want full cut‑flow arrays)
CUTFLOW_TEMPLATE = {ct: np.zeros(CUT_NUM, dtype=float) for ct in CUTFLOW_TYPE}  # noqa: E501

# ────────────────────────────────
# Helpers
# ────────────────────────────────
TAGGER_REGEX = re.compile(
    r"\[Tagger\]\s*:.*?cut type\s*:\s*(\w+),\s*cut\s*:\s*"  # Cut‑type
    r"([\w\s|<>.=!'\\-]+?),?\s*yields\s*:?\s*([\d.]+)"        # Cut‑name & yields
)


def append_output(section: str) -> None:
    """Append ``section`` to OUTPUT_FILE, creating folders if needed."""
    Path(OUTPUT_FILE).parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "a", encoding="utf-8") as fh:
        fh.write(section + "\n")


# ────────────────────────────────
# Main Driver
# ────────────────────────────────

if __name__ == "__main__":
    t0 = time.time()

    append_output("📌📌📌 Start of cut‑flow summary 📌📌📌")

    for year, datasets in YEAR_TO_DATASETS.items():
        print(f"\n▶ Year {year}: processing {len(datasets)} dataset(s)…")

        # One dict per year, accumulates yields across all datasets
        yields_dict: dict[str, dict[str, float]] = {}

        for dataset in datasets:
            parquet_file = f"{EOS_PATH}{DATASET_TYP}/{dataset}_{year}/merged_nominal.parquet"
            weight = 1.0
            try:
                df = pd.read_parquet(parquet_file)
                if {"weight_central", "weight_central_initial"}.issubset(df.columns):
                    weight = (
                        # 乘 2.0 是因為沒有 skimmed root files ，填上ㄧ樣的 root file 導致 sum of weights 變兩倍，scale1fb 小了兩倍
                        2.0* df["weight_central"].to_numpy(dtype="float64").sum()
                        / df["weight_central_initial"].to_numpy(dtype="float64").sum()
                    )
                del df  # free memory
            except Exception as exc:
                print(f"  [WARN] {dataset}_{year}: cannot read parquet – {exc!s}; using weight = 1")

            # -------- Parse *.out logs --------
            log_dataset_dir = Path(LOG_PATH) / DATASET_TYP / f"{dataset}_{year}"
            if not log_dataset_dir.is_dir():
                print(f"  [WARN] {log_dataset_dir} does not exist – skipping")
                continue

            for subdir in log_dataset_dir.iterdir():
                if not subdir.is_dir():
                    continue

                # Choose the largest .out file (by numeric prefix)
                out_files = [p for p in subdir.glob("*.out")]
                if not out_files:
                    continue
                try:
                    largest = max(out_files, key=lambda p: float(p.stem))
                except ValueError:
                    # non‑numeric names – fall back to latest mtime
                    largest = max(out_files, key=lambda p: p.stat().st_mtime)

                with open(largest, "r", encoding="utf-8", errors="ignore") as fh:
                    segments = fh.read().split("DEBUG")

                if len(segments) < CUT_NUM * 3 + 40:
                    continue  # suspiciously short – skip

                for seg in segments:
                    seg = re.sub(r"tagger\.py:\d+", "", seg)
                    seg = re.sub(r"analysis\.py:\d+", "", seg)
                    seg = " ".join(seg.split())  # normalise spaces

                    if "[Tagger]" not in seg or "WARNING" in seg or "[[INFO]]" in seg:
                        continue

                    m = TAGGER_REGEX.match(seg)
                    if not m:
                        continue
                    cut_type, cut_label, yields = m.group(1).strip(), m.group(2).strip(), float(m.group(3))

                    # init sub‑dict if needed
                    yields_dict.setdefault(cut_type, {})

                    scale = weight if "_w" in cut_type else 1.0  # weight only for *_w types
                    yields_dict[cut_type][cut_label] = yields_dict[cut_type].get(cut_label, 0.0) + yields * scale

        # ───── Finished all datasets for the year ─────
        print(f"   ↪ finished aggregation for {year}")

        # ⇣ Write human‑readable summary to file ⇣
        lines = [f"\n📌 Year {year}"]
        for ct, cuts in yields_dict.items():
            lines.append(f"Cut Type: {ct}")
            if ct in CUT_NAME:
                for idx, (cut, yield_val) in enumerate(cuts.items()):
                    pretty = CUT_NAME[ct][idx] if idx < len(CUT_NAME[ct]) else cut
                    lines.append(f"     {pretty:30} & {yield_val:.0f} \\" )
            else:
                for cut, yield_val in cuts.items():
                    lines.append(f"     {cut:30} & {yield_val:.0f} \\" )
        append_output("\n".join(lines))

    # ────────────────────────────────
    # Done
    # ────────────────────────────────
    elapsed = time.time() - t0
    print(f"\n✔ All done – total runtime {elapsed:.1f} s")
