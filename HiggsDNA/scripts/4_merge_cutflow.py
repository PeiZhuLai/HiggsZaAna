#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge per-sample cutflow JSON files into analysis-note LaTeX tables.

Input files are expected to come from 3_collect_cutflow.py:
  cutflow_<dataset_type>_<dataset>_<era>.json

Outputs:
  - Signal: two tables per era, mA=1-7 and mA=8-30.
  - Background: one table per era.
  - Data: one table per era.

Each sample has two columns: event count and N-1 efficiency.
Counts are rounded to integers; efficiencies are shown as one decimal percent.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


DEFAULT_INPUT_DIR = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list")
DEFAULT_OUTPUT_SUBDIR = "merged_cutflow_latex"

YEARS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
DATASET_TYPES = ["Sig_MC", "Bkg_MC", "Data"]

SIGNAL_MASS_GROUPS = [
    ("mA1to7", [1, 2, 3, 4, 5, 6, 7], r"$m_{a}=1$--$7$"),
    ("mA8to30", [8, 9, 10, 15, 20, 25, 30], r"$m_{a}=8$--$30$"),
]

CUT_ORDER = [
    "all",
    "N_lep_sel",
    "trig_cut",
    "lep_pt_cut",
    "ele_sip3d_cut",
    "has_z_cand",
    "g_kin_cut",
    "has_2g_cand",
    "sel_h_1",
    "sel_h_2",
    "sel_h",
    "event",
    "all cuts",
]

LABELS_COMMON = {
    "all": "Initial Events",
    "N_lep_sel": r"$N_{l}\geq 2$",
    "trig_cut": "Trigger",
    "lep_pt_cut": r"Lepton Trigger $p_T$ Cut",
    "ele_sip3d_cut": "Electron SIP3D Cut",
    "has_z_cand": r"$m_{ll} > 50\,\mathrm{GeV}$",
    "g_kin_cut": "Photon Kinematic Cuts",
    "has_2g_cand": r"$N_{\gamma}\geq 2$",
    "sel_h_1": r"$m_{ll} + m_{ll\gamma\gamma} > 185\,\mathrm{GeV}$",
    "sel_h_2": r"$95\,\mathrm{GeV} < m_{ll\gamma\gamma} < 180\,\mathrm{GeV}$",
    "sel_h": r"$95\,\mathrm{GeV} < m_{ll\gamma\gamma} < 180\,\mathrm{GeV}$",
    "event": "Event Filtering",
    "all cuts": "Total Baseline Events",
}


@dataclass(frozen=True)
class SampleCutflow:
    dataset_type: str
    dataset: str
    year: str
    path: Path
    cutflows: Dict[str, Dict[str, float]]


@dataclass(frozen=True)
class TableColumn:
    key: str
    header: str
    cut_type: str
    cuts: Dict[str, float]


_LATEX_ESCAPE_MAP = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}


def latex_escape_text(s: str) -> str:
    return "".join(_LATEX_ESCAPE_MAP.get(ch, ch) for ch in str(s))


def latex_texttt(s: str) -> str:
    return rf"\texttt{{{latex_escape_text(s)}}}"


def latex_label_id(s: str) -> str:
    label = re.sub(r"[^A-Za-z0-9]+", "-", str(s)).strip("-").lower()
    return label or "table"


def label_for(cut_type: str, cut_key: str) -> str:
    if cut_key == "N_lep_sel":
        if "_ele" in cut_type:
            return r"$N_{e}\geq 2$"
        if "_mu" in cut_type:
            return r"$N_{\mu}\geq 2$"
    if cut_key == "has_z_cand":
        if "_ele" in cut_type:
            return r"$m_{ee} > 50\,\mathrm{GeV}$"
        if "_mu" in cut_type:
            return r"$m_{\mu\mu} > 50\,\mathrm{GeV}$"
    return LABELS_COMMON.get(cut_key, latex_escape_text(cut_key))


def parse_cutflow_filename(path: Path) -> Optional[Tuple[str, str, str]]:
    stem = path.stem
    if not stem.startswith("cutflow_"):
        return None

    for dataset_type in DATASET_TYPES:
        prefix = f"cutflow_{dataset_type}_"
        if not stem.startswith(prefix):
            continue

        rest = stem[len(prefix):]
        for year in YEARS:
            suffix = f"_{year}"
            if rest.endswith(suffix):
                dataset = rest[:-len(suffix)]
                if dataset:
                    return dataset_type, dataset, year
    return None


def load_samples(input_dir: Path) -> List[SampleCutflow]:
    samples: List[SampleCutflow] = []
    for path in sorted(input_dir.glob("cutflow_*.json")):
        parsed = parse_cutflow_filename(path)
        if parsed is None:
            continue

        dataset_type, dataset, year = parsed
        try:
            with path.open("r", encoding="utf-8") as f:
                payload = json.load(f)
        except Exception as exc:
            print(f"[SKIP] cannot read {path}: {exc}")
            continue

        cutflows = payload.get("cutflows")
        if not isinstance(cutflows, dict):
            print(f"[SKIP] missing cutflows block: {path}")
            continue

        samples.append(
            SampleCutflow(
                dataset_type=dataset_type,
                dataset=dataset,
                year=year,
                path=path,
                cutflows=cutflows,
            )
        )
    return samples


def default_input_dir() -> Path:
    if DEFAULT_INPUT_DIR.exists() and any(DEFAULT_INPUT_DIR.glob("cutflow_*.json")):
        return DEFAULT_INPUT_DIR
    return Path(__file__).resolve().parent


def mass_from_dataset(dataset: str) -> Optional[int]:
    match = re.fullmatch(r"mA_M(\d+)", dataset)
    return int(match.group(1)) if match else None


def choose_cutflow(sample: SampleCutflow, preferred_cut_type: str) -> Optional[TableColumn]:
    cut_type = preferred_cut_type
    if cut_type not in sample.cutflows:
        if cut_type.endswith("_w") and cut_type[:-2] in sample.cutflows:
            cut_type = cut_type[:-2]
        elif "zgammas" in sample.cutflows:
            cut_type = "zgammas"
        else:
            return None

    cuts = sample.cutflows.get(cut_type)
    if not isinstance(cuts, dict):
        return None

    if sample.dataset_type == "Sig_MC":
        mass = mass_from_dataset(sample.dataset)
        header = rf"$m_a={mass}$" if mass is not None else latex_texttt(sample.dataset)
    else:
        header = latex_texttt(sample.dataset)

    return TableColumn(
        key=f"{sample.dataset_type}_{sample.dataset}_{sample.year}",
        header=header,
        cut_type=cut_type,
        cuts={str(k): float(v) for k, v in cuts.items()},
    )


def ordered_cut_keys(columns: Sequence[TableColumn]) -> List[str]:
    keys: List[str] = []
    for cut_key in CUT_ORDER:
        if any(cut_key in col.cuts for col in columns):
            keys.append(cut_key)

    seen = set(keys)
    for col in columns:
        for cut_key in col.cuts:
            if cut_key not in seen:
                keys.append(cut_key)
                seen.add(cut_key)
    return keys


def n_minus_one_efficiencies(cuts: Dict[str, float], rows: Sequence[str]) -> Dict[str, Optional[float]]:
    out: Dict[str, Optional[float]] = {}
    previous_value: Optional[float] = None

    for cut_key in rows:
        if cut_key not in cuts:
            out[cut_key] = None
            continue

        value = float(cuts[cut_key])
        if previous_value is None:
            out[cut_key] = 100.0
        elif previous_value == 0.0:
            out[cut_key] = None
        else:
            out[cut_key] = 100.0 * value / previous_value

        previous_value = value

    return out


def overall_efficiency(cuts: Dict[str, float], rows: Sequence[str]) -> Optional[float]:
    initial_key = "all" if "all" in cuts else next((key for key in rows if key in cuts), None)
    final_key = "all cuts" if "all cuts" in cuts else next((key for key in reversed(rows) if key in cuts), None)
    if initial_key is None or final_key is None:
        return None

    initial = float(cuts[initial_key])
    if initial == 0.0:
        return None
    return 100.0 * float(cuts[final_key]) / initial


def format_count(value: float) -> str:
    return f"{value:.0f}"


def format_eff(value: Optional[float]) -> str:
    return "--" if value is None else f"{value:.1f}"


def render_table(*, caption: str, label: str, columns: Sequence[TableColumn]) -> str:
    rows = ordered_cut_keys(columns)
    col_spec = "l" + "rr" * len(columns)

    lines: List[str] = [
        r"% Requires \usepackage{graphicx} for \resizebox.",
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{{caption}}}",
        rf"\label{{tab:{latex_label_id(label)}}}",
        r"\resizebox{\textwidth}{!}{%",
        rf"\begin{{tabular}}{{{col_spec}}}",
        r"\hline",
    ]

    top_header = ["Selection"]
    for col in columns:
        top_header.append(rf"\multicolumn{{2}}{{c}}{{{col.header}}}")
    lines.append(" & ".join(top_header) + r" \\")

    sub_header = [""]
    for _ in columns:
        sub_header.extend(["Events", r"N-1 eff. (\%)"])
    lines.append(" & ".join(sub_header) + r" \\")
    lines.append(r"\hline")

    eff_by_column = [n_minus_one_efficiencies(col.cuts, rows) for col in columns]
    for cut_key in rows:
        row = [label_for(columns[0].cut_type, cut_key)]
        for col, effs in zip(columns, eff_by_column):
            if cut_key not in col.cuts:
                row.extend(["--", "--"])
                continue
            row.extend([format_count(col.cuts[cut_key]), format_eff(effs[cut_key])])
        lines.append(" & ".join(row) + r" \\")

    overall_row = ["Overall efficiency"]
    for col in columns:
        overall_row.extend(["--", format_eff(overall_efficiency(col.cuts, rows))])

    lines.extend([
        r"\hline",
        " & ".join(overall_row) + r" \\",
        r"\hline",
        r"\end{tabular}%",
        r"}",
        r"\end{table}",
        "",
    ])
    return "\n".join(lines)


def samples_for(samples: Sequence[SampleCutflow], dataset_type: str, year: str) -> List[SampleCutflow]:
    return [s for s in samples if s.dataset_type == dataset_type and s.year == year]


def signal_group_columns(
    samples: Sequence[SampleCutflow],
    year: str,
    masses: Sequence[int],
    cut_type: str,
) -> List[TableColumn]:
    by_mass = {
        mass_from_dataset(sample.dataset): sample
        for sample in samples_for(samples, "Sig_MC", year)
        if mass_from_dataset(sample.dataset) is not None
    }

    columns: List[TableColumn] = []
    for mass in masses:
        sample = by_mass.get(mass)
        if sample is None:
            continue
        col = choose_cutflow(sample, cut_type)
        if col is not None:
            columns.append(col)
    return columns


def background_columns(samples: Sequence[SampleCutflow], year: str, cut_type: str) -> List[TableColumn]:
    columns: List[TableColumn] = []
    for sample in sorted(samples_for(samples, "Bkg_MC", year), key=lambda s: s.dataset):
        col = choose_cutflow(sample, cut_type)
        if col is not None:
            columns.append(col)
    return columns


def data_columns(samples: Sequence[SampleCutflow], year: str, cut_type: str) -> List[TableColumn]:
    columns: List[TableColumn] = []
    for sample in sorted(samples_for(samples, "Data", year), key=lambda s: s.dataset):
        col = choose_cutflow(sample, cut_type)
        if col is not None:
            columns.append(col)
    return columns


def write_table(out_dir: Path, filename: str, table: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / filename
    path.write_text(table + "\n", encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge cutflow JSON files into LaTeX tables.")
    parser.add_argument("--input-dir", default=None, help="Directory containing cutflow_*.json files.")
    parser.add_argument("--output-dir", default=None, help="Directory for merged LaTeX tables.")
    parser.add_argument("--sig-cut-type", default="zgammas", help="Cutflow key for signal tables.")
    parser.add_argument("--bkg-cut-type", default="zgammas_w", help="Cutflow key for background tables.")
    parser.add_argument("--data-cut-type", default="zgammas", help="Cutflow key for data tables.")
    args = parser.parse_args()

    input_dir = Path(args.input_dir) if args.input_dir else default_input_dir()
    output_dir = Path(args.output_dir) if args.output_dir else input_dir / DEFAULT_OUTPUT_SUBDIR

    samples = load_samples(input_dir)
    if not samples:
        print(f"[ERROR] no cutflow_*.json files found in {input_dir}")
        return 1

    all_tables: List[str] = []
    written: List[Path] = []

    for year in YEARS:
        for group_name, masses, mass_caption in SIGNAL_MASS_GROUPS:
            columns = signal_group_columns(samples, year, masses, args.sig_cut_type)
            if columns:
                table = render_table(
                    caption=rf"Signal cutflow for {latex_escape_text(year)}, {mass_caption}",
                    label=f"cutflow-sig-mc-{year}-{group_name}",
                    columns=columns,
                )
                written.append(write_table(output_dir, f"cutflow_Sig_MC_{year}_{group_name}.tex", table))
                all_tables.append(table)

        columns = background_columns(samples, year, args.bkg_cut_type)
        if columns:
            table = render_table(
                caption=rf"Background cutflow for {latex_escape_text(year)}",
                label=f"cutflow-bkg-mc-{year}",
                columns=columns,
            )
            written.append(write_table(output_dir, f"cutflow_Bkg_MC_{year}.tex", table))
            all_tables.append(table)

        columns = data_columns(samples, year, args.data_cut_type)
        if columns:
            table = render_table(
                caption=rf"Data cutflow for {latex_escape_text(year)}",
                label=f"cutflow-data-{year}",
                columns=columns,
            )
            written.append(write_table(output_dir, f"cutflow_Data_{year}.tex", table))
            all_tables.append(table)

    if all_tables:
        written.append(write_table(output_dir, "cutflow_merged_all.tex", "\n".join(all_tables)))

    for path in written:
        print(f"[OK] wrote {path}")
    print(f"Done. tables={len(written) - (1 if all_tables else 0)} input_dir={input_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
