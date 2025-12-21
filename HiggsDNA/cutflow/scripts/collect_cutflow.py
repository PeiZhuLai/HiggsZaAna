#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collect_cutflow.py

Collects cutflow counters from HiggsDNA tagger .out files under:
  <base_dir>/<dataset_type>/<dataset>_<year>/job_*/<N>.out

Key features (requested):
- If a job directory contains >=2 .out files, pick the one with the largest numeric name (e.g. 10633893.0.out).
- Parse "CutFlow JSON {...}" payloads even when the log wraps the JSON across lines (including inside strings).
- Support the default 6 CUTFLOW_TYPES plus:
  - 15 photonID scenarios: 5 kinds × 3 WPs (tight/medium/loose)
  - 1 sip3d scenario: zgammas_eleip3d (and optionally its _w variant if present)

Output:
- Prints a LaTeX-friendly cutflow table rows to stdout.
- Writes the same content to an output text file.

Usage examples:
  python collect_cutflow.py --dataset_type Sig_MC --dataset mA_M1 --year 2022postEE
  python collect_cutflow.py --path "/afs/.../eos_logs/Sig_MC/mA_M1_2022postEE"   # direct directory
  python collect_cutflow.py --out "/afs/.../job_1/*.out"                         # direct glob

Notes:
- This script only needs the .out files. No parquet is required.
"""

import os
import numpy as np
import pandas as pd
import re
from pdb import set_trace
import time
import glob
import json
from dataclasses import dataclass
from collections import defaultdict, OrderedDict
from typing import Dict, List, Tuple, Optional, Iterator, Iterable

start_time = time.time()

baseDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs"

dataType = ["Data", "Bkg_MC", "Sig_MC"]

outputDir = "../cutflow_list"

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
name_Data_2022 = ["DYJetsToLL"]
name_Data_2023 = ["DYJetsToLL"]
name_Data_2024 = ["DYJetsToLL"]

# -----------------------------
# Defaults
# -----------------------------
DEFAULT_BASE_DIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs"

# The 6 "normal" types you listed
BASE_CUTFLOW_TYPES = [
    "zgammas",
    "zgammas_ele",
    "zgammas_mu",
    "zgammas_w",
    "zgammas_ele_w",
    "zgammas_mu_w",
]

# photon ID scenarios (5 kinds × 3 wps = 15)
_PHID_WPS = ["tight", "medium", "loose"]
_PHID_KINDS = ["custom", "custom_extend", "sieie", "PFECalIso", "official"]
PHID_CUTFLOW_TYPES = [
    f"zgammas_phid_{kind}_{wp}"
    for wp in _PHID_WPS
    for kind in _PHID_KINDS
]

# sip3d scenario (your code uses ele_ip3d_cut key; cut_type is zgammas_eleip3d)
SIP3D_CUTFLOW_TYPES = [
    "zgammas_eleip3d",
    "zgammas_eleip3d_w",  # keep in case it exists
]

DEFAULT_CUTFLOW_TYPES = BASE_CUTFLOW_TYPES + PHID_CUTFLOW_TYPES + SIP3D_CUTFLOW_TYPES

# A minimal mapping for nicer printing
LABELS_COMMON = {
    "all": "Initial Events",
    "N_lep_sel": r"$N_{l}\geq 2$",
    "trig_cut": "Trigger",
    "lep_pt_cut": r"Lepton trigger $p_T$ cut",
    "has_z_cand": r"$m_{ll} > 50\,\mathrm{GeV}$",
    "g_kin_cut": "Photon kinematic cuts",
    "has_2g_cand": r"$N_{\gamma}\geq 2$",
    "sel_h_1": r"$m_{ll} + m_{ll\gamma\gamma} > 185\,\mathrm{GeV}$",
    "sel_h_2": r"$95\,\mathrm{GeV} < m_{ll\gamma\gamma} < 180\,\mathrm{GeV}$",
    "event": "Event filtering",
    "all cuts": "All cuts",
}


# A minimal mapping for nicer printing
LABELS_COMMON_ele_ip3d = {
    "all": "Initial Events",
    "N_lep_sel": r"$N_{l}\geq 2$",
    "trig_cut": "Trigger",
    "lep_pt_cut": r"Lepton trigger $p_T$ cut",
    "ele_ip3d_cut": r"Electron sip3d cut",
    "has_z_cand": r"$m_{ll} > 50\,\mathrm{GeV}$",
    "g_kin_cut": "Photon kinematic cuts",
    "has_2g_cand": r"$N_{\gamma}\geq 2$",
    "sel_h_1": r"$m_{ll} + m_{ll\gamma\gamma} > 185\,\mathrm{GeV}$",
    "sel_h_2": r"$95\,\mathrm{GeV} < m_{ll\gamma\gamma} < 180\,\mathrm{GeV}$",
    "event": "Event filtering",
    "all cuts": "All cuts",
}
# -----------------------------
# JSON de-wrapping helpers
# -----------------------------
def _dewrap_json(s: str) -> str:
    """
    HiggsDNA logger prints:  logger.debug("CutFlow JSON %s", json.dumps(payload, separators=(",", ":")))

    In some .out files the JSON is line-wrapped by the logging/terminal layer, even *inside* JSON strings.
    That introduces literal newlines and indentation inside strings => invalid JSON.

    Strategy:
    - Remove all whitespace outside strings.
    - Inside strings: remove line breaks and the indentation after them, and also strip trailing spaces before
      the line break (these are wrap artifacts).
    This reconstructs a valid JSON string in practice for these logs.
    """
    out: List[str] = []
    in_str = False
    esc = False
    i = 0
    n = len(s)

    while i < n:
        ch = s[i]
        if in_str:
            if ch in "\r\n":
                # remove trailing wrap-padding spaces just before line break
                while out and out[-1] == " ":
                    out.pop()
                i += 1
                # skip indentation after the line break
                while i < n and s[i] in " \t":
                    i += 1
                continue

            if esc:
                out.append(ch)
                esc = False
                i += 1
                continue

            if ch == "\\":
                out.append(ch)
                esc = True
                i += 1
                continue

            if ch == '"':
                in_str = False
                out.append(ch)
                i += 1
                continue

            out.append(ch)
            i += 1
        else:
            if ch in " \t\r\n":
                i += 1
                continue
            if ch == '"':
                in_str = True
                out.append(ch)
                i += 1
                continue
            out.append(ch)
            i += 1

    return "".join(out)


def _extract_braced(text: str, start: int) -> Tuple[Optional[str], int]:
    """Extract a {...} JSON object starting at index start (which must be '{')."""
    depth = 0
    in_str = False
    esc = False
    for i in range(start, len(text)):
        ch = text[i]
        if in_str:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
        else:
            if ch == '"':
                in_str = True
            elif ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    return text[start : i + 1], i + 1
    return None, start + 1


def iter_cutflow_payloads_from_text(text: str, syst: Optional[str] = None) -> Iterator[Dict]:
    """
    Yield dict payloads for each 'CutFlow JSON {..}' (or raw '{..}' with 'tagger' key).
    Handles wrapped JSON via _dewrap_json.

    If syst is not None, only yield payloads with payload["syst"] == syst.
    """
    idx = 0
    # Prefer the explicit prefix
    while True:
        pos = text.find("CutFlow JSON", idx)
        if pos == -1:
            break
        brace = text.find("{", pos)
        if brace == -1:
            idx = pos + 10
            continue
        obj, end = _extract_braced(text, brace)
        if not obj:
            idx = brace + 1
            continue

        cleaned = _dewrap_json(obj)
        try:
            payload = json.loads(cleaned)
            if (syst is None) or (payload.get("syst") == syst):
                yield payload
        except Exception:
            pass
        idx = end

    # Fallback: sometimes the prefix may be missing, but payload starts with {"tagger":...}
    idx = 0
    while True:
        pos = text.find('{"tagger"', idx)
        if pos == -1:
            break
        obj, end = _extract_braced(text, pos)
        if not obj:
            idx = pos + 1
            continue
        cleaned = _dewrap_json(obj)
        try:
            payload = json.loads(cleaned)
            if (syst is None) or (payload.get("syst") == syst):
                yield payload
        except Exception:
            pass
        idx = end


# -----------------------------
# File selection
# -----------------------------
_NUM_OUT_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)\.out$")


def pick_best_out_file(out_files: List[str]) -> Optional[str]:
    """
    If multiple *.out exist:
    - Prefer those whose filename ends with a number like 10633893.0.out; pick the largest numeric one.
    - Else fall back to newest mtime.
    """
    if not out_files:
        return None

    scored: List[Tuple[float, str]] = []
    for fp in out_files:
        m = _NUM_OUT_RE.search(os.path.basename(fp))
        if m:
            try:
                scored.append((float(m.group(1)), fp))
            except Exception:
                pass
    if scored:
        scored.sort(key=lambda x: x[0])
        return scored[-1][1]

    # Fallback: newest modified
    out_files.sort(key=lambda x: os.path.getmtime(x))
    return out_files[-1]


def iter_job_out_files(sample_dir: str) -> Iterator[str]:
    """Yield the chosen .out file for each job_* directory under sample_dir."""
    job_dirs = sorted(
        d for d in glob.glob(os.path.join(sample_dir, "job_*"))
        if os.path.isdir(d)
    )
    for jd in job_dirs:
        outs = glob.glob(os.path.join(jd, "*.out"))
        best = pick_best_out_file(outs)
        if best:
            yield best


# -----------------------------
# Cutflow collection
# -----------------------------
@dataclass
class CutflowResult:
    cutflows: Dict[str, "OrderedDict[str, float]"]
    n_jobs_used: int
    out_files_used: List[str]
    per_file_all: Dict[str, Dict[str, float]]  # file -> {cut_type -> all}


def collect_cutflows(
    out_files: Iterable[str],
    cutflow_types: List[str],
    syst: str = "nominal",
) -> CutflowResult:
    """
    Parse out files and sum 'cuts' across jobs for selected cutflow types and syst.
    """
    wanted = set(cutflow_types)

    accum: Dict[str, OrderedDict] = defaultdict(OrderedDict)
    used_files: List[str] = []
    per_file_all: Dict[str, Dict[str, float]] = OrderedDict()
    n_jobs = 0

    for fp in out_files:
        try:
            text = open(fp, "r", errors="ignore").read()
        except Exception:
            continue

        used_files.append(fp)
        n_jobs += 1
        per_file_all.setdefault(fp, OrderedDict())

        # NEW: collect last payload per cut_type in this file to avoid double counting
        last_cuts_by_ct: Dict[str, OrderedDict] = {}

        for payload in iter_cutflow_payloads_from_text(text, syst=syst):
            ct = payload.get("cut_type")
            if ct not in wanted:
                continue
            cuts = payload.get("cuts")
            if not isinstance(cuts, dict):
                continue

            fixed_cuts = OrderedDict()
            for k, v in cuts.items():
                kk = "all cuts" if k == "allcuts" else k
                fixed_cuts[kk] = float(v)

            # record per-file "all"
            if "all" in fixed_cuts:
                per_file_all[fp][ct] = float(fixed_cuts["all"])
            elif "all cuts" in fixed_cuts:
                per_file_all[fp][ct] = float(fixed_cuts["all cuts"])

            # keep only the last one per cut_type for this file
            last_cuts_by_ct[ct] = fixed_cuts

        # NEW: apply once per (file, cut_type)
        for ct, fixed_cuts in last_cuts_by_ct.items():
            if not accum[ct]:
                for k in fixed_cuts.keys():
                    accum[ct][k] = 0.0
            else:
                for k in fixed_cuts.keys():
                    if k not in accum[ct]:
                        accum[ct][k] = 0.0
            for k, v in fixed_cuts.items():
                accum[ct][k] += v

    return CutflowResult(
        cutflows=dict(accum),
        n_jobs_used=n_jobs,
        out_files_used=used_files,
        per_file_all=per_file_all,
    )


# -----------------------------
# Trigger-eff (cut_type_prefix) collection
# -----------------------------
@dataclass
class TrigEffResult:
    by_prefix: Dict[str, Dict]   # prefix -> {"bins": OrderedDict[ptbin]->dict, "overall": dict}
    n_jobs_used: int
    out_files_used: List[str]


def _safe_float(x, default: float = 0.0) -> float:
    try:
        return float(x)
    except Exception:
        return default


def collect_trigeff(
    out_files: Iterable[str],
    syst: str = "nominal",
) -> TrigEffResult:
    """
    Collect payloads like:
      {"syst":"nominal","cut_type_prefix":"...","bins":[{ptbin,in_bin,pass_trigger,eff},...],
       "overall":{"in_bin":...,"pass_trigger":...,"eff":...}}

    We sum in_bin/pass_trigger across jobs and recompute eff = pass/in.
    """
    acc: Dict[str, Dict] = {}
    used_files: List[str] = []
    n_jobs = 0

    for fp in out_files:
        try:
            text = open(fp, "r", errors="ignore").read()
        except Exception:
            continue

        used_files.append(fp)
        n_jobs += 1

        for payload in iter_cutflow_payloads_from_text(text, syst=syst):
            prefix = payload.get("cut_type_prefix")
            if not prefix:
                continue

            bins = payload.get("bins", [])
            overall = payload.get("overall", {})

            if prefix not in acc:
                acc[prefix] = {
                    "bins": OrderedDict(),   # ptbin -> {"in_bin": float, "pass_trigger": float}
                    "overall": {"in_bin": 0.0, "pass_trigger": 0.0},
                }

            # overall
            acc[prefix]["overall"]["in_bin"] += _safe_float(overall.get("in_bin"))
            acc[prefix]["overall"]["pass_trigger"] += _safe_float(overall.get("pass_trigger"))

            # bins
            if isinstance(bins, list):
                for b in bins:
                    if not isinstance(b, dict):
                        continue
                    ptbin = b.get("ptbin")
                    if not ptbin:
                        continue
                    if ptbin not in acc[prefix]["bins"]:
                        acc[prefix]["bins"][ptbin] = {"in_bin": 0.0, "pass_trigger": 0.0}
                    acc[prefix]["bins"][ptbin]["in_bin"] += _safe_float(b.get("in_bin"))
                    acc[prefix]["bins"][ptbin]["pass_trigger"] += _safe_float(b.get("pass_trigger"))

    # recompute eff fields for rendering convenience
    out: Dict[str, Dict] = OrderedDict()
    for prefix, d in sorted(acc.items(), key=lambda kv: kv[0]):
        o_in = d["overall"]["in_bin"]
        o_pass = d["overall"]["pass_trigger"]
        o_eff = (o_pass / o_in) if o_in > 0 else 0.0

        bins_od = OrderedDict()
        for ptbin, bb in d["bins"].items():
            b_in = bb["in_bin"]
            b_pass = bb["pass_trigger"]
            b_eff = (b_pass / b_in) if b_in > 0 else 0.0
            bins_od[ptbin] = {"in_bin": b_in, "pass_trigger": b_pass, "eff": b_eff}

        out[prefix] = {
            "overall": {"in_bin": o_in, "pass_trigger": o_pass, "eff": o_eff},
            "bins": bins_od,
        }

    return TrigEffResult(by_prefix=out, n_jobs_used=n_jobs, out_files_used=used_files)


def render_trigeff(result: TrigEffResult) -> str:
    """
    LaTeX-friendly:
      TrigEff: <cut_type_prefix>
        overall & in & pass & eff \\
        <ptbin> & in & pass & eff \\
    """
    lines: List[str] = []
    if not result.by_prefix:
        return ""

    lines.append("Trigger efficiency payloads:")
    for prefix, d in result.by_prefix.items():
        o = d["overall"]
        lines.append(f"TrigEff: {prefix}")
        lines.append(f"  overall{'':30} & {o['in_bin']:.0f} & {o['pass_trigger']:.0f} & {o['eff']:.6f} \\\\")
        for ptbin, bb in d["bins"].items():
            lines.append(f"  {ptbin:35} & {bb['in_bin']:.0f} & {bb['pass_trigger']:.0f} & {bb['eff']:.6f} \\\\")
        lines.append("")
    return "\n".join(lines)


# -----------------------------
# Pretty printing
# -----------------------------
def label_for(ct: str, cut_key: str) -> str:
    # Provide a bit more specific labels for ele/mu channels if needed
    if cut_key == "N_lep_sel":
        if "_ele" in ct:
            return r"$N_{e}\geq 2$"
        if "_mu" in ct:
            return r"$N_{\mu}\geq 2$"
    if cut_key == "has_z_cand":
        if "_ele" in ct:
            return r"$m_{ee} > 50\,\mathrm{GeV}$"
        if "_mu" in ct:
            return r"$m_{\mu\mu} > 50\,\mathrm{GeV}$"

    # sip3d variants should use the dedicated label mapping (includes ele_ip3d_cut)
    if ct in set(SIP3D_CUTFLOW_TYPES):
        return LABELS_COMMON_ele_ip3d.get(cut_key, cut_key)

    return LABELS_COMMON.get(cut_key, cut_key)


def format_value(ct: str, v: float) -> str:
    # weighted cutflows typically are non-integers; keep higher precision
    if "_w" in ct:
        return f"{v:.6f}"
    # otherwise print as integer-like
    return f"{v:.0f}"


def render_cutflows(result: CutflowResult) -> str:
    """
    Render in a LaTeX-friendly style:
      Cut Type: <ct>
          <label> & <value> \\
    """
    lines: List[str] = []
    lines.append(f"# jobs used: {result.n_jobs_used}")
    for fp in result.out_files_used:
        lines.append(f"#   {fp}")
    lines.append("")

    for ct, cuts in result.cutflows.items():
        lines.append(f"Cut Type: {ct}")
        for cut_key, v in cuts.items():
            lines.append(f"  {label_for(ct, cut_key):35} & {format_value(ct, v)} \\\\")
        lines.append("")

    return "\n".join(lines)


# --- NEW: JSON rendering (raw cut keys, no label replacement) ---
def render_cutflows_json_dict(
    result: CutflowResult,
    syst: str = "nominal",
) -> Dict:
    """
    JSON-friendly structure with raw cut keys.
    {
      "meta": {...},
      "cutflows": {cut_type: {cut_key: value, ...}, ...}
    }
    """
    cutflows_out: Dict[str, Dict[str, float]] = OrderedDict()
    for ct, cuts in result.cutflows.items():
        cutflows_out[ct] = OrderedDict((k, float(v)) for k, v in cuts.items())

    return OrderedDict(
        meta=OrderedDict(
            syst=syst,
            n_jobs_used=result.n_jobs_used,
            out_files_used=list(result.out_files_used),
        ),
        cutflows=cutflows_out,
    )


# -----------------------------
# CLI
# -----------------------------
def build_sample_dir(base_dir: str, dataset_type: str, dataset: str, year: str) -> str:
    return os.path.join(base_dir, dataset_type, f"{dataset}_{year}")


def iter_all_samples() -> Iterator[Tuple[str, str, str]]:
    """
    Auto-loop all samples based on the lists defined at top of file.
    Yields: (dataset_type, dataset, year)
    """
    # Sig
    for y in year_sig_2022:
        for ds in name_sig_2022:
            yield ("Sig_MC", ds, y)
    for y in year_sig_2023:
        for ds in name_sig_2023:
            yield ("Sig_MC", ds, y)
    for y in year_sig_2024:
        for ds in name_sig_2024:
            yield ("Sig_MC", ds, y)

    # Bkg (DYG)
    for y in year_DYG_2022:
        for ds in name_DYG_2022:
            yield ("Bkg_MC", ds, y)
    for y in year_DYG_2023:
        for ds in name_DYG_2023:
            yield ("Bkg_MC", ds, y)
    for y in year_DYG_2024:
        for ds in name_DYG_2024:
            yield ("Bkg_MC", ds, y)

    # Bkg (DYJets)
    for y in years_DYJet_2022:
        for ds in name_DYJet_2022:
            yield ("Bkg_MC", ds, y)
    for y in year_DYJet_2023:
        for ds in name_DYJet_2023:
            yield ("Bkg_MC", ds, y)
    for y in years_DYJet_2024:
        for ds in name_DYJet_2024:
            yield ("Bkg_MC", ds, y)

    # Data
    for y in years_Data_2022:
        for ds in name_Data_2022:
            yield ("Data", ds, y)
    for y in year_Data_2023:
        for ds in name_Data_2023:
            yield ("Data", ds, y)
    for y in years_Data_2024:
        for ds in name_Data_2024:
            yield ("Data", ds, y)


def main() -> int:
    base_dir = DEFAULT_BASE_DIR  # use your fixed eos_logs base
    cutflow_types = DEFAULT_CUTFLOW_TYPES
    syst = "nominal"

    out_dir = outputDir
    os.makedirs(out_dir, exist_ok=True)

    n_ok = 0
    n_skip = 0

    for dataset_type, dataset, year in iter_all_samples():
        sample_dir = build_sample_dir(base_dir, dataset_type, dataset, year)
        if not os.path.isdir(sample_dir):
            n_skip += 1
            print(f"[SKIP] missing dir: {sample_dir}")
            continue

        out_files = list(iter_job_out_files(sample_dir))
        if not out_files:
            n_skip += 1
            print(f"[SKIP] no job_*/.out under: {sample_dir}")
            continue

        result = collect_cutflows(out_files, cutflow_types=cutflow_types, syst=syst)
        text = render_cutflows(result)

        # NEW: trigger-eff payloads
        trigeff = collect_trigeff(out_files, syst=syst)
        trigeff_text = render_trigeff(trigeff)
        if (trigeff_text):
            text = text + "\n" + trigeff_text

        out_path = os.path.join(out_dir, f"cutflow_{dataset_type}_{dataset}_{year}.txt")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(text + "\n")

        # --- NEW: write JSON output (raw cut keys like N_lep_sel) ---
        out_json_path = os.path.join(out_dir, f"cutflow_{dataset_type}_{dataset}_{year}.json")
        payload = render_cutflows_json_dict(result, syst=syst)
        with open(out_json_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)

        n_ok += 1
        print(f"[OK] {dataset_type} {dataset} {year}  jobs={result.n_jobs_used}  -> {out_path}")

    print(f"\nDone. OK={n_ok}, SKIP={n_skip}, elapsed={time.time()-start_time:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())