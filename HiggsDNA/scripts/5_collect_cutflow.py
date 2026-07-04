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
  - muon isolation/sip3d studies: zgammas_mu_muiso04 and zgammas_mu_munosip3d

Output:
- Writes directly usable LaTeX cutflow tables to an output text file.
- Writes the same information to a JSON sidecar file with metadata.

Usage examples:
  python collect_cutflow.py --dataset_type Sig_MC --dataset mA_M1 --year 2022postEE
  python collect_cutflow.py --path "/afs/.../eos_logs/Sig_MC/mA_M1_2022postEE"   # direct directory
  python collect_cutflow.py --out "/afs/.../job_1/*.out"                         # direct glob

Notes:
- This script only needs the .out files. No parquet is required.
"""

import os
import re
import time
import glob
import json
from dataclasses import dataclass
from collections import defaultdict, OrderedDict
from typing import Dict, List, Tuple, Optional, Iterator, Iterable, Set

start_time = time.time()

baseDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs"

dataType = ["Data", "Bkg_MC", "Sig_MC"]

outputDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list"

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
# pb
xs_sig_2022 = 0.1
xs_sig_2022 = 0.1
xs_sig_2022 = 0.1
# pb
xs_DYG_2022 = [124, 2.088] 
xs_DYG_2023 = [126.6]
xs_DYG_2024 = [126.6]
xs_DYJet_2022 = [6688.0]
xs_DYJet_2023 = [6688.0]
xs_DYJet_2024 = [2213.0,2220.66,2230.33]
pfTofb = 1000.0
# fb
lumi_2022preEE = 7.9804
lumi_2022postEE = 26.6717
lumi_2023preBPix = 17.794
lumi_2023postBPix = 9.451
lumi_2024 = 108.96

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

# --- NEW: also include weighted variants for the 15 PHID scenarios ---
PHID_CUTFLOW_TYPES_W = [f"{ct}_w" for ct in PHID_CUTFLOW_TYPES]

# sip3d scenario (your code uses ele_ip3d_cut key; cut_type is zgammas_eleip3d)
SIP3D_CUTFLOW_TYPES = [
    "zgammas_ele_elesip3d",
    "zgammas_ele_elesip3d_w",  # keep in case it exists
]

MUON_STUDY_CUTFLOW_TYPES = [
    "zgammas_mu_muiso04",
    "zgammas_mu_muiso04_w",
    "zgammas_mu_munosip3d",
    "zgammas_mu_munosip3d_w",
]

# --- CHANGED: include PHID *_w types in default list ---
DEFAULT_CUTFLOW_TYPES = (
    BASE_CUTFLOW_TYPES
    + SIP3D_CUTFLOW_TYPES
    + MUON_STUDY_CUTFLOW_TYPES
    + PHID_CUTFLOW_TYPES
)

# A minimal mapping for nicer printing
LABELS_COMMON = {
    "all": "Initial Events",
    "N_lep_sel": r"$N_{l}\geq 2$",
    "trig_cut": "Trigger",
    "lep_pt_cut": r"Lepton Trigger $p_T$ Cut",
    "has_z_cand": r"$m_{ll} > 50\,\mathrm{GeV}$",
    "g_kin_cut": "Photon Kinematic Cuts",
    "has_2g_cand": r"$N_{\gamma}\geq 2$",
    "sel_h": r"$95\,\mathrm{GeV} < m_{ll\gamma\gamma} < 180\,\mathrm{GeV}$",
    "sel_h_1": r"$m_{ll} + m_{ll\gamma\gamma} > 185\,\mathrm{GeV}$",
    "sel_h_2": r"$95\,\mathrm{GeV} < m_{ll\gamma\gamma} < 180\,\mathrm{GeV}$",
    "event": "Event Filtering",
    "all cuts": "Total Baseline Events",
}


# A minimal mapping for nicer printing
LABELS_COMMON_ele_ip3d = {
    "all": "Initial Events",
    "N_lep_sel": r"$N_{l}\geq 2$",
    "trig_cut": "Trigger",
    "lep_pt_cut": r"Lepton Trigger $p_T$ Cut",
    "ele_ip3d_cut": r"Electron SIP3D Cut",
    "has_z_cand": r"$m_{ll} > 50\,\mathrm{GeV}$",
    "g_kin_cut": "Photon Kinematic Cuts",
    "has_2g_cand": r"$N_{\gamma}\geq 2$",
    "sel_h": r"$95 < m_{ll\gamma\gamma} < 180\,\mathrm{GeV}$",
    "sel_h_1": r"$m_{ll} + m_{ll\gamma\gamma} > 185\,\mathrm{GeV}$",
    "sel_h_2": r"$95 < m_{ll\gamma\gamma} < 180\,\mathrm{GeV}$",
    "event": "Event Filtering",
    "ele_sip3d_cut": "Electron SIP3D Cut",
    "all cuts": "Total Baseline Events",
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

# --- NEW: GeneratorWeightSum parsing ---
_GENW_SUM_RE = re.compile(r"Generator_weight:\s*([+-]?\d+(?:\.\d+)?)")

def parse_generator_weight_sum(text: str) -> Optional[float]:
    """
    Parse lines like:
      INFO [AnalysisManager : GeneratorWeightSum] Sum of ...
               Generator_weight: 183430976.0
    Returns float or None if not found / not parseable.
    """
    pos = text.find("AnalysisManager : GeneratorWeightSum")
    if pos == -1:
        return None
    # search in a limited window after the tag (to avoid accidental matches)
    window = text[pos : min(len(text), pos + 2000)]
    m = _GENW_SUM_RE.search(window)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


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
    # --- NEW ---
    gen_weight_sum_by_file: Dict[str, Optional[float]]
    gen_weight_sum_total: float


# --- NEW: help function for "output-only" weighted variants ---
def _is_phid_cutflow_type(ct: str) -> bool:
    return ct.startswith("zgammas_phid_") and (not ct.endswith("_w"))

def collect_cutflows(
    out_files: Iterable[str],
    cutflow_types: List[str],
    syst: str = "nominal",
) -> CutflowResult:
    """
    Parse out files and sum 'cuts' across jobs for selected cutflow types and syst.

    NEW: For any cut_type containing '_w', divide all cut values by the per-job GeneratorWeightSum.

    CHANGED: PHID cutflows are only present without '_w' in input logs, but we additionally
    produce an output-only '<phid>_w' variant by normalizing the same cuts with GeneratorWeightSum.
    """
    wanted = set(cutflow_types)

    accum: Dict[str, OrderedDict] = defaultdict(OrderedDict)
    used_files: List[str] = []
    per_file_all: Dict[str, Dict[str, float]] = OrderedDict()
    n_jobs = 0

    gen_weight_sum_by_file: Dict[str, Optional[float]] = OrderedDict()
    gen_weight_sum_total: float = 0.0

    for fp in out_files:
        try:
            text = open(fp, "r", errors="ignore").read()
        except Exception:
            continue

        used_files.append(fp)
        n_jobs += 1
        per_file_all.setdefault(fp, OrderedDict())

        # --- NEW ---
        genw_sum = parse_generator_weight_sum(text)
        gen_weight_sum_by_file[fp] = genw_sum
        if genw_sum is not None and genw_sum != 0.0:
            gen_weight_sum_total += genw_sum

        # collect last payload per cut_type in this file to avoid double counting
        last_cuts_by_ct: Dict[str, OrderedDict] = {}

        # --- NEW: also remember last PHID cuts for output-only weighted rendering ---
        last_phid_cuts: Dict[str, OrderedDict] = {}

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

            last_cuts_by_ct[ct] = fixed_cuts

            # --- NEW: capture PHID base type for synthetic '<ct>_w' output ---
            if _is_phid_cutflow_type(ct):
                last_phid_cuts[ct] = fixed_cuts

        # apply once per (file, cut_type)
        for ct, fixed_cuts in last_cuts_by_ct.items():
            if not accum[ct]:
                for k in fixed_cuts.keys():
                    accum[ct][k] = 0.0
            else:
                for k in fixed_cuts.keys():
                    if k not in accum[ct]:
                        accum[ct][k] = 0.0

            norm = 1.0
            if "_w" in ct:
                if genw_sum is None or genw_sum == 0.0:
                    norm = 1.0
                else:
                    norm = genw_sum

            for k, v in fixed_cuts.items():
                accum[ct][k] += (v / norm)

        # --- NEW: add output-only PHID weighted variants: '<phid>_w' ---
        for ct, fixed_cuts in last_phid_cuts.items():
            ct_w = f"{ct}_w"

            if not accum[ct_w]:
                for k in fixed_cuts.keys():
                    accum[ct_w][k] = 0.0
            else:
                for k in fixed_cuts.keys():
                    if k not in accum[ct_w]:
                        accum[ct_w][k] = 0.0

            if genw_sum is None or genw_sum == 0.0:
                norm = 1.0
            else:
                norm = genw_sum

            for k, v in fixed_cuts.items():
                accum[ct_w][k] += (v / norm)

    return CutflowResult(
        cutflows=dict(accum),
        n_jobs_used=n_jobs,
        out_files_used=used_files,
        per_file_all=per_file_all,
        gen_weight_sum_by_file=gen_weight_sum_by_file,
        gen_weight_sum_total=gen_weight_sum_total,
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


# --- NEW: discover which syst values have trigger-eff payloads in given out files ---
def discover_trigeff_systs(out_files: Iterable[str]) -> Set[str]:
    """
    Scan all CutFlow JSON payloads and return the set of payload["syst"] values
    that appear together with a 'cut_type_prefix' field (trigger-eff payloads).
    """
    systs: Set[str] = set()
    for fp in out_files:
        try:
            text = open(fp, "r", errors="ignore").read()
        except Exception:
            continue
        for payload in iter_cutflow_payloads_from_text(text, syst=None):
            if payload.get("cut_type_prefix"):
                s = payload.get("syst")
                if isinstance(s, str) and s:
                    systs.add(s)
    return systs


# --- NEW: JSON rendering for trigeff ---
def render_trigeff_json_dict(result: TrigEffResult) -> Dict:
    """
    JSON-friendly structure:
    { prefix: { "overall": {...}, "bins": {ptbin: {...}} } }
    """
    out: Dict[str, Dict] = OrderedDict()
    for prefix, d in result.by_prefix.items():
        out[prefix] = OrderedDict(
            overall=OrderedDict(
                in_bin=float(d["overall"]["in_bin"]),
                pass_trigger=float(d["overall"]["pass_trigger"]),
                eff=float(d["overall"]["eff"]),
            ),
            bins=OrderedDict(
                (ptbin, OrderedDict(
                    in_bin=float(bb["in_bin"]),
                    pass_trigger=float(bb["pass_trigger"]),
                    eff=float(bb["eff"]),
                ))
                for ptbin, bb in d["bins"].items()
            ),
        )
    return out


def render_trigeff(result: TrigEffResult, sample_label: Optional[str] = None) -> str:
    """
    Render trigger-efficiency payloads as LaTeX tables.
    """
    lines: List[str] = []
    if not result.by_prefix:
        return ""

    for prefix, d in result.by_prefix.items():
        o = d["overall"]
        caption = (
            rf"Trigger efficiency for {_latex_texttt(sample_label)}: {_latex_texttt(prefix)}"
            if sample_label
            else rf"Trigger efficiency: {_latex_texttt(prefix)}"
        )
        label_parts = [p for p in (sample_label, prefix) if p]
        lines.extend([
            r"\begin{table}[htbp]",
            r"\centering",
            rf"\caption{{{caption}}}",
            rf"\label{{tab:trigeff-{_latex_label_id('-'.join(label_parts))}}}",
            r"\begin{tabular}{lrrr}",
            r"\hline",
            r"$p_T$ bin & In bin & Pass trigger & Efficiency \\",
            r"\hline",
            f"Overall & {o['in_bin']:.0f} & {o['pass_trigger']:.0f} & {o['eff']:.6f} \\\\",
        ])
        for ptbin, bb in d["bins"].items():
            lines.append(
                f"{_latex_escape_text(ptbin)} & {bb['in_bin']:.0f} & {bb['pass_trigger']:.0f} & {bb['eff']:.6f} \\\\"
            )
        lines.extend([
            r"\hline",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ])
    return "\n".join(lines)


# -----------------------------
# Pretty printing
# -----------------------------
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


def _latex_escape_text(s: str) -> str:
    return "".join(_LATEX_ESCAPE_MAP.get(ch, ch) for ch in str(s))


def _latex_texttt(s: str) -> str:
    return rf"\texttt{{{_latex_escape_text(s)}}}"


def _latex_label_id(s: str) -> str:
    label = re.sub(r"[^A-Za-z0-9]+", "-", str(s)).strip("-").lower()
    return label or "table"


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
        return LABELS_COMMON_ele_ip3d.get(cut_key, _latex_escape_text(cut_key))

    return LABELS_COMMON.get(cut_key, _latex_escape_text(cut_key))


def format_value(ct: str, v: float) -> str:
    # weighted cutflows typically are non-integers; keep higher precision
    if "_w" in ct:
        return f"{v:.6f}"
    # otherwise print as integer-like
    return f"{v:.0f}"


def render_cutflows(
    result: CutflowResult,
    xs_pb: Optional[float] = None,
    lumi_fb: Optional[float] = None,
    sample_label: Optional[str] = None,
) -> str:
    """
    Render cutflows as directly usable LaTeX tables.

    The xs/lumi arguments are kept for call-site compatibility; text output
    intentionally mirrors the same unscaled values that were previously written.
    """
    lines: List[str] = []

    for ct, cuts in result.cutflows.items():
        caption = (
            rf"Cutflow for {_latex_texttt(sample_label)}: {_latex_texttt(ct)}"
            if sample_label
            else rf"Cutflow: {_latex_texttt(ct)}"
        )
        label_parts = [p for p in (sample_label, ct) if p]
        lines.extend([
            r"\begin{table}[htbp]",
            r"\centering",
            rf"\caption{{{caption}}}",
            rf"\label{{tab:cutflow-{_latex_label_id('-'.join(label_parts))}}}",
            r"\begin{tabular}{lr}",
            r"\hline",
            r"Selection & Events \\",
            r"\hline",
        ])
        for cut_key, v in cuts.items():
            vv = v
            lines.append(f"{label_for(ct, cut_key)} & {format_value(ct, vv)} \\\\")
        lines.extend([
            r"\hline",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ])

    return "\n".join(lines)


def render_cutflows_json_dict(
    result: CutflowResult,
    syst: str = "nominal",
    xs_pb: Optional[float] = None,
    lumi_fb: Optional[float] = None,
) -> Dict:
    """
    JSON-friendly structure with raw cut keys.

    NEW: for any ct containing '_w', additionally scale by xs * pfTofb * lumi.
    """
    w_scale = 1.0
    if xs_pb is not None and lumi_fb is not None:
        w_scale = float(xs_pb) * float(pfTofb) * float(lumi_fb)

    cutflows_out: Dict[str, Dict[str, float]] = OrderedDict()
    for ct, cuts in result.cutflows.items():
        if "_w" in ct:
            cutflows_out[ct] = OrderedDict((k, float(v) * w_scale) for k, v in cuts.items())
        else:
            cutflows_out[ct] = OrderedDict((k, float(v)) for k, v in cuts.items())

    return OrderedDict(
        meta=OrderedDict(
            syst=syst,
            n_jobs_used=result.n_jobs_used,
            out_files_used=list(result.out_files_used),
            gen_weight_sum_total=float(result.gen_weight_sum_total),
            gen_weight_sum_by_file=OrderedDict(
                (fp, (None if v is None else float(v))) for fp, v in result.gen_weight_sum_by_file.items()
            ),
            # --- NEW ---
            xs_pb=(None if xs_pb is None else float(xs_pb)),
            lumi_fb=(None if lumi_fb is None else float(lumi_fb)),
            pfTofb=float(pfTofb),
            w_scale=float(w_scale),
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

        # --- NEW: determine xs and lumi for scaling *_w outputs ---
        xs_pb: Optional[float] = None
        lumi_fb: Optional[float] = None

        if dataset_type in ("Sig_MC", "Bkg_MC"):
            # xs (pb)
            if dataset_type == "Sig_MC":
                xs_pb = float(xs_sig_2022)  # NOTE: your code currently uses same value for all years
            else:
                if dataset in name_DYG_2022:
                    xs_pb = float(xs_DYG_2022[name_DYG_2022.index(dataset)])
                elif dataset in name_DYG_2023:
                    xs_pb = float(xs_DYG_2023[name_DYG_2023.index(dataset)])
                elif dataset in name_DYG_2024:
                    xs_pb = float(xs_DYG_2024[name_DYG_2024.index(dataset)])
                elif dataset in name_DYJet_2022:
                    xs_pb = float(xs_DYJet_2022[name_DYJet_2022.index(dataset)])
                elif dataset in name_DYJet_2023:
                    xs_pb = float(xs_DYJet_2023[name_DYJet_2023.index(dataset)])
                elif dataset in name_DYJet_2024:
                    xs_pb = float(xs_DYJet_2024[name_DYJet_2024.index(dataset)])

            # lumi (fb^-1)
            if year == "2022preEE":
                lumi_fb = float(lumi_2022preEE)
            elif year == "2022postEE":
                lumi_fb = float(lumi_2022postEE)
            elif year == "2023preBPix":
                lumi_fb = float(lumi_2023preBPix)
            elif year == "2023postBPix":
                lumi_fb = float(lumi_2023postBPix)
            elif year == "2024":
                lumi_fb = float(lumi_2024)

        result = collect_cutflows(out_files, cutflow_types=cutflow_types, syst=syst)
        sample_label = f"{dataset_type}_{dataset}_{year}"
        text = render_cutflows(result, xs_pb=xs_pb, lumi_fb=lumi_fb, sample_label=sample_label)

        trigeff_systs = discover_trigeff_systs(out_files)
        trigeff_json = None
        if "nominal" in trigeff_systs:
            trigeff = collect_trigeff(out_files, syst="nominal")
            trigeff_text = render_trigeff(trigeff, sample_label=sample_label)
            if trigeff_text:
                text = text + "\n" + trigeff_text
            trigeff_json = render_trigeff_json_dict(trigeff)

        out_path = os.path.join(out_dir, f"cutflow_{dataset_type}_{dataset}_{year}.txt")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(text + "\n")

        out_json_path = os.path.join(out_dir, f"cutflow_{dataset_type}_{dataset}_{year}.json")
        payload = render_cutflows_json_dict(result, syst=syst, xs_pb=xs_pb, lumi_fb=lumi_fb)
        if trigeff_json is not None:
            payload["trigeff_nominal"] = trigeff_json
        payload["meta"]["trigeff_systs_present"] = sorted(trigeff_systs)

        with open(out_json_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)

        n_ok += 1
        print(f"[OK] {dataset_type} {dataset} {year}  jobs={result.n_jobs_used}  -> {out_path}")

    elapsed = time.time() - start_time
    print(f"\nDone. OK={n_ok}, SKIP={n_skip}, elapsed={time.strftime('%H:%M:%S', time.gmtime(elapsed))}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
