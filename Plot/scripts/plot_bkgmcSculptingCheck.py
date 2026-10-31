import argparse
import json
import re
import sys
import time
from array import array
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import ROOT
import uproot

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_DIR = SCRIPT_DIR.parent
LIB_DIR = PLOT_DIR / "lib"
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

# Allow importing the sideband reweighter helper that lives next to the BDT.
_HZAMVA_SCRIPTS_DIR = SCRIPT_DIR.parent.parent / "HZaMVA" / "scripts"
if str(_HZAMVA_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_HZAMVA_SCRIPTS_DIR))

from ROOT import TH1F, THStack, TCanvas, TLegend, TPad, TColor, gROOT

from Analyzer_Configs import Analyzer_Config
from Plot_Configs import Plot_Config
from Plot_Helper import LoadNtuples, MakeRatioPlot, Get_StatUnc, Draw_unc, SaveCanvPic
import CMS_lumi
import tdrstyle

try:
    from sideband_reweight import (
        DERIVED_VARS,
        VAR_ALIASES,
        _mass_from_event,
        _mass_from_row_order,
        load_sideband_reweighter,
    )
except Exception as _exc:  # pragma: no cover - best-effort import
    load_sideband_reweighter = None
    print(f"[Warning] Cannot import sideband_reweight helper: {_exc}")


DEFAULT_MVA_CUT_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
DEFAULT_OUTPUT_DIR = "/eos/home-p/pelai/HZa/output_plots/bkgmcScupltingCheck"
DEFAULT_CONTROL_SAMPLE_DIR = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_control"
DEFAULT_SIGEFF_ELE_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json"
DEFAULT_SIGEFF_MU_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json"
LOCAL_MVA_CUT_CANDIDATES = [
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/HiggsZaAna/Plot/output/MVAcut_points_run3.json",
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/optimize_alias/MVAcut_points_run3.json",
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/selection_alias/MVAcut_points_run3.json",
]
LOCAL_SIGEFF_ELE_CANDIDATES = [
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/HiggsZaAna/Plot/output/sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json",
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/optimize_alias/sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json",
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/selection_alias/sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json",
]
LOCAL_SIGEFF_MU_CANDIDATES = [
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/HiggsZaAna/Plot/output/sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json",
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/optimize_alias/sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json",
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/selection_alias/sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json",
]
LOCAL_OUTPUT_DIR = PLOT_DIR / "plots" / "bkgmcScupltingCheck"

TARGET_MASSES = list(range(1, 31))
SIGNAL_MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
BLIND_LOW = 115.0
BLIND_HIGH = 135.0
H_M_NBINS = 85
H_M_XMIN = 95.0
H_M_XMAX = 180.0
H_M_BIN_COARSE_WIDTH = 5.0
H_M_BIN_COARSE_XMAX = 180.0
H_M_BIN_COARSE10_WIDTH = 10.0
H_M_BIN_COARSE10_XMAX = 185.0   # 95..185 -> 9 bins of 10 GeV; edges align with SR (115/135) and peak (125)
H_M_BIN_SIGNAL_OVERLAY_WIDTH = 2.0
H_M_BIN_SIGNAL_OVERLAY_XMAX = 181.0
BDT_SHAPE_NBINS = 10
BDT_SHAPE_SCORE_MIN = 0.0
BDT_SHAPE_SCORE_MAX = 1.0
BDT_SHAPE_MASS_BIN_WIDTH = 5.0
BDT_SHAPE_LUMI_TEXT = "172.13 fb^{-1} (13.6 TeV)"
BDT_SHAPE_CMS_TEXT_SIZE = 0.050
BDT_SHAPE_PRELIM_TEXT_SIZE = 0.042
BDT_SHAPE_LUMI_TEXT_SIZE = 0.044
SIGNAL_DRAW_SCALE = 0.5
SIGNAL_DRAW_SCALE_NEAREST_BIN5 = 0.2
LUMI_MAP = {
    "2022preEE": 7.99,
    "2022postEE": 26.68,
    "2023preBPix": 17.96,
    "2023postBPix": 9.68,
    "2024": 109.82,
}
BKG_LABELS = {"DYJetsToLL": "Z + jets", "DYGto2LG": "Z + #gamma"}
_EFF_CACHE: Dict[str, dict] = {}
DATA_MARKER_SIZE = 1.8
BDT_SHAPE_COLORS = [
    ROOT.kBlack,
    ROOT.kBlue + 1,
    ROOT.kAzure + 7,
    ROOT.kCyan + 2,
    ROOT.kTeal + 3,
    ROOT.kGreen + 2,
    ROOT.kSpring + 5,
    ROOT.kOrange + 1,
    ROOT.kRed + 1,
    ROOT.kMagenta,
]
BDT_SHAPE_LINE_STYLES = [1, 2, 7, 9]
# Per-bin overrides (BDT bin index -> line style). Force the highest BDT bin
# (0.9 < BDT <= 1.0, idx=9) to a solid line so it does not get confused with idx=8.
BDT_SHAPE_LINE_STYLE_OVERRIDES = {
    9: 1,
}
SIGNAL_COLOR_HEX = {
    "M1":  "#2563EB",
    "M2":  "#14B8A6",
    "M3":  "#10B981",
    "M4":  "#22C55E",
    "M5":  "#65A30D",
    "M6":  "#84CC16",
    "M7":  "#A3E635",
    "M8":  "#D9F99D",
    "M9":  "#FDE047",
    "M10": "#CA8A04",
    "M15": "#DC2626",
    "M20": "#DB2777",
    "M25": "#6D28D9",
    "M30": "#111827",
}


def _to_int(value) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return int(round(value))
    if isinstance(value, str):
        match = re.search(r"-?\d+", value)
        return int(match.group(0)) if match else None
    return None


def _to_float(value) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        match = re.search(r"-?\d+(?:\.\d+)?", value)
        return float(match.group(0)) if match else None
    return None


def _parse_ma(name: str) -> Optional[int]:
    match = re.search(r"[Mm](\d+)", name)
    return int(match.group(1)) if match else None


def _parse_mva_cuts(path: str) -> Dict[int, float]:
    def _entry_to_pair(entry: dict):
        mass_keys = ("mA", "ma", "mass", "massA")
        cut_keys = ("MVAcut", "mvaCut", "mva_cut", "cut", "bdtCut", "bdt_cut", "best_MVAcut")

        ma = None
        for key in mass_keys:
            if key in entry:
                ma = _to_int(entry[key])
                break

        if ma is None:
            for key in ("sample", "name", "label", "title"):
                if key in entry and isinstance(entry[key], str):
                    ma = _parse_ma(entry[key])
                    if ma is not None:
                        break

        cut = None
        for key in cut_keys:
            if key in entry:
                cut = _to_float(entry[key])
                break

        if cut is None:
            for key in ("best", "opt", "result", "payload"):
                if key in entry and isinstance(entry[key], dict):
                    for cut_key in cut_keys:
                        if cut_key in entry[key]:
                            cut = _to_float(entry[key][cut_key])
                            break
                if cut is not None:
                    break

        return ma, cut

    with open(path, "r") as handle:
        data = json.load(handle)

    out: Dict[int, float] = {}

    if isinstance(data, list):
        for entry in data:
            if not isinstance(entry, dict):
                continue
            ma, cut = _entry_to_pair(entry)
            if ma is not None and cut is not None:
                out[ma] = cut

    if not out and isinstance(data, dict):
        for key in ("results", "points", "entries"):
            entries = data.get(key)
            if not isinstance(entries, list):
                continue
            for entry in entries:
                if not isinstance(entry, dict):
                    continue
                ma, cut = _entry_to_pair(entry)
                if ma is not None and cut is not None:
                    out[ma] = cut
            if out:
                break

    if not out and isinstance(data, dict):
        for key, value in data.items():
            ma_key = _parse_ma(str(key)) or _to_int(key)
            if isinstance(value, dict):
                ma, cut = _entry_to_pair({"mA": ma_key, **value})
            else:
                ma, cut = ma_key, _to_float(value)
            if ma is not None and cut is not None:
                out[ma] = cut

    return out


def _complete_mva_cuts(cuts: Dict[int, float], target_masses: List[int]) -> Dict[int, float]:
    if not cuts:
        raise ValueError("No valid MVA cut points were parsed from JSON.")

    points = sorted((int(mass), float(cut)) for mass, cut in cuts.items())
    completed = dict(cuts)
    interpolated = []

    for mass in target_masses:
        if mass in completed:
            continue

        # Nearest simulated mass, as the analysis does (apply_bdt_data.py::find_nearest_cut);
        # a linear interpolation drew some panels 0.001-0.004 away from the adopted cut (2026-09-27).
        completed[mass] = min(points, key=lambda p: abs(p[0] - mass))[1]

        interpolated.append((mass, completed[mass]))

    if interpolated:
        msg = ", ".join(f"{mass}:{cut:.4f}" for mass, cut in interpolated)
        print(f"[Info] Nearest-anchor MVA cuts for missing masses -> {msg}")

    return completed


def _resolve_mva_cut_json(user_path: str) -> Path:
    candidates = [Path(user_path), Path(DEFAULT_MVA_CUT_JSON), *LOCAL_MVA_CUT_CANDIDATES]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "Cannot find MVA cut JSON. Tried:\n  - " + "\n  - ".join(str(path) for path in candidates)
    )


def _resolve_output_dir(user_path: str) -> Path:
    candidates = [Path(user_path), Path(DEFAULT_OUTPUT_DIR), LOCAL_OUTPUT_DIR]
    tried = []
    for candidate in candidates:
        if candidate in tried:
            continue
        tried.append(candidate)
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            return candidate
        except Exception as exc:
            print(f"[Warning] Cannot create output dir {candidate}: {exc}")
    raise OSError("Failed to create any output directory.")


def _resolve_sig_eff_json(channel_mode: str) -> Path:
    if channel_mode == "ele":
        candidates = [Path(DEFAULT_SIGEFF_ELE_JSON), *LOCAL_SIGEFF_ELE_CANDIDATES]
    elif channel_mode == "mu":
        candidates = [Path(DEFAULT_SIGEFF_MU_JSON), *LOCAL_SIGEFF_MU_CANDIDATES]
    else:
        raise ValueError(f"Unsupported efficiency channel: {channel_mode}")

    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"Cannot find signal efficiency JSON for channel={channel_mode}. Tried:\n  - "
        + "\n  - ".join(str(path) for path in candidates)
    )


def _load_sig_eff_payload(channel_mode: str) -> dict:
    if channel_mode not in _EFF_CACHE:
        with open(_resolve_sig_eff_json(channel_mode), "r") as handle:
            _EFF_CACHE[channel_mode] = json.load(handle)
    return _EFF_CACHE[channel_mode]


def _signal_efficiency_for_mass(mass: int, channel_mode: str, years: List[str]) -> float:
    def _weighted(channel_key: str) -> float:
        payload = _load_sig_eff_payload(channel_key)
        total = 0.0
        for year in years:
            total += LUMI_MAP[year] * float(payload["values"][year][str(int(mass))])
        return total

    if channel_mode == "ele":
        return _weighted("ele")
    if channel_mode == "mu":
        return _weighted("mu")
    if channel_mode == "inclusive":
        return _weighted("ele") + _weighted("mu")
    raise ValueError(f"Unsupported plot channel mode: {channel_mode}")


def _resolve_mva_branch_for_mass(chain, mass_tag: str) -> Optional[str]:
    if not chain or not chain.GetListOfBranches():
        return None

    branches = chain.GetListOfBranches()
    candidates = [
        f"MVA_Score_mA_{mass_tag}",
        f"MVA_Score_{mass_tag}",
        f"MVA_score_{mass_tag}",
        f"BDTScore_{mass_tag}",
        f"BDT_{mass_tag}",
        "MVA_Score",
    ]
    for name in candidates:
        if branches.FindObject(name):
            return name
    return None


def _signal_sample_name(ma: int) -> Optional[str]:
    return f"M{ma}" if ma in SIGNAL_MASSES else None


def _signal_color(sample_name: Optional[str]) -> int:
    if not sample_name:
        return ROOT.kMagenta + 1
    color_hex = SIGNAL_COLOR_HEX.get(sample_name)
    if color_hex:
        return TColor.GetColor(color_hex)
    return ROOT.kMagenta + 1


def _nearest_anchor_mass(mass: int) -> int:
    return min(SIGNAL_MASSES, key=lambda anchor: abs(anchor - int(mass)))


def _find_bracketing_anchors(mass: int) -> Tuple[int, int]:
    target = int(mass)
    anchors = sorted(SIGNAL_MASSES)
    if target <= anchors[0]:
        return anchors[0], anchors[0]
    if target >= anchors[-1]:
        return anchors[-1], anchors[-1]
    for left, right in zip(anchors[:-1], anchors[1:]):
        if left <= target <= right:
            return left, right
    raise ValueError(f"Cannot find bracketing anchors for mA={mass}")


def _resolve_signal_components(mass: int, shape_mode: str = "mixture") -> List[Dict[str, float]]:
    target = int(mass)
    if shape_mode == "nearest":
        return [{"anchor_mass": _nearest_anchor_mass(target), "shape_weight": 1.0}]
    if target in SIGNAL_MASSES:
        return [{"anchor_mass": target, "shape_weight": 1.0}]

    left, right = _find_bracketing_anchors(target)
    if left == right:
        return [{"anchor_mass": left, "shape_weight": 1.0}]

    weight_left = float(right - target) / float(right - left)
    weight_right = float(target - left) / float(right - left)
    return [
        {"anchor_mass": left, "shape_weight": weight_left},
        {"anchor_mass": right, "shape_weight": weight_right},
    ]


def _build_signal_overlay(
    mass: int,
    histos: Dict[str, TH1F],
    all_histos: Dict[int, Dict[str, TH1F]],
    plot_cfg: Plot_Config,
    channel_mode: str,
    years: List[str],
    shape_mode: str = "mixture",
    signal_scale: float = SIGNAL_DRAW_SCALE,
) -> Tuple[Optional[TH1F], Optional[str]]:
    components = _resolve_signal_components(mass, shape_mode=shape_mode)
    signal_hist = None

    for idx, component in enumerate(components):
        anchor_mass = int(component["anchor_mass"])
        weight = float(component["shape_weight"])
        sample = _signal_sample_name(anchor_mass)
        if not sample or sample not in histos:
            continue

        source_hist = histos[sample]
        if source_hist.GetEntries() <= 0 and source_hist.Integral() <= 0.0:
            continue

        temp_hist = source_hist.Clone(f"h_sig_component_mA_{mass}_{sample}_{idx}")
        temp_hist.SetDirectory(0)
        temp_hist.Scale(weight)

        if signal_hist is None:
            signal_hist = temp_hist.Clone(f"h_sig_draw_mA_{mass}")
            signal_hist.SetDirectory(0)
        else:
            signal_hist.Add(temp_hist)

    if signal_hist is None:
        return None, None

    nearest_anchor = _nearest_anchor_mass(mass)
    nearest_sample = _signal_sample_name(nearest_anchor)
    signal_hist.SetLineColor(_signal_color(nearest_sample))
    signal_hist.SetLineWidth(4)
    signal_hist.SetFillStyle(0)
    signal_hist.SetFillColor(0)

    ref_hist = all_histos.get(nearest_anchor, {}).get(nearest_sample) if nearest_sample else None
    ref_yield = ref_hist.Integral() if ref_hist else 0.0
    eff_target = _signal_efficiency_for_mass(mass, channel_mode, years)
    eff_nearest = _signal_efficiency_for_mass(nearest_anchor, channel_mode, years)
    eff_scale = (eff_target / eff_nearest) if eff_nearest > 0.0 else 0.0
    target_yield = ref_yield * eff_scale * signal_scale
    current_yield = signal_hist.Integral()
    if current_yield > 0.0:
        signal_hist.Scale(target_yield / current_yield)
    else:
        signal_hist.Scale(0.0)

    if shape_mode == "mixture" and len(components) > 1:
        signal_hist.SetLineStyle(2)
    else:
        signal_hist.SetLineStyle(1)
    legend_label = f"m_{{a}} = {mass} GeV #times {signal_scale:.1f}"

    return signal_hist, legend_label


def _build_uniform_bin_edges(xmin: float, xmax: float, step: float) -> List[float]:
    n_steps = int(round((xmax - xmin) / step))
    return [xmin + step * idx for idx in range(n_steps + 1)]


def _build_n_uniform_edges(xmin: float, xmax: float, nbins: int) -> List[float]:
    width = (xmax - xmin) / float(nbins)
    return [xmin + width * idx for idx in range(nbins + 1)]


def _format_elapsed_hms(elapsed_seconds: float) -> str:
    total_seconds = int(round(elapsed_seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}h {minutes:02d}m {seconds:02d}s"


def _sideband_integral(hist, signal_low: float = BLIND_LOW, signal_high: float = BLIND_HIGH) -> float:
    axis = hist.GetXaxis()
    nbins = axis.GetNbins()
    low_bin = axis.FindFixBin(signal_low)
    high_bin = axis.FindFixBin(signal_high)

    left = hist.Integral(1, low_bin - 1) if low_bin > 1 else 0.0
    right = hist.Integral(high_bin + 1, nbins) if high_bin < nbins else 0.0
    return left + right


def _scale_bkg_to_data_sideband(histos: Dict[str, TH1F], analyzer_cfg: Analyzer_Config) -> float:
    data_hist = histos.get("Data")
    if not data_hist:
        print("[Scale] Data histogram missing. Skip sideband normalization.")
        return 1.0

    data_sb = _sideband_integral(data_hist)
    if data_sb <= 0.0:
        print("[Scale] Data sideband yield <= 0. Skip sideband normalization.")
        return 1.0

    bkg_sb = 0.0
    for sample in analyzer_cfg.bkg_names:
        if sample in histos:
            bkg_sb += _sideband_integral(histos[sample])

    if bkg_sb <= 0.0:
        print("[Scale] Background sideband yield <= 0. Skip sideband normalization.")
        return 1.0

    scale = data_sb / bkg_sb
    for sample in analyzer_cfg.bkg_names:
        if sample in histos:
            histos[sample].Scale(scale)

    return scale


def _clone_total_bkg(histos: Dict[str, TH1F], analyzer_cfg: Analyzer_Config, name: str):
    total = None
    for sample in analyzer_cfg.bkg_names:
        hist = histos.get(sample)
        if not hist:
            continue
        if total is None:
            total = hist.Clone(name)
            total.SetDirectory(0)
        else:
            total.Add(hist)
    return total


def _style_histograms(histos: Dict[str, TH1F], plot_cfg: Plot_Config, analyzer_cfg: Analyzer_Config) -> None:
    for sample, hist in histos.items():
        if sample == "Data":
            hist.SetMarkerStyle(20)
            hist.SetMarkerSize(DATA_MARKER_SIZE)
            hist.SetLineColor(ROOT.kBlack)
            hist.SetLineWidth(2)
        elif sample in analyzer_cfg.bkg_names:
            hist.SetFillColor(plot_cfg.colors[sample])
            hist.SetLineColor(plot_cfg.colors[sample])
            hist.SetLineWidth(1)
        elif sample in analyzer_cfg.sig_names:
            hist.SetLineColor(_signal_color(sample))
            hist.SetLineWidth(4)
            hist.SetFillStyle(0)
            hist.SetFillColor(0)


def _make_stack(histos: Dict[str, TH1F], analyzer_cfg: Analyzer_Config, name: str) -> THStack:
    stack = THStack(f"h_stack_{name}", name)
    for sample in analyzer_cfg.bkg_names:
        if sample in histos:
            stack.Add(histos[sample])
    return stack


def _configure_needed_branches(chain, needed_branches: List[str]) -> None:
    if not chain:
        return
    try:
        chain.SetBranchStatus("*", 0)
        for branch in needed_branches:
            if branch:
                chain.SetBranchStatus(branch, 1)
    except Exception as exc:
        print(f"[Warning] Failed to restrict branches: {exc}")


def _draw_mass_plot(
    mass: int,
    histos: Dict[str, TH1F],
    all_histos: Dict[int, Dict[str, TH1F]],
    analyzer_cfg: Analyzer_Config,
    plot_cfg: Plot_Config,
    output_dir: Path,
    logy: bool = False,
    show_signal: bool = True,
    name_suffix: str = "",
    channel_mode: str = "inclusive",
    signal_shape_mode: str = "mixture",
    signal_scale: float = SIGNAL_DRAW_SCALE,
    signal_source_histos: Optional[Dict[str, TH1F]] = None,
    signal_source_all_histos: Optional[Dict[int, Dict[str, TH1F]]] = None,
    selection_label: Optional[str] = None,
) -> None:
    draw_histos: Dict[str, TH1F] = {}
    draw_tag = name_suffix if name_suffix else "_nominal"
    for sample, hist in histos.items():
        cloned = hist.Clone(f"{hist.GetName()}_draw_{mass}{draw_tag}_{sample}")
        cloned.SetDirectory(0)
        draw_histos[sample] = cloned

    _style_histograms(draw_histos, plot_cfg, analyzer_cfg)

    scale = _scale_bkg_to_data_sideband(draw_histos, analyzer_cfg)
    if abs(scale - 1.0) > 1e-12:
        print(f"[mA={mass:02d}] Applied sideband scale factor = {scale:.4f}")

    bkg_total = _clone_total_bkg(draw_histos, analyzer_cfg, f"h_bkg_total_mA_{mass}")
    if bkg_total is None or bkg_total.Integral() <= 0.0:
        print(f"[mA={mass:02d}] Background is empty after MVA cut. Skip plot.")
        return

    stack = _make_stack(draw_histos, analyzer_cfg, f"mH_mA_{mass}")
    ratio = MakeRatioPlot(draw_histos["Data"], bkg_total, f"H_m_mA_{mass}")
    ratio.SetMarkerSize(DATA_MARKER_SIZE)
    stat_abs, stat_norm = Get_StatUnc(bkg_total)

    signal_hist = None
    signal_legend_label = None
    if show_signal:
        signal_input_histos = signal_source_histos if signal_source_histos is not None else draw_histos
        signal_input_all_histos = signal_source_all_histos if signal_source_all_histos is not None else all_histos
        signal_hist, signal_legend_label = _build_signal_overlay(
            mass,
            signal_input_histos,
            signal_input_all_histos,
            plot_cfg,
            channel_mode,
            getattr(analyzer_cfg, "years_sig", []),
            shape_mode=signal_shape_mode,
            signal_scale=signal_scale,
        )
        if signal_hist:
            print(
                f"[mA={mass:02d}] signal draw integral={signal_hist.Integral():.3f}"
                + (f" label={signal_legend_label}" if signal_legend_label else "")
            )

    canvas = TCanvas(f"canv_mA_{mass}{'_log' if logy else ''}", "", 800, 800)
    canvas.SetBottomMargin(0.012)
    canvas.cd()

    upper_pad = TPad(f"upperpad_mA_{mass}", "", 0, 0.25, 1, 1)
    upper_pad.SetLeftMargin(0.16)
    upper_pad.SetRightMargin(0.05)
    upper_pad.SetBottomMargin(0.19)
    upper_pad.SetTopMargin(0.085)
    upper_pad.Draw()
    upper_pad.cd()
    upper_pad.SetTickx(1)
    upper_pad.SetTicky(1)

    peak_ref = bkg_total.GetMaximum()
    if signal_hist:
        peak_ref = max(peak_ref, signal_hist.GetMaximum())
    if peak_ref <= 0.0:
        peak_ref = 1.0
    ymax = peak_ref * 1.8

    draw_histos["Data"].SetMinimum(1e-2 if logy else 0.0)
    draw_histos["Data"].SetMaximum(ymax)
    draw_histos["Data"].GetXaxis().SetRangeUser(H_M_XMIN, H_M_XMAX)
    if logy:
        upper_pad.SetLogy()

    draw_histos["Data"].Draw("PE")
    draw_histos["Data"].GetXaxis().SetLabelSize(0.0)
    draw_histos["Data"].GetXaxis().SetTitleOffset(0.95)
    draw_histos["Data"].GetYaxis().SetTitle(f"Events / {draw_histos['Data'].GetBinWidth(1):.2f} GeV")
    draw_histos["Data"].GetYaxis().SetTitleSize(0.07)
    draw_histos["Data"].GetYaxis().SetLabelSize(0.055)
    draw_histos["Data"].GetYaxis().SetTitleFont(42)
    draw_histos["Data"].GetYaxis().SetTitleOffset(1.15)

    stack.SetMinimum(1e-2 if logy else 0.0)
    stack.SetMaximum(ymax)
    stack.Draw("HIST SAME")

    Draw_unc(stat_abs, TColor.GetColor("#404040"), alpha=0.90, fill_style=3354)

    if signal_hist:
        signal_hist.Draw("HIST SAME")

    draw_histos["Data"].Draw("PE SAME")
    draw_histos["Data"].Draw("AXIS SAME")

    legend = TLegend(0.58, 0.60, 0.95, 0.87)
    ROOT.SetOwnership(legend, False)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetFillColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.AddEntry(draw_histos["Data"], "Data", "PE")
    for sample in analyzer_cfg.bkg_names:
        if sample in draw_histos:
            legend.AddEntry(draw_histos[sample], BKG_LABELS.get(sample, sample), "f")
    legend.AddEntry(stat_abs, "Stat. Unc.", "f")
    if signal_hist and signal_legend_label:
        legend.AddEntry(signal_hist, signal_legend_label, "l")
    legend.Draw("SAME")

    tag = ROOT.TLatex()
    tag.SetNDC()
    tag.SetTextFont(42)
    tag.SetTextSize(0.045)
    tag.DrawLatex(0.20, 0.84, f"After BDT Cut m_{{a}} = {mass} GeV")
    if selection_label:
        tag.SetTextSize(0.038)
        tag.DrawLatex(0.20, 0.79, selection_label)

    CMS_lumi.cmsText = "CMS"
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.cmsTextSize = 0.95
    CMS_lumi.CMSText_posX = -0.03
    CMS_lumi.outOfFrame = True
    CMS_lumi.lumiText_posX = -0.000
    CMS_lumi.CMS_lumi(canvas, 5, 0, plot_cfg.year)

    canvas.cd()
    lower_pad = TPad(f"lowerpad_mA_{mass}", "", 0, 0, 1, 0.35)
    lower_pad.SetTopMargin(0.00001)
    lower_pad.SetBottomMargin(0.36)
    lower_pad.SetLeftMargin(0.16)
    lower_pad.SetRightMargin(0.05)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    lower_pad.SetTickx(1)
    lower_pad.SetTicky(1)

    ratio.GetXaxis().SetLimits(H_M_XMIN, H_M_XMAX)
    ratio.GetXaxis().SetTitle(plot_cfg.var_title_map["H_m"])
    ratio.GetXaxis().SetTitleOffset(1.0)
    ratio.Draw("APZ SAME")
    Draw_unc(stat_norm, TColor.GetColor("#404040"), alpha=0.90, fill_style=3354)
    ratio.Draw("PZ SAME")

    save_name = f"bkgmcScupltingCheck_mA{mass:02d}{name_suffix}"
    if logy:
        save_name += "_log"
    SaveCanvPic(canvas, str(output_dir), save_name)


def _book_histograms(
    analyzer_cfg: Analyzer_Config,
    sample_names: Optional[List[str]] = None,
    bin_edges: Optional[List[float]] = None,
    nbins: int = H_M_NBINS,
    xmin: float = H_M_XMIN,
    xmax: float = H_M_XMAX,
) -> Dict[int, Dict[str, TH1F]]:
    histos: Dict[int, Dict[str, TH1F]] = {}
    samples = sample_names or analyzer_cfg.samp_names
    root_edges = array("d", bin_edges) if bin_edges is not None else None
    for mass in TARGET_MASSES:
        histos[mass] = {}
        for sample in samples:
            if root_edges is not None:
                hist = TH1F(
                    f"H_m_mA_{mass}_{sample}",
                    f"H_m_mA_{mass}_{sample}",
                    len(root_edges) - 1,
                    root_edges,
                )
            else:
                hist = TH1F(
                    f"H_m_mA_{mass}_{sample}",
                    f"H_m_mA_{mass}_{sample}",
                    nbins,
                    xmin,
                    xmax,
                )
            hist.Sumw2()
            hist.SetDirectory(0)
            histos[mass][sample] = hist
    return histos


def _book_bdt_shape_histograms(
    n_bdt_bins: int = BDT_SHAPE_NBINS,
    score_min: float = BDT_SHAPE_SCORE_MIN,
    score_max: float = BDT_SHAPE_SCORE_MAX,
    mass_bin_width: float = BDT_SHAPE_MASS_BIN_WIDTH,
    xmin: float = H_M_XMIN,
    xmax: float = H_M_XMAX,
) -> Tuple[Dict[int, List[TH1F]], List[float]]:
    score_edges = _build_n_uniform_edges(score_min, score_max, n_bdt_bins)
    mass_nbins = int(round((xmax - xmin) / mass_bin_width))
    histos: Dict[int, List[TH1F]] = {}
    for mass in TARGET_MASSES:
        histos[mass] = []
        for idx in range(n_bdt_bins):
            hist = TH1F(
                f"h_bkg_mass_shape_mA_{mass}_bdt_{idx}",
                f"h_bkg_mass_shape_mA_{mass}_bdt_{idx}",
                mass_nbins,
                xmin,
                xmax,
            )
            hist.Sumw2()
            hist.SetDirectory(0)
            histos[mass].append(hist)
    return histos, score_edges


def _draw_bkg_mass_shapes_by_bdt(
    mass: int,
    histos: List[TH1F],
    score_edges: List[float],
    plot_cfg: Plot_Config,
    output_dir: Path,
    logy: bool = False,
    selection_label: Optional[str] = None,
) -> None:
    draw_histos = []
    ymax = 0.0

    for idx, hist in enumerate(histos):
        if hist.Integral() <= 0.0:
            continue
        draw_hist = hist.Clone(f"{hist.GetName()}_norm")
        draw_hist.SetDirectory(0)
        norm = draw_hist.Integral("width")
        if norm > 0.0:
            draw_hist.Scale(1.0 / norm)

        color = BDT_SHAPE_COLORS[idx % len(BDT_SHAPE_COLORS)]
        line_style = BDT_SHAPE_LINE_STYLE_OVERRIDES.get(
            idx, BDT_SHAPE_LINE_STYLES[idx % len(BDT_SHAPE_LINE_STYLES)]
        )
        draw_hist.SetLineColor(color)
        draw_hist.SetLineStyle(line_style)
        draw_hist.SetLineWidth(3)
        draw_hist.SetMarkerSize(0)
        draw_hist.SetFillStyle(0)
        ymax = max(ymax, draw_hist.GetMaximum())
        draw_histos.append((idx, draw_hist))

    if not draw_histos:
        print(f"[mA={mass:02d}] No background entries for BDT-binned mass-shape plot.")
        return

    y_axis_max = 1.6 * max(ymax, 1e-6)
    y_axis_min = 1e-5 if logy else 0.0

    canvas = TCanvas(f"canv_bkg_mass_shapes_by_bdt_mA_{mass}", "", 800, 600)
    canvas.SetLeftMargin(0.16)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.14)
    canvas.SetTopMargin(0.08)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    if logy:
        canvas.SetLogy()

    frame = TH1F(
        f"frame_bkg_mass_shapes_by_bdt_mA_{mass}",
        "",
        1,
        H_M_XMIN,
        H_M_XMAX,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(y_axis_min)
    frame.SetMaximum(y_axis_max)
    frame.GetXaxis().SetTitle(plot_cfg.var_title_map["H_m"])
    frame.GetXaxis().SetTitleSize(0.055)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.25)
    frame.GetYaxis().SetTitle(f"A.U. / {draw_hist.GetXaxis().GetBinWidth(1):.2f} GeV")
    frame.GetYaxis().SetTitleSize(0.055)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.55)
    frame.Draw("AXIS")

    legend = TLegend(0.69, 0.45, 0.99, 0.84)
    ROOT.SetOwnership(legend, False)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetFillColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.032)

    for idx, hist in draw_histos:
        hist.Draw("HIST SAME")
        legend.AddEntry(hist, f"{score_edges[idx]:.1f} < BDT #leq {score_edges[idx + 1]:.1f}", "l")

    vertical_lines = []
    for x_value in (BLIND_LOW, BLIND_HIGH):
        line = ROOT.TLine(x_value, y_axis_min, x_value, y_axis_max)
        line.SetLineColor(ROOT.kBlue)
        line.SetLineStyle(2)
        line.SetLineWidth(3)
        line.Draw("SAME")
        vertical_lines.append(line)

    legend.Draw("SAME")

    tag = ROOT.TLatex()
    tag.SetNDC()
    tag.SetTextFont(42)
    tag.SetTextSize(0.035)
    tag.DrawLatex(0.765, 0.86, f"m_{{a}} = {mass} GeV")
    if selection_label:
        tag.SetTextSize(0.030)
        tag.DrawLatex(0.705, 0.81, selection_label)

    cms_label = ROOT.TLatex()
    cms_label.SetNDC()
    cms_label.SetTextAlign(11)
    cms_label.SetTextFont(61)
    cms_label.SetTextSize(BDT_SHAPE_CMS_TEXT_SIZE)
    cms_label.DrawLatex(0.16, 0.935, "CMS")

    prelim_label = ROOT.TLatex()
    prelim_label.SetNDC()
    prelim_label.SetTextAlign(11)
    prelim_label.SetTextFont(52)
    prelim_label.SetTextSize(BDT_SHAPE_PRELIM_TEXT_SIZE)
    prelim_label.DrawLatex(0.24, 0.935, "Preliminary")

    lumi_label = ROOT.TLatex()
    lumi_label.SetNDC()
    lumi_label.SetTextAlign(31)
    lumi_label.SetTextFont(42)
    lumi_label.SetTextSize(BDT_SHAPE_LUMI_TEXT_SIZE)
    lumi_label.DrawLatex(0.95, 0.935, BDT_SHAPE_LUMI_TEXT)

    save_name = f"bkg_mass_shapes_by_bdt_mA{mass:02d}"
    if logy:
        save_name += "_log"
    SaveCanvPic(canvas, str(output_dir), save_name)


def _chain_files(chain) -> List[str]:
    return [str(element.GetTitle()) for element in chain.GetListOfFiles()]


def _fill_n(hist, x, w) -> None:
    """TH1::FillN == TH1::Fill called once per entry in the same order (same float/Sumw2/stat sums)."""
    if len(x):
        hist.FillN(len(x), np.ascontiguousarray(x, dtype=np.float64), np.ascontiguousarray(w, dtype=np.float64))


def _bdt_shape_bin_indices(scores: np.ndarray, score_edges: List[float]) -> np.ndarray:
    """BDT-score bin index per entry (last bin closed on the right); -1 outside [min, max]."""
    n_bins = len(score_edges) - 1
    width = score_edges[1] - score_edges[0]
    out = np.full(len(scores), -1, dtype=np.int64)
    in_range = ~((scores < score_edges[0]) | (scores > score_edges[-1]))
    at_top = scores == score_edges[-1]
    other = in_range & ~at_top
    if np.isnan(scores[other]).any():
        raise ValueError("cannot convert float NaN to integer")
    out[other] = ((scores[other] - score_edges[0]) / width).astype(np.int64)
    out[at_top] = n_bins - 1
    out[(out < 0) | (out >= n_bins)] = -1
    return out


def _sideband_weights_as_seen_by_chain(sideband_reweighter, chain, enabled: set) -> np.ndarray:
    """Per-entry sideband weights equal to sideband_reweighter.weight_for_object(chain, row_index=i).

    The per-event loop read the chain after SetBranchStatus("*", 0): a branch that exists but was
    not re-enabled is seen by getattr as 0.0, a missing one as NaN (factor 1).
    The frame below reproduces exactly those values, then weights_for_dataframe does the reweighting.
    """
    step_vars = list(dict.fromkeys(
        step["var"] for iteration in sideband_reweighter.iterations for step in iteration.get("steps", [])
    ))
    frames = []
    offset = 0
    for path in _chain_files(chain):
        with uproot.open(path) as handle:
            tree = handle[chain.GetName()]
            n_rows = tree.num_entries
            names = set(tree.keys())

            def _seen(name, dtype=np.float64):
                if name not in names:
                    return None
                if name not in enabled:
                    return np.zeros(n_rows, dtype=dtype)
                return tree[name].array(library="np").astype(dtype)

            def _first(logical, dtype=np.float64):
                for candidate in VAR_ALIASES.get(logical, (logical,)):
                    values = _seen(candidate, dtype)
                    if values is not None:
                        return values
                return None

            def _finite(values):
                if values is None:
                    return np.full(n_rows, np.nan)
                return np.where(np.isfinite(values), values, np.nan)

            columns = {}
            for var in step_vars:
                if var == "param":
                    param = _first("param")
                    if param is not None:
                        columns[var] = _finite(param)
                        continue
                    h_m, alp_m = _first("H_m"), _first("ALP_m")
                    if h_m is None or alp_m is None:
                        columns[var] = np.full(n_rows, np.nan)
                        continue
                    events = _first("event", np.int64)
                    if sideband_reweighter.param_method == "event_hash" and events is not None:
                        mass_hyp = _mass_from_event(events, sideband_reweighter.seed, sideband_reweighter.mass_hypotheses)
                    else:
                        mass_hyp = _mass_from_row_order(
                            offset + n_rows, sideband_reweighter.seed, sideband_reweighter.mass_hypotheses
                        )[offset:]
                    bad = ~np.isfinite(h_m) | (h_m == 0.0) | ~np.isfinite(alp_m)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        columns[var] = np.where(bad, np.nan, (alp_m - mass_hyp) / h_m)
                elif var in DERIVED_VARS:
                    if var == "pho_pt_asym":
                        pt1, pt2 = _finite(_first("pho1Pt_oHm")), _finite(_first("pho2Pt_oHm"))
                        with np.errstate(divide="ignore", invalid="ignore"):
                            columns[var] = np.where(
                                np.isfinite(pt1) & np.isfinite(pt2), (pt1 - pt2) / (pt1 + pt2 + 1e-6), np.nan
                            )
                    else:
                        columns[var] = np.full(n_rows, np.nan)
                else:
                    columns[var] = _finite(_first(var))
            frames.append(pd.DataFrame(columns, index=np.arange(offset, offset + n_rows)))
            offset += n_rows
    if not frames:
        return np.ones(0)
    return sideband_reweighter.weights_for_dataframe(pd.concat(frames))


def _read_chain_columns(chain, defaults: Dict[str, float]) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Concatenate the chain's files in chain order; a branch missing in a file takes its getattr default
    (NaN, flagged in the returned 'present' mask, when the per-event loop would have skipped it)."""
    parts: Dict[str, List[np.ndarray]] = {name: [] for name in defaults}
    present: Dict[str, List[np.ndarray]] = {name: [] for name in defaults}
    for path in _chain_files(chain):
        with uproot.open(path) as handle:
            tree = handle[chain.GetName()]
            n_rows = tree.num_entries
            names = set(tree.keys())
            to_read = [name for name in defaults if name in names]
            arrays = tree.arrays(to_read, library="np") if to_read else {}
            for name, default in defaults.items():
                if name in arrays:
                    parts[name].append(np.asarray(arrays[name]).astype(np.float64))
                    present[name].append(np.ones(n_rows, dtype=bool))
                else:
                    parts[name].append(np.full(n_rows, np.nan if default is None else default, dtype=np.float64))
                    present[name].append(np.full(n_rows, default is not None, dtype=bool))
    columns = {name: (np.concatenate(chunks) if chunks else np.zeros(0)) for name, chunks in parts.items()}
    masks = {name: (np.concatenate(chunks) if chunks else np.zeros(0, dtype=bool)) for name, chunks in present.items()}
    return columns, masks


def _fill_histograms(
    ntuples,
    analyzer_cfg: Analyzer_Config,
    mva_cuts: Dict[int, float],
    histos: Dict[int, Dict[str, TH1F]],
    extra_histos: Optional[List[Dict[int, Dict[str, TH1F]]]],
    bdt_shape_histos: Optional[Dict[int, List[TH1F]]],
    bdt_shape_edges: Optional[List[float]],
    blind: bool,
    only_ele: bool,
    only_mu: bool,
    sideband_reweighter=None,
) -> None:
    """Columnar version of the former per-event GetEvent loop (identical histograms).

    Each histogram is filled with TH1::FillN over the selected entries in chain order, which performs the
    same per-entry float32 bin-content / double Sumw2 accumulation as the former TH1::Fill calls.
    """
    mva_branch_map: Dict[str, Dict[int, Optional[str]]] = {}
    enabled_branches: Dict[str, set] = {}
    extra_histos = extra_histos or []
    for sample in analyzer_cfg.samp_names:
        masses_to_use = TARGET_MASSES if (
            sample in analyzer_cfg.bkg_names or sample == "Data" or sample in analyzer_cfg.sig_names
        ) else []

        branch_map = {}
        for mass in masses_to_use:
            branch_map[mass] = _resolve_mva_branch_for_mass(ntuples[sample], f"M{mass}")
        mva_branch_map[sample] = branch_map

        needed_branches = ["H_m", "weight", "z_mumu", "z_ee"]
        needed_branches.extend(sorted({branch for branch in branch_map.values() if branch}))
        _configure_needed_branches(ntuples[sample], needed_branches)
        enabled_branches[sample] = set(needed_branches)

        missing = [mass for mass, branch in branch_map.items() if branch is None]
        if missing:
            print(f"[Warning] sample={sample} missing MVA branch for mA={missing}")

    # When sideband reweight is enabled, also need the input vars for the reweight steps.
    if sideband_reweighter is not None:
        for sample in analyzer_cfg.bkg_names:
            chain = ntuples.get(sample)
            if not chain:
                continue
            try:
                rwgt_branches = set()
                for iteration in getattr(sideband_reweighter, "iterations", []):
                    for step in iteration.get("steps", []):
                        rwgt_branches.add(step.get("var"))
                # param and ALP_m (with H_m, event) build the mass-hypothesis variable; before 2026-10-31 they
                # were left disabled, so the reweight saw param = 0 (nominal) or (0 - m_hyp)/H_m (control)
                rwgt_branches.update(("pho1ECALIso", "pho2ECALIso", "pho1PIso_noCorr", "pho2PIso_noCorr",
                                      "ALP_lead_photon_ecalPFClusterIso", "ALP_sublead_photon_ecalPFClusterIso",
                                      "H_mass", "ALP_mass", "ALP_m", "event"))
                for branch in list(rwgt_branches):
                    if not branch:
                        continue
                    if chain.GetListOfBranches() and chain.GetListOfBranches().FindObject(branch):
                        chain.SetBranchStatus(branch, 1)
                        enabled_branches[sample].add(branch)
            except Exception as _exc_rw:
                print(f"[Warning] sideband reweight branch enable failed for {sample}: {_exc_rw}")

    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample]
        entries = ntup.GetEntries()
        print(f"[Fill] sample={sample:>10s} entries={entries}")

        masses_to_fill = list(mva_branch_map.get(sample, {}).keys())
        if not masses_to_fill:
            continue

        active = [
            (mass, mva_cuts.get(mass), mva_branch_map[sample].get(mass))
            for mass in masses_to_fill
        ]
        active = [(mass, cut, branch) for mass, cut, branch in active if cut is not None and branch is not None]
        defaults = {"H_m": -999.0, "z_mumu": 0.0, "z_ee": 0.0}
        if sample != "Data":
            defaults["weight"] = 1.0
        defaults.update({branch: None for _, _, branch in active})
        columns, present = _read_chain_columns(ntup, defaults)

        h_mass = columns["H_m"]
        keep = ~((h_mass < H_M_XMIN) | (h_mass > H_M_XMAX))
        if only_ele:
            keep &= ~(np.abs(columns["z_mumu"].astype(np.int64)) == 1)
        if only_mu:
            keep &= ~(np.abs(columns["z_ee"].astype(np.int64)) == 1)

        if sample == "Data":
            weight = np.ones(len(h_mass))
        else:
            weight = columns["weight"]
            if sideband_reweighter is not None and sample in analyzer_cfg.bkg_names:
                weight = weight * _sideband_weights_as_seen_by_chain(
                    sideband_reweighter, ntup, enabled_branches[sample]
                )

        sel = np.nonzero(keep)[0]
        h_sel, w_sel = h_mass[sel], weight[sel]
        blind_sel = (BLIND_LOW <= h_sel) & (h_sel <= BLIND_HIGH)
        fill_bdt_shape = (
            bdt_shape_histos is not None and bdt_shape_edges is not None and sample in analyzer_cfg.bkg_names
        )

        for mass, cut, branch in active:
            has_score = present[branch][sel]
            score = columns[branch][sel]

            if fill_bdt_shape:
                bdt_idx = _bdt_shape_bin_indices(np.where(has_score, score, -1.0), bdt_shape_edges)
                bdt_idx[~has_score] = -1
                for idx, hist in enumerate(bdt_shape_histos[mass]):
                    in_bin = bdt_idx == idx
                    _fill_n(hist, h_sel[in_bin], w_sel[in_bin])

            passed = has_score & ~(score < cut)
            if sample == "Data" and blind:
                passed &= ~blind_sel

            x_pass, w_pass = h_sel[passed], w_sel[passed]
            _fill_n(histos[mass][sample], x_pass, w_pass)
            for hist_set in extra_histos:
                hist = hist_set.get(mass, {}).get(sample)
                if hist is not None:
                    _fill_n(hist, x_pass, w_pass)


def _channel_mode_from_args(only_ele: bool, only_mu: bool) -> str:
    if only_ele:
        return "ele"
    if only_mu:
        return "mu"
    return "inclusive"


def _build_output_dir_map(base_output_dir: Path, skip_bdt_shape_plots: bool) -> Dict[str, Path]:
    output_dirs = {
        "signal_mixture": base_output_dir / "signalShapeMixture",
        "bkg_only": base_output_dir / "bkgOnly",
        "bkg_only_bin5": base_output_dir / "bkgOnly_bin5GeV",
        "signal_nearest_bin5": base_output_dir / "signalNearest_bin5GeV",
        "signal_nearest_bin10": base_output_dir / "signalNearest_bin10GeV",
        "signal_nearest_bin1": base_output_dir / "signalNearest_bin1GeV",
        "bdt_mass_shapes": base_output_dir / "bdtMassShapes",
    }
    output_dirs["signal_mixture"].mkdir(parents=True, exist_ok=True)
    output_dirs["bkg_only"].mkdir(parents=True, exist_ok=True)
    output_dirs["bkg_only_bin5"].mkdir(parents=True, exist_ok=True)
    output_dirs["signal_nearest_bin5"].mkdir(parents=True, exist_ok=True)
    output_dirs["signal_nearest_bin10"].mkdir(parents=True, exist_ok=True)
    output_dirs["signal_nearest_bin1"].mkdir(parents=True, exist_ok=True)
    if not skip_bdt_shape_plots:
        output_dirs["bdt_mass_shapes"].mkdir(parents=True, exist_ok=True)
    return output_dirs


def _print_output_dir_map(label: str, base_output_dir: Path, output_dirs: Dict[str, Path], skip_bdt_shape_plots: bool) -> None:
    print(f"[{label}] Plot directory: {base_output_dir}")
    print(f"[{label}] Signal mixture directory: {output_dirs['signal_mixture']}")
    print(f"[{label}] Bkg-only directory: {output_dirs['bkg_only']}")
    print(f"[{label}] Bkg-only 5 GeV directory: {output_dirs['bkg_only_bin5']}")
    print(f"[{label}] Signal nearest 5 GeV directory: {output_dirs['signal_nearest_bin5']}")
    print(f"[{label}] Signal nearest 10 GeV directory: {output_dirs['signal_nearest_bin10']}")
    if not skip_bdt_shape_plots:
        print(f"[{label}] BDT-binned bkg mass-shape directory: {output_dirs['bdt_mass_shapes']}")


def _build_analyzer_cfg(year: str, plot_output_dir: Path, sample_loc_override: Optional[str] = None) -> Analyzer_Config:
    analyzer_cfg = Analyzer_Config("inclusive", year, 0, True)
    analyzer_cfg.mva = True
    analyzer_cfg.plot_output_path = str(plot_output_dir)
    if sample_loc_override:
        analyzer_cfg.sample_loc = sample_loc_override
    return analyzer_cfg


def _run_plot_suite(
    suite_label: str,
    year: str,
    sample_loc_override: Optional[str],
    base_output_dir: Path,
    mva_cuts: Dict[int, float],
    blind: bool,
    only_ele: bool,
    only_mu: bool,
    logy: bool,
    skip_bdt_shape_plots: bool,
    bdt_shape_bins: int,
    selection_label: Optional[str] = None,
    sideband_reweighter=None,
) -> None:
    output_dirs = _build_output_dir_map(base_output_dir, skip_bdt_shape_plots)
    _print_output_dir_map(suite_label, base_output_dir, output_dirs, skip_bdt_shape_plots)

    analyzer_cfg = _build_analyzer_cfg(year, base_output_dir, sample_loc_override=sample_loc_override)
    plot_cfg = Plot_Config(analyzer_cfg, year)
    print(f"[{suite_label}] Input sample directory: {analyzer_cfg.sample_loc}")

    if selection_label and "control" in selection_label.lower():
        print(f"[{suite_label}] Reusing nominal signal-efficiency JSON for interpolated signal overlays.")

    ntuples = LoadNtuples(analyzer_cfg)
    channel_mode = _channel_mode_from_args(only_ele=only_ele, only_mu=only_mu)

    histos = _book_histograms(analyzer_cfg)
    histos_bkg_only_bin5 = _book_histograms(
        analyzer_cfg,
        sample_names=analyzer_cfg.bkg_names + ["Data"],
        bin_edges=_build_uniform_bin_edges(H_M_XMIN, H_M_BIN_COARSE_XMAX, H_M_BIN_COARSE_WIDTH),
    )
    histos_signal_nearest_bin5 = _book_histograms(
        analyzer_cfg,
        bin_edges=_build_uniform_bin_edges(H_M_XMIN, H_M_BIN_COARSE_XMAX, H_M_BIN_COARSE_WIDTH),
    )
    histos_signal_nearest_bin10 = _book_histograms(
        analyzer_cfg,
        bin_edges=_build_uniform_bin_edges(H_M_XMIN, H_M_BIN_COARSE10_XMAX, H_M_BIN_COARSE10_WIDTH),
    )
    histos_signal_nearest_bin1 = _book_histograms(
        analyzer_cfg,
        bin_edges=_build_uniform_bin_edges(H_M_XMIN, H_M_XMAX, 1.0),
    )
    histos_signal_nearest_overlay_2gev = _book_histograms(
        analyzer_cfg,
        sample_names=analyzer_cfg.sig_names,
        bin_edges=_build_uniform_bin_edges(
            H_M_XMIN,
            H_M_BIN_SIGNAL_OVERLAY_XMAX,
            H_M_BIN_SIGNAL_OVERLAY_WIDTH,
        ),
    )

    if skip_bdt_shape_plots:
        bdt_shape_histos = None
        bdt_shape_edges = None
    else:
        bdt_shape_histos, bdt_shape_edges = _book_bdt_shape_histograms(n_bdt_bins=bdt_shape_bins)

    _fill_histograms(
        ntuples=ntuples,
        analyzer_cfg=analyzer_cfg,
        mva_cuts=mva_cuts,
        histos=histos,
        extra_histos=[histos_bkg_only_bin5, histos_signal_nearest_bin5, histos_signal_nearest_bin10, histos_signal_nearest_bin1, histos_signal_nearest_overlay_2gev],
        bdt_shape_histos=bdt_shape_histos,
        bdt_shape_edges=bdt_shape_edges,
        blind=blind,
        only_ele=only_ele,
        only_mu=only_mu,
        sideband_reweighter=sideband_reweighter,
    )

    for mass in TARGET_MASSES:
        print(
            f"[{suite_label}][mA={mass:02d}] yields: "
            f"Data={histos[mass]['Data'].Integral():.3f}, "
            f"DYJetsToLL={histos[mass]['DYJetsToLL'].Integral():.3f}, "
            f"DYGto2LG={histos[mass]['DYGto2LG'].Integral():.3f}"
        )
        signal_sample = _signal_sample_name(mass)
        if signal_sample:
            print(
                f"[{suite_label}][mA={mass:02d}] signal "
                f"{signal_sample}={histos[mass][signal_sample].Integral():.3f}"
            )

        _draw_mass_plot(
            mass=mass,
            histos=histos[mass],
            all_histos=histos,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dirs["signal_mixture"],
            logy=logy,
            show_signal=True,
            name_suffix="",
            channel_mode=channel_mode,
            signal_shape_mode="mixture",
            signal_scale=SIGNAL_DRAW_SCALE,
            selection_label=selection_label,
        )
        _draw_mass_plot(
            mass=mass,
            histos=histos[mass],
            all_histos=histos,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dirs["bkg_only"],
            logy=logy,
            show_signal=False,
            name_suffix="_bkgOnly",
            channel_mode=channel_mode,
            selection_label=selection_label,
        )
        _draw_mass_plot(
            mass=mass,
            histos=histos_bkg_only_bin5[mass],
            all_histos=histos_bkg_only_bin5,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dirs["bkg_only_bin5"],
            logy=logy,
            show_signal=False,
            name_suffix="_bkgOnly_bin5GeV",
            channel_mode=channel_mode,
            selection_label=selection_label,
        )
        _draw_mass_plot(
            mass=mass,
            histos=histos_signal_nearest_bin5[mass],
            all_histos=histos_signal_nearest_bin5,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dirs["signal_nearest_bin5"],
            logy=logy,
            show_signal=True,
            name_suffix="_signalNearest_bin5GeV",
            channel_mode=channel_mode,
            signal_shape_mode="nearest",
            signal_scale=SIGNAL_DRAW_SCALE_NEAREST_BIN5,
            signal_source_histos=histos_signal_nearest_overlay_2gev[mass],
            signal_source_all_histos=histos_signal_nearest_overlay_2gev,
            selection_label=selection_label,
        )
        _draw_mass_plot(
            mass=mass,
            histos=histos_signal_nearest_bin10[mass],
            all_histos=histos_signal_nearest_bin10,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dirs["signal_nearest_bin10"],
            logy=logy,
            show_signal=True,
            name_suffix="_signalNearest_bin10GeV",
            channel_mode=channel_mode,
            signal_shape_mode="nearest",
            signal_scale=SIGNAL_DRAW_SCALE_NEAREST_BIN5,
            signal_source_histos=histos_signal_nearest_overlay_2gev[mass],
            signal_source_all_histos=histos_signal_nearest_overlay_2gev,
            selection_label=selection_label,
        )
        _draw_mass_plot(
            mass=mass,
            histos=histos_signal_nearest_bin1[mass],
            all_histos=histos_signal_nearest_bin1,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dirs["signal_nearest_bin1"],
            logy=logy,
            show_signal=True,
            name_suffix="_signalNearest_bin1GeV",
            channel_mode=channel_mode,
            signal_shape_mode="nearest",
            signal_scale=SIGNAL_DRAW_SCALE_NEAREST_BIN5,
            signal_source_histos=histos_signal_nearest_overlay_2gev[mass],
            signal_source_all_histos=histos_signal_nearest_overlay_2gev,
            selection_label=selection_label,
        )
        if bdt_shape_histos is not None and bdt_shape_edges is not None:
            _draw_bkg_mass_shapes_by_bdt(
                mass=mass,
                histos=bdt_shape_histos[mass],
                score_edges=bdt_shape_edges,
                plot_cfg=plot_cfg,
                output_dir=output_dirs["bdt_mass_shapes"],
                logy=logy,
                selection_label=selection_label,
            )


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Plot background mH sculpting after MVA cuts for mA = 1..30.")
    parser.add_argument("-y", "--Year", dest="year", default="run3", help="Only run3 is supported.")
    parser.add_argument("--mva-cut-json", default=DEFAULT_MVA_CUT_JSON, help="Path to MVA cut JSON.")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Directory for output plots.")
    parser.add_argument(
        "--control-sample-dir",
        default=DEFAULT_CONTROL_SAMPLE_DIR,
        help="Control-sample ntuple directory. If found, an extra control-sample plot set is saved.",
    )
    parser.add_argument(
        "--skip-control-sample-plots",
        action="store_true",
        default=False,
        help="Do not save the extra control-sample plot set.",
    )
    parser.add_argument("--unblind", dest="blind", action="store_false", default=True, help="Show data inside 115-135 GeV.")
    parser.add_argument("-ln", "--ln", dest="ln", action="store_true", default=False, help="Save logY plots.")
    parser.add_argument("--ele", dest="ele", action="store_true", default=False, help="Electron channel only.")
    parser.add_argument("--mu", dest="mu", action="store_true", default=False, help="Muon channel only.")
    parser.add_argument(
        "--skip-bdt-shape-plots",
        action="store_true",
        default=False,
        help="Do not save normalized background H_m shapes in BDT-score bins.",
    )
    parser.add_argument(
        "--bdt-shape-bins",
        type=int,
        default=BDT_SHAPE_NBINS,
        help="Number of uniform BDT-score bins for bkg_mass_shapes_by_bdt plots.",
    )
    parser.add_argument(
        "--use-sideband-reweight",
        dest="use_sideband_reweight",
        action="store_true",
        default=False,
        help="Apply the iterative sideband reweighting JSON to MC background weights.",
    )
    parser.add_argument(
        "--sideband-reweight-json",
        dest="sideband_reweight_json",
        default=None,
        help="Path to sideband reweight JSON. Falls back to HZA_SIDEBAND_REWEIGHT_JSON or the default.",
    )
    args = parser.parse_args()

    if args.year != "run3":
        raise ValueError("plot_bkgmcScupltingCheck.py currently supports only --Year run3.")
    if args.ele and args.mu:
        raise ValueError("Choose at most one of --ele or --mu.")
    if args.bdt_shape_bins <= 0:
        raise ValueError("--bdt-shape-bins must be positive.")

    gROOT.SetBatch(True)
    tdrstyle.setTDRStyle()

    mva_cut_path = _resolve_mva_cut_json(args.mva_cut_json)
    output_dir = _resolve_output_dir(args.output_dir)
    print(f"[Input] MVA cut JSON: {mva_cut_path}")
    mva_cuts = _complete_mva_cuts(_parse_mva_cuts(str(mva_cut_path)), TARGET_MASSES)

    sideband_reweighter = None
    if args.use_sideband_reweight:
        if load_sideband_reweighter is None:
            print("[Warning] sideband_reweight helper not available; skip reweighting.")
        else:
            if args.sideband_reweight_json:
                sideband_reweighter = load_sideband_reweighter(args.sideband_reweight_json, required=True)
            else:
                sideband_reweighter = load_sideband_reweighter()
            if sideband_reweighter is not None:
                print(f"[Reweight] Loaded sideband reweight JSON: {sideband_reweighter.source_path}")
            else:
                print("[Warning] Sideband reweight JSON not found; running without reweighting.")

    _run_plot_suite(
        suite_label="Nominal",
        year=args.year,
        sample_loc_override=None,
        base_output_dir=output_dir,
        mva_cuts=mva_cuts,
        blind=args.blind,
        only_ele=args.ele,
        only_mu=args.mu,
        logy=args.ln,
        skip_bdt_shape_plots=args.skip_bdt_shape_plots,
        bdt_shape_bins=args.bdt_shape_bins,
        selection_label=None,
        sideband_reweighter=sideband_reweighter,
    )

    if not args.skip_control_sample_plots:
        control_sample_dir = Path(args.control_sample_dir).expanduser()
        if control_sample_dir.is_dir():
            _run_plot_suite(
                suite_label="Control",
                year=args.year,
                sample_loc_override=str(control_sample_dir),
                base_output_dir=output_dir / "controlSample",
                mva_cuts=mva_cuts,
                blind=args.blind,
                only_ele=args.ele,
                only_mu=args.mu,
                logy=args.ln,
                skip_bdt_shape_plots=args.skip_bdt_shape_plots,
                bdt_shape_bins=args.bdt_shape_bins,
                selection_label="Control sample",
                sideband_reweighter=sideband_reweighter,
            )
        else:
            print(f"[Warning] Control sample directory not found. Skip extra control plots: {control_sample_dir}")

    elapsed = time.time() - start_time
    print(f"Total runtime: {_format_elapsed_hms(elapsed)}")
    print("Done")


if __name__ == "__main__":
    main()
