import argparse
import json
import re
import sys
import time
from array import array
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import ROOT

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_DIR = SCRIPT_DIR.parent
LIB_DIR = PLOT_DIR / "lib"
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from ROOT import TH1F, THStack, TCanvas, TLegend, TPad, TColor, gROOT

from Analyzer_Configs import Analyzer_Config
from Plot_Configs import Plot_Config
from Plot_Helper import LoadNtuples, MakeRatioPlot, Get_StatUnc, Draw_unc, SaveCanvPic
import CMS_lumi
import tdrstyle


DEFAULT_MVA_CUT_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
DEFAULT_OUTPUT_DIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/bkgmcScupltingCheck"
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
H_M_BIN_SIGNAL_OVERLAY_WIDTH = 2.0
H_M_BIN_SIGNAL_OVERLAY_XMAX = 181.0
BDT_SHAPE_NBINS = 10
BDT_SHAPE_SCORE_MIN = 0.0
BDT_SHAPE_SCORE_MAX = 1.0
BDT_SHAPE_LUMI_TEXT = "170.8 fb^{-1} (13.6 TeV)"
BDT_SHAPE_CMS_TEXT_SIZE = 0.050
BDT_SHAPE_PRELIM_TEXT_SIZE = 0.042
BDT_SHAPE_LUMI_TEXT_SIZE = 0.044
SIGNAL_DRAW_SCALE = 0.5
SIGNAL_DRAW_SCALE_NEAREST_BIN5 = 0.2
LUMI_MAP = {
    "2022preEE": 7.98,
    "2022postEE": 26.70,
    "2023preBPix": 17.79,
    "2023postBPix": 9.45,
    "2024": 108.95,
}
BKG_LABELS = {"DYJetsToLL": "Z + jets", "DYGto2LG": "Z + #gamma"}
_EFF_CACHE: Dict[str, dict] = {}
DATA_MARKER_SIZE = 1.3
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
    ROOT.kMagenta + 1,
]
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

        if mass <= points[0][0]:
            completed[mass] = points[0][1]
        elif mass >= points[-1][0]:
            completed[mass] = points[-1][1]
        else:
            for (m1, c1), (m2, c2) in zip(points[:-1], points[1:]):
                if m1 <= mass <= m2:
                    fraction = (mass - m1) / float(m2 - m1) if m2 > m1 else 0.0
                    completed[mass] = c1 + fraction * (c2 - c1)
                    break

        interpolated.append((mass, completed[mass]))

    if interpolated:
        msg = ", ".join(f"{mass}:{cut:.4f}" for mass, cut in interpolated)
        print(f"[Info] Interpolated MVA cuts for missing masses -> {msg}")

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
    ymax = peak_ref * 1.4

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
    nbins: int = 50,
    xmin: float = H_M_XMIN,
    xmax: float = H_M_XMAX,
) -> Tuple[Dict[int, List[TH1F]], List[float]]:
    score_edges = _build_n_uniform_edges(score_min, score_max, n_bdt_bins)
    histos: Dict[int, List[TH1F]] = {}
    for mass in TARGET_MASSES:
        histos[mass] = []
        for idx in range(n_bdt_bins):
            hist = TH1F(
                f"h_bkg_mass_shape_mA_{mass}_bdt_{idx}",
                f"h_bkg_mass_shape_mA_{mass}_bdt_{idx}",
                nbins,
                xmin,
                xmax,
            )
            hist.Sumw2()
            hist.SetDirectory(0)
            histos[mass].append(hist)
    return histos, score_edges


def _bdt_shape_bin_index(score: float, score_edges: List[float]) -> Optional[int]:
    if score < score_edges[0] or score > score_edges[-1]:
        return None
    if score == score_edges[-1]:
        return len(score_edges) - 2

    width = score_edges[1] - score_edges[0]
    idx = int((score - score_edges[0]) / width)
    if idx < 0 or idx >= len(score_edges) - 1:
        return None
    return idx


def _fill_bdt_shape_histogram(
    sample: str,
    mass: int,
    score: float,
    h_mass: float,
    weight: float,
    analyzer_cfg: Analyzer_Config,
    bdt_shape_histos: Optional[Dict[int, List[TH1F]]],
    bdt_shape_edges: Optional[List[float]],
) -> None:
    if bdt_shape_histos is None or bdt_shape_edges is None:
        return
    if sample not in analyzer_cfg.bkg_names:
        return

    idx = _bdt_shape_bin_index(float(score), bdt_shape_edges)
    if idx is None:
        return

    bdt_shape_histos[mass][idx].Fill(h_mass, weight)


def _draw_bkg_mass_shapes_by_bdt(
    mass: int,
    histos: List[TH1F],
    score_edges: List[float],
    plot_cfg: Plot_Config,
    output_dir: Path,
    logy: bool = False,
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
        draw_hist.SetLineColor(color)
        draw_hist.SetLineWidth(3)
        draw_hist.SetMarkerSize(0)
        draw_hist.SetFillStyle(0)
        ymax = max(ymax, draw_hist.GetMaximum())
        draw_histos.append((idx, draw_hist))

    if not draw_histos:
        print(f"[mA={mass:02d}] No background entries for BDT-binned mass-shape plot.")
        return

    y_axis_max = 1.5 * max(ymax, 1e-6)
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
    frame.GetXaxis().SetTitleOffset(1.34)
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
) -> None:
    mva_branch_map: Dict[str, Dict[int, Optional[str]]] = {}
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

        missing = [mass for mass, branch in branch_map.items() if branch is None]
        if missing:
            print(f"[Warning] sample={sample} missing MVA branch for mA={missing}")

    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample]
        entries = ntup.GetEntries()
        print(f"[Fill] sample={sample:>10s} entries={entries}")

        masses_to_fill = list(mva_branch_map.get(sample, {}).keys())
        if not masses_to_fill:
            continue

        for i_evt in range(entries):
            ntup.GetEvent(i_evt)
            # For short test to tune the style
            if i_evt == 1000:
                break
            if i_evt % 100000 == 1:
                print(f"  looking at event {i_evt}")

            h_mass = getattr(ntup, "H_m", -999.0)
            if h_mass < H_M_XMIN or h_mass > H_M_XMAX:
                continue

            z_mumu = abs(int(getattr(ntup, "z_mumu", 0)))
            z_ee = abs(int(getattr(ntup, "z_ee", 0)))
            if only_ele and z_mumu == 1:
                continue
            if only_mu and z_ee == 1:
                continue

            weight = float(getattr(ntup, "weight", 1.0))
            if sample == "Data":
                weight = 1.0

            for mass in masses_to_fill:
                cut = mva_cuts.get(mass)
                branch = mva_branch_map[sample].get(mass)
                if cut is None or branch is None:
                    continue

                score = getattr(ntup, branch, None)
                if score is None:
                    continue

                _fill_bdt_shape_histogram(
                    sample=sample,
                    mass=mass,
                    score=float(score),
                    h_mass=h_mass,
                    weight=weight,
                    analyzer_cfg=analyzer_cfg,
                    bdt_shape_histos=bdt_shape_histos,
                    bdt_shape_edges=bdt_shape_edges,
                )

                if score < cut:
                    continue

                if sample == "Data" and blind and BLIND_LOW <= h_mass <= BLIND_HIGH:
                    continue

                histos[mass][sample].Fill(h_mass, weight)
                for hist_set in extra_histos:
                    hist = hist_set.get(mass, {}).get(sample)
                    if hist is not None:
                        hist.Fill(h_mass, weight)


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Plot background mH sculpting after MVA cuts for mA = 1..30.")
    parser.add_argument("-y", "--Year", dest="year", default="run3", help="Only run3 is supported.")
    parser.add_argument("--mva-cut-json", default=DEFAULT_MVA_CUT_JSON, help="Path to MVA cut JSON.")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Directory for output plots.")
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
    output_dir_signal_mixture = output_dir / "signalShapeMixture"
    output_dir_bkg_only = output_dir / "bkgOnly"
    output_dir_bkg_only_bin5 = output_dir / "bkgOnly_bin5GeV"
    output_dir_signal_nearest_bin5 = output_dir / "signalNearest_bin5GeV"
    output_dir_bdt_mass_shapes = output_dir / "bdtMassShapes"
    output_dir_signal_mixture.mkdir(parents=True, exist_ok=True)
    output_dir_bkg_only.mkdir(parents=True, exist_ok=True)
    output_dir_bkg_only_bin5.mkdir(parents=True, exist_ok=True)
    output_dir_signal_nearest_bin5.mkdir(parents=True, exist_ok=True)
    if not args.skip_bdt_shape_plots:
        output_dir_bdt_mass_shapes.mkdir(parents=True, exist_ok=True)
    print(f"[Input] MVA cut JSON: {mva_cut_path}")
    print(f"[Output] Plot directory: {output_dir}")
    print(f"[Output] Signal mixture directory: {output_dir_signal_mixture}")
    print(f"[Output] Bkg-only directory: {output_dir_bkg_only}")
    print(f"[Output] Bkg-only 5 GeV directory: {output_dir_bkg_only_bin5}")
    print(f"[Output] Signal nearest 5 GeV directory: {output_dir_signal_nearest_bin5}")
    if not args.skip_bdt_shape_plots:
        print(f"[Output] BDT-binned bkg mass-shape directory: {output_dir_bdt_mass_shapes}")

    mva_cuts = _complete_mva_cuts(_parse_mva_cuts(str(mva_cut_path)), TARGET_MASSES)

    analyzer_cfg = Analyzer_Config("inclusive", args.year, 0, True)
    analyzer_cfg.mva = True
    analyzer_cfg.plot_output_path = str(output_dir)
    plot_cfg = Plot_Config(analyzer_cfg, args.year)
    ntuples = LoadNtuples(analyzer_cfg)
    if args.ele:
        channel_mode = "ele"
    elif args.mu:
        channel_mode = "mu"
    else:
        channel_mode = "inclusive"

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
    histos_signal_nearest_overlay_2gev = _book_histograms(
        analyzer_cfg,
        sample_names=analyzer_cfg.sig_names,
        bin_edges=_build_uniform_bin_edges(
            H_M_XMIN,
            H_M_BIN_SIGNAL_OVERLAY_XMAX,
            H_M_BIN_SIGNAL_OVERLAY_WIDTH,
        ),
    )
    if args.skip_bdt_shape_plots:
        bdt_shape_histos = None
        bdt_shape_edges = None
    else:
        bdt_shape_histos, bdt_shape_edges = _book_bdt_shape_histograms(n_bdt_bins=args.bdt_shape_bins)
    _fill_histograms(
        ntuples=ntuples,
        analyzer_cfg=analyzer_cfg,
        mva_cuts=mva_cuts,
        histos=histos,
        extra_histos=[histos_bkg_only_bin5, histos_signal_nearest_bin5, histos_signal_nearest_overlay_2gev],
        bdt_shape_histos=bdt_shape_histos,
        bdt_shape_edges=bdt_shape_edges,
        blind=args.blind,
        only_ele=args.ele,
        only_mu=args.mu,
    )

    for mass in TARGET_MASSES:
        print(
            f"[mA={mass:02d}] yields: "
            f"Data={histos[mass]['Data'].Integral():.3f}, "
            f"DYJetsToLL={histos[mass]['DYJetsToLL'].Integral():.3f}, "
            f"DYGto2LG={histos[mass]['DYGto2LG'].Integral():.3f}"
        )
        signal_sample = _signal_sample_name(mass)
        if signal_sample:
            print(f"[mA={mass:02d}] signal {signal_sample}={histos[mass][signal_sample].Integral():.3f}")
        _draw_mass_plot(
            mass=mass,
            histos=histos[mass],
            all_histos=histos,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dir_signal_mixture,
            logy=args.ln,
            show_signal=True,
            name_suffix="",
            channel_mode=channel_mode,
            signal_shape_mode="mixture",
            signal_scale=SIGNAL_DRAW_SCALE,
        )
        _draw_mass_plot(
            mass=mass,
            histos=histos[mass],
            all_histos=histos,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dir_bkg_only,
            logy=args.ln,
            show_signal=False,
            name_suffix="_bkgOnly",
            channel_mode=channel_mode,
        )
        _draw_mass_plot(
            mass=mass,
            histos=histos_bkg_only_bin5[mass],
            all_histos=histos_bkg_only_bin5,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dir_bkg_only_bin5,
            logy=args.ln,
            show_signal=False,
            name_suffix="_bkgOnly_bin5GeV",
            channel_mode=channel_mode,
        )
        _draw_mass_plot(
            mass=mass,
            histos=histos_signal_nearest_bin5[mass],
            all_histos=histos_signal_nearest_bin5,
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dir_signal_nearest_bin5,
            logy=args.ln,
            show_signal=True,
            name_suffix="_signalNearest_bin5GeV",
            channel_mode=channel_mode,
            signal_shape_mode="nearest",
            signal_scale=SIGNAL_DRAW_SCALE_NEAREST_BIN5,
            signal_source_histos=histos_signal_nearest_overlay_2gev[mass],
            signal_source_all_histos=histos_signal_nearest_overlay_2gev,
        )
        if bdt_shape_histos is not None and bdt_shape_edges is not None:
            _draw_bkg_mass_shapes_by_bdt(
                mass=mass,
                histos=bdt_shape_histos[mass],
                score_edges=bdt_shape_edges,
                plot_cfg=plot_cfg,
                output_dir=output_dir_bdt_mass_shapes,
                logy=args.ln,
            )

    elapsed = time.time() - start_time
    print(f"Total runtime: {_format_elapsed_hms(elapsed)}")
    print("Done")


if __name__ == "__main__":
    main()
