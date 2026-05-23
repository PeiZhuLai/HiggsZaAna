####################################################
####################################################

import os
import sys
import numpy as np
import gc
import time  # NEW
import errno
import shutil

sys.path.insert(0, '%s/lib' % os.getcwd())
PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
HZA_MVA_SCRIPTS_DIR = os.path.join(PROJECT_DIR, 'HZaMVA', 'scripts')
if HZA_MVA_SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, HZA_MVA_SCRIPTS_DIR)

from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend, Total_Unc, ScaleBkgToData
from Analyzer_Helper import getMassSigma
import Analyzer_Configs as AC
import Plot_Configs     as PC

from Analyzer_ALP import PIso2D

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle
import copy 
import random
from sideband_reweight import load_sideband_reweighter

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="run3", help="which year's datasetes")
parser.add_argument("--region", dest="region", type=int, default=0, help="0 for full region, 1 for signal region, 2 for sideband region")
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument('--cut', dest='cut', action='store_true', default=False, help='apply mva')
parser.add_argument("--cutVal", dest="cutVal", type=float, default=0.0, help="mva cut value")
parser.add_argument("--mA", dest="mA", default="M5", help="ALP mass")
parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')
parser.add_argument('-ln', '--ln', dest='ln', action='store_true', default=False, help='log plot?')
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
parser.add_argument('--useSidebandReweight', dest='use_sideband_reweight', action='store_true', default=False, help='use sideband-reweighted background weights')
parser.add_argument('--sidebandReweightJson', dest='sideband_reweight_json', default=None, help='sideband reweight JSON path; defaults to HZA_SIDEBAND_REWEIGHT_JSON or HZaMVA/reweights/sideband_run3_iterative.json')
parser.add_argument('--sidebandReweightUnc', dest='sideband_reweight_unc', action='store_true', default=False, help='add sideband reweight uncertainty to the MC error band')
parser.add_argument('--noSidebandReweightUnc', dest='sideband_reweight_unc', action='store_false', help=argparse.SUPPRESS)
parser.add_argument('--outputTag', dest='output_tag', default=None, help='append a tag to the output ROOT file and plot directory')
parser.add_argument('--samples', dest='samples', default=None, help="comma-separated sample names to run; aliases: all, data, bkg, sig; source filters are supported as sample@year or DYGto2LG@subsample:year")
parser.add_argument('--histOnly', dest='hist_only', action='store_true', default=False, help='write raw_plots/sys_dir only; skip stack and PDF drawing')
parser.add_argument('--skipSystematics', dest='skip_systematics', action='store_true', default=False, help='write nominal histograms and copy them into sys_dir instead of filling per-event systematic variations')
parser.add_argument('--optimizeBranches', dest='optimize_branches', action='store_true', default=False, help='disable unused TTree branches before the event loop')
parser.add_argument('--maxEvents', dest='max_events', type=int, default=-1, help='maximum events per sample; <=0 runs all events')
# 新增：MVA 偵錯輸出控制
parser.add_argument('--mvaDebug', dest='mva_debug', action='store_true', default=False, help='print per-mA fill info')
parser.add_argument('--mvaDebugN', dest='mva_debug_n', type=int, default=15, help='max debug prints per (sample,mA)')
args = parser.parse_args()

def _append_unique(items, value):
    if value not in items:
        items.append(value)

def _parse_source_selectors(sample, source_spec, analyzer_cfg):
    selectors = []
    raw_tokens = str(source_spec).replace("+", ",").split(",")

    if sample in analyzer_cfg.sig_names:
        allowed_years = set(getattr(analyzer_cfg, "years_sig", []))
        for token in raw_tokens:
            year = token.strip()
            if not year:
                continue
            if year not in allowed_years:
                raise ValueError(
                    "Unknown source year '%s' for sample '%s'. Allowed years are: %s"
                    % (year, sample, ", ".join(sorted(allowed_years)))
                )
            _append_unique(selectors, year)
        return selectors

    if sample in ("Data", "DYJetsToLL"):
        allowed_years = set(getattr(analyzer_cfg, "years_dyll", []))
        for token in raw_tokens:
            year = token.strip()
            if not year:
                continue
            if year not in allowed_years:
                raise ValueError(
                    "Unknown source year '%s' for sample '%s'. Allowed years are: %s"
                    % (year, sample, ", ".join(sorted(allowed_years)))
                )
            _append_unique(selectors, year)
        return selectors

    if sample == "DYGto2LG":
        allowed_sources = set()
        for subsample in getattr(analyzer_cfg, "bkg_2022", []):
            for year in getattr(analyzer_cfg, "years_22", []):
                allowed_sources.add((subsample, year))
        for subsample in getattr(analyzer_cfg, "bkg_2023", []):
            for year in getattr(analyzer_cfg, "years_23", []):
                allowed_sources.add((subsample, year))

        for token in raw_tokens:
            selector = token.strip()
            if not selector:
                continue
            if ":" not in selector:
                raise ValueError(
                    "DYGto2LG source selector must be '<subsample>:<year>', got '%s'"
                    % selector
                )
            subsample, year = [part.strip() for part in selector.split(":", 1)]
            source = (subsample, year)
            if source not in allowed_sources:
                allowed_text = ", ".join(
                    "%s:%s" % (allowed_subsample, allowed_year)
                    for allowed_subsample, allowed_year in sorted(allowed_sources)
                )
                raise ValueError(
                    "Unknown source '%s:%s' for sample '%s'. Allowed sources are: %s"
                    % (subsample, year, sample, allowed_text)
                )
            _append_unique(selectors, source)
        return selectors

    raise ValueError("Sample source selectors are not supported for sample '%s'" % sample)

def _format_source_selectors(source_selectors):
    parts = []
    for item in source_selectors:
        if isinstance(item, tuple):
            parts.append("%s:%s" % item)
        else:
            parts.append(str(item))
    return ", ".join(parts)

def _parse_sample_filter(sample_filter, analyzer_cfg):
    if sample_filter is None:
        return None, {}

    allowed = list(analyzer_cfg.samp_names)
    selected = []
    source_filters = {}
    for token in str(sample_filter).replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue

        sample_token, has_source_filter, source_spec = token.partition("@")
        sample_token = sample_token.strip()
        source_spec = source_spec.strip() if has_source_filter else None
        token_lower = sample_token.lower()
        if token_lower in ("all", "*"):
            if source_spec:
                raise ValueError("Sample aliases like '%s' cannot use a source filter." % sample_token)
            return None, {}
        if token_lower in ("data",):
            if source_spec:
                raise ValueError("Sample aliases like '%s' cannot use a source filter." % sample_token)
            expanded = ["Data"]
        elif token_lower in ("bkg", "bkgs", "background", "backgrounds"):
            if source_spec:
                raise ValueError("Sample aliases like '%s' cannot use a source filter." % sample_token)
            expanded = list(analyzer_cfg.bkg_names)
        elif token_lower in ("sig", "sigs", "signal", "signals"):
            if source_spec:
                raise ValueError("Sample aliases like '%s' cannot use a source filter." % sample_token)
            expanded = list(analyzer_cfg.sig_names)
        else:
            expanded = [sample_token]

        for sample in expanded:
            if sample not in allowed:
                raise ValueError(
                    "Unknown sample '%s'. Allowed samples are: %s"
                    % (sample, ", ".join(allowed))
                )
            if sample not in selected:
                selected.append(sample)
            if source_spec:
                sample_sources = source_filters.setdefault(sample, [])
                for source in _parse_source_selectors(sample, source_spec, analyzer_cfg):
                    _append_unique(sample_sources, source)
            else:
                source_filters.pop(sample, None)

    return (selected if selected else None), source_filters

SIDEBAND_REWEIGHTER = None
SIDEBAND_REWEIGHT_UNC_SYS_NAMES = ("sideband_rwgt_reweight_up", "sideband_rwgt_reweight_down")
SIDEBAND_REWEIGHT_UNCERTAINTY_FRACTION = 1.0
if args.use_sideband_reweight:
    if args.sideband_reweight_json:
        SIDEBAND_REWEIGHTER = load_sideband_reweighter(args.sideband_reweight_json, required=True)
    else:
        SIDEBAND_REWEIGHTER = load_sideband_reweighter()
    if SIDEBAND_REWEIGHTER is not None:
        print("[SidebandReweight] Loaded JSON:", SIDEBAND_REWEIGHTER.source_path)
    else:
        print("[SidebandReweight] JSON not found; will use ROOT weight_sideband_rwgt branch when available.")
    if not args.sideband_reweight_unc:
        print("[SidebandReweight] Reweight uncertainty is disabled. Use --sidebandReweightUnc to enable it.")

def _lookup_sideband_factor(value, edges, factors):
    try:
        value = float(value)
    except Exception:
        return 1.0
    if not np.isfinite(value):
        return 1.0
    idx = np.searchsorted(np.asarray(edges, dtype=float), value, side="right") - 1
    idx = int(np.clip(idx, 0, len(factors) - 1))
    out = float(factors[idx])
    return out if np.isfinite(out) else 1.0

def _sideband_step_weight(reweighter, ntup, step, row_index=None):
    value = reweighter._object_value(ntup, step["var"], row_index=row_index)
    factor = _lookup_sideband_factor(value, step["edges"], step["clipped_factors"])
    norm = step.get("normalization_scale", 1.0)
    try:
        norm = float(norm)
    except Exception:
        norm = 1.0
    if not np.isfinite(norm):
        norm = 1.0
    return factor * norm

def _sideband_reweight_for_object(reweighter, ntup, row_index=None, varied_var=None, variation_scale=1.0):
    weight = 1.0
    for iteration in getattr(reweighter, "iterations", []):
        for step in iteration.get("steps", []):
            step_weight = _sideband_step_weight(reweighter, ntup, step, row_index=row_index)
            if varied_var is not None and step.get("var") == varied_var:
                step_weight = 1.0 + variation_scale * (step_weight - 1.0)
            weight *= step_weight
    return weight if np.isfinite(weight) else 1.0

def get_sideband_reweight_uncertainty_weights(ntup, sample, analyzer_cfg, central_weight, row_index=None):
    if (
        not args.use_sideband_reweight
        or not args.sideband_reweight_unc
        or SIDEBAND_REWEIGHTER is None
        or not is_background_sample(sample, analyzer_cfg)
    ):
        return {}

    nominal = _sideband_reweight_for_object(SIDEBAND_REWEIGHTER, ntup, row_index=row_index)
    if nominal == 0.0 or not np.isfinite(nominal):
        return {
            SIDEBAND_REWEIGHT_UNC_SYS_NAMES[0]: central_weight,
            SIDEBAND_REWEIGHT_UNC_SYS_NAMES[1]: central_weight,
        }

    ratios = [1.0]
    for var in getattr(SIDEBAND_REWEIGHTER, "reweight_vars", []):
        for scale in (
            1.0 + SIDEBAND_REWEIGHT_UNCERTAINTY_FRACTION,
            1.0 - SIDEBAND_REWEIGHT_UNCERTAINTY_FRACTION,
        ):
            varied = _sideband_reweight_for_object(
                SIDEBAND_REWEIGHTER,
                ntup,
                row_index=row_index,
                varied_var=var,
                variation_scale=scale,
            )
            ratio = varied / nominal
            if np.isfinite(ratio) and ratio > 0.0:
                ratios.append(ratio)

    return {
        SIDEBAND_REWEIGHT_UNC_SYS_NAMES[0]: central_weight * max(ratios),
        SIDEBAND_REWEIGHT_UNC_SYS_NAMES[1]: central_weight * min(ratios),
    }

def is_background_sample(sample, analyzer_cfg):
    return sample in analyzer_cfg.bkg_names

def get_event_weight(ntup, sample, analyzer_cfg, row_index=None):
    weight = ntup.weight
    if not args.use_sideband_reweight or not is_background_sample(sample, analyzer_cfg):
        return weight

    try:
        return ntup.weight_sideband_rwgt
    except Exception:
        pass

    if SIDEBAND_REWEIGHTER is not None:
        return weight * SIDEBAND_REWEIGHTER.weight_for_object(ntup, row_index=row_index)

    return weight

SYSTEMATIC_WEIGHT_WARNINGS = set()

def _systematic_central_branch(sys_name):
    if sys_name.endswith("_up"):
        return sys_name[:-3] + "_central"
    if sys_name.endswith("_down"):
        return sys_name[:-5] + "_central"
    return None

def _collect_systematic_branch_names(sys_names):
    names = []
    for sys_name in sys_names:
        if sys_name in SIDEBAND_REWEIGHT_UNC_SYS_NAMES:
            continue
        _append_unique(names, sys_name)
        central_branch = _systematic_central_branch(sys_name)
        if central_branch:
            _append_unique(names, central_branch)
    return names

def _warn_systematic_fallback(sample, sys_name, message):
    key = (sample, sys_name, message)
    if key in SYSTEMATIC_WEIGHT_WARNINGS:
        return
    SYSTEMATIC_WEIGHT_WARNINGS.add(key)
    print("[Systematics][WARN] sample=%s sys=%s: %s" % (sample, sys_name, message))

def _safe_systematic_weight(ntup, sample, sys_name, nominal_weight):
    central_branch = _systematic_central_branch(sys_name)
    if central_branch is None:
        return nominal_weight

    try:
        varied = float(getattr(ntup, sys_name))
    except Exception:
        _warn_systematic_fallback(sample, sys_name, "missing varied branch; use nominal weight")
        return nominal_weight

    try:
        central = float(getattr(ntup, central_branch))
    except Exception:
        _warn_systematic_fallback(sample, sys_name, "missing central branch; use nominal weight")
        return nominal_weight

    if not np.isfinite(varied):
        _warn_systematic_fallback(sample, sys_name, "varied branch is not finite; use nominal weight")
        return nominal_weight

    if not np.isfinite(central) or abs(central) < 1e-12:
        _warn_systematic_fallback(sample, sys_name, "central branch is zero/non-finite; use nominal weight")
        return nominal_weight

    ratio = varied / central
    if not np.isfinite(ratio):
        _warn_systematic_fallback(sample, sys_name, "varied/central ratio is not finite; use nominal weight")
        return nominal_weight

    return nominal_weight * ratio


gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


# load the model from disk
# model = pickle.load(open(BDT_filename, 'rb'))

pdfName_map = {'pho1Pt': '1_pho1Pt', 'pho1R9': '2_pho1R9', 'pho1IetaIeta55': '3_pho1IetaIeta55', 
'pho2Pt': '4_pho2Pt', 'pho2R9': '5_pho2R9', 'pho2IetaIeta55': '6_pho2IetaIeta55',
'pho1ECALIso': '7_pho1ECALIso', 'pho2ECALIso': '8_pho2ECALIso', 'ALP_calculatedPhotonIso': '9_ALP_calculatedPhotonIso',
'var_dR_Za': '10_var_dR_Za', 'var_dR_g1g2': '11_var_dR_g1g2', 'var_dR_g1Z': '12_var_dR_g1Z',
'var_PtaOverMh': '13_var_PtaOverMh', 'H_pt': '14_H_pt', 'param': '15_param', 'ALP_m': '16_ALP_m'}

# NEW: 輸出檔名依 pdfName_map 排序；沒對應就用原本 var_name
def _pdf_output_name(var_name: str) -> str:
    return pdfName_map.get(var_name, var_name)

MVA_LARGER_NBINS = 20
MVA_LARGER_XMIN = 0.0
MVA_LARGER_XMAX = 1.0
MVA_FULL_RANGE_BLIND_BINS = 2
MVA_FULL_RANGE_DATA_BLIND_THRESHOLD = (
    MVA_LARGER_XMAX
    - MVA_FULL_RANGE_BLIND_BINS * (MVA_LARGER_XMAX - MVA_LARGER_XMIN) / MVA_LARGER_NBINS
)
SIGMA_REGION_SCALES = {
    "1sigma": 1.0,
    "1P5sigma": 1.5,
    "2sigma": 2.0,
    "3sigma": 3.0,
}

def _is_full_range_mva_larger_var(var_name, target_masses):
    return var_name in [f"mvaVal_larger_{mass}" for mass in target_masses]

def _sigma_region_info(var_name, target_masses):
    for mass in target_masses:
        for region, scale in SIGMA_REGION_SCALES.items():
            if var_name in (
                f"mvaVal_{region}_{mass}",
                f"mvaVal_larger_{region}_{mass}",
            ):
                return mass, scale
    return None

def _passes_sigma_region(h_mass, sigma_low, sigma_hig, mass, scale):
    return (
        h_mass < 125.0 + sigma_hig[mass] * scale
        and h_mass > 125.0 + sigma_low[mass] * scale
    )

def _should_fill_var_for_event(var_name, h_mass, sigma_low, sigma_hig, target_masses):
    sigma_info = _sigma_region_info(var_name, target_masses)
    if sigma_info is None:
        return True
    mass, scale = sigma_info
    return _passes_sigma_region(h_mass, sigma_low, sigma_hig, mass, scale)

def _visible_integral(hist):
    if not hist:
        return 0.0
    return float(hist.Integral(1, hist.GetNbinsX()))

def _cdf_edges_from_reference_hist(reference_hist):
    if not reference_hist or _visible_integral(reference_hist) <= 0.0:
        return None

    nbins = reference_hist.GetNbinsX()
    bin_weights = [
        max(0.0, float(reference_hist.GetBinContent(i_bin)))
        for i_bin in range(1, nbins + 1)
    ]
    total = sum(bin_weights)
    if total <= 0.0:
        return None

    cdf_edges = [0.0]
    running = 0.0
    for weight in bin_weights:
        running += weight
        cdf_edges.append(min(1.0, max(0.0, running / total)))
    cdf_edges[-1] = 1.0
    return cdf_edges

def _make_cdf_transformed_hist(hist, name, cdf_edges):
    transformed_hist = hist.Clone(name)
    transformed_hist.Reset("ICES")
    transformed_hist.SetMaximum(-1111)
    transformed_hist.SetMinimum(-1111)
    transformed_hist.SetDirectory(0)

    if not cdf_edges:
        return transformed_hist

    out_axis = transformed_hist.GetXaxis()
    out_nbins = transformed_hist.GetNbinsX()
    for i_bin in range(1, hist.GetNbinsX() + 1):
        content = float(hist.GetBinContent(i_bin))
        error = float(hist.GetBinError(i_bin))
        if content == 0.0 and error == 0.0:
            continue

        u_low = min(1.0, max(0.0, cdf_edges[i_bin - 1]))
        u_high = min(1.0, max(0.0, cdf_edges[i_bin]))
        if u_high <= u_low:
            out_bin = transformed_hist.FindFixBin(u_high)
            out_bin = min(out_nbins, max(1, out_bin))
            transformed_hist.SetBinContent(out_bin, transformed_hist.GetBinContent(out_bin) + content)
            transformed_hist.SetBinError(out_bin, np.sqrt(transformed_hist.GetBinError(out_bin) ** 2 + error ** 2))
            continue

        width = u_high - u_low
        for out_bin in range(1, out_nbins + 1):
            bin_low = out_axis.GetBinLowEdge(out_bin)
            bin_high = out_axis.GetBinUpEdge(out_bin)
            overlap = max(0.0, min(u_high, bin_high) - max(u_low, bin_low))
            if overlap <= 0.0:
                continue
            frac = overlap / width
            transformed_hist.SetBinContent(out_bin, transformed_hist.GetBinContent(out_bin) + content * frac)
            transformed_hist.SetBinError(out_bin, np.sqrt(transformed_hist.GetBinError(out_bin) ** 2 + (error * frac) ** 2))

    transformed_hist.SetEntries(hist.GetEntries())
    return transformed_hist

def _build_cdf_histo_maps(var_name, histos_var, histos_sys_var, analyzer_cfg):
    """
    Build normal histograms after a CDF-based score transform.
    The reference CDF is taken from the matching signal mass when available,
    matching the usual BDT score_t convention.
    """
    mass_tag = var_name.split("_")[-1]
    reference_hist = histos_var.get(mass_tag) or histos_var.get("Data")
    cdf_edges = _cdf_edges_from_reference_hist(reference_hist)

    histos_cdf = {}
    histos_sys_cdf = {}

    for sample in analyzer_cfg.samp_names:
        histos_cdf[sample] = _make_cdf_transformed_hist(histos_var[sample], f"{var_name}_{sample}_cdf", cdf_edges)
        histos_sys_cdf[sample] = {}

        for sys_name in analyzer_cfg.sys_names:
            histos_sys_cdf[sample][sys_name] = _make_cdf_transformed_hist(
                histos_sys_var[sample][sys_name],
                f"{var_name}_{sample}_{sys_name}_cdf",
                cdf_edges,
            )

    return histos_cdf, histos_sys_cdf

def SideBandScaleBkgToData(histos, histos_sys, analyzer_cfg, signal_low=115., signal_high=135.):
    """
    只針對 histos['H_m']：
      將背景總和在 sideband (H_m < signal_low 或 H_m > signal_high) 的 yield
      正規化到資料在同一 sideband 的 yield。
    """
    if 'H_m' not in histos:
        print('[SideBandScale] H_m not found, skip.')
        return 1.0

    # 資料 sample key 可能大小寫不一致
    data_key = None
    for k in histos['H_m'].keys():
        if k.lower() == 'data':
            data_key = k
            break
    if data_key is None:
        print('[SideBandScale] Data hist not found, skip.')
        return 1.0

    axis = histos['H_m'][data_key].GetXaxis()
    nb = axis.GetNbins()
    # 找到 signal window 所對應的 bin
    sig_low_bin = axis.FindFixBin(signal_low)
    sig_high_bin = axis.FindFixBin(signal_high)

    def sideband_integral(h):
        # 左側
        left = 0.0
        if sig_low_bin > 1:
            left = h.Integral(1, sig_low_bin - 1)
        # 右側
        right = 0.0
        if sig_high_bin < nb:
            right = h.Integral(sig_high_bin + 1, nb)
        return left + right

    data_sb = sideband_integral(histos['H_m'][data_key])

    bkg_samples = [
        s for s in analyzer_cfg.samp_names
        if (s.lower() != 'data') and (s not in analyzer_cfg.sig_names)
    ]

    bkg_sb = 0.0
    for s in bkg_samples:
        bkg_sb += sideband_integral(histos['H_m'][s])

    if bkg_sb <= 0:
        print('[SideBandScale] Background sideband yield = 0, skip scaling.')
        return 1.0

    scale_factor = data_sb / bkg_sb
    # 套用縮放到背景 (包含 systematics)
    for s in bkg_samples:
        histos['H_m'][s].Scale(scale_factor)
        if histos_sys and 'H_m' in histos_sys:
            for sys_name in analyzer_cfg.sys_names:
                histos_sys['H_m'][s][sys_name].Scale(scale_factor)

    print(f'[SideBandScale] H_m sideband Data={data_sb:.3f}  Bkg(before)={bkg_sb:.3f}  Scale={scale_factor:.4f}')
    return scale_factor


# 新增：偵測指定 mA 的 MVA 分支名稱（回傳第一個存在的分支，否則預設為 "MVA_Score"）
def _resolve_mva_branch_for_mass(chain, mass_tag):
    """
    run3 更新：同一個 era 檔內同時有多個 mA 的 MVA 分數：
      - 優先：MVA_Score_mA_<M>  (e.g. MVA_Score_mA_M4)
    其餘命名仍保留回退相容。
    """
    if not chain or not chain.GetListOfBranches():
        return "MVA_Score"
    br_list = chain.GetListOfBranches()
    candidates = [
        f"MVA_Score_mA_{mass_tag}",   # NEW primary
        f"MVA_Score_{mass_tag}",
        f"MVA_score_{mass_tag}",
        f"BDTScore_{mass_tag}",
        f"BDT_{mass_tag}",
        "MVA_Score",
    ]
    for name in candidates:
        if br_list.FindObject(name):
            return name
    return "MVA_Score"

def _validate_and_publish_output(out_file, tmp_output_path, final_output_path):
    out_file.Close()

    check_file = TFile.Open(tmp_output_path, "READ")
    if not check_file or check_file.IsZombie():
        raise RuntimeError("[Output][ERROR] temporary ROOT file is unreadable: %s" % tmp_output_path)

    keys = check_file.GetListOfKeys()
    n_keys = keys.GetEntries() if keys else 0
    check_file.Close()
    if n_keys <= 0:
        raise RuntimeError("[Output][ERROR] temporary ROOT file has no keys: %s" % tmp_output_path)

    try:
        os.replace(tmp_output_path, final_output_path)
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise
        shutil.move(tmp_output_path, final_output_path)
    print("[Output] Published ROOT file: %s (%d top-level keys)" % (final_output_path, n_keys))

def _systematic_branch_names(sys_names):
    return _collect_systematic_branch_names(sys_names)

def _enable_used_branches(chain, sample, analyzer_cfg, mva_branches):
    if not chain:
        return

    needed = {
        "weight",
        "H_m",
        "Z_mass",
        "ALP_m",
        "H_pt",
        "pho1Pt",
        "ALP_lead_photon_pt",
        "ALP_lead_photon_eta",
        "ALP_lead_photon_phi",
        "ALP_lead_photon_r9",
        "ALP_lead_photon_sieie",
        "ALP_lead_photon_ecalPFClusterIso",
        "ALP_lead_photon_chiso",
        "ALP_lead_photon_hcalPFClusterIso",
        "ALP_lead_photon_hoe_PUcorr",
        "ALP_sublead_photon_pt",
        "ALP_sublead_photon_eta",
        "ALP_sublead_photon_phi",
        "ALP_sublead_photon_r9",
        "ALP_sublead_photon_sieie",
        "ALP_sublead_photon_ecalPFClusterIso",
        "ALP_sublead_photon_chiso",
        "ALP_sublead_photon_hcalPFClusterIso",
        "ALP_sublead_photon_hoe_PUcorr",
        "ALP_calculatedPhotonIso",
        "var_dR_Za",
        "var_dR_g1g2",
        "var_dR_g1Z",
        "var_PtaOverMh",
        "var_PtaOverMa",
        "var_Pta",
        "var_MhMa",
        "var_MhMZ",
    }

    if args.ele:
        needed.add("z_mumu")
    if args.mu:
        needed.add("z_ee")
    if args.mva:
        needed.update(branch for branch in mva_branches if branch)
    if args.use_sideband_reweight and is_background_sample(sample, analyzer_cfg):
        needed.update((
            "weight_sideband_rwgt",
            "event",
            "param",
            "H_mass",
            "ALP_mass",
            "pho1ECALIso",
            "pho1PIso_noCorr",
            "pho2ECALIso",
            "pho2PIso_noCorr",
        ))
    if not args.skip_systematics and is_background_sample(sample, analyzer_cfg):
        needed.update(_systematic_branch_names(analyzer_cfg.sys_names))

    branch_list = chain.GetListOfBranches()
    chain.SetBranchStatus("*", 0)
    enabled = 0
    for name in sorted(needed):
        if branch_list and branch_list.FindObject(name):
            chain.SetBranchStatus(name, 1)
            enabled += 1

    print("[Branches] sample=%s enabled %d/%d requested branches" % (sample, enabled, len(needed)))


def main():
    start_time = time.time()  # NEW

    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year, args.region, args.mva)

    analyzer_cfg.mva = bool(args.mva)
    analyzer_cfg.mva_alp_mass = str(args.mA) if args.mva else "M1"
    all_sig_names = list(analyzer_cfg.sig_names)
    if args.skip_systematics:
        print("[Systematics] --skipSystematics enabled: sys_dir will use nominal histogram clones.")
    if args.use_sideband_reweight and args.sideband_reweight_unc and SIDEBAND_REWEIGHTER is not None:
        for sys_name in SIDEBAND_REWEIGHT_UNC_SYS_NAMES:
            if sys_name not in analyzer_cfg.sys_names:
                analyzer_cfg.sys_names.append(sys_name)
        print("[SidebandReweight] Add uncertainty pair to MC error band:", SIDEBAND_REWEIGHT_UNC_SYS_NAMES)
    elif args.use_sideband_reweight:
        analyzer_cfg.sys_names = [
            sys_name
            for sys_name in analyzer_cfg.sys_names
            if sys_name not in SIDEBAND_REWEIGHT_UNC_SYS_NAMES
        ]
        print("[SidebandReweight] Apply sideband reweighting; sideband reweight uncertainty pair is not added.")
    if args.output_tag:
        output_tag = os.path.basename(str(args.output_tag).strip().strip('/'))
        if output_tag:
            root_base, root_ext = os.path.splitext(analyzer_cfg.root_output_name)
            analyzer_cfg.root_output_name = f"{root_base}_{output_tag}{root_ext or '.root'}"
            analyzer_cfg.plot_output_path = os.path.join(analyzer_cfg.plot_output_path, output_tag)

    # 在 mva + run3/run3_NFlow 下，自動按 self.sig_names 跑每個 mA；其他情況維持單一目標質量
    if args.mva and analyzer_cfg.year in ['run3', 'run3_NFlow']:
        target_masses = list(all_sig_names)
    else:
        target_masses = [analyzer_cfg.mva_alp_mass] if args.mva else []

    if not os.path.exists(analyzer_cfg.out_dir):
        os.makedirs(analyzer_cfg.out_dir)

    if not os.path.exists(analyzer_cfg.plot_output_path):
        os.makedirs(analyzer_cfg.plot_output_path)

    final_output_path = os.path.join(analyzer_cfg.out_dir, analyzer_cfg.root_output_name)
    tmp_output_dir = os.environ.get("ALP_OUTPUT_TMP_DIR", analyzer_cfg.out_dir)
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)
    tmp_output_path = os.path.join(
        tmp_output_dir,
        ".%s.%d.tmp" % (analyzer_cfg.root_output_name, os.getpid()),
    )
    if os.path.exists(tmp_output_path):
        os.remove(tmp_output_path)

    print("[Output] Writing temporary ROOT file: %s" % tmp_output_path)
    print("[Output] Final ROOT file: %s" % final_output_path)
    out_file = TFile(tmp_output_path, "RECREATE")
    if not out_file or out_file.IsZombie():
        raise RuntimeError("[Output][ERROR] cannot create temporary ROOT file: %s" % tmp_output_path)
    # model = pickle.load(open(analyzer_cfg.BDT_filename, 'rb'))

    
    if args.cut: 
        analyzer_cfg.sig_names = [args.mA]
        # 修正大小寫，與其他地方一致使用 'Data'
        analyzer_cfg.samp_names = analyzer_cfg.bkg_names + analyzer_cfg.sig_names + ['Data']

    sample_source_filters = {}
    if args.samples:
        try:
            selected_samples, sample_source_filters = _parse_sample_filter(args.samples, analyzer_cfg)
        except ValueError as exc:
            print("[Samples][ERROR]", exc)
            sys.exit(1)

        if selected_samples:
            if not args.hist_only:
                print("[Samples] --samples is set; enabling --histOnly because partial sample files cannot draw complete stacks.")
                args.hist_only = True
            selected_set = set(selected_samples)
            analyzer_cfg.samp_names = [s for s in analyzer_cfg.samp_names if s in selected_set]
            if not analyzer_cfg.samp_names:
                print("[Samples][ERROR] No samples selected after applying --samples=%s" % args.samples)
                sys.exit(1)
            print("[Samples] Running selected samples:", analyzer_cfg.samp_names)
            for sample in analyzer_cfg.samp_names:
                if sample in sample_source_filters:
                    print(
                        "[Samples] sample=%s source filters: %s"
                        % (sample, _format_source_selectors(sample_source_filters[sample]))
                    )

    analyzer_cfg.run3_source_filters = sample_source_filters

    analyzer_cfg.Print_Config()

    # NEW: 一律走 LoadNtuples（不再 split_by_mass；檔案本身就含多 mA 的 MVA branch）
    ntuples = LoadNtuples(analyzer_cfg)

    # NEW: 印出使用的樣本/背景清單，並檢查 ntuples 實際載入到哪些 sample
    print("\n[Samples] analyzer_cfg lists:")
    print(f"  bkg_names({len(analyzer_cfg.bkg_names)}): {analyzer_cfg.bkg_names}")
    print(f"  sig_names({len(analyzer_cfg.sig_names)}): {analyzer_cfg.sig_names}")
    print(f"  samp_names({len(analyzer_cfg.samp_names)}): {analyzer_cfg.samp_names}")
    print(f"  sys_names({len(analyzer_cfg.sys_names)}): {analyzer_cfg.sys_names}")

    loaded_keys = sorted(list(ntuples.keys())) if ntuples else []
    print("\n[Samples] ntuples loaded keys:")
    print(f"  loaded({len(loaded_keys)}): {loaded_keys}")

    cfg_set = set(analyzer_cfg.samp_names)
    loaded_set = set(loaded_keys)
    missing_in_ntuples = sorted(list(cfg_set - loaded_set))
    extra_in_ntuples = sorted(list(loaded_set - cfg_set))
    if missing_in_ntuples:
        print("\n[Warning][Samples] In cfg.samp_names but NOT loaded by LoadNtuples:")
        print(f"  missing({len(missing_in_ntuples)}): {missing_in_ntuples}")
    if extra_in_ntuples:
        print("\n[Info][Samples] Loaded by LoadNtuples but NOT in cfg.samp_names:")
        print(f"  extra({len(extra_in_ntuples)}): {extra_in_ntuples}")

    # 另外特別看 bkg 是否有缺
    bkg_expected = set([s for s in analyzer_cfg.bkg_names if s.lower() != "data"])
    bkg_loaded = set([s for s in loaded_keys if (s.lower() != "data") and (s in analyzer_cfg.bkg_names)])
    bkg_missing = sorted(list(bkg_expected - bkg_loaded))
    print("\n[Samples] bkg cross-check:")
    print(f"  expected_bkg({len(bkg_expected)}): {sorted(list(bkg_expected))}")
    print(f"  loaded_bkg({len(bkg_loaded)}): {sorted(list(bkg_loaded))}")
    if bkg_missing:
        print(f"  [Warning] missing_bkg({len(bkg_missing)}): {bkg_missing}")

    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)
    var_names = list(plot_cfg.var_title_map.keys())

    if args.mva:
        # 僅針對選定質量新增 MVA 變數
        var_mva = []
        for ALP_mass in target_masses:
            var_names.append('mvaVal_'+ALP_mass)
            var_names.append('mvaVal_larger_'+ALP_mass)
            for r in ['1sigma', '1P5sigma', '2sigma', '3sigma']:
                var_mva.append('mvaVal_'+r+'_'+ALP_mass)
                var_mva.append('mvaVal_larger_'+r+'_'+ALP_mass)
        var_names = var_names + var_mva

    # get mass region
    sigma_low, sigma_hig = getMassSigma(analyzer_cfg)

    # 為每個 (sample, mA) 解析一次 MVA 分支（從同一條 chain 上找 MVA_Score_mA_<mass>）
    mva_branch_map = {}  # dict[sample][mass] -> branch_name
    if args.mva:
        for s in analyzer_cfg.samp_names:
            chain = ntuples.get(s, None)
            mva_branch_map[s] = {}
            for ALP_mass in target_masses:
                br = _resolve_mva_branch_for_mass(chain, ALP_mass)
                mva_branch_map[s][ALP_mass] = br
                print(f"[MVA] sample={s:>12s} mass={ALP_mass:>3s} uses branch '{br}'")

    if args.optimize_branches:
        for s in analyzer_cfg.samp_names:
            _enable_used_branches(ntuples.get(s, None), s, analyzer_cfg, mva_branch_map.get(s, {}).values())

    ### declare histograms

    histos = {}
    histos_sys = {}

    for var_name in var_names:
        histos[var_name] = {}

    for sample in analyzer_cfg.samp_names:
        histos['pho1Pt'][sample]    = TH1F('pho1Pt'    + '_' + sample, 'pho1Pt'    + '_' + sample, 42,  8., 50.)
        histos['pho1eta'][sample]    = TH1F('pho1eta'    + '_' + sample, 'pho1eta'    + '_' + sample, 30,  -3., 3.)
        histos['pho1phi'][sample]    = TH1F('pho1phi'    + '_' + sample, 'pho1phi'    + '_' + sample, 40,  -4., 4.)
        histos['pho1R9'][sample]    = TH1F('pho1R9'    + '_' + sample, 'pho1R9'    + '_' + sample, 25,  0.1, 1.)
        histos['pho1IetaIeta55'][sample]    = TH1F('pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 25,  0., 0.08)
        histos['pho1ECALIso'][sample]    = TH1F('pho1ECALIso'    + '_' + sample, 'pho1ECALIso'    + '_' + sample, 50, 0., 50.)
        histos['pho1CIso'][sample]    = TH1F('pho1CIso'    + '_' + sample, 'pho1CIso'    + '_' + sample, 35, 0., 0.35)
        histos['pho1HCALIso'][sample]  = TH1F('pho1HCALIso'    + '_' + sample, 'pho1HCALIso'    + '_' + sample, 25, 0., 4.0)
        histos['pho1HOE'][sample]    = TH1F('pho1HOE'    + '_' + sample, 'pho1HOE'    + '_' + sample, 25, 0., 0.01)
        histos['pho2Pt'][sample]    = TH1F('pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 22,  8., 30.)
        histos['pho2eta'][sample]    = TH1F('pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 30,  -3., 3.)
        histos['pho2phi'][sample]    = TH1F('pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 40,  -4., 4.)
        histos['pho2R9'][sample]    = TH1F('pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 25,  0.1, 1.)
        histos['pho2IetaIeta55'][sample]    = TH1F('pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 25,  0., 0.08)
        histos['pho2ECALIso'][sample]    = TH1F('pho2ECALIso'    + '_' + sample, 'pho2ECALIso'    + '_' + sample, 50, 0., 50.)
        histos['pho2CIso'][sample]    = TH1F('pho2CIso'    + '_' + sample, 'pho2CIso'    + '_' + sample, 35, 0., 0.35)
        histos['pho2HCALIso'][sample]    = TH1F('pho2HCALIso'    + '_' + sample, 'pho2HCALIso'    + '_' + sample, 25, 0., 4.0)
        histos['pho2HOE'][sample]    = TH1F('pho2HOE'    + '_' + sample, 'pho2HOE'    + '_' + sample, 25, 0., 0.1)
        histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 80,  50., 130.)
        histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 85,  95., 180.)
        histos['H_pt'][sample]    = TH1F('H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 80,  0., 160.)
        histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 40, 0., 40.)
        histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 25, 0., 5)
        histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 25, 0., 100.)
        histos['var_dR_Za'][sample] = TH1F('var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 35, 0., 7.)
        histos['var_dR_g1Z'][sample] = TH1F('var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 35, 0., 7)
        histos['var_PtaOverMh'][sample] = TH1F('var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 25, 0., 0.75)
        histos['var_Pta'][sample] = TH1F('var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 60, 0., 60.)
        histos['var_MhMa'][sample] = TH1F('var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 100, 100., 200.)
        histos['var_MhMZ'][sample] = TH1F('var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 130, 180., 310.)
        histos['ALP_calculatedPhotonIso'][sample] = TH1F('ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 25, 0., 125.)
        histos['param'][sample] = TH1F('param' + '_' + sample, 'param' + '_' + sample, 25, -0.3, 0.6)

        if args.mva:
            # 僅為選定質量建立 MVA 相關直方圖
            for ALP_mass in target_masses:
                histos['mvaVal_'+ALP_mass][sample]           = TH1F('mvaVal_'+ALP_mass           + '_' + sample, 'mvaVal_'+ALP_mass    + '_' + sample, 240, -0.1, 1.1)
                histos['mvaVal_1sigma_'+ALP_mass][sample]    = TH1F('mvaVal_1sigma_'+ALP_mass    + '_' + sample, 'mvaVal_1sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_1P5sigma_'+ALP_mass][sample]  = TH1F('mvaVal_1P5sigma_'+ALP_mass  + '_' + sample, 'mvaVal_1P5sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_2sigma_'+ALP_mass][sample]    = TH1F('mvaVal_2sigma_'+ALP_mass    + '_' + sample, 'mvaVal_2sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_3sigma_'+ALP_mass][sample]    = TH1F('mvaVal_3sigma_'+ALP_mass    + '_' + sample, 'mvaVal_3sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)

                histos['mvaVal_larger_'+ALP_mass][sample]           = TH1F('mvaVal_larger_'+ALP_mass           + '_' + sample, 'mvaVal_larger_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)
                histos['mvaVal_larger_1sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_1sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_1sigma_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)
                histos['mvaVal_larger_1P5sigma_'+ALP_mass][sample]  = TH1F('mvaVal_larger_1P5sigma_'+ALP_mass  + '_' + sample, 'mvaVal_larger_1P5sigma_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)
                histos['mvaVal_larger_2sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_2sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_2sigma_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)
                histos['mvaVal_larger_3sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_3sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_3sigma_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)


    for var_name in var_names:
        histos_sys[var_name] = {}
        for sample in analyzer_cfg.samp_names:
            histos_sys[var_name][sample] = {}
            for sys in analyzer_cfg.sys_names:
                if not args.skip_systematics:
                    histos_sys[var_name][sample][sys] = copy.deepcopy(histos[var_name][sample])
                    histos_sys[var_name][sample][sys].SetNameTitle(var_name+'_'+sample+'_'+sys, var_name+'_'+sample+'_'+sys)

    ### loop over samples and events
    mass_list = {'M1':1.0, 'M2':2.0, 'M3':3.0, 'M4':4.0, 'M5':5.0, 'M6':6.0, 'M7':7.0, 'M8':8.0, 'M9':9.0, 'M10':10.0, 'M15':15.0, 'M20':20.0, 'M25':25.0, 'M30':30.0}
    mass_values = list(mass_list.values())
    # 新增：控制每個 (sample, mA) 的偵錯輸出次數
    debug_printed = {}

    # NEW: 移除/停用 split-by-mass 的事件迴圈，統一用原本單一鏈流程
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        fill_systematics = (not args.skip_systematics) and is_background_sample(sample, analyzer_cfg)
        print('\n\nOn sample: %s' %sample)
        print('total events: %d' %ntup.GetEntries())
        if not args.skip_systematics and not fill_systematics:
            print("[Systematics] sample=%s is not a background sample; sys_dir will use nominal clones." % sample)

        for iEvt in range( ntup.GetEntries() ):
    
            ntup.GetEvent(iEvt)
            if args.max_events > 0 and iEvt >= args.max_events:
                break

            if (iEvt % 100000 == 1):
                print("looking at event %d" %iEvt)

            
            if args.ele:
                if abs(ntup.z_mumu) == 1: 
                    continue
            if args.mu:
                if abs(ntup.z_ee) == 1: 
                    continue
            
            if (ntup.H_m > -90):
                if ntup.H_m>180. or ntup.H_m<95.: continue

                if  args.region == 1 and (ntup.H_m>135. or ntup.H_m<115.): continue
                if  args.region == 2 and (ntup.H_m<135. and ntup.H_m>115.): continue

                # weight = ntup.factor * ntup.pho1SFs * ntup.pho2SFs
                weight = get_event_weight(ntup, sample, analyzer_cfg, row_index=iEvt)
                sideband_rwgt_sys_weights = {}

                MVA_value = {}
                if args.mva:
                    # 以各 mA 對應分支取得 MVA 分數（若分支不存在，回退 "MVA_Score"，再不行設為 -1）
                    for ALP_mass in target_masses:
                        br = mva_branch_map.get(sample, {}).get(ALP_mass, "MVA_Score")
                        try:
                            MVA_value[ALP_mass] = getattr(ntup, br)
                        except Exception:
                            # 兩層保險：先試指定分支，再試通用分支，最後給 -1.0
                            try:
                                MVA_value[ALP_mass] = getattr(ntup, "MVA_Score")
                            except Exception:
                                MVA_value[ALP_mass] = -1.0

                        # 新增：每個 (sample, mA) 前 args.mva_debug_n 筆偵錯印出
                        if args.mva_debug:
                            key = (sample, ALP_mass)
                            if debug_printed.get(key, 0) < args.mva_debug_n:
                                print(f"[MVA-FILL] sample={sample:>10s} mA={ALP_mass:>3s} branch={br:>20s} score={MVA_value[ALP_mass]:+.4f} H_m={ntup.H_m:6.2f} w={weight: .3e}")
                                debug_printed[key] = debug_printed.get(key, 0) + 1

                    if args.cut:
                        if MVA_value[args.mA] < analyzer_cfg.mvaCut[args.mA]: continue

                    for ALP_mass in target_masses:
                        if args.blind:
                            if not (sample == 'Data' and MVA_value[ALP_mass] >= MVA_FULL_RANGE_DATA_BLIND_THRESHOLD):
                                histos['mvaVal_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                                histos['mvaVal_larger_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                        else:
                            histos['mvaVal_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )

                        if _passes_sigma_region(ntup.H_m, sigma_low, sigma_hig, ALP_mass, 1.0):
                            histos['mvaVal_1sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_1sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )

                        if _passes_sigma_region(ntup.H_m, sigma_low, sigma_hig, ALP_mass, 1.5):
                            histos['mvaVal_1P5sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_1P5sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )

                        if _passes_sigma_region(ntup.H_m, sigma_low, sigma_hig, ALP_mass, 2.0):
                            histos['mvaVal_2sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_2sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )

                        if _passes_sigma_region(ntup.H_m, sigma_low, sigma_hig, ALP_mass, 3.0):
                            histos['mvaVal_3sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_3sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )

                histos['pho1Pt'][sample].Fill( ntup.ALP_lead_photon_pt, weight )
                histos['pho1eta'][sample].Fill( ntup.ALP_lead_photon_eta, weight )
                histos['pho1phi'][sample].Fill( ntup.ALP_lead_photon_phi, weight )
                histos['pho1R9'][sample].Fill( ntup.ALP_lead_photon_r9, weight )
                histos['pho1IetaIeta55'][sample].Fill( ntup.ALP_lead_photon_sieie, weight )
                histos['pho1ECALIso'][sample].Fill( ntup.ALP_lead_photon_ecalPFClusterIso, weight )
                histos['pho2Pt'][sample].Fill( ntup.ALP_sublead_photon_pt, weight )
                histos['pho2eta'][sample].Fill( ntup.ALP_sublead_photon_eta, weight )
                histos['pho2phi'][sample].Fill( ntup.ALP_sublead_photon_phi, weight )
                histos['pho2R9'][sample].Fill( ntup.ALP_sublead_photon_r9, weight )
                histos['pho2IetaIeta55'][sample].Fill( ntup.ALP_sublead_photon_sieie, weight )
                histos['pho2ECALIso'][sample].Fill( ntup.ALP_sublead_photon_ecalPFClusterIso, weight )

                histos['pho1CIso'][sample].Fill( ntup.ALP_lead_photon_chiso, weight)
                histos['pho1HCALIso'][sample].Fill( ntup.ALP_lead_photon_hcalPFClusterIso, weight)
                histos['pho1HOE'][sample].Fill( ntup.ALP_lead_photon_hoe_PUcorr, weight)
                histos['pho2CIso'][sample].Fill( ntup.ALP_sublead_photon_chiso, weight)
                histos['pho2HCALIso'][sample].Fill( ntup.ALP_sublead_photon_hcalPFClusterIso, weight)
                histos['pho2HOE'][sample].Fill( ntup.ALP_sublead_photon_hoe_PUcorr, weight)

                if args.blind:
                    if not (sample == 'Data' and (ntup.H_m<135. and ntup.H_m>115.)): 
                        histos['H_m'][sample].Fill( ntup.H_m, weight )
                else:        
                    histos['H_m'][sample].Fill( ntup.H_m, weight )

                histos['H_pt'][sample].Fill( ntup.H_pt, weight )
                histos['ALP_m'][sample].Fill( ntup.ALP_m, weight )
                histos['Z_m'][sample].Fill( ntup.Z_mass, weight )

                histos['var_dR_Za'][sample].Fill( ntup.var_dR_Za, weight )
                histos['var_dR_g1g2'][sample].Fill( ntup.var_dR_g1g2, weight )
                histos['var_dR_g1Z'][sample].Fill( ntup.var_dR_g1Z, weight )
                histos['var_PtaOverMa'][sample].Fill( ntup.var_PtaOverMa, weight )
                histos['var_PtaOverMh'][sample].Fill( ntup.var_PtaOverMh, weight )
                histos['var_Pta'][sample].Fill( ntup.var_Pta, weight )
                histos['var_MhMa'][sample].Fill( ntup.var_MhMa, weight )
                histos['var_MhMZ'][sample].Fill( ntup.var_MhMZ, weight )
                histos['ALP_calculatedPhotonIso'][sample].Fill( ntup.ALP_calculatedPhotonIso, weight )

                if sample in analyzer_cfg.sig_names:
                    param_val = (ntup.ALP_m - mass_list[sample])/ntup.H_m
                else:
                    # mass_random = random.choice(search_mA_list)
                    mass_random = random.choice(mass_values)
                    param_val = (ntup.ALP_m - mass_random)/ntup.H_m

                histos['param'][sample].Fill( param_val, weight )

                if args.sideband_reweight_unc and fill_systematics:
                    sideband_rwgt_sys_weights = get_sideband_reweight_uncertainty_weights(
                        ntup,
                        sample,
                        analyzer_cfg,
                        weight,
                        row_index=iEvt,
                    )

                if fill_systematics:
                    var_map = {'Z_m':ntup.Z_mass, 'H_m':ntup.H_m, 'ALP_m':ntup.ALP_m,'pho1Pt':ntup.pho1Pt, 'pho1eta':ntup.ALP_lead_photon_eta, 'pho1phi':ntup.ALP_lead_photon_phi, 'pho1R9':ntup.ALP_lead_photon_r9, 'pho1IetaIeta':ntup.ALP_lead_photon_sieie, 'pho1IetaIeta55':ntup.ALP_lead_photon_sieie,'pho1ECALIso':ntup.ALP_lead_photon_ecalPFClusterIso, 'pho1CIso':ntup.ALP_lead_photon_chiso, 'pho1HCALIso':ntup.ALP_lead_photon_hcalPFClusterIso, 'pho1HOE':ntup.ALP_lead_photon_hoe_PUcorr, 'pho2Pt':ntup.ALP_sublead_photon_pt, 'pho2eta':ntup.ALP_sublead_photon_eta, 'pho2phi':ntup.ALP_sublead_photon_phi, 'pho2R9':ntup.ALP_sublead_photon_r9, 'pho2IetaIeta':ntup.ALP_sublead_photon_sieie, 'pho2IetaIeta55':ntup.ALP_sublead_photon_sieie,'pho2ECALIso':ntup.ALP_sublead_photon_ecalPFClusterIso, 'pho2CIso':ntup.ALP_sublead_photon_chiso, 'pho2HCALIso':ntup.ALP_sublead_photon_hcalPFClusterIso, 'pho2HOE':ntup.ALP_sublead_photon_hoe_PUcorr,'ALP_calculatedPhotonIso':ntup.ALP_calculatedPhotonIso, 'var_dR_Za':ntup.var_dR_Za, 'var_dR_g1g2':ntup.var_dR_g1g2, 'var_dR_g1Z':ntup.var_dR_g1Z, 'var_PtaOverMh':ntup.var_PtaOverMh, 'var_Pta':ntup.var_Pta, 'var_MhMZ':ntup.var_MhMZ, 'H_pt':ntup.H_pt, 'var_PtaOverMa':ntup.var_PtaOverMa, 'var_MhMa':ntup.var_MhMa, 'param':param_val}
                    if args.mva:
                        for ALP_mass in target_masses:
                            var_map['mvaVal_'+ALP_mass] = MVA_value[ALP_mass]
                            var_map['mvaVal_larger_'+ALP_mass] = MVA_value[ALP_mass]
                            for r in ['1sigma', '1P5sigma', '2sigma', '3sigma']:
                                var_map['mvaVal_'+r+'_'+ALP_mass] = MVA_value[ALP_mass]
                                var_map['mvaVal_larger_'+r+'_'+ALP_mass] = MVA_value[ALP_mass]
                    for sys_name in analyzer_cfg.sys_names:
                        if sys_name in sideband_rwgt_sys_weights:
                            weight_sys = sideband_rwgt_sys_weights[sys_name]
                        else:
                            weight_sys = _safe_systematic_weight(ntup, sample, sys_name, weight)

                        for var in var_names:
                            if _should_fill_var_for_event(var, ntup.H_m, sigma_low, sigma_hig, target_masses):
                                histos_sys[var][sample][sys_name].Fill(var_map[var], weight_sys)
                



        ## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names

    # 新增：MVA 直方圖摘要（每個 mA x 每個樣本）
    if args.mva:
        print("\n[MVA] Per-mass histogram summary:")
        for ALP_mass in target_masses:
            hname = 'mvaVal_' + ALP_mass
            if hname not in histos:
                continue
            for s in analyzer_cfg.samp_names:
                h = histos[hname].get(s)
                if not h:
                    continue
                print(f"  mass={ALP_mass:>3s} sample={s:>10s} Entries={h.GetEntries():8.0f}  Integral={h.Integral():.3f}")

    for var in ['H_m', 'pho1Pt', 'Z_m']:
        entries = histos[var][sample].GetEntries()
        print(f"[{sample}] {var} histogram entries: {entries}")


    ### save raw histograms
    raw_dir = out_file.mkdir('raw_plots')
    raw_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            plot_cfg.SetHistStyles(histos[var_name][sample], sample)
            histos[var_name][sample].Write()

    #### save sys histograms
    sys_dir = out_file.mkdir('sys_dir')
    sys_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            for sys in analyzer_cfg.sys_names:
                if args.skip_systematics or not is_background_sample(sample, analyzer_cfg):
                    histos_sys[var_name][sample][sys] = histos[var_name][sample].Clone(var_name+'_'+sample+'_'+sys)
                    histos_sys[var_name][sample][sys].SetDirectory(0)
                plot_cfg.SetHistStyles(histos_sys[var_name][sample][sys], sample)
                histos_sys[var_name][sample][sys].Write()

    if args.hist_only:
        _validate_and_publish_output(out_file, tmp_output_path, final_output_path)
        elapsed = time.time() - start_time
        print("[histOnly] Wrote raw_plots/sys_dir only; skipped stack/PDF drawing.")
        print(f"Total runtime: {time.strftime('%H:%M:%S', time.gmtime(elapsed))}")
        print('Done')
        return

    ### save stack plots and make ratio plots
    out_file.cd()
    lumi_label = MakeLumiLabel(plot_cfg.lumi)
    cms_label  = MakeCMSDASLabel()

    scaled_sig = {}
    for var_name in var_names:
        # 依變數選擇縮放策略：H_m 用 sideband，其他維持原本全區間的 ScaleBkgToData
        if var_name == 'H_m':
            scale_factor = SideBandScaleBkgToData(histos, histos_sys, analyzer_cfg, signal_low=115., signal_high=135.)
        else:
            scale_factor = ScaleBkgToData(histos[var_name], analyzer_cfg, histos_sys.get(var_name))
        if scale_factor != 1.0:
            print(f"[Info] Applied background scaling ({var_name}) = {scale_factor:.4f}")

        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        
        if stacks['all'].GetStack().GetEntries() == 0:
            stack_entry = stacks['all'].GetStack().GetEntries()
            print(f"stack_entry: {stack_entry}")
            print(f"[Warning] Stack for {var_name} is empty. Skipping drawing.")
            continue  # 空 stack 直接跳過後續

        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['Data'], stacks['all'].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

        #### uncertainty graph
        total_unc = Total_Unc(stacks['bkg'], histos_sys[var_name], analyzer_cfg)

        if args.ln:
            canv_log = CreateCanvas(var_name+'_log')
            DrawOnCanv(canv_log, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, args.cut, args.mA, logY=True)

            canv_log.Write()
            SaveCanvPic(canv_log, analyzer_cfg.plot_output_path, _pdf_output_name(var_name) + '_log')
        else:
            canv = CreateCanvas(var_name)
            DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, args.cut, args.mA, logY=False)

            canv.Write()
            SaveCanvPic(canv, analyzer_cfg.plot_output_path, _pdf_output_name(var_name))

        if args.mva and _is_full_range_mva_larger_var(var_name, target_masses):
            histos_cdf, histos_sys_cdf = _build_cdf_histo_maps(
                var_name,
                histos[var_name],
                histos_sys[var_name],
                analyzer_cfg,
            )
            cdf_var_name = var_name + "_cdf"
            stacks_cdf = MakeStack(histos_cdf, analyzer_cfg, cdf_var_name)

            if stacks_cdf['all'].GetStack().GetEntries() == 0:
                stack_entry = stacks_cdf['all'].GetStack().GetEntries()
                print(f"stack_entry: {stack_entry}")
                print(f"[Warning] CDF stack for {var_name} is empty. Skipping drawing.")
                continue

            scaled_sig_cdf = {}
            for sample in analyzer_cfg.sig_names:
                scaled_sig_cdf[sample] = ScaleSignal(plot_cfg, stacks_cdf[sample], histos_cdf[sample], cdf_var_name)

            ratio_plot_cdf = MakeRatioPlot(histos_cdf['Data'], stacks_cdf['all'].GetStack().Last(), cdf_var_name)
            legend_cdf = MakeLegend(plot_cfg, histos_cdf, scaled_sig_cdf)
            total_unc_cdf = Total_Unc(stacks_cdf['bkg'], histos_sys_cdf, analyzer_cfg)
            cdf_y_axis_title = "Entries / %.3f" % histos_cdf['Data'].GetBinWidth(1)

            if args.ln:
                canv_cdf = CreateCanvas(cdf_var_name + "_log")
                DrawOnCanv(
                    canv_cdf,
                    cdf_var_name + "_log",
                    plot_cfg,
                    stacks_cdf,
                    histos_cdf,
                    scaled_sig_cdf,
                    ratio_plot_cdf,
                    legend_cdf,
                    lumi_label,
                    cms_label,
                    total_unc_cdf,
                    args.cut,
                    args.mA,
                    logY=True,
                    y_axis_title=cdf_y_axis_title,
                    axis_var_name=var_name,
                    x_axis_title="Transformed BDT Score",
                    log_y_max_scale=1.5e4,
                )
                canv_cdf.Write()
                SaveCanvPic(canv_cdf, analyzer_cfg.plot_output_path, _pdf_output_name(var_name) + "_cdf_log")
            else:
                canv_cdf = CreateCanvas(cdf_var_name)
                DrawOnCanv(
                    canv_cdf,
                    cdf_var_name,
                    plot_cfg,
                    stacks_cdf,
                    histos_cdf,
                    scaled_sig_cdf,
                    ratio_plot_cdf,
                    legend_cdf,
                    lumi_label,
                    cms_label,
                    total_unc_cdf,
                    args.cut,
                    args.mA,
                    logY=False,
                    y_axis_title=cdf_y_axis_title,
                    axis_var_name=var_name,
                    x_axis_title="Transformed BDT Score",
                )
                canv_cdf.Write()
                SaveCanvPic(canv_cdf, analyzer_cfg.plot_output_path, _pdf_output_name(var_name) + "_cdf")


    
    print('\n\n')
    #CountYield(analyzer_cfg, histos['ALP_m'])
    _validate_and_publish_output(out_file, tmp_output_path, final_output_path)

    elapsed = time.time() - start_time  # NEW
    print(f"Total runtime: {time.strftime('%H:%M:%S', time.gmtime(elapsed))}")  # NEW
    print('Done')


main()
