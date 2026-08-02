#!/usr/bin/env python3
"""Fast histogram backend for 1_prepare_dataVmc.py.

This module keeps the ROOT output contract of 1_prepare_dataVmc.py:

  raw_plots/<var>_<sample>
  sys_dir/<var>_<sample>_<systematic>

The slow script fills those histograms in a Python event loop.  This backend
books the same fills as ROOT RDataFrame actions so the event loop runs in C++.
It is intentionally restricted to the histogram-production mode used by the
Condor workflow.  Unsupported options raise FastBackendUnsupported so the
caller can fall back to the original loop implementation.
"""

from __future__ import annotations

import errno
import os
import shutil
import sys
import time
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir))
PROJECT_DIR = os.path.abspath(os.path.join(PLOT_DIR, os.pardir))

for path in (
    os.path.join(PLOT_DIR, "lib"),
    os.path.join(PROJECT_DIR, "HZaMVA", "scripts"),
):
    if path not in sys.path:
        sys.path.insert(0, path)

from root_compat import import_pyroot

ROOT = import_pyroot()

import Analyzer_Configs as AC
import Plot_Configs as PC
from Analyzer_Helper import getMassSigma
from Plot_Helper import LoadNtuples


ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gInterpreter.Declare("#include <cmath>")


class FastBackendUnsupported(RuntimeError):
    """Raised when the fast backend should defer to the original loop."""


@dataclass(frozen=True)
class HistSpec:
    nominal_col: str
    sys_col: str
    nbins: int
    xmin: float
    xmax: float
    extra_cut: str = "1"


MASS_VALUES = {
    "M1": 1.0,
    "M2": 2.0,
    "M3": 3.0,
    "M4": 4.0,
    "M5": 5.0,
    "M6": 6.0,
    "M7": 7.0,
    "M8": 8.0,
    "M9": 9.0,
    "M10": 10.0,
    "M15": 15.0,
    "M20": 20.0,
    "M25": 25.0,
    "M30": 30.0,
}

MASS_VALUES_FOR_PARAM = list(MASS_VALUES.values())
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
SIDEBAND_REWEIGHT_UNC_SYS_NAMES = ("sideband_rwgt_reweight_up", "sideband_rwgt_reweight_down")


def _clean_tag(tag: object) -> Optional[str]:
    if tag is None:
        return None
    tag = os.path.basename(str(tag).strip().strip("/"))
    return tag or None


def _apply_output_tag(analyzer_cfg: AC.Analyzer_Config, output_tag: object) -> None:
    tag = _clean_tag(output_tag)
    if not tag:
        return
    root_base, root_ext = os.path.splitext(analyzer_cfg.root_output_name)
    analyzer_cfg.root_output_name = f"{root_base}_{tag}{root_ext or '.root'}"
    analyzer_cfg.plot_output_path = os.path.join(analyzer_cfg.plot_output_path, tag)


def _validate_and_publish_output(out_file, tmp_output_path: str, final_output_path: str) -> None:
    out_file.Close()

    check_file = ROOT.TFile.Open(tmp_output_path, "READ")
    if not check_file or check_file.IsZombie():
        raise RuntimeError(f"[FastBackend][ERROR] unreadable temporary ROOT file: {tmp_output_path}")

    keys = check_file.GetListOfKeys()
    n_keys = keys.GetEntries() if keys else 0
    check_file.Close()
    if n_keys <= 0:
        raise RuntimeError(f"[FastBackend][ERROR] temporary ROOT file has no keys: {tmp_output_path}")

    try:
        os.replace(tmp_output_path, final_output_path)
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise
        shutil.move(tmp_output_path, final_output_path)
    print(f"[FastBackend] Published ROOT file: {final_output_path} ({n_keys} top-level keys)")


def _has_branch(chain, branch: str) -> bool:
    branches = chain.GetListOfBranches() if chain else None
    return bool(branches and branches.FindObject(branch))


def _is_background_sample(sample: str, analyzer_cfg: AC.Analyzer_Config) -> bool:
    return sample in analyzer_cfg.bkg_names


def _is_mc_sample(sample: str, analyzer_cfg: AC.Analyzer_Config) -> bool:
    return sample in analyzer_cfg.bkg_names or sample in analyzer_cfg.sig_names


def _systematic_central_branch(sys_name: str) -> Optional[str]:
    if sys_name.endswith("_up"):
        return sys_name[:-3] + "_central"
    if sys_name.endswith("_down"):
        return sys_name[:-5] + "_central"
    return None


def _safe_sys_weight_expr(nominal_col: str, sys_name: str, central_branch: str) -> str:
    return (
        f"(std::isfinite(static_cast<double>({sys_name})) && "
        f"std::isfinite(static_cast<double>({central_branch})) && "
        f"std::abs(static_cast<double>({central_branch})) >= 1e-12 && "
        f"std::isfinite(static_cast<double>({sys_name}) / static_cast<double>({central_branch}))"
        f" ? ({nominal_col}) * static_cast<double>({sys_name}) / static_cast<double>({central_branch})"
        f" : ({nominal_col}))"
    )


def _resolve_mva_branch_for_mass(chain, mass_tag: str) -> Optional[str]:
    candidates = (
        f"MVA_Score_mA_{mass_tag}",
        f"MVA_Score_{mass_tag}",
        f"MVA_score_{mass_tag}",
        f"BDTScore_{mass_tag}",
        f"BDT_{mass_tag}",
        "MVA_Score",
    )
    for name in candidates:
        if _has_branch(chain, name):
            return name
    return None


def _mva_col(mass_tag: str) -> str:
    return f"rdf_mva_score_{mass_tag}"


def _sys_weight_col(sys_name: str) -> str:
    return "rdf_weight_" + "".join(ch if ch.isalnum() else "_" for ch in sys_name)


def _param_cycle_expr() -> str:
    n_masses = len(MASS_VALUES_FOR_PARAM)
    terms = [
        f"((rdfentry_ % {n_masses}) == {idx}) * {mass:.8g}"
        for idx, mass in enumerate(MASS_VALUES_FOR_PARAM)
    ]
    return "(" + " + ".join(terms) + ")"


def _region_cut(region: int) -> str:
    base = "(H_m > -90.0 && H_m >= 95.0 && H_m <= 180.0)"
    if region == 1:
        return f"({base} && H_m <= 135.0 && H_m >= 115.0)"
    if region == 2:
        return f"({base} && !(H_m < 135.0 && H_m > 115.0))"
    return base


def _channel_cut(args) -> str:
    cuts: List[str] = []
    if args.ele:
        cuts.append("std::abs(z_mumu) != 1")
    if args.mu:
        cuts.append("std::abs(z_ee) != 1")
    return " && ".join(cuts) if cuts else "1"


def _sigma_cut(sigma_low: Mapping[str, float], sigma_hig: Mapping[str, float], mass: str, scale: float) -> str:
    low = 125.0 + float(sigma_low[mass]) * scale
    high = 125.0 + float(sigma_hig[mass]) * scale
    return f"(H_m < {high:.8g} && H_m > {low:.8g})"


def _materialize_extra_cut(
    extra_cut: str,
    var_name: str,
    sigma_low: Mapping[str, float],
    sigma_hig: Mapping[str, float],
) -> str:
    if not extra_cut.startswith("__SIGMA_"):
        return extra_cut
    mass = var_name.split("_")[-1]
    scale_key = extra_cut[len("__SIGMA_") : -2]
    if scale_key == "1P5":
        scale = 1.5
    else:
        scale = float(scale_key)
    return _sigma_cut(sigma_low, sigma_hig, mass, scale)


def _is_full_range_mva_var(var_name: str, target_masses: Sequence[str]) -> bool:
    return var_name in [f"mvaVal_{mass}" for mass in target_masses] or var_name in [
        f"mvaVal_larger_{mass}" for mass in target_masses
    ]


def _blind_cut(var_name: str, value_col: str, target_masses: Sequence[str]) -> str:
    if var_name == "H_m":
        return "!(H_m < 135.0 && H_m > 115.0)"
    if _is_full_range_mva_var(var_name, target_masses):
        return f"({value_col} < {MVA_FULL_RANGE_DATA_BLIND_THRESHOLD:.8g})"
    return "1"


def _build_var_names(plot_cfg: PC.Plot_Config, target_masses: Sequence[str], use_mva: bool) -> List[str]:
    var_names = list(plot_cfg.var_title_map.keys())
    if not use_mva:
        return var_names

    mva_region_vars: List[str] = []
    for mass in target_masses:
        var_names.append(f"mvaVal_{mass}")
        var_names.append(f"mvaVal_larger_{mass}")
        for region in ("1sigma", "1P5sigma", "2sigma", "3sigma"):
            mva_region_vars.append(f"mvaVal_{region}_{mass}")
            mva_region_vars.append(f"mvaVal_larger_{region}_{mass}")
    var_names.extend(mva_region_vars)
    return var_names


def _build_hist_specs(var_names: Iterable[str], target_masses: Sequence[str]) -> Dict[str, HistSpec]:
    base_specs = {
        "pho1Pt": HistSpec("ALP_lead_photon_pt", "pho1Pt", 42, 8.0, 50.0),
        # 2026-08-02: complete the 06-03 var rename to p_T/m_llgg ratios. var_title_map uses
        # the *_oHm keys; the branches exist in the scored ntuples. Binning mirrors the legacy
        # 1_prepare_dataVmc.py so data_vmc_fast (prepare) and 2_plot_dataVmc (draw) agree.
        "pho1Pt_oHm": HistSpec("pho1Pt_oHm", "pho1Pt_oHm", 30, 0.0, 0.6),
        "pho1eta": HistSpec("ALP_lead_photon_eta", "ALP_lead_photon_eta", 30, -3.0, 3.0),
        "pho1phi": HistSpec("ALP_lead_photon_phi", "ALP_lead_photon_phi", 40, -4.0, 4.0),
        "pho1R9": HistSpec("ALP_lead_photon_r9", "ALP_lead_photon_r9", 25, 0.1, 1.0),
        "pho1IetaIeta55": HistSpec("ALP_lead_photon_sieie", "ALP_lead_photon_sieie", 25, 0.0, 0.08),
        "pho1ECALIso": HistSpec(
            "ALP_lead_photon_ecalPFClusterIso",
            "ALP_lead_photon_ecalPFClusterIso",
            50,
            0.0,
            50.0,
        ),
        "pho1CIso": HistSpec("ALP_lead_photon_chiso", "ALP_lead_photon_chiso", 35, 0.0, 0.35),
        "pho1HCALIso": HistSpec(
            "ALP_lead_photon_hcalPFClusterIso",
            "ALP_lead_photon_hcalPFClusterIso",
            25,
            0.0,
            4.0,
        ),
        "pho1HOE": HistSpec("ALP_lead_photon_hoe_PUcorr", "ALP_lead_photon_hoe_PUcorr", 25, 0.0, 0.01),
        "pho2Pt": HistSpec("ALP_sublead_photon_pt", "ALP_sublead_photon_pt", 22, 8.0, 30.0),
        "pho2Pt_oHm": HistSpec("pho2Pt_oHm", "pho2Pt_oHm", 25, 0.0, 0.5),
        "pho2eta": HistSpec("ALP_sublead_photon_eta", "ALP_sublead_photon_eta", 30, -3.0, 3.0),
        "pho2phi": HistSpec("ALP_sublead_photon_phi", "ALP_sublead_photon_phi", 40, -4.0, 4.0),
        "pho2R9": HistSpec("ALP_sublead_photon_r9", "ALP_sublead_photon_r9", 25, 0.1, 1.0),
        "pho2IetaIeta55": HistSpec("ALP_sublead_photon_sieie", "ALP_sublead_photon_sieie", 25, 0.0, 0.08),
        "pho2ECALIso": HistSpec(
            "ALP_sublead_photon_ecalPFClusterIso",
            "ALP_sublead_photon_ecalPFClusterIso",
            50,
            0.0,
            50.0,
        ),
        "pho2CIso": HistSpec("ALP_sublead_photon_chiso", "ALP_sublead_photon_chiso", 35, 0.0, 0.35),
        "pho2HCALIso": HistSpec(
            "ALP_sublead_photon_hcalPFClusterIso",
            "ALP_sublead_photon_hcalPFClusterIso",
            25,
            0.0,
            4.0,
        ),
        "pho2HOE": HistSpec("ALP_sublead_photon_hoe_PUcorr", "ALP_sublead_photon_hoe_PUcorr", 25, 0.0, 0.1),
        "Z_m": HistSpec("Z_mass", "Z_mass", 80, 50.0, 130.0),
        "H_m": HistSpec("H_m", "H_m", 85, 95.0, 180.0),
        "H_pt": HistSpec("H_pt", "H_pt", 80, 0.0, 160.0),
        "H_pt_oHm": HistSpec("H_pt_oHm", "H_pt_oHm", 50, 0.0, 2.5),
        "ALP_m": HistSpec("ALP_m", "ALP_m", 40, 0.0, 40.0),
        "var_dR_g1g2": HistSpec("var_dR_g1g2", "var_dR_g1g2", 25, 0.0, 5.0),
        "var_PtaOverMa": HistSpec("var_PtaOverMa", "var_PtaOverMa", 25, 0.0, 100.0),
        "var_dR_Za": HistSpec("var_dR_Za", "var_dR_Za", 35, 0.0, 7.0),
        "var_dR_g1Z": HistSpec("var_dR_g1Z", "var_dR_g1Z", 35, 0.0, 7.0),
        "var_PtaOverMh": HistSpec("var_PtaOverMh", "var_PtaOverMh", 25, 0.0, 0.75),
        "var_Pta": HistSpec("var_Pta", "var_Pta", 60, 0.0, 60.0),
        "var_MhMa": HistSpec("var_MhMa", "var_MhMa", 100, 100.0, 200.0),
        "var_MhMZ": HistSpec("var_MhMZ", "var_MhMZ", 130, 180.0, 310.0),
        "ALP_calculatedPhotonIso": HistSpec(
            "ALP_calculatedPhotonIso",
            "ALP_calculatedPhotonIso",
            25,
            0.0,
            125.0,
        ),
        "param": HistSpec("rdf_param_value", "rdf_param_value", 25, -0.3, 0.6),
        # Derived BDT feature (data/MC comparison)
        "pho_pt_asym": HistSpec("rdf_pho_pt_asym", "rdf_pho_pt_asym", 40, 0.0, 1.0),
    }

    specs = {var: base_specs[var] for var in var_names if var in base_specs}
    for mass in target_masses:
        mva_col = _mva_col(mass)
        specs[f"mvaVal_{mass}"] = HistSpec(mva_col, mva_col, 240, -0.1, 1.1)
        specs[f"mvaVal_larger_{mass}"] = HistSpec(mva_col, mva_col, MVA_LARGER_NBINS, 0.0, 1.0)
        for region, scale in SIGMA_REGION_SCALES.items():
            extra_cut = "__SIGMA_%s__" % ("1P5" if scale == 1.5 else int(scale))
            specs[f"mvaVal_{region}_{mass}"] = HistSpec(mva_col, mva_col, 240, -0.1, 1.1, extra_cut)
            specs[f"mvaVal_larger_{region}_{mass}"] = HistSpec(
                mva_col,
                mva_col,
                MVA_LARGER_NBINS,
                0.0,
                1.0,
                extra_cut,
            )
    return specs


def _book_hist(df, name: str, spec: HistSpec, value_col: str, weight_col: str):
    model = ROOT.RDF.TH1DModel(name, name, spec.nbins, spec.xmin, spec.xmax)
    return df.Histo1D(model, value_col, weight_col)


def _to_th1f(result_ptr, name: str, spec: HistSpec):
    source = result_ptr.GetValue()
    source.SetName(source.GetName() + "__rdf_source")
    old = ROOT.gROOT.FindObject(name)
    if old:
        old.Delete()
    hist = ROOT.TH1F(name, name, spec.nbins, spec.xmin, spec.xmax)
    hist.Sumw2()
    for i_bin in range(0, spec.nbins + 2):
        hist.SetBinContent(i_bin, source.GetBinContent(i_bin))
        hist.SetBinError(i_bin, source.GetBinError(i_bin))
    hist.SetEntries(source.GetEntries())
    hist.SetDirectory(0)
    return hist


def _empty_hist(name: str, spec: HistSpec):
    old = ROOT.gROOT.FindObject(name)
    if old:
        old.Delete()
    hist = ROOT.TH1F(name, name, spec.nbins, spec.xmin, spec.xmax)
    hist.Sumw2()
    hist.SetDirectory(0)
    return hist


def _clone_hist(hist, name: str):
    clone = hist.Clone(name)
    clone.SetNameTitle(name, name)
    clone.SetDirectory(0)
    return clone


def _column_is_available(chain, defined_cols: Iterable[str], col: str) -> bool:
    return col in defined_cols or _has_branch(chain, col)


def _require_branch(chain, sample: str, branch: str, purpose: str) -> None:
    if not _has_branch(chain, branch):
        raise FastBackendUnsupported(f"sample={sample} is missing branch '{branch}' needed for {purpose}")


def _nominal_weight_expr(chain, sample: str, analyzer_cfg: AC.Analyzer_Config, args) -> str:
    if args.use_sideband_reweight and _is_mc_sample(sample, analyzer_cfg):
        if _has_branch(chain, "weight_sideband_rwgt"):
            return "weight_sideband_rwgt"
        raise FastBackendUnsupported(
            "fast backend needs branch weight_sideband_rwgt for --useSidebandReweight; "
            "falling back to the original JSON-aware loop"
        )
    if _has_branch(chain, "weight"):
        return "weight"
    return "1.0"


def _prepare_dataframe(
    chain,
    sample: str,
    analyzer_cfg: AC.Analyzer_Config,
    target_masses: Sequence[str],
    args,
) -> Tuple[object, set, Dict[str, Optional[str]]]:
    _require_branch(chain, sample, "H_m", "region selection")
    if args.ele:
        _require_branch(chain, sample, "z_mumu", "--ele channel veto")
    if args.mu:
        _require_branch(chain, sample, "z_ee", "--mu channel veto")

    df = ROOT.RDataFrame(chain)
    if args.max_events and args.max_events > 0:
        df = df.Range(int(args.max_events))

    defined_cols = {"rdf_nominal_weight"}
    df = df.Define("rdf_nominal_weight", _nominal_weight_expr(chain, sample, analyzer_cfg, args))

    if sample in MASS_VALUES:
        param_expr = f"(ALP_m - {MASS_VALUES[sample]:.8g}) / H_m"
    else:
        param_expr = f"(ALP_m - {_param_cycle_expr()}) / H_m"
    df = df.Define("rdf_param_value", param_expr)
    defined_cols.add("rdf_param_value")

    # Derived BDT feature for data/MC plots (only pho_pt_asym kept after corr study).
    # Lead pT >= sublead pT by construction, so the signed form lies in [0, 1].
    df = df.Define(
        "rdf_pho_pt_asym",
        "((double)ALP_lead_photon_pt - (double)ALP_sublead_photon_pt) "
        "/ ((double)ALP_lead_photon_pt + (double)ALP_sublead_photon_pt + 1e-6)",
    )
    defined_cols.add("rdf_pho_pt_asym")

    mva_branches: Dict[str, Optional[str]] = {}
    for mass in target_masses:
        branch = _resolve_mva_branch_for_mass(chain, mass)
        mva_branches[mass] = branch
        expr = branch if branch else "-1.0"
        df = df.Define(_mva_col(mass), expr)
        defined_cols.add(_mva_col(mass))
        print(
            "[FastBackend][MVA] sample=%12s mass=%3s branch=%s"
            % (sample, mass, branch if branch else "<missing, fill -1>")
        )

    cuts = [_region_cut(args.region), _channel_cut(args)]
    if args.cut:
        if args.mA not in target_masses:
            raise FastBackendUnsupported(f"--cut requested {args.mA}, but it is not in target_masses")
        cuts.append(f"({_mva_col(args.mA)} >= {float(analyzer_cfg.mvaCut[args.mA]):.8g})")

    df = df.Filter(" && ".join(f"({cut})" for cut in cuts if cut and cut != "1"))
    return df, defined_cols, mva_branches


def _book_sample_actions(
    df,
    chain,
    sample: str,
    analyzer_cfg: AC.Analyzer_Config,
    plot_cfg: PC.Plot_Config,
    var_names: Sequence[str],
    hist_specs: Mapping[str, HistSpec],
    sigma_low: Mapping[str, float],
    sigma_hig: Mapping[str, float],
    target_masses: Sequence[str],
    defined_cols: set,
    args,
):
    nominal_results = {}

    for var_name in var_names:
        spec = hist_specs[var_name]
        value_col = spec.nominal_col
        if not _column_is_available(chain, defined_cols, value_col):
            raise FastBackendUnsupported(f"sample={sample} is missing column '{value_col}' for {var_name}")

        cut_parts = [
            _materialize_extra_cut(spec.extra_cut, var_name, sigma_low, sigma_hig),
        ]
        if args.blind and sample == "Data":
            cut_parts.append(_blind_cut(var_name, value_col, target_masses))

        hist_df = df
        cut_expr = " && ".join(f"({cut})" for cut in cut_parts if cut and cut != "1")
        if cut_expr:
            hist_df = hist_df.Filter(cut_expr)

        nominal_results[var_name] = _book_hist(
            hist_df,
            f"{var_name}_{sample}",
            spec,
            value_col,
            "rdf_nominal_weight",
        )

    sys_results = {}
    if not args.skip_systematics and _is_background_sample(sample, analyzer_cfg):
        df_sys = df
        for sys_name in analyzer_cfg.sys_names:
            weight_col = _sys_weight_col(sys_name)
            central_branch = _systematic_central_branch(sys_name)
            if central_branch and _has_branch(chain, sys_name) and _has_branch(chain, central_branch):
                weight_expr = _safe_sys_weight_expr("rdf_nominal_weight", sys_name, central_branch)
            else:
                print(f"[FastBackend][Sys] sample={sample} sys={sys_name}: missing varied/central branch; use nominal")
                weight_expr = "rdf_nominal_weight"
            df_sys = df_sys.Define(weight_col, weight_expr)
            defined_cols.add(weight_col)

        for var_name in var_names:
            spec = hist_specs[var_name]
            value_col = spec.sys_col
            if not _column_is_available(chain, defined_cols, value_col):
                value_col = spec.nominal_col
            if not _column_is_available(chain, defined_cols, value_col):
                raise FastBackendUnsupported(f"sample={sample} is missing sys column for {var_name}")

            cut_expr = _materialize_extra_cut(spec.extra_cut, var_name, sigma_low, sigma_hig)
            hist_df = df_sys.Filter(cut_expr) if cut_expr and cut_expr != "1" else df_sys

            sys_results[var_name] = {}
            for sys_name in analyzer_cfg.sys_names:
                sys_results[var_name][sys_name] = _book_hist(
                    hist_df,
                    f"{var_name}_{sample}_{sys_name}",
                    spec,
                    value_col,
                    _sys_weight_col(sys_name),
                )

    return nominal_results, sys_results


def _configure_analyzer(args, parse_sample_filter: Callable, format_source_selectors: Callable):
    analyzer_cfg = AC.Analyzer_Config("inclusive", args.year, args.region, args.mva)
    analyzer_cfg.mva = bool(args.mva)
    analyzer_cfg.mva_alp_mass = str(args.mA) if args.mva else "M1"
    all_sig_names = list(analyzer_cfg.sig_names)

    if args.use_sideband_reweight and args.sideband_reweight_unc:
        raise FastBackendUnsupported("sideband reweight uncertainty uses Python JSON callbacks; use loop backend")
    if args.use_sideband_reweight:
        analyzer_cfg.sys_names = [
            sys_name
            for sys_name in analyzer_cfg.sys_names
            if sys_name not in SIDEBAND_REWEIGHT_UNC_SYS_NAMES
        ]

    _apply_output_tag(analyzer_cfg, args.output_tag)

    if args.mva and analyzer_cfg.year in ["run3", "run3_NFlow"]:
        target_masses = list(all_sig_names)
    else:
        target_masses = [analyzer_cfg.mva_alp_mass] if args.mva else []

    if args.cut:
        analyzer_cfg.sig_names = [args.mA]
        analyzer_cfg.samp_names = analyzer_cfg.bkg_names + analyzer_cfg.sig_names + ["Data"]
        if args.mva:
            target_masses = [args.mA]

    sample_source_filters = {}
    if args.samples:
        selected_samples, sample_source_filters = parse_sample_filter(args.samples, analyzer_cfg)
        if selected_samples:
            selected_set = set(selected_samples)
            analyzer_cfg.samp_names = [s for s in analyzer_cfg.samp_names if s in selected_set]
            if not analyzer_cfg.samp_names:
                raise RuntimeError(f"No samples selected after applying --samples={args.samples}")
            print("[FastBackend][Samples] Running selected samples:", analyzer_cfg.samp_names)
            for sample in analyzer_cfg.samp_names:
                if sample in sample_source_filters:
                    print(
                        "[FastBackend][Samples] sample=%s source filters: %s"
                        % (sample, format_source_selectors(sample_source_filters[sample]))
                    )

    analyzer_cfg.run3_source_filters = sample_source_filters
    return analyzer_cfg, target_masses


def run_fast_prepare(args, parse_sample_filter: Callable, format_source_selectors: Callable) -> None:
    if not args.hist_only:
        raise FastBackendUnsupported("fast backend currently writes histogram ROOT files only; use --histOnly")

    threads = int(os.environ.get("ALP_RDF_THREADS", "0") or "0")
    if threads > 0:
        ROOT.ROOT.EnableImplicitMT(threads)
        print(f"[FastBackend] Enabled ROOT implicit MT with {threads} threads")

    start_time = time.time()
    analyzer_cfg, target_masses = _configure_analyzer(args, parse_sample_filter, format_source_selectors)

    os.makedirs(analyzer_cfg.out_dir, exist_ok=True)
    final_output_path = os.path.join(analyzer_cfg.out_dir, analyzer_cfg.root_output_name)
    tmp_output_dir = os.environ.get("ALP_OUTPUT_TMP_DIR", analyzer_cfg.out_dir)
    os.makedirs(tmp_output_dir, exist_ok=True)
    tmp_output_path = os.path.join(tmp_output_dir, ".%s.%d.tmp" % (analyzer_cfg.root_output_name, os.getpid()))
    if os.path.exists(tmp_output_path):
        os.remove(tmp_output_path)

    print("[FastBackend] Writing temporary ROOT file:", tmp_output_path)
    print("[FastBackend] Final ROOT file:", final_output_path)
    analyzer_cfg.Print_Config()

    ntuples = LoadNtuples(analyzer_cfg)
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)
    var_names = _build_var_names(plot_cfg, target_masses, args.mva)
    hist_specs = _build_hist_specs(var_names, target_masses)
    var_names = [var for var in var_names if var in hist_specs]
    sigma_low, sigma_hig = getMassSigma(analyzer_cfg)

    histos = {var_name: {} for var_name in var_names}
    histos_sys = {var_name: {} for var_name in var_names}

    for sample in analyzer_cfg.samp_names:
        chain = ntuples[sample]
        print(f"\n[FastBackend][Sample] {sample}: entries={chain.GetEntries()}")
        df, defined_cols, _ = _prepare_dataframe(chain, sample, analyzer_cfg, target_masses, args)
        nominal_results, sys_results = _book_sample_actions(
            df,
            chain,
            sample,
            analyzer_cfg,
            plot_cfg,
            var_names,
            hist_specs,
            sigma_low,
            sigma_hig,
            target_masses,
            defined_cols,
            args,
        )

        for var_name in var_names:
            hist = _to_th1f(nominal_results[var_name], f"{var_name}_{sample}", hist_specs[var_name])
            plot_cfg.SetHistStyles(hist, sample)
            histos[var_name][sample] = hist

        for var_name in var_names:
            histos_sys[var_name][sample] = {}
            for sys_name in analyzer_cfg.sys_names:
                if var_name in sys_results and sys_name in sys_results[var_name]:
                    hist_sys = _to_th1f(
                        sys_results[var_name][sys_name],
                        f"{var_name}_{sample}_{sys_name}",
                        hist_specs[var_name],
                    )
                else:
                    hist_sys = _clone_hist(histos[var_name][sample], f"{var_name}_{sample}_{sys_name}")
                plot_cfg.SetHistStyles(hist_sys, sample)
                histos_sys[var_name][sample][sys_name] = hist_sys

        if args.mva:
            print("[FastBackend][MVA] Per-mass histogram summary for sample:", sample)
            for mass in target_masses:
                hname = f"mvaVal_{mass}"
                if hname in histos and sample in histos[hname]:
                    hist = histos[hname][sample]
                    print(
                        "  mass=%3s sample=%12s Entries=%8.0f Integral=%.3f"
                        % (mass, sample, hist.GetEntries(), hist.Integral())
                    )

    out_file = ROOT.TFile(tmp_output_path, "RECREATE")
    if not out_file or out_file.IsZombie():
        raise RuntimeError(f"[FastBackend][ERROR] cannot create temporary ROOT file: {tmp_output_path}")

    raw_dir = out_file.mkdir("raw_plots")
    raw_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            histos[var_name][sample].Write()

    sys_dir = out_file.mkdir("sys_dir")
    sys_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            for sys_name in analyzer_cfg.sys_names:
                histos_sys[var_name][sample][sys_name].Write()

    _validate_and_publish_output(out_file, tmp_output_path, final_output_path)
    elapsed = time.time() - start_time
    print("[FastBackend] Wrote raw_plots/sys_dir only; skipped stack/PDF drawing.")
    print(f"Total runtime: {time.strftime('%H:%M:%S', time.gmtime(elapsed))}")
    print("Done")
