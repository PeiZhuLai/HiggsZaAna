#!/usr/bin/env python3
#
# Faster version of plot_variable_dataVmc.py.
#
# The slow script fills every histogram in a Python event loop.  This version
# uses TChain.Draw/TTreeFormula so the event scan runs inside ROOT.

import argparse
import os
import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import ROOT

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, os.path.join(PLOT_DIR, "lib"))
sys.path.insert(0, "%s/lib" % os.getcwd())

import Analyzer_Configs as AC
import Plot_Configs as PC
from Analyzer_Helper import getMassSigma
from Plot_Helper import (
    CreateCanvas,
    Draw_unc,
    Get_StatUnc,
    LoadNtuples,
    MakeRatioPlot,
    MakeStack,
    SaveCanvPic,
    ScaleBkgToData,
    ScaleSignal,
    Total_Unc,
)


ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)


PDF_NAME_MAP = {
    "pho1Pt": "1_pho1Pt",
    "pho1R9": "2_pho1R9",
    "pho1IetaIeta55": "3_pho1IetaIeta55",
    "pho2Pt": "4_pho2Pt",
    "pho2R9": "5_pho2R9",
    "pho2IetaIeta55": "6_pho2IetaIeta55",
    "pho1ECALIso": "7_pho1ECALIso",
    "pho2ECALIso": "8_pho2ECALIso",
    "ALP_calculatedPhotonIso": "9_ALP_calculatedPhotonIso",
    "var_dR_Za": "10_var_dR_Za",
    "var_dR_g1g2": "11_var_dR_g1g2",
    "var_dR_g1Z": "12_var_dR_g1Z",
    "var_PtaOverMh": "13_var_PtaOverMh",
    "H_pt": "14_H_pt",
    "param": "15_param",
    "ALP_m": "16_ALP_m",
}


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

search_mA_list = [float(m) for m in range(1, 31)]


@dataclass(frozen=True)
class HistSpec:
    expr: str
    nbins: int
    xmin: float
    xmax: float
    extra_cut: str = "1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fast data/MC variable plotter")
    parser.add_argument("-y", "--Year", dest="year", default="run3")
    parser.add_argument("--region", dest="region", type=int, default=0)
    parser.add_argument("-m", "--mva", dest="mva", action="store_true", default=False)
    parser.add_argument("--cut", dest="cut", action="store_true", default=False)
    parser.add_argument("--cutVal", dest="cutVal", type=float, default=0.0)
    parser.add_argument("--mA", dest="mA", default="M5")
    parser.add_argument("-b", "--blind", dest="blind", action="store_true", default=False)
    parser.add_argument("-ln", "--ln", dest="ln", action="store_true", default=False)
    parser.add_argument("--ele", dest="ele", action="store_true", default=False)
    parser.add_argument("--mu", dest="mu", action="store_true", default=False)
    parser.add_argument("--input-dir", default=None, help="Override Analyzer_Configs sample_loc")
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Plot output directory. Default: <Plot>/plots/fast_plot_<region>",
    )
    parser.add_argument(
        "--skip-sys",
        action="store_true",
        help="Reuse central histograms for the uncertainty band instead of drawing all systematic variations.",
    )
    return parser.parse_args()


def output_name(var_name: str) -> str:
    return PDF_NAME_MAP.get(var_name, var_name)


def set_no_title_or_stats(obj) -> None:
    if obj:
        obj.SetTitle("")
        if hasattr(obj, "SetStats"):
            obj.SetStats(False)


def draw_cms_labels(canvas: ROOT.TCanvas, lumi: str) -> None:
    canvas.cd()
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextFont(42)

    label.SetTextAlign(13)
    label.SetTextSize(0.050)
    label.DrawLatex(0.16, 0.935, "#bf{CMS} #it{Preliminary}")

    label.SetTextAlign(31)
    label.SetTextSize(0.040)
    label.DrawLatex(0.95, 0.935, f"{float(lumi):.2f} fb^{{-1}} (13.6 TeV)")


def draw_fast_on_canvas(
    canvas: ROOT.TCanvas,
    var_name: str,
    plot_cfg: PC.Plot_Config,
    stacks,
    histos,
    scaled_sig,
    ratio_plot,
    total_unc,
    bdt_cut: bool,
    mA: str,
    log_y: bool,
) -> None:
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    canvas.SetBottomMargin(0.012)
    canvas.cd()

    left_margin = 0.16
    right_margin = 0.05
    bottom_margin = 0.19
    top_margin = 0.085

    upper_pad = ROOT.TPad("upperpad_" + var_name, "upperpad_" + var_name, 0, 0.25, 1, 1)
    upper_pad.SetLeftMargin(left_margin)
    upper_pad.SetRightMargin(right_margin)
    upper_pad.SetBottomMargin(bottom_margin)
    upper_pad.SetTopMargin(top_margin)
    upper_pad.Draw()
    upper_pad.cd()
    upper_pad.SetTickx(1)
    upper_pad.SetTicky(1)
    canvas.SetTickx()
    canvas.SetTicky()

    for hist in histos.values():
        set_no_title_or_stats(hist)
    for stack in stacks.values():
        set_no_title_or_stats(stack)

    h_max = max(histos["Data"].GetMaximum(), stacks["all"].GetMaximum(), stacks["sig"].GetMaximum())
    if log_y:
        upper_pad.SetLogy()
        histos["Data"].SetMinimum(1e-2)
        stacks["all"].SetMinimum(1e-2)
        histos["Data"].SetMaximum(h_max * 1.1e4 if h_max > 0 else 1.0)
        stacks["all"].SetMaximum(h_max * 1.1e4 if h_max > 0 else 1.0)
    else:
        histos["Data"].SetMaximum(h_max * 1.4 if h_max > 0 else 1.0)
        stacks["all"].SetMaximum(h_max * 1.4 if h_max > 0 else 1.0)

    histos["Data"].Draw("PE")
    histos["Data"].GetXaxis().SetLabelSize(0)
    histos["Data"].GetXaxis().SetTitleOffset(0.95)

    if var_name in ["H_m", "ALP_m", "Z_m"]:
        histos["Data"].GetYaxis().SetTitle("Events / %.2f GeV" % histos["Data"].GetBinWidth(1))
    else:
        histos["Data"].GetYaxis().SetTitle("Events / %.3f" % histos["Data"].GetBinWidth(1))
    histos["Data"].GetYaxis().SetTitleSize(0.07)
    histos["Data"].GetYaxis().SetLabelSize(0.055)
    histos["Data"].GetYaxis().SetTitleFont(42)
    histos["Data"].GetYaxis().SetTitleOffset(1.15)

    stacks["all"].Draw("HISTSAME")
    if var_name.split("_")[-1] in plot_cfg.ana_cfg.sig_names:
        scaled_sig[var_name.split("_")[-1]].Draw("HISTSAME")
    elif bdt_cut:
        scaled_sig[mA].Draw("HISTSAME")
    else:
        line_styles = [1, 2, 5, 7]
        for sample, style in zip(["M1", "M10", "M20", "M30"], line_styles):
            if sample in scaled_sig:
                scaled_sig[sample].SetLineStyle(style)
                scaled_sig[sample].SetLineWidth(4)
                scaled_sig[sample].Draw("HISTSAME")

    histos["Data"].SetMarkerStyle(20)
    histos["Data"].SetLineWidth(2)
    histos["Data"].GetXaxis().SetTickLength(0.04)

    stat_err, stat_err_norm = Get_StatUnc(stacks["bkg"].GetStack().Last())
    total_abs = total_unc[0]
    total_norm = total_unc[1]

    Draw_unc(total_abs, ROOT.TColor.GetColor("#A0A0A0"), alpha=0.3, fill_style=1001)
    Draw_unc(stat_err, ROOT.TColor.GetColor("#404040"), alpha=0.90, fill_style=3354)

    histos["Data"].Draw("SAMEPE")
    histos["Data"].Draw("AXIS SAME")

    if var_name.split("_")[-1] in plot_cfg.ana_cfg.sig_names:
        legend_1 = ROOT.TLegend(0.66, 0.59, 0.97, 0.86)
        ROOT.SetOwnership(legend_1, False)
        legend_1.AddEntry(histos["Data"], "Data", "PE")
        bkg_labels = {"DYGto2LG": "Z + #gamma", "DYJetsToLL": "Z + jets"}
        for sample_bkg in plot_cfg.ana_cfg.bkg_names:
            legend_1.AddEntry(histos[sample_bkg], bkg_labels.get(sample_bkg, sample_bkg), "f")
        legend_1.AddEntry(total_abs, "Total Unc.", "f")
        legend_1.AddEntry(stat_err, "Stat. Unc.", "f")

        legend_2 = ROOT.TLegend(0.40, 0.80, 0.66, 0.86)
        ROOT.SetOwnership(legend_2, False)
        legend_2.AddEntry(
            scaled_sig[var_name.split("_")[-1]],
            "m_{a} = %s GeV" % var_name.split("_")[-1].lstrip("M"),
            "l",
        )
        legends = (legend_1, legend_2)
    else:
        legend_1 = ROOT.TLegend(0.29, 0.816, 0.50, 0.87)
        ROOT.SetOwnership(legend_1, False)
        legend_1.AddEntry(histos["Data"], "Data", "PE")

        legend_2 = ROOT.TLegend(0.44, 0.654, 0.73, 0.87)
        ROOT.SetOwnership(legend_2, False)
        bkg_labels = {"DYGto2LG": "Z + #gamma", "DYJetsToLL": "Z + jets"}
        for sample_bkg in plot_cfg.ana_cfg.bkg_names:
            legend_2.AddEntry(histos[sample_bkg], bkg_labels.get(sample_bkg, sample_bkg), "f")
        legend_2.AddEntry(total_abs, "Total Unc.", "f")
        legend_2.AddEntry(stat_err, "Stat. Unc.", "f")

        legend_3 = ROOT.TLegend(0.67, 0.654, 0.97, 0.87)
        ROOT.SetOwnership(legend_3, False)
        if bdt_cut:
            legend_3.AddEntry(scaled_sig[mA], "m_{a} = %s GeV" % mA.lstrip("M"), "l")
        else:
            for sample in ["M1", "M10", "M20", "M30"]:
                if sample in scaled_sig:
                    legend_3.AddEntry(scaled_sig[sample], "m_{a} = %s GeV" % sample.lstrip("M"), "l")
        legends = (legend_1, legend_2, legend_3)

    for legend in legends:
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetFillColor(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.045)
        legend.Draw("SAME")

    draw_cms_labels(canvas, plot_cfg.lumi)

    canvas.cd()
    lower_pad = ROOT.TPad("lowerpad_" + var_name, "lowerpad_" + var_name, 0, 0, 1, 0.35)
    lower_pad.SetTopMargin(0.00001)
    lower_pad.SetBottomMargin(0.36)
    lower_pad.SetLeftMargin(left_margin)
    lower_pad.SetRightMargin(right_margin)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    lower_pad.SetTickx(1)
    lower_pad.SetTicky(1)

    if var_name.split("_")[-1] in plot_cfg.ana_cfg.sig_names:
        ratio_plot.GetXaxis().SetTitle("BDT Score")
    else:
        ratio_plot.GetXaxis().SetTitle(plot_cfg.var_title_map[var_name])
    ratio_plot.GetXaxis().SetTitleOffset(1.0)
    ratio_plot.Draw("APZ SAME")

    Draw_unc(total_norm, ROOT.TColor.GetColor("#A0A0A0"), alpha=0.3, fill_style=1001)
    Draw_unc(stat_err_norm, ROOT.TColor.GetColor("#404040"), alpha=0.90, fill_style=3354)
    ratio_plot.Draw("SAMEPZ")


def has_branch(chain: ROOT.TChain, branch: str) -> bool:
    branches = chain.GetListOfBranches()
    return bool(branches and branches.FindObject(branch))


def first_existing_branch(chain: ROOT.TChain, candidates: Sequence[str]) -> Optional[str]:
    for branch in candidates:
        if has_branch(chain, branch):
            return branch
    return None


def common_branch(ntuples: Dict[str, ROOT.TChain], samples: Iterable[str], candidates: Sequence[str]) -> str:
    for branch in candidates:
        if all(has_branch(ntuples[sample], branch) for sample in samples if sample in ntuples):
            return branch
    return candidates[0]


def resolve_mva_branch(chain: ROOT.TChain, mass_tag: str) -> str:
    return first_existing_branch(
        chain,
        (
            f"MVA_Score_mA_{mass_tag}",
            f"MVA_Score_{mass_tag}",
            f"MVA_score_{mass_tag}",
            f"BDTScore_{mass_tag}",
            f"BDT_{mass_tag}",
            "MVA_Score",
        ),
    ) or "MVA_Score"


def region_cut(region: int, mass_branch: str) -> str:
    base = f"({mass_branch} > 95.0 && {mass_branch} < 180.0)"
    if region == 1:
        return f"({base} && {mass_branch} > 115.0 && {mass_branch} < 135.0)"
    if region == 2:
        return f"({base} && !({mass_branch} > 115.0 && {mass_branch} < 135.0))"
    return base


def channel_cut(args: argparse.Namespace) -> str:
    cuts = []
    if args.ele:
        cuts.append("abs(z_mumu) != 1")
    if args.mu:
        cuts.append("abs(z_ee) != 1")
    return " && ".join(cuts) if cuts else "1"


def safe_div_weight(base_weight: str, varied_branch: str, central_branch: str) -> str:
    return f"({base_weight}) * ({varied_branch}) / ({central_branch})"


def sys_central_branch(sys_name: str) -> str:
    if sys_name.endswith("_up"):
        return sys_name[:-3] + "_central"
    if sys_name.endswith("_down"):
        return sys_name[:-5] + "_central"
    return sys_name


def param_cycle_expr() -> str:
    terms = [f"((Entry$ % {len(search_mA_list)}) == {idx}) * {mass:.8g}" for idx, mass in enumerate(search_mA_list)]
    return "(" + " + ".join(terms) + ")"


def build_hist_specs(
    var_names: Sequence[str],
    target_masses: Sequence[str],
    mva_branch_map: Dict[str, Dict[str, str]],
    mass_branch: str,
) -> Dict[str, HistSpec]:
    base_specs = {
        "pho1Pt": HistSpec("ALP_lead_photon_pt", 42, 8.0, 50.0),
        "pho1eta": HistSpec("ALP_lead_photon_eta", 30, -3.0, 3.0),
        "pho1phi": HistSpec("ALP_lead_photon_phi", 40, -4.0, 4.0),
        "pho1R9": HistSpec("ALP_lead_photon_r9", 25, 0.1, 1.0),
        "pho1IetaIeta55": HistSpec("ALP_lead_photon_sieie", 25, 0.0, 0.08),
        "pho1ECALIso": HistSpec("ALP_lead_photon_ecalPFClusterIso", 50, 0.0, 50.0),
        "pho1CIso": HistSpec("ALP_lead_photon_chiso", 35, 0.0, 0.35),
        "pho1HCALIso": HistSpec("ALP_lead_photon_hcalPFClusterIso", 25, 0.0, 4.0),
        "pho1HOE": HistSpec("ALP_lead_photon_hoe_PUcorr", 25, 0.0, 0.01),
        "pho2Pt": HistSpec("ALP_sublead_photon_pt", 22, 8.0, 30.0),
        "pho2eta": HistSpec("ALP_sublead_photon_eta", 30, -3.0, 3.0),
        "pho2phi": HistSpec("ALP_sublead_photon_phi", 40, -4.0, 4.0),
        "pho2R9": HistSpec("ALP_sublead_photon_r9", 25, 0.1, 1.0),
        "pho2IetaIeta55": HistSpec("ALP_sublead_photon_sieie", 25, 0.0, 0.08),
        "pho2ECALIso": HistSpec("ALP_sublead_photon_ecalPFClusterIso", 50, 0.0, 50.0),
        "pho2CIso": HistSpec("ALP_sublead_photon_chiso", 35, 0.0, 0.35),
        "pho2HCALIso": HistSpec("ALP_sublead_photon_hcalPFClusterIso", 25, 0.0, 4.0),
        "pho2HOE": HistSpec("ALP_sublead_photon_hoe_PUcorr", 25, 0.0, 0.1),
        "Z_m": HistSpec("Z_mass", 80, 50.0, 130.0),
        "H_m": HistSpec(mass_branch, 85, 95.0, 180.0),
        "H_pt": HistSpec("H_pt", 80, 0.0, 160.0),
        "ALP_m": HistSpec("ALP_m", 40, 0.0, 40.0),
        "var_dR_g1g2": HistSpec("var_dR_g1g2", 25, 0.0, 5.0),
        "var_PtaOverMa": HistSpec("var_PtaOverMa", 25, 0.0, 100.0),
        "var_dR_Za": HistSpec("var_dR_Za", 35, 0.0, 7.0),
        "var_dR_g1Z": HistSpec("var_dR_g1Z", 35, 0.0, 7.0),
        "var_PtaOverMh": HistSpec("var_PtaOverMh", 25, 0.0, 0.75),
        "var_Pta": HistSpec("var_Pta", 60, 0.0, 60.0),
        "var_MhMa": HistSpec("var_MhMa", 100, 100.0, 200.0),
        "var_MhMZ": HistSpec("var_MhMZ", 130, 180.0, 310.0),
        "ALP_calculatedPhotonIso": HistSpec("ALP_calculatedPhotonIso", 25, 0.0, 125.0),
        # The original script used a random ALP mass for background/data in param.
        # For a deterministic fast TTreeFormula version, background/data cycle
        # through the same mass hypotheses by Entry$.
        "param": HistSpec(f"(ALP_m - {param_cycle_expr()}) / {mass_branch}", 25, -0.3, 0.6),
    }

    specs = {var: base_specs[var] for var in var_names if var in base_specs}
    for mass in target_masses:
        mva_expr = mva_branch_map.get("__default__", {}).get(mass, f"MVA_Score_mA_{mass}")
        mass_windows = {
            "": "1",
            "1sigma": "__SIGMA_1__",
            "1P5sigma": "__SIGMA_1P5__",
            "2sigma": "__SIGMA_2__",
            "3sigma": "__SIGMA_3__",
        }
        for suffix, placeholder in mass_windows.items():
            prefix = f"mvaVal_{suffix + '_' if suffix else ''}{mass}"
            specs[prefix] = HistSpec(mva_expr, 240, -0.1, 1.1, placeholder)
            prefix_larger = f"mvaVal_larger_{suffix + '_' if suffix else ''}{mass}"
            specs[prefix_larger] = HistSpec(mva_expr, 20, 0.0, 1.0, placeholder)
    return specs


def sigma_cut(mass_branch: str, sigma_low: dict, sigma_hig: dict, mass: str, scale: float) -> str:
    low = 125.0 + sigma_low[mass] * scale
    high = 125.0 + sigma_hig[mass] * scale
    return f"({mass_branch} > {low:.8g} && {mass_branch} < {high:.8g})"


def materialize_extra_cut(extra_cut: str, var_name: str, mass_branch: str, sigma_low: dict, sigma_hig: dict) -> str:
    mass = var_name.split("_")[-1]
    if extra_cut == "__SIGMA_1__":
        return sigma_cut(mass_branch, sigma_low, sigma_hig, mass, 1.0)
    if extra_cut == "__SIGMA_1P5__":
        return sigma_cut(mass_branch, sigma_low, sigma_hig, mass, 1.5)
    if extra_cut == "__SIGMA_2__":
        return sigma_cut(mass_branch, sigma_low, sigma_hig, mass, 2.0)
    if extra_cut == "__SIGMA_3__":
        return sigma_cut(mass_branch, sigma_low, sigma_hig, mass, 3.0)
    return extra_cut


def make_hist(name: str, spec: HistSpec) -> ROOT.TH1F:
    old = ROOT.gROOT.FindObject(name)
    if old:
        old.Delete()
    hist = ROOT.TH1F(name, name, spec.nbins, spec.xmin, spec.xmax)
    hist.Sumw2()
    return hist


def add_overflow(hist: ROOT.TH1) -> None:
    nbins = hist.GetNbinsX()
    first = hist.GetBinContent(1) + hist.GetBinContent(0)
    first_err = (hist.GetBinError(1) ** 2 + hist.GetBinError(0) ** 2) ** 0.5
    last = hist.GetBinContent(nbins) + hist.GetBinContent(nbins + 1)
    last_err = (hist.GetBinError(nbins) ** 2 + hist.GetBinError(nbins + 1) ** 2) ** 0.5
    hist.SetBinContent(1, first)
    hist.SetBinError(1, first_err)
    hist.SetBinContent(nbins, last)
    hist.SetBinError(nbins, last_err)
    hist.SetBinContent(0, 0.0)
    hist.SetBinError(0, 0.0)
    hist.SetBinContent(nbins + 1, 0.0)
    hist.SetBinError(nbins + 1, 0.0)


def draw_into_hist(chain: ROOT.TChain, hist: ROOT.TH1, expr: str, weight: str, cut: str) -> None:
    draw_cut = f"({weight}) * ({cut})"
    chain.Draw(f"{expr}>>{hist.GetName()}", draw_cut, "goff")
    hist.SetDirectory(0)
    add_overflow(hist)


def branch_requirements(expr: str) -> List[str]:
    # This is intentionally conservative: it only catches the direct branch
    # expressions used here.  TTreeFormula will report complex expression issues.
    cleaned = expr.strip()
    if cleaned.replace("_", "").replace("mA", "").isalnum():
        return [cleaned]
    return []


def draw_sample_hist(
    chain: ROOT.TChain,
    sample: str,
    var_name: str,
    spec: HistSpec,
    base_weight: str,
    base_cut: str,
    mass_branch: str,
    sigma_low: dict,
    sigma_hig: dict,
    blind: bool,
    mva_branch_map: Dict[str, Dict[str, str]],
) -> ROOT.TH1F:
    sample_spec = spec
    if var_name == "param" and sample in MASS_VALUES:
        sample_spec = HistSpec(f"(ALP_m - {MASS_VALUES[sample]:.8g}) / {mass_branch}", spec.nbins, spec.xmin, spec.xmax)

    expr = sample_spec.expr
    if var_name.startswith("mvaVal"):
        mass = var_name.split("_")[-1]
        expr = mva_branch_map[sample][mass]

    extra_cut = materialize_extra_cut(sample_spec.extra_cut, var_name, mass_branch, sigma_low, sigma_hig)
    cut = f"({base_cut}) && ({extra_cut})"
    if blind and sample == "Data":
        if var_name == "H_m":
            cut = f"({cut}) && !({mass_branch} > 115.0 && {mass_branch} < 135.0)"
        if var_name.startswith("mvaVal"):
            cut = f"({cut}) && !(({expr}) > 0.95)"

    hist = make_hist(f"{var_name}_{sample}", sample_spec)
    draw_weight = "1.0" if sample == "Data" else base_weight
    draw_into_hist(chain, hist, expr, draw_weight, cut)
    return hist


def draw_sys_hist(
    chain: ROOT.TChain,
    sample: str,
    var_name: str,
    spec: HistSpec,
    sys_name: str,
    base_weight: str,
    base_cut: str,
    mass_branch: str,
    sigma_low: dict,
    sigma_hig: dict,
    mva_branch_map: Dict[str, Dict[str, str]],
) -> ROOT.TH1F:
    sample_spec = spec
    if var_name == "param" and sample in MASS_VALUES:
        sample_spec = HistSpec(f"(ALP_m - {MASS_VALUES[sample]:.8g}) / {mass_branch}", spec.nbins, spec.xmin, spec.xmax)

    expr = sample_spec.expr
    if var_name.startswith("mvaVal"):
        mass = var_name.split("_")[-1]
        expr = mva_branch_map[sample][mass]

    hist = make_hist(f"{var_name}_{sample}_{sys_name}", sample_spec)
    extra_cut = materialize_extra_cut(sample_spec.extra_cut, var_name, mass_branch, sigma_low, sigma_hig)
    cut = f"({base_cut}) && ({extra_cut})"

    if sample == "Data":
        draw_into_hist(chain, hist, expr, "1.0", cut)
        return hist

    central = sys_central_branch(sys_name)
    if has_branch(chain, sys_name) and has_branch(chain, central):
        sys_weight = safe_div_weight(base_weight, sys_name, central)
    else:
        print(f"[Sys] {sample}/{sys_name}: missing branch, using central weight")
        sys_weight = base_weight
    draw_into_hist(chain, hist, expr, sys_weight, cut)
    return hist


def sideband_scale_bkg_to_data(histos, histos_sys, analyzer_cfg, signal_low=115.0, signal_high=135.0):
    data_hist = histos.get("H_m", {}).get("Data")
    if not data_hist:
        return 1.0

    axis = data_hist.GetXaxis()
    nbins = axis.GetNbins()
    low_bin = axis.FindFixBin(signal_low)
    high_bin = axis.FindFixBin(signal_high)

    def sideband_integral(hist):
        left = hist.Integral(1, low_bin - 1) if low_bin > 1 else 0.0
        right = hist.Integral(high_bin + 1, nbins) if high_bin < nbins else 0.0
        return left + right

    data_sb = sideband_integral(data_hist)
    bkg_sb = sum(sideband_integral(histos["H_m"][sample]) for sample in analyzer_cfg.bkg_names)
    if data_sb <= 0 or bkg_sb <= 0:
        return 1.0

    scale = data_sb / bkg_sb
    for sample in analyzer_cfg.bkg_names:
        histos["H_m"][sample].Scale(scale)
        if histos_sys and "H_m" in histos_sys:
            for sys_name in analyzer_cfg.sys_names:
                histos_sys["H_m"][sample][sys_name].Scale(scale)
    print(f"[SideBandScale] H_m Data={data_sb:.3f} Bkg={bkg_sb:.3f} Scale={scale:.4f}")
    return scale


def clone_sys_from_central(histos, analyzer_cfg):
    histos_sys = {}
    for var_name, sample_hists in histos.items():
        histos_sys[var_name] = {}
        for sample, hist in sample_hists.items():
            histos_sys[var_name][sample] = {}
            for sys_name in analyzer_cfg.sys_names:
                clone = hist.Clone(f"{var_name}_{sample}_{sys_name}")
                clone.SetDirectory(0)
                histos_sys[var_name][sample][sys_name] = clone
    return histos_sys


def main() -> None:
    args = parse_args()
    start_time = time.time()

    analyzer_cfg = AC.Analyzer_Config("inclusive", args.year, args.region, args.mva)
    analyzer_cfg.mva = bool(args.mva)
    analyzer_cfg.mva_alp_mass = str(args.mA) if args.mva else "M1"
    if args.input_dir:
        analyzer_cfg.sample_loc = args.input_dir

    if args.cut:
        analyzer_cfg.sig_names = [args.mA]
        analyzer_cfg.samp_names = analyzer_cfg.bkg_names + analyzer_cfg.sig_names + ["Data"]

    target_masses = list(analyzer_cfg.sig_names) if args.mva and args.year == "run3" else ([args.mA] if args.mva else [])
    plot_output_path = args.out_dir or os.path.join(PLOT_DIR, "plots", f"fast_plot_{analyzer_cfg.out_region_name}")
    os.makedirs(plot_output_path, exist_ok=True)
    root_output_path = os.path.join(PLOT_DIR, "plots", f"ALP_fast_plot_{args.year}_{analyzer_cfg.out_region_name}.root")
    os.makedirs(os.path.dirname(root_output_path), exist_ok=True)

    print(f"[Config] input={analyzer_cfg.sample_loc}")
    print(f"[Config] output={plot_output_path}")
    analyzer_cfg.Print_Config()

    ntuples = LoadNtuples(analyzer_cfg)
    mass_branch = common_branch(ntuples, analyzer_cfg.samp_names, ("H_m", "H_mass"))
    base_weight = common_branch(ntuples, [s for s in analyzer_cfg.samp_names if s != "Data"], ("weight", "weight_central"))
    base_cut = f"({region_cut(args.region, mass_branch)}) && ({channel_cut(args)})"
    print(f"[Config] mass_branch={mass_branch}, weight_branch={base_weight}")

    mva_branch_map = {"__default__": {}}
    for sample in analyzer_cfg.samp_names:
        mva_branch_map[sample] = {}
        for mass in target_masses:
            mva_branch_map[sample][mass] = resolve_mva_branch(ntuples[sample], mass)
            mva_branch_map["__default__"][mass] = mva_branch_map[sample][mass]
            print(f"[MVA] sample={sample:>12s} mass={mass:>3s} branch={mva_branch_map[sample][mass]}")

    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)
    var_names = list(plot_cfg.var_title_map.keys())
    if args.mva:
        mva_vars = []
        for mass in target_masses:
            var_names.append(f"mvaVal_{mass}")
            var_names.append(f"mvaVal_larger_{mass}")
            for region in ("1sigma", "1P5sigma", "2sigma", "3sigma"):
                mva_vars.append(f"mvaVal_{region}_{mass}")
                mva_vars.append(f"mvaVal_larger_{region}_{mass}")
        var_names.extend(mva_vars)

    sigma_low, sigma_hig = getMassSigma(analyzer_cfg)
    hist_specs = build_hist_specs(var_names, target_masses, mva_branch_map, mass_branch)
    var_names = [var for var in var_names if var in hist_specs]

    histos = {var_name: {} for var_name in var_names}
    for sample in analyzer_cfg.samp_names:
        print(f"\n[Sample] {sample}: entries={ntuples[sample].GetEntries()}")
        for var_name in var_names:
            spec = hist_specs[var_name]
            requirements_expr = spec.expr
            if var_name.startswith("mvaVal"):
                requirements_expr = mva_branch_map[sample][var_name.split("_")[-1]]
            missing = branch_requirements(requirements_expr)
            if any(not has_branch(ntuples[sample], branch) for branch in missing):
                print(f"[Skip] {sample}/{var_name}: missing {missing}")
                histos[var_name][sample] = make_hist(f"{var_name}_{sample}", spec)
                continue
            histos[var_name][sample] = draw_sample_hist(
                ntuples[sample],
                sample,
                var_name,
                spec,
                base_weight,
                base_cut,
                mass_branch,
                sigma_low,
                sigma_hig,
                args.blind,
                mva_branch_map,
            )
            plot_cfg.SetHistStyles(histos[var_name][sample], sample)

    if args.skip_sys:
        histos_sys = clone_sys_from_central(histos, analyzer_cfg)
    else:
        histos_sys = {var_name: {sample: {} for sample in analyzer_cfg.samp_names} for var_name in var_names}
        for var_name in var_names:
            print(f"\n[Sys] {var_name}")
            for sample in analyzer_cfg.samp_names:
                for sys_name in analyzer_cfg.sys_names:
                    histos_sys[var_name][sample][sys_name] = draw_sys_hist(
                        ntuples[sample],
                        sample,
                        var_name,
                        hist_specs[var_name],
                        sys_name,
                        base_weight,
                        base_cut,
                        mass_branch,
                        sigma_low,
                        sigma_hig,
                        mva_branch_map,
                    )
                    plot_cfg.SetHistStyles(histos_sys[var_name][sample][sys_name], sample)

    out_file = ROOT.TFile(root_output_path, "RECREATE")
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

    out_file.cd()
    scaled_sig = {}

    for var_name in var_names:
        if var_name == "H_m":
            scale_factor = sideband_scale_bkg_to_data(histos, histos_sys, analyzer_cfg)
        else:
            scale_factor = ScaleBkgToData(histos[var_name], analyzer_cfg, histos_sys.get(var_name))
        if scale_factor != 1.0:
            print(f"[Info] Applied background scaling ({var_name}) = {scale_factor:.4f}")

        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        if not stacks["all"].GetStack() or stacks["all"].GetStack().GetEntries() == 0:
            print(f"[Warning] Stack for {var_name} is empty. Skipping drawing.")
            continue

        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]["Data"], stacks["all"].GetStack().Last(), var_name)
        total_unc = Total_Unc(stacks["bkg"], histos_sys[var_name], analyzer_cfg)

        if args.ln:
            canvas = CreateCanvas(var_name + "_log")
            draw_fast_on_canvas(
                canvas,
                var_name,
                plot_cfg,
                stacks,
                histos[var_name],
                scaled_sig,
                ratio_plot,
                total_unc,
                args.cut,
                args.mA,
                log_y=True,
            )
            canvas.Write()
            SaveCanvPic(canvas, plot_output_path, output_name(var_name) + "_log")
        else:
            canvas = CreateCanvas(var_name)
            draw_fast_on_canvas(
                canvas,
                var_name,
                plot_cfg,
                stacks,
                histos[var_name],
                scaled_sig,
                ratio_plot,
                total_unc,
                args.cut,
                args.mA,
                log_y=False,
            )
            canvas.Write()
            SaveCanvPic(canvas, plot_output_path, output_name(var_name))

    out_file.Close()
    elapsed = time.time() - start_time
    print(f"Total runtime: {time.strftime('%H:%M:%S', time.gmtime(elapsed))}")
    print("Done")


if __name__ == "__main__":
    main()
