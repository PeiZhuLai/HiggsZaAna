####################################################
####################################################

import argparse
import os
import sys
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir))
PROJECT_DIR = os.path.abspath(os.path.join(PLOT_DIR, os.pardir))

sys.path.insert(0, os.path.join(PLOT_DIR, "lib"))

from ROOT import *  # noqa: F401,F403
import ROOT

from Plot_Helper import (
    MakeStack,
    CreateCanvas,
    DrawOnCanv,
    SaveCanvPic,
    MakeLumiLabel,
    MakeCMSDASLabel,
    ScaleSignal,
    MakeRatioPlot,
    MakeLegend,
    Total_Unc,
    ScaleBkgToData,
)
import Analyzer_Configs as AC
import Plot_Configs as PC
import tdrstyle


parser = argparse.ArgumentParser(description="Draw data/MC plots from plot_variable_dataVmc.py ROOT histograms")
parser.add_argument("-y", "--Year", dest="year", default="run3", help="which year's datasets")
parser.add_argument("--region", dest="region", type=int, default=0, help="0 for full region, 1 for signal region, 2 for sideband region")
parser.add_argument("-m", "--mva", dest="mva", action="store_true", default=False, help="use mva variables")
parser.add_argument("--cut", dest="cut", action="store_true", default=False, help="apply mva cut style")
parser.add_argument("--mA", dest="mA", default="M5", help="ALP mass used for cut-style drawing")
parser.add_argument("-b", "--blind", dest="blind", action="store_true", default=False, help="kept for command compatibility; histograms are already blinded if needed")
parser.add_argument("-ln", "--ln", dest="ln", action="store_true", default=False, help="draw log-y plots")
parser.add_argument("--inputTag", dest="input_tag", default=None, help="tag appended to the input ROOT file name")
parser.add_argument("--outputTag", dest="output_tag", default=None, help="tag appended to the output plot directory; defaults to inputTag")
parser.add_argument("--inputRoot", dest="input_root", default=None, help="explicit input ROOT file path")
parser.add_argument("--outputRoot", dest="output_root", default=None, help="explicit output ROOT file path for canvases")
parser.add_argument("--plotDir", dest="plot_dir", default=None, help="explicit output directory for PDFs")
parser.add_argument("--noCdf", dest="draw_cdf", action="store_false", default=True, help="do not draw CDF-transformed larger-MVA plots")
args = parser.parse_args()


SIDEBAND_REWEIGHT_UNC_SYS_NAMES = ("sideband_rwgt_reweight_up", "sideband_rwgt_reweight_down")

pdfName_map = {
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
    "pho_pt_asym": "15_pho_pt_asym",
    "param": "16_param",
    "ALP_m": "17_ALP_m",
}

MVA_LARGER_NBINS = 20
MVA_LARGER_XMIN = 0.0
MVA_LARGER_XMAX = 1.0
SIGMA_REGION_SCALES = {
    "1sigma": 1.0,
    "1P5sigma": 1.5,
    "2sigma": 2.0,
    "3sigma": 3.0,
}


def _clean_tag(tag):
    if tag is None:
        return None
    tag = os.path.basename(str(tag).strip().strip("/"))
    return tag or None


def _add_tag(root_name, tag):
    tag = _clean_tag(tag)
    if not tag:
        return root_name
    root_base, root_ext = os.path.splitext(root_name)
    return f"{root_base}_{tag}{root_ext or '.root'}"


def _append_suffix(root_path, suffix):
    root_dir = os.path.dirname(root_path)
    root_name = os.path.basename(root_path)
    root_base, root_ext = os.path.splitext(root_name)
    return os.path.join(root_dir, f"{root_base}_{suffix}{root_ext or '.root'}")


def _pdf_output_name(var_name):
    return pdfName_map.get(var_name, var_name)


def _is_full_range_mva_larger_var(var_name, target_masses):
    return var_name in [f"mvaVal_larger_{mass}" for mass in target_masses]


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


def SideBandScaleBkgToData(histos, histos_sys, analyzer_cfg, signal_low=115.0, signal_high=135.0):
    if "H_m" not in histos:
        print("[SideBandScale] H_m not found, skip.")
        return 1.0

    data_key = None
    for key in histos["H_m"].keys():
        if key.lower() == "data":
            data_key = key
            break
    if data_key is None:
        print("[SideBandScale] Data hist not found, skip.")
        return 1.0

    axis = histos["H_m"][data_key].GetXaxis()
    nbins = axis.GetNbins()
    sig_low_bin = axis.FindFixBin(signal_low)
    sig_high_bin = axis.FindFixBin(signal_high)

    def sideband_integral(hist):
        left = hist.Integral(1, sig_low_bin - 1) if sig_low_bin > 1 else 0.0
        right = hist.Integral(sig_high_bin + 1, nbins) if sig_high_bin < nbins else 0.0
        return left + right

    data_sb = sideband_integral(histos["H_m"][data_key])
    bkg_samples = [
        sample for sample in analyzer_cfg.samp_names
        if (sample.lower() != "data") and (sample not in analyzer_cfg.sig_names)
    ]
    bkg_sb = sum(sideband_integral(histos["H_m"][sample]) for sample in bkg_samples)

    if bkg_sb <= 0:
        print("[SideBandScale] Background sideband yield = 0, skip scaling.")
        return 1.0

    scale_factor = data_sb / bkg_sb
    for sample in bkg_samples:
        histos["H_m"][sample].Scale(scale_factor)
        if histos_sys and "H_m" in histos_sys:
            for sys_name in analyzer_cfg.sys_names:
                histos_sys["H_m"][sample][sys_name].Scale(scale_factor)

    print(f"[SideBandScale] H_m sideband Data={data_sb:.3f} Bkg(before)={bkg_sb:.3f} Scale={scale_factor:.4f}")
    return scale_factor


def _build_var_names(plot_cfg, analyzer_cfg, target_masses):
    var_names = list(plot_cfg.var_title_map.keys())
    if analyzer_cfg.mva:
        var_mva = []
        for alp_mass in target_masses:
            var_names.append("mvaVal_" + alp_mass)
            var_names.append("mvaVal_larger_" + alp_mass)
            for region in ["1sigma", "1P5sigma", "2sigma", "3sigma"]:
                var_mva.append("mvaVal_" + region + "_" + alp_mass)
                var_mva.append("mvaVal_larger_" + region + "_" + alp_mass)
        var_names = var_names + var_mva
    return var_names


def _open_root(path):
    root_file = TFile.Open(path, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Cannot open input ROOT file: %s" % path)
    return root_file


def _clone_hist(directory, name):
    obj = directory.Get(name) if directory else None
    if not obj:
        raise RuntimeError("Missing histogram '%s' in directory '%s'" % (name, directory.GetName() if directory else "None"))
    hist = obj.Clone(name)
    hist.SetDirectory(0)
    return hist


def _sys_key_names(root_file):
    sys_dir = root_file.Get("sys_dir")
    if not sys_dir:
        return []
    return [key.GetName() for key in sys_dir.GetListOfKeys()]


def _append_available_extra_systematics(root_file, analyzer_cfg):
    key_names = _sys_key_names(root_file)
    for sys_name in SIDEBAND_REWEIGHT_UNC_SYS_NAMES:
        if sys_name in analyzer_cfg.sys_names:
            continue
        if any(name.endswith("_" + sys_name) for name in key_names):
            analyzer_cfg.sys_names.append(sys_name)


def _load_histo_maps(root_file, var_names, analyzer_cfg, plot_cfg):
    raw_dir = root_file.Get("raw_plots")
    sys_dir = root_file.Get("sys_dir")
    if not raw_dir:
        raise RuntimeError("Input ROOT file has no raw_plots directory.")
    if not sys_dir:
        raise RuntimeError("Input ROOT file has no sys_dir directory.")

    histos = {}
    histos_sys = {}
    missing_sideband_sys = []
    for var_name in var_names:
        histos[var_name] = {}
        histos_sys[var_name] = {}
        for sample in analyzer_cfg.samp_names:
            hist = _clone_hist(raw_dir, f"{var_name}_{sample}")
            plot_cfg.SetHistStyles(hist, sample)
            histos[var_name][sample] = hist

            histos_sys[var_name][sample] = {}
            for sys_name in analyzer_cfg.sys_names:
                hist_sys_name = f"{var_name}_{sample}_{sys_name}"
                hist_sys_obj = sys_dir.Get(hist_sys_name)
                if hist_sys_obj:
                    hist_sys = hist_sys_obj.Clone(hist_sys_name)
                    hist_sys.SetDirectory(0)
                elif sys_name in SIDEBAND_REWEIGHT_UNC_SYS_NAMES:
                    hist_sys = hist.Clone(hist_sys_name)
                    hist_sys.SetDirectory(0)
                    missing_sideband_sys.append(hist_sys_name)
                else:
                    hist_sys = _clone_hist(sys_dir, hist_sys_name)
                plot_cfg.SetHistStyles(hist_sys, sample)
                histos_sys[var_name][sample][sys_name] = hist_sys

    if missing_sideband_sys:
        print(
            "[Warning] Missing {} sideband reweight uncertainty histograms; "
            "using nominal shapes for those entries.".format(len(missing_sideband_sys))
        )

    return histos, histos_sys


def _draw_all(histos, histos_sys, var_names, target_masses, analyzer_cfg, plot_cfg, output_root):
    lumi_label = MakeLumiLabel(plot_cfg.lumi)
    cms_label = MakeCMSDASLabel()

    for var_name in var_names:
        if var_name == "H_m":
            scale_factor = SideBandScaleBkgToData(histos, histos_sys, analyzer_cfg, signal_low=115.0, signal_high=135.0)
        else:
            scale_factor = ScaleBkgToData(histos[var_name], analyzer_cfg, histos_sys.get(var_name))
        if scale_factor != 1.0:
            print(f"[Info] Applied background scaling ({var_name}) = {scale_factor:.4f}")

        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        if stacks["all"].GetStack().GetEntries() == 0:
            stack_entry = stacks["all"].GetStack().GetEntries()
            print(f"stack_entry: {stack_entry}")
            print(f"[Warning] Stack for {var_name} is empty. Skipping drawing.")
            continue

        scaled_sig = {}
        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]["Data"], stacks["all"].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)
        total_unc = Total_Unc(stacks["bkg"], histos_sys[var_name], analyzer_cfg)

        if args.ln:
            canv = CreateCanvas(var_name + "_log")
            DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, args.cut, args.mA, logY=True)
            output_root.cd()
            canv.Write()
            SaveCanvPic(canv, analyzer_cfg.plot_output_path, _pdf_output_name(var_name) + "_log")
        else:
            canv = CreateCanvas(var_name)
            DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, args.cut, args.mA, logY=False)
            output_root.cd()
            canv.Write()
            SaveCanvPic(canv, analyzer_cfg.plot_output_path, _pdf_output_name(var_name))

        if args.draw_cdf and analyzer_cfg.mva and _is_full_range_mva_larger_var(var_name, target_masses):
            histos_cdf, histos_sys_cdf = _build_cdf_histo_maps(var_name, histos[var_name], histos_sys[var_name], analyzer_cfg)
            cdf_var_name = var_name + "_cdf"
            stacks_cdf = MakeStack(histos_cdf, analyzer_cfg, cdf_var_name)

            if stacks_cdf["all"].GetStack().GetEntries() == 0:
                stack_entry = stacks_cdf["all"].GetStack().GetEntries()
                print(f"stack_entry: {stack_entry}")
                print(f"[Warning] CDF stack for {var_name} is empty. Skipping drawing.")
                continue

            scaled_sig_cdf = {}
            for sample in analyzer_cfg.sig_names:
                scaled_sig_cdf[sample] = ScaleSignal(plot_cfg, stacks_cdf[sample], histos_cdf[sample], cdf_var_name)

            ratio_plot_cdf = MakeRatioPlot(histos_cdf["Data"], stacks_cdf["all"].GetStack().Last(), cdf_var_name)
            legend_cdf = MakeLegend(plot_cfg, histos_cdf, scaled_sig_cdf)
            total_unc_cdf = Total_Unc(stacks_cdf["bkg"], histos_sys_cdf, analyzer_cfg)
            cdf_y_axis_title = "Entries / %.3f" % histos_cdf["Data"].GetBinWidth(1)

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
                output_root.cd()
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
                output_root.cd()
                canv_cdf.Write()
                SaveCanvPic(canv_cdf, analyzer_cfg.plot_output_path, _pdf_output_name(var_name) + "_cdf")


def main():
    start_time = time.time()

    gROOT.SetBatch(True)
    tdrstyle.setTDRStyle()

    analyzer_cfg = AC.Analyzer_Config("inclusive", args.year, args.region, args.mva)
    analyzer_cfg.mva = bool(args.mva)
    analyzer_cfg.mva_alp_mass = str(args.mA) if args.mva else "M1"
    all_sig_names = list(analyzer_cfg.sig_names)

    if args.cut:
        analyzer_cfg.sig_names = [args.mA]
        analyzer_cfg.samp_names = analyzer_cfg.bkg_names + analyzer_cfg.sig_names + ["Data"]

    target_masses = list(all_sig_names) if args.mva and analyzer_cfg.year in ["run3", "run3_NFlow"] else ([analyzer_cfg.mva_alp_mass] if args.mva else [])

    base_root_name = analyzer_cfg.root_output_name
    input_tag = _clean_tag(args.input_tag)
    input_root_path = args.input_root or os.path.join(analyzer_cfg.out_dir, _add_tag(base_root_name, input_tag))

    plot_tag = _clean_tag(args.output_tag) or input_tag
    if args.plot_dir:
        analyzer_cfg.plot_output_path = args.plot_dir
    elif plot_tag:
        analyzer_cfg.plot_output_path = os.path.join(analyzer_cfg.plot_output_path, plot_tag)
    os.makedirs(analyzer_cfg.plot_output_path, exist_ok=True)

    output_root_path = args.output_root or _append_suffix(input_root_path, "plots")
    os.makedirs(os.path.dirname(output_root_path) or ".", exist_ok=True)

    print("[Input] ROOT:", input_root_path)
    print("[Output] plot dir:", analyzer_cfg.plot_output_path)
    print("[Output] canvas ROOT:", output_root_path)

    input_file = _open_root(input_root_path)
    _append_available_extra_systematics(input_file, analyzer_cfg)

    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)
    var_names = _build_var_names(plot_cfg, analyzer_cfg, target_masses)
    histos, histos_sys = _load_histo_maps(input_file, var_names, analyzer_cfg, plot_cfg)
    input_file.Close()

    output_root = TFile(output_root_path, "RECREATE")
    _draw_all(histos, histos_sys, var_names, target_masses, analyzer_cfg, plot_cfg, output_root)
    output_root.Close()

    elapsed = time.time() - start_time
    print(f"Total runtime: {time.strftime('%H:%M:%S', time.gmtime(elapsed))}")
    print("Done")


main()
