import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

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
LOCAL_MVA_CUT_CANDIDATES = [
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/HiggsZaAna/Plot/output/MVAcut_points_run3.json",
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/optimize_alias/MVAcut_points_run3.json",
    PLOT_DIR.parent.parent / "AN-25-172/figure_ALP/run3/selection_alias/MVAcut_points_run3.json",
]
LOCAL_OUTPUT_DIR = PLOT_DIR / "plots" / "bkgmcScupltingCheck"

TARGET_MASSES = list(range(1, 31))
SIGNAL_MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
BLIND_LOW = 115.0
BLIND_HIGH = 135.0
H_M_NBINS = 85
H_M_XMIN = 95.0
H_M_XMAX = 180.0
BKG_LABELS = {"DYJetsToLL": "Z + jets", "DYGto2LG": "Z + #gamma"}


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
            hist.SetMarkerSize(1.0)
            hist.SetLineColor(ROOT.kBlack)
            hist.SetLineWidth(2)
        elif sample in analyzer_cfg.bkg_names:
            hist.SetFillColor(plot_cfg.colors[sample])
            hist.SetLineColor(plot_cfg.colors[sample])
            hist.SetLineWidth(1)
        elif sample in analyzer_cfg.sig_names:
            hist.SetLineColor(plot_cfg.colors[sample])
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
    analyzer_cfg: Analyzer_Config,
    plot_cfg: Plot_Config,
    output_dir: Path,
    logy: bool = False,
) -> None:
    _style_histograms(histos, plot_cfg, analyzer_cfg)

    scale = _scale_bkg_to_data_sideband(histos, analyzer_cfg)
    if abs(scale - 1.0) > 1e-12:
        print(f"[mA={mass:02d}] Applied sideband scale factor = {scale:.4f}")

    bkg_total = _clone_total_bkg(histos, analyzer_cfg, f"h_bkg_total_mA_{mass}")
    if bkg_total is None or bkg_total.Integral() <= 0.0:
        print(f"[mA={mass:02d}] Background is empty after MVA cut. Skip plot.")
        return

    stack = _make_stack(histos, analyzer_cfg, f"mH_mA_{mass}")
    ratio = MakeRatioPlot(histos["Data"], bkg_total, f"H_m_mA_{mass}")
    stat_abs, stat_norm = Get_StatUnc(bkg_total)

    signal_sample = _signal_sample_name(mass)
    signal_hist = None
    if signal_sample and signal_sample in histos and histos[signal_sample].Integral() > 0.0:
        signal_hist = histos[signal_sample].Clone(f"h_sig_draw_{signal_sample}_mA_{mass}")
        signal_hist.SetDirectory(0)
        signal_hist.SetLineColor(plot_cfg.colors[signal_sample])
        signal_hist.SetLineWidth(4)
        signal_hist.SetFillStyle(0)
        signal_hist.SetFillColor(0)

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

    ymax = max(histos["Data"].GetMaximum(), bkg_total.GetMaximum(), signal_hist.GetMaximum() if signal_hist else 0.0)
    if ymax <= 0.0:
        ymax = 1.0

    histos["Data"].SetMinimum(1e-2 if logy else 0.0)
    histos["Data"].SetMaximum(ymax * (1.0e4 if logy else 1.4))
    if logy:
        upper_pad.SetLogy()

    histos["Data"].Draw("PE")
    histos["Data"].GetXaxis().SetLabelSize(0.0)
    histos["Data"].GetXaxis().SetTitleOffset(0.95)
    histos["Data"].GetYaxis().SetTitle(f"Events / {histos['Data'].GetBinWidth(1):.2f} GeV")
    histos["Data"].GetYaxis().SetTitleSize(0.07)
    histos["Data"].GetYaxis().SetLabelSize(0.055)
    histos["Data"].GetYaxis().SetTitleFont(42)
    histos["Data"].GetYaxis().SetTitleOffset(1.15)

    stack.SetMinimum(1e-2 if logy else 0.0)
    stack.SetMaximum(ymax * (1.0e4 if logy else 1.4))
    stack.Draw("HIST SAME")

    Draw_unc(stat_abs, TColor.GetColor("#404040"), alpha=0.90, fill_style=3354)

    if signal_hist:
        signal_hist.Draw("HIST SAME")

    histos["Data"].Draw("PE SAME")
    histos["Data"].Draw("AXIS SAME")

    legend = TLegend(0.58, 0.60, 0.95, 0.87)
    ROOT.SetOwnership(legend, False)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetFillColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.AddEntry(histos["Data"], "Data", "PE")
    for sample in analyzer_cfg.bkg_names:
        legend.AddEntry(histos[sample], BKG_LABELS.get(sample, sample), "f")
    legend.AddEntry(stat_abs, "Stat. Unc.", "f")
    if signal_hist:
        legend.AddEntry(signal_hist, f"m_{{a}} = {mass} GeV", "l")
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

    ratio.GetXaxis().SetTitle(plot_cfg.var_title_map["H_m"])
    ratio.GetXaxis().SetTitleOffset(1.0)
    ratio.Draw("APZ SAME")
    Draw_unc(stat_norm, TColor.GetColor("#404040"), alpha=0.90, fill_style=3354)
    ratio.Draw("PZ SAME")

    save_name = f"bkgmcScupltingCheck_mA{mass:02d}"
    if logy:
        save_name += "_log"
    SaveCanvPic(canvas, str(output_dir), save_name)


def _book_histograms(analyzer_cfg: Analyzer_Config) -> Dict[int, Dict[str, TH1F]]:
    histos: Dict[int, Dict[str, TH1F]] = {}
    for mass in TARGET_MASSES:
        histos[mass] = {}
        for sample in analyzer_cfg.samp_names:
            hist = TH1F(
                f"H_m_mA_{mass}_{sample}",
                f"H_m_mA_{mass}_{sample}",
                H_M_NBINS,
                H_M_XMIN,
                H_M_XMAX,
            )
            hist.Sumw2()
            hist.SetDirectory(0)
            histos[mass][sample] = hist
    return histos


def _fill_histograms(
    ntuples,
    analyzer_cfg: Analyzer_Config,
    mva_cuts: Dict[int, float],
    histos: Dict[int, Dict[str, TH1F]],
    blind: bool,
    only_ele: bool,
    only_mu: bool,
) -> None:
    mva_branch_map: Dict[str, Dict[int, Optional[str]]] = {}
    for sample in analyzer_cfg.samp_names:
        masses_to_use = TARGET_MASSES if sample in analyzer_cfg.bkg_names or sample == "Data" else []
        signal_mass = _parse_ma(sample)
        if signal_mass in TARGET_MASSES:
            masses_to_use = [signal_mass]

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
                if score is None or score < cut:
                    continue

                if sample == "Data" and blind and BLIND_LOW <= h_mass <= BLIND_HIGH:
                    continue

                histos[mass][sample].Fill(h_mass, weight)


def main():
    parser = argparse.ArgumentParser(description="Plot background mH sculpting after MVA cuts for mA = 1..30.")
    parser.add_argument("-y", "--Year", dest="year", default="run3", help="Only run3 is supported.")
    parser.add_argument("--mva-cut-json", default=DEFAULT_MVA_CUT_JSON, help="Path to MVA cut JSON.")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Directory for output plots.")
    parser.add_argument("--unblind", dest="blind", action="store_false", default=True, help="Show data inside 115-135 GeV.")
    parser.add_argument("-ln", "--ln", dest="ln", action="store_true", default=False, help="Save logY plots.")
    parser.add_argument("--ele", dest="ele", action="store_true", default=False, help="Electron channel only.")
    parser.add_argument("--mu", dest="mu", action="store_true", default=False, help="Muon channel only.")
    args = parser.parse_args()

    if args.year != "run3":
        raise ValueError("plot_bkgmcScupltingCheck.py currently supports only --Year run3.")
    if args.ele and args.mu:
        raise ValueError("Choose at most one of --ele or --mu.")

    gROOT.SetBatch(True)
    tdrstyle.setTDRStyle()

    mva_cut_path = _resolve_mva_cut_json(args.mva_cut_json)
    output_dir = _resolve_output_dir(args.output_dir)
    print(f"[Input] MVA cut JSON: {mva_cut_path}")
    print(f"[Output] Plot directory: {output_dir}")

    mva_cuts = _complete_mva_cuts(_parse_mva_cuts(str(mva_cut_path)), TARGET_MASSES)

    analyzer_cfg = Analyzer_Config("inclusive", args.year, 0, True)
    analyzer_cfg.mva = True
    analyzer_cfg.plot_output_path = str(output_dir)
    plot_cfg = Plot_Config(analyzer_cfg, args.year)
    ntuples = LoadNtuples(analyzer_cfg)

    histos = _book_histograms(analyzer_cfg)
    _fill_histograms(
        ntuples=ntuples,
        analyzer_cfg=analyzer_cfg,
        mva_cuts=mva_cuts,
        histos=histos,
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
            analyzer_cfg=analyzer_cfg,
            plot_cfg=plot_cfg,
            output_dir=output_dir,
            logy=args.ln,
        )

    print("Done")


if __name__ == "__main__":
    main()
