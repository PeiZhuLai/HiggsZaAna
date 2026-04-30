#!/usr/bin/env python3
#
# Fast sideband data/MC scale-factor validation plots.
#
# This script intentionally uses TChain.Draw instead of a Python event loop.

import argparse
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import ROOT

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, os.path.join(PLOT_DIR, "lib"))
sys.path.insert(0, "%s/lib" % os.getcwd())
import Analyzer_Configs as AC


ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)


@dataclass(frozen=True)
class PlotSpec:
    name: str
    expr: str
    title: str
    nbins: int
    xmin: float
    xmax: float
    selection: str
    flavor: str


@dataclass(frozen=True)
class WeightSpec:
    name: str
    after_weight: str
    before_weight: str
    after_label: str
    before_label: str
    denominator_branches: Tuple[str, ...]
    plot_flavors: Tuple[str, ...]


def sideband_cut(mass_branch: str) -> str:
    return (
        f"({mass_branch} > 95.0 && {mass_branch} < 180.0 && "
        f"!({mass_branch} > 115.0 && {mass_branch} < 135.0))"
    )

PLOTS: Tuple[PlotSpec, ...] = (
    PlotSpec(
        "lead_electron_pt",
        "Z_lead_lepton_pt",
        "Lead electron p_{T} [GeV]",
        50,
        0.0,
        200.0,
        "z_ee == 1",
        "electron",
    ),
    PlotSpec(
        "sublead_electron_pt",
        "Z_sublead_lepton_pt",
        "Sublead electron p_{T} [GeV]",
        50,
        0.0,
        120.0,
        "z_ee == 1",
        "electron",
    ),
    PlotSpec(
        "lead_muon_pt",
        "Z_lead_lepton_pt",
        "Lead muon p_{T} [GeV]",
        50,
        0.0,
        200.0,
        "z_mumu == 1",
        "muon",
    ),
    PlotSpec(
        "sublead_muon_pt",
        "Z_sublead_lepton_pt",
        "Sublead muon p_{T} [GeV]",
        50,
        0.0,
        120.0,
        "z_mumu == 1",
        "muon",
    ),
    PlotSpec(
        "lead_photon_pt",
        "ALP_lead_photon_pt",
        "Lead photon p_{T} [GeV]",
        50,
        0.0,
        200.0,
        "1",
        "photon",
    ),
    PlotSpec(
        "sublead_photon_pt",
        "ALP_sublead_photon_pt",
        "Sublead photon p_{T} [GeV]",
        50,
        0.0,
        120.0,
        "1",
        "photon",
    ),
)


WEIGHTS: Tuple[WeightSpec, ...] = (
    WeightSpec(
        "trigger_sf",
        "weight_central",
        "weight_central / weight_hlt_sf_central",
        "Weighted After Trigger SFs",
        "Weighted Before Trigger SFs",
        ("weight_hlt_sf_central",),
        ("electron", "muon"),
    ),
    WeightSpec(
        "electron_sf",
        "weight_central",
        (
            "weight_central / "
            "(weight_electron_wplid_sf_SelectedElectron_central * "
            "weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_central)"
        ),
        "Weighted Electron SFs",
        "Weighted Before Electron SFs",
        (
            "weight_electron_wplid_sf_SelectedElectron_central",
            "weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_central",
        ),
        ("electron",),
    ),
    WeightSpec(
        "muon_sf",
        "weight_central",
        (
            "weight_central / "
            "(weight_muon_looseid_sf_SelectedMuon_central * "
            "weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_central)"
        ),
        "Weighted Muon SFs",
        "Weighted Before Muon SFs",
        (
            "weight_muon_looseid_sf_SelectedMuon_central",
            "weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_central",
        ),
        ("muon",),
    ),
    WeightSpec(
        "photon_id_sf",
        "weight_central",
        "weight_central / weight_photon_id_sf_SelectedPhoton_central",
        "Weighted Photon ID SFs",
        "Weighted Before Photon ID SFs",
        ("weight_photon_id_sf_SelectedPhoton_central",),
        ("photon",),
    ),
    WeightSpec(
        "photon_csev_sf",
        "weight_central",
        "weight_central / weight_photon_csev_sf_SelectedPhoton_central",
        "Weighted Photon CSEV SFs",
        "Weighted Before Photon CSEV SFs",
        ("weight_photon_csev_sf_SelectedPhoton_central",),
        ("photon",),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fast sideband SF validation plots")
    parser.add_argument("-y", "--Year", dest="year", default="run3")
    parser.add_argument("--input-dir", default=None, help="Override Analyzer_Configs sample_loc")
    parser.add_argument("--out-dir", default=None, help="Output directory")
    parser.add_argument("--eras", default="all", help="Comma list of eras, or 'all'")
    parser.add_argument("--normalize", action="store_true", help="Normalize each histogram to unit area")
    parser.add_argument("-ln", "--ln", "--log", dest="log_y", action="store_true", help="Use log scale")
    parser.add_argument("-m", "--mva", action="store_true", help="Accepted for compatibility; ignored")
    parser.add_argument("-b", "--blind", action="store_true", help="Accepted for compatibility; ignored")
    return parser.parse_args()


def set_style() -> None:
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetTitleFillColor(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetPadColor(ROOT.kWhite)
    ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
    ROOT.gStyle.SetLegendBorderSize(0)


def add_file(chain: ROOT.TChain, path: str) -> bool:
    if os.path.exists(path):
        chain.Add(path)
        return True
    print(f"[MISS] {path}")
    return False


def dyg_subsamples_for_era(cfg: AC.Analyzer_Config, era: str) -> Sequence[str]:
    if era.startswith("2022"):
        return cfg.bkg_2022
    return cfg.bkg_2023


def build_chains(
    cfg: AC.Analyzer_Config,
    base: str,
    eras: Sequence[str],
) -> Dict[str, ROOT.TChain]:
    chains = {
        "Data": ROOT.TChain("inclusive", "chain_Data"),
        "MC": ROOT.TChain("inclusive", "chain_MC"),
    }
    added = {"Data": 0, "MC": 0}

    for era in eras:
        if add_file(chains["Data"], os.path.join(base, "Data", f"{era}.root")):
            added["Data"] += 1

        if add_file(chains["MC"], os.path.join(base, "DYJetsToLL", f"{era}.root")):
            added["MC"] += 1

        for subsample in dyg_subsamples_for_era(cfg, era):
            if add_file(chains["MC"], os.path.join(base, subsample, f"{era}.root")):
                added["MC"] += 1

    for sample, nfiles in added.items():
        print(f"[Files] {sample}: added {nfiles} files, entries={chains[sample].GetEntries()}")
    return chains


def has_branch(chain: ROOT.TChain, branch: str) -> bool:
    branches = chain.GetListOfBranches()
    return bool(branches and branches.FindObject(branch))


def missing_branches(chain: ROOT.TChain, branches: Iterable[str]) -> List[str]:
    return [branch for branch in branches if not has_branch(chain, branch)]


def resolve_mass_branch(chains: Dict[str, ROOT.TChain]) -> Optional[str]:
    for branch in ("H_m", "H_mass"):
        if all(has_branch(chain, branch) for chain in chains.values()):
            return branch
    return None


def denominator_cut(branches: Sequence[str]) -> str:
    if not branches:
        return "1"
    return " && ".join(f"{branch} != 0" for branch in branches)


def draw_hist(
    chain: ROOT.TChain,
    hist_name: str,
    plot: PlotSpec,
    weight: str,
    mass_branch: str,
    extra_selection: str = "1",
    data: bool = False,
) -> ROOT.TH1D:
    old = ROOT.gROOT.FindObject(hist_name)
    if old:
        old.Delete()

    hist = ROOT.TH1D(hist_name, hist_name, plot.nbins, plot.xmin, plot.xmax)
    hist.Sumw2()
    selection = f"({sideband_cut(mass_branch)}) && ({plot.selection}) && ({extra_selection})"
    draw_weight = "1.0" if data else weight
    cut = f"({draw_weight}) * ({selection})"
    chain.Draw(f"{plot.expr}>>{hist_name}", cut, "goff")
    hist.SetDirectory(0)
    add_overflow(hist)
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


def scale_to_unit(hist: ROOT.TH1) -> None:
    integral = hist.Integral()
    if integral > 0:
        hist.Scale(1.0 / integral)


def make_ratio(data: ROOT.TH1, mc: ROOT.TH1, name: str) -> Optional[ROOT.TH1]:
    if mc.Integral() <= 0:
        return None
    ratio = data.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(mc)
    return ratio


def decorate_histograms(data: ROOT.TH1, after: ROOT.TH1, before: ROOT.TH1) -> None:
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.9)
    data.SetLineColor(ROOT.kBlack)
    data.SetMarkerColor(ROOT.kBlack)

    after.SetLineColor(ROOT.kRed + 1)
    after.SetLineWidth(2)
    after.SetFillStyle(0)
    after.SetMarkerColor(ROOT.kRed + 1)

    before.SetLineColor(ROOT.kAzure + 1)
    before.SetLineWidth(2)
    before.SetLineStyle(2)
    before.SetFillStyle(0)
    before.SetMarkerColor(ROOT.kAzure + 1)


def draw_plot(
    data: ROOT.TH1,
    after: ROOT.TH1,
    before: ROOT.TH1,
    plot: PlotSpec,
    weight: WeightSpec,
    era_label: str,
    out_path: str,
    log_y: bool,
) -> None:
    decorate_histograms(data, after, before)

    canvas = ROOT.TCanvas(f"c_{era_label}_{weight.name}_{plot.name}", "", 800, 800)
    upper = ROOT.TPad("upper", "upper", 0.0, 0.30, 1.0, 1.0)
    lower = ROOT.TPad("lower", "lower", 0.0, 0.0, 1.0, 0.30)
    upper.SetBottomMargin(0.03)
    upper.SetLeftMargin(0.13)
    upper.SetRightMargin(0.05)
    lower.SetTopMargin(0.03)
    lower.SetBottomMargin(0.35)
    lower.SetLeftMargin(0.13)
    lower.SetRightMargin(0.05)
    if log_y:
        upper.SetLogy()

    upper.Draw()
    lower.Draw()

    upper.cd()
    ymax = max(data.GetMaximum(), after.GetMaximum(), before.GetMaximum())
    if log_y:
        positive = [h.GetMinimum(0.0) for h in (data, after, before) if h.GetMinimum(0.0) > 0]
        ymin = min(positive) * 0.5 if positive else 0.01
        after.SetMinimum(ymin)
        after.SetMaximum(ymax * 20.0)
    else:
        after.SetMinimum(0.0)
        after.SetMaximum(ymax * 1.45 if ymax > 0 else 1.0)

    after.GetXaxis().SetLabelSize(0)
    after.GetYaxis().SetTitle("Events")
    after.GetYaxis().SetTitleOffset(1.15)
    after.Draw("hist")
    before.Draw("hist same")
    data.Draw("E1 same")

    legend = ROOT.TLegend(0.52, 0.70, 0.90, 0.89)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.AddEntry(data, "Data sideband", "lep")
    legend.AddEntry(after, weight.after_label, "l")
    legend.AddEntry(before, weight.before_label, "l")
    legend.Draw()

    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.04)
    label.DrawLatex(0.16, 0.86, era_label)
    label.DrawLatex(0.16, 0.81, "Sideband: 95 < m_{ll#gamma#gamma} < 180 GeV, excluding 115-135 GeV")

    lower.cd()
    ratio_after = make_ratio(data, after, f"ratio_after_{era_label}_{weight.name}_{plot.name}")
    ratio_before = make_ratio(data, before, f"ratio_before_{era_label}_{weight.name}_{plot.name}")
    frame = ROOT.TH1D(
        f"frame_{era_label}_{weight.name}_{plot.name}",
        "",
        plot.nbins,
        plot.xmin,
        plot.xmax,
    )
    frame.SetDirectory(0)
    frame.GetYaxis().SetTitle("Data / MC")
    frame.GetXaxis().SetTitle(plot.title)
    frame.GetYaxis().SetRangeUser(0.45, 1.55)
    frame.GetYaxis().SetNdivisions(505)
    frame.GetYaxis().SetTitleSize(0.10)
    frame.GetYaxis().SetTitleOffset(0.52)
    frame.GetYaxis().SetLabelSize(0.085)
    frame.GetXaxis().SetTitleSize(0.11)
    frame.GetXaxis().SetTitleOffset(1.10)
    frame.GetXaxis().SetLabelSize(0.09)
    frame.Draw("axis")

    line = ROOT.TLine(plot.xmin, 1.0, plot.xmax, 1.0)
    line.SetLineColor(ROOT.kGray + 2)
    line.SetLineStyle(2)
    line.Draw("same")

    if ratio_after:
        ratio_after.SetMarkerStyle(20)
        ratio_after.SetMarkerSize(0.75)
        ratio_after.SetMarkerColor(ROOT.kRed + 1)
        ratio_after.SetLineColor(ROOT.kRed + 1)
        ratio_after.Draw("E1 same")
    if ratio_before:
        ratio_before.SetMarkerStyle(24)
        ratio_before.SetMarkerSize(0.75)
        ratio_before.SetMarkerColor(ROOT.kAzure + 1)
        ratio_before.SetLineColor(ROOT.kAzure + 1)
        ratio_before.Draw("E1 same")

    canvas.cd()
    canvas.SaveAs(out_path + ".pdf")
    canvas.SaveAs(out_path + ".png")
    canvas.Close()


def resolve_eras(cfg: AC.Analyzer_Config, requested: str) -> List[Tuple[str, Sequence[str]]]:
    all_eras = list(cfg.years_dyll)
    if requested.strip().lower() == "all":
        return [(era, [era]) for era in all_eras] + [("run3_all", all_eras)]

    groups = []
    for era in [item.strip() for item in requested.split(",") if item.strip()]:
        if era == "run3_all":
            groups.append((era, all_eras))
        else:
            groups.append((era, [era]))
    return groups


def main() -> None:
    args = parse_args()
    set_style()

    cfg = AC.Analyzer_Config("inclusive", args.year, 2, False)
    if args.year != "run3":
        print("[Warning] This script was built for run3 layout; proceeding with configured sample_loc.")

    input_dir = args.input_dir or cfg.sample_loc
    out_dir = args.out_dir or os.path.join(cfg.out_dir, "SF_validation")
    era_groups = resolve_eras(cfg, args.eras)

    os.makedirs(out_dir, exist_ok=True)
    print(f"[Config] input_dir={input_dir}")
    print(f"[Config] out_dir={out_dir}")

    for era_label, eras in era_groups:
        print(f"\n[Era] {era_label}: {', '.join(eras)}")
        chains = build_chains(cfg, input_dir, eras)
        if chains["Data"].GetEntries() <= 0 or chains["MC"].GetEntries() <= 0:
            print(f"[Skip] {era_label}: empty Data or MC chain")
            continue
        mass_branch = resolve_mass_branch(chains)
        if not mass_branch:
            print(f"[Skip] {era_label}: no common H_m/H_mass branch in Data and MC")
            continue
        print(f"[Mass] using {mass_branch}")

        for weight in WEIGHTS:
            for plot in PLOTS:
                if plot.flavor not in weight.plot_flavors:
                    continue

                required_data = [plot.expr, mass_branch]
                if plot.selection != "1":
                    required_data.extend(token.strip() for token in plot.selection.split("==")[:1])
                required_mc = required_data + ["weight_central"] + list(weight.denominator_branches)

                missing_data = missing_branches(chains["Data"], required_data)
                missing_mc = missing_branches(chains["MC"], required_mc)
                if missing_data or missing_mc:
                    print(
                        f"[Skip] {era_label}/{weight.name}/{plot.name}: "
                        f"missing Data={missing_data}, MC={missing_mc}"
                    )
                    continue

                out_subdir = os.path.join(out_dir, era_label, weight.name)
                os.makedirs(out_subdir, exist_ok=True)
                tag = f"{era_label}_{weight.name}_{plot.name}"

                print(f"[Draw] {tag}")
                data_hist = draw_hist(
                    chains["Data"],
                    f"h_data_{tag}",
                    plot,
                    "1.0",
                    mass_branch,
                    data=True,
                )
                after_hist = draw_hist(
                    chains["MC"],
                    f"h_after_{tag}",
                    plot,
                    weight.after_weight,
                    mass_branch,
                )
                before_hist = draw_hist(
                    chains["MC"],
                    f"h_before_{tag}",
                    plot,
                    weight.before_weight,
                    mass_branch,
                    denominator_cut(weight.denominator_branches),
                )

                if args.normalize:
                    for hist in (data_hist, after_hist, before_hist):
                        scale_to_unit(hist)

                out_path = os.path.join(out_subdir, plot.name)
                draw_plot(data_hist, after_hist, before_hist, plot, weight, era_label, out_path, args.log_y)

    print("\nDone")


if __name__ == "__main__":
    main()
