####################################################
# plot_validate_zmmg.py
#
# Data/MC agreement validation for the privately selected
# Z -> mu mu gamma tag-and-probe samples, focused on the photon ID
# variables. The plotting style (two-pad Data/MC + ratio, CMS labels,
# green Z->mumugamma stack, MC stat-err band) follows 2_plot_dataVmc.py
# and reproduces the figures in
#   doc/HZa/Papers/202601_Zmmg_PT10.pdf
#
# Input ntuples (tree: tnpPhoIDs/fitter_tree):
#   Data   : run3_tnp_zmmg/Data_tnp_zmmg/<year>.root
#   MC NLO : run3_tnp_zmmg/DYJetsTo2Mu/<year>.root      (amcatnlo, default)
#   MC LO  : run3_tnp_zmmg/DYJetsToLL_MLM/<year>.root
####################################################

import argparse
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir))
sys.path.insert(0, os.path.join(PLOT_DIR, "lib"))

import ROOT  # noqa: E402
import tdrstyle  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

TREE_NAME = "tnpPhoIDs/fitter_tree"

# Base directory holding the per-sample subfolders.
SAMPLE_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_tnp_zmmg"

# MC sample folder per "--mc" choice.
MC_DIRS = {
    "nlo": "DYJetsTo2Mu",       # amcatnlo DYto2Mu (default, 2024 only)
    "lo": "DYJetsToLL_MLM",     # LO MadGraph MLM
}

DATA_DIR = "Data_tnp_zmmg"

# Integrated luminosity (fb^-1) used only for the on-plot label.
LUMI_MAP = {
    "2022preEE": 8.0,
    "2022postEE": 27.0,
    "2023preBPix": 18.0,
    "2023postBPix": 10.0,
    "2024": 109.0,
    "run3": 172.0,
}

# Colours matching the reference figures.
COL_MC_FILL = "#ffa90e"     # Z->mumugamma stack -> COLOR B (kP10Yellow)
COL_MC_LINE = "#e76300"     # darker edge (kP10Orange) for the yellow stack
COL_STAT = "#000000"        # MC Stat. Err. band -> black, hatched (hash)

MC_WEIGHT = "totWeight"     # genWeight x SFs (amcatnlo: signed); Data has totWeight==1

# Eta-region selections. EB/EE follow the H->gg fiducial gap.
REGIONS = {
    "incl": {"cut": "1", "label": ""},
    "EB": {"cut": "ph_isEB==1", "label": "|#eta| < 1.4442"},
    "EE": {"cut": "ph_isEE==1", "label": "|#eta| > 1.566"},
}

# -------------------------------------------------------------------------
# Variable definitions.
#   key      : output file stem
#   branch   : tree branch (or expression)
#   title    : x-axis title
#   nbins/xlo/xhi : binning
#   regions  : which eta regions to draw
#   logy     : log-scale main pad
#   group    : "phoid" (photon ID / shower-shape / isolation) or "kin"
# -------------------------------------------------------------------------
VARIABLES = [
    # ---- Photon ID variables (EB/EE separated, as in the reference) ----
    dict(key="phoIDMVA", branch="ph_mvaID", title="Photon ID MVA score",
         nbins=100, xlo=-1.0, xhi=1.0, regions=["EB", "EE"], logy=False, group="phoid"),
    dict(key="hoe_PUcorr", branch="ph_hoe_PUcorr", title="PU corr. H/E",
         nbins=50, xlo=0.0, xhi=0.05, regions=["EB", "EE"], logy=True, group="phoid"),
    dict(key="pfChargedIso", branch="ph_pfChargedIso", title="Relative pf charged hadron Iso (GeV)",
         nbins=50, xlo=0.0, xhi=1.0, regions=["EB", "EE"], logy=True, group="phoid"),
    dict(key="r9", branch="ph_r9", title="R_{9}",
         nbins=60, xlo=0.4, xhi=1.0, regions=["EB", "EE"], logy=False, group="phoid"),
    dict(key="sieie", branch="ph_sieie", title="#sigma_{i#eta i#eta}",
         nbins=60, xlo=0.0, xhi=0.06, regions=["EB", "EE"], logy=True, group="phoid"),
    dict(key="s4", branch="ph_s4", title="S_{4}",
         nbins=50, xlo=0.4, xhi=1.0, regions=["EB", "EE"], logy=False, group="phoid"),
    dict(key="etaWidth", branch="ph_etaWidth", title="#eta width",
         nbins=50, xlo=0.0, xhi=0.03, regions=["EB", "EE"], logy=False, group="phoid"),
    dict(key="phiWidth", branch="ph_phiWidth", title="#phi width",
         nbins=50, xlo=0.0, xhi=0.15, regions=["EB", "EE"], logy=False, group="phoid"),
    dict(key="ecalIso", branch="ph_ecalIso", title="ECAL Iso (GeV)",
         nbins=50, xlo=0.0, xhi=10.0, regions=["EB", "EE"], logy=True, group="phoid"),
    # ---- Kinematics / selection context (inclusive) ----
    dict(key="gammaPt", branch="ph_et", title="p_{T} of #gamma (GeV)",
         nbins=45, xlo=10.0, xhi=100.0, regions=["incl"], logy=False, group="kin"),
    dict(key="gammaEta", branch="ph_sc_eta", title="#eta_{#gamma}",
         nbins=50, xlo=-2.5, xhi=2.5, regions=["incl"], logy=False, group="kin"),
    dict(key="gammaPhi", branch="ph_phi", title="#phi_{#gamma}",
         nbins=32, xlo=-3.2, xhi=3.2, regions=["incl"], logy=False, group="kin"),
    dict(key="mMuMu", branch="tag_mumu_mass", title="m_{#mu#mu} (GeV)",
         nbins=60, xlo=30.0, xhi=90.0, regions=["incl"], logy=False, group="kin"),
    dict(key="mMuMuGamma", branch="pair_mass", title="m_{#mu#mu#gamma} (GeV)",
         nbins=40, xlo=80.0, xhi=100.0, regions=["incl"], logy=False, group="kin"),
    dict(key="dRGammaNearMu", branch="ph_lepNearDR", title="#DeltaR(#gamma, near #mu)",
         nbins=40, xlo=0.4, xhi=0.8, regions=["incl"], logy=False, group="kin"),
]


def _resolve_input(year, mc):
    data_path = os.path.join(SAMPLE_BASE, DATA_DIR, f"{year}.root")
    mc_path = os.path.join(SAMPLE_BASE, MC_DIRS[mc], f"{year}.root")
    for label, path in (("Data", data_path), ("MC", mc_path)):
        if not os.path.exists(path):
            raise RuntimeError(f"Missing {label} input: {path}")
    return data_path, mc_path


def _fill_hist(tree, name, var, region_cut, weight=None, truth_match=False):
    # Book in the current directory so TTree::Draw can locate it by name.
    h = ROOT.TH1F(name, name, var["nbins"], var["xlo"], var["xhi"])
    h.Sumw2()

    cut = f"({region_cut})"
    if truth_match:
        cut = f"({cut}) && (ph_genPartFlav==1)"
    if weight:
        cut = f"({weight})*({cut})"

    tree.Draw(f"{var['branch']}>>{name}", cut, "goff")
    h.SetDirectory(0)
    # Include the overflow/underflow into the edge bins for a fair area norm.
    nb = h.GetNbinsX()
    h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    h.SetBinContent(nb, h.GetBinContent(nb) + h.GetBinContent(nb + 1))
    h.SetBinContent(0, 0.0)
    h.SetBinContent(nb + 1, 0.0)
    return h


def _stat_band(h_mc):
    """MC stat-error band as a filled TGraphAsymmErrors (absolute scale)."""
    g = ROOT.TGraphAsymmErrors()
    g.SetName(h_mc.GetName() + "_statband")
    for i in range(1, h_mc.GetNbinsX() + 1):
        x = h_mc.GetBinCenter(i)
        y = h_mc.GetBinContent(i)
        ex = 0.5 * h_mc.GetBinWidth(i)
        ey = h_mc.GetBinError(i)
        n = g.GetN()
        g.SetPoint(n, x, y)
        g.SetPointError(n, ex, ex, ey, ey)
    g.SetFillColorAlpha(ROOT.TColor.GetColor(COL_STAT), 1.0)
    g.SetLineColor(ROOT.TColor.GetColor(COL_STAT))   # hatch lines drawn in fill colour
    g.SetFillStyle(3354)                              # hatched (hash / cross-hatch) uncertainty
    g.SetLineWidth(0)
    g.SetMarkerStyle(0)
    return g


def _ratio_band(h_mc):
    """MC stat-error band normalised to the MC central value (for the ratio pad)."""
    g = ROOT.TGraphAsymmErrors()
    g.SetName(h_mc.GetName() + "_ratioband")
    for i in range(1, h_mc.GetNbinsX() + 1):
        x = h_mc.GetBinCenter(i)
        y = h_mc.GetBinContent(i)
        ex = 0.5 * h_mc.GetBinWidth(i)
        rel = (h_mc.GetBinError(i) / y) if y > 0 else 0.0
        n = g.GetN()
        g.SetPoint(n, x, 1.0)
        g.SetPointError(n, ex, ex, rel, rel)
    g.SetFillColorAlpha(ROOT.TColor.GetColor(COL_STAT), 1.0)
    g.SetLineColor(ROOT.TColor.GetColor(COL_STAT))   # hatch lines drawn in fill colour
    g.SetFillStyle(3354)                              # hatched (hash / cross-hatch) uncertainty
    g.SetLineWidth(0)
    g.SetMarkerStyle(0)
    return g


def _make_lumi_label(lumi):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.05)
    tex.SetTextAlign(31)
    tex.DrawLatex(0.95, 0.940, f"{lumi:.0f} fb^{{-1}} (13.6 TeV)")
    return tex


def _make_cms_label():
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(61)
    tex.SetTextSize(0.058)
    tex.SetTextAlign(11)
    tex.DrawLatex(0.16, 0.94, "CMS")
    sub = ROOT.TLatex()
    sub.SetNDC()
    sub.SetTextFont(52)
    sub.SetTextSize(0.045)
    sub.SetTextAlign(11)
    sub.DrawLatex(0.255, 0.94, "Preliminary")
    return tex, sub


def _draw_one(var, region_key, h_data, h_mc, lumi, out_dir):
    region = REGIONS[region_key]
    # Shift the y-axis exponent (x10^N) to the left so it doesn't collide with the CMS label.
    ROOT.TGaxis.SetExponentOffset(-0.075, 0.0, "y")

    # Area-normalise MC to data (MC normalized to #data).
    data_int = h_data.Integral()
    mc_int = h_mc.Integral()
    if mc_int > 0 and data_int > 0:
        h_mc.Scale(data_int / mc_int)

    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(1.0)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetLineWidth(1)

    h_mc.SetFillColor(ROOT.TColor.GetColor(COL_MC_FILL))
    h_mc.SetLineColor(ROOT.TColor.GetColor(COL_MC_LINE))
    h_mc.SetLineWidth(1)

    stat_band = _stat_band(h_mc)
    ratio_band = _ratio_band(h_mc)

    # Data / MC ratio points.
    ratio = ROOT.TGraphAsymmErrors()
    ratio.Divide(h_data, h_mc, "pois")
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.0)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetLineColor(ROOT.kBlack)

    canv = ROOT.TCanvas("c_" + var["key"] + "_" + region_key, "", 800, 800)
    canv.cd()

    L, R, B, T = 0.16, 0.05, 0.19, 0.085

    upper = ROOT.TPad("upper", "upper", 0, 0.31, 1, 1)
    upper.SetLeftMargin(L)
    upper.SetRightMargin(R)
    upper.SetBottomMargin(0.03)
    upper.SetTopMargin(T)
    upper.SetTickx(1)
    upper.SetTicky(1)
    if var["logy"]:
        upper.SetLogy()
    upper.Draw()
    upper.cd()

    hmax = max(h_data.GetMaximum(), h_mc.GetMaximum())
    if var["logy"]:
        h_mc.SetMinimum(0.5)
        h_mc.SetMaximum(hmax * 50.0)
    else:
        h_mc.SetMinimum(0.0)
        h_mc.SetMaximum(hmax * 1.45)

    h_mc.GetYaxis().SetTitle("Events / %.3g" % h_mc.GetBinWidth(1))
    h_mc.GetYaxis().SetTitleSize(0.07)
    h_mc.GetYaxis().SetLabelSize(0.055)
    h_mc.GetYaxis().SetTitleFont(42)
    h_mc.GetYaxis().SetTitleOffset(1.25)
    h_mc.GetXaxis().SetLabelSize(0)
    h_mc.Draw("HIST")
    stat_band.Draw("2 SAME")
    h_data.Draw("PE SAME")
    h_data.Draw("AXIS SAME")

    legend = ROOT.TLegend(0.62, 0.66, 0.92, 0.86)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.045)
    legend.AddEntry(h_data, "Data", "PE")
    legend.AddEntry(h_mc, "Z #rightarrow #mu#mu#gamma", "f")
    legend.AddEntry(stat_band, "MC Stat. Uncert.", "f")
    legend.Draw("SAME")

    if region["label"]:
        reg_tex = ROOT.TLatex()
        reg_tex.SetNDC()
        reg_tex.SetTextFont(42)
        reg_tex.SetTextSize(0.045)
        reg_tex.SetTextAlign(13)
        reg_tex.DrawLatex(0.20, 0.83, region["label"])

    _make_cms_label()
    _make_lumi_label(lumi)

    # ---- ratio pad ----
    canv.cd()
    lower = ROOT.TPad("lower", "lower", 0, 0, 1, 0.29)
    lower.SetTopMargin(0.00001)
    lower.SetBottomMargin(0.36)
    lower.SetLeftMargin(L)
    lower.SetRightMargin(R)
    lower.SetTickx(1)
    lower.SetTicky(1)
    lower.SetGridy()
    lower.Draw()
    lower.cd()

    frame = ROOT.TH1F("frame_" + var["key"] + "_" + region_key, "",
                      var["nbins"], var["xlo"], var["xhi"])
    frame.SetMinimum(0.5)
    frame.SetMaximum(1.5)
    frame.SetStats(0)
    frame.GetYaxis().SetTitle("Data / MC")
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetNdivisions(505)
    frame.GetYaxis().SetTitleSize(0.13)
    frame.GetYaxis().SetTitleOffset(0.58)
    frame.GetYaxis().SetLabelSize(0.11)
    frame.GetXaxis().SetTitle(var["title"])
    frame.GetXaxis().SetTitleSize(0.15)
    frame.GetXaxis().SetTitleOffset(1.0)
    frame.GetXaxis().SetLabelSize(0.12)
    frame.GetXaxis().SetTickLength(0.07)
    frame.Draw("AXIS")

    ratio_band.Draw("2 SAME")

    line = ROOT.TLine(var["xlo"], 1.0, var["xhi"], 1.0)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kBlack)
    line.Draw("SAME")

    ratio.Draw("PZ SAME")
    frame.Draw("AXIS SAME")

    canv.cd()
    out_name = f"{var['group']}_{var['key']}"
    if region_key != "incl":
        out_name += f"_{region_key}"
    canv.SaveAs(os.path.join(out_dir, out_name + ".pdf"))
    canv.Close()


def main():
    parser = argparse.ArgumentParser(description="Z->mumugamma data/MC photon-ID validation plots")
    parser.add_argument("-y", "--year", default="2024",
                        help="era tag of the input ROOT files (default: 2024)")
    parser.add_argument("--mc", choices=["nlo", "lo"], default="nlo",
                        help="MC sample: nlo=DYJetsTo2Mu (amcatnlo, default), lo=DYJetsToLL_MLM")
    parser.add_argument("--outDir", default=None,
                        help="output directory (default: Plot/plots/validate_zmmg/<year>_<mc>)")
    parser.add_argument("--truthMatch", action="store_true",
                        help="require ph_genPartFlav==1 (prompt photon) on MC")
    parser.add_argument("--group", choices=["all", "phoid", "kin"], default="all",
                        help="which variable group to draw")
    args = parser.parse_args()

    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetErrorX(0.5)
    # Use a x10^n multiplier (as in the reference figures) so wide event
    # counts do not collide with the y-axis title.
    ROOT.TGaxis.SetMaxDigits(3)

    data_path, mc_path = _resolve_input(args.year, args.mc)
    lumi = LUMI_MAP.get(args.year, 0.0)

    out_dir = args.outDir or os.path.join(
        PLOT_DIR, "plots", "validate_zmmg", f"{args.year}_{args.mc}")
    os.makedirs(out_dir, exist_ok=True)

    print(f"[Input] Data : {data_path}")
    print(f"[Input] MC   : {mc_path}  ({args.mc})")
    print(f"[Output] dir : {out_dir}")
    print(f"[Config] lumi label = {lumi:.0f} fb^-1, truthMatch = {args.truthMatch}")

    f_data = ROOT.TFile.Open(data_path)
    f_mc = ROOT.TFile.Open(mc_path)
    t_data = f_data.Get(TREE_NAME)
    t_mc = f_mc.Get(TREE_NAME)
    if not t_data or not t_mc:
        raise RuntimeError("Could not read tree '%s' from inputs." % TREE_NAME)

    n_drawn = 0
    for var in VARIABLES:
        if args.group != "all" and var["group"] != args.group:
            continue
        for region_key in var["regions"]:
            region_cut = REGIONS[region_key]["cut"]
            tag = f"{var['key']}_{region_key}"
            h_data = _fill_hist(t_data, "hdata_" + tag, var, region_cut)
            h_mc = _fill_hist(t_mc, "hmc_" + tag, var, region_cut,
                              weight=MC_WEIGHT, truth_match=args.truthMatch)

            if h_data.Integral() <= 0 or h_mc.Integral() <= 0:
                print(f"[Warning] empty histogram for {tag}; skip.")
                continue

            _draw_one(var, region_key, h_data, h_mc, lumi, out_dir)
            n_drawn += 1
            print(f"[Done] {tag}")

    f_data.Close()
    f_mc.Close()
    print(f"\nDrew {n_drawn} plots into {out_dir}")


if __name__ == "__main__":
    main()
