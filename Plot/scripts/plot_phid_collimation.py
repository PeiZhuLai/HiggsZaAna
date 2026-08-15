#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Signal-MC characterization of the custom photon-ID efficiency vs diphoton
collimation, for the L3-convener question (AN-25-172 L319 / Sec 4.3):
  "Does the MC signal efficiency match data in this (collimated) regime?"

Idea: each HZa signal sample is a fixed ALP mass mA, and mA is a clean label for
the diphoton opening angle dR(gg) (small mA => merged, large mA => resolved).
We read the per-(mA,era) cutflow JSONs and measure the custom photon-ID
efficiency as a function of mA for 5 candidate photon-ID definitions, showing
that the analysis' custom ID (H/E + chIso + HCal-iso, NO shower-shape) is the
only one whose efficiency is flat across the merged->resolved transition, while
the shower-shape / ECAL-iso based IDs degrade by 10-21% in the merged regime.

Inputs  : Plot/output/cutflow_list/cutflow_Sig_MC_mA_M<ma>_<era>.json (70 files)
          review_comments/_gendR_map.json  (gen dR(gg) vs mA, from parquet)
Outputs : Plot/plots/photonid_collimation/{phid_eff_vs_mA_custom,
          phid_eff_vs_mA_scenarios, mA_to_dRgg_map}.{pdf,png}
Env     : conda higgs-zg-plot_python  (PyROOT 6.26)
"""
import os, json, array
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# ---------------------------------------------------------------- config
CUTFLOW_DIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list"
OUT_DIR     = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/photonid_collimation"
GENDR_JSON  = os.path.join(OUT_DIR, "_gendR_map.json")
os.makedirs(OUT_DIR, exist_ok=True)

MA    = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
ERAS  = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
# per-era lumi (fb-1) for labels only; combined = sum of eras used
LUMI  = {"2022preEE": 7.99, "2022postEE": 26.68, "2023preBPix": 17.96,
         "2023postBPix": 9.68, "2024": 109.82}
LUMI_RUN3 = sum(LUMI.values())

# exact definitions (higgs_dna/selections/photon_selections.py:352-360):
#   custom        = H/E & chIso & HCaliso          <-- nominal analysis ID
#   custom_extend = custom & sieie
#   official      = custom & sieie & ECAL-iso       (full EGamma-style tight)
#   sieie         = sieie cut ALONE
#   PFECalIso     = ECAL-iso cut ALONE
# Legend order: nested trio first (the collimation story), single cuts after.
KINDS = ["custom", "custom_extend", "official", "sieie", "PFECalIso"]
KIND_LABEL = {
    "custom":        "custom = H/E+chIso+HCaliso [nominal]",
    "custom_extend": "custom + #sigma_{i#etai#eta}",
    "official":      "custom + #sigma_{i#etai#eta} + ECAL-iso",
    "sieie":         "#sigma_{i#etai#eta} only",
    "PFECalIso":     "ECAL-iso only",
}
KIND_COLOR = {
    "custom":        ROOT.kRed + 1,
    "custom_extend": ROOT.kOrange + 7,
    "official":      ROOT.kMagenta + 2,
    "sieie":         ROOT.kGreen + 2,
    "PFECalIso":     ROOT.kAzure + 1,
}
KIND_MARKER = {"custom": 20, "custom_extend": 21, "official": 33,
               "sieie": 24, "PFECalIso": 25}
KIND_LSTYLE = {"custom": 1, "custom_extend": 1, "official": 1,
               "sieie": 2, "PFECalIso": 2}

ERA_COLOR = {"2022preEE": ROOT.kAzure + 1, "2022postEE": ROOT.kGreen + 2,
             "2023preBPix": ROOT.kOrange + 7, "2023postBPix": ROOT.kMagenta + 2,
             "2024": ROOT.kGray + 2}
ERA_MARKER = {"2022preEE": 24, "2022postEE": 25, "2023preBPix": 26,
              "2023postBPix": 32, "2024": 27}


# ---------------------------------------------------------------- read cutflows
def load():
    """Return nested dict of summed weighted counts: D[era][ma][field]."""
    D = {e: {m: {} for m in MA} for e in ERAS}
    for ma in MA:
        for era in ERAS:
            fn = os.path.join(CUTFLOW_DIR, f"cutflow_Sig_MC_mA_M{ma}_{era}.json")
            with open(fn) as f:
                cf = json.load(f)["cutflows"]
            z = cf["zgammas_w"]
            zu = cf["zgammas"]  # unweighted, for binomial stat errors
            d = {
                "ph_eta":        z["ph_eta"],
                "ph_id_hoe":     z["ph_id_hoe"],
                "ph_id_chiso":   z["ph_id_chiso"],
                "ph_id_hcaliso": z["ph_id_hcaliso"],
                "ph_eta_u":        zu["ph_eta"],
                "ph_id_hcaliso_u": zu["ph_id_hcaliso"],
            }
            for k in KINDS:
                c = cf[f"zgammas_phid_{k}_tight_w"]
                d[f"{k}_den"] = c["g_kin_cut"]
                d[f"{k}_num"] = c["has_2g_cand"]
            D[era][ma] = d
    return D


def eff_curve_scenario(D, kind, eras):
    """Run3(-combined over `eras`) ID-scenario efficiency vs mA."""
    ys = []
    for ma in MA:
        num = sum(D[e][ma][f"{kind}_num"] for e in eras)
        den = sum(D[e][ma][f"{kind}_den"] for e in eras)
        ys.append(num / den if den else 0.0)
    return ys


def eff_curve_nominal(D, era_or_eras):
    """Nominal custom-ID eff = ph_id_hcaliso / ph_eta.
    Central value from weighted counts; statistical error bar from the
    (unweighted) event counts as a binomial uncertainty sqrt(p(1-p)/N)."""
    import math
    eras = era_or_eras if isinstance(era_or_eras, list) else [era_or_eras]
    ys, errs = [], []
    for ma in MA:
        num = sum(D[e][ma]["ph_id_hcaliso"] for e in eras)
        den = sum(D[e][ma]["ph_eta"] for e in eras)
        ys.append(num / den if den else 0.0)
        nu = sum(D[e][ma]["ph_eta_u"] for e in eras)
        ku = sum(D[e][ma]["ph_id_hcaliso_u"] for e in eras)
        pu = ku / nu if nu else 0.0
        errs.append(math.sqrt(pu * (1 - pu) / nu) if nu else 0.0)
    return ys, errs


# ---------------------------------------------------------------- style helpers
def cms_lumi(pad, lumi_fb):
    pad.cd()
    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.045)
    lat.SetTextAlign(11)
    lat.DrawLatex(pad.GetLeftMargin(), 0.945, "#bf{CMS} #it{Preliminary}")
    lat.SetTextAlign(31)
    lat.DrawLatex(1 - pad.GetRightMargin(), 0.945,
                  f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")


def frame(pad, xlo, xhi, ylo, yhi, xt, yt):
    h = pad.DrawFrame(xlo, ylo, xhi, yhi)
    h.GetXaxis().SetTitle(xt); h.GetYaxis().SetTitle(yt)
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetTitleSize(0.055); ax.SetLabelSize(0.050)
    h.GetXaxis().SetTitleOffset(1.15); h.GetYaxis().SetTitleOffset(1.20)
    return h


def tgraph(xs, ys, color, marker, lstyle=1, msize=1.7):
    g = ROOT.TGraph(len(xs), array.array('d', xs), array.array('d', ys))
    g.SetLineColor(color); g.SetMarkerColor(color)
    g.SetMarkerStyle(marker); g.SetLineWidth(3); g.SetLineStyle(lstyle)
    g.SetMarkerSize(msize)
    return g


def tgraph_err(xs, ys, eys, color, marker, lstyle=1, msize=1.7, lw=3):
    n = len(xs)
    g = ROOT.TGraphErrors(n, array.array('d', xs), array.array('d', ys),
                          array.array('d', [0.0] * n), array.array('d', eys))
    g.SetLineColor(color); g.SetMarkerColor(color)
    g.SetMarkerStyle(marker); g.SetLineWidth(lw); g.SetLineStyle(lstyle)
    g.SetMarkerSize(msize)
    return g


def save(c, name):
    for ext in ("pdf", "png"):
        c.SaveAs(os.path.join(OUT_DIR, f"{name}.{ext}"))


# ---------------------------------------------------------------- plots
def plot_scenarios(D, gendr):
    """THE key plot: 5 ID scenarios vs mA, Run3-combined."""
    c = ROOT.TCanvas("cscen", "", 900, 700)
    c.SetLeftMargin(0.16); c.SetRightMargin(0.05)
    c.SetTopMargin(0.08);  c.SetBottomMargin(0.13)
    frame(c, 0, 31, 0.40, 1.06, "m_{a} [GeV]",
          "photon-ID efficiency (#geq2#gamma, per event)")
    graphs = []
    leg = ROOT.TLegend(0.30, 0.15, 0.95, 0.39)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42)
    leg.SetTextSize(0.040)
    for k in KINDS:
        ys = eff_curve_scenario(D, k, ERAS)
        g = tgraph(MA, ys, KIND_COLOR[k], KIND_MARKER[k], lstyle=KIND_LSTYLE[k])
        g.Draw("PL SAME"); graphs.append(g)
        leg.AddEntry(g, KIND_LABEL[k], "pl")
    leg.Draw()
    cms_lumi(c, LUMI_RUN3)
    save(c, "phid_eff_vs_mA_scenarios")


def plot_custom_per_era(D):
    c = ROOT.TCanvas("cera", "", 900, 700)
    c.SetLeftMargin(0.16); c.SetRightMargin(0.05)
    c.SetTopMargin(0.08);  c.SetBottomMargin(0.13)
    frame(c, 0, 31, 0.70, 0.97, "m_{a} [GeV]",
          "custom photon-ID efficiency (nominal)")
    graphs = []
    leg = ROOT.TLegend(0.49, 0.62, 0.94, 0.91)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42)
    leg.SetTextSize(0.040)
    for era in ERAS:
        ys, eys = eff_curve_nominal(D, era)
        g = tgraph_err(MA, ys, eys, ERA_COLOR[era], ERA_MARKER[era], msize=1.6)
        g.Draw("PL SAME"); graphs.append(g)
        leg.AddEntry(g, era, "pl")
    ycomb, eycomb = eff_curve_nominal(D, ERAS)
    gc = tgraph_err(MA, ycomb, eycomb, ROOT.kBlue + 2, 20, msize=2.0)
    gc.Draw("PL SAME"); graphs.append(gc)
    leg.AddEntry(gc, "Run 3 combined", "pl")
    leg.Draw()
    cms_lumi(c, LUMI_RUN3)
    save(c, "phid_eff_vs_mA_custom")


def plot_mapping(gendr):
    """mA -> collimation: mean gen dR(gg) + merged/transition/resolved frac."""
    c = ROOT.TCanvas("cmap", "", 900, 700)
    c.SetLeftMargin(0.16); c.SetRightMargin(0.155)
    c.SetTopMargin(0.08);  c.SetBottomMargin(0.13)
    frame(c, 0, 31, 0.0, 1.02, "m_{a} [GeV]",
          "fraction of signal events (gen)")
    # stacked-ish fractions as lines
    fr = {reg: [] for reg in ("lt01", "0103", "gt03")}
    dr_mean = []
    for ma in MA:
        mean, med, f1, f2, f3 = gendr[str(ma)]
        fr["lt01"].append(f1); fr["0103"].append(f2); fr["gt03"].append(f3)
        dr_mean.append(mean)
    regs = [("lt01", ROOT.kRed + 1, 20, "merged  #Delta R(#gamma,#gamma)<0.1"),
            ("0103", ROOT.kOrange + 7, 21, "transition 0.1#minus0.3"),
            ("gt03", ROOT.kAzure + 1, 22, "resolved  >0.3")]
    graphs = []
    leg = ROOT.TLegend(0.44, 0.15, 0.82, 0.40)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42)
    leg.SetTextSize(0.040)
    for key, col, mk, lab in regs:
        g = tgraph(MA, fr[key], col, mk)
        g.Draw("PL SAME"); graphs.append(g)
        leg.AddEntry(g, lab, "pl")
    leg.Draw()
    # overlay mean gen dR on a right axis. The frame top is FRAME_TOP (not 1.0),
    # so the curve is scaled by FRAME_TOP/drmax and the right axis maps the same
    # world range [0, FRAME_TOP] -> [0, drmax], making the two exactly consistent.
    FRAME_TOP = 1.02
    drmax = 2.2
    g_dr = tgraph(MA, [d * FRAME_TOP / drmax for d in dr_mean],
                  ROOT.kGray + 3, 34, lstyle=2)
    g_dr.Draw("PL SAME"); graphs.append(g_dr)
    raxis = ROOT.TGaxis(31, 0.0, 31, FRAME_TOP, 0.0, drmax, 510, "+L")
    raxis.SetTitle("mean #Delta R(#gamma,#gamma)  (gen)")
    raxis.SetLineColor(ROOT.kGray + 3); raxis.SetLabelColor(ROOT.kGray + 3)
    raxis.SetTitleColor(ROOT.kGray + 3)
    raxis.SetTitleFont(42); raxis.SetLabelFont(42)
    raxis.SetTitleSize(0.055); raxis.SetLabelSize(0.050)
    raxis.SetTitleOffset(1.15); raxis.Draw()
    leg.AddEntry(g_dr, "mean #Delta R(#gamma,#gamma) (right)", "pl")
    cms_lumi(c, LUMI_RUN3)
    save(c, "mA_to_dRgg_map")


def main():
    D = load()
    with open(GENDR_JSON) as f:
        gendr = json.load(f)
    plot_scenarios(D, gendr)
    plot_custom_per_era(D)
    plot_mapping(gendr)
    # dump a small numeric summary next to plots
    summ = {"mA": MA,
            "eff_custom_run3": eff_curve_scenario(D, "custom", ERAS),
            "eff_sieie_run3": eff_curve_scenario(D, "sieie", ERAS),
            "eff_official_run3": eff_curve_scenario(D, "official", ERAS),
            "eff_PFECalIso_run3": eff_curve_scenario(D, "PFECalIso", ERAS),
            "eff_custom_extend_run3": eff_curve_scenario(D, "custom_extend", ERAS),
            "gen_dRgg_mean": [gendr[str(m)][0] for m in MA]}
    with open(os.path.join(OUT_DIR, "phid_collimation_summary.json"), "w") as f:
        json.dump(summ, f, indent=2)
    print("wrote plots + summary to", OUT_DIR)


if __name__ == "__main__":
    main()
