#!/usr/bin/env python3
"""eta -> gamma gamma control sample: is there an eta peak in the merged-photon mass of
Z + X data, and does MC put it at the same place with the same width?

Leading MLPhoton (EXO-22-022 model), events outside the merged signal selection
(extract_mlphoton_eta_control.py). For each selection (diphotonScore and pT cut) the
mass in [0.30, 1.00] GeV is fitted with Gaussian + 3rd-order polynomial (binned
likelihood for data, chi2 with Sumw2 errors for weighted MC), and once more with the
background alone.
The peak position is limited to 0.45-0.65 GeV (eta = 0.548), so the significance
sqrt(delta chi2) is a local one inside that range.

Plots: data with fit, MC with fit, data vs MC shape (area-normalized), per selection.
Table: mu, sigma, N_sig, local significance for data and MC; mu_data / mu_MC and
sigma_data / sigma_MC are the regressor's data/MC scale and resolution.

PyROOT, LCG_104, hzg_style. Lumi from the HZgamma constants.py (2022 + 2023 + 2024).
Input : <out>/eta_control_{data,mc}.npz
Output: <out>/eta_<sel>_{data,mc,shape}.pdf/.png ; <tables>/eta_control_fits.txt
"""
import argparse
import math
import os
import sys

import numpy as np
import ROOT

if ROOT.gROOT.GetVersion().startswith("6.34"):
    sys.exit("ROOT 6.34 drops histograms from PDF output; use LCG_104")

HZG = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3"
sys.path.insert(0, f"{HZG}/plot/scripts")
from hzg_style import (  # noqa: E402
    HEADROOM, LINE_WIDTH, PAIR_COLORS, canvas_single, cms_header, house_legend,
    panel_label, save, style_axes, _keep,
)

OUT = "/eos/home-p/pelai/HZa_merged/output_plots/MLPhoton_eta_control"
TABLES = "/eos/home-p/pelai/HZa_merged/output_tables/MLPhoton_eta_control"
XLO, XHI, NBIN = 0.30, 1.00, 35
MU_RANGE = (0.45, 0.65)
M_ETA = 0.548
SELECTIONS = [(0.5, 15.0), (0.9, 15.0), (0.5, 25.0)]   # (diphotonScore >, pT >)
XTITLE = "m_{#gamma#gamma}^{reco} (merged photon) [GeV]"
TITLE = "Z + X, merged-photon candidate"

def load_lumi():
    """LUMI from the HZgamma constants.py. The module cannot be imported (its ERA dict
    uses sets as keys), so only the `LUMI = {...}` assignment is parsed, as
    plot_dataVmc_root_disjoint.py does."""
    import ast
    tree = ast.parse(open(f"{HZG}/higgs_dna/metaconditions/corrections/constants.py").read())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "LUMI" for t in node.targets):
            return ast.literal_eval(node.value)
    raise SystemExit("no LUMI in constants.py")


LUMI = load_lumi()

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.TH1.AddDirectory(False)


def hist(name, m, w):
    h = ROOT.TH1D(name, "", NBIN, XLO, XHI)
    h.Sumw2()
    h.FillN(len(m), m.astype(np.float64), w.astype(np.float64))
    _keep.append(h)
    return h


def fit(h, weighted):
    """Gaussian + pol3, then pol3 alone; returns (sb function, b function, result dict)."""
    # weighted (MC): chi2 with Sumw2 errors -- the weighted likelihood ("WL") crashes
    # ROOT 6.28 in DoWeightMinimization on this sample (negative weights)
    opt = "SQ0" if weighted else "LSQ0"
    c = 0.5 * (XLO + XHI)
    poly = "[3]+[4]*(x-%g)+[5]*(x-%g)^2+[6]*(x-%g)^3" % (c, c, c)
    sb = ROOT.TF1(f"sb_{h.GetName()}", "[0]*exp(-0.5*((x-[1])/[2])^2)+" + poly, XLO, XHI)
    peak = h.GetBinContent(h.FindBin(M_ETA))
    sb.SetParameters(0.2 * peak, M_ETA, 0.06, h.Integral() / NBIN, 0, 0, 0)
    sb.SetParLimits(0, 0.0, 10 * max(peak, 1e-9))
    sb.SetParLimits(1, *MU_RANGE)
    sb.SetParLimits(2, 0.02, 0.15)
    r_sb = h.Fit(sb, opt)
    b = ROOT.TF1(f"b_{h.GetName()}", poly.replace("[3]", "[0]").replace("[4]", "[1]")
                 .replace("[5]", "[2]").replace("[6]", "[3]"), XLO, XHI)
    b.SetParameters(*[sb.GetParameter(i) for i in range(3, 7)])
    r_b = h.Fit(b, opt)
    dchi2 = max(r_b.Chi2() - r_sb.Chi2(), 0.0)
    bw = (XHI - XLO) / NBIN
    nsig = sb.GetParameter(0) * sb.GetParameter(2) * math.sqrt(2 * math.pi) / bw
    res = {"mu": sb.GetParameter(1), "dmu": sb.GetParError(1),
           "sigma": sb.GetParameter(2), "dsigma": sb.GetParError(2),
           "nsig": nsig, "z": math.sqrt(dchi2), "status": int(r_sb.Status()),
           "n": h.Integral()}
    # draw the background component of the s+b fit, not the separate b-only fit
    bdraw = ROOT.TF1(f"bd_{h.GetName()}", poly.replace("[3]", "[0]").replace("[4]", "[1]")
                     .replace("[5]", "[2]").replace("[6]", "[3]"), XLO, XHI)
    bdraw.SetParameters(*[sb.GetParameter(i) for i in range(3, 7)])
    _keep.extend([sb, b, bdraw])
    return sb, bdraw, res


def draw_fit(name, h, sb, b, res, notes, lumi, is_data, out):
    c = canvas_single(name)
    label = "Data" if is_data else "Z + jets MC (2024)"
    leg = house_legend(c, [label, "Gaussian + pol3", "pol3 component"],
                       clear_label=TITLE, clear_notes=notes, size=0.038)
    h.SetMarkerStyle(20)
    h.SetMarkerSize(1.0)
    h.SetLineColor(ROOT.kBlack)
    h.SetMarkerColor(ROOT.kBlack)
    style_axes(h, c, XTITLE, "Events / %.0f MeV" % (1000 * (XHI - XLO) / NBIN))
    h.SetMinimum(0.0)
    h.SetMaximum(HEADROOM * max(h.GetMaximum(), 1e-9))
    h.Draw("E1")
    for f, style in ((sb, 1), (b, 7)):
        f.SetLineColor(PAIR_COLORS[1])
        f.SetLineStyle(style)
        f.SetLineWidth(LINE_WIDTH)
        f.SetNpx(400)
        f.Draw("SAME")
    h.Draw("E1 SAME")
    leg.AddEntry(h, label, "pe")
    leg.AddEntry(sb, "Gaussian + pol3", "l")
    leg.AddEntry(b, "pol3 component", "l")
    leg.Draw()
    cms_header(c, lumi if is_data else None, "13.6", extra="Preliminary" if is_data else "Simulation")
    panel_label(c, TITLE, notes)
    c.RedrawAxis()
    save(c, os.path.join(out, name))


def draw_shape(name, hd, hm, notes, lumi, out):
    c = canvas_single(name)
    leg = house_legend(c, ["Data", "Z + jets MC (2024)"], clear_label=TITLE,
                       clear_notes=notes)
    hd, hm = hd.Clone(name + "_d"), hm.Clone(name + "_m")
    _keep.extend([hd, hm])
    for h, col in ((hd, PAIR_COLORS[0]), (hm, PAIR_COLORS[1])):
        if h.Integral() > 0:
            h.Scale(1.0 / h.Integral())
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetLineWidth(LINE_WIDTH)
    hd.SetMarkerStyle(20)
    style_axes(hm, c, XTITLE, "A.U. / %.0f MeV" % (1000 * (XHI - XLO) / NBIN))
    hm.SetMinimum(0.0)
    hm.SetMaximum(HEADROOM * max(hm.GetMaximum(), hd.GetMaximum()))
    hm.Draw("HIST E")
    hd.Draw("E1 SAME")
    leg.AddEntry(hd, "Data", "pe")
    leg.AddEntry(hm, "Z + jets MC (2024)", "l")
    line = ROOT.TLine(M_ETA, 0.0, M_ETA, 1.05 * max(hm.GetMaximum(), hd.GetMaximum()) / HEADROOM)
    line.SetLineStyle(7)
    line.SetLineColor(ROOT.kGray + 2)
    line.SetLineWidth(2)
    line.Draw()
    _keep.append(line)
    leg.Draw()
    cms_header(c, lumi, "13.6", extra="Preliminary")
    panel_label(c, TITLE, notes)
    c.RedrawAxis()
    save(c, os.path.join(out, name))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--tables", default=TABLES)
    args = ap.parse_args()
    d = np.load(os.path.join(args.out, "eta_control_data.npz"))
    m = np.load(os.path.join(args.out, "eta_control_mc.npz"))
    lumi = LUMI["2022"] + LUMI["2023"] + LUMI["2024"]

    lines = ["# eta -> gg control: Gaussian + pol3 in [%.2f, %.2f] GeV, mu limited to %s; "
             "z = local sqrt(delta chi2) vs pol3 alone" % (XLO, XHI, MU_RANGE),
             f"{'selection':22s} {'sample':5s} {'N':>9s} {'mu':>14s} {'sigma':>14s} "
             f"{'N_sig':>8s} {'z':>5s} st"]
    for ds_cut, pt_cut in SELECTIONS:
        tag = f"ds{ds_cut:g}_pt{pt_cut:g}".replace(".", "p")
        sel_txt = f"score > {ds_cut:g}, p_{{T}} > {pt_cut:g} GeV"
        results = {}
        hists = {}
        for key, z, weighted in (("data", d, False), ("mc", m, True)):
            s = (z["ds"] > ds_cut) & (z["pt"] > pt_cut)
            w = np.ones(s.sum()) if not weighted else z["weight"][s]
            h = hist(f"h_{tag}_{key}", z["mass"][s], w)
            sb, b, res = fit(h, weighted)
            results[key], hists[key] = res, h
            notes = [sel_txt, "#mu = %.3f #pm %.3f, #sigma = %.3f GeV, %.1f#sigma"
                     % (res["mu"], res["dmu"], res["sigma"], res["z"])]
            draw_fit(f"eta_{tag}_{key}", h, sb, b, res, notes, lumi, key == "data", args.out)
            lines.append(f"{sel_txt.replace('_{T}', 'T'):22s} {key:5s} {res['n']:9.0f} "
                         f"{res['mu']:7.4f}+-{res['dmu']:.4f} {res['sigma']:7.4f}+-{res['dsigma']:.4f} "
                         f"{res['nsig']:8.0f} {res['z']:5.1f} {res['status']}")
        rd, rm = results["data"], results["mc"]
        lines.append(f"{'':22s} data/MC: mu {rd['mu'] / rm['mu']:.3f}  sigma {rd['sigma'] / rm['sigma']:.3f}")
        draw_shape(f"eta_{tag}_shape", hists["data"], hists["mc"], [sel_txt, "area normalized"],
                   lumi, args.out)

    os.makedirs(args.tables, exist_ok=True)
    txt = os.path.join(args.tables, "eta_control_fits.txt")
    with open(txt, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote", txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
