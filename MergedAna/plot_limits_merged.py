#!/usr/bin/env python3
"""Plot the merged low-mA expected upper limit (in units of the xs=0.1 pb
reference signal strength r) vs m_a, with +/-1,2 sigma bands.  PyROOT, CMS
Preliminary.  Run in `higgs-alp-ana`.
"""
import argparse
import csv
import json
import os

import ROOT

ROOT.gROOT.SetBatch(True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "output"))
    ap.add_argument("--plotdir", default=os.path.join(os.path.dirname(__file__),
                                                      "../Plot/plots/merged_lowma"))
    ap.add_argument("--lumi", type=float, default=109.82)
    args = ap.parse_args()
    os.makedirs(args.plotdir, exist_ok=True)

    yields = {r["tag"]: r["ma"] for r in
              json.load(open(os.path.join(args.outdir, "signal_yields.json")))}
    csv_path = os.path.join(args.outdir, "limits", "limits.csv")
    rows = list(csv.DictReader(open(csv_path)))
    rows = [r for r in rows if r["tag"] in yields]
    rows.sort(key=lambda r: yields[r["tag"]])
    if not rows:
        raise SystemExit(f"[plot] no rows in {csv_path}")

    n = len(rows)
    g_med = ROOT.TGraph(n)
    g_1s = ROOT.TGraphAsymmErrors(n)
    g_2s = ROOT.TGraphAsymmErrors(n)
    def fnum(x, default):
        try:
            v = float(x)
            return v if v == v else default  # NaN -> default
        except (TypeError, ValueError):
            return default

    for i, r in enumerate(rows):
        ma = yields[r["tag"]]
        med = fnum(r["exp_med"], 0.0)
        m1, p1 = fnum(r["exp_m1"], med), fnum(r["exp_p1"], med)
        m2, p2 = fnum(r["exp_m2"], m1), fnum(r["exp_p2"], p1)
        g_med.SetPoint(i, ma, med)
        g_1s.SetPoint(i, ma, med)
        g_1s.SetPointError(i, 0, 0, med - m1, p1 - med)
        g_2s.SetPoint(i, ma, med)
        g_2s.SetPointError(i, 0, 0, med - m2, p2 - med)

    c = ROOT.TCanvas("c", "", 800, 700)
    c.SetLogy()
    c.SetLeftMargin(0.13)
    g_2s.SetFillColor(ROOT.kOrange)
    g_2s.SetTitle("")
    g_2s.GetXaxis().SetTitle("m_{a} [GeV]")
    g_2s.GetYaxis().SetTitle("95% CL exp. limit on r (#sigma_{ref}=0.1 pb)")
    g_2s.Draw("A3")
    g_1s.SetFillColor(ROOT.kGreen + 1)
    g_1s.Draw("3 SAME")
    g_med.SetLineStyle(2)
    g_med.SetLineWidth(2)
    g_med.Draw("L SAME")

    leg = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.AddEntry(g_med, "Median expected", "l")
    leg.AddEntry(g_1s, "68% expected", "f")
    leg.AddEntry(g_2s, "95% expected", "f")
    leg.Draw()

    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(42)
    l.SetTextSize(0.045)
    l.DrawLatex(0.13, 0.92, "#bf{CMS} #it{Preliminary}")
    l.SetTextSize(0.040)
    l.DrawLatex(0.62, 0.92, f"{args.lumi:.1f} fb^{{-1}} (13.6 TeV)")

    src_file = os.path.join(args.outdir, "bkg_source.txt")
    if os.path.exists(src_file) and open(src_file).read().strip() == "pseudo_mc":
        w = ROOT.TLatex()
        w.SetNDC()
        w.SetTextColor(ROOT.kRed + 1)
        w.SetTextFont(42)
        w.SetTextSize(0.032)
        w.DrawLatex(0.15, 0.18, "MC pseudo-data (closure only)")

    for ext in ("png", "pdf"):
        c.SaveAs(f"{args.plotdir}/limit_vs_mA_merged.{ext}")
    print(f"[plot] wrote {args.plotdir}/limit_vs_mA_merged.(png,pdf)")


if __name__ == "__main__":
    main()
