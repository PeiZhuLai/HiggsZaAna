#!/usr/bin/env python3
"""
V2 study: overlay the trigger-efficiency curves obtained with the three denominators
(production = no requirement on the other leg; N2 = >= 2 selected same-flavor leptons;
OL = N2 + other leg above its double-lepton threshold) for one era, m_a = 1 GeV.
All three come from the same dedicated run, i.e. identical events.

Output: <out>/<era>/denomCompare_<Flavor>_<leg>_<trig>_mA_M1_<era>.pdf
"""
import argparse
import json
import os
import re
from array import array

import ROOT

CONSTANTS = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py"


def load_lumi():
    """Read the LUMI dict from the authoritative constants.py without importing it
    (its ERA dict uses set keys and cannot be exec'd)."""
    import ast
    txt = open(CONSTANTS).read()
    i = txt.index("LUMI = {")
    j = txt.index("{", i)
    depth = 0
    for k in range(j, len(txt)):
        depth += txt[k] == "{"
        depth -= txt[k] == "}"
        if depth == 0:
            break
    body = "\n".join(l.split("#")[0] for l in txt[j:k + 1].splitlines())
    return ast.literal_eval(body)


LUMI = load_lumi()


def graph(bins):
    xs, ys, el, eh = [], [], [], []
    for name, rec in bins.items():
        m = re.match(r"pt(\d+)to(\d+|Inf)", name)
        if not m or rec["in_bin"] <= 0:
            continue
        lo = float(m.group(1))
        hi = lo + 2.0 if m.group(2) == "Inf" else float(m.group(2))
        n, k = int(round(rec["in_bin"])), int(round(rec["pass_trigger"]))
        e = k / n
        xs.append(0.5 * (lo + hi))
        ys.append(100 * e)
        el.append(100 * (e - ROOT.TEfficiency.ClopperPearson(n, k, 0.682689492137086, False)))
        eh.append(100 * (ROOT.TEfficiency.ClopperPearson(n, k, 0.682689492137086, True) - e))
    z = [0.0] * len(xs)
    return ROOT.TGraphAsymmErrors(len(xs), array("d", xs), array("d", ys), array("d", z), array("d", z),
                                  array("d", el), array("d", eh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--era", default="2024")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    blk = json.load(open(os.path.join(a.in_dir, f"cutflow_Sig_MC_mA_M1_{a.era}.json")))["trigeff_nominal"]
    colors = {"": "#5790fc", "N2": "#f89c20", "OL": "#e42536"}
    marks = {"": 24, "N2": 26, "OL": 20}
    thr = {"ele": (25, 15), "mu": (20, 10)}
    fl_name = {"ele": ("Electron", "e"), "mu": ("Muon", "#mu")}
    for lep, trigs in (("ele", ("double_ele", "OR_ele")), ("mu", ("double_mu", "OR_mu"))):
        for leg in ("lead", "sublead"):
            for trig in trigs:
                c = ROOT.TCanvas("c", "", 800, 600)
                c.SetMargin(0.13, 0.04, 0.15, 0.08)
                c.SetTicks(1, 1)
                fr = c.DrawFrame(8, 0, 102, 125)
                FL, sym = fl_name[lep]
                fr.GetXaxis().SetTitle(f"{'Lead' if leg == 'lead' else 'Sublead'} {FL} p_{{T}} [GeV]")
                fr.GetYaxis().SetTitle("Trigger Efficiency (%)")
                for ax in (fr.GetXaxis(), fr.GetYaxis()):
                    ax.SetTitleSize(0.055)
                    ax.SetLabelSize(0.050)
                    ax.SetTitleOffset(1.2)
                other = "sublead" if leg == "lead" else "lead"
                othr = thr[lep][1] if leg == "lead" else thr[lep][0]
                labels = {
                    "": "no requirement on the other leg (current AN)",
                    "N2": f"N_{{{sym}}} #geq 2",
                    "OL": f"N_{{{sym}}} #geq 2, {other} {sym} p_{{T}} > {othr} GeV",
                }
                leg_box = ROOT.TLegend(0.17, 0.70, 0.93, 0.89)
                leg_box.SetBorderSize(0)
                leg_box.SetFillStyle(0)
                leg_box.SetTextFont(42)
                leg_box.SetTextSize(0.036)
                keep = []
                for suf in ("", "N2", "OL"):
                    key = f"trigeff_{lep}_{leg}{suf}_{trig}"
                    if key not in blk:
                        continue
                    g = graph(blk[key]["bins"])
                    col = ROOT.TColor.GetColor(colors[suf])
                    g.SetMarkerColor(col)
                    g.SetLineColor(col)
                    g.SetMarkerStyle(marks[suf])
                    g.SetMarkerSize(1.1)
                    g.Draw("P SAME")
                    leg_box.AddEntry(g, labels[suf], "lp")
                    keep.append(g)
                leg_box.Draw()
                line = ROOT.TLine(8, 100, 102, 100)
                line.SetLineStyle(2)
                line.Draw()
                lat = ROOT.TLatex()
                lat.SetNDC()
                lat.SetTextFont(42)
                lat.SetTextSize(0.045)
                lat.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Simulation}")
                lat.SetTextAlign(31)
                lat.SetTextSize(0.040)
                lat.DrawLatex(0.96, 0.93, f"{LUMI[a.era]:.2f} fb^{{-1}} (13.6 TeV)")
                lat.SetTextAlign(11)
                lat.SetTextSize(0.038)
                tlabel = ("Single- OR Double-" if trig.startswith("OR") else "Double-") + FL + " Trigger"
                lat.DrawLatex(0.52, 0.27, f"m_{{a}} = 1 GeV, {a.era}")
                lat.DrawLatex(0.52, 0.21, tlabel)
                os.makedirs(os.path.join(a.out, a.era), exist_ok=True)
                c.SaveAs(os.path.join(a.out, a.era, f"denomCompare_{FL}_{leg}_{trig}_mA_M1_{a.era}.pdf"))
                c.Close()


if __name__ == "__main__":
    main()
