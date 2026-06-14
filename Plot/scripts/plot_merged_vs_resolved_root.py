#!/usr/bin/env python3
"""
PyROOT (re)plots for the merged-vs-resolved category study, replacing the
earlier matplotlib versions in Plot/plots/merged_vs_resolved/.

Produces:
  genDR_vs_mA.{png,pdf}            -- true dR(g1,g2) vs m_a, median + 10-90% band
  significance_vs_mA.{png,pdf}     -- Z_resolved / Z_merged / Z_combined
  significance_vs_mA_narrow.{png,pdf} -- + merged narrow-MergedML projection

Run in higgs-alp-ana (working ROOT).  CMS Preliminary / 109.82 fb^-1 (13.6 TeV).
"""
from __future__ import annotations

import array
import glob
import json
import os
import sys

import numpy as np
import pyarrow.parquet as pq
import ROOT

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from _cms_style import make_legend  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

EOS = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"
OUTDIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/merged_vs_resolved"
JDIR = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output"
LUMI = 109.82
DR_MERGE, DR_RESOLVE = 0.07, 0.30


def _a(xs):
    return array.array("d", [float(v) for v in xs])


def _save(c, stub):
    for ext in ("png", "pdf"):
        c.SaveAs(f"{OUTDIR}/{stub}.{ext}")
        print("saved", f"{OUTDIR}/{stub}.{ext}")


def _header(label="Preliminary"):
    # CMS-standard label: bold "CMS" + italic suffix (font 42 + #bf/#it inline),
    # matching plot_MVASigEffVmA.py. Returns TLatex to keep it alive.
    lat = ROOT.TLatex()
    lat.SetNDC(True)
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.SetTextAlign(11)
    lat.DrawLatex(0.13, 0.93, "#bf{CMS} #it{%s}" % label)
    lat.SetTextAlign(31)
    lat.DrawLatex(0.95, 0.93, "%.1f fb^{-1} (13.6 TeV)" % LUMI)
    return lat


# ----------------------------------------------------------------- gen dR ----
RES_BASE = f"{EOS}/parquet_DNA/Sig_MC"
MRG_BASES = [f"{EOS}/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all",
             f"{EOS}/parquet_merged_DNA_tmp/Sig_MC_MLNANO",
             f"{EOS}/parquet_merged_DNA_tmp/Sig_MC_MLNANO_M1"]
POINTS = [(f"M0p{i}", i / 10.0, "mrg") for i in range(1, 10)] + \
         [(f"M{m}", float(m), "res") for m in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]]


def _path(tag, src):
    if src == "res":
        p = f"{RES_BASE}/mA_{tag}_2024/merged_nominal.parquet"
        return p if glob.glob(p) else None
    for b in MRG_BASES:
        p = f"{b}/mA_MLNANO_{tag}_2024/merged_nominal.parquet"
        if glob.glob(p):
            return p
    return None


def _gendr(path):
    t = pq.read_table(path, columns=["GenALPPhoton1_eta", "GenALPPhoton1_phi",
                                     "GenALPPhoton2_eta", "GenALPPhoton2_phi"]).to_pandas()
    e1, p1 = t.GenALPPhoton1_eta.to_numpy(), t.GenALPPhoton1_phi.to_numpy()
    e2, p2 = t.GenALPPhoton2_eta.to_numpy(), t.GenALPPhoton2_phi.to_numpy()
    acc = (np.abs(e1) < 2.5) & (np.abs(e2) < 2.5)
    dp = np.arctan2(np.sin(p1 - p2), np.cos(p1 - p2))
    dr = np.sqrt((e1 - e2) ** 2 + dp ** 2)
    return dr[acc & np.isfinite(dr)]


def plot_gendr():
    ma, med, p10, p90 = [], [], [], []
    for tag, m, src in POINTS:
        p = _path(tag, src)
        if not p:
            continue
        dr = _gendr(p)
        if len(dr) == 0:
            continue
        ma.append(m); med.append(np.median(dr))
        p10.append(np.percentile(dr, 10)); p90.append(np.percentile(dr, 90))

    c = ROOT.TCanvas("c_gendr", "", 900, 700)
    c.SetLogx(); c.SetLogy()
    c.SetLeftMargin(0.13); c.SetRightMargin(0.05); c.SetTopMargin(0.08); c.SetBottomMargin(0.13)

    frame = ROOT.TH2F("fr_gendr", ";m_{a} [GeV];gen #DeltaR(#gamma_{1},#gamma_{2})  (|#eta|<2.5)",
                      100, 0.08, 40, 100, 3e-3, 5)
    frame.GetXaxis().SetTitleSize(0.05); frame.GetXaxis().SetLabelSize(0.045); frame.GetXaxis().SetTitleOffset(1.15)
    frame.GetYaxis().SetTitleSize(0.05); frame.GetYaxis().SetLabelSize(0.045); frame.GetYaxis().SetTitleOffset(1.20)
    frame.Draw()

    # regime shading
    boxes = []
    for y1, y2, col in [(3e-3, DR_MERGE, ROOT.kRed - 10),
                        (DR_MERGE, DR_RESOLVE, ROOT.kOrange - 9),
                        (DR_RESOLVE, 5, ROOT.kGreen - 10)]:
        b = ROOT.TBox(0.08, y1, 40, y2); b.SetFillColorAlpha(col, 0.35); b.SetLineWidth(0)
        b.Draw("same"); boxes.append(b)
    ROOT.gPad.RedrawAxis()

    lines = []
    for yv, col in [(DR_MERGE, ROOT.kRed + 1), (DR_RESOLVE, ROOT.kGreen + 2)]:
        ln = ROOT.TLine(0.08, yv, 40, yv); ln.SetLineColor(col); ln.SetLineStyle(2); ln.SetLineWidth(2)
        ln.Draw("same"); lines.append(ln)

    n = len(ma)
    band = ROOT.TGraph(2 * n)
    for i in range(n):
        band.SetPoint(i, ma[i], p90[i])
        band.SetPoint(2 * n - 1 - i, ma[i], p10[i])
    band.SetFillColorAlpha(ROOT.kAzure + 1, 0.30); band.SetLineWidth(0)
    band.Draw("F same")

    g = ROOT.TGraph(n, _a(ma), _a(med))
    g.SetLineColor(ROOT.kBlue + 2); g.SetLineWidth(3)
    g.SetMarkerColor(ROOT.kBlue + 2); g.SetMarkerStyle(20); g.SetMarkerSize(1.1)
    g.Draw("LP same")

    # regime text
    txts = []
    for x, y, s, col in [(0.10, 0.03, "MERGED", ROOT.kRed + 2),
                         (0.10, 0.14, "transition (both)", ROOT.kOrange + 3),
                         (0.10, 2.2, "RESOLVED", ROOT.kGreen + 3)]:
        t = ROOT.TLatex(x, y, s); t.SetTextSize(0.035); t.SetTextColor(col); t.Draw("same"); txts.append(t)

    leg = make_legend(0.55, 0.18, 0.93, 0.30)
    leg.AddEntry(g, "median #DeltaR(#gamma#gamma)", "lp")
    leg.AddEntry(band, "10-90% of #DeltaR(#gamma#gamma)", "f")
    leg.Draw()
    hdr = _header()
    c.Update()
    _save(c, "genDR_vs_mA")
    return c, frame, boxes, lines, band, g, txts, leg, hdr


# -------------------------------------------------------------- significance --
def _series(rows, key, xkey="ma"):
    xs, ys = [], []
    for r in rows:
        if r.get(key) is not None:
            xs.append(r[xkey]); ys.append(r[key])
    return xs, ys


def _graph(xs, ys, color, mstyle, lstyle=1):
    g = ROOT.TGraph(len(xs), _a(xs), _a(ys))
    g.SetLineColor(color); g.SetLineWidth(3); g.SetLineStyle(lstyle)
    g.SetMarkerColor(color); g.SetMarkerStyle(mstyle); g.SetMarkerSize(1.2)
    return g


def plot_significance(narrow=False, logy=True):
    d = json.load(open(f"{JDIR}/significance_merged_vs_resolved.json"))
    xr, yr = _series(d, "Z_resolved")
    xm, ym = _series(d, "Z_merged")
    xc, yc = _series(d, "Z_combined")

    suf = ("n" if narrow else "") + ("" if logy else "L")
    ymin, ymax = (0.1, 60) if logy else (0.0, 32.0)
    c = ROOT.TCanvas("c_sig" + suf, "", 900, 700)
    c.SetLogx()
    if logy:
        c.SetLogy()
    c.SetLeftMargin(0.13); c.SetRightMargin(0.05); c.SetTopMargin(0.08); c.SetBottomMargin(0.13)

    frame = ROOT.TH2F("fr_sig" + suf,
                      ";m_{a} [GeV];Expected Significance",
                      100, 0.08, 40, 100, ymin, ymax)
    frame.GetXaxis().SetTitleSize(0.05); frame.GetXaxis().SetLabelSize(0.045); frame.GetXaxis().SetTitleOffset(1.15)
    frame.GetYaxis().SetTitleSize(0.05); frame.GetYaxis().SetLabelSize(0.045); frame.GetYaxis().SetTitleOffset(1.20)
    frame.Draw()

    gm = _graph(xm, ym, ROOT.kRed + 1, 21)
    gr = _graph(xr, yr, ROOT.kAzure + 2, 20)
    gc = _graph(xc, yc, ROOT.kGreen + 2, 22, lstyle=2)
    gm.Draw("LP same"); gr.Draw("LP same"); gc.Draw("LP same")

    objs = [gm, gr, gc]
    # log: lower-right is empty; linear: upper-left is empty.
    leg = (make_legend(0.50, 0.16, 0.93, 0.34, text_size=0.034) if logy
           else make_legend(0.16, 0.63, 0.60, 0.88, text_size=0.034))
    leg.AddEntry(gm, "Merged (m_{H} window only)", "lp")
    gn = None
    if narrow:
        pj = json.load(open(f"{JDIR}/significance_merged_narrow_projection.json"))
        xn, yn = _series(pj, "Z_narrow")
        gn = _graph(xn, yn, ROOT.kRed + 3, 5, lstyle=2)
        gn.Draw("LP same"); objs.append(gn)
        leg.AddEntry(gn, "Merged + MergedML m_{a} window", "lp")
    leg.AddEntry(gr, "Resolved (m_{a} peak, 2#gamma)", "lp")
    leg.AddEntry(gc, "Resolved + Merged (quad.)", "lp")
    leg.Draw()

    lv = ROOT.TLine(1.0, ymin, 1.0, ymax); lv.SetLineStyle(3); lv.SetLineColor(ROOT.kGray + 2); lv.Draw("same")
    hdr = _header()
    c.Update()
    stub = "significance_vs_mA_narrow" if narrow else "significance_vs_mA"
    if not logy:
        stub += "_liny"
    _save(c, stub)
    return c, frame, objs, leg, lv, hdr


if __name__ == "__main__":
    keep = []
    keep.append(plot_gendr())
    keep.append(plot_significance(narrow=False, logy=True))
    keep.append(plot_significance(narrow=True, logy=True))
    keep.append(plot_significance(narrow=False, logy=False))
    keep.append(plot_significance(narrow=True, logy=False))
    print("done")
