#!/usr/bin/env python3
"""v4 的 diphotonScore cut 掃描：收緊 classifier，峰能不能變得跟舊模型一樣準、還留多少 signal。

baseline 是 tagger 的 diphotonScore > 0.5（pass_allcuts_merged_ML，兩版相同）。對每個 m_a、
每個 cut 值算 v4 的
  f20   = 重建質量落在真值 ±20% 內的比例（峰的品質）
  N_win = N × f20（窗內 signal 數，同一批 2024 signal MC）
再除以舊模型（lepclean_old，baseline cut）的同一個量：
  f20 比 = 1 表示峰的品質追上舊模型；N_win 比 > 1 表示窗內 signal 比舊模型多。
只看 signal：背景（hadron fake）也隨 cut 變，結論要等 v4 背景產物算 S/sqrt(B)。

⚠️ v4 的 parquet 是 8 月的 tagger（沒有 lepton cleaning），舊模型是現行 tagger，
兩者的事件集不完全可比；lepclean 版 v4 產出後重跑本腳本。

PyROOT（LCG_104），版面走 hzg_style。
Input : <plots>/<tag>_mass.npz（extract_mlphoton_mass_versions.py）
Output: <plots>/score_scan_{f20,nwin}_ratio.pdf/.png、<tables>/score_scan.txt
"""
import argparse
import os
import sys

import numpy as np
import ROOT

if ROOT.gROOT.GetVersion().startswith("6.34"):
    sys.exit("ROOT 6.34 drops histograms from PDF output; use LCG_104")

HZG = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3"
sys.path.insert(0, f"{HZG}/plot/scripts")
from hzg_style import (  # noqa: E402
    LINE_WIDTH, MARKERS, PETROFF, canvas_single, cms_header, color,
    house_legend, panel_label, save, style_axes, _keep,
)

PLOTS = "/eos/home-p/pelai/HZa_merged/output_plots/MLPhoton_mass_versions"
TABLES = "/eos/home-p/pelai/HZa_merged/output_tables/MLPhoton_mass_versions"
TAGS = ["M0p1", "M0p2", "M0p3", "M0p4", "M0p5", "M0p6", "M0p7", "M0p8", "M0p9", "M1"]
CUTS = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99]
TITLE = "H #rightarrow Za, a #rightarrow #gamma#gamma"

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def m_of(tag):
    return float(tag[1:].replace("p", "."))


def f20(v, m):
    return float(np.mean(np.abs(v - m) < 0.2 * m)) if len(v) else 0.0


def sigma_eff(v, frac=0.683):
    if len(v) < 10:
        return float("nan")
    v = np.sort(v)
    k = int(np.ceil(frac * len(v)))
    return float((v[k - 1:] - v[:len(v) - k + 1]).min()) / 2.0


def graph_plot(name, curves, ytitle, notes, out):
    """curves: [(label, [(cut, y)])]; y from 0 up to where the highest point clears the legend"""
    c = canvas_single(name)
    leg = house_legend(c, [lab for lab, _ in curves], ncols=3, clear_label=TITLE,
                       clear_notes=notes, size=0.040)
    ymin, ymax_data = 0.0, max(y for _, p in curves for _, y in p)
    frame_h = 1.0 - c.GetTopMargin() - c.GetBottomMargin()
    below_leg = (leg.GetY1() - c.GetBottomMargin()) / frame_h - 0.04   # GetY1NDC is 0 until drawn
    ymax = ymin + (ymax_data - ymin) / below_leg
    frame = c.DrawFrame(CUTS[0] - 0.02, ymin, 1.0, ymax)
    style_axes(frame, c, "v4 diphotonScore cut", ytitle)
    one = ROOT.TLine(CUTS[0] - 0.02, 1.0, 1.0, 1.0)
    one.SetLineStyle(7)
    one.SetLineColor(ROOT.kGray + 2)
    one.SetLineWidth(2)
    one.Draw()
    _keep.extend([frame, one])
    for i, (lab, pts) in enumerate(curves):
        g = ROOT.TGraph(len(pts))
        for j, (x, y) in enumerate(pts):
            g.SetPoint(j, x, y)
        col = color(PETROFF[i % len(PETROFF)])
        g.SetLineColor(col)
        g.SetMarkerColor(col)
        g.SetLineWidth(LINE_WIDTH)
        g.SetMarkerStyle(MARKERS[i % len(MARKERS)])
        g.SetMarkerSize(1.2)
        g.Draw("LP SAME")
        leg.AddEntry(g, lab, "lp")
        _keep.append(g)
    leg.Draw()
    cms_header(c, None, "13.6")
    panel_label(c, TITLE, notes)
    c.RedrawAxis()
    save(c, os.path.join(out, name))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", nargs="+", default=TAGS)
    ap.add_argument("--plots", default=PLOTS)
    ap.add_argument("--tables", default=TABLES)
    ap.add_argument("--new", default="v4", help="npz key of the version scanned")
    ap.add_argument("--ref", default="old", help="npz key of the reference version")
    args = ap.parse_args()

    lines = [f"# {args.new} diphotonScore cut scan vs {args.ref} at its baseline cut (> 0.5); "
             f"f20 = fraction within +-20% of m_a, Nwin = N x f20; 2024 signal MC",
             f"{'tag':5s} {'cut':>5s} {'N':>7s} {'eff':>6s} {'sig_eff':>7s} {'f20':>6s} "
             f"{'Nwin':>7s} | {'N_ref':>7s} {'f20_ref':>7s} {'Nwin_ref':>8s} | "
             f"{'f20/ref':>7s} {'Nwin/ref':>8s}"]
    curves_f, curves_n = [], []
    for tag in args.tags:
        path = os.path.join(args.plots, f"{tag}_mass.npz")
        if not os.path.exists(path):
            continue
        z = np.load(path)
        if args.new not in z or args.ref not in z:
            print(f"{tag}: needs {args.new} and {args.ref}")
            continue
        m = m_of(tag)
        mv, sv = z[args.new], z[f"{args.new}_score"]
        ref = z[args.ref]
        f_ref, n_ref = f20(ref, m), len(ref)
        pf, pn = [], []
        for cut in CUTS:
            v = mv[sv > cut]
            f = f20(v, m)
            lines.append(f"{tag:5s} {cut:5.2f} {len(v):7d} {len(v) / len(mv):6.3f} "
                         f"{sigma_eff(v):7.4f} {f:6.3f} {len(v) * f:7.0f} | {n_ref:7d} "
                         f"{f_ref:7.3f} {n_ref * f_ref:8.0f} | {f / f_ref:7.3f} "
                         f"{len(v) * f / (n_ref * f_ref):8.3f}")
            pf.append((cut, f / f_ref))
            pn.append((cut, len(v) * f / (n_ref * f_ref)))
        lab = f"{m:g} GeV"
        curves_f.append((lab, pf))
        curves_n.append((lab, pn))

    os.makedirs(args.tables, exist_ok=True)
    txt = os.path.join(args.tables, "score_scan.txt")
    with open(txt, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote", txt)

    notes = ["v4 / EXO-22-022 model, 2024", "curves: m_{a}"]
    graph_plot("score_scan_f20_ratio", curves_f, "f(#pm20%) ratio", notes, args.plots)
    graph_plot("score_scan_nwin_ratio", curves_n, "N_{window} ratio", notes, args.plots)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
