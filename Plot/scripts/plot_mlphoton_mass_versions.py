#!/usr/bin/env python3
"""重建 merged-photon 質量分布：舊模型 (EXO-22-022) / v3 / v4。

三種圖（--mode，可多選）：
  versions  每個 m_a 一張，三版各自的選擇疊圖（面積歸一）
  same      每個 m_a 一張，兩版都選到、且是同一顆 cluster（dR < 0.02）的事件，舊 vs v4
            —— 去掉 classifier 收的事件不同，只比 regressor
  all       每個版本一張，所有 m_a 疊在同一張圖

「舊模型」是 lepclean_old（現行分析鏈用的那份），見 extract_mlphoton_mass_versions.py。
真值 m_a 用虛線標出 —— 看峰有沒有往真值靠，而不只是變窄。

PyROOT，版面走共用的 hzg_style（800x600、房規字級、legend 左緣由字寬算）。
⚠️ 出 PDF 一律用 LCG_104(ROOT 6.28)：6.34 會把 histogram 丟掉只留座標軸
（rc=0 無警告，見 memory ref_root634_pdf_drops_histograms）。

Input : <out>/<tag>_mass.npz  (extract_mlphoton_mass_versions.py，hza_ana 抽出)
Output: <out>/<tag>_mass_versions、<tag>_mass_sameCluster、all_masses_<ver> 的 .pdf/.png，
        stdout 印 median / σ_eff / 窗內比例表
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
    HEADROOM, HEADROOM_MANY, LINE_WIDTH, PAIR_COLORS, PETROFF, canvas_single, cms_header,
    color, house_legend, panel_label, save, style_axes, _keep,
)

OUT = "/eos/home-p/pelai/HZa_merged/output_plots/MLPhoton_mass_versions"
TAGS = ["M0p1", "M0p2", "M0p3", "M0p4", "M0p5", "M0p6", "M0p7", "M0p8", "M0p9", "M1"]
VERSIONS = [
    ("old", "EXO-22-022 model", 1),
    ("v3",  "retrained v3",     2),
    ("v4",  "retrained v4",     1),
]
LABEL = dict((k, lab) for k, lab, _ in VERSIONS)
TITLE = "H #rightarrow Za, a #rightarrow #gamma#gamma"
XTITLE = "m_{#gamma#gamma}^{reco} (merged photon) [GeV]"
SAME_DR = 0.02          # same cluster: leading MLPhoton of both versions within this dR

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.TH1.AddDirectory(False)


def m_of(tag):
    return float(tag[1:].replace("p", "."))


def x_range(m_true):
    # 20 MeV bins up to 1.2 GeV; the heavier points need room for the old model's tail
    return (0.0, 1.2, 60) if m_true <= 0.4 else (0.0, 2.0, 50)


def sigma_eff(v, frac=0.683):
    v = np.sort(v)
    k = int(np.ceil(frac * len(v)))
    return float((v[k - 1:] - v[:len(v) - k + 1]).min()) / 2.0


def stats(v, m_true):
    q16, q50, q84 = np.quantile(v, [0.16, 0.5, 0.84])
    return len(v), q50, sigma_eff(v), float(np.mean(np.abs(v - m_true) < 0.2 * m_true))


def load(tag, out):
    path = os.path.join(out, f"{tag}_mass.npz")
    if not os.path.exists(path):
        print(f"{tag}: MISSING {path}")
        return None
    return np.load(path)


def hist(v, name, xlo, xhi, nbin, col, style):
    h = ROOT.TH1D(name, "", nbin, xlo, xhi)
    h.FillN(len(v), np.asarray(v, dtype=np.float64), np.ones(len(v)))
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())          # 面積歸一：比形狀不比 yield
    h.SetLineColor(col)
    h.SetLineWidth(LINE_WIDTH)
    h.SetLineStyle(style)
    _keep.append(h)
    return h


def draw_overlay(cname, series, xlo, xhi, nbin, notes, m_true=None, headroom=HEADROOM, ncols=1):
    """series: [(values, legend label, color, line style)]; returns the canvas."""
    c = canvas_single(cname)
    labels = [s[1] for s in series]
    leg = house_legend(c, labels, ncols=ncols, clear_label=TITLE, clear_notes=notes,
                       size=0.045 if ncols == 1 else 0.040)
    hists = []
    for i, (v, label, col, style) in enumerate(series):
        h = hist(v, f"h_{cname}_{i}", xlo, xhi, nbin, col, style)
        hists.append(h)
        leg.AddEntry(h, label, "l")
    hmax = max(h.GetMaximum() for h in hists)
    first = hists[0]
    style_axes(first, c, XTITLE, "A.U. / %.0f MeV" % (1000.0 * (xhi - xlo) / nbin))
    first.SetMinimum(0.0)
    first.SetMaximum(headroom * hmax)
    first.Draw("HIST")
    for h in hists[1:]:
        h.Draw("HIST SAME")
    if m_true is not None:
        # 真值：停在最高峰上方一點，不穿過 legend
        line = ROOT.TLine(m_true, 0.0, m_true, 1.05 * hmax)
        line.SetLineColor(ROOT.kGray + 2)
        line.SetLineStyle(7)
        line.SetLineWidth(LINE_WIDTH)
        line.Draw()
        _keep.append(line)
    leg.Draw()
    cms_header(c, None, "13.6")
    panel_label(c, TITLE, notes)
    c.RedrawAxis()
    return c


def plot_versions(tag, z, out, rows):
    m_true = m_of(tag)
    keys = [k for k, _, _ in VERSIONS if k in z]
    if not keys:
        return
    series = [(z[k], LABEL[k], color(PETROFF[i]), style)
              for i, (k, _, style) in enumerate(VERSIONS) if k in z]
    c = draw_overlay(f"c_{tag}", series, *x_range(m_true),
                     [f"m_{{a}} = {m_true:g} GeV, 2024", "dashed: true m_{a}"], m_true)
    save(c, os.path.join(out, f"{tag}_mass_versions"))
    rows += [(tag, "own selection", LABEL[k]) + stats(z[k], m_true) for k in keys]


def plot_same(tag, z, out, rows):
    m_true = m_of(tag)
    if "old" not in z or "v4" not in z:
        print(f"{tag}: same-cluster needs old and v4")
        return
    _, io, iv = np.intersect1d(z["old_key"], z["v4_key"], return_indices=True)
    dphi = np.abs((z["old_phi"][io] - z["v4_phi"][iv] + np.pi) % (2 * np.pi) - np.pi)
    same = np.hypot(z["old_eta"][io] - z["v4_eta"][iv], dphi) < SAME_DR
    a, b = z["old"][io][same], z["v4"][iv][same]
    series = [(a, LABEL["old"], PAIR_COLORS[0], 1), (b, LABEL["v4"], PAIR_COLORS[1], 1)]
    notes = [f"m_{{a}} = {m_true:g} GeV, 2024", f"same events + cluster: {same.sum()}"]
    c = draw_overlay(f"cs_{tag}", series, *x_range(m_true), notes, m_true)
    save(c, os.path.join(out, f"{tag}_mass_sameCluster"))
    rows += [(tag, "same cluster", LABEL["old"]) + stats(a, m_true),
             (tag, "same cluster", LABEL["v4"]) + stats(b, m_true)]


def plot_all(zs, out):
    for key, label, _ in VERSIONS:
        series = [(z[key], f"m_{{a}} = {m_of(tag):g} GeV", color(PETROFF[i % len(PETROFF)]), 1)
                  for i, (tag, z) in enumerate(zs) if key in z]
        if not series:
            continue
        c = draw_overlay(f"ca_{key}", series, 0.0, 1.5, 50, [label, "2024"],
                         headroom=HEADROOM_MANY, ncols=2)
        save(c, os.path.join(out, f"all_masses_{key}"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", nargs="+", default=TAGS)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--mode", nargs="+", default=["versions", "same", "all"],
                    choices=["versions", "same", "all"])
    args = ap.parse_args()
    zs = [(tag, z) for tag in args.tags for z in [load(tag, args.out)] if z is not None]
    rows = []
    for tag, z in zs:
        if "versions" in args.mode:
            plot_versions(tag, z, args.out, rows)
        if "same" in args.mode:
            plot_same(tag, z, args.out, rows)
    if "all" in args.mode:
        plot_all(zs, args.out)
    print()
    print(f"{'tag':5s} {'events':>13s} {'version':>17s} {'N':>7s} {'median':>7s} "
          f"{'sig_eff':>7s} {'f(+-20%)':>8s}")
    for r in rows:
        print("%-5s %13s %17s %7d %7.4f %7.4f %8.3f" % r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
