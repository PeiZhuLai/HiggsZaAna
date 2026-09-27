#!/usr/bin/env python3
"""M0p1 的重建 merged-photon 質量分布：舊模型 / v3 / v4 疊圖。

回答的問題：ROI 窗變窄，是「解析度真的變好」還是「先驗收縮」？
真值 m_a = 0.1 GeV 用虛線標出 —— 看峰有沒有往真值靠，而不只是變窄。

PyROOT。⚠️ 出 PDF 一律用 LCG_104(ROOT 6.28)：6.34 會把 histogram 丟掉只留座標軸
（rc=0 無警告，見 memory ref_root634_pdf_drops_histograms）。

Input : <out>/m0p1_mass.npz  (由 hza_ana 抽出，避開 LCG 沒有 pyarrow 的問題)
Output: <out>/m0p1_mass_versions.pdf / .png
"""
import os
import numpy as np
import ROOT

OUT = "/eos/home-p/pelai/HZa/output_plots/MLPhoton_M0p1"
M_TRUE = 0.1
XLO, XHI, NBIN = 0.0, 1.2, 60
LEFT_MARGIN = 0.215

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

z = np.load(os.path.join(OUT, "m0p1_mass.npz"))
VERSIONS = [
    ("old", "EXO-22-022 model", ROOT.kBlack, 1),
    ("v3",  "retrained v3",     ROOT.kAzure + 2, 2),
    ("v4",  "retrained v4",     ROOT.kRed + 1, 1),
]

c = ROOT.TCanvas("c", "c", 900, 750)
c.SetLeftMargin(LEFT_MARGIN)
c.SetRightMargin(0.05)
c.SetBottomMargin(0.14)
c.SetTopMargin(0.085)
c.SetTickx(1)
c.SetTicky(1)

hists = []
hmax = 0.0
for key, label, color, style in VERSIONS:
    if key not in z:
        continue
    v = z[key]
    h = ROOT.TH1D(f"h_{key}", "", NBIN, XLO, XHI)
    for x in v:
        h.Fill(float(x))
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())          # 面積歸一：比形狀不比 yield
    h.SetLineColor(color)
    h.SetLineWidth(3)                        # 一律 3
    h.SetLineStyle(style)
    hmax = max(hmax, h.GetMaximum())
    hists.append((h, label, v))

first = True
for h, _, _ in hists:
    h.SetMaximum(hmax * 1.45)
    h.GetXaxis().SetTitle("m_{#gamma#gamma}^{reco} (merged photon) [GeV]")
    h.GetYaxis().SetTitle("A.U. / %.0f MeV" % (1000.0 * (XHI - XLO) / NBIN))
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetTitleSize(0.055)
        ax.SetLabelSize(0.050)
    h.GetYaxis().SetTitleOffset(1.55)
    h.GetXaxis().SetTitleOffset(1.15)
    h.Draw("HIST" if first else "HISTSAME")
    first = False

# 真值：所有版本都應該往這裡收斂
line = ROOT.TLine(M_TRUE, 0.0, M_TRUE, hmax * 1.45 * 0.58)   # 停在 legend 底邊以下
line.SetLineColor(ROOT.kGray + 2)
line.SetLineStyle(7)
line.SetLineWidth(3)
line.Draw()

lat = ROOT.TLatex()
lat.SetNDC(True)
lat.SetTextFont(42)
lat.SetTextAlign(13)
lat.SetTextSize(0.050)
lat.DrawLatex(LEFT_MARGIN + 0.005, 0.935, "#bf{CMS} #it{Simulation}")
lat.SetTextAlign(31)
lat.SetTextSize(0.040)
lat.DrawLatex(0.95, 0.935, "13.6 TeV")

lat.SetTextAlign(13)
# Split across lines: as one line this runs past x=0.58 and collides with the
# legend (the legend x range is fixed by the house style).
lat.SetTextSize(0.038)
lat.DrawLatex(0.245, 0.80, "H #rightarrow Za, a #rightarrow #gamma#gamma")
lat.DrawLatex(0.245, 0.745, "m_{a} = 0.1 GeV, 2024")
lat.SetTextSize(0.034)
lat.SetTextColor(ROOT.kGray + 3)
lat.DrawLatex(0.245, 0.692, "dashed: true m_{a}")
lat.SetTextColor(ROOT.kBlack)

leg = ROOT.TLegend(0.58, 0.68, 0.93, 0.68 + 0.075 * len(hists))
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.045)
for h, label, v in hists:
    q16, q50, q84 = np.quantile(v, [0.16, 0.5, 0.84])
    leg.AddEntry(h, f"{label}", "l")
leg.Draw()

c.RedrawAxis()
os.makedirs(OUT, exist_ok=True)
pdf = os.path.join(OUT, "m0p1_mass_versions.pdf")
png = os.path.join(OUT, "m0p1_mass_versions.png")
c.SaveAs(pdf)
c.SaveAs(png)
print(f"wrote {pdf}")
print(f"wrote {png}")

print()
print(f"{'version':>22s} {'N':>7s} {'median':>8s} {'q16':>7s} {'q84':>7s} {'width':>7s} {'|med-0.1|':>10s}")
for h, label, v in hists:
    q16, q50, q84 = np.quantile(v, [0.16, 0.5, 0.84])
    print(f"{label:>22s} {len(v):>7d} {q50:>8.4f} {q16:>7.4f} {q84:>7.4f} "
          f"{q84-q16:>7.4f} {abs(q50-M_TRUE):>10.4f}")
