# ==== New imports/config for plotting ====
import os
from pathlib import Path
# 先安全匯入 PyROOT，避免 Cling 解析到重複 LLVM 參數
import sys
import ROOT
import numpy as np
import re
import ast
from typing import Dict, List, Tuple, Optional
import uproot
import awkward as ak
# 加回：SciPy 嘗試導入 (供 _quadratic_curve 使用)
try:
    from scipy.interpolate import make_interp_spline
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False
import json  # 新增：讀取 MVAcut 的 JSON
import argparse  # 新增：CLI 參數

baseDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/cutflow/cutflow_list"

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/preselectSigEffVmA"

# 你指定要畫的 cuts（分母皆為 all）
CUTS_TO_PLOT = ["event"]
# JSON 裡 cutflow 的 key（分 channel）
CUTFLOW_KEY_ELE = "zgammas_ele_w"
CUTFLOW_KEY_MU  = "zgammas_mu_w"

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

palette_hex_mu = [
    "#540D6E", "#EE4266", "#FFB640",
    "#3BCEAC", "#086788",
]

palette_hex_ele = [
    "#5C7AFF", "#242038", "#59D2FE",
    "#44E5E7", "#417B5A",
]

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, '2024':108.95,
            'combined_run3':170.84 }
# 移除：YEAR 未定義時會直接噴錯；改為在 _plot_year() 依 year 動態取用
# LUMI_FB = float(lumiMap[YEAR])

def _lumi_fb_for_year(year: str) -> Optional[float]:
    v = lumiMap.get(year)
    return float(v) if v is not None else None

def _parse_year_from_filename(fn: str) -> Optional[str]:
    # 例如：cutflow_Sig_MC_mA_M1_2022preEE.json
    for y in YEAR_ORDER:
        if y in fn:
            return y
    return None

def _parse_ma_from_filename(fn: str) -> Optional[int]:
    # 例如：...mA_M15_...
    m = re.search(r"mA_M(\d+)", fn)
    return int(m.group(1)) if m else None

def _load_cutflow_point(path: str, channel_key: str) -> Optional[Tuple[str,int,Dict[str,float]]]:
    fn = os.path.basename(path)
    year = _parse_year_from_filename(fn)
    ma = _parse_ma_from_filename(fn)
    if year is None or ma is None:
        return None
    try:
        with open(path, "r") as f:
            data = json.load(f)
    except Exception:
        return None
    cf = (data.get("cutflows") or {}).get(channel_key) or {}
    allv = cf.get("all")
    if allv is None or float(allv) <= 0:
        return None
    effs: Dict[str, float] = {}
    for c in CUTS_TO_PLOT:
        v = cf.get(c)
        if v is None:
            continue
        # 直接存成百分比
        effs[c] = (float(v) / float(allv)) * 100.0
    if not effs:
        return None
    return year, ma, effs

def _g_from_xy(xs: List[int], ys: List[float]) -> ROOT.TGraph:
    from array import array as carray
    x = carray("d", [float(v) for v in xs])
    y = carray("d", [float(v) for v in ys])
    return ROOT.TGraph(len(xs), x, y)

def _hex_to_rgb01(hexstr: str) -> Tuple[int, int, int]:
    s = hexstr.strip()
    if s.startswith("#"):
        s = s[1:]
    if len(s) != 6:
        raise ValueError(f"Invalid hex color: {hexstr}")
    r = int(s[0:2], 16)
    g = int(s[2:4], 16)
    b = int(s[4:6], 16)
    return r, g, b

def _root_color(hexstr: str, *, fallback: int = 1) -> int:
    """
    ROOT 有些環境/版本對 TColor::GetColor('#RRGGBB') 支援不完整；
    若取不到顏色（常見返回 0/負值），改用 (r,g,b) 介面。
    """
    try:
        ci = int(ROOT.TColor.GetColor(hexstr))
    except Exception:
        ci = 0
    if ci > 0:
        return ci
    try:
        r, g, b = _hex_to_rgb01(hexstr)
        ci2 = int(ROOT.TColor.GetColor(r, g, b))
        return ci2 if ci2 > 0 else fallback
    except Exception:
        return fallback

def _plot_year(year: str, points: Dict[int, Dict[str,float]], out_dir: Path) -> None:
    mas = sorted(points.keys())
    if not mas:
        return

    # 每個 cut 一條線
    graphs: Dict[str, ROOT.TGraph] = {}
    for c in CUTS_TO_PLOT:
        xs, ys = [], []
        for ma in mas:
            if c in points[ma]:
                xs.append(ma)
                ys.append(points[ma][c])
        if len(xs) >= 1:
            graphs[c] = _g_from_xy(xs, ys)

    if not graphs:
        return

    # 樣式
    color_map = {
        "trig_cut": _root_color("#3185FC"),
        "has_z_cand": _root_color("#FF8000"),
        "has_2g_cand": _root_color("#33BA91"),
        "event": _root_color("#E3170A"),
    }

    marker_map = {"trig_cut": 20, "has_z_cand": 22, "has_2g_cand": 21, "event": 33}
    marker_size_map = {
        "trig_cut": 1.4,      # 第一
        "has_z_cand": 1.6,    # 第二
        "has_2g_cand": 1.4,   # 第三（改成不同）
        "event": 2.0,         # 第四（改成不同）
    }

    legend_map = {
        "trig_cut": "Trigger",
        "has_z_cand": "Trigger + Z Candidate",
        "has_2g_cand": "Trigger + Z Candidate + a Candidate",
        "event": "Trigger + Z Candidate + a Candidate + Za Candidate",
    }

    c1 = ROOT.TCanvas(f"c_cutflow_{year}", "", 800, 600)
    c1.SetMargin(0.13, 0.04, 0.13, 0.08)
    c1.SetTickx()
    c1.SetTicky()

    # frame：用第一條 graph 建軸
    first_key = next(iter(graphs.keys()))
    g0 = graphs[first_key]
    g0.SetTitle("")
    g0.GetXaxis().SetTitle("m_{a} [GeV]")
    g0.GetYaxis().SetTitle("Preselection Signal Efficiency (%)")
    g0.GetXaxis().SetTitleOffset(1.1)
    g0.GetYaxis().SetTitleOffset(1.15)
    g0.GetXaxis().SetLimits(0, 31)
    g0.GetXaxis().SetTitleSize(0.055)
    g0.GetYaxis().SetTitleSize(0.055)
    g0.GetXaxis().SetLabelSize(0.05)
    g0.GetYaxis().SetLabelSize(0.05)
    g0.SetMinimum(0.0)
    g0.SetMaximum(47)
    g0.SetLineColor(color_map.get(first_key, 1))
    g0.SetMarkerColor(color_map.get(first_key, 1))
    g0.SetMarkerStyle(marker_map.get(first_key, 20))
    g0.SetMarkerSize(marker_size_map.get(first_key, 1.4))
    g0.SetLineWidth(3)
    g0.Draw("ALPC")

    for k, g in graphs.items():
        if k == first_key:
            continue
        col = color_map.get(k, 1)
        g.SetLineColor(col)
        g.SetMarkerColor(col)
        g.SetMarkerStyle(marker_map.get(k, 20))
        g.SetMarkerSize(marker_size_map.get(k, 1.4))
        g.SetLineWidth(3)
        g.Draw("LPC SAME")

    leg = ROOT.TLegend(0.16, 0.43, 0.41, 0.68)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.AddEntry(0, year, "")   # 左對齊 header
    for k in CUTS_TO_PLOT:
        if k in graphs:
            label = legend_map.get(k, k)
            leg.AddEntry(graphs[k], label, "lp")
    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Preliminary}")

    # 右上角 lumi（每一年各自顯示）
    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)  # right-aligned
        lat_lumi.SetTextSize(0.040)
        lat_lumi.DrawLatex(0.96, 0.93, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = out_dir / f"preselectSigEffVmA_{CUTFLOW_KEY}_{year}.pdf"
    c1.SaveAs(str(out_pdf))

def _plot_channel(channel: str, store: Dict[str, Dict[int, Dict[str, float]]], palette_hex: List[str], out_dir: Path) -> None:
    # 只畫 CUTS_TO_PLOT 的第一個（你目前是 ["event"]）
    cut = CUTS_TO_PLOT[0]

    # 準備每個 year 的 TGraph
    graphs: Dict[str, ROOT.TGraph] = {}
    for year in YEAR_ORDER:
        points = store.get(year)
        if not points:
            continue
        mas = sorted(points.keys())
        xs, ys = [], []
        for ma in mas:
            v = points[ma].get(cut)
            if v is None:
                continue
            xs.append(ma)
            ys.append(v)
        if len(xs) >= 1:
            graphs[year] = _g_from_xy(xs, ys)

    if not graphs:
        return

    # ---- user requested styles ----
    line_syles = [1, 2, 3]
    marker_styles = [20, 21, 23, 33, 34, 47, 29]
    marker_size = [1.4, 1.3, 1.6, 1.9, 1.7, 1.6, 1.9]

    # canvas / axes：用 YEAR_ORDER 第一個存在的 year 當 frame
    first_year = next((y for y in YEAR_ORDER if y in graphs), None)
    if first_year is None:
        return
    g0 = graphs[first_year]

    c1 = ROOT.TCanvas(f"c_{channel}_{cut}", "", 800, 600)
    c1.SetMargin(0.13, 0.04, 0.13, 0.08)
    c1.SetTickx()
    c1.SetTicky()

    g0.SetTitle("")
    g0.GetXaxis().SetTitle("m_{a} [GeV]")
    g0.GetYaxis().SetTitle(f"Preselection Signal Efficiency (%)")
    g0.GetXaxis().SetTitleOffset(1.1)
    g0.GetYaxis().SetTitleOffset(1.15)
    g0.GetXaxis().SetLimits(0, 31)
    g0.GetXaxis().SetTitleSize(0.055)
    g0.GetYaxis().SetTitleSize(0.055)
    g0.GetXaxis().SetLabelSize(0.05)
    g0.GetYaxis().SetLabelSize(0.05)
    g0.SetMinimum(1.2)
    g0.SetMaximum(8)

    # 依 YEAR_ORDER 上色（palette 對應 YEAR_ORDER index）
    def _col_for_year(y: str) -> int:
        idx = YEAR_ORDER.index(y) if y in YEAR_ORDER else 0
        hexstr = palette_hex[idx % len(palette_hex)]
        return _root_color(hexstr)

    # 依 years_present index 指定 style（line/marker/size）
    years_present = [y for y in YEAR_ORDER if y in graphs]

    def _apply_style(g: ROOT.TGraph, y: str) -> None:
        i = years_present.index(y) if y in years_present else 0
        col = _col_for_year(y)
        g.SetLineColor(col)
        g.SetMarkerColor(col)
        g.SetLineStyle(line_syles[i % len(line_syles)])
        g.SetMarkerStyle(marker_styles[i % len(marker_styles)])
        g.SetMarkerSize(marker_size[i % len(marker_size)])
        g.SetLineWidth(3)

    _apply_style(g0, first_year)
    g0.Draw("ALPC")

    for year in YEAR_ORDER:
        g = graphs.get(year)
        if not g or year == first_year:
            continue
        _apply_style(g, year)
        g.Draw("LPC SAME")

    # --- legend (channel-dependent header + dynamic geometry at top-right) ---
    channel_label = "Muon" if channel.lower().startswith("mu") else "Electron"
    n_lines = 1 + len(years_present)  # header + years

    # geometry: keep same top-right anchor, vary height with lines
    x1, y1 = 0.95, 0.85
    x0 = 0.63
    line_h = 0.052
    pad = 0.018
    h_leg = pad + n_lines * line_h
    y0 = max(0.14, y1 - h_leg)

    leg = ROOT.TLegend(x0, y0, x1, y1)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)

    if channel_label == "Muon":
        header = "H #rightarrow Za #rightarrow #mu^{+}#mu^{-} + 2#gamma"
    else:
        header = "H #rightarrow Za #rightarrow e^{+}e^{-} + 2#gamma"
    leg.SetHeader(header, "L")

    for year in years_present:
        leg.AddEntry(graphs[year], year, "lp")
    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Simulation}")
    x_lumi = 1.0 - c1.GetRightMargin() - 0.005
    y_lumi = 0.92 + 0.01
    lat.SetTextAlign(31)
    lat.DrawLatexNDC(x_lumi, y_lumi, "13.6 TeV")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = out_dir / f"preselectSigEffVmA_{channel}_{cut}.pdf"
    c1.SaveAs(str(out_pdf))

def main():
    parser = argparse.ArgumentParser(description="Plot cutflow efficiency vs mA.")
    args = parser.parse_args()

    in_dir = Path(baseDir)
    out_dir = Path(outDir)

    # channel store: year -> ma -> {cut:eff}
    store_ele: Dict[str, Dict[int, Dict[str, float]]] = {}
    store_mu: Dict[str, Dict[int, Dict[str, float]]] = {}

    for p in sorted(in_dir.glob("cutflow_*.json")):
        got_ele = _load_cutflow_point(str(p), CUTFLOW_KEY_ELE)
        if got_ele is not None:
            year, ma, effs = got_ele
            store_ele.setdefault(year, {})[ma] = effs

        got_mu = _load_cutflow_point(str(p), CUTFLOW_KEY_MU)
        if got_mu is not None:
            year, ma, effs = got_mu
            store_mu.setdefault(year, {})[ma] = effs

    _plot_channel("electron", store_ele, palette_hex_ele, out_dir)
    _plot_channel("muon", store_mu, palette_hex_mu, out_dir)

if __name__ == "__main__":
    main()
