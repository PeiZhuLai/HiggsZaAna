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

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/cutflowVmA"

# 你指定要畫的 cuts（分母皆為 all）
CUTS_TO_PLOT = ["trig_cut", "has_z_cand", "has_2g_cand", "event"]
# JSON 裡 cutflow 的 key（用這個最通用；若你想改 ele/mu 就換 key）
CUTFLOW_KEY = "zgammas_w"

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, '2024':108.95,
            'combined_run3':61.89 }
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

def _load_cutflow_point(path: str) -> Optional[Tuple[str,int,Dict[str,float]]]:
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
    cf = (data.get("cutflows") or {}).get(CUTFLOW_KEY) or {}
    allv = cf.get("all")
    if allv is None or float(allv) <= 0:
        return None
    effs: Dict[str, float] = {}
    for c in CUTS_TO_PLOT:
        v = cf.get(c)
        if v is None:
            continue
        effs[c] = float(v) / float(allv)
    # 至少要有一個 cut 才收
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
    g0.GetYaxis().SetTitle("Signal Efficiency")
    g0.GetXaxis().SetTitleOffset(1.1)
    g0.GetYaxis().SetTitleOffset(1.15)
    g0.GetXaxis().SetLimits(0, 31)
    g0.GetXaxis().SetTitleSize(0.055)
    g0.GetYaxis().SetTitleSize(0.055)
    g0.GetXaxis().SetLabelSize(0.05)
    g0.GetYaxis().SetLabelSize(0.05)
    g0.SetMinimum(0.0)
    g0.SetMaximum(0.47)
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
    out_png = out_dir / f"cutflowVmA_{CUTFLOW_KEY}_{year}.png"
    c1.SaveAs(str(out_png))
    c1.SaveAs(str(out_png.with_suffix(".pdf")))

def main():
    parser = argparse.ArgumentParser(description="Plot cutflow efficiency vs mA per year.")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022preEE). Use 'all' to plot all years.")
    args = parser.parse_args()

    in_dir = Path(baseDir)
    out_dir = Path(outDir)

    # year -> ma -> {cut:eff}
    store: Dict[str, Dict[int, Dict[str,float]]] = {}

    for p in sorted(in_dir.glob("cutflow_*.json")):
        got = _load_cutflow_point(str(p))
        if got is None:
            continue
        year, ma, effs = got
        store.setdefault(year, {})[ma] = effs

    # 依序輸出每個年份一張（或只輸出指定年份）
    if args.year != "all":
        if args.year in store:
            _plot_year(args.year, store[args.year], out_dir)
        return

    for year in YEAR_ORDER:
        if year in store:
            _plot_year(year, store[year], out_dir)

if __name__ == "__main__":
    main()
