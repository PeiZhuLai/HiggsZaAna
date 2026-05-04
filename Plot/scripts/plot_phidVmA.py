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
import ctypes  # 新增：取代 ROOT.Double，避免某些 PyROOT 環境沒有 ROOT.Double

# NEW: globally disable ROOT stat/fit boxes (Entries/Mean/Std Dev)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

baseDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list"

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/phidVmA"

# 你指定要畫的 cuts（分母皆為 all）
CUTS_TO_PLOT = ["event"]

# JSON 裡 cutflow 的 key：改成 15 種 photon ID scenarios
# _PHID_WPS = ["tight", "medium", "loose"]
# FIX: 固定順序（tight/medium/loose），避免不同地方迭代順序不一致造成誤判
_PHID_WPS = ["loose", "medium", "tight"]
_PHID_KINDS = ["custom", "custom_extend", "sieie", "PFECalIso", "official"]
PHID_CUTFLOW_TYPES = [f"zgammas_phid_{kind}_{wp}" for wp in _PHID_WPS for kind in _PHID_KINDS]

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

def _load_cutflow_point(path: str) -> Optional[Tuple[str, int, Dict[str, Dict[str, float]]]]:
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

    cutflows = data.get("cutflows") or {}

    # scenario -> {cut: eff}
    effs_by_scenario: Dict[str, Dict[str, float]] = {}
    for scenario in PHID_CUTFLOW_TYPES:
        cf = cutflows.get(scenario) or {}
        allv = cf.get("all")
        if allv is None or float(allv) <= 0:
            continue

        effs: Dict[str, float] = {}
        for c in CUTS_TO_PLOT:
            v = cf.get(c)
            if v is None and c == "event":
                v = cf.get("all cuts")
            if v is None:
                continue
            effs[c] = float(v) / float(allv)

        if effs:
            effs_by_scenario[scenario] = effs

    if not effs_by_scenario:
        return None
    return year, ma, effs_by_scenario

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

def _plot_year(year: str, points: Dict[int, Dict[str, Dict[str, float]]], out_dir: Path) -> None:
    mas = sorted(points.keys())
    if not mas:
        return

    # scenario -> graph(event eff vs mA)
    scenario_graphs: Dict[str, ROOT.TGraph] = {}
    for scenario in PHID_CUTFLOW_TYPES:
        xs, ys = [], []
        for ma in mas:
            per_ma = points.get(ma) or {}
            effs = per_ma.get(scenario) or {}
            if "event" in effs:
                xs.append(ma)
                ys.append(effs["event"])
        if len(xs) >= 1:
            scenario_graphs[scenario] = _g_from_xy(xs, ys)

    if not scenario_graphs:
        return

    # kind 用顏色
    kind_colors = {
        "custom": _root_color("#3185FC"),
        "custom_extend": _root_color("#33BA91"),
        "sieie": _root_color("#FF8000"),
        "PFECalIso": _root_color("#7D4E57"),
        "official": _root_color("#E3170A"),
    }

    legend_map = {
        "custom": "H/E + CHIso + HIso (Custom)",
        "custom_extend": "H/E + CHIso + HIso + #sigma_{i#etai#eta}",
        "sieie": "#sigma_{i#etai#eta}",
        "PFECalIso": "EIso",
        "official": "H/E + CHIso + HIso + #sigma_{i#etai#eta} + EIso (Official)",
    }

    # wp 用「黑色線型」對應 tight/medium/loose
    # 原本 medium=7 在某些 ROOT 版本/輸出（特別是 PDF）會看起來接近實線
    # FIX: 改用更常見且明顯的線型：tight=1(實線), medium=2(虛線), loose=3(點線)
    wp_linestyles = {"tight": 1, "medium": 2, "loose": 3}  # solid / dashed / dotted
    # marker 只用來讓點更清楚（同色），不再用來區分 wp
    kind_markers = {"custom": 20, "custom_extend": 21, "sieie": 22, "PFECalIso": 33, "official": 29}
    marker_size = 1.2

    # NEW: legend positions (x1, y1, x2, y2) — tweak only these 4 lines
    LEG_POS_LEG      = (0.16, 0.60, 0.52, 0.90)  # linear: kinds
    LEG_POS_LEG_WP   = (0.57, 0.70, 0.82, 0.85)  # linear: wps
    LEG_POS_LEG_LOG  = (0.38, 0.18, 0.78, 0.48)  # logy:  kinds
    LEG_POS_WP_LOG   = (0.78, 0.28, 0.98, 0.43)  # logy:  wps

    def _kind_wp_from_scenario(s: str) -> Tuple[str, str]:
        # "zgammas_phid_{kind}_{wp}"
        m = re.match(r"^zgammas_phid_(.+)_(tight|medium|loose)$", s)
        if not m:
            return ("unknown", "unknown")
        return (m.group(1), m.group(2))

    def _label_for_scenario(s: str) -> str:
        kind, wp = _kind_wp_from_scenario(s)
        return f"{kind} ({wp})"

    c1 = ROOT.TCanvas(f"c_phid_event_{year}", "", 900, 700)
    # NEW (safety): ensure stat box stays off on this canvas
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    c1.SetMargin(0.13, 0.04, 0.13, 0.08)
    c1.SetTickx()
    c1.SetTicky()

    # NEW: keepalive list for ROOT objects attached to pads (avoid PyROOT GC/dtor order issues)
    _keepalive = []

    # 用第一條 graph 取 axis 物件，但不要用它來「代表第一條線」
    first_key = next(iter(scenario_graphs.keys()))
    g_axis = scenario_graphs[first_key]

    # --- FIX: use TH1 frame (robust) instead of TGraph clone as dummy frame ---
    frame_name = f"frame_phid_event_{year}"
    x_min = 0
    x_max = 31
    h_frame = ROOT.TH1F(frame_name, "", 1, x_min, x_max)
    h_frame.SetDirectory(0)
    _keepalive.append(h_frame)

    h_frame.SetTitle("")
    h_frame.GetXaxis().SetTitle("m_{a} [GeV]")
    h_frame.GetYaxis().SetTitle("Signal Efficiency")
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.15)
    h_frame.GetXaxis().SetTitleSize(0.055)
    h_frame.GetYaxis().SetTitleSize(0.055)
    h_frame.GetXaxis().SetLabelSize(0.05)
    h_frame.GetYaxis().SetLabelSize(0.05)
    h_frame.SetMinimum(0.008)
    h_frame.SetMaximum(0.09)

    h_frame.Draw("AXIS")

    # --- draw graphs one-by-one (no MultiGraph) ---
    for scenario, g in scenario_graphs.items():
        kind, wp = _kind_wp_from_scenario(scenario)
        col = kind_colors.get(kind, 1)

        g.SetFillStyle(0)
        g.SetLineColor(col)
        g.SetMarkerColor(col)
        g.SetLineStyle(int(wp_linestyles.get(wp, 1)))
        g.SetLineWidth(3)
        g.SetMarkerStyle(kind_markers.get(kind, 20))
        g.SetMarkerSize(marker_size)

        g.Draw("LP SAME")

    # --- legends: 5 kinds (colors) + 3 wps (black line styles) ---
    leg = ROOT.TLegend(*LEG_POS_LEG)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.033)

    # 重要：保存 legend 內引用的 dummy 物件，避免被 Python GC 回收造成 ROOT crash
    _legend_keepalive = []

    # ROOT 不同版本對 AddEntry(nullptr, ...) 穩定性不同，改用空字串較保險
    leg.AddEntry("", f"{year}", "")

    # 5 kinds legend (use dummy lines)
    for kind in _PHID_KINDS:
        if kind not in kind_colors:
            continue
        dummy = ROOT.TGraph(1)
        dummy.SetLineColor(kind_colors[kind])
        dummy.SetLineStyle(1)
        dummy.SetLineWidth(3)
        dummy.SetMarkerStyle(kind_markers.get(kind, 20))
        dummy.SetMarkerColor(kind_colors[kind])
        dummy.SetMarkerSize(marker_size)
        _legend_keepalive.append(dummy)
        leg.AddEntry(dummy, legend_map[kind], "lp")

    leg.Draw()

    # NEW: 獨立宣告一個新的 legend 給 3 wps（可自行調整位置）
    leg_wp = ROOT.TLegend(*LEG_POS_LEG_WP)
    leg_wp.SetBorderSize(0)
    leg_wp.SetFillStyle(0)
    leg_wp.SetTextFont(42)
    leg_wp.SetTextSize(0.033)

    wp_label = {"tight": "Tight", "medium": "Medium", "loose": "Loose"}
    for wp in _PHID_WPS:
        dummy_wp = ROOT.TGraph(1)
        dummy_wp.SetLineColor(ROOT.kBlack)
        dummy_wp.SetLineStyle(int(wp_linestyles.get(wp, 1)))
        dummy_wp.SetLineWidth(3)
        dummy_wp.SetMarkerStyle(0)
        _legend_keepalive.append(dummy_wp)
        leg_wp.AddEntry(dummy_wp, wp_label.get(wp, wp), "l")

    leg_wp.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Preliminary}")

    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetTextSize(0.040)
        lat_lumi.DrawLatex(0.96, 0.93, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"cutflowVmA_phid_event_{year}.png"

    # 讓 ROOT 先完成 paint（在某些環境可降低 SaveAs 時的 legend paint 問題）
    c1.Modified()
    c1.Update()

    c1.SaveAs(str(out_png))
    c1.SaveAs(str(out_png.with_suffix(".pdf")))

    # ========== log(y) version ==========
    # 收集所有正的 y 值，避免 log scale 下出現 0/負值導致問題
    _all_pos_ys = []
    for _g in scenario_graphs.values():
        for i in range(int(_g.GetN())):
            x = ctypes.c_double(0.0)
            y = ctypes.c_double(0.0)
            _g.GetPoint(i, x, y)
            if float(y.value) > 0.0:
                _all_pos_ys.append(float(y.value))

    if _all_pos_ys:
        y_min_pos = min(_all_pos_ys)
        y_max_pos = max(_all_pos_ys)

        c2 = ROOT.TCanvas(f"c_phid_event_{year}_logy", "", 900, 700)
        # NEW (safety): ensure stat box stays off on this canvas
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
        c2.SetMargin(0.13, 0.04, 0.13, 0.08)
        c2.SetTickx()
        c2.SetTicky()
        c2.SetLogy()

        # --- FIX: TH1 frame for logy too ---
        h_frame_log = ROOT.TH1F(f"{frame_name}_log", "", 1, x_min, x_max)
        h_frame_log.SetDirectory(0)
        _keepalive.append(h_frame_log)

        h_frame_log.SetTitle("")
        h_frame_log.SetMinimum(0.01)
        h_frame_log.SetMaximum(0.1)
        h_frame_log.GetXaxis().SetTitle("m_{a} [GeV]")
        h_frame_log.GetYaxis().SetTitle("Signal Efficiency")
        h_frame_log.GetXaxis().SetTitleOffset(1.1)
        h_frame_log.GetYaxis().SetTitleOffset(1.15)
        h_frame_log.GetXaxis().SetTitleSize(0.055)
        h_frame_log.GetYaxis().SetTitleSize(0.055)
        h_frame_log.GetXaxis().SetLabelSize(0.05)
        h_frame_log.GetYaxis().SetLabelSize(0.05)

        h_frame_log.Draw("AXIS")

        # --- draw graphs one-by-one (logy) ---
        for scenario, g in scenario_graphs.items():
            g.Draw("LP SAME")

        # --- legends: 5 kinds (colors) + 3 wps (black line styles) ---
        leg_log = ROOT.TLegend(*LEG_POS_LEG_LOG)
        leg_log.SetBorderSize(0)
        leg_log.SetFillStyle(0)
        leg_log.SetTextFont(42)
        leg_log.SetTextSize(0.032)
        leg_log.AddEntry("", f"{year}", "")

        for kind in _PHID_KINDS:
            if kind not in kind_colors:
                continue
            dummy = ROOT.TGraph(1)
            dummy.SetLineColor(kind_colors[kind])
            dummy.SetLineStyle(1)
            dummy.SetLineWidth(3)
            dummy.SetMarkerStyle(kind_markers.get(kind, 20))
            dummy.SetMarkerColor(kind_colors[kind])
            dummy.SetMarkerSize(marker_size)
            _legend_keepalive.append(dummy)
            leg_log.AddEntry(dummy, legend_map[kind], "lp")

        leg_wp_log = ROOT.TLegend(*LEG_POS_WP_LOG)
        leg_wp_log.SetBorderSize(0)
        leg_wp_log.SetFillStyle(0)
        leg_wp_log.SetTextFont(42)
        leg_wp_log.SetTextSize(0.032)

        wp_label = {"tight": "Tight", "medium": "Medium", "loose": "Loose"}
        for wp in _PHID_WPS:
            dummy_wp = ROOT.TGraph(1)
            dummy_wp.SetLineColor(ROOT.kBlack)
            dummy_wp.SetLineStyle(int(wp_linestyles.get(wp, 1)))
            dummy_wp.SetLineWidth(3)
            dummy_wp.SetMarkerStyle(0)
            _legend_keepalive.append(dummy_wp)
            leg_wp_log.AddEntry(dummy_wp, wp_label.get(wp, wp), "l")

        leg_log.Draw()
        leg_wp_log.Draw()

        # FIX: 不用 lat.Clone()；在某些 ROOT 版本 clone 後屬性/NDC 會不見，導致不顯示
        lat2 = ROOT.TLatex()
        lat2.SetNDC()
        lat2.SetTextFont(42)
        lat2.SetTextSize(0.045)
        lat2.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Preliminary}")

        if lumi_fb is not None:
            lat_lumi2 = ROOT.TLatex()
            lat_lumi2.SetNDC()
            lat_lumi2.SetTextFont(42)
            lat_lumi2.SetTextAlign(31)
            lat_lumi2.SetTextSize(0.040)
            lat_lumi2.DrawLatex(0.96, 0.93, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

        out_png_log = out_dir / f"cutflowVmA_phid_event_{year}_logy.png"
        c2.Modified()
        c2.Update()
        c2.SaveAs(str(out_png_log))
        c2.SaveAs(str(out_png_log.with_suffix(".pdf")))

        # --- cleanup: remove only the frames we own ---
        try:
            prim2 = c2.GetListOfPrimitives()
            if prim2:
                prim2.Remove(h_frame_log)
        except Exception:
            pass
        try:
            c2.Close()
        except Exception:
            pass
        del h_frame_log, c2

    # --- cleanup: remove only the frame ---
    try:
        prim1 = c1.GetListOfPrimitives()
        if prim1:
            prim1.Remove(h_frame)
    except Exception:
        pass
    try:
        c1.Close()
    except Exception:
        pass
    del h_frame, c1

    # keep refs alive until after canvases are closed
    _ = (_keepalive, _legend_keepalive)

def main():
    parser = argparse.ArgumentParser(description="Plot PHID event efficiency vs mA per year.")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022preEE). Use 'all' to plot all years.")
    args = parser.parse_args()

    in_dir = Path(baseDir)
    out_dir = Path(outDir)

    # year -> ma -> scenario -> {cut:eff}
    store: Dict[str, Dict[int, Dict[str, Dict[str, float]]]] = {}

    for p in sorted(in_dir.glob("cutflow_*.json")):
        got = _load_cutflow_point(str(p))
        if got is None:
            continue
        year, ma, effs_by_scenario = got
        store.setdefault(year, {})[ma] = effs_by_scenario

    if args.year != "all":
        if args.year in store:
            _plot_year(args.year, store[args.year], out_dir)
        return

    for year in YEAR_ORDER:
        if year in store:
            _plot_year(year, store[year], out_dir)

if __name__ == "__main__":
    main()
