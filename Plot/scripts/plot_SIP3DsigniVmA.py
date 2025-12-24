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
import math   # NEW: significance formula

# NEW: globally disable ROOT stat/fit boxes (Entries/Mean/Std Dev)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

baseDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/cutflow/cutflow_list"

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/SIP3DsigniVmA"

# 你指定要畫的 cuts（分母皆為 all）
CUTS_TO_PLOT = ["event"]

# JSON 裡 cutflow 的 key：改成只比較這兩條
SCENARIOS_TO_PLOT = ["zgammas_ele_w", "zgammas_ele_eleip3d_w"]

legend_map = { "zgammas_ele_w" : "e SIP3D Removal",
               "zgammas_ele_eleip3d_w" : "e SIP3D"}

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, '2024':108.95,
            'combined_run3':61.89 }
# 移除：YEAR 未定義時會直接噴錯；改為在 _plot_year() 依 year 動態取用
# LUMI_FB = float(lumiMap[YEAR])

# Name of Bkg Sample
name_DYG_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
name_DYG_2023 = ["DYGto2LG_10to100"]
name_DYG_2024 = ["DYGto2LG_10to100"]
name_DYJet_2022 = ["DYJetsToLL"]
name_DYJet_2023 = ["DYJetsToLL"]
name_DYJet_2024 = ["DYJetsTo2E","DYJetsTo2Mu","DYJetsTo2Tau"]
# Year of Bkg
year_DYG_2022 = ["2022preEE", "2022postEE"]
year_DYG_2023 = ["2023preBPix", "2023postBPix"]
year_DYG_2024 = ["2024"]
years_DYJet_2022 = ["2022preEE","2022postEE"]
year_DYJet_2023  = ["2023preBPix", "2023postBPix"]
years_DYJet_2024 = ["2024"]

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

    # scenario -> {cut: count}
    effs_by_scenario: Dict[str, Dict[str, float]] = {}
    for scenario in SCENARIOS_TO_PLOT:
        cf = cutflows.get(scenario) or {}
        allv = cf.get("all")
        if allv is None or float(allv) <= 0:
            continue

        effs: Dict[str, float] = {}

        # NEW: cache "all" and "all cuts" (or "event" fallback) for efficiency
        effs["all"] = float(allv)
        allcuts = cf.get("all cuts")
        if allcuts is None:
            allcuts = cf.get("event")  # fallback (your old usage)
        if allcuts is not None:
            effs["all cuts"] = float(allcuts)

        for c in CUTS_TO_PLOT:
            v = cf.get(c)
            if v is None:
                continue
            effs[c] = float(v)  # store count

        if effs:
            effs_by_scenario[scenario] = effs

    if not effs_by_scenario:
        return None
    return year, ma, effs_by_scenario

def _load_bkg_event_count_by_year(in_dir: Path) -> Dict[str, float]:
    """
    讀入每年 bkg 的 'all cuts 之後' 數目（這裡以 cutflow JSON 的 event 欄位為準）。
    回傳: year -> bkg_event_count
    假設 bkg cutflow JSON 檔名含有 year 字串（同 signal）。
    """
    out: Dict[str, float] = {}
    for p in sorted(in_dir.glob("cutflow_*.json")):
        fn = p.name
        year = _parse_year_from_filename(fn)
        if year is None:
            continue
        # 只挑 bkg 檔（檔名規則不確定，先用常見關鍵字；不符合就略過）
        low = fn.lower()
        if not any(k in low for k in ["bkg", "background", "dyg", "dyjet", "dyjets"]):
            continue
        try:
            with open(str(p), "r") as f:
                data = json.load(f)
        except Exception:
            continue
        cutflows = data.get("cutflows") or {}

        # 對齊 signal JSON 的格式：優先用 "zgammas" 這條當作 bkg 代表
        b_evt = None
        cf0 = cutflows.get("zgammas")
        if isinstance(cf0, dict):
            if cf0.get("event") is not None:
                b_evt = float(cf0["event"])
            elif cf0.get("all cuts") is not None:
                b_evt = float(cf0["all cuts"])
        # fallback：真的沒有 zgammas 才掃描
        if b_evt is None:
            for _, cf in cutflows.items():
                if not isinstance(cf, dict):
                    continue
                evt = cf.get("event") if cf.get("event") is None else cf.get("all cuts")
                if evt is None:
                    continue
                try:
                    b_evt = float(evt)
                    break
                except Exception:
                    continue

        if b_evt is None:
            continue
        out[year] = out.get(year, 0.0) + b_evt
    return out

def _load_bkg_eff_by_year(in_dir: Path) -> Dict[str, float]:
    """
    回傳 year -> eff_bkg，其中 eff_bkg = (all cuts)/(all)。
    bkg 會跨檔累加 all / all cuts，再做比值。
    """
    sum_all: Dict[str, float] = {}
    sum_allcuts: Dict[str, float] = {}

    for p in sorted(in_dir.glob("cutflow_*.json")):
        fn = p.name
        year = _parse_year_from_filename(fn)
        if year is None:
            continue
        low = fn.lower()
        if not any(k in low for k in ["bkg", "background", "dyg", "dyjet", "dyjets"]):
            continue
        try:
            with open(str(p), "r") as f:
                data = json.load(f)
        except Exception:
            continue

        cutflows = data.get("cutflows") or {}

        # 以 zgammas 為優先代表；沒有就掃第一個 dict
        cf0 = cutflows.get("zgammas") if isinstance(cutflows.get("zgammas"), dict) else None
        if cf0 is None:
            for _, cf in cutflows.items():
                if isinstance(cf, dict):
                    cf0 = cf
                    break
        if cf0 is None:
            continue

        a = cf0.get("all")
        ac = cf0.get("all cuts")
        if ac is None:
            ac = cf0.get("event")  # fallback
        try:
            if a is None or float(a) <= 0:
                continue
            if ac is None:
                continue
            sum_all[year] = sum_all.get(year, 0.0) + float(a)
            sum_allcuts[year] = sum_allcuts.get(year, 0.0) + float(ac)
        except Exception:
            continue

    out: Dict[str, float] = {}
    for y, a in sum_all.items():
        ac = sum_allcuts.get(y, 0.0)
        out[y] = (ac / a) if a > 0 else 0.0
    return out

def _significance_asimov_like(s: float, b: float) -> float:
    """開根號版本：sqrt( (2*(s+b)*log(1+s/b)) - 2*s )"""
    if s <= 0.0:
        return 0.0
    if b <= 0.0:
        return 0.0
    val = (2.0 * (s + b) * math.log(1.0 + (s / b))) - (2.0 * s)
    # 數值誤差保護：避免極小負值導致 sqrt 出錯
    if val <= 0.0:
        return 0.0
    return math.sqrt(val)

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
    # NEW: bkg event count for this year (precomputed in main, injected via global)
    b_evt = float(_BKG_EVENT_BY_YEAR.get(year, 0.0))

    # scenario -> graph(significance vs mA)
    scenario_graphs: Dict[str, ROOT.TGraph] = {}
    for scenario in SCENARIOS_TO_PLOT:
        xs, ys = [], []
        for ma in mas:
            per_ma = points.get(ma) or {}
            effs = per_ma.get(scenario) or {}
            if "event" in effs:
                xs.append(ma)
                # 現在 effs["event"] 已是「count」
                s_evt = float(effs["event"])
                ys.append(_significance_asimov_like(s_evt, b_evt))
        if len(xs) >= 1:
            scenario_graphs[scenario] = _g_from_xy(xs, ys)

    if not scenario_graphs:
        return

    # 兩條線的外觀（固定顏色）
    scenario_style = {
        "zgammas_w":            {"color": _root_color("#1f77b4"), "ls": 1, "ms": 20, "label": "e SIP3D < 4 Removal"},
        "zgammas_eleip3d_w":    {"color": _root_color("#d62728"), "ls": 2, "ms": 21, "label": "e SIP3D < 4"},
    }

    c1 = ROOT.TCanvas(f"c_w_vs_eleip3d_{year}", "", 900, 700)
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
    h_frame.GetYaxis().SetTitle("Significance (AMS)")
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.15)
    h_frame.GetXaxis().SetTitleSize(0.055)
    h_frame.GetYaxis().SetTitleSize(0.055)
    h_frame.GetXaxis().SetLabelSize(0.05)
    h_frame.GetYaxis().SetLabelSize(0.05)
    h_frame.SetMinimum(10.0)
    h_frame.SetMaximum(200.0)  # 先給個安全範圍；可再依資料自動調

    h_frame.Draw("AXIS")

    # --- draw graphs one-by-one ---
    for scenario in SCENARIOS_TO_PLOT:
        g = scenario_graphs.get(scenario)
        if not g:
            continue
        st = scenario_style.get(scenario, {"color": 1, "ls": 1, "ms": 20, "label": scenario})

        g.SetFillStyle(0)
        g.SetLineColor(int(st["color"]))
        g.SetMarkerColor(int(st["color"]))
        g.SetLineStyle(int(st["ls"]))
        g.SetLineWidth(3)
        g.SetMarkerStyle(int(st["ms"]))
        g.SetMarkerSize(1.2)
        g.Draw("LP SAME")

    # --- legend: just two scenarios ---
    leg = ROOT.TLegend(0.16, 0.72, 0.52, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)

    leg.AddEntry("", f"{year}", "")
    for scenario in SCENARIOS_TO_PLOT:
        g = scenario_graphs.get(scenario)
        if not g:
            continue
        st = scenario_style.get(scenario, {"label": scenario})
        leg.AddEntry(g, st["label"], "lp")
    leg.Draw()

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
    out_png = out_dir / f"SIP3D_significanceVmA_{year}.pdf"

    # 讓 ROOT 先完成 paint（在某些環境可降低 SaveAs 時的 legend paint 問題）
    c1.Modified()
    c1.Update()

    c1.SaveAs(str(out_png))

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
    _ = _keepalive

    # =========================
    # NEW: plot eff_sig / eff_bkg vs mA
    # =========================
    eff_bkg = float(_BKG_EFF_BY_YEAR.get(year, 0.0))

    scenario_ratio_graphs: Dict[str, ROOT.TGraph] = {}
    for scenario in SCENARIOS_TO_PLOT:
        xs, ys = [], []
        for ma in mas:
            per_ma = points.get(ma) or {}
            effs = per_ma.get(scenario) or {}
            a = float(effs.get("all", 0.0))
            ac = float(effs.get("all cuts", 0.0))
            if a <= 0.0 or ac <= 0.0:
                continue
            eff_sig = ac / a
            if eff_bkg <= 0.0:
                continue
            xs.append(ma)
            ys.append(eff_sig / eff_bkg)
        if xs:
            scenario_ratio_graphs[scenario] = _g_from_xy(xs, ys)

    if scenario_ratio_graphs:
        c2 = ROOT.TCanvas(f"c_effSigOverEffBkg_{year}", "", 900, 700)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
        c2.SetMargin(0.13, 0.04, 0.13, 0.08)
        c2.SetTickx()
        c2.SetTicky()

        _keepalive2 = []
        frame2 = ROOT.TH1F(f"frame_effSigOverEffBkg_{year}", "", 1, 0, 31)
        frame2.SetDirectory(0)
        _keepalive2.append(frame2)

        frame2.SetTitle("")
        frame2.GetXaxis().SetTitle("m_{a} [GeV]")
        frame2.GetYaxis().SetTitle("#varepsilon_{sig} / #varepsilon_{bkg}   (#varepsilon = all cuts / all)")
        frame2.GetXaxis().SetTitleOffset(1.1)
        frame2.GetYaxis().SetTitleOffset(1.25)
        frame2.GetXaxis().SetTitleSize(0.055)
        frame2.GetYaxis().SetTitleSize(0.045)
        frame2.GetXaxis().SetLabelSize(0.05)
        frame2.GetYaxis().SetLabelSize(0.05)
        frame2.SetMinimum(80.0)
        frame2.SetMaximum(220.0)  # 視你的結果可再調

        frame2.Draw("AXIS")

        for scenario in SCENARIOS_TO_PLOT:
            g = scenario_ratio_graphs.get(scenario)
            if not g:
                continue
            st = scenario_style.get(scenario, {"color": 1, "ls": 1, "ms": 20, "label": scenario})
            g.SetFillStyle(0)
            g.SetLineColor(int(st["color"]))
            g.SetMarkerColor(int(st["color"]))
            g.SetLineStyle(int(st["ls"]))
            g.SetLineWidth(3)
            g.SetMarkerStyle(int(st["ms"]))
            g.SetMarkerSize(1.2)
            g.Draw("LP SAME")

        leg2 = ROOT.TLegend(0.16, 0.72, 0.62, 0.88)
        leg2.SetBorderSize(0)
        leg2.SetFillStyle(0)
        leg2.SetTextFont(42)
        leg2.SetTextSize(0.045)
        leg2.AddEntry("", f"{year}   (eff_bkg={eff_bkg:.3g})", "")
        for scenario in SCENARIOS_TO_PLOT:
            g = scenario_ratio_graphs.get(scenario)
            if not g:
                continue
            st = scenario_style.get(scenario, {"label": scenario})
            leg2.AddEntry(g, st["label"], "lp")
        leg2.Draw()

        lat2 = ROOT.TLatex()
        lat2.SetNDC()
        lat2.SetTextFont(42)
        lat2.SetTextSize(0.045)
        lat2.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Preliminary}")

        lumi_fb = _lumi_fb_for_year(year)
        if lumi_fb is not None:
            lat_lumi2 = ROOT.TLatex()
            lat_lumi2.SetNDC()
            lat_lumi2.SetTextFont(42)
            lat_lumi2.SetTextAlign(31)
            lat_lumi2.SetTextSize(0.040)
            lat_lumi2.DrawLatex(0.96, 0.93, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

        out_dir.mkdir(parents=True, exist_ok=True)
        out2 = out_dir / f"SIP3D_effSigOverEffBkgVmA_{year}.pdf"
        c2.Modified()
        c2.Update()
        c2.SaveAs(str(out2))

        try:
            prim2 = c2.GetListOfPrimitives()
            if prim2:
                prim2.Remove(frame2)
        except Exception:
            pass
        try:
            c2.Close()
        except Exception:
            pass
        del frame2, c2
        _ = _keepalive2

def main():
    parser = argparse.ArgumentParser(description="Plot PHID event efficiency vs mA per year.")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022preEE). Use 'all' to plot all years.")
    args = parser.parse_args()

    in_dir = Path(baseDir)
    out_dir = Path(outDir)

    # NEW: load bkg counts once (year -> bkg event count)
    global _BKG_EVENT_BY_YEAR
    _BKG_EVENT_BY_YEAR = _load_bkg_event_count_by_year(in_dir)

    # NEW: bkg efficiency per year (all cuts / all)
    global _BKG_EFF_BY_YEAR
    _BKG_EFF_BY_YEAR = _load_bkg_eff_by_year(in_dir)

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
