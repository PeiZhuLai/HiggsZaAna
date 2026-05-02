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

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/eachphidVmA"

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

    # tight/medium/loose 用黑色線型
    wp_linestyles = {"tight": 1, "medium": 2, "loose": 3}  # solid / dashed / dotted
    wp_label = {"tight": "Tight", "medium": "Medium", "loose": "Loose"}

    # 同一 kind 內三條線：同 marker，靠線型區分
    marker_style = 20
    marker_size = 1.2

    # NEW: legend positions（每張圖只有 WP legend）
    LEG_POS         = (0.12, 0.80, 0.55, 0.90)
    LEG_POS_WP      = (0.12, 0.65, 0.38, 0.80)

    # NEW: per-photonID y-range (ymin, ymax)
    y_range = {"custom": (0.04, 0.09), "custom_extend": (0.02, 0.09), "sieie": (0.02, 0.08),
               "PFECalIso": (0.006, 0.082), "official": (0.008, 0.08)}
    _ymin_default, _ymax_default = (0.008, 0.09)

    def _kind_wp_from_scenario(s: str) -> Tuple[str, str]:
        m = re.match(r"^zgammas_phid_(.+)_(tight|medium|loose)$", s)
        if not m:
            return ("unknown", "unknown")
        return (m.group(1), m.group(2))

    def _scenario(kind: str, wp: str) -> str:
        return f"zgammas_phid_{kind}_{wp}"

    # NEW: allow a few common key variants (to survive JSON naming inconsistencies)
    def _scenario_candidates(kind: str, wp: str) -> List[str]:
        wp_l = wp.lower()
        wp_u = wp.upper()
        return [
            f"zgammas_phid_{kind}_{wp_l}",      # expected
            f"zgammas_phid_{kind}_{wp_u}",      # wp upper
            f"zgammas_phid_{wp_l}_{kind}",      # swapped order (seen in some legacy outputs)
            f"zgammas_phid_{wp_u}_{kind}",      # swapped + upper
        ]

    def _make_graph(kind: str, wp: str) -> Optional[ROOT.TGraph]:
        xs, ys = [], []
        missing_mas = []

        for ma in mas:
            per_ma = points.get(ma) or {}

            effs = None
            picked_key = None
            for k in _scenario_candidates(kind, wp):
                cand = per_ma.get(k)
                if cand and ("event" in cand):
                    effs = cand
                    picked_key = k
                    break

            if effs is None:
                missing_mas.append(ma)
                continue

            xs.append(ma)
            ys.append(float(effs["event"]))

        # NEW: debug summary for custom (so you can see why only tight appears)
        if kind == "custom":
            scen_keys = set()
            # show a small sample of keys to avoid huge spam
            for ma in mas[:3]:
                scen_keys |= set((points.get(ma) or {}).keys())
            print(
                f"[DEBUG] year={year} kind={kind} wp={wp} "
                f"npoints={len(xs)} missing={len(missing_mas)} "
                f"example_keys={sorted([k for k in scen_keys if 'zgammas_phid' in k])[:20]}"
            )
            if missing_mas:
                print(f"[DEBUG] year={year} kind={kind} wp={wp} missing_mA={missing_mas}")

        if not xs:
            return None
        return _g_from_xy(xs, ys)

    def _draw_one(kind: str) -> None:
        col = kind_colors.get(kind, 1)

        c = ROOT.TCanvas(f"c_phid_{kind}_event_{year}", "", 900, 700)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
        c.SetMargin(0.13, 0.04, 0.13, 0.08)
        c.SetTickx()
        c.SetTicky()

        _keepalive = []
        _legend_keepalive = []

        # frame
        x_min, x_max = 0, 31
        frame_name = f"frame_phid_{kind}_event_{year}"
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

        ymin, ymax = y_range.get(kind, (_ymin_default, _ymax_default))
        h_frame.SetMinimum(float(ymin))
        h_frame.SetMaximum(float(ymax))

        h_frame.Draw("AXIS")

        # graphs: exactly 3 WPs
        graphs = {}
        for wp in _PHID_WPS:  # ["loose","medium","tight"]
            g = _make_graph(kind, wp)
            if g is None:
                continue
            g.SetFillStyle(0)
            g.SetLineColor(col)
            g.SetMarkerColor(col)
            g.SetLineStyle(int(wp_linestyles.get(wp, 1)))
            g.SetLineWidth(3)
            g.SetMarkerStyle(marker_style)
            g.SetMarkerSize(marker_size)
            g.Draw("LP SAME")
            graphs[wp] = g
            _keepalive.append(g)

        if not graphs:
            try:
                c.Close()
            except Exception:
                pass
            return

        # labels
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

        # title-ish (year + kind description) 用 legend 無框文字
        leg_title = ROOT.TLegend(*LEG_POS)
        leg_title.SetBorderSize(0)
        leg_title.SetFillStyle(0)
        leg_title.SetTextFont(42)
        leg_title.SetTextSize(0.045)
        leg_title.AddEntry("", f"{year}", "")
        leg_title.AddEntry("", legend_map.get(kind, kind), "")
        leg_title.Draw()

        # WP legend (改成使用該 kind 的顏色；用線型區分 WP)
        leg_wp = ROOT.TLegend(*LEG_POS_WP)
        leg_wp.SetBorderSize(0)
        leg_wp.SetFillStyle(0)
        leg_wp.SetTextFont(42)
        leg_wp.SetTextSize(0.033)

        for wp in _PHID_WPS:
            dummy_wp = ROOT.TGraph(1)
            dummy_wp.SetLineColor(col)  # was ROOT.kBlack
            dummy_wp.SetLineStyle(int(wp_linestyles.get(wp, 1)))
            dummy_wp.SetLineWidth(3)
            dummy_wp.SetMarkerStyle(0)
            _legend_keepalive.append(dummy_wp)
            leg_wp.AddEntry(dummy_wp, wp_label.get(wp, wp), "l")
        leg_wp.Draw()

        out_dir.mkdir(parents=True, exist_ok=True)
        fname = f"cutflowVmA_phid_{kind}_event_{year}"
        base = out_dir / fname
        c.Modified()
        c.Update()
        c.SaveAs(str(base.with_suffix(".png")))
        c.SaveAs(str(base.with_suffix(".pdf")))

        # cleanup
        try:
            prim = c.GetListOfPrimitives()
            if prim:
                prim.Remove(h_frame)
        except Exception:
            pass
        try:
            c.Close()
        except Exception:
            pass
        _ = (_keepalive, _legend_keepalive)

    # per-kind plots: each kind -> 1 plot (3 WPs)
    for kind in _PHID_KINDS:
        _draw_one(kind)

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
