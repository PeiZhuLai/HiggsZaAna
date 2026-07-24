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

# ---------------------------------------------------------------------------
# 主 cutflowVmA（inclusive，cut_type = zgammas_w）
# 順序 = e selection / e+mu selection / e+mu+trigger / Z candidate /
#        a-candidate(photon) 逐項分解 / a candidate / Za candidate。
# 說明：內部 key `trig_cut` 是累積(含 e/mu 選取)，故 legend 標為 "e+#mu+trigger"；
#      末點用 "all cuts"（tagger 不再輸出舊的 "event" key）。
#      ph_lepovlp == has_2g_cand（a candidate），故主圖不另畫 has_2g_cand 避免重複點。
# ---------------------------------------------------------------------------
CUTS_TO_PLOT = [
    "e_sel", "emu_sel", "trig_cut", "has_z_cand",
    "ph_ge2", "ph_pt", "ph_eta", "ph_id_hoe", "ph_id_chiso", "ph_id_hcaliso",
    "ph_eveto", "ph_lepovlp", "all cuts",
]
CUTFLOW_KEY = "zgammas_w"

MAIN_NOTE = "Z #rightarrow (ee or #mu#mu) 67% #times electron acceptance (70%) = 47%"

# ---------------------------------------------------------------------------
# dR(gamma,gamma)-binned cutflowVmA（同一張 cutflowVmA，但依 a-candidate 兩光子
# dR 分三互斥區間；cut_type = zgammas_<dr_key>_w）。
# dR 只在 a-candidate（pair 建立）後才有定義，故 pair 之前的步驟三個 bin 相同；
# has_2g_cand / sel_h 為該 dR bin 的 a-candidate / Za-candidate。
# ---------------------------------------------------------------------------
DR_BINS = [
    ("dr_lt_0p1",  "dR(#gamma,#gamma) < 0.1"),
    ("dr_0p1_0p3", "0.1 #leq dR(#gamma,#gamma) < 0.3"),
    ("dr_gt_0p3",  "dR(#gamma,#gamma) #geq 0.3"),
]
DR_CUTS_TO_PLOT = [
    "e_sel", "emu_sel", "trig_cut", "has_z_cand",
    "ph_ge2", "ph_pt", "ph_eta", "ph_id_hoe", "ph_id_chiso", "ph_id_hcaliso",
    "ph_eveto", "ph_lepovlp", "has_2g_cand", "sel_h",
]

# 共用配色/marker（涵蓋主圖與 dR 圖的所有 key）
COLOR_MAP = {
    "e_sel": "#8338EC", "emu_sel": "#3A86FF", "trig_cut": "#3185FC", "has_z_cand": "#FF8000",
    "ph_ge2": "#FFBE0B", "ph_pt": "#FB5607", "ph_eta": "#8AC926", "ph_id_hoe": "#1982C4",
    "ph_id_chiso": "#6A4C93", "ph_id_hcaliso": "#33BA91", "ph_eveto": "#C1121F",
    "ph_lepovlp": "#00A896", "has_2g_cand": "#F72585", "all cuts": "#E3170A", "sel_h": "#E3170A",
}
MARKER_MAP = {
    "e_sel": 24, "emu_sel": 25, "trig_cut": 20, "has_z_cand": 22,
    "ph_ge2": 26, "ph_pt": 32, "ph_eta": 27, "ph_id_hoe": 28,
    "ph_id_chiso": 30, "ph_id_hcaliso": 21, "ph_eveto": 34, "ph_lepovlp": 23,
    "has_2g_cand": 29, "all cuts": 33, "sel_h": 33,
}
_BIG_MARKER_KEYS = {"all cuts", "sel_h"}  # 最終 Za candidate 用大 marker

# legend labels（主圖）
LEGEND_MAIN = {
    "e_sel": "e selection", "emu_sel": "e + #mu selection", "trig_cut": "e + #mu + trigger",
    "has_z_cand": "+ Z candidate (m_{ll} > 50)",
    "ph_ge2": "+ N_{#gamma} #geq 2", "ph_pt": "+ 2#gamma p_{T} > 10 GeV",
    "ph_eta": "+ 2#gamma #eta acceptance", "ph_id_hoe": "+ 2#gamma H/E",
    "ph_id_chiso": "+ 2#gamma PF ch. iso", "ph_id_hcaliso": "+ 2#gamma PF HCal iso",
    "ph_eveto": "+ 2#gamma e-veto", "ph_lepovlp": "+ 2#gamma lep/#gamma overlap (a cand.)",
    "all cuts": "+ Za candidate",
}
# legend labels（dR 圖：多了 has_2g_cand / sel_h 為該 dR bin 的 a/Za candidate）
LEGEND_DR = dict(LEGEND_MAIN)
LEGEND_DR.update({
    "ph_lepovlp": "+ 2#gamma lep/#gamma overlap (a cand.)",
    "has_2g_cand": "+ a cand. in this #DeltaR bin",
    "sel_h": "+ Za cand. in this #DeltaR bin",
})

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = {
    "2016" : 36.31,
    "2016preVFP" : 19.51,
    "2016postVFP" : 16.80,
    "2017" : 41.48,
    "2018" : 59.83,
    "2022" : 34.6521,
    "2022preEE": 7.99,
    "2022postEE": 26.68,
    "2022": 34.67,
    "2023preBPix": 17.96,
    "2023postBPix": 9.68,
    "2023" : 27.64,
    "2024": 109.82,
    "2025": 110.67,
    "combined_run3" : 172.13
}

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

def _load_cutflow_point(path: str, cutflow_key: str, cuts: List[str]) -> Optional[Tuple[str, int, Dict[str, float]]]:
    """讀單一 signal JSON，回傳 (year, ma, {cut: eff%})。分母為該 cut_type 的 'all'。"""
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
    cf = (data.get("cutflows") or {}).get(cutflow_key) or {}
    allv = cf.get("all")
    if allv is None or float(allv) <= 0:
        return None
    effs: Dict[str, float] = {}
    for c in cuts:
        v = cf.get(c)
        if v is None:
            continue
        effs[c] = (float(v) / float(allv)) * 100.0
    if not effs:
        return None
    return year, ma, effs

def _build_store(in_dir: Path, cutflow_key: str, cuts: List[str]) -> Dict[str, Dict[int, Dict[str, float]]]:
    """掃 baseDir 的 signal JSON → year -> ma -> {cut: eff%}。"""
    store: Dict[str, Dict[int, Dict[str, float]]] = {}
    for p in sorted(in_dir.glob("cutflow_*.json")):
        got = _load_cutflow_point(str(p), cutflow_key, cuts)
        if got is None:
            continue
        year, ma, effs = got
        store.setdefault(year, {})[ma] = effs
    return store

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

def _draw_cutflow(
    year: str,
    points: Dict[int, Dict[str, float]],
    out_dir: Path,
    *,
    cuts: List[str],
    legend_map: Dict[str, str],
    out_name: str,
    top_lines: List[Tuple[float, float, float, str]],
) -> None:
    """通用 cutflowVmA 繪圖（效率 vs m_a，每個 cut 一條線）。主圖與 dR 圖共用。

    top_lines: [(x_ndc, y_ndc, textsize, text), ...] 放在 legend 上方（era + 說明/dR bin）。
    """
    mas = sorted(points.keys())
    if not mas:
        return

    graphs: Dict[str, ROOT.TGraph] = {}
    for c in cuts:
        xs, ys = [], []
        for ma in mas:
            if c in points[ma]:
                xs.append(ma)
                ys.append(points[ma][c])
        if len(xs) >= 1:
            graphs[c] = _g_from_xy(xs, ys)
    if not graphs:
        return

    color_map = {k: _root_color(COLOR_MAP.get(k, "#000000")) for k in cuts}

    c1 = ROOT.TCanvas(f"c_cutflow_{out_name}", "", 800, 600)
    c1.SetMargin(0.13, 0.04, 0.13, 0.08)
    c1.SetTickx()
    c1.SetTicky()

    def _msize(k):
        return 2.3 if k in _BIG_MARKER_KEYS else 1.7

    # 依 cuts 順序畫（第一個建軸）
    first_key = next((k for k in cuts if k in graphs), None)
    if first_key is None:
        return
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
    g0.SetMaximum(100)  # headroom for the 2-column legend at the top
    g0.SetLineColor(color_map.get(first_key, 1))
    g0.SetMarkerColor(color_map.get(first_key, 1))
    g0.SetMarkerStyle(MARKER_MAP.get(first_key, 20))
    g0.SetMarkerSize(_msize(first_key))
    g0.SetLineWidth(3)
    g0.Draw("ALPC")

    for k in cuts:
        if k == first_key or k not in graphs:
            continue
        g = graphs[k]
        col = color_map.get(k, 1)
        g.SetLineColor(col)
        g.SetMarkerColor(col)
        g.SetMarkerStyle(MARKER_MAP.get(k, 20))
        g.SetMarkerSize(_msize(k))
        g.SetLineWidth(3)
        g.Draw("LPC SAME")

    # 拆細後線變多：用 2 欄 legend 放在圖頂（y 上限拉到 100 留白），避免重疊。
    leg = ROOT.TLegend(0.15, 0.55, 0.83, 0.835)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.033)
    for k in cuts:
        if k in graphs:
            leg.AddEntry(graphs[k], legend_map.get(k, k), "lp")
    leg.Draw()

    # 頂端註記：每筆 (x, y, size, text, align)，align 為 ROOT TLatex 對齊碼
    # (11=左下, 31=右下)。era 用左對齊放左上；說明/dR bin 用右對齊放右上角。
    _keep_lat = []
    for (x, y, size, text, align) in top_lines:
        lt = ROOT.TLatex()
        lt.SetNDC(); lt.SetTextFont(42); lt.SetTextSize(size)
        lt.SetTextAlign(align)
        lt.DrawLatex(x, y, text)
        _keep_lat.append(lt)

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Simulation}")

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
    c1.SaveAs(str(out_dir / f"{out_name}.pdf"))

def _plot_main_year(year: str, points: Dict[int, Dict[str, float]], out_dir: Path) -> None:
    _draw_cutflow(
        year, points, out_dir,
        cuts=CUTS_TO_PLOT,
        legend_map=LEGEND_MAIN,
        out_name=f"cutflowVmA_{CUTFLOW_KEY}_{year}",
        top_lines=[
            (0.15, 0.855, 0.042, year, 11),          # era：左上
            (0.955, 0.858, 0.024, MAIN_NOTE, 31),    # 說明：右上角、靠右對齊
        ],
    )

def _plot_dr_year(year: str, dr_key: str, dr_label: str,
                  points: Dict[int, Dict[str, float]], out_dir: Path) -> None:
    _draw_cutflow(
        year, points, out_dir,
        cuts=DR_CUTS_TO_PLOT,
        legend_map=LEGEND_DR,
        out_name=f"cutflow_by_dr_{dr_key}_{year}",
        top_lines=[
            (0.15, 0.855, 0.042, year, 11),          # era：左上
            (0.955, 0.858, 0.032, dr_label, 31),     # dR bin：右上角、靠右對齊
        ],
    )

def main():
    parser = argparse.ArgumentParser(description="Plot cutflow efficiency vs mA per year (main + dR-binned).")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022preEE). Use 'all' to plot all years.")
    parser.add_argument(
        "--mode", default="all", choices=["main", "by-dr", "all"],
        help="main: inclusive cutflowVmA; by-dr: 依 dR(gamma,gamma) 三區間; all: 兩者都畫(預設)。",
    )
    parser.add_argument(
        "--unweighted", action="store_true",
        help="用未加權 cut_type（zgammas / zgammas_dr_*，預設用 *_w 加權）。",
    )
    args = parser.parse_args()

    in_dir = Path(baseDir)
    out_dir = Path(outDir)
    years = YEAR_ORDER if args.year == "all" else [args.year]
    wsfx = "" if args.unweighted else "_w"

    # --- 主 cutflowVmA ---
    if args.mode in ("main", "all"):
        key = "zgammas" if args.unweighted else "zgammas_w"
        global CUTFLOW_KEY
        CUTFLOW_KEY = key  # 讓輸出檔名沿用 cutflowVmA_<key>_<year>
        store = _build_store(in_dir, key, CUTS_TO_PLOT)
        for year in years:
            if year in store:
                _plot_main_year(year, store[year], out_dir)

    # --- dR-binned cutflowVmA ---
    if args.mode in ("by-dr", "all"):
        for dr_key, dr_label in DR_BINS:
            key = f"zgammas_{dr_key}{wsfx}"
            store = _build_store(in_dir, key, DR_CUTS_TO_PLOT)
            for year in years:
                if year in store:
                    _plot_dr_year(year, dr_key, dr_label, store[year], out_dir)

if __name__ == "__main__":
    main()
