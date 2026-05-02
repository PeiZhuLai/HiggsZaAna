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

baseDir = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal"

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/dREff"

# Define the list of mA values
mAs = ["mA_M1", "mA_M2", "mA_M3", "mA_M4", "mA_M5", "mA_M6", "mA_M7", "mA_M8", "mA_M9", "mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]


YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, '2024':108.95,
            'combined_run3':61.89 }


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
                # FIX: 原本條件寫反，會導致永遠選錯欄位
                evt = cf.get("event") if cf.get("event") is not None else cf.get("all cuts")
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

# NEW: var_dR_g1g2 cut parsing utilities
_VAR_DR_RE = re.compile(r"var[_ ]?dR[_ ]?g1g2[^0-9]*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)

def _is_bkg_filename(fn: str) -> bool:
    low = fn.lower()
    return any(k in low for k in ["bkg", "background", "dyg", "dyjet", "dyjets"])

def _parse_var_dR_cut_from_cutname(cutname: str) -> Optional[float]:
    """
    從 cut 名稱中抓 var_dR_g1g2 的 cut 值；支援像:
      "var_dR_g1g2 < 0.6", "var_dR_g1g2_0p6", "cut_var_dR_g1g2_1.2" ...
    回傳 float cut；抓不到回 None
    """
    if not cutname:
        return None
    m = _VAR_DR_RE.search(str(cutname))
    if not m:
        return None
    try:
        v = float(m.group(1))
    except Exception:
        return None
    # restrict to requested range
    if v < 0.0 or v > 3.0:
        return None
    return v

def _load_eff_scan_from_root_by_year_ma(
    base_dir: Path,
    *,
    years: List[str],
    ma_dirs: List[str],
    tree_name: str = "inclusive",
    dr_branch: str = "var_dR_g1g2",
    dr_min: float = 0.0,
    dr_max: float = 4.0,
    dr_step: float = 0.05,
) -> Dict[str, Dict[int, Dict[float, float]]]:
    """
    直接從 ROOT 檔讀取 inclusive tree，掃描 cut: var_dR_g1g2 > x 的效率：
      eff(x) = N_pass(x) / N_all
    回傳: year -> ma(int) -> {cut_value: eff_pass_over_all}

    ROOT 路徑假設：
      base_dir / mA_M{ma} / {year}.root
    """
    n_steps = int(round((dr_max - dr_min) / dr_step))
    cuts = [round(dr_min + i * dr_step, 10) for i in range(n_steps + 1)]

    out: Dict[str, Dict[int, Dict[float, float]]] = {y: {} for y in years}

    for ma_tag in ma_dirs:
        m = re.match(r"mA_M(\d+)$", ma_tag)
        if not m:
            continue
        ma = int(m.group(1))
        ma_path = base_dir / ma_tag
        if not ma_path.is_dir():
            continue

        for year in years:
            fpath = ma_path / f"{year}.root"
            if not fpath.exists():
                continue

            try:
                with uproot.open(str(fpath)) as f:
                    t = f[tree_name]
                    keys = set(t.keys())
                    if dr_branch not in keys:
                        continue
                    arrs = t.arrays([dr_branch], library="ak")
            except Exception:
                continue

            try:
                dr = arrs[dr_branch]
            except Exception:
                continue
            if dr is None:
                continue

            # flatten + numeric mask
            dr = ak.to_numpy(ak.flatten(dr, axis=None))
            if dr.size == 0:
                continue
            msk = np.isfinite(dr) & (dr >= 0.0)
            dr = dr[msk]
            if dr.size == 0:
                continue

            # N_all after cleaning
            n_all = float(dr.size)

            # sort once to compute N_pass efficiently
            order = np.argsort(dr, kind="mergesort")
            drs = dr[order]

            ym = out.setdefault(year, {}).setdefault(ma, {})
            for x in cuts:
                # pass: dr > x
                # N_le = count(dr <= x) = searchsorted(right)
                n_le = int(np.searchsorted(drs, x, side="right"))
                n_pass = int(drs.size - n_le)
                eff = (float(n_pass) / n_all) if n_all > 0 else 0.0
                ym[x] = eff  # overwrite is fine (one file per year/mA)

    out = {y: m for y, m in out.items() if m}
    return out

def _plot_sumw_vs_var_dR_year(
    year: str,
    by_ma: Dict[int, Dict[float, float]],
    out_dir: Path,
    *,
    scenario_label: str = "efficiency",
    out_name: Optional[str] = None,  # NEW
) -> None:
    mas = sorted(by_ma.keys())
    if not mas:
        return

    c = ROOT.TCanvas(f"c_sumw_vs_var_dR_{year}", "", 800, 600)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    c.SetMargin(0.13, 0.04, 0.14, 0.08)
    c.SetTickx()
    c.SetTicky()
    c.cd()  # NEW: ensure this canvas is the current pad

    _keepalive = []

    frame = ROOT.TH1F(f"frame_sumw_vs_var_dR_{year}", "", 1, 0.0, 4.0)
    frame.SetDirectory(0)
    _keepalive.append(frame)
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("#DeltaR(#gamma_{1}, #gamma_{2}) Cut")
    frame.GetYaxis().SetTitle(scenario_label)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetYaxis().SetTitleOffset(1.15)
    frame.GetXaxis().SetTitleSize(0.055)
    frame.GetYaxis().SetTitleSize(0.055)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelOffset(0.01)
    frame.SetMaximum(1.1)
    frame.SetMinimum(0.0)
    frame.Draw()  # CHANGED: 用預設 Draw()，避免只畫 AXIS 導致 grid 不更新/不顯示

    # NEW: grid 要設在 pad 上，且在 Draw 之後做最穩
    if ROOT.gPad:
        ROOT.gPad.SetGridx(1)
        ROOT.gPad.SetGridy(1)

        # 某些 ROOT 版本的 TVirtualPad 沒有 SetGridColor/Style/Width；做相容保護
        try:
            if hasattr(ROOT.gPad, "SetGridColor"):
                ROOT.gPad.SetGridColor(ROOT.kGray + 1)
            if hasattr(ROOT.gPad, "SetGridStyle"):
                ROOT.gPad.SetGridStyle(3)
            if hasattr(ROOT.gPad, "SetGridWidth"):
                ROOT.gPad.SetGridWidth(1)
        except Exception:
            pass

    # # ------------------------------------------------------------------
    # # 3 gradients × 4 colors + 2 neutral = 14 colors (publication-safe)
    # # ------------------------------------------------------------------
    #     palette = [
    #     # Cold (blue → green)
    #     ROOT.TColor.GetColor("#A6CEE3"),
    #     ROOT.TColor.GetColor("#1F78B4"),
    #     ROOT.TColor.GetColor("#33A02C"),
    #     ROOT.TColor.GetColor("#0B5D3A"),

    #     # Warm (orange → red)
    #     ROOT.TColor.GetColor("#FDBF6F"),
    #     ROOT.TColor.GetColor("#FF7F00"),
    #     ROOT.TColor.GetColor("#E31A1C"),
    #     ROOT.TColor.GetColor("#99000D"),

    #     # Neutral (purple → blue-gray)
    #     ROOT.TColor.GetColor("#CAB2D6"),
    #     ROOT.TColor.GetColor("#6A3D9A"),
    #     ROOT.TColor.GetColor("#5A6E8C"),
    #     ROOT.TColor.GetColor("#2F3E55"),

    #     # Neutral dark (gray → black)
    #     ROOT.TColor.GetColor("#BDBDBD"),
    #     ROOT.TColor.GetColor("#7F7F7F"),
    #     ROOT.TColor.GetColor("#4A4A4A"),
    #     ROOT.TColor.GetColor("#1A1A1A"),
    # ]

    palette = [
        # Cold & dark (背景用)
        ROOT.TColor.GetColor("#0B3C5D"),  # deep navy
        ROOT.TColor.GetColor("#1F618D"),
        ROOT.TColor.GetColor("#2874A6"),
        ROOT.TColor.GetColor("#3498DB"),
        ROOT.TColor.GetColor("#5DADE2"),

        # Green → cyan
        ROOT.TColor.GetColor("#117864"),
        ROOT.TColor.GetColor("#17A589"),
        ROOT.TColor.GetColor("#48C9B0"),

        # Purple / pink
        ROOT.TColor.GetColor("#7D3C98"),
        ROOT.TColor.GetColor("#C0392B"),  # dark red jump
        ROOT.TColor.GetColor("#E74C3C"),

        # Nuclear highlights (後畫王者)
        ROOT.TColor.GetColor("#FF6F00"),  # neon orange
        ROOT.TColor.GetColor("#FF1744"),  # hot red
        ROOT.TColor.GetColor("#000000"),  # final boss
    ]



    line_syles = [1, 2, 3]
    marker_styles = [20,  21,  23,  33,  34,  47,   29]
    marker_size   = [1.0, 0.9, 1.0, 1.5, 1.3, 1.2, 1.5]

    graphs: List[Tuple[int, ROOT.TGraph]] = []
    for i, ma in enumerate(mas):
        xy = by_ma.get(ma) or {}
        xs = sorted(xy.keys())
        if len(xs) < 1:
            continue
        ys = [float(xy[x]) for x in xs]
        g = _g_from_xy([0]*len(xs), ys)  # placeholder, will overwrite X below
        # overwrite X with float cuts
        from array import array as carray
        xarr = carray("d", [float(x) for x in xs])
        yarr = carray("d", [float(y) for y in ys])
        g = ROOT.TGraph(len(xs), xarr, yarr)

        col = palette[i % len(palette)]
        g.SetLineColor(col)
        g.SetMarkerColor(col)
        g.SetLineWidth(3)
        g.SetLineStyle(3)
        g.SetLineStyle(line_syles[i % len(line_syles)])
        g.SetMarkerStyle(marker_styles[i % len(marker_styles)])
        g.SetMarkerSize(marker_size[i % len(marker_size)])
        g.Draw("LP SAME")
        graphs.append((ma, g))

        _keepalive.append(g)

    # legend
    leg = ROOT.TLegend(0.745, 0.24, 0.99, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.033)
    leg.AddEntry("", f"{year}", "")
    for ma, g in graphs:
        leg.AddEntry(g, f"m_{{a}} = {ma} GeV", "lp")
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
    out_pdf = out_dir / (out_name or f"efficiency_vs_var_dR_g1g2_{year}.pdf")  # CHANGED
    c.Modified()
    c.Update()
    if ROOT.gPad:
        ROOT.gPad.RedrawAxis()  # NEW: 確保軸/格線最終刷新在正確層
    c.SaveAs(str(out_pdf))

    try:
        prim = c.GetListOfPrimitives()
        if prim:
            prim.Remove(frame)
    except Exception:
        pass
    try:
        c.Close()
    except Exception:
        pass
    del frame, c
    _ = _keepalive

def main():
    parser = argparse.ArgumentParser(description="Plot PHID event efficiency vs mA per year.")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022preEE). Use 'all' to plot all years.")
    args = parser.parse_args()

    in_dir = Path(baseDir)
    out_dir = Path(outDir)

    # NEW: load bkg counts once (year -> bkg event count)
    global _BKG_EVENT_BY_YEAR
    _BKG_EVENT_BY_YEAR = _load_bkg_event_count_by_year(in_dir)


    # year -> ma -> scenario -> {cut:eff}
    store: Dict[str, Dict[int, Dict[str, Dict[str, float]]]] = {}

    # NEW: efficiency vs var_dR_g1g2 cut (0~3), 14 lines (mA) per year
    eff_scan = _load_eff_scan_from_root_by_year_ma(
        Path(baseDir),
        years=YEAR_ORDER if args.year == "all" else [args.year],
        ma_dirs=mAs,
        tree_name="inclusive",
        dr_branch="var_dR_g1g2",
        dr_min=0.0,
        dr_max=4.0,
        dr_step=0.05,
    )

    if args.year != "all":
        by_ma = eff_scan.get(args.year)
        if by_ma:
            _plot_sumw_vs_var_dR_year(
                args.year,
                by_ma,
                out_dir,
                scenario_label="Efficiency / 0.05",
                out_name=f"EffVdRg1g2_{args.year}.pdf",
            )
        return

    for year in YEAR_ORDER:
        by_ma = eff_scan.get(year)
        if by_ma:
            _plot_sumw_vs_var_dR_year(
                year,
                by_ma,
                out_dir,
                scenario_label="Efficiency / 0.05",
                out_name=f"EffVdRg1g2_{year}.pdf",
            )

if __name__ == "__main__":
    main()
