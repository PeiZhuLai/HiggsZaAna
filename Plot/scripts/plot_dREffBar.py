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

baseDir = "/eos/home-p/pelai/HZa/root_P2Root/run3"

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/dREffBar"

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

def _eff_at_cut(xy: Dict[float, float], cut: float) -> Optional[float]:
    """Get eff(dr > cut) from scan dict; choose nearest available cut key."""
    if not xy:
        return None
    try:
        keys = list(xy.keys())
        k = min(keys, key=lambda a: abs(float(a) - float(cut)))
        return float(xy[k])
    except Exception:
        return None

def _fractions_from_eff_scan(xy: Dict[float, float]) -> Optional[Tuple[float, float, float]]:
    """
    Convert eff(dr > x) into fractions in 3 regions:
      f1 = P(dR < 0.1)
      f2 = P(0.1 < dR < 0.3)
      f3 = P(dR > 0.3)
    """
    e01 = _eff_at_cut(xy, 0.1)
    e03 = _eff_at_cut(xy, 0.3)
    if e01 is None or e03 is None:
        return None
    # clamp for safety
    e01 = max(0.0, min(1.0, e01))
    e03 = max(0.0, min(1.0, e03))
    f1 = 1.0 - e01
    f3 = e03
    f2 = e01 - e03
    # numeric cleanup
    f1 = max(0.0, min(1.0, f1))
    f2 = max(0.0, min(1.0, f2))
    f3 = max(0.0, min(1.0, f3))
    s = f1 + f2 + f3
    if s > 0:
        f1, f2, f3 = f1 / s, f2 / s, f3 / s
    return f1, f2, f3

def _plot_dr_fraction_bar_year(
    year: str,
    by_ma: Dict[int, Dict[float, float]],
    out_dir: Path,
    *,
    out_name: Optional[str] = None,
) -> None:
    mas = sorted(by_ma.keys())
    if not mas:
        return

    # Build fractions table
    rows: List[Tuple[int, float, float, float]] = []
    for ma in mas:
        fr = _fractions_from_eff_scan(by_ma.get(ma) or {})
        if fr is None:
            continue
        f1, f2, f3 = fr
        rows.append((ma, f1, f2, f3))
    if not rows:
        return

    n = len(rows)

    c = ROOT.TCanvas(f"c_drFracBar_{year}", "", 800, 800)
    ROOT.gStyle.SetOptStat(0)
    # 留更大的 top margin 讓 legend 放在「圖筐外」(margin 區)
    c.SetMargin(0.13, 0.06, 0.12, 0.10) #left, right, bottom, top
    c.SetTickx()
    c.SetTicky()
    c.cd()

    # frame: x in percent, y is category index
    frame = ROOT.TH2F(f"frame_drFracBar_{year}", "", 1, 0.0, 100.0, n, 0.0, float(n))
    frame.SetDirectory(0)
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("Fraction (%)")
    frame.GetYaxis().SetTitle("m_{a} (GeV)")
    frame.GetXaxis().SetTitleOffset(1.05)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetLabelSize(0.060)

    # y labels as mA
    for i, (ma, _, _, _) in enumerate(rows, start=1):
        frame.GetYaxis().SetBinLabel(i, f"{ma}")
    frame.Draw("")

    if ROOT.gPad:
        ROOT.gPad.SetGridx(1)
        ROOT.gPad.SetGridy(0)

    # colors (3 categories)
    col1 = _root_color("#FF688B", fallback=ROOT.kBlue + 1)   # dR < 0.1
    col2 = _root_color("#FFD557", fallback=ROOT.kGreen + 2)  # 0.1 < dR < 0.3
    col3 = _root_color("#9CD1FF", fallback=ROOT.kRed + 1)    # dR > 0.3

    # 水平 legend：放在上方圖筐外 (top margin 區域)，用 NDC 座標
    leg = ROOT.TLegend(0.12, 0.90, 0.94, 0.95)  # (x1,y1,x2,y2) in NDC
    leg.SetNColumns(3)  # 水平排列
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    # dummy boxes for legend
    b1 = ROOT.TBox(0, 0, 0, 0); b1.SetFillColor(col1); b1.SetLineColor(col1)
    b2 = ROOT.TBox(0, 0, 0, 0); b2.SetFillColor(col2); b2.SetLineColor(col2)
    b3 = ROOT.TBox(0, 0, 0, 0); b3.SetFillColor(col3); b3.SetLineColor(col3)
    leg.AddEntry(b1, "#DeltaR(#gamma,#gamma) < 0.1", "f")
    leg.AddEntry(b2, "0.1 < #DeltaR(#gamma,#gamma) < 0.3", "f")
    leg.AddEntry(b3, "#DeltaR(#gamma,#gamma) > 0.3", "f")

    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextSize(0.032)

    # draw stacked horizontal bars with TBox
    _keep = [frame, b1, b2, b3, leg]

    bar_h = 0.70  # fraction of bin height
    for idx, (ma, f1, f2, f3) in enumerate(rows):
        y0 = float(idx) + (1.0 - bar_h) / 2.0
        y1 = float(idx) + 1.0 - (1.0 - bar_h) / 2.0

        p1 = 100.0 * f1
        p2 = 100.0 * f2
        p3 = 100.0 * f3

        x0 = 0.0
        x1 = x0 + p1
        x2 = x1 + p2
        x3 = x2 + p3

        bx1 = ROOT.TBox(x0, y0, x1, y1); bx1.SetFillColor(col1); bx1.SetLineColor(col1)
        bx2 = ROOT.TBox(x1, y0, x2, y1); bx2.SetFillColor(col2); bx2.SetLineColor(col2)
        bx3 = ROOT.TBox(x2, y0, x3, y1); bx3.SetFillColor(col3); bx3.SetLineColor(col3)
        bx1.Draw("SAME"); bx2.Draw("SAME"); bx3.Draw("SAME")
        _keep += [bx1, bx2, bx3]

        # percent labels on bars (only if segment is wide enough)
        def _draw_pct(xl: float, xr: float, yy0: float, yy1: float, val: float):
            w = xr - xl
            if w < 6.0:  # too narrow, skip to avoid clutter
                return
            latex.SetTextAlign(22)  # center
            latex.DrawLatex((xl + xr) / 2.0, (yy0 + yy1) / 2.0 - 0.02, f"{val:.1f}%")

        _draw_pct(x0, x1, y0, y1, p1)
        _draw_pct(x1, x2, y0, y1, p2)
        _draw_pct(x2, x3, y0, y1, p3)

    # CMS label（往上移到 top margin 區，避免和 legend 重疊）
    lat2 = ROOT.TLatex()
    lat2.SetNDC()
    lat2.SetTextFont(42)
    lat2.SetTextSize(0.045)
    lat2.DrawLatex(c.GetLeftMargin(), 0.955, "#bf{CMS} #it{Simulation}")
    _keep.append(lat2)

    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetTextSize(0.040)
        lat_lumi.DrawLatex(0.94, 0.955, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")
        _keep.append(lat_lumi)

    leg.Draw()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = out_dir / (out_name or f"dR_fraction_bar_{year}.pdf")
    c.Modified()
    c.Update()
    if ROOT.gPad:
        ROOT.gPad.RedrawAxis()
    c.SaveAs(str(out_pdf))

    try:
        c.Close()
    except Exception:
        pass
    _ = _keep

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
            _plot_dr_fraction_bar_year(
                args.year,
                by_ma,
                out_dir,
                out_name=f"dRFracBar_{args.year}.pdf",
            )
        return

    for year in YEAR_ORDER:
        by_ma = eff_scan.get(year)
        if by_ma:
            _plot_dr_fraction_bar_year(
                year,
                by_ma,
                out_dir,
                out_name=f"dRFracBar_{year}.pdf",
            )

if __name__ == "__main__":
    main()
