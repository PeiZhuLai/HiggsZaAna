import os
from pathlib import Path
import sys
import ROOT
import numpy as np
import re
from typing import Dict, List, Tuple, Optional
import uproot
import awkward as ak
import argparse  # CLI 參數

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

baseDir = "/eos/home-p/pelai/HZa/root_P2Root/run3"

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/sigGenInfo"

# Define the list of mA values
mAs = ["mA_M1", "mA_M2", "mA_M3", "mA_M4", "mA_M5", "mA_M6", "mA_M7", "mA_M8", "mA_M9", "mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]

# NEW: split MA choices
ALL_MA_DIRS = mAs
PT_MA_DIRS = ["mA_M1", "mA_M5", "mA_M10", "mA_M20", "mA_M30"]

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
# YEAR_ORDER = ["2022preEE"]

# NEW: run3 組成用的子年份清單
RUN3 = YEAR_ORDER[:]  # ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, '2024':108.95,
            'Run3':170.84 }


def _lumi_fb_for_year(year: str) -> Optional[float]:
    v = lumiMap.get(year)
    return float(v) if v is not None else None

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

def _safe_np(a) -> np.ndarray:
    try:
        return ak.to_numpy(a)
    except Exception:
        return np.asarray([])

def _load_arrays_for_year_mas(
    base_dir: Path,
    *,
    year: str,
    ma_dirs: List[str],
    branches: List[str],
    tree_name: str = "inclusive",
) -> Dict[str, np.ndarray]:
    """
    讀取指定 year + 多個 mA 檔案，將 branches 串接起來（缺 branch/檔就跳過）。
    回傳: branch -> numpy array (flatten 後)
    """
    out: Dict[str, List[np.ndarray]] = {b: [] for b in branches}

    for ma_tag in ma_dirs:
        ma_path = base_dir / ma_tag
        if not ma_path.is_dir():
            continue
        fpath = ma_path / f"{year}.root"
        if not fpath.exists():
            continue

        try:
            with uproot.open(str(fpath)) as f:
                t = f[tree_name]
                keys = set(t.keys())
                want = [b for b in branches if b in keys]
                if not want:
                    continue
                arrs = t.arrays(want, library="ak")
        except Exception:
            continue

        for b in branches:
            if b not in arrs.fields:
                continue
            a = arrs[b]
            # 以 "event-level scalar" 為主：flatten 全部
            a = ak.flatten(a, axis=None)
            anp = _safe_np(a)
            if anp.size:
                out[b].append(anp)

    merged: Dict[str, np.ndarray] = {}
    for b, chunks in out.items():
        if chunks:
            merged[b] = np.concatenate(chunks, axis=0)
    return merged

def _load_arrays_for_year_ma(
    base_dir: Path,
    *,
    year: str,
    ma_dir: str,
    branches: List[str],
    tree_name: str = "inclusive",
) -> Dict[str, np.ndarray]:
    """
    讀取指定 year + 單一 mA 檔案。
    回傳: branch -> numpy array (flatten 後)
    """
    ma_path = base_dir / ma_dir
    if not ma_path.is_dir():
        return {}
    fpath = ma_path / f"{year}.root"
    if not fpath.exists():
        return {}

    try:
        with uproot.open(str(fpath)) as f:
            t = f[tree_name]
            keys = set(t.keys())
            want = [b for b in branches if b in keys]
            if not want:
                return {}
            arrs = t.arrays(want, library="ak")
    except Exception:
        return {}

    out: Dict[str, np.ndarray] = {}
    for b in want:
        a = ak.flatten(arrs[b], axis=None)
        anp = _safe_np(a)
        if anp.size:
            out[b] = anp
    return out

# NEW: 合併多個年份（同一個 mA）成單一 arrays dict
def _merge_arrays_dicts(dicts: List[Dict[str, np.ndarray]]) -> Dict[str, np.ndarray]:
    merged: Dict[str, List[np.ndarray]] = {}
    for d in dicts:
        for k, v in d.items():
            if v is None or getattr(v, "size", 0) == 0:
                continue
            merged.setdefault(k, []).append(v)
    out: Dict[str, np.ndarray] = {}
    for k, chunks in merged.items():
        if chunks:
            out[k] = np.concatenate(chunks, axis=0)
    return out

# NEW: 讀取多年份 + 多 mA，將 branches 串接起來（缺檔就跳過）
def _load_arrays_for_years_mas(
    base_dir: Path,
    *,
    years: List[str],
    ma_dirs: List[str],
    branches: List[str],
    tree_name: str = "inclusive",
) -> Dict[str, np.ndarray]:
    per_year = [
        _load_arrays_for_year_mas(base_dir, year=y, ma_dirs=ma_dirs, branches=branches, tree_name=tree_name)
        for y in years
    ]
    return _merge_arrays_dicts(per_year)

# NEW: 讀取多年份 + 單一 mA
def _load_arrays_for_years_ma(
    base_dir: Path,
    *,
    years: List[str],
    ma_dir: str,
    branches: List[str],
    tree_name: str = "inclusive",
) -> Dict[str, np.ndarray]:
    per_year = [
        _load_arrays_for_year_ma(base_dir, year=y, ma_dir=ma_dir, branches=branches, tree_name=tree_name)
        for y in years
    ]
    return _merge_arrays_dicts(per_year)

def _ma_value_from_dir(ma_dir: str) -> Optional[int]:
    m = re.match(r"mA_M(\d+)$", str(ma_dir))
    return int(m.group(1)) if m else None

def _colors_many(n: int) -> List[int]:
    # 簡單的可辨識 palette；不夠就循環
    # hexes = [
    #     "#1B4F72", "#2874A6", "#3498DB", "#5DADE2",
    #     "#117864", "#17A589", "#48C9B0",
    #     "#7D3C98", "#8E44AD", "#A569BD",
    #     "#AF601A", "#D68910", "#E67E22",
    #     "#7F7F7F",
    # ]

    hexes = [
        # Cold & dark (背景用)
        "#0B3C5D",  # deep navy
        "#1F618D",
        "#2874A6",
        "#3498DB",
        "#5DADE2",

        # Green → cyan
        "#117864",
        "#17A589",
        "#48C9B0",

        # Purple / pink
        "#7D3C98",
        "#C0392B",  # dark red jump
        "#E74C3C",

        # Nuclear highlights (後畫王者)
        "#FF6F00",  # neon orange
        "#FF1744",  # hot red
        "#000000",  # final boss
    ]

    cols = [_root_color(hx, fallback=ROOT.kBlack) for hx in hexes]
    if n <= len(cols):
        return cols[:n]
    # 不夠就循環
    return [cols[i % len(cols)] for i in range(n)]

def _markers_many(n: int) -> List[int]:
    ms = [20, 21, 22, 23, 24, 25, 26, 27, 28, 33, 34]
    if n <= len(ms):
        return ms[:n]
    return [ms[i % len(ms)] for i in range(n)]

def _linestyles_many(n: int) -> List[int]:
    ls = [1,2,3]
    if n <= len(ls):
        return ls[:n]
    return [ls[i % len(ls)] for i in range(n)]

def _make_hist_from_np(
    name: str,
    x: np.ndarray,
    w: Optional[np.ndarray],
    *,
    nbins: int,
    xlow: float,
    xhigh: float,
) -> ROOT.TH1F:
    h = ROOT.TH1F(name, "", nbins, xlow, xhigh)
    h.SetDirectory(0)

    if x is None or x.size == 0:
        return h

    # NEW: 對齊 x / w 長度，避免不同 branch 合併後長度不一致造成 boolean mask mismatch
    if w is not None and getattr(w, "size", 0):
        nmin = min(int(x.size), int(w.size))
        if nmin <= 0:
            return h
        x = x[:nmin]
        w = w[:nmin]

    # NEW: finite mask（x 必須 finite；w 若存在也必須 finite）
    m = np.isfinite(x) & (x >= xlow) & (x <= xhigh)
    if w is not None and getattr(w, "size", 0):
        m = m & np.isfinite(w)

    x = x[m]
    if w is not None and getattr(w, "size", 0):
        w = w[m]
    else:
        w = None

    bins = np.linspace(xlow, xhigh, nbins + 1, dtype=float)
    counts, _ = np.histogram(x, bins=bins, weights=w)

    for i in range(nbins):
        h.SetBinContent(i + 1, float(counts[i]))
    return h

def _root_style_line(h, color: int, mstyle: int, lstyle: int = 1) -> None:
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetLineWidth(3)
    h.SetLineStyle(lstyle)  # NEW
    h.SetMarkerStyle(mstyle)
    h.SetMarkerSize(0.9)

def _plot_overlay_hists(
    *,
    year: str,
    title: str,
    xlabel: str,
    ylabel: str,
    out_path: Path,
    hists: List[Tuple[str, ROOT.TH1F]],
    normalize: bool = False,
    legend_map: Optional[Dict[str, str]] = None,  # NEW: map branch->label for legend
    subtitle_lines: Optional[List[str]] = None,    # NEW: extra lines (e.g. mA info)
    norm_left_margin: Optional[float] = None,       # NEW
    norm_y_title_offset: Optional[float] = None,    # NEW
    cms_label_x: Optional[float] = None,            # NEW: CMS label x (NDC)
    cms_label_y: Optional[float] = None,            # NEW: CMS label y (NDC)
) -> None:
    if not hists:
        return

    c = ROOT.TCanvas(f"c_{out_path.stem}_{year}", "", 800, 600)

    # CHANGED: normalize 圖可獨立調整 left margin（不影響非 normalize）
    left, right, bottom, top = 0.13, 0.04, 0.14, 0.08
    if normalize and norm_left_margin is not None:
        left = float(norm_left_margin)
    c.SetMargin(left, right, bottom, top)

    c.SetTickx()
    c.SetTicky()
    c.cd()

    # y-max
    ymax = 0.0
    for _, h in hists:
        if normalize and h.Integral() > 0:
            h.Scale(1.0 / h.Integral())
        ymax = max(ymax, h.GetMaximum())

    frame = ROOT.TH1F(
        f"frame_{out_path.stem}_{year}",
        "",
        1,
        hists[0][1].GetXaxis().GetXmin(),
        hists[0][1].GetXaxis().GetXmax(),
    )
    frame.SetDirectory(0)
    frame.SetTitle("")
    frame.GetXaxis().SetTitle(xlabel)
    frame.GetXaxis().SetTitleOffset(1.1)

    # CHANGED: normalize 圖可獨立調整 y axis title offset
    ytitle_off = 1.15
    if normalize and norm_y_title_offset is not None:
        ytitle_off = float(norm_y_title_offset)
    frame.GetYaxis().SetTitleOffset(ytitle_off)

    frame.GetXaxis().SetTitleSize(0.055)
    frame.GetYaxis().SetTitleSize(0.055)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelOffset(0.01)

    # CHANGED: y-title append bin width
    xax = hists[0][1].GetXaxis()
    b1 = max(1, min(xax.GetNbins(), xax.FindBin(0.5 * (xax.GetXmin() + xax.GetXmax()))))
    bw = float(xax.GetBinWidth(b1)) if xax.GetNbins() > 0 else 1.0
    if bw <= 0:
        bw = 1.0
    if "GeV" in xlabel:
        frame.GetYaxis().SetTitle(f"{ylabel} / {bw:.2f} GeV")
    else:
        frame.GetYaxis().SetTitle(f"{ylabel} / {bw:.2f}")

    frame.SetMaximum(1.25 * ymax if ymax > 0 else 1.0)
    frame.SetMinimum(0.0)
    frame.Draw()

    # if ROOT.gPad:
    #     ROOT.gPad.SetGridx(1)
    #     ROOT.gPad.SetGridy(1)

    # --- CHANGED: legend sizing in one place, auto-scale by number of entries ---
    do_legend = len(hists) >= 2

    # fixed header lines always shown in legend
    header_lines = [f"{year}", r"H#rightarrowZa#rightarrowll#gamma#gamma"]

    # count lines (header + optionally entries)
    n_lines = len(header_lines) + (len(hists) if do_legend else 0)

    # geometry: keep same top-right anchor, vary height with lines
    x1, y1 = 0.98, 0.88
    x0 = 0.68
    line_h = 0.042     # per-line height in NDC (tune if needed)
    pad = 0.018        # extra padding
    h_leg = pad + n_lines * line_h
    y0 = max(0.14, y1 - h_leg)

    leg = ROOT.TLegend(x0, y0, x1, y1)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.039)

    for s in header_lines:
        leg.AddEntry("", s, "")

    for i, (lab, h) in enumerate(hists):
        h.Draw("HIST SAME")
        if do_legend:
            pretty = (legend_map or {}).get(lab, lab)
            leg.AddEntry(h, pretty, "l")

    # draw legend always (even single hist -> header only), size auto-adjusted
    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)

    # NEW: CMS Preliminary 位置跟著 left margin（normalize 調過 margin 時自動右移）
    cms_x = 0.13 if cms_label_x is None else float(cms_label_x)
    cms_y = 0.93 if cms_label_y is None else float(cms_label_y)
    if cms_label_x is None and normalize and norm_left_margin is not None:
        cms_x = float(norm_left_margin)

    lat.DrawLatex(cms_x, cms_y, "#bf{CMS} #it{Preliminary}")

    # NEW: 右上角次標題（例如 mA 資訊）
    if subtitle_lines:
        lat_sub = ROOT.TLatex()
        lat_sub.SetNDC()
        lat_sub.SetTextFont(42)
        lat_sub.SetTextAlign(13)   # right aligned
        lat_sub.SetTextSize(0.04)
        x0, y0, dy = 0.18, 0.88, 0.042
        for i, line in enumerate(subtitle_lines):
            if line:
                lat_sub.DrawLatex(x0, y0 - i * dy, line)

    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetTextSize(0.040)
        lat_lumi.DrawLatex(0.96, 0.93, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    c.Modified()
    c.Update()
    if ROOT.gPad:
        ROOT.gPad.RedrawAxis()
    c.SaveAs(str(out_path))
    try:
        c.Close()
    except Exception:
        pass
    del frame, c

def _plot_gen_distributions_for_year(
    *,
    base_dir: Path,
    out_dir: Path,
    year: str,
    ma_dirs: List[str],
    tree_name: str = "inclusive",
    ma_dirs_pt: Optional[List[str]] = None,  # NEW: pt-only MA selection
) -> None:
    # 你指定要畫的分布
    single_vars = [
        ("GenHzaHiggs_pt", "Gen Higgs P_{T} [GeV]", 120, 0.0, 300.0),
        ("GenHzaZ_pt",     "Gen Z P_{T} [GeV]",     100, 0.0, 250.0),
        ("GenHzaALP_pt",   "Gen a P_{T} [GeV]",     60, 0.0, 150.0),
        ("GenHza_dR_ZALP", "#DeltaR(Z, ALP)",       60, 0.0, 6.0),
        ("GenHza_dR_ll",   "#DeltaR(l, l)",         60, 0.0, 6.0),
        ("GenALP_dR_gg",   "#DeltaR(#gamma, #gamma)",100, 0.0, 3.0),
    ]

    legend_map={
        "GenALPLeadPho_pt": "Lead", 
        "GenALPSubleadPho_pt": "Sublead", 
        "GenHzaZLeadLep_pt": "Lead", 
        "GenHzaZSubleadLep_pt": "Sublead"
    }

    # CHANGED: pt overlay 改成兩張圖（Lead/Sublead），各自疊 mA 線
    overlay_pairs = [
        (["GenALPLeadPho_pt"],    "Lead Photon P_{T} [GeV]",    50, 0.0, 50.0, "GenALP_leadPho_pt_overlay_byMA"),
        (["GenALPSubleadPho_pt"], "Sublead Photon P_{T} [GeV]", 50, 0.0, 50.0, "GenALP_subleadPho_pt_overlay_byMA"),
    ]

    # CHANGED: 原本 lep_overlay 連續賦值會被覆蓋；改成 lead/sublead 各自的 spec（含對應 pdg branch）
    lep_pt_specs = [
        {
            "pt_branch": "GenHzaZLeadLep_pt",
            "pdg_branch": "GenHzaZLeadLep_pdgId",
            "xlabel": "Lead Lepton P_{T} [GeV]",
            "nb": 80,
            "xl": 0.0,
            "xh": 160.0,
            "tag_suffix": "leadLep",
        },
        {
            "pt_branch": "GenHzaZSubleadLep_pt",
            "pdg_branch": "GenHzaZSubleadLep_pdgId",
            "xlabel": "Sublead Lepton P_{T} [GeV]",
            "nb": 80,
            "xl": 0.0,
            "xh": 160.0,
            "tag_suffix": "subleadLep",
        },
    ]

    # 一次把需要的 branch 全列出來讀
    # CHANGED: pdg_branch 不再用單一變數，改把 lead/sublead 的 pdg branch 都加入
    need = {"weight"}
    for spec in lep_pt_specs:
        need.add(spec["pdg_branch"])
        need.add(spec["pt_branch"])
    for b, _, _, _, _ in single_vars:
        need.add(b)
    for bs, _, _, _, _, _ in overlay_pairs:
        need.update(bs)

    # --- CHANGED: general/pt arrays 讀取支援 Run3 ---
    is_combined = (year == "Run3")
    years_to_merge = RUN3 if is_combined else [year]

    ma_dirs_pt = ma_dirs_pt or ma_dirs

    arr = _load_arrays_for_years_mas(
        base_dir,
        years=years_to_merge,
        ma_dirs=ma_dirs,
        branches=sorted(need),
        tree_name=tree_name,
    )
    arr_pt = _load_arrays_for_years_mas(
        base_dir,
        years=years_to_merge,
        ma_dirs=ma_dirs_pt,
        branches=sorted(need),
        tree_name=tree_name,
    )

    w = arr.get("weight")
    w_pt = arr_pt.get("weight")

    # NEW: mA subtitle text separately (all vs pt)
    def _ma_text_from_dirs(dirs: List[str]) -> str:
        vals = sorted({int(m.group(1)) for d in dirs if (m := re.match(r"mA_M(\d+)$", str(d)))})
        if len(vals) == 1:
            return f"m_{{a}} = {vals[0]} GeV"
        if len(vals) >= 2:
            return f"m_{{a}} = {vals[0]}-{vals[-1]} GeV"
        return ""

    ma_text = _ma_text_from_dirs(ma_dirs)
    ma_text_pt = _ma_text_from_dirs(ma_dirs_pt)

    # colors/markers
    cols = [
        _root_color("#E31A1C", fallback=ROOT.kRed+1),
        _root_color("#1F78B4", fallback=ROOT.kBlue+1),
        _root_color("#33A02C", fallback=ROOT.kGreen+2),
        _root_color("#FF7F00", fallback=ROOT.kOrange+7),
    ]
    mstyles = [20, 21, 23, 33]

    # NEW: for "overlay by MA" plots
    cols_many = _colors_many(len(ma_dirs_pt))
    mstyles_many = _markers_many(len(ma_dirs_pt))
    ls_many = _linestyles_many(len(ma_dirs_pt))  # NEW

    # NEW: for single_vars overlay across ALL mA (14 lines)
    cols_many_all = _colors_many(len(ma_dirs))
    mstyles_many_all = _markers_many(len(ma_dirs))
    ls_many_all = _linestyles_many(len(ma_dirs))  # NEW

    # --- CHANGED: 單變數改成「每個 mA 一條線」(總共 14 條) ---
    for (b, xlabel, nb, xl, xh) in single_vars:
        hs = []
        legend_map_by_ma: Dict[str, str] = {}

        for i, ma_tag in enumerate(ma_dirs):
            arr_ma = _load_arrays_for_years_ma(
                base_dir,
                years=years_to_merge,   # CHANGED
                ma_dir=ma_tag,
                branches=["weight", b],
                tree_name=tree_name,
            )
            x = arr_ma.get(b)
            w_ma = arr_ma.get("weight")
            if x is None or x.size == 0:
                continue

            ma_val = _ma_value_from_dir(ma_tag)
            key = ma_tag
            legend_map_by_ma[key] = (f"m_{{a}} = {ma_val} GeV" if ma_val is not None else str(ma_tag))

            h = _make_hist_from_np(f"h_{b}_{year}_{ma_tag}", x, w_ma, nbins=nb, xlow=xl, xhigh=xh)

            # guard: palette 長度不足時循環使用
            ci = cols_many_all[i % len(cols_many_all)] if cols_many_all else ROOT.kBlack
            mi = mstyles_many_all[i % len(mstyles_many_all)] if mstyles_many_all else 20
            li = ls_many_all[i % len(ls_many_all)] if ls_many_all else 1
            _root_style_line(h, ci, mi, li)

            hs.append((key, h))

        if hs:
            # 原本：未歸一化
            _plot_overlay_hists(
                year=year,
                title=b,
                xlabel=xlabel,
                ylabel="Events",
                out_path=out_dir / year / f"{b}_{year}.pdf",
                hists=hs,
                normalize=False,
                legend_map=legend_map_by_ma,
                subtitle_lines=None,
            )
            # NEW：歸一化到 1
            _plot_overlay_hists(
                year=year,
                title=b,
                xlabel=xlabel,
                ylabel="Events",
                out_path=out_dir / year / f"{b}_{year}_norm.pdf",
                hists=hs,
                normalize=True,
                legend_map=legend_map_by_ma,
                subtitle_lines=None,
                norm_left_margin=0.15,      
                norm_y_title_offset=1.25,
            )

    # ------------------------------------------------------------------
    # CHANGED: photon lead/sublead pT 一律畫 14 條線（ALL mA），不再只用 pt subset
    # ------------------------------------------------------------------
    for (bs, xlabel, nb, xl, xh, tag) in overlay_pairs:
        for b in bs:
            hs = []
            legend_map_by_ma: Dict[str, str] = {}

            for i, ma_tag in enumerate(ma_dirs):  # CHANGED: ma_dirs (ALL)
                arr_ma = _load_arrays_for_years_ma(
                    base_dir,
                    years=years_to_merge,  # CHANGED
                    ma_dir=ma_tag,
                    branches=["weight", b],
                    tree_name=tree_name,
                )
                x = arr_ma.get(b)
                w_ma = arr_ma.get("weight")
                if x is None or x.size == 0:
                    continue

                ma_val = _ma_value_from_dir(ma_tag)
                key = ma_tag
                legend_map_by_ma[key] = (f"m_{{a}} = {ma_val} GeV" if ma_val is not None else str(ma_tag))

                h = _make_hist_from_np(f"h_{b}_{year}_{ma_tag}", x, w_ma, nbins=nb, xlow=xl, xhigh=xh)
                _root_style_line(h, cols_many_all[i], mstyles_many_all[i], ls_many_all[i])  # CHANGED: all palette
                hs.append((key, h))

            if hs:
                # 原本：未歸一化
                _plot_overlay_hists(
                    year=year,
                    title=tag,
                    xlabel=xlabel,
                    ylabel="Events",
                    out_path=out_dir / year / f"{tag}_{year}.pdf",
                    hists=hs,
                    normalize=False,
                    legend_map=legend_map_by_ma,
                    subtitle_lines=None,
                )
                # NEW：歸一化到 1
                _plot_overlay_hists(
                    year=year,
                    title=tag,
                    xlabel=xlabel,
                    ylabel="Events",
                    out_path=out_dir / year / f"{tag}_{year}_norm.pdf",
                    hists=hs,
                    normalize=True,
                    legend_map=legend_map_by_ma,
                    subtitle_lines=None,
                    norm_left_margin=0.16,      # NEW
                    norm_y_title_offset=1.45,   # NEW
                )

    # ------------------------------------------------------------------
    # CHANGED: electron/muon Lead/Sublead pT 各自畫 14 條線（ALL mA）
    # ------------------------------------------------------------------
    for spec in lep_pt_specs:
        pt_branch = spec["pt_branch"]
        pdg_branch = spec["pdg_branch"]
        nb, xl, xh = spec["nb"], spec["xl"], spec["xh"]
        tag_suffix = spec["tag_suffix"]

        lep_specs = [
            ("electron", 11, pt_branch, f"GenHza_e_{tag_suffix}_pt_overlay_byMA", f"{'Lead' if tag_suffix=='leadLep' else 'Sublead'} Electron P_{{T}} [GeV]"),
            ("muon",     13, pt_branch, f"GenHza_mu_{tag_suffix}_pt_overlay_byMA", f"{'Lead' if tag_suffix=='leadLep' else 'Sublead'} Muon P_{{T}} [GeV]"),
        ]

        for flav_name, abs_pdg, b_lep, out_tag, flav_xlabel in lep_specs:
            hs = []
            legend_map_by_ma: Dict[str, str] = {}

            for i, ma_tag in enumerate(ma_dirs):
                arr_ma = _load_arrays_for_years_ma(
                    base_dir,
                    years=years_to_merge,  # CHANGED
                    ma_dir=ma_tag,
                    branches=["weight", pdg_branch, b_lep],
                    tree_name=tree_name,
                )
                x = arr_ma.get(b_lep)
                w_ma = arr_ma.get("weight")
                pdg_ma = arr_ma.get(pdg_branch)
                if x is None or w_ma is None or pdg_ma is None:
                    continue
                if x.size == 0 or w_ma.size == 0 or pdg_ma.size == 0:
                    continue

                nmin = min(int(x.size), int(w_ma.size), int(pdg_ma.size))
                if nmin <= 0:
                    continue
                x = x[:nmin]
                w0 = w_ma[:nmin]
                pdg0 = pdg_ma[:nmin]

                absid = np.abs(pdg0.astype(np.int64, copy=False))
                msk = absid == int(abs_pdg)

                xx = x[msk]
                ww = w0[msk]
                if xx.size == 0:
                    continue

                ma_val = _ma_value_from_dir(ma_tag)
                key = ma_tag
                legend_map_by_ma[key] = (f"m_{{a}} = {ma_val} GeV" if ma_val is not None else str(ma_tag))

                h = _make_hist_from_np(
                    f"h_{b_lep}_{year}_{flav_name}_{ma_tag}",
                    xx,
                    ww,
                    nbins=nb,
                    xlow=xl,
                    xhigh=xh,
                )
                _root_style_line(h, cols_many_all[i], mstyles_many_all[i], ls_many_all[i])
                hs.append((key, h))

            if hs:
                _plot_overlay_hists(
                    year=year,
                    title=out_tag,
                    xlabel=flav_xlabel,  # CHANGED
                    ylabel="Events",
                    out_path=out_dir / year / f"{out_tag}_{year}.pdf",
                    hists=hs,
                    normalize=False,
                    legend_map=legend_map_by_ma,
                    subtitle_lines=None,
                )
                _plot_overlay_hists(
                    year=year,
                    title=out_tag,
                    xlabel=flav_xlabel,  # CHANGED
                    ylabel="Events",
                    out_path=out_dir / year / f"{out_tag}_{year}_norm.pdf",
                    hists=hs,
                    normalize=True,
                    legend_map=legend_map_by_ma,
                    subtitle_lines=None,
                    norm_left_margin=0.16,
                    norm_y_title_offset=1.45,
                )

    # ------------------------------------------------------------------
    # NEW: GenHza_lep_pt_overlay_byFlav_{electron,muon}_{year}
    #      只有兩條線：Lead lepton pT vs Sublead lepton pT
    #      且「所有 mA」與（year==Run3 時）「所有 era」都合併在一起
    # ------------------------------------------------------------------
    byflav_nb, byflav_xl, byflav_xh = 80, 0.0, 160.0
    byflav_specs = [
        {
            "flav": "electron",
            "abs_pdg": 11,
            "out_tag": "GenHza_lep_pt_overlay_byFlav_electron",
            "xlabel": "Gen e P_{T} [GeV]",
        },
        {
            "flav": "muon",
            "abs_pdg": 13,
            "out_tag": "GenHza_lep_pt_overlay_byFlav_muon",
            "xlabel": "Gen #mu P_{T} [GeV]",
        },
    ]

    lead_pt = arr.get("GenHzaZLeadLep_pt")
    sub_pt  = arr.get("GenHzaZSubleadLep_pt")
    lead_id = arr.get("GenHzaZLeadLep_pdgId")
    sub_id  = arr.get("GenHzaZSubleadLep_pdgId")
    w_all   = arr.get("weight")

    if (lead_pt is not None and sub_pt is not None and lead_id is not None and sub_id is not None and w_all is not None
        and lead_pt.size and sub_pt.size and lead_id.size and sub_id.size and w_all.size):

        def _align3(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
            n = min(int(a.size), int(b.size), int(c.size))
            return a[:n], b[:n], c[:n]

        # align each (pt,id,w) triplet separately to avoid mismatch
        lead_pt0, lead_id0, w0 = _align3(lead_pt, lead_id, w_all)
        sub_pt0,  sub_id0,  w1 = _align3(sub_pt,  sub_id,  w_all)

        for sp in byflav_specs:
            abs_pdg = int(sp["abs_pdg"])

            m_lead = (np.abs(lead_id0.astype(np.int64, copy=False)) == abs_pdg)
            m_sub  = (np.abs(sub_id0.astype(np.int64, copy=False))  == abs_pdg)

            x_lead, w_lead = lead_pt0[m_lead], w0[m_lead]
            x_sub,  w_sub  = sub_pt0[m_sub],  w1[m_sub]

            h_lead = _make_hist_from_np(
                f"h_{sp['out_tag']}_{year}_lead",
                x_lead, w_lead,
                nbins=byflav_nb, xlow=byflav_xl, xhigh=byflav_xh,
            )
            h_sub = _make_hist_from_np(
                f"h_{sp['out_tag']}_{year}_sublead",
                x_sub, w_sub,
                nbins=byflav_nb, xlow=byflav_xl, xhigh=byflav_xh,
            )

            # only two lines
            _root_style_line(h_lead, _root_color("#E31A1C", fallback=ROOT.kRed+1), 20, 1)
            _root_style_line(h_sub,  _root_color("#1F78B4", fallback=ROOT.kBlue+1), 21, 1)

            hs2 = [("lead", h_lead), ("sublead", h_sub)]
            leg2 = {"lead": "Lead", "sublead": "Sublead"}

            _plot_overlay_hists(
                year=year,
                title=sp["out_tag"],
                xlabel=sp["xlabel"],
                ylabel="Events",
                out_path=out_dir / year / f"{sp['out_tag']}_{year}.pdf",
                hists=hs2,
                normalize=False,
                legend_map=leg2,
                subtitle_lines=None,
            )
            _plot_overlay_hists(
                year=year,
                title=sp["out_tag"],
                xlabel=sp["xlabel"],
                ylabel="Events",
                out_path=out_dir / year / f"{sp['out_tag']}_{year}_norm.pdf",
                hists=hs2,
                normalize=True,
                legend_map=leg2,
                subtitle_lines=None,
                norm_left_margin=0.16,
                norm_y_title_offset=1.45,
            )

    # ------------------------------------------------------------------
    # NEW: GenHza_pho_pt_overlay_{year}.pdf
    #      只有兩條線：Lead photon pT vs Sublead photon pT
    #      且（year==Run3 時）「所有 mA + 所有 era」合併在一起
    # ------------------------------------------------------------------
    pho_nb, pho_xl, pho_xh = 50, 0.0, 50.0
    lead_pho = arr.get("GenALPLeadPho_pt")
    sub_pho  = arr.get("GenALPSubleadPho_pt")
    w_all_pho = arr.get("weight")

    if (lead_pho is not None and sub_pho is not None and w_all_pho is not None
        and lead_pho.size and sub_pho.size and w_all_pho.size):

        n0 = min(int(lead_pho.size), int(w_all_pho.size))
        n1 = min(int(sub_pho.size),  int(w_all_pho.size))

        x_lead, w_lead = lead_pho[:n0], w_all_pho[:n0]
        x_sub,  w_sub  = sub_pho[:n1],  w_all_pho[:n1]

        h_lead = _make_hist_from_np(
            f"h_GenHza_pho_pt_overlay_{year}_lead",
            x_lead, w_lead,
            nbins=pho_nb, xlow=pho_xl, xhigh=pho_xh,
        )
        h_sub = _make_hist_from_np(
            f"h_GenHza_pho_pt_overlay_{year}_sublead",
            x_sub, w_sub,
            nbins=pho_nb, xlow=pho_xl, xhigh=pho_xh,
        )

        _root_style_line(h_lead, _root_color("#E31A1C", fallback=ROOT.kRed+1), 20, 1)
        _root_style_line(h_sub,  _root_color("#1F78B4", fallback=ROOT.kBlue+1), 21, 1)

        hs_pho = [("lead", h_lead), ("sublead", h_sub)]
        leg_pho = {"lead": "Lead", "sublead": "Sublead"}

        out_tag = "GenHza_pho_pt_overlay"
        _plot_overlay_hists(
            year=year,
            title=out_tag,
            xlabel="Gen #gamma P_{T} [GeV]",
            ylabel="Events",
            out_path=out_dir / year / f"{out_tag}_{year}.pdf",
            hists=hs_pho,
            normalize=False,
            legend_map=leg_pho,
            subtitle_lines=None,
        )
        _plot_overlay_hists(
            year=year,
            title=out_tag,
            xlabel="Gen #gamma P_{T} [GeV]",
            ylabel="Events",
            out_path=out_dir / year / f"{out_tag}_{year}_norm.pdf",
            hists=hs_pho,
            normalize=True,
            legend_map=leg_pho,
            subtitle_lines=None,
            norm_left_margin=0.16,
            norm_y_title_offset=1.45,
        )

def main():
    parser = argparse.ArgumentParser(description="Plot PHID event efficiency vs mA per year.")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022preEE). Use 'all' to plot all years.")
    args = parser.parse_args()

    out_dir = Path(outDir)

    # CHANGED: year=all 時，除了逐年也多跑 run3
    years = YEAR_ORDER if args.year == "all" else [args.year]
    if args.year == "all":
        years = years + ["Run3"]

    for y in years:
        _plot_gen_distributions_for_year(
            base_dir=Path(baseDir),
            out_dir=Path(outDir),
            year=y,
            ma_dirs=ALL_MA_DIRS,
            ma_dirs_pt=PT_MA_DIRS,
            tree_name="inclusive",
        )

if __name__ == "__main__":
    main()
