# ==== Photon electron-veto (CSEV) study: significance & eff_sig/eff_bkg vs mA ====
# 比照 plot_SIP3DsigniVmA.py，研究「拿掉 photon e-veto」對 cutflow 的影響。
# 兩條 scenario：
#   zgammas_ph_eveto_w    : 保留 e-veto（baseline，== nominal has_2gamma_cand）
#   zgammas_ph_noeveto_w  : 拿掉 e-veto
# 與 SIP3D 圖不同之處：背景 b / eff_bkg 為「scenario-specific」，因為 e-veto 主要
# 作用是壓 Z->ee 電子假光子背景（DYJets），拿掉後背景會變多，significance/eff_bkg
# 才能反映此效應。s/b 加權慣例與 SIP3D 一致（signal 用 lumi-scaled _w，bkg 用非加權
# raw MC 計數的 scenario 對應鍵）。
import os
from pathlib import Path
import sys
import ROOT
import numpy as np
import re
from typing import Dict, List, Tuple, Optional
import json
import argparse
import ctypes
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

baseDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list"
outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/evetosigniVmA"

CUTS_TO_PLOT = ["event"]

# 每個 scenario(weighted signal 鍵) 對應的「背景 scenario 鍵」(非加權 raw MC 計數)。
BKG_SCENARIO_FOR = {
    "zgammas_ph_eveto_w":   "zgammas_ph_eveto",
    "zgammas_ph_noeveto_w": "zgammas_ph_noeveto",
}

PLOT_GROUPS = [
    {
        "key": "ph_eveto",
        "out_prefix": "PhEveto",
        "scenarios": ["zgammas_ph_eveto_w", "zgammas_ph_noeveto_w"],
        "legend": {
            "zgammas_ph_eveto_w":   "Applied #gamma e-veto",
            "zgammas_ph_noeveto_w": "Removed #gamma e-veto",
        },
        "legend_order": ["zgammas_ph_eveto_w", "zgammas_ph_noeveto_w"],
        "style": {
            "zgammas_ph_eveto_w":   {"color": "#1f77b4", "ls": 1, "ms": 20},
            "zgammas_ph_noeveto_w": {"color": "#d62728", "ls": 2, "ms": 21},
        },
    },
]

SCENARIOS_TO_PLOT = sorted({
    scenario
    for group in PLOT_GROUPS
    for scenario in group["scenarios"]
})

legend_map = {
    scenario: label
    for group in PLOT_GROUPS
    for scenario, label in group["legend"].items()
}

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = {
    "2016": 36.31, "2016preVFP": 19.51, "2016postVFP": 16.80,
    "2017": 41.48, "2018": 59.83,
    "2022preEE": 7.99, "2022postEE": 26.68, "2022": 34.67,
    "2023preBPix": 17.96, "2023postBPix": 9.68, "2023": 27.64,
    "2024": 109.82, "2025": 110.67, "combined_run3": 172.13,
}


def _lumi_fb_for_year(year: str) -> Optional[float]:
    v = lumiMap.get(year)
    return float(v) if v is not None else None


def _parse_year_from_filename(fn: str) -> Optional[str]:
    for y in YEAR_ORDER:
        if y in fn:
            return y
    return None


def _parse_ma_from_filename(fn: str) -> Optional[int]:
    m = re.search(r"mA_M(\d+)", fn)
    return int(m.group(1)) if m else None


def _load_cutflow_point(path: str) -> Optional[Tuple[str, int, Dict[str, Dict[str, float]]]]:
    """讀 signal cutflow JSON → (year, ma, {scenario: {all, all cuts, event}})。"""
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
    effs_by_scenario: Dict[str, Dict[str, float]] = {}
    for scenario in SCENARIOS_TO_PLOT:
        cf = cutflows.get(scenario) or {}
        allv = cf.get("all")
        if allv is None or float(allv) <= 0:
            continue
        effs: Dict[str, float] = {"all": float(allv)}
        allcuts = cf.get("all cuts")
        if allcuts is None:
            allcuts = cf.get("event")
        if allcuts is not None:
            effs["all cuts"] = float(allcuts)
        for c in CUTS_TO_PLOT:
            v = cf.get(c)
            if v is None and c == "event":
                v = cf.get("all cuts")
            if v is None:
                continue
            effs[c] = float(v)
        if effs:
            effs_by_scenario[scenario] = effs

    if not effs_by_scenario:
        return None
    return year, ma, effs_by_scenario


def _is_bkg_file(fn: str) -> bool:
    low = fn.lower()
    return any(k in low for k in ["bkg", "background", "dyg", "dyjet", "dyjets"])


def _load_bkg_scenarios_by_year(in_dir: Path, bkg_scenarios: List[str]) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    回傳 scenario -> year -> {'all':a, 'allcuts':ac, 'event':ac}，跨 bkg 檔累加。
    使用「非加權 raw MC 計數」(與 SIP3D 圖一致)；scenario 為背景 cutflow 鍵。
    """
    sums: Dict[str, Dict[str, Dict[str, float]]] = {s: {} for s in bkg_scenarios}
    for p in sorted(in_dir.glob("cutflow_*.json")):
        fn = p.name
        year = _parse_year_from_filename(fn)
        if year is None or not _is_bkg_file(fn):
            continue
        try:
            with open(str(p), "r") as f:
                data = json.load(f)
        except Exception:
            continue
        cutflows = data.get("cutflows") or {}
        for s in bkg_scenarios:
            cf = cutflows.get(s)
            if not isinstance(cf, dict):
                continue
            a = cf.get("all")
            ac = cf.get("all cuts")
            if ac is None:
                ac = cf.get("event")
            try:
                if a is None or float(a) <= 0 or ac is None:
                    continue
                d = sums[s].setdefault(year, {"all": 0.0, "allcuts": 0.0})
                d["all"] += float(a)
                d["allcuts"] += float(ac)
            except Exception:
                continue
    return sums


def _bkg_event(scn_w: str, year: str) -> float:
    bscn = BKG_SCENARIO_FOR.get(scn_w)
    if bscn is None:
        return 0.0
    d = _BKG_BY_SCENARIO_YEAR.get(bscn, {}).get(year)
    return float(d["allcuts"]) if d else 0.0


def _bkg_eff(scn_w: str, year: str) -> float:
    bscn = BKG_SCENARIO_FOR.get(scn_w)
    if bscn is None:
        return 0.0
    d = _BKG_BY_SCENARIO_YEAR.get(bscn, {}).get(year)
    if not d or d["all"] <= 0:
        return 0.0
    return float(d["allcuts"]) / float(d["all"])


def _significance_asimov_like(s: float, b: float) -> float:
    if s <= 0.0 or b <= 0.0:
        return 0.0
    val = (2.0 * (s + b) * math.log(1.0 + (s / b))) - (2.0 * s)
    if val <= 0.0:
        return 0.0
    return math.sqrt(val)


def _g_from_xy(xs: List[int], ys: List[float]) -> ROOT.TGraph:
    from array import array as carray
    x = carray("d", [float(v) for v in xs])
    y = carray("d", [float(v) for v in ys])
    return ROOT.TGraph(len(xs), x, y)


def _hex_to_rgb01(hexstr: str) -> Tuple[int, int, int]:
    s = hexstr.strip().lstrip("#")
    if len(s) != 6:
        raise ValueError(f"Invalid hex color: {hexstr}")
    return int(s[0:2], 16), int(s[2:4], 16), int(s[4:6], 16)


def _root_color(hexstr: str, *, fallback: int = 1) -> int:
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


def _positive_ys_from_graph(g: ROOT.TGraph) -> List[float]:
    ys: List[float] = []
    for i in range(int(g.GetN())):
        x = ctypes.c_double(0.0)
        y = ctypes.c_double(0.0)
        g.GetPoint(i, x, y)
        yv = float(y.value)
        if yv > 0.0:
            ys.append(yv)
    return ys


def _plot_group_year(year: str, points: Dict[int, Dict[str, Dict[str, float]]], out_dir: Path, group: Dict) -> None:
    mas = sorted(points.keys())
    if not mas:
        return
    group_key = str(group["key"])
    out_prefix = str(group["out_prefix"])
    scenarios_to_plot = list(group["scenarios"])
    legend_order = list(group.get("legend_order", scenarios_to_plot))
    group_legend = dict(group.get("legend", {}))
    group_style_cfg = dict(group.get("style", {}))

    scenario_style = {
        scenario: {
            "color": _root_color(str(st.get("color", "#000000"))),
            "ls": int(st.get("ls", 1)),
            "ms": int(st.get("ms", 20)),
        }
        for scenario, st in group_style_cfg.items()
    }

    # ---- significance vs mA (per-scenario bkg) ----
    scenario_graphs: Dict[str, ROOT.TGraph] = {}
    for scenario in scenarios_to_plot:
        b_evt = _bkg_event(scenario, year)
        xs, ys = [], []
        for ma in mas:
            effs = (points.get(ma) or {}).get(scenario) or {}
            if "event" in effs:
                xs.append(ma)
                ys.append(_significance_asimov_like(float(effs["event"]), b_evt))
        if xs:
            scenario_graphs[scenario] = _g_from_xy(xs, ys)
        print(f"[DEBUG] group={group_key} year={year} scenario={scenario} bkg_event={b_evt:.6g}")

    if not scenario_graphs:
        return

    c1 = ROOT.TCanvas(f"c_{group_key}_significance_{year}", "", 800, 600)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    c1.SetMargin(0.13, 0.04, 0.13, 0.08)
    c1.SetTickx()
    c1.SetTicky()

    _keepalive = []
    h_frame = ROOT.TH1F(f"frame_{group_key}_significance_{year}", "", 1, 0, 31)
    h_frame.SetDirectory(0)
    _keepalive.append(h_frame)

    def _auto_yrange(ref: Optional[ROOT.TGraph], lo_default, hi_default) -> Tuple[float, float]:
        if ref is None:
            return lo_default, hi_default
        ys = _positive_ys_from_graph(ref)
        if not ys:
            return lo_default, hi_default
        ymin, ymax = min(ys), max(ys)
        if ymax <= ymin:
            ymax = ymin * 2.0 if ymin > 0 else 1.0
        return max(1e-6, ymin * 0.1), ymax * 1.7

    _ref_graph = scenario_graphs.get(scenarios_to_plot[0])
    y_min_frame, y_max_frame = _auto_yrange(_ref_graph, 10.0, 200.0)

    h_frame.SetTitle("")
    h_frame.GetXaxis().SetTitle("m_{a} [GeV]")
    h_frame.GetYaxis().SetTitle("Significance (AMS)")
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.15)
    h_frame.GetXaxis().SetTitleSize(0.055)
    h_frame.GetYaxis().SetTitleSize(0.055)
    h_frame.GetXaxis().SetLabelSize(0.05)
    h_frame.GetYaxis().SetLabelSize(0.05)
    h_frame.SetMinimum(y_min_frame)
    h_frame.SetMaximum(y_max_frame)
    h_frame.Draw("AXIS")

    for scenario in scenarios_to_plot:
        g = scenario_graphs.get(scenario)
        if not g:
            continue
        st = scenario_style.get(scenario, {"color": 1, "ls": 1, "ms": 20})
        g.SetFillStyle(0)
        g.SetLineColor(int(st["color"]))
        g.SetMarkerColor(int(st["color"]))
        g.SetLineStyle(int(st["ls"]))
        g.SetLineWidth(3)
        g.SetMarkerStyle(int(st["ms"]))
        g.SetMarkerSize(1.2)
        g.Draw("LP SAME")

    leg = ROOT.TLegend(0.16, 0.72, 0.58, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.AddEntry("", f"{year}", "")
    for scenario in legend_order:
        g = scenario_graphs.get(scenario)
        if not g:
            continue
        leg.AddEntry(g, group_legend.get(scenario, legend_map.get(scenario, scenario)), "lp")
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
    out_pdf = out_dir / f"{out_prefix}_significanceVmA_{year}.pdf"
    c1.Modified()
    c1.Update()
    c1.SaveAs(str(out_pdf))
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
    _ = _keepalive

    # ---- eff_sig / eff_bkg vs mA (per-scenario bkg eff) ----
    scenario_ratio_graphs: Dict[str, ROOT.TGraph] = {}
    for scenario in scenarios_to_plot:
        eff_bkg = _bkg_eff(scenario, year)
        print(f"[DEBUG] group={group_key} year={year} scenario={scenario} bkg_eff(allcuts/all)={eff_bkg:.6g}")
        if eff_bkg <= 0.0:
            continue
        xs, ys = [], []
        for ma in mas:
            effs = (points.get(ma) or {}).get(scenario) or {}
            a = float(effs.get("all", 0.0))
            ac = float(effs.get("all cuts", 0.0))
            if a <= 0.0 or ac <= 0.0:
                continue
            eff_sig = ac / a
            ratio = eff_sig / eff_bkg
            print(
                f"[DEBUG] group={group_key} year={year} mA={ma:>2d} scenario={scenario} "
                f"sig_all={a:.6g} sig_allcuts={ac:.6g} sig_eff={eff_sig:.6g} bkg_eff={eff_bkg:.6g} "
                f"(sig_eff/bkg_eff)={ratio:.6g}"
            )
            xs.append(ma)
            ys.append(ratio)
        if xs:
            scenario_ratio_graphs[scenario] = _g_from_xy(xs, ys)

    if not scenario_ratio_graphs:
        return

    c2 = ROOT.TCanvas(f"c_effSigOverEffBkg_{group_key}_{year}", "", 800, 600)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    c2.SetMargin(0.13, 0.04, 0.13, 0.08)
    c2.SetTickx()
    c2.SetTicky()

    _keepalive2 = []
    frame2 = ROOT.TH1F(f"frame_effSigOverEffBkg_{group_key}_{year}", "", 1, 0, 31)
    frame2.SetDirectory(0)
    _keepalive2.append(frame2)

    def _auto_yrange2(ref: Optional[ROOT.TGraph]) -> Tuple[float, float]:
        lo_default, hi_default = 1.0, 300.0
        if ref is None:
            return lo_default, hi_default
        ys = _positive_ys_from_graph(ref)
        if not ys:
            return lo_default, hi_default
        ymin, ymax = min(ys), max(ys)
        if ymax <= ymin:
            ymax = ymin * 2.0 if ymin > 0 else 1.0
        return max(1e-6, ymin * 0.8), ymax * 1.2

    # 用兩條 scenario 的整體範圍決定 y 軸（scenario 間差異大）
    _all_ys = []
    for g in scenario_ratio_graphs.values():
        _all_ys += _positive_ys_from_graph(g)
    if _all_ys:
        y2_min, y2_max = max(1e-6, min(_all_ys) * 0.8), max(_all_ys) * 1.2
    else:
        y2_min, y2_max = 1.0, 300.0

    frame2.SetTitle("")
    frame2.GetXaxis().SetTitle("m_{a} [GeV]")
    frame2.GetYaxis().SetTitle("#varepsilon_{sig} / #varepsilon_{bkg}")
    frame2.GetXaxis().SetTitleOffset(1.1)
    frame2.GetYaxis().SetTitleOffset(1.25)
    frame2.GetXaxis().SetTitleSize(0.055)
    frame2.GetYaxis().SetTitleSize(0.045)
    frame2.GetXaxis().SetLabelSize(0.05)
    frame2.GetYaxis().SetLabelSize(0.05)
    frame2.SetMinimum(y2_min)
    frame2.SetMaximum(y2_max)
    frame2.Draw("AXIS")

    for scenario in scenarios_to_plot:
        g = scenario_ratio_graphs.get(scenario)
        if not g:
            continue
        st = scenario_style.get(scenario, {"color": 1, "ls": 1, "ms": 20})
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
    leg2.AddEntry("", f"{year}", "")
    for scenario in legend_order:
        g = scenario_ratio_graphs.get(scenario)
        if not g:
            continue
        leg2.AddEntry(g, group_legend.get(scenario, legend_map.get(scenario, scenario)), "lp")
    leg2.Draw()

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

    out_dir.mkdir(parents=True, exist_ok=True)
    out2 = out_dir / f"{out_prefix}_effSigOverEffBkgVmA_{year}.pdf"
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


def _plot_year(year: str, points: Dict[int, Dict[str, Dict[str, float]]], out_dir: Path) -> None:
    for group in PLOT_GROUPS:
        _plot_group_year(year, points, out_dir, group)


def main():
    parser = argparse.ArgumentParser(description="Plot photon e-veto significance / eff_sig/eff_bkg vs mA per year.")
    parser.add_argument("--year", default="all", help="Year to plot, or 'all'.")
    args = parser.parse_args()

    in_dir = Path(baseDir)
    out_dir = Path(outDir)

    # per-scenario bkg (non-weighted raw counts) by year
    global _BKG_BY_SCENARIO_YEAR
    _BKG_BY_SCENARIO_YEAR = _load_bkg_scenarios_by_year(in_dir, list(BKG_SCENARIO_FOR.values()))

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
