import os
from pathlib import Path
import re
import json
import argparse
from typing import Dict, List, Tuple, Optional

import ROOT

baseDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/cutflow/cutflow_list"

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/trigEffCompareVlepPt"

rootDir = "/eos/home-p-pelai/HZa/root_P2Root/run3"

treeName = "inclusive"

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = {
    '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
    '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45,'2024':108.95,
    'combined_run3':170.84
}

TRIG_GROUPS = {
    "ele_lead": ("trigeff_ele_lead_OR_ele", "trigeff_ele_lead_double_ele", "Electron lead"),
    "ele_sublead": ("trigeff_ele_sublead_OR_ele", "trigeff_ele_sublead_double_ele", "Electron sublead"),
    "mu_lead": ("trigeff_mu_lead_OR_mu", "trigeff_mu_lead_double_mu", "Muon lead"),
    "mu_sublead": ("trigeff_mu_sublead_OR_mu", "trigeff_mu_lead_double_mu", "Muon sublead"),
}

PT_BIN_ORDER = [
    "pt8to10","pt10to12","pt12to14","pt14to16","pt16to18","pt18to20","pt20to22","pt22to24",
    "pt24to26","pt26to28","pt28to30","pt30to32","pt32to34","pt34to36","pt36to38","pt38to40",
    "pt40to42","pt42to44","pt44to46","pt46to48","pt48to50","pt50to52","pt52to54","pt54to56",
    "pt56to58","pt58to60","pt60to62","pt62to64","pt64to66","pt66to68","pt68to70","pt70to72",
    "pt72to74","pt74to76","pt76to78","pt78to80","pt80to82","pt82to84","pt84to86","pt86to88",
    "pt88to90","pt90to92","pt92to94","pt94to96","pt96to98","pt98to100","pt100toInf",
]

# --- NEW: consistent x-range and binning for all pads/hists ---
PT_XMIN = 8.0
PT_XMAX = 102.0
PT_BIN_W = 2.0
ONE_SIGMA_CL = 0.682689492137086

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

def _ptbin_center(bin_name: str) -> float:
    if bin_name.endswith("toInf"):
        m = re.match(r"pt(\d+)toInf", bin_name)
        lo = float(m.group(1)) if m else 100.0
        # For display, place the open-ended last bin at the center of a 2 GeV bin.
        # This aligns points with our fixed histogram binning [8,102] with width=2.
        return lo + 0.5 * PT_BIN_W
    m = re.match(r"pt(\d+)to(\d+)", bin_name)
    if not m:
        return float("nan")
    lo, hi = float(m.group(1)), float(m.group(2))
    return 0.5 * (lo + hi)

def _edges_from_pt_bins() -> List[float]:
    # Force fixed 2 GeV binning in [8,102] for all histograms/pads.
    n = int(round((PT_XMAX - PT_XMIN) / PT_BIN_W))
    edges = [PT_XMIN + i * PT_BIN_W for i in range(n + 1)]
    return edges

def _make_hist(name: str, edges: List[float]) -> ROOT.TH1D:
    from array import array as carray
    a = carray("d", edges)
    h = ROOT.TH1D(name, "", len(edges) - 1, a)
    h.Sumw2()
    return h

def _graph_cumulative_eff_from_hist(h: ROOT.TH1) -> ROOT.TGraph:
    tot = h.Integral(0, h.GetNbinsX() + 1)
    from array import array as carray
    xs: List[float] = []
    ys: List[float] = []
    if tot <= 0:
        return ROOT.TGraph(0, carray("d", []), carray("d", []))
    nb = h.GetNbinsX()
    for i in range(1, nb + 1):
        # Use bin-center as "threshold" x-value for a fixed-binning histogram
        x = h.GetXaxis().GetBinCenter(i)
        num = h.Integral(i, nb + 1)
        xs.append(float(x))
        ys.append(100.0 * float(num) / float(tot))
    return ROOT.TGraph(len(xs), carray("d", xs), carray("d", ys))

def _binom_eff_errors(pass_count: float, total_count: float) -> Tuple[float, float, float]:
    if total_count <= 0.0:
        return 0.0, 0.0, 0.0

    passed = max(0.0, min(float(pass_count), float(total_count)))
    total = float(total_count)
    eff = passed / total

    passed_i = int(round(passed))
    total_i = int(round(total))
    if abs(passed - passed_i) < 1e-6 and abs(total - total_i) < 1e-6 and total_i > 0:
        low = float(ROOT.TEfficiency.ClopperPearson(total_i, passed_i, ONE_SIGMA_CL, False))
        high = float(ROOT.TEfficiency.ClopperPearson(total_i, passed_i, ONE_SIGMA_CL, True))
        return eff, max(0.0, eff - low), max(0.0, high - eff)

    stat = (eff * max(0.0, 1.0 - eff) / total) ** 0.5
    return eff, stat, stat

def _g_from_bins(eff_by_bin: Dict[str, float]) -> ROOT.TGraphAsymmErrors:
    from array import array as carray
    xs: List[float] = []
    ys: List[float] = []
    exl: List[float] = []
    exh: List[float] = []
    eyl: List[float] = []
    eyh: List[float] = []
    for b in PT_BIN_ORDER:
        if b not in eff_by_bin:
            continue
        xs.append(_ptbin_center(b))
        ys.append(float(eff_by_bin[b]) * 100.0)
        exl.append(0.0)  # keep only vertical error bars
        exh.append(0.0)
        eyl.append(0.0)
        eyh.append(0.0)
    x = carray("d", xs)
    y = carray("d", ys)
    return ROOT.TGraphAsymmErrors(
        len(xs),
        x,
        y,
        carray("d", exl),
        carray("d", exh),
        carray("d", eyl),
        carray("d", eyh),
    )

def _g_from_bin_records(bins_obj: Dict[str, Dict[str, float]]) -> ROOT.TGraphAsymmErrors:
    from array import array as carray

    xs: List[float] = []
    ys: List[float] = []
    exl: List[float] = []
    exh: List[float] = []
    eyl: List[float] = []
    eyh: List[float] = []

    for b in PT_BIN_ORDER:
        rec = bins_obj.get(b)
        if not isinstance(rec, dict):
            continue

        if "in_bin" in rec and "pass_trigger" in rec:
            eff, err_low, err_high = _binom_eff_errors(
                float(rec["pass_trigger"]),
                float(rec["in_bin"]),
            )
        elif "eff" in rec:
            eff = float(rec["eff"])
            err_low = 0.0
            err_high = 0.0
        else:
            continue

        xs.append(_ptbin_center(b))
        ys.append(100.0 * eff)
        exl.append(0.0)
        exh.append(0.0)
        eyl.append(100.0 * err_low)
        eyh.append(100.0 * err_high)

    x = carray("d", xs)
    y = carray("d", ys)
    return ROOT.TGraphAsymmErrors(
        len(xs),
        x,
        y,
        carray("d", exl),
        carray("d", exh),
        carray("d", eyl),
        carray("d", eyh),
    )

def _make_eff_graph(
    *,
    eff_by_bin: Dict[str, float],
    bins_obj: Optional[Dict[str, Dict[str, float]]] = None,
) -> ROOT.TGraphAsymmErrors:
    if bins_obj:
        g = _g_from_bin_records(bins_obj)
        if g.GetN() > 0:
            return g
    return _g_from_bins(eff_by_bin)

# --- NEW: graph of pass_trigger per pT-bin (y = counts) ---
def _g_pass_trigger_from_bins(bins_obj: Dict[str, Dict[str, float]]) -> ROOT.TGraph:
    from array import array as carray
    xs: List[float] = []
    ys: List[float] = []
    for b in PT_BIN_ORDER:
        rec = bins_obj.get(b)
        if not isinstance(rec, dict):
            continue
        if "pass_trigger" not in rec:
            continue
        xs.append(_ptbin_center(b))
        ys.append(float(rec["pass_trigger"]))
    return ROOT.TGraph(len(xs), carray("d", xs), carray("d", ys))

def _hist_from_pass_trigger_bins(
    *,
    name: str,
    edges: List[float],
    bins_obj: Dict[str, Dict[str, float]],
) -> ROOT.TH1D:
    """
    Fill TH1 by pT-bin using pass_trigger as bin content.
    JSON schema expected: bins_obj[ptBin] contains {"pass_trigger": <N>, ...}
    """
    h = _make_hist(name, edges)

    # Fill by mapping each JSON pt-bin to histogram bin via FindBin(center).
    # This is robust when the JSON bins do NOT start at 8 (e.g. sublead), while our histogram is fixed [8,102], step=2.
    for b, rec in (bins_obj or {}).items():
        if not isinstance(rec, dict) or ("pass_trigger" not in rec):
            continue

        m = re.match(r"pt(\d+)to(\d+|Inf)$", str(b))
        if not m:
            continue
        lo = float(m.group(1))
        hi_s = m.group(2)

        # pt100toInf etc: center is outside [8,102] -> ignore (also consistent with dropping last bin later)
        if hi_s == "Inf":
            continue
        hi = float(hi_s)

        center = 0.5 * (lo + hi)
        if not (PT_XMIN <= center < PT_XMAX):
            continue

        ib = int(h.GetXaxis().FindBin(center))
        if ib < 1 or ib > h.GetNbinsX():
            continue

        val = float(rec["pass_trigger"])
        h.SetBinContent(ib, val)
        h.SetBinError(ib, 0.0)

    return h

def _load_pt_hists_from_root(
    *,
    root_path: Path,
    tree_name: str,
    edges: List[float],
) -> Optional[Dict[str, ROOT.TH1D]]:
    if not root_path.exists():
        return None
    f = ROOT.TFile.Open(str(root_path), "READ")
    if not f or f.IsZombie():
        return None
    t = f.Get(tree_name)
    if not t:
        f.Close()
        return None

    h: Dict[str, ROOT.TH1D] = {}
    h["ee_lead"] = _make_hist("h_ee_lead", edges)
    h["ee_sublead"] = _make_hist("h_ee_sublead", edges)
    h["mumu_lead"] = _make_hist("h_mumu_lead", edges)
    h["mumu_sublead"] = _make_hist("h_mumu_sublead", edges)

    t.Draw("Z_lead_lepton_pt>>h_ee_lead", "z_ee==1", "goff")
    t.Draw("Z_sublead_lepton_pt>>h_ee_sublead", "z_ee==1", "goff")
    t.Draw("Z_lead_lepton_pt>>h_mumu_lead", "z_mumu==1", "goff")
    t.Draw("Z_sublead_lepton_pt>>h_mumu_sublead", "z_mumu==1", "goff")

    for hh in h.values():
        hh.SetDirectory(0)
    f.Close()
    return h

def _load_trigeff_point(path: str, syst: str) -> Optional[Tuple[str, int, Dict[str, Dict[str, float]]]]:
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

    teff_block = (data.get(f"trigeff_{syst}") or {})
    if not teff_block:
        return None

    out: Dict[str, Dict[str, float]] = {}
    for trigkey in {k for pair in TRIG_GROUPS.values() for k in pair[:2]}:
        obj = teff_block.get(trigkey) or {}
        bins = obj.get("bins") or {}
        effs: Dict[str, float] = {}
        for b, rec in bins.items():
            if not isinstance(rec, dict):
                continue
            if "eff" in rec:
                effs[b] = float(rec["eff"])
        if effs:
            out[trigkey] = effs

    return (year, ma, out) if out else None

# --- NEW: load pass_trigger bins for turn-on from cutflow JSON ---
def _load_trig_bins_for_turnon(
    path: str,
    syst: str,
) -> Optional[Tuple[str, int, Dict[str, Dict[str, Dict[str, float]]]]]:
    """
    Returns (year, ma, byTrigKeyBins) where byTrigKeyBins[trigkey] == bins dict from JSON.
    bins dict keeps per-pt-bin records (needs pass_trigger).
    """
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

    teff_block = (data.get(f"trigeff_{syst}") or {})
    if not teff_block:
        return None

    out: Dict[str, Dict[str, Dict[str, float]]] = {}
    for trigkey in {k for pair in TRIG_GROUPS.values() for k in pair[:2]}:
        obj = teff_block.get(trigkey) or {}
        bins = obj.get("bins") or {}
        # keep full per-bin record, but require pass_trigger to exist in at least one bin
        if isinstance(bins, dict) and any(isinstance(v, dict) and ("pass_trigger" in v) for v in bins.values()):
            out[trigkey] = bins

    return (year, ma, out) if out else None

def _apply_style(g: ROOT.TGraph, color: int, mstyle: int, msize: float = 1.2, lw: int = 3) -> None:
    g.SetLineColor(color)
    g.SetMarkerColor(color)
    g.SetMarkerStyle(mstyle)
    g.SetMarkerSize(msize)
    g.SetLineWidth(lw)

def _cms_label(year: str) -> None:
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Simulation}")

    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)  # align top-right
        lat_lumi.SetTextSize(0.040)
        lat_lumi.DrawLatex(0.86, 0.93, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

def _plot_or_vs_double(
    *,
    year_label: str,
    ma: int,
    title_right: str,
    eff_or: Dict[str, float],
    eff_double: Dict[str, float],
    # optional: per-pt-bin records (must contain pass_trigger) for right-axis histograms
    bins_or: Optional[Dict[str, Dict[str, float]]] = None,
    bins_double: Optional[Dict[str, Dict[str, float]]] = None,
    out_path: Path,
    # --- NEW: axis controls ---
    ymax_left: Optional[float] = None,
    ymax_right: Optional[float] = None,
    # --- NEW: per-channel x-axis title override ---
    x_title: Optional[str] = None,
) -> None:
    edges = _edges_from_pt_bins()
    xmin = PT_XMIN
    xmax = PT_XMAX

    # Build pass_trigger histograms (counts) for right axis (optional)
    hpass_or: Optional[ROOT.TH1D] = None
    hpass_db: Optional[ROOT.TH1D] = None
    max_pass = 0.0

    if bins_or:
        hpass_or = _hist_from_pass_trigger_bins(
            name=f"hpass_or_{year_label}_mA{ma}_{title_right}".replace(" ", "_"),
            edges=edges,
            bins_obj=bins_or,
        )
        hpass_or.SetDirectory(0)

    if bins_double:
        hpass_db = _hist_from_pass_trigger_bins(
            name=f"hpass_db_{year_label}_mA{ma}_{title_right}".replace(" ", "_"),
            edges=edges,
            bins_obj=bins_double,
        )
        hpass_db.SetDirectory(0)

    # --- NEW: drop last bin + normalize both histograms (right axis) ---
    def _drop_last_bin(h: Optional[ROOT.TH1D]) -> None:
        if not h:
            return
        nb = h.GetNbinsX()
        h.SetBinContent(nb, 0.0)
        h.SetBinError(nb, 0.0)

    def _normalize(h: Optional[ROOT.TH1D]) -> None:
        if not h:
            return
        integ = float(h.Integral(1, h.GetNbinsX()))
        if integ > 0.0:
            h.Scale(1.0 / integ)

    _drop_last_bin(hpass_or)
    _drop_last_bin(hpass_db)
    _normalize(hpass_or)
    _normalize(hpass_db)

    # recompute max after modifications
    if hpass_or:
        max_pass = max(max_pass, float(hpass_or.GetMaximum()))
    if hpass_db:
        max_pass = max(max_pass, float(hpass_db.GetMaximum()))

    has_pass = (max_pass > 0.0)

    c = ROOT.TCanvas(f"c_{year_label}_mA{ma}_{title_right}".replace(" ", "_"), "", 800, 600)
    c.SetTickx()
    c.SetTicky(0)

    # Two transparent-overlaid pads:
    #   pad1: efficiency curves (left y-axis)
    #   pad2: pass_trigger histograms + right y-axis (only if available)
    lmargin = 0.13
    rmargin = 0.14 if has_pass else 0.04

    pad1 = ROOT.TPad("pad1", "", 0.0, 0.0, 1.0, 1.0)
    pad1.SetMargin(lmargin, rmargin, 0.15, 0.08) # left, right, bottom, top
    pad1.SetTickx()
    pad1.SetTicky(0)
    pad1.Draw()
    pad1.cd()

    g_or = _make_eff_graph(eff_by_bin=eff_or, bins_obj=bins_or)
    g_db = _make_eff_graph(eff_by_bin=eff_double, bins_obj=bins_double)

    g_or.SetTitle("")
    g_or.GetXaxis().SetTitle(x_title if x_title else "p_{T} [GeV]")

    # IMPORTANT: For TGraph, the x-axis range defaults to the min/max of existing points.
    # Sublead may miss the first bin(s) (e.g. pt8to10), which makes the drawn axis start at 10/12 GeV.
    # Force a fixed range for ALL channels.
    g_or.GetXaxis().SetLimits(PT_XMIN, PT_XMAX)

    g_or.GetYaxis().SetTitle("Trigger Efficiency (%)")
    g_or.SetMinimum(0.0)
    g_or.SetMaximum(float(ymax_left) if ymax_left is not None else 105.0)

    # --- NEW: axis font sizes ---
    g_or.GetXaxis().SetTitleSize(0.055)
    g_or.GetYaxis().SetTitleSize(0.055)
    g_or.GetXaxis().SetLabelSize(0.05)
    g_or.GetYaxis().SetLabelSize(0.05)

    g_or.GetXaxis().SetTitleOffset(1.2)
    g_or.GetYaxis().SetTitleOffset(1.2)
    g_or.GetXaxis().SetLabelOffset(0.02)

    _apply_style(g_or, ROOT.TColor.GetColor("#e42536"), 20, 1.2, 2)
    _apply_style(g_db, ROOT.TColor.GetColor("#5790fc"), 22, 1.2, 2)

    # draw points with vertical error bars and no connecting line
    g_or.Draw("AP")          # axes + points
    g_db.Draw("P SAME")      # points only (reuse axes)

    leg = ROOT.TLegend(0.24, 0.70, 0.66, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.AddEntry("", f"m_{{a}} = {ma} GeV", "")
    # --- NEW: legend text by channel (Electron/Muon) ---
    _tr = (title_right or "").lower()
    if "electron" in _tr:
        or_label = "Single- OR Double-Electron Trigger"
        db_label = "Double-Electron Trigger"
    elif "muon" in _tr:
        or_label = "Single- OR Double-Muon Trigger"
        db_label = "Double-Muon Trigger"
    else:
        # fallback (keep generic if unexpected title_right)
        or_label = "Single- OR Double-Lepton Trigger"
        db_label = "Double-Lepton Trigger"

    leg.AddEntry(g_or, or_label, "lp")
    leg.AddEntry(g_db, db_label, "lp")

    lat2 = ROOT.TLatex()
    # lat2.SetNDC()
    # lat2.SetTextFont(42)
    # lat2.SetTextSize(0.045)
    # lat2.SetTextAlign(13)
    # lat2.DrawLatex(0.26, 0.70, f"m_{{a}} = {ma} GeV")
    # lat2.DrawLatex(0.16, 0.63, title_right)

    _cms_label(year_label)

    keep: List[object] = [pad1, g_or, g_db, leg, lat2]

    # --- right-axis overlay for pass_trigger counts ---
    if has_pass:
        c.cd()
        pad2 = ROOT.TPad("pad2", "", 0.0, 0.0, 1.0, 1.0)
        pad2.SetFillStyle(4000)
        pad2.SetFrameFillStyle(4000)
        pad2.SetMargin(0.13, rmargin, 0.15, 0.08)
        pad2.SetTicks(0, 0)
        pad2.SetTickx()
        pad2.SetTicky(0) 
        pad2.Draw()
        pad2.cd()

        ymax = 1.4 * max_pass
        if ymax_right is not None:
            ymax = float(ymax_right)

        frame2 = pad2.DrawFrame(xmin, 0.0, xmax, ymax)
        frame2.SetTitle("")
        frame2.GetXaxis().SetLabelSize(0)
        frame2.GetXaxis().SetTitle("")
        frame2.GetYaxis().SetLabelSize(0)
        frame2.GetYaxis().SetTitle("")
        frame2.GetYaxis().SetTickLength(0)
        frame2.GetXaxis().SetTickLength(0)

        keep.extend([pad2, frame2])

        if hpass_or:
            hpass_or.SetLineColor(ROOT.TColor.GetColor("#e42536"))
            hpass_or.SetLineWidth(2)
            hpass_or.SetLineStyle(1)
            hpass_or.SetFillStyle(0)
            hpass_or.Draw("HIST SAME")
            keep.append(hpass_or)

        if hpass_db:
            hpass_db.SetLineColor(ROOT.TColor.GetColor("#5790fc"))
            hpass_db.SetLineWidth(2)
            hpass_db.SetLineStyle(2)
            hpass_db.SetFillStyle(0)
            hpass_db.Draw("HIST SAME")
            keep.append(hpass_db)

        # IMPORTANT: TGaxis expects coordinates in *user space* (not NDC) unless using special options.
        # Use pad2 user coordinates so the axis is really on the right edge.
        pad2.Update()
        x_right = pad2.GetUxmax()
        y_low   = pad2.GetUymin()
        y_high  = pad2.GetUymax()
        axis = ROOT.TGaxis(
            x_right, y_low,
            x_right, y_high,
            0.0, ymax,
            510,
            "+L",
        )
        axis.SetTitle("A.U. / 2 GeV")  # was: pass_trigger (counts)
        axis.SetTitleFont(42)
        axis.SetLabelFont(42)
        axis.SetTitleSize(0.055)  # was 0.038
        axis.SetLabelSize(0.05)   # was 0.032
        axis.SetTitleOffset(1.2)
        axis.Draw()

        # redraw legend on top
        pad1.cd()
        leg.Draw()

        keep.append(axis)
    else:
        leg.Draw()

    c._keep = keep

    out_path.parent.mkdir(parents=True, exist_ok=True)
    c.Modified()
    c.Update()
    c.SaveAs(str(out_path))
    ROOT.gSystem.ProcessEvents()
    c.Close()

def _plot_turnon_overlay(
    *,
    title_right: str,
    year_to_graph: Dict[str, ROOT.TGraph],
    out_path: Path,
    label_left: str,
    # --- NEW ---
    year_to_pass: Optional[Dict[str, ROOT.TGraph]] = None,
    # --- NEW: axis controls ---
    ymax_left: Optional[float] = None,
    ymax_right: Optional[float] = None,
) -> None:
    c = ROOT.TCanvas(f"c_turnon_{label_left}_{title_right}".replace(" ", "_"), "", 950, 650)
    c.SetMargin(0.13, 0.04, 0.15, 0.08)
    c.SetTickx()
    c.SetTicky(0)

    # --- NEW: secondary axis for pass_trigger ---
    pad = ROOT.TPad("pad", "", 0.0, 0.0, 1.0, 1.0)
    pad.SetMargin(0.13, 0.10, 0.15, 0.08)  # a bit more right margin for right axis
    pad.SetTickx()
    pad.SetTicky(0)
    pad.Draw()
    pad.cd()

    year_colors = {
        "2022preEE": ROOT.kAzure + 1,
        "2022postEE": ROOT.kBlue + 2,
        "2023preBPix": ROOT.kGreen + 2,
        "2023postBPix": ROOT.kOrange + 7,
        "2024": ROOT.kRed + 1,
        "run3": ROOT.kBlack,
    }

    leg = ROOT.TLegend(0.14, 0.60, 0.55, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.034)

    first = True
    keep = [leg, pad]

    # Draw cumulative-eff graphs (left axis)
    for y in YEAR_ORDER + ["run3"]:
        if y not in year_to_graph:
            continue
        g = year_to_graph[y]
        keep.append(g)
        col = year_colors.get(y, ROOT.kBlack)
        msty = 20 if y != "run3" else 21
        _apply_style(g, col, msty, 0.9, 3)
        if y == "run3":
            g.SetLineWidth(4)

        if first:
            g.SetTitle("")
            g.GetXaxis().SetTitle("p_{T} threshold [GeV]")
            g.GetYaxis().SetTitle("Cumulative efficiency: N(p_{T}>thr)/N(all) (%)")
            g.SetMinimum(0.0)
            g.SetMaximum(float(ymax_left) if ymax_left is not None else 105.0)

            # --- NEW: axis font sizes ---
            g.GetXaxis().SetTitleSize(0.055)
            g.GetYaxis().SetTitleSize(0.055)
            g.GetXaxis().SetLabelSize(0.05)
            g.GetYaxis().SetLabelSize(0.05)

            g.Draw("ALP")
            first = False
        else:
            g.Draw("LP SAME")
        leg.AddEntry(g, f"{y} (cum eff)", "lp")

    # --- NEW: overlay pass_trigger (right axis) if provided ---
    max_pass = 0.0
    if year_to_pass:
        for y, gpass in year_to_pass.items():
            n = gpass.GetN()
            for i in range(n):
                x = ROOT.Double(0.0)
                v = ROOT.Double(0.0)
                gpass.GetPoint(i, x, v)
                if float(v) > max_pass:
                    max_pass = float(v)

        if max_pass > 0:
            y2max = 1.05 * max_pass
            if ymax_right is not None:
                y2max = float(ymax_right)

            frame2 = pad.DrawFrame(
                PT_XMIN,
                0.0,
                PT_XMAX,
                y2max,
            )
            frame2.SetTitle("")
            frame2.GetXaxis().SetLabelSize(0)  # keep bottom axis from first frame
            frame2.GetXaxis().SetTitle("")
            frame2.GetYaxis().SetLabelSize(0)  # hide left labels for this frame
            frame2.GetYaxis().SetTitle("")
            keep.append(frame2)

            # draw pass graphs on top
            for y in YEAR_ORDER + ["run3"]:
                if y not in year_to_pass:
                    continue
                gpass = year_to_pass[y]
                keep.append(gpass)
                col = year_colors.get(y, ROOT.kBlack)
                msty = 24 if y != "run3" else 25
                _apply_style(gpass, col, msty, 0.7, 2)
                gpass.SetLineStyle(3)
                gpass.Draw("LP SAME")
                leg.AddEntry(gpass, f"{y} pass_trigger", "lp")

            # right axis
            # IMPORTANT: TGaxis expects coordinates in *user space* (not NDC).
            pad.Update()
            x_right = pad.GetUxmax()
            y_low   = pad.GetUymin()
            y_high  = pad.GetUymax()
            axis = ROOT.TGaxis(
                x_right, y_low,
                x_right, y_high,
                0.0, y2max,
                510,
                "+L",
            )
            axis.SetTitle("pass_trigger (counts)")
            axis.SetTitleFont(42)
            axis.SetLabelFont(42)
            axis.SetTitleSize(0.055)  # was 0.038
            axis.SetLabelSize(0.05)   # was 0.032
            axis.SetTitleOffset(1.1)
            axis.Draw()
            keep.append(axis)

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.042)
    lat.DrawLatex(0.16, 0.52, label_left)
    lat.DrawLatex(0.16, 0.47, title_right)
    keep.append(lat)

    lat_cms = ROOT.TLatex()
    lat_cms.SetNDC()
    lat_cms.SetTextFont(42)
    lat_cms.SetTextSize(0.045)
    lat_cms.DrawLatex(0.13, 0.93, "#bf{CMS} #it{Simulation}")
    keep.append(lat_cms)

    leg.Draw()
    c._keep = keep

    out_path.parent.mkdir(parents=True, exist_ok=True)
    c.Modified()
    c.Update()
    c.SaveAs(str(out_path))
    ROOT.gSystem.ProcessEvents()
    c.Close()

def _plot_combined_years_overlay(
    *,
    ma: int,
    title_right: str,
    by_year_or: Dict[str, Dict[str, float]],
    by_year_double: Dict[str, Dict[str, float]],
    by_year_or_bins: Optional[Dict[str, Dict[str, Dict[str, float]]]] = None,
    by_year_double_bins: Optional[Dict[str, Dict[str, Dict[str, float]]]] = None,
    out_path: Path,
) -> None:
    """
    Overlay multiple years in one plot:
      - solid: OR
      - dashed: Double
    Input dict maps: year -> {ptBinName -> eff}
    """
    c = ROOT.TCanvas(f"c_combined_mA{ma}_{title_right}".replace(" ", "_"), "", 950, 650)
    c.SetMargin(0.13, 0.04, 0.15, 0.08)
    c.SetTickx()
    c.SetTicky(0)

    year_colors = {
        "2022preEE": ROOT.kAzure + 1,
        "2022postEE": ROOT.kBlue + 2,
        "2023preBPix": ROOT.kGreen + 2,
        "2023postBPix": ROOT.kOrange + 7,
        "2024": ROOT.kRed + 1,
    }

    leg = ROOT.TLegend(0.14, 0.60, 0.62, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.034)

    first = True
    keep: List[object] = [leg, c]

    # draw in canonical order
    for y in YEAR_ORDER:
        if y not in by_year_or:
            continue

        g_or = _make_eff_graph(
            eff_by_bin=by_year_or[y],
            bins_obj=(by_year_or_bins or {}).get(y),
        )
        keep.append(g_or)
        col = year_colors.get(y, ROOT.kBlack)
        _apply_style(g_or, col, 20, 0.9, 2)

        if first:
            g_or.SetTitle("")
            g_or.GetXaxis().SetTitle("p_{T} (GeV)")
            g_or.GetYaxis().SetTitle("Trigger efficiency (%)")
            g_or.SetMinimum(0.0)
            g_or.SetMaximum(105.0)
            g_or.GetXaxis().SetLimits(PT_XMIN, PT_XMAX)

            # --- NEW: axis font sizes ---
            g_or.GetXaxis().SetTitleSize(0.055)
            g_or.GetYaxis().SetTitleSize(0.055)
            g_or.GetXaxis().SetLabelSize(0.05)
            g_or.GetYaxis().SetLabelSize(0.05)

            g_or.Draw("ALP")
            first = False
        else:
            g_or.Draw("LP SAME")

        leg.AddEntry(g_or, f"{y} Sigle- OR Double-Trigger", "lp")

        if y in by_year_double:
            g_db = _make_eff_graph(
                eff_by_bin=by_year_double[y],
                bins_obj=(by_year_double_bins or {}).get(y),
            )
            keep.append(g_db)
            _apply_style(g_db, col, 22, 0.9, 2)
            g_db.SetLineStyle(2)
            g_db.Draw("LP SAME")
            leg.AddEntry(g_db, f"{y} Double Trigger", "lp")

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.042)
    lat.DrawLatex(0.16, 0.52, f"m_{{a}} = {ma} GeV")
    lat.DrawLatex(0.16, 0.47, title_right)
    keep.append(lat)

    # use a generic CMS label for combined overlay
    _cms_label("combined_run3")

    leg.Draw()
    c._keep = keep

    out_path.parent.mkdir(parents=True, exist_ok=True)
    c.Modified()
    c.Update()
    c.SaveAs(str(out_path))
    ROOT.gSystem.ProcessEvents()
    c.Close()

def main() -> None:
    parser = argparse.ArgumentParser(description="Compare trigger efficiencies (OR vs Double) per mA and year.")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022postEE). Use 'all' for all years.")
    parser.add_argument("--ma", default="all", help="mA to plot (e.g. 1, 15). Use 'all' for all mA found.")
    parser.add_argument("--syst", default="nominal", help="Trigger syst key suffix: trigeff_<syst> (default: nominal).")
    parser.add_argument("--out", default=outDir, help="Output directory.")
    parser.add_argument(
        "--mode",
        default="json",
        choices=["json", "jsonpt", "rootpt"],
        help="json: existing cutflow JSON trigeff (eff vs pT bin center); jsonpt: build pT turn-on from JSON pass_trigger bins; rootpt: legacy ROOT-tree mode.",
    )
    parser.add_argument("--root-base", default=rootDir, help="Base dir containing mA_M*/{year}.root (default: rootDir).")
    parser.add_argument("--tree", default=treeName, help="TTree name (default: inclusive).")
    # --- NEW: y-axis controls ---
    parser.add_argument("--ymax-left", type=float, default=160, help="Left y-axis maximum (e.g. 105).")
    parser.add_argument("--ymax-right", type=float, default=0.2, help="Right y-axis maximum (pass_trigger axis).")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetEndErrorSize(6)

    out_dir = Path(args.out)

    # --- NEW: jsonpt mode (no ROOT files) ---
    if args.mode == "jsonpt":
        in_dir = Path(baseDir)
        edges = _edges_from_pt_bins()

        # store_bins[year][ma][trigkey] = bins_obj (per-pt record containing pass_trigger)
        store_bins: Dict[str, Dict[int, Dict[str, Dict[str, Dict[str, float]]]]] = {}

        for p in sorted(in_dir.glob("cutflow_*.json")):
            got = _load_trig_bins_for_turnon(str(p), args.syst)
            if got is None:
                continue
            year, ma, by_trig_bins = got
            store_bins.setdefault(year, {})[ma] = by_trig_bins

        years = [args.year] if args.year != "all" else YEAR_ORDER
        mas_all = sorted({ma for y in store_bins for ma in store_bins[y].keys()})
        mas = [int(args.ma)] if args.ma != "all" else mas_all

        # Define which trigger key to use for pT-shape; use OR keys (pass_trigger should be present)
        turnon_defs = [
            ("ee", "lead", "Electron lead", TRIG_GROUPS["ele_lead"][0]),
            ("ee", "sublead", "Electron sublead", TRIG_GROUPS["ele_sublead"][0]),
            ("mumu", "lead", "Muon lead", TRIG_GROUPS["mu_lead"][0]),
            ("mumu", "sublead", "Muon sublead", TRIG_GROUPS["mu_lead"][0]),
        ]

        for ma in mas:
            for (ch, leg, label, trigkey) in turnon_defs:
                year_to_graph: Dict[str, ROOT.TGraph] = {}
                # --- NEW ---
                year_to_pass: Dict[str, ROOT.TGraph] = {}

                # also build a "run3" summed hist across selected years
                h_run3: Optional[ROOT.TH1D] = None
                # --- NEW ---
                hpass_run3: Optional[ROOT.TH1D] = None

                for y in years:
                    if y not in store_bins or ma not in store_bins[y]:
                        continue
                    by_trig = store_bins[y][ma]
                    bins_obj = by_trig.get(trigkey)
                    if not bins_obj:
                        continue

                    h = _hist_from_pass_trigger_bins(
                        name=f"h_{ch}_{leg}_{y}_mA{ma}",
                        edges=edges,
                        bins_obj=bins_obj,
                    )
                    g = _graph_cumulative_eff_from_hist(h)
                    year_to_graph[y] = g

                    # --- NEW: per-bin pass_trigger graph (counts) ---
                    year_to_pass[y] = _g_pass_trigger_from_bins(bins_obj)

                    if h_run3 is None:
                        h_run3 = h.Clone(f"h_{ch}_{leg}_run3_mA{ma}")
                        h_run3.SetDirectory(0)
                    else:
                        h_run3.Add(h)

                    # --- NEW: also sum pass_trigger counts for run3 ---
                    if hpass_run3 is None:
                        hpass_run3 = h.Clone(f"hpass_{ch}_{leg}_run3_mA{ma}")
                        hpass_run3.SetDirectory(0)
                    else:
                        hpass_run3.Add(h)

                if h_run3 is not None:
                    year_to_graph["run3"] = _graph_cumulative_eff_from_hist(h_run3)
                # --- NEW ---
                if hpass_run3 is not None:
                    # reuse the same helper but build graph from summed hist content per bin
                    # easiest: convert hist -> bins_obj-like graph by reading bin centers and contents
                    from array import array as carray
                    xs: List[float] = []
                    ys: List[float] = []
                    for i in range(1, hpass_run3.GetNbinsX() + 1):
                        xs.append(hpass_run3.GetXaxis().GetBinCenter(i))
                        ys.append(hpass_run3.GetBinContent(i))
                    year_to_pass["run3"] = ROOT.TGraph(len(xs), carray("d", xs), carray("d", ys))

                if not year_to_graph:
                    continue

                out_pdf = out_dir / "jsonpt" / f"mA_M{ma}" / f"turnon_{ch}_{leg}_mA_M{ma}.pdf"
                _plot_turnon_overlay(
                    title_right=f"m_{{a}} = {ma} GeV",
                    year_to_graph=year_to_graph,
                    out_path=out_pdf,
                    label_left=label,
                    year_to_pass=year_to_pass if year_to_pass else None,
                    ymax_left=args.ymax_left,
                    ymax_right=args.ymax_right,
                )
        return

    if args.mode == "rootpt":
        years = [args.year] if args.year != "all" else YEAR_ORDER
        mas = [int(args.ma)] if args.ma != "all" else []
        if not mas:
            rb = Path(args.root_base)
            for d in sorted(rb.glob("mA_M*")):
                m = re.match(r"mA_M(\d+)", d.name)
                if m:
                    mas.append(int(m.group(1)))

        edges = _edges_from_pt_bins()

        for ma in mas:
            ma_dir = Path(args.root_base) / f"mA_M{ma}"

            combos = [
                ("ee", "lead", "Electron lead", "ee_lead"),
                ("ee", "sublead", "Electron sublead", "ee_sublead"),
                ("mumu", "lead", "Muon lead", "mumu_lead"),
                ("mumu", "sublead", "Muon sublead", "mumu_sublead"),
            ]

            for (ch, leg, label, key) in combos:
                year_to_graph: Dict[str, ROOT.TGraph] = {}

                for y in years:
                    rp = ma_dir / f"{y}.root"
                    hists = _load_pt_hists_from_root(root_path=rp, tree_name=args.tree, edges=edges)
                    if not hists:
                        continue
                    g = _graph_cumulative_eff_from_hist(hists[key])
                    year_to_graph[y] = g

                rp_run3 = ma_dir / "run3.root"
                hists_run3 = _load_pt_hists_from_root(root_path=rp_run3, tree_name=args.tree, edges=edges)
                if hists_run3:
                    year_to_graph["run3"] = _graph_cumulative_eff_from_hist(hists_run3[key])

                if not year_to_graph:
                    continue

                out_pdf = out_dir / "rootpt" / f"mA_M{ma}" / f"turnon_{ch}_{leg}_mA_M{ma}.pdf"
                _plot_turnon_overlay(
                    title_right=f"m_{{a}} = {ma} GeV",
                    year_to_graph=year_to_graph,
                    out_path=out_pdf,
                    label_left=label,
                    ymax_left=args.ymax_left,
                    ymax_right=args.ymax_right,
                )
        return

    # --- default json mode (eff vs pT bin center) ---
    in_dir = Path(baseDir)

    store: Dict[str, Dict[int, Dict[str, Dict[str, float]]]] = {}
    # also keep per-pt-bin records for pass_trigger overlay
    store_bins: Dict[str, Dict[int, Dict[str, Dict[str, Dict[str, float]]]]] = {}

    for p in sorted(in_dir.glob("cutflow_*.json")):
        # Prefer the bins-loader (contains both eff + pass_trigger); fall back to eff-only loader.
        got_bins = _load_trig_bins_for_turnon(str(p), args.syst)
        if got_bins is not None:
            year, ma, by_trig_bins = got_bins
            store_bins.setdefault(year, {})[ma] = by_trig_bins

            # derive eff dictionary from per-bin records (so we don't need to read JSON twice)
            teffs: Dict[str, Dict[str, float]] = {}
            for trigkey, bins_obj in by_trig_bins.items():
                effs: Dict[str, float] = {}
                if isinstance(bins_obj, dict):
                    for b, rec in bins_obj.items():
                        if isinstance(rec, dict) and ("eff" in rec):
                            try:
                                effs[b] = float(rec["eff"])
                            except Exception:
                                pass
                if effs:
                    teffs[trigkey] = effs
            if teffs:
                store.setdefault(year, {})[ma] = teffs
            continue

        got = _load_trigeff_point(str(p), args.syst)
        if got is None:
            continue
        year, ma, teffs = got
        store.setdefault(year, {})[ma] = teffs

    years = [args.year] if args.year != "all" else YEAR_ORDER
    mas_all = sorted({ma for y in store for ma in store[y].keys()})
    mas = [int(args.ma)] if args.ma != "all" else mas_all

    for ma in mas:
        # per-year plots
        for year in years:
            if year not in store or ma not in store[year]:
                continue

            te = store[year][ma]
            for trig_group_key, (k_or, k_db, label) in TRIG_GROUPS.items():
                if k_or not in te or k_db not in te:
                    continue

                # --- NEW: x-axis title by channel ---
                x_title = None
                if trig_group_key == "ele_lead":
                    x_title = "Lead Electron p_{T} [GeV]"
                elif trig_group_key == "ele_sublead":
                    x_title = "Sublead Electron p_{T} [GeV]"
                elif trig_group_key == "mu_lead":
                    x_title = "Lead Muon p_{T} [GeV]"
                elif trig_group_key == "mu_sublead":
                    x_title = "Sublead Muon p_{T} [GeV]"

                out_pdf = out_dir / year / f"mA_M{ma}" / f"trigeffCompare_{label.replace(' ','_')}_mA_M{ma}_{year}.pdf"
                _plot_or_vs_double(
                    year_label=year,
                    ma=ma,
                    title_right=label,
                    eff_or=te[k_or],
                    eff_double=te[k_db],
                    bins_or=store_bins.get(year, {}).get(ma, {}).get(k_or),
                    bins_double=store_bins.get(year, {}).get(ma, {}).get(k_db),
                    out_path=out_pdf,
                    ymax_left=args.ymax_left,
                    ymax_right=args.ymax_right,
                    x_title=x_title,
                )

        # combined overlay across years (eff only)
        by_year = {y: store[y][ma] for y in YEAR_ORDER if y in store and ma in store[y]}
        if len(by_year) >= 1:
            for _, (k_or, k_db, label) in TRIG_GROUPS.items():
                by_year_or: Dict[str, Dict[str, float]] = {}
                by_year_db: Dict[str, Dict[str, float]] = {}
                by_year_or_bins: Dict[str, Dict[str, Dict[str, float]]] = {}
                by_year_db_bins: Dict[str, Dict[str, Dict[str, float]]] = {}
                for y in by_year:
                    te = by_year[y]
                    if k_or in te and k_db in te:
                        by_year_or[y] = te[k_or]
                        by_year_db[y] = te[k_db]
                        if k_or in store_bins.get(y, {}).get(ma, {}):
                            by_year_or_bins[y] = store_bins[y][ma][k_or]
                        if k_db in store_bins.get(y, {}).get(ma, {}):
                            by_year_db_bins[y] = store_bins[y][ma][k_db]
                if not by_year_or:
                    continue

                out_pdf = out_dir / "combined" / f"mA_M{ma}" / f"trigeffCompare_{label.replace(' ','_')}_mA_M{ma}_combined.pdf"
                _plot_combined_years_overlay(
                    ma=ma,
                    title_right=label,
                    by_year_or=by_year_or,
                    by_year_double=by_year_db,
                    by_year_or_bins=by_year_or_bins if by_year_or_bins else None,
                    by_year_double_bins=by_year_db_bins if by_year_db_bins else None,
                    out_path=out_pdf,
                )



if __name__ == "__main__":
    main()
