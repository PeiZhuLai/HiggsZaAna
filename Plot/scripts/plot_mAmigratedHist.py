# ==== Imports/config for plotting ====
from pathlib import Path
import ROOT
import numpy as np
import re
from typing import Dict, List, Tuple, Optional
import uproot
import awkward as ak
import json
import argparse
import traceback
import os
import uuid

# NEW: run ROOT in batch to avoid GUI issues / silent pad problems on batch nodes
ROOT.gROOT.SetBatch(True)

# NEW: globally disable ROOT stat/fit boxes (Entries/Mean/Std Dev)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

# RENAME: clarify these are input/output roots for plotting
# FIX: eos home path should be /eos/home-<letter>/<user>/... (NOT /eos/home-<letter>-<user>/...)
input_root_dir = "/eos/home-p/pelai/HZa/root_P2Root/run3_mergedBDT"
output_plot_dir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/mAmigratedHist"

# Define the list of mA values (these are "true mass" directories)
signal_ma_dirs = ["mA_M1", "mA_M2", "mA_M3", "mA_M4", "mA_M5", "mA_M6", "mA_M7", "mA_M8", "mA_M9", "mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, '2024':108.95,
            'Run3':171.41 }

ALT_WEIGHT_CANDS = ["weight", "genWeight", "eventWeight"]

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

# NEW: y-axis label with bin width
def _ylabel_with_binwidth(hist: ROOT.TH1, base: str, *, unit: str = "GeV") -> str:
    """
    e.g.
      base='Events' -> 'Events / (0.67 GeV)'
      base='A.U.' -> 'A.U. / (0.67 GeV)'
    """
    try:
        bw = float(hist.GetXaxis().GetBinWidth(1))
        if bw > 0:
            # keep compact but readable
            bw_str = f"{bw:.3g}"
            return f"{base} / ({bw_str} {unit})"
    except Exception:
        pass
    return str(base)

def _lumi_fb_for_year(year: str) -> Optional[float]:
    v = lumiMap.get(year)
    return float(v) if v is not None else None

def _to_int(v) -> Optional[int]:
    if v is None:
        return None
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, float):
        return int(round(v))
    if isinstance(v, str):
        m = re.search(r"-?\d+", v)
        return int(m.group(0)) if m else None
    return None

def _to_float(v) -> Optional[float]:
    if v is None:
        return None
    if isinstance(v, (int, float, np.integer, np.floating)):
        return float(v)
    if isinstance(v, str):
        m = re.search(r"-?\d+(?:\.\d+)?", v)
        return float(m.group(0)) if m else None
    return None

def _parse_ma(name: str) -> Optional[int]:
    m = re.search(r"[Mm](\d+)", str(name))
    return int(m.group(1)) if m else None

def _parse_mva_cuts(path: str) -> Dict[int, float]:
    def _entry_to_pair(entry: dict) -> Tuple[Optional[int], Optional[float]]:
        mass_keys = ("mA", "ma", "mass", "massA")
        cut_keys  = ("MVAcut", "mvaCut", "mva_cut", "cut", "bdtCut", "bdt_cut", "best_MVAcut")
        ma = None
        for k in mass_keys:
            if k in entry:
                ma = _to_int(entry[k]); break
        if ma is None:
            for k in ("sample", "name", "label", "title"):
                if k in entry and isinstance(entry[k], str):
                    ma = _parse_ma(entry[k])
                    if ma is not None:
                        break
        cut = None
        for k in cut_keys:
            if k in entry:
                cut = _to_float(entry[k]); break
        if cut is None:
            for k in ("best", "opt", "result", "payload"):
                if k in entry and isinstance(entry[k], dict):
                    for ck in cut_keys:
                        if ck in entry[k]:
                            cut = _to_float(entry[k][ck]); break
                if cut is not None:
                    break
        return ma, cut

    try:
        with open(path, "r") as f:
            data = json.load(f)
    except Exception:
        return {}

    out: Dict[int, float] = {}
    if isinstance(data, list):
        for ent in data:
            if not isinstance(ent, dict): 
                continue
            ma, cut = _entry_to_pair(ent)
            if ma is not None and cut is not None:
                out[ma] = cut

    if not out and isinstance(data, dict):
        for key in ("results", "points", "entries"):
            arr = data.get(key)
            if isinstance(arr, list):
                for ent in arr:
                    if not isinstance(ent, dict):
                        continue
                    ma, cut = _entry_to_pair(ent)
                    if ma is not None and cut is not None:
                        out[ma] = cut
                if out:
                    break

    if not out and isinstance(data, dict):
        for v in data.values():
            if isinstance(v, dict):
                for key in ("results", "points", "entries"):
                    arr = v.get(key)
                    if isinstance(arr, list):
                        for ent in arr:
                            if not isinstance(ent, dict):
                                continue
                            ma, cut = _entry_to_pair(ent)
                            if ma is not None and cut is not None:
                                out[ma] = cut
                        if out:
                            break
            if out:
                break
    return out

def _pick_weight_branch(t) -> Optional[str]:
    keys = set(map(str, t.keys()))
    for k in ALT_WEIGHT_CANDS:
        if k in keys:
            return k
    return None

def _tree_pick(f) -> Optional[str]:
    # mergedBDT 可能用 test 或 inclusive；兩者都試
    for cand in ("test", "inclusive"):
        try:
            if cand in f:
                return cand
        except Exception:
            pass
    return None

# NEW: pick score prefix with controllable preference/override
def _pick_score_prefix(keys: set, prefer: str = "auto") -> Optional[str]:
    """
    prefer:
      - 'auto': prefer MVA then BDT
      - 'mva' : only MVA
      - 'bdt' : only BDT
    """
    prefer = (prefer or "auto").strip().lower()
    has_mva = any(k.startswith("MVA_Score_mA_M") for k in keys)
    has_bdt = any(k.startswith("BDT_Score_mA_M") for k in keys)

    if prefer == "mva":
        return "MVA_Score_mA_M" if has_mva else None
    if prefer == "bdt":
        return "BDT_Score_mA_M" if has_bdt else None

    return "MVA_Score_mA_M"

def _true_ma_from_dir(ma_tag: str) -> Optional[int]:
    m = re.match(r"mA_M(\d+)$", str(ma_tag))
    return int(m.group(1)) if m else None

def _accumulate_pred_pass_sumw_by_year(
    base_dir: Path,
    *,
    years: List[str],
    true_ma_dirs: List[str],
    pred_mas: List[int],
    mva_cuts: Dict[int, float],
    debug: bool = False,
    max_files: Optional[int] = None,
    score_prefer: str = "auto",
    strict: bool = False,  # NEW
) -> Dict[str, Dict[int, Dict[int, float]]]:
    """
    回傳: year -> pred_ma -> true_ma -> sumw_pass
    (sumw_pass: 在 true 檔案內，用 branch BDT_Score_mA_M{pred} >= cut(pred) 的加權通過數)
    """
    out: Dict[str, Dict[int, Dict[int, float]]] = {y: {p: {} for p in pred_mas} for y in years}
    n_opened = 0

    for true_tag in true_ma_dirs:
        true_ma = _true_ma_from_dir(true_tag)
        if true_ma is None:
            if debug:
                print(f"[DEBUG] skip true_tag (cannot parse true mA): {true_tag}")
            continue
        d_true = base_dir / true_tag
        if not d_true.is_dir():
            if debug:
                print(f"[DEBUG] missing dir: {d_true}")
            continue

        for year in years:
            fp = d_true / f"{year}.root"
            if not fp.exists():
                if debug:
                    print(f"[DEBUG] missing file: {fp}")
                continue

            if max_files is not None and n_opened >= max_files:
                if debug:
                    print(f"[DEBUG] max_files reached ({max_files}), stop reading more files.")
                break

            try:
                with uproot.open(str(fp)) as f:
                    n_opened += 1
                    tname = _tree_pick(f)
                    if not tname:
                        if debug:
                            print(f"[DEBUG] no tree found in: {fp} (tried test/inclusive)")
                        if strict:
                            raise RuntimeError(f"No tree found in {fp} (tried test/inclusive)")
                        continue
                    t = f[tname]
                    keys = set(map(str, t.keys()))
                    wname = _pick_weight_branch(t)
                    if debug:
                        print(f"[DEBUG] {fp}: tree={tname} score_prefix={_pick_score_prefix(keys, prefer=score_prefer)} weight={wname}")

                    pred_here = []
                    need = []
                    miss_cut = 0
                    miss_branch = 0

                    score_prefix = _pick_score_prefix(keys, prefer=score_prefer)
                    if score_prefix is None:
                        if debug:
                            print(f"[DEBUG] {fp}: no MVA_Score_mA_M* nor BDT_Score_mA_M* branches found.")
                        if strict:
                            raise RuntimeError(f"No score branches found in {fp} tree={tname}")
                        continue

                    for p in pred_mas:
                        if p not in mva_cuts:
                            miss_cut += 1
                            continue
                        b = f"{score_prefix}{p}"
                        if b in keys:
                            pred_here.append(p)
                            need.append(b)
                        else:
                            miss_branch += 1

                    if not pred_here:
                        if debug:
                            print(f"[DEBUG] {fp}: pred_here empty (miss_cut={miss_cut}, miss_branch={miss_branch})")
                        if strict:
                            raise RuntimeError(f"No usable pred branches in {fp} (check cuts/branches/prefix)")
                        continue

                    if wname:
                        need.append(wname)

                    any_positive = False
                    for arrs in t.iterate(need, library="ak", step_size="200 MB"):
                        w = ak.values_astype(arrs[wname], np.float64) if wname else ak.ones_like(arrs[need[0]], dtype=np.float64)

                        for p in pred_here:
                            b = f"{score_prefix}{p}"
                            cut = float(mva_cuts.get(p, 999.0))
                            sc = arrs[b]
                            mask = sc >= cut
                            sw = float(ak.sum(w[mask]))
                            if sw == 0.0:
                                continue
                            any_positive = True
                            cur = out[year][p].get(true_ma, 0.0)
                            out[year][p][true_ma] = cur + sw

                    if debug and not any_positive:
                        print(f"[DEBUG] {fp}: all pred branches have 0 weighted pass (after cuts).")

            except Exception as e:
                if debug:
                    print(f"[DEBUG] exception while reading {fp}: {e}")
                    print(traceback.format_exc())
                if strict:
                    raise
                continue

        if max_files is not None and n_opened >= max_files:
            break

    # prune empties
    out2: Dict[str, Dict[int, Dict[int, float]]] = {}
    for y, by_pred in out.items():
        by_pred2 = {p: d for p, d in by_pred.items() if d}
        if by_pred2:
            out2[y] = by_pred2
    if debug:
        n_year = len(out2)
        n_pred = sum(len(v) for v in out2.values())
        print(f"[DEBUG] pred_pass summary: years_with_data={n_year}, total_pred_with_data={n_pred}")
    return out2

def _root_color(hex_color: str, fallback: int = 1) -> int:
    """Return ROOT color index for a hex string (e.g. '#FF00AA')."""
    try:
        return int(ROOT.TColor.GetColor(str(hex_color)))
    except Exception:
        return int(fallback)

def _merge_pred_pass_years_to_run3(
    pred_pass_by_year: Dict[str, Dict[int, Dict[int, float]]],
    years: List[str],
    *,
    run3_label: str = "Run3",
) -> Dict[str, Dict[int, Dict[int, float]]]:
    """
    pred_pass_by_year: year -> pred -> true -> sumw
    回傳: {Run3: pred -> true -> sumw} (跨年相加)
    """
    merged: Dict[int, Dict[int, float]] = {}
    for y in years:
        by_pred = pred_pass_by_year.get(y) or {}
        for p, true_map in by_pred.items():
            mp = merged.setdefault(int(p), {})
            for tma, sw in (true_map or {}).items():
                mp[int(tma)] = mp.get(int(tma), 0.0) + float(sw)

    # prune empties
    merged = {p: m for p, m in merged.items() if m}
    return {run3_label: merged} if merged else {}

def _eos_path_default(base_rel: str, *, eos_user: Optional[str] = None, eos_username: Optional[str] = None) -> str:
    """
    Build EOS user path:
      /eos/home-<eos_user>/<eos_username>/<base_rel>

    Defaults:
      eos_user: env EOS_USER or 'p'
      eos_username: env EOS_USERNAME or 'pelai'
    """
    eu = (eos_user or os.environ.get("EOS_USER") or "p").strip()
    un = (eos_username or os.environ.get("EOS_USERNAME") or "pelai").strip()
    base_rel = str(base_rel).lstrip("/")
    return f"/eos/home-{eu}/{un}/{base_rel}"

def _root_style_line(h, color: int, mstyle: int, lstyle: int = 1) -> None:
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetLineWidth(3)
    h.SetLineStyle(lstyle)
    h.SetMarkerStyle(mstyle)
    h.SetMarkerSize(0.9)

def _plot_single_hist_sigstyle(
    *,
    year: str,
    xlabel: str,
    ylabel: str,
    out_path: Path,
    hist: ROOT.TH1F,
    legend_label: str,
    normalize: bool = False,
    strict: bool = False,  # NEW
) -> None:
    if not hist:
        return

    cname = f"c_{_safe_root_name(out_path.stem)}_{_safe_root_name(year)}_{uuid.uuid4().hex[:8]}"
    c = ROOT.TCanvas(cname, "", 800, 600)
    left, right, bottom, top = 0.13, 0.04, 0.14, 0.08
    if normalize:
        left = 0.16
    c.SetMargin(left, right, bottom, top)
    c.SetTickx(); c.SetTicky()
    c.cd()

    h = hist
    if normalize and h.Integral() > 0:
        h.Scale(1.0 / h.Integral())

    xmax = h.GetXaxis().GetXmax()
    xmin = h.GetXaxis().GetXmin()

    frame = ROOT.TH1F(f"frame_{out_path.stem}_{year}", "", 1, xmin, xmax)
    frame.SetDirectory(0)
    frame.SetTitle("")
    frame.GetXaxis().SetTitle(xlabel)
    # NEW: include bin width in y-axis title
    frame.GetYaxis().SetTitle(_ylabel_with_binwidth(h, ylabel, unit="GeV"))

    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetTitleSize(0.055)
    frame.GetYaxis().SetTitleSize(0.055)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)

    ytitle_off = 1.15 if not normalize else 1.45
    frame.GetYaxis().SetTitleOffset(ytitle_off)

    ymax = h.GetMaximum()
    frame.SetMaximum(1.25 * ymax if ymax > 0 else 1.0)
    frame.SetMinimum(0.0)
    frame.Draw()

    h.Draw("HIST SAME")

    header_lines = [f"{year}", r"H#rightarrowZa#rightarrowll#gamma#gamma"]
    x1, y1 = 0.98, 0.88
    # NEW: 3 columns legend; widen a bit for 14 masses
    x0, y0 = 0.65, 0.70
    leg = ROOT.TLegend(x0, y0, x1, y1)
    # leg.SetNColumns(3)  # NEW
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    for s in header_lines:
        leg.AddEntry("", s, "")
    leg.AddEntry(h, legend_label, "l")
    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(left, 0.93, "#bf{CMS} #it{Simulation}")

    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetTextSize(0.040)
        lat_lumi.DrawLatex(0.96, 0.93, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    c.Modified(); c.Update()
    if ROOT.gPad:
        ROOT.gPad.RedrawAxis()

    # NEW: SaveAs can fail silently; verify output exists
    out_str = str(out_path)
    c.SaveAs(out_str)
    if not out_path.exists():
        msg = f"[WARN] SaveAs did not create file: {out_str}"
        if strict:
            raise RuntimeError(msg)
        print(msg)

    try:
        c.Close()
    except Exception:
        pass
    del frame, c

def _plot_alp_m_pred_true_from_tree_year(
    year: str,
    h_pred: ROOT.TH1F,
    h_true: ROOT.TH1F,
    out_dir: Path,
    *,
    out_tag: str = "mA_",
    debug: bool = False,  # NEW
    strict: bool = False, # NEW
) -> None:
    if not h_pred and not h_true:
        return

    if debug:
        try:
            print(f"[DEBUG] ALP_m plotting year={year} predInt={h_pred.Integral() if h_pred else None} trueInt={h_true.Integral() if h_true else None}")
            print(f"[DEBUG] out_dir={out_dir} (exists={out_dir.exists()})")
        except Exception:
            pass

    col_pred = _root_color("#E31A1C", fallback=ROOT.kRed + 1)
    col_true = _root_color("#1F78B4", fallback=ROOT.kBlue + 1)
    _root_style_line(h_pred, col_pred, 20, 1)
    _root_style_line(h_true, col_true, 21, 1)

    out_dir.mkdir(parents=True, exist_ok=True)

    _plot_single_hist_sigstyle(
        year=year,
        xlabel="m_{#gamma#gamma} (GeV)",
        ylabel="Events",
        out_path=out_dir / f"{out_tag}_pred_{year}.pdf",
        hist=h_pred,
        legend_label="m^{Predict}_{#gamma#gamma}",
        normalize=False,
        strict=strict,  # NEW
    )
    _plot_single_hist_sigstyle(
        year=year,
        xlabel="m_{#gamma#gamma} (GeV)",
        ylabel="A.U.",
        out_path=out_dir / f"{out_tag}_pred_{year}_norm.pdf",
        hist=h_pred,
        legend_label="m^{Predict}_{#gamma#gamma}",
        normalize=True,
        strict=strict,  # NEW
    )
    _plot_single_hist_sigstyle(
        year=year,
        xlabel="m_{#gamma#gamma} (GeV)",
        ylabel="Events",
        out_path=out_dir / f"{out_tag}_true_{year}.pdf",
        hist=h_true,
        legend_label="m^{True}_{#gamma#gamma}",
        normalize=False,
        strict=strict,  # NEW
    )
    _plot_single_hist_sigstyle(
        year=year,
        xlabel="m_{#gamma#gamma} (GeV)",
        ylabel="A.U.",
        out_path=out_dir / f"{out_tag}_true_{year}_norm.pdf",
        hist=h_true,
        legend_label="m^{True}_{#gamma#gamma}",
        normalize=True,
        strict=strict,  # NEW
    )

def _scan_inputs(
    base_dir: Path,
    *,
    years: List[str],
    true_ma_dirs: List[str],
    pred_mas: List[int],
    mva_cuts: Dict[int, float],
    score_prefer: str = "auto",
    max_files: int = 10,
) -> None:
    """
    --dry-run helper:
      - 掃描 ROOT 檔案：是否存在、tree 名稱(test/inclusive)、weight branch、score prefix/branches
      - 統計缺檔/缺tree/缺cut/缺branch
    """
    base_dir = Path(base_dir)
    n_opened = 0

    cnt_missing_dir = 0
    cnt_missing_file = 0
    cnt_no_tree = 0
    cnt_no_score_prefix = 0
    cnt_open_fail = 0

    miss_cut_total = 0
    miss_branch_total = 0
    ok_pred_total = 0

    print(f"[SCAN] base_dir={base_dir}")
    print(f"[SCAN] years={years}")
    print(f"[SCAN] score_prefix_prefer={score_prefer} max_files={max_files}")
    print(f"[SCAN] pred_mas(n={len(pred_mas)}) example={pred_mas[:10]}")

    for true_tag in true_ma_dirs:
        d_true = base_dir / true_tag
        if not d_true.is_dir():
            cnt_missing_dir += 1
            continue

        for year in years:
            if n_opened >= int(max_files):
                break

            fp = d_true / f"{year}.root"
            if not fp.exists():
                cnt_missing_file += 1
                continue

            try:
                with uproot.open(str(fp)) as f:
                    n_opened += 1
                    tname = _tree_pick(f)
                    if not tname:
                        cnt_no_tree += 1
                        print(f"[SCAN] {fp}: no tree found (tried test/inclusive)")
                        continue

                    t = f[tname]
                    keys = set(map(str, t.keys()))
                    wname = _pick_weight_branch(t)
                    score_prefix = _pick_score_prefix(keys, prefer=score_prefer)
                    if score_prefix is None:
                        cnt_no_score_prefix += 1
                        print(f"[SCAN] {fp}: tree={tname} NO score prefix (no MVA_Score_mA_M*/BDT_Score_mA_M*)")
                        continue

                    miss_cut = 0
                    miss_branch = 0
                    ok_pred = 0
                    for p in pred_mas:
                        if p not in mva_cuts:
                            miss_cut += 1
                            continue
                        b = f"{score_prefix}{p}"
                        if b in keys:
                            ok_pred += 1
                        else:
                            miss_branch += 1

                    miss_cut_total += miss_cut
                    miss_branch_total += miss_branch
                    ok_pred_total += ok_pred

                    print(
                        f"[SCAN] {fp}: tree={tname} weight={wname} "
                        f"score_prefix={score_prefix} ok_pred={ok_pred} miss_branch={miss_branch} miss_cut={miss_cut}"
                    )

            except Exception as e:
                cnt_open_fail += 1
                print(f"[SCAN] {fp}: open/read failed: {e}")

        if n_opened >= int(max_files):
            break

    print("[SCAN] ---- summary ----")
    print(f"[SCAN] opened_files={n_opened}")
    print(f"[SCAN] missing_dirs={cnt_missing_dir} missing_files={cnt_missing_file}")
    print(f"[SCAN] no_tree={cnt_no_tree} no_score_prefix={cnt_no_score_prefix} open_fail={cnt_open_fail}")
    print(f"[SCAN] ok_pred_total={ok_pred_total} miss_branch_total={miss_branch_total} miss_cut_total={miss_cut_total}")

def _safe_root_name(s: str) -> str:
    """Make a string safe for ROOT object names."""
    s = re.sub(r"[^A-Za-z0-9_]+", "_", str(s))
    return s[:200] if len(s) > 200 else s

def _new_th1f(name: str, title: str, bins: int, x_min: float, x_max: float) -> ROOT.TH1F:
    h = ROOT.TH1F(_safe_root_name(name), str(title), int(bins), float(x_min), float(x_max))
    h.SetDirectory(0)
    h.Sumw2()
    return h

def _accumulate_alp_m_hists_by_year(
    base_dir: Path,
    *,
    years: List[str],
    true_ma_dirs: List[str],
    pred_mas: List[int],
    mva_cuts: Dict[int, float],
    alp_m_branch: str = "ALP_m",
    bins: int = 60,
    x_min: float = 0.0,
    x_max: float = 40.0,
    debug: bool = False,
    max_files: Optional[int] = None,
    score_prefer: str = "auto",
    strict: bool = False,
) -> Dict[str, Dict[str, ROOT.TH1F]]:
    """
    回傳: year -> {"pred": TH1F, "true": TH1F}

    pred: 混合 true mA；event 若通過任一 pred cut (OR) 則填入 ALP_m
    true: 只取 diagonal (pred=true) 的 cut (分支存在且 cut 存在才做) 填入 ALP_m
    """
    base_dir = Path(base_dir)
    out: Dict[str, Dict[str, ROOT.TH1F]] = {}
    for y in years:
        out[y] = {
            "pred": _new_th1f(f"h_alp_m_pred_{y}", "", bins, x_min, x_max),
            "true": _new_th1f(f"h_alp_m_true_{y}", "", bins, x_min, x_max),
        }

    n_opened = 0

    for true_tag in true_ma_dirs:
        true_ma = _true_ma_from_dir(true_tag)
        if true_ma is None:
            if debug:
                print(f"[DEBUG] ALP_m: skip true_tag (cannot parse true mA): {true_tag}")
            continue

        d_true = base_dir / true_tag
        if not d_true.is_dir():
            if debug:
                print(f"[DEBUG] ALP_m: missing dir: {d_true}")
            continue

        for year in years:
            fp = d_true / f"{year}.root"
            if not fp.exists():
                if debug:
                    print(f"[DEBUG] ALP_m: missing file: {fp}")
                continue

            if max_files is not None and n_opened >= max_files:
                if debug:
                    print(f"[DEBUG] ALP_m: max_files reached ({max_files}), stop reading more files.")
                break

            try:
                with uproot.open(str(fp)) as f:
                    n_opened += 1
                    tname = _tree_pick(f)
                    if not tname:
                        if debug:
                            print(f"[DEBUG] ALP_m: no tree found in: {fp} (tried test/inclusive)")
                        if strict:
                            raise RuntimeError(f"ALP_m: No tree found in {fp} (tried test/inclusive)")
                        continue

                    t = f[tname]
                    keys = set(map(str, t.keys()))

                    if alp_m_branch not in keys:
                        if debug:
                            print(f"[DEBUG] ALP_m: missing branch {alp_m_branch} in {fp} tree={tname}")
                        if strict:
                            raise RuntimeError(f"ALP_m: missing branch {alp_m_branch} in {fp} tree={tname}")
                        continue

                    wname = _pick_weight_branch(t)
                    score_prefix = _pick_score_prefix(keys, prefer=score_prefer)
                    if score_prefix is None:
                        if debug:
                            print(f"[DEBUG] ALP_m: no score prefix in {fp} tree={tname}")
                        if strict:
                            raise RuntimeError(f"ALP_m: No score branches found in {fp} tree={tname}")
                        continue

                    # pred OR-mask branches (only those that exist and have cuts)
                    pred_here: List[int] = []
                    pred_branches: List[str] = []
                    miss_cut = 0
                    miss_branch = 0
                    for p in pred_mas:
                        if p not in mva_cuts:
                            miss_cut += 1
                            continue
                        b = f"{score_prefix}{p}"
                        if b in keys:
                            pred_here.append(int(p))
                            pred_branches.append(b)
                        else:
                            miss_branch += 1

                    # diagonal true-branch if available
                    true_branch = None
                    true_cut = None
                    if true_ma in mva_cuts:
                        btrue = f"{score_prefix}{true_ma}"
                        if btrue in keys:
                            true_branch = btrue
                            true_cut = float(mva_cuts[true_ma])

                    if not pred_here and true_branch is None:
                        if debug:
                            print(f"[DEBUG] ALP_m: no usable pred branches AND no diagonal true branch in {fp} (miss_cut={miss_cut}, miss_branch={miss_branch})")
                        if strict:
                            raise RuntimeError(f"ALP_m: No usable score branches in {fp} (check cuts/branches/prefix)")
                        continue

                    if debug:
                        print(
                            f"[DEBUG] ALP_m: {fp} tree={tname} weight={wname} "
                            f"score_prefix={score_prefix} pred_here={len(pred_here)} miss_branch={miss_branch} miss_cut={miss_cut} "
                            f"diag={'yes' if true_branch else 'no'}"
                        )

                    need = [alp_m_branch] + pred_branches
                    if true_branch and true_branch not in need:
                        need.append(true_branch)
                    if wname and wname not in need:
                        need.append(wname)

                    h_pred = out[year]["pred"]
                    h_true = out[year]["true"]

                    for arrs in t.iterate(need, library="ak", step_size="200 MB"):
                        alp = arrs[alp_m_branch]
                        w = ak.values_astype(arrs[wname], np.float64) if wname else ak.ones_like(alp, dtype=np.float64)

                        # pred: OR over available predictions
                        if pred_here:
                            mask_pred = ak.zeros_like(alp, dtype=bool)
                            for p in pred_here:
                                b = f"{score_prefix}{p}"
                                cut = float(mva_cuts.get(p, 999.0))
                                mask_pred = mask_pred | (arrs[b] >= cut)

                            alp_sel = alp[mask_pred]
                            w_sel = w[mask_pred]
                            for xv, wv in zip(ak.to_numpy(alp_sel), ak.to_numpy(w_sel)):
                                h_pred.Fill(float(xv), float(wv))

                        # true: diagonal only
                        if true_branch is not None and true_cut is not None:
                            mask_true = arrs[true_branch] >= true_cut
                            alp_sel = alp[mask_true]
                            w_sel = w[mask_true]
                            for xv, wv in zip(ak.to_numpy(alp_sel), ak.to_numpy(w_sel)):
                                h_true.Fill(float(xv), float(wv))

            except Exception as e:
                if debug:
                    print(f"[DEBUG] ALP_m: exception while reading {fp}: {e}")
                    print(traceback.format_exc())
                if strict:
                    raise
                continue

        if max_files is not None and n_opened >= max_files:
            break

    # prune years with no entries at all
    out2: Dict[str, Dict[str, ROOT.TH1F]] = {}
    for y, d in out.items():
        hp, ht = d["pred"], d["true"]
        if (hp.GetEntries() > 0) or (ht.GetEntries() > 0):
            out2[y] = d
    return out2

def _merge_alp_m_hists_to_run3(
    alp_h_by_year: Dict[str, Dict[str, ROOT.TH1F]],
    years: List[str],
    *,
    run3_label: str = "Run3",
) -> Optional[Dict[str, ROOT.TH1F]]:
    """
    回傳: {"pred": TH1F, "true": TH1F}，把指定 years 的 hist 相加。
    若沒有任何輸入則回傳 None。
    """
    h_pred_sum = None
    h_true_sum = None

    for y in years:
        d = alp_h_by_year.get(y)
        if not d:
            continue

        hp = d.get("pred")
        ht = d.get("true")
        if hp:
            if h_pred_sum is None:
                h_pred_sum = hp.Clone(_safe_root_name(f"h_alp_m_pred_{run3_label}"))
                h_pred_sum.SetDirectory(0)
            else:
                h_pred_sum.Add(hp)
        if ht:
            if h_true_sum is None:
                h_true_sum = ht.Clone(_safe_root_name(f"h_alp_m_true_{run3_label}"))
                h_true_sum.SetDirectory(0)
            else:
                h_true_sum.Add(ht)

    if h_pred_sum is None and h_true_sum is None:
        return None
    return {"pred": h_pred_sum, "true": h_true_sum}

def _palette_color(i: int, fallback: int = 1) -> int:
    try:
        if not palette:
            return int(fallback)
        return int(palette[int(i) % len(palette)])
    except Exception:
        return int(fallback)

def _plot_multi_hists_sigstyle(
    *,
    year: str,
    xlabel: str,
    ylabel: str,
    out_path: Path,
    hists: List[Tuple[int, ROOT.TH1F]],  # [(ma, hist), ...]
    normalize: bool = False,
    strict: bool = False,
) -> None:
    hists = [(ma, h) for (ma, h) in (hists or []) if h]
    if not hists:
        return

    cname = f"c_multi_{_safe_root_name(out_path.stem)}_{_safe_root_name(year)}_{uuid.uuid4().hex[:8]}"
    c = ROOT.TCanvas(cname, "", 900, 700)
    left, right, bottom, top = 0.14, 0.04, 0.14, 0.08
    if normalize:
        left = 0.16
    c.SetMargin(left, right, bottom, top)
    c.SetTickx(); c.SetTicky()
    c.cd()

    # normalize (clone to avoid touching originals)
    draw_hists: List[Tuple[int, ROOT.TH1F]] = []
    for ma, h in hists:
        hh = h
        if normalize:
            hh = h.Clone(_safe_root_name(f"{h.GetName()}_{year}_norm_{uuid.uuid4().hex[:6]}"))
            hh.SetDirectory(0)
            if hh.Integral() > 0:
                hh.Scale(1.0 / hh.Integral())
        draw_hists.append((ma, hh))

    # frame range from first hist; y max over all
    xax = draw_hists[0][1].GetXaxis()
    xmin, xmax = xax.GetXmin(), xax.GetXmax()
    ymax = 0.0
    for _, h in draw_hists:
        ymax = max(ymax, float(h.GetMaximum()))

    frame = ROOT.TH1F(f"frame_{out_path.stem}_{year}_{uuid.uuid4().hex[:6]}", "", 1, float(xmin), float(xmax))
    frame.SetDirectory(0)
    frame.SetTitle("")
    frame.GetXaxis().SetTitle(xlabel)
    # NEW: include bin width in y-axis title (use first drawn hist as reference)
    frame.GetYaxis().SetTitle(_ylabel_with_binwidth(draw_hists[0][1], ylabel, unit="GeV"))
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetTitleSize(0.055)
    frame.GetYaxis().SetTitleSize(0.055)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.15 if not normalize else 1.45)
    frame.SetMaximum(1.30 * ymax if ymax > 0 else 1.0)
    frame.SetMinimum(0.0)
    frame.Draw()

    # legend
    header_lines = [f"{year}", r"H#rightarrowZa#rightarrowll#gamma#gamma"]
    n_lines = len(header_lines)
    x1, y1 = 0.93, 0.88
    # NEW: 3 columns legend; widen a bit for 14 masses
    x0, y0 = 0.35, 0.55
    leg = ROOT.TLegend(x0, y0, x1, y1)
    leg.SetNColumns(3)  # NEW
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.040)
    for s in header_lines:
        leg.AddEntry("", s, "")
    for ma, h in draw_hists:
        leg.AddEntry(h, f"mA={ma}", "l")

    # draw all
    for ma, h in draw_hists:
        h.Draw("HIST SAME")
    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(left, 0.93, "#bf{CMS} #it{Simulation}")

    lumi_fb = _lumi_fb_for_year(year if year != "Run3" else "Run3")
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetTextSize(0.040)
        lat_lumi.DrawLatex(0.96, 0.93, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    c.Modified(); c.Update()
    if ROOT.gPad:
        ROOT.gPad.RedrawAxis()

    out_str = str(out_path)
    c.SaveAs(out_str)
    if not out_path.exists():
        msg = f"[WARN] SaveAs did not create file: {out_str}"
        if strict:
            raise RuntimeError(msg)
        print(msg)

    try:
        c.Close()
    except Exception:
        pass
    del frame, c

def _accumulate_alp_m_hists_by_true_ma_year(
    base_dir: Path,
    *,
    year: str,
    true_ma_dirs: List[str],
    mva_cuts: Dict[int, float],
    alp_m_branch: str = "ALP_m",
    bins: int = 60,
    x_min: float = 0.0,
    x_max: float = 40.0,
    debug: bool = False,
    max_files: Optional[int] = None,
    score_prefer: str = "auto",
    strict: bool = False,
) -> Dict[int, ROOT.TH1F]:
    """
    回傳: true_ma -> TH1F
    每個質量點(=true_ma)各自一個 histogram；只用 diagonal cut：score(true_ma) >= cut(true_ma)
    """
    base_dir = Path(base_dir)
    out: Dict[int, ROOT.TH1F] = {}
    n_opened = 0

    for true_tag in true_ma_dirs:
        true_ma = _true_ma_from_dir(true_tag)
        if true_ma is None:
            continue
        fp = base_dir / true_tag / f"{year}.root"
        if not fp.exists():
            if debug:
                print(f"[DEBUG] multi(ALP_m): missing file: {fp}")
            continue

        if max_files is not None and n_opened >= max_files:
            break

        if true_ma not in mva_cuts:
            if debug:
                print(f"[DEBUG] multi(ALP_m): missing cut for mA={true_ma} (skip)")
            continue

        try:
            with uproot.open(str(fp)) as f:
                n_opened += 1
                tname = _tree_pick(f)
                if not tname:
                    if strict:
                        raise RuntimeError(f"multi(ALP_m): No tree found in {fp} (tried test/inclusive)")
                    continue
                t = f[tname]
                keys = set(map(str, t.keys()))
                if alp_m_branch not in keys:
                    if strict:
                        raise RuntimeError(f"multi(ALP_m): missing branch {alp_m_branch} in {fp} tree={tname}")
                    continue

                score_prefix = _pick_score_prefix(keys, prefer=score_prefer)
                if score_prefix is None:
                    if strict:
                        raise RuntimeError(f"multi(ALP_m): No score prefix in {fp} tree={tname}")
                    continue

                score_branch = f"{score_prefix}{true_ma}"
                if score_branch not in keys:
                    if debug:
                        print(f"[DEBUG] multi(ALP_m): missing diag score branch {score_branch} in {fp}")
                    if strict:
                        raise RuntimeError(f"multi(ALP_m): missing diag score branch {score_branch} in {fp}")
                    continue

                wname = _pick_weight_branch(t)
                need = [alp_m_branch, score_branch]
                if wname:
                    need.append(wname)

                h = out.get(true_ma)
                if h is None:
                    h = _new_th1f(f"h_alp_m_trueMa{true_ma}_{year}", "", bins, x_min, x_max)
                    # style here (color from palette by index in true_ma_dirs ordering is set later by caller if needed)
                    out[int(true_ma)] = h

                cut = float(mva_cuts[true_ma])

                for arrs in t.iterate(need, library="ak", step_size="200 MB"):
                    alp = arrs[alp_m_branch]
                    w = ak.values_astype(arrs[wname], np.float64) if wname else ak.ones_like(alp, dtype=np.float64)
                    mask = arrs[score_branch] >= cut
                    alp_sel = alp[mask]
                    w_sel = w[mask]
                    for xv, wv in zip(ak.to_numpy(alp_sel), ak.to_numpy(w_sel)):
                        h.Fill(float(xv), float(wv))

        except Exception as e:
            if debug:
                print(f"[DEBUG] multi(ALP_m): exception while reading {fp}: {e}")
                print(traceback.format_exc())
            if strict:
                raise
            continue

    # prune truly empty
    out = {ma: h for ma, h in out.items() if h and (h.GetEntries() > 0)}
    return out

def _accumulate_alp_m_hists_predtrue_by_true_ma_year(
    base_dir: Path,
    *,
    year: str,
    true_ma_dirs: List[str],
    pred_mas: List[int],
    mva_cuts: Dict[int, float],
    alp_m_branch: str = "ALP_m",
    bins: int = 60,
    x_min: float = 0.0,
    x_max: float = 40.0,
    debug: bool = False,
    max_files: Optional[int] = None,
    score_prefer: str = "auto",
    strict: bool = False,
) -> Dict[int, Dict[str, ROOT.TH1F]]:
    """
    回傳: true_ma -> {"pred": TH1F, "true": TH1F}

    pred: 固定在同一個 true_ma 檔案內，事件通過「任一 pred cut (OR)」後填 ALP_m
    true: 固定在同一個 true_ma 檔案內，只做 diagonal cut(score(true_ma) >= cut(true_ma)) 後填 ALP_m
    """
    base_dir = Path(base_dir)
    out: Dict[int, Dict[str, ROOT.TH1F]] = {}
    n_opened = 0

    for true_tag in true_ma_dirs:
        true_ma = _true_ma_from_dir(true_tag)
        if true_ma is None:
            continue

        fp = base_dir / true_tag / f"{year}.root"
        if not fp.exists():
            if debug:
                print(f"[DEBUG] predtrue14(ALP_m): missing file: {fp}")
            continue

        if max_files is not None and n_opened >= max_files:
            break

        try:
            with uproot.open(str(fp)) as f:
                n_opened += 1
                tname = _tree_pick(f)
                if not tname:
                    if strict:
                        raise RuntimeError(f"predtrue14(ALP_m): No tree found in {fp} (tried test/inclusive)")
                    continue

                t = f[tname]
                keys = set(map(str, t.keys()))
                if alp_m_branch not in keys:
                    if strict:
                        raise RuntimeError(f"predtrue14(ALP_m): missing branch {alp_m_branch} in {fp} tree={tname}")
                    continue

                score_prefix = _pick_score_prefix(keys, prefer=score_prefer)
                if score_prefix is None:
                    if strict:
                        raise RuntimeError(f"predtrue14(ALP_m): No score prefix in {fp} tree={tname}")
                    continue

                wname = _pick_weight_branch(t)

                # pred OR-mask branches (must have cut AND branch exists)
                pred_here: List[int] = []
                pred_branches: List[str] = []
                miss_cut = 0
                miss_branch = 0
                for p in pred_mas:
                    if p not in mva_cuts:
                        miss_cut += 1
                        continue
                    b = f"{score_prefix}{p}"
                    if b in keys:
                        pred_here.append(int(p))
                        pred_branches.append(b)
                    else:
                        miss_branch += 1

                # diagonal true branch (optional)
                true_branch = None
                true_cut = None
                if true_ma in mva_cuts:
                    btrue = f"{score_prefix}{true_ma}"
                    if btrue in keys:
                        true_branch = btrue
                        true_cut = float(mva_cuts[true_ma])

                if not pred_here and true_branch is None:
                    if debug:
                        print(f"[DEBUG] predtrue14(ALP_m): no usable pred AND no diag in {fp} (miss_cut={miss_cut}, miss_branch={miss_branch})")
                    if strict:
                        raise RuntimeError(f"predtrue14(ALP_m): No usable score branches in {fp}")
                    continue

                if debug:
                    print(
                        f"[DEBUG] predtrue14(ALP_m): {fp} tree={tname} weight={wname} score_prefix={score_prefix} "
                        f"pred_here={len(pred_here)} miss_branch={miss_branch} miss_cut={miss_cut} diag={'yes' if true_branch else 'no'}"
                    )

                d = out.get(int(true_ma))
                if d is None:
                    d = {
                        "pred": _new_th1f(f"h_alp_m_pred_trueMa{true_ma}_{year}", "", bins, x_min, x_max),
                        "true": _new_th1f(f"h_alp_m_true_trueMa{true_ma}_{year}", "", bins, x_min, x_max),
                    }
                    out[int(true_ma)] = d

                need = [alp_m_branch] + pred_branches
                if true_branch and true_branch not in need:
                    need.append(true_branch)
                if wname and wname not in need:
                    need.append(wname)

                h_pred = d["pred"]
                h_true = d["true"]

                for arrs in t.iterate(need, library="ak", step_size="200 MB"):
                    alp = arrs[alp_m_branch]
                    w = ak.values_astype(arrs[wname], np.float64) if wname else ak.ones_like(alp, dtype=np.float64)

                    # pred: OR over all available pred branches
                    if pred_here:
                        mask_pred = ak.zeros_like(alp, dtype=bool)
                        for p in pred_here:
                            b = f"{score_prefix}{p}"
                            cut = float(mva_cuts.get(p, 999.0))
                            mask_pred = mask_pred | (arrs[b] >= cut)
                        alp_sel = alp[mask_pred]
                        w_sel = w[mask_pred]
                        for xv, wv in zip(ak.to_numpy(alp_sel), ak.to_numpy(w_sel)):
                            h_pred.Fill(float(xv), float(wv))

                    # true: diagonal only
                    if true_branch is not None and true_cut is not None:
                        mask_true = arrs[true_branch] >= true_cut
                        alp_sel = alp[mask_true]
                        w_sel = w[mask_true]
                        for xv, wv in zip(ak.to_numpy(alp_sel), ak.to_numpy(w_sel)):
                            h_true.Fill(float(xv), float(wv))

        except Exception as e:
            if debug:
                print(f"[DEBUG] predtrue14(ALP_m): exception while reading {fp}: {e}")
                print(traceback.format_exc())
            if strict:
                raise
            continue

    # prune empties (both pred/true empty)
    out2: Dict[int, Dict[str, ROOT.TH1F]] = {}
    for ma, d in out.items():
        hp, ht = d.get("pred"), d.get("true")
        if (hp and hp.GetEntries() > 0) or (ht and ht.GetEntries() > 0):
            out2[int(ma)] = d
    return out2

def main():
    parser = argparse.ArgumentParser(description="Plot mA migration (4-cat) and ΔR fraction bars per year.")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022preEE). Use 'all' to plot all years.")
    parser.add_argument("--pred-mas", default="same", help="Predicted mA list: 'same' or comma list like '1,2,3,...,30'.")
    parser.add_argument("--debug", action="store_true", help="Print diagnostics (why files/trees/branches are skipped).")
    parser.add_argument("--max-files", type=int, default=None, help="Debug helper: only read first N ROOT files.")
    parser.add_argument(
        "--score-prefix",
        default="auto",
        choices=["auto", "mva", "bdt"],
        help="Which score branch prefix to use for migration cuts: auto(prefer MVA), mva, bdt.",
    )
    parser.add_argument("--strict", action="store_true", help="Fail fast: raise on missing tree/branch or SaveAs failure.")  # NEW
    parser.add_argument("--dry-run", action="store_true", help="Only scan ROOT contents (tree/branches), do not plot.")      # NEW
    parser.add_argument("--scan-files", type=int, default=10, help="--dry-run: max files to scan.")                         # NEW
    parser.add_argument("--eos-user", default=None, help="EOS area letter for /eos/home-<letter>/... (default: $EOS_USER or 'p').")
    parser.add_argument("--eos-username", default=None, help="EOS username for /eos/home-*/<username>/... (default: $EOS_USERNAME or 'pelai').")
    parser.add_argument("--alp-m-bins", type=int, default=330, help="ALP_m histogram bins.")
    parser.add_argument("--alp-m-xmin", type=float, default=0.0, help="ALP_m histogram xmin.")
    parser.add_argument("--alp-m-xmax", type=float, default=33.0, help="ALP_m histogram xmax.")

    args = parser.parse_args()

    global input_root_dir
    if input_root_dir.startswith("/eos/home-") and "/HZa/root_P2Root/run3_mergedBDT" in input_root_dir:
        input_root_dir = _eos_path_default(
            "HZa/root_P2Root/run3_mergedBDT",
            eos_user=args.eos_user,
            eos_username=args.eos_username,
        )

    in_dir = Path(input_root_dir)
    out_dir = Path(output_plot_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    if args.debug:
        print(f"[DEBUG] ROOT batch={ROOT.gROOT.IsBatch()} version={ROOT.gROOT.GetVersion()}")
        print(f"[DEBUG] output_plot_dir ensured: {out_dir} (exists={out_dir.exists()})")
        print(f"[DEBUG] input_root_dir: {in_dir} (exists={in_dir.exists()})")

    mva_cuts = _parse_mva_cuts(optimized_BDT_Cut)
    if not mva_cuts:
        print(f"[錯誤] 從 JSON 讀不到任何 MVA cut：{optimized_BDT_Cut}")
        return
    if args.debug:
        print(f"[DEBUG] loaded mva_cuts: n={len(mva_cuts)} example={list(sorted(mva_cuts.items()))[:5]}")

    true_mas_int = [v for v in (_true_ma_from_dir(t) for t in signal_ma_dirs) if v is not None]
    if args.pred_mas.strip().lower() == "same":
        pred_mas = sorted(set(true_mas_int))
    else:
        try:
            pred_mas = sorted(set(int(x.strip()) for x in args.pred_mas.split(",") if x.strip()))
        except Exception:
            pred_mas = sorted(set(true_mas_int))

    years = YEAR_ORDER if args.year == "all" else [args.year]

    if args.dry_run:
        _scan_inputs(
            Path(input_root_dir),
            years=years,
            true_ma_dirs=signal_ma_dirs,
            pred_mas=pred_mas,
            mva_cuts=mva_cuts,
            score_prefer=args.score_prefix,
            max_files=args.scan_files,
        )
        return

    pred_pass = _accumulate_pred_pass_sumw_by_year(
        Path(input_root_dir),
        years=years,
        true_ma_dirs=signal_ma_dirs,
        pred_mas=pred_mas,
        mva_cuts=mva_cuts,
        debug=args.debug,
        max_files=args.max_files,
        score_prefer=args.score_prefix,
        strict=args.strict,  # NEW
    )

    if args.debug and not pred_pass:
        print("[DEBUG] pred_pass is empty: likely原因=找不到tree/branch、mva cut缺、或全都沒通過cut(加權為0)。")

    pred_run3 = _merge_pred_pass_years_to_run3(pred_pass, years, run3_label="Run3")
    by_pred_run3 = pred_run3.get("Run3") or {}

    # always make ALP_m
    if True:
        alp_h_by_year = _accumulate_alp_m_hists_by_year(
            Path(input_root_dir),
            years=years,
            true_ma_dirs=signal_ma_dirs,
            pred_mas=pred_mas,
            mva_cuts=mva_cuts,
            alp_m_branch="ALP_m",
            bins=args.alp_m_bins,
            x_min=args.alp_m_xmin,
            x_max=args.alp_m_xmax,
            debug=args.debug,
            max_files=args.max_files,
            score_prefer=args.score_prefix,
            strict=args.strict,  # NEW
        )

        if args.debug:
            print(f"[DEBUG] alp_h_by_year years={list(alp_h_by_year.keys())}")

        for y in years:
            d = alp_h_by_year.get(y)
            if not d:
                if args.debug:
                    print(f"[DEBUG] skip ALP_m plot for {y}: no hist in alp_h_by_year")
                continue
            if args.debug:
                print(f"[DEBUG] {y}: h_pred.Integral={d['pred'].Integral()} h_true.Integral={d['true'].Integral()}")

            _plot_alp_m_pred_true_from_tree_year(
                y,
                d["pred"],
                d["true"],
                out_dir,
                out_tag="mA",
                debug=args.debug,
                strict=args.strict,  # NEW
            )

        run3 = _merge_alp_m_hists_to_run3(alp_h_by_year, years, run3_label="Run3")
        if run3:
            if args.debug:
                print(f"[DEBUG] Run3: h_pred.Integral={run3['pred'].Integral()} h_true.Integral={run3['true'].Integral()}")

            _plot_alp_m_pred_true_from_tree_year(
                "Run3",
                run3["pred"],
                run3["true"],
                out_dir,
                out_tag="mA",
                debug=args.debug,
                strict=args.strict,  # NEW
            )

        # NEW: pred/true 也各自疊 14 個 hists（每個 true_ma 一條）
        ordered_true_mas = [v for v in (_true_ma_from_dir(t) for t in signal_ma_dirs) if v is not None]
        ma_to_idx = {ma: i for i, ma in enumerate(ordered_true_mas)}

        for y in years:
            by_ma = _accumulate_alp_m_hists_predtrue_by_true_ma_year(
                Path(input_root_dir),
                year=y,
                true_ma_dirs=signal_ma_dirs,
                pred_mas=pred_mas,
                mva_cuts=mva_cuts,
                alp_m_branch="ALP_m",
                bins=args.alp_m_bins,
                x_min=args.alp_m_xmin,
                x_max=args.alp_m_xmax,
                debug=args.debug,
                max_files=args.max_files,
                score_prefer=args.score_prefix,
                strict=args.strict,
            )
            if not by_ma:
                if args.debug:
                    print(f"[DEBUG] predtrue14(ALP_m): no hists for year={y}")
                continue

            pred_list: List[Tuple[int, ROOT.TH1F]] = []
            true_list: List[Tuple[int, ROOT.TH1F]] = []

            for ma in sorted(by_ma.keys(), key=lambda m: ma_to_idx.get(m, 10**9)):
                idx = ma_to_idx.get(ma, 0)
                col = _palette_color(idx, fallback=ROOT.kBlack)
                ls = line_syles[idx % len(line_syles)] if line_syles else 1

                hp = by_ma[ma].get("pred")
                if hp:
                    _root_style_line(hp, col, 20, ls)
                    pred_list.append((int(ma), hp))

                ht = by_ma[ma].get("true")
                if ht:
                    _root_style_line(ht, col, 20, ls)
                    true_list.append((int(ma), ht))

            _plot_multi_hists_sigstyle(
                year=y,
                xlabel="m_{#gamma#gamma} (GeV)",
                ylabel="Events",
                out_path=out_dir / f"mA_pred_14masses_{y}.pdf",
                hists=pred_list,
                normalize=False,
                strict=args.strict,
            )
            _plot_multi_hists_sigstyle(
                year=y,
                xlabel="m_{#gamma#gamma} (GeV)",
                ylabel="A.U.",
                out_path=out_dir / f"mA_pred_14masses_{y}_norm.pdf",
                hists=pred_list,
                normalize=True,
                strict=args.strict,
            )
            _plot_multi_hists_sigstyle(
                year=y,
                xlabel="m_{#gamma#gamma} (GeV)",
                ylabel="Events",
                out_path=out_dir / f"mA_true_14masses_{y}.pdf",
                hists=true_list,
                normalize=False,
                strict=args.strict,
            )
            _plot_multi_hists_sigstyle(
                year=y,
                xlabel="m_{#gamma#gamma} (GeV)",
                ylabel="A.U.",
                out_path=out_dir / f"mA_true_14masses_{y}_norm.pdf",
                hists=true_list,
                normalize=True,
                strict=args.strict,
            )

        # NEW: Run3 合併（同一個 true_ma 的 pred/true 各自跨年相加）
        if years:
            run3_sum_pred: Dict[int, ROOT.TH1F] = {}
            run3_sum_true: Dict[int, ROOT.TH1F] = {}

            for y in years:
                by_ma = _accumulate_alp_m_hists_predtrue_by_true_ma_year(
                    Path(input_root_dir),
                    year=y,
                    true_ma_dirs=signal_ma_dirs,
                    pred_mas=pred_mas,
                    mva_cuts=mva_cuts,
                    alp_m_branch="ALP_m",
                    bins=args.alp_m_bins,
                    x_min=args.alp_m_xmin,
                    x_max=args.alp_m_xmax,
                    debug=args.debug,
                    max_files=args.max_files,
                    score_prefer=args.score_prefix,
                    strict=args.strict,
                )
                for ma, d in (by_ma or {}).items():
                    hp, ht = d.get("pred"), d.get("true")

                    if hp:
                        if ma not in run3_sum_pred:
                            hh = hp.Clone(_safe_root_name(f"h_alp_m_pred_trueMa{ma}_Run3"))
                            hh.SetDirectory(0)
                            run3_sum_pred[int(ma)] = hh
                        else:
                            run3_sum_pred[int(ma)].Add(hp)

                    if ht:
                        if ma not in run3_sum_true:
                            hh = ht.Clone(_safe_root_name(f"h_alp_m_true_trueMa{ma}_Run3"))
                            hh.SetDirectory(0)
                            run3_sum_true[int(ma)] = hh
                        else:
                            run3_sum_true[int(ma)].Add(ht)

            pred_list: List[Tuple[int, ROOT.TH1F]] = []
            true_list: List[Tuple[int, ROOT.TH1F]] = []
            for ma in sorted(set(run3_sum_pred.keys()) | set(run3_sum_true.keys()), key=lambda m: ma_to_idx.get(m, 10**9)):
                idx = ma_to_idx.get(ma, 0)
                col = _palette_color(idx, fallback=ROOT.kBlack)
                ls = line_syles[idx % len(line_syles)] if line_syles else 1

                if ma in run3_sum_pred:
                    _root_style_line(run3_sum_pred[ma], col, 20, ls)
                    pred_list.append((int(ma), run3_sum_pred[ma]))
                if ma in run3_sum_true:
                    _root_style_line(run3_sum_true[ma], col, 20, ls)
                    true_list.append((int(ma), run3_sum_true[ma]))

            _plot_multi_hists_sigstyle(
                year="Run3",
                xlabel="m_{#gamma#gamma} (GeV)",
                ylabel="Events",
                out_path=out_dir / "mA_pred_14masses_Run3.pdf",
                hists=pred_list,
                normalize=False,
                strict=args.strict,
            )
            _plot_multi_hists_sigstyle(
                year="Run3",
                xlabel="m_{#gamma#gamma} (GeV)",
                ylabel="A.U.",
                out_path=out_dir / "mA_pred_14masses_Run3_norm.pdf",
                hists=pred_list,
                normalize=True,
                strict=args.strict,
            )
            _plot_multi_hists_sigstyle(
                year="Run3",
                xlabel="m_{#gamma#gamma} (GeV)",
                ylabel="Events",
                out_path=out_dir / "mA_true_14masses_Run3.pdf",
                hists=true_list,
                normalize=False,
                strict=args.strict,
            )
            _plot_multi_hists_sigstyle(
                year="Run3",
                xlabel="m_{#gamma#gamma} (GeV)",
                ylabel="A.U.",
                out_path=out_dir / "mA_true_14masses_Run3_norm.pdf",
                hists=true_list,
                normalize=True,
                strict=args.strict,
            )

if __name__ == "__main__":
    main()
