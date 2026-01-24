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
import traceback  # NEW
import os  # NEW

# NEW: globally disable ROOT stat/fit boxes (Entries/Mean/Std Dev)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

# RENAME: clarify these are input/output roots for plotting
# FIX: eos home path should be /eos/home-<letter>/<user>/... (NOT /eos/home-<letter>-<user>/...)
input_root_dir = "/eos/home-p/pelai/HZa/root_P2Root/run3_mergedBDT"
output_plot_dir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/mAmigratedMatrix"

# Define the list of mA values (these are "true mass" directories)
signal_ma_dirs = ["mA_M1", "mA_M2", "mA_M3", "mA_M4", "mA_M5", "mA_M6", "mA_M7", "mA_M8", "mA_M9", "mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]

YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, '2024':108.95,
            'Run3':171.41 }

ALT_WEIGHT_CANDS = ["weight", "genWeight", "eventWeight"]

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

    # auto
    if has_mva:
        return "MVA_Score_mA_M"
    if has_bdt:
        return "BDT_Score_mA_M"
    return None

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
    debug: bool = False,          # NEW
    max_files: Optional[int] = None,  # NEW
    score_prefer: str = "auto",   # NEW
) -> Dict[str, Dict[int, Dict[int, float]]]:
    """
    回傳: year -> pred_ma -> true_ma -> sumw_pass
    (sumw_pass: 在 true 檔案內，用 branch BDT_Score_mA_M{pred} >= cut(pred) 的加權通過數)
    """
    out: Dict[str, Dict[int, Dict[int, float]]] = {y: {p: {} for p in pred_mas} for y in years}
    n_opened = 0  # NEW

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
                        continue
                    t = f[tname]
                    keys = set(map(str, t.keys()))
                    wname = _pick_weight_branch(t)
                    if debug and not wname:
                        print(f"[DEBUG] no weight branch found in {fp}. Will use ones.")

                    pred_here = []
                    need = []
                    miss_cut = 0
                    miss_branch = 0

                    # CHANGED: decide which score prefix exists in this tree (prefer MVA, since your files are MVA_*)
                    score_prefix = _pick_score_prefix(keys, prefer=score_prefer)

                    if debug and score_prefix is None:
                        print(f"[DEBUG] {fp}: no MVA_Score_mA_M* nor BDT_Score_mA_M* branches found.")
                        # keep going; pred_here will be empty

                    for p in pred_mas:
                        if p not in mva_cuts:
                            miss_cut += 1
                            continue
                        if score_prefix is None:
                            miss_branch += 1
                            continue
                        b = f"{score_prefix}{p}"
                        if b in keys:
                            pred_here.append(p)
                            need.append(b)
                        else:
                            miss_branch += 1

                    if debug:
                        print(f"[DEBUG] {fp} tree={tname} keys={len(keys)} pred_here={len(pred_here)} "
                              f"miss_cut={miss_cut} miss_branch={miss_branch} score_prefix={score_prefix} prefer={score_prefer}")

                    if not pred_here:
                        continue

                    if wname:
                        need.append(wname)

                    any_positive = False  # NEW
                    for arrs in t.iterate(need, library="ak", step_size="200 MB"):
                        w = ak.values_astype(arrs[wname], np.float64) if wname else ak.ones_like(arrs[need[0]], dtype=np.float64)

                        for p in pred_here:
                            b = f"{score_prefix}{p}"  # NEW: must match what we put into `need`
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

def _pred4cats_from_true_map(pred_ma: int, true_map: Dict[int, float]) -> Optional[Tuple[float, float, float, float]]:
    """
    4 categories (fractions, sum to 1):
      c0:  (pred-true)=0
      cP1: (pred-true)=+1  -> true = pred-1
      cM1: (pred-true)=-1  -> true = pred+1
      cGT: |pred-true|>1
    """
    if not true_map:
        return None
    tot = float(sum(true_map.values()))
    if tot <= 0:
        return None

    c0 = float(true_map.get(pred_ma, 0.0))
    cP1 = float(true_map.get(pred_ma - 1, 0.0))  # pred-true=+1
    cM1 = float(true_map.get(pred_ma + 1, 0.0))  # pred-true=-1
    cGT = max(0.0, tot - (c0 + cP1 + cM1))
    return c0 / tot, cP1 / tot, cM1 / tot, cGT / tot

# RENAME: make name match purpose (mA migration 4-cat stacked horizontal bar)
def _plot_ma_migration_4cat_bar_year(
    year: str,
    pred_to_true: Dict[int, Dict[int, float]],
    out_dir: Path,
    *,
    out_name: Optional[str] = None,
) -> None:
    """
    每個 pred mass 一根 bar，顯示 4 類比例：
      Δ=0, Δ=+1, Δ=-1, |Δ|>1
    """
    pred_mas = sorted(pred_to_true.keys())
    rows: List[Tuple[int, float, float, float, float]] = []
    for p in pred_mas:
        fr = _pred4cats_from_true_map(p, pred_to_true.get(p) or {})
        if fr is None:
            continue
        c0, cP1, cM1, cGT = fr
        rows.append((p, c0, cP1, cM1, cGT))
    if not rows:
        return

    n = len(rows)
    c = ROOT.TCanvas(f"c_mig4_{year}", "", 800, 800)
    ROOT.gStyle.SetOptStat(0)
    leftMargin = 0.14
    c.SetMargin(leftMargin, 0.06, 0.12, 0.10) # left, right, bottom, top
    c.SetTickx(); c.SetTicky()
    c.cd()

    frame = ROOT.TH2F(f"frame_mig4_{year}", "", 1, 0.0, 100.0, n, 0.0, float(n))
    frame.SetDirectory(0)
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("Fraction (%)")
    frame.GetYaxis().SetTitle("m_{a}^{Predict} (GeV)")
    frame.GetXaxis().SetTitleOffset(1.05)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetLabelSize(0.055)

    for i, (p, _, _, _, _) in enumerate(rows, start=1):
        frame.GetYaxis().SetBinLabel(i, f"{p}")
    frame.Draw("")
    if ROOT.gPad:
        ROOT.gPad.SetGridx(1)

    # 4 colors
    col0  = _root_color("#5C7AFF", fallback=ROOT.kBlue + 1)   # Δ=0
    colP1 = _root_color("#FFD557", fallback=ROOT.kOrange + 7) # Δ=+1
    colM1 = _root_color("#FF688B", fallback=ROOT.kRed + 1)    # Δ=-1
    colGT = _root_color("#9CD1FF", fallback=ROOT.kAzure + 1)  # |Δ|>1

    leg = ROOT.TLegend(leftMargin, 0.90, 0.97, 0.96)
    leg.SetNColumns(4)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.SetTextFont(42); leg.SetTextSize(0.030)
    b0 = ROOT.TBox(0,0,0,0);  b0.SetFillColor(col0);  b0.SetLineColor(col0)
    b1 = ROOT.TBox(0,0,0,0);  b1.SetFillColor(colP1); b1.SetLineColor(colP1)
    b2 = ROOT.TBox(0,0,0,0);  b2.SetFillColor(colM1); b2.SetLineColor(colM1)
    b3 = ROOT.TBox(0,0,0,0);  b3.SetFillColor(colGT); b3.SetLineColor(colGT)
    leg.AddEntry(b0, "#Deltam = 0", "f")
    leg.AddEntry(b1, "#Deltam = -1", "f")
    leg.AddEntry(b2, "#Deltam = +1", "f")
    leg.AddEntry(b3, "|#Deltam| > 1", "f")

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    # latex.DrawLatex(c.GetLeftMargin(), 0.855, "#Deltam = m_{a}^{Predict} - m_{a}^{True}")

    # NEW: percentage labels (in user coordinates, not NDC)
    pct = ROOT.TLatex()
    pct.SetTextFont(42)
    pct.SetTextSize(0.030)

    def _draw_pct(xl: float, xr: float, yy0: float, yy1: float, val: float):
        w = xr - xl
        if w < 6.0:  # too narrow in "percent units" (0~100)
            return
        pct.SetTextAlign(22)  # centered
        pct.DrawLatex((xl + xr) / 2.0, (yy0 + yy1) / 2.0 - 0.02, f"{val:.1f}%")

    _keep = [frame, leg, b0, b1, b2, b3, pct]
    bar_h = 0.70
    for idx, (p, c0f, cP1f, cM1f, cGTf) in enumerate(rows):
        y0 = float(idx) + (1.0 - bar_h) / 2.0
        y1 = float(idx) + 1.0 - (1.0 - bar_h) / 2.0

        p0 = 100.0 * c0f
        p1v = 100.0 * cP1f
        p2v = 100.0 * cM1f
        p3v = 100.0 * cGTf

        x0 = 0.0
        x1 = x0 + p0
        x2 = x1 + p1v
        x3 = x2 + p2v
        x4 = x3 + p3v

        bx0 = ROOT.TBox(x0, y0, x1, y1); bx0.SetFillColor(col0);  bx0.SetLineColor(col0)
        bx1 = ROOT.TBox(x1, y0, x2, y1); bx1.SetFillColor(colP1); bx1.SetLineColor(colP1)
        bx2 = ROOT.TBox(x2, y0, x3, y1); bx2.SetFillColor(colM1); bx2.SetLineColor(colM1)
        bx3 = ROOT.TBox(x3, y0, x4, y1); bx3.SetFillColor(colGT); bx3.SetLineColor(colGT)
        bx0.Draw("SAME"); bx1.Draw("SAME"); bx2.Draw("SAME"); bx3.Draw("SAME")
        _keep += [bx0, bx1, bx2, bx3]

        # NEW: draw percent labels on segments
        _draw_pct(x0, x1, y0, y1, p0)
        _draw_pct(x1, x2, y0, y1, p1v)
        _draw_pct(x2, x3, y0, y1, p2v)
        _draw_pct(x3, x4, y0, y1, p3v)

    # CMS line
    lat2 = ROOT.TLatex()
    lat2.SetNDC()
    lat2.SetTextFont(42)
    lat2.SetTextSize(0.04)
    lat2.DrawLatex(c.GetLeftMargin(), 0.965, "#bf{CMS} #it{Simulation}")
    _keep.append(lat2)

    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetTextSize(0.038)
        lat_lumi.DrawLatex(0.94, 0.965, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")
        _keep.append(lat_lumi)

    leg.Draw()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = out_dir / (out_name or f"mA_migration_4cat_{year}.pdf")
    c.Modified(); c.Update()
    if ROOT.gPad:
        ROOT.gPad.RedrawAxis()
    c.SaveAs(str(out_pdf))
    try:
        c.Close()
    except Exception:
        pass
    _ = _keep

def _root_color(hex_color: str, fallback: int = 1) -> int:
    """Return ROOT color index for a hex string (e.g. '#FF00AA')."""
    try:
        return int(ROOT.TColor.GetColor(str(hex_color)))
    except Exception:
        return int(fallback)

# RENAME: name reflects it returns ΔR fractions (not "eff scan")
def _load_dr_fraction_scan_by_year_ma(
    base_dir: Path,
    *,
    years: List[str],
    ma_dirs: List[str],
    tree_name: str = "inclusive",
    dr_branch: str = "var_dR_g1g2",
    dr_min: float = 0.0,
    dr_max: float = 4.0,
    dr_step: float = 0.05,
    debug: bool = False,               # NEW
    max_files: Optional[int] = None,   # NEW
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    回傳：
      year -> ma_dir -> { "bin_edges": [...], "fractions": [...] }

    fractions 是依 dr 區間的分段比例（sum=1），用於畫 stacked bar。
    bin_edges: [dr_min, ..., dr_max]，bins = len(edges)-1
    """
    base_dir = Path(base_dir)
    edges = np.arange(dr_min, dr_max + 1e-9, dr_step, dtype=np.float64)
    if len(edges) < 2:
        return {}

    out: Dict[str, Dict[str, Dict[str, float]]] = {}
    n_opened = 0  # NEW

    for year in years:
        for ma_tag in ma_dirs:
            fp = base_dir / ma_tag / f"{year}.root"
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
                    tname = tree_name if tree_name in f else _tree_pick(f)
                    if not tname:
                        if debug:
                            print(f"[DEBUG] no tree found in: {fp} (wanted '{tree_name}', fallback test/inclusive)")
                        continue
                    t = f[tname]
                    keys = set(map(str, t.keys()))
                    if dr_branch not in keys:
                        if debug:
                            print(f"[DEBUG] missing branch '{dr_branch}' in {fp} tree={tname}.")
                        continue

                    wname = _pick_weight_branch(t)
                    need = [dr_branch] + ([wname] if wname else [])

                    counts = np.zeros(len(edges) - 1, dtype=np.float64)
                    total = 0.0

                    for arrs in t.iterate(need, library="ak", step_size="200 MB"):
                        dr = ak.values_astype(arrs[dr_branch], np.float64)
                        w = ak.values_astype(arrs[wname], np.float64) if wname else ak.ones_like(dr, dtype=np.float64)

                        # 清掉 NaN/inf（避免用 ak.isfinite：舊版 awkward 沒這個 API）
                        dr_np = ak.to_numpy(dr)
                        w_np  = ak.to_numpy(w)

                        fin = np.isfinite(dr_np) & np.isfinite(w_np)
                        if not np.any(fin):
                            continue
                        dr_np = dr_np[fin]
                        w_np  = w_np[fin]

                        # 範圍外丟棄（避免最後 overflow）
                        inrng = (dr_np >= edges[0]) & (dr_np < edges[-1])
                        if not np.any(inrng):
                            continue
                        dr_np = dr_np[inrng]
                        w_np  = w_np[inrng]

                        hist, _ = np.histogram(dr_np, bins=edges, weights=w_np)
                        counts += hist
                        total += float(np.sum(w_np))

                    if total <= 0.0:
                        if debug:
                            print(f"[DEBUG] total weight <= 0 for {fp} (branch={dr_branch}, weight={wname}).")
                        continue

                    fracs = (counts / total).tolist()
                    payload = {
                        "bin_edges": edges.tolist(),
                        "fractions": fracs,
                        "totalW": float(total),  # NEW: allow Run3 weighted merge
                    }
                    out.setdefault(year, {})[ma_tag] = payload

                    if debug:
                        print(f"[DEBUG] loaded dR fractions: year={year} ma_tag={ma_tag} totalW={total:.3g}")

            except Exception as e:
                if debug:
                    print(f"[DEBUG] exception while reading {fp}: {e}")
                    print(traceback.format_exc())
                continue

        if max_files is not None and n_opened >= max_files:
            break

    if debug:
        tot = sum(len(v) for v in out.values())
        print(f"[DEBUG] dr_frac_scan summary: years_with_data={len(out)}, total_ma_entries={tot}")
    return out

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

def _merge_dr_fracs_years_to_run3(
    dr_frac_by_year: Dict[str, Dict[str, Dict[str, float]]],
    years: List[str],
    *,
    run3_label: str = "Run3",
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    dr_frac_by_year: year -> ma_tag -> {bin_edges, fractions}
    用每個 (year, ma_tag) 的 total weight 當權重做加權平均：
      weight = totalW = sum(counts) = sum(frac)*totalW -> 我們只有 frac，因此用 totalW 需另外估。
    由於 loader 內部已用 counts/totalW 的 fractions，但沒有存 totalW，
    這裡改用「以該 ma_tag 在該 year 的分佈總量 proxy=1」是不合理。
    所以：在 loader payload 內補存 totalW（最小侵入：如果沒 totalW 就退化成等權重平均）。
    """
    # year -> ma_tag -> payload
    out: Dict[str, Dict[str, Dict[str, float]]] = {}

    # collect union ma_tags
    ma_tags = set()
    for y in years:
        ma_tags |= set((dr_frac_by_year.get(y) or {}).keys())

    by_ma: Dict[str, Dict[str, float]] = {}
    for ma_tag in sorted(ma_tags):
        # pick reference edges from first available year
        ref_edges = None
        nbins = None
        for y in years:
            payload = (dr_frac_by_year.get(y) or {}).get(ma_tag)
            if payload and isinstance(payload.get("bin_edges"), list) and isinstance(payload.get("fractions"), list):
                ref_edges = payload["bin_edges"]
                nbins = len(payload["fractions"])
                break
        if ref_edges is None or nbins is None:
            continue

        num = np.zeros(nbins, dtype=np.float64)
        den = 0.0

        for y in years:
            payload = (dr_frac_by_year.get(y) or {}).get(ma_tag)
            if not payload:
                continue
            edges = payload.get("bin_edges")
            fracs = payload.get("fractions")
            if not (isinstance(edges, list) and isinstance(fracs, list)):
                continue
            if edges != ref_edges or len(fracs) != nbins:
                continue

            w = payload.get("totalW")
            w = float(w) if isinstance(w, (int, float, np.integer, np.floating)) else 1.0  # fallback: equal weight
            if w <= 0:
                continue

            num += w * np.asarray(fracs, dtype=np.float64)
            den += w

        if den <= 0:
            continue

        fr_run3 = (num / den).tolist()
        by_ma[ma_tag] = {
            "bin_edges": ref_edges,
            "fractions": fr_run3,
            "totalW": den,  # optional: for debug/consistency
        }

    if by_ma:
        out[run3_label] = by_ma
    return out

def _plot_dr_fraction_bar_year(
    year: str,
    by_ma: Dict[str, Dict[str, float]],
    out_dir: Path,
    *,
    out_name: str,
) -> None:
    """
    每個 mA 一根水平 stacked bar：
      x: fraction (%), y: mA
    顏色依 dR bin 變化（從低到高）。
    """
    # rows: (ma_int, ma_tag, fracs, edges)
    rows = []
    for ma_tag, payload in (by_ma or {}).items():
        ma_int = _true_ma_from_dir(ma_tag)
        if ma_int is None:
            continue
        edges = payload.get("bin_edges")
        fracs = payload.get("fractions")
        if not isinstance(edges, list) or not isinstance(fracs, list):
            continue
        if len(edges) != len(fracs) + 1:
            continue
        if sum(fracs) <= 0:
            continue
        rows.append((ma_int, ma_tag, fracs, edges))

    if not rows:
        return

    rows.sort(key=lambda t: t[0])
    n = len(rows)

    # colors for bins (repeat if too many)
    palette = [
        "#5C7AFF", "#59D2FE", "#3BCEAC", "#FFD557", "#FFB640",
        "#FF688B", "#C77DFF", "#9CD1FF", "#417B5A", "#242038",
    ]

    c = ROOT.TCanvas(f"c_drfrac_{year}", "", 800, 800)
    ROOT.gStyle.SetOptStat(0)
    c.SetMargin(0.15, 0.06, 0.12, 0.12) # left, right, bottom, top
    c.SetTickx(); c.SetTicky()
    c.cd()

    frame = ROOT.TH2F(f"frame_drfrac_{year}", "", 1, 0.0, 100.0, n, 0.0, float(n))
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

    for i, (ma_int, _, _, _) in enumerate(rows, start=1):
        frame.GetYaxis().SetBinLabel(i, f"{ma_int}")
    frame.Draw("")
    if ROOT.gPad:
        ROOT.gPad.SetGridx(1)

    # legend: show a few representative bins (avoid ultra-wide legend)
    # use first row edges as reference
    _, _, fr0, edges0 = rows[0]
    nbins = len(fr0)
    leg = ROOT.TLegend(0.12, 0.915, 0.94, 0.94)
    leg.SetNColumns(min(5, max(1, nbins)))
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    leg.SetTextFont(42); leg.SetTextSize(0.05)

    dummy_boxes = []
    for ib in range(nbins):
        lo = float(edges0[ib]); hi = float(edges0[ib + 1])
        col = _root_color(palette[ib % len(palette)], fallback=ROOT.kAzure + 1)
        b = ROOT.TBox(0, 0, 0, 0)
        b.SetFillColor(col); b.SetLineColor(col)
        dummy_boxes.append(b)
        leg.AddEntry(b, f"{lo:.2f} ≤ ΔR < {hi:.2f}", "f")

    # NEW: percentage labels on bar segments (user coordinates, like dREffBar)
    pct = ROOT.TLatex()
    pct.SetTextFont(42)
    pct.SetTextSize(0.028)

    def _draw_pct(xl: float, xr: float, yy0: float, yy1: float, val: float):
        w = xr - xl
        if w < 6.0:  # too narrow in "percent units" (0~100)
            return
        pct.SetTextAlign(22)  # center
        pct.DrawLatex((xl + xr) / 2.0, (yy0 + yy1) / 2.0 - 0.02, f"{val:.1f}%")

    keep = [frame, leg, pct] + dummy_boxes

    bar_h = 0.70
    for idx, (_, _, fracs, edges) in enumerate(rows):
        y0 = float(idx) + (1.0 - bar_h) / 2.0
        y1 = float(idx) + 1.0 - (1.0 - bar_h) / 2.0

        x = 0.0
        for ib, frac in enumerate(fracs):
            w = 100.0 * float(frac)
            if w <= 0:
                continue
            col = _root_color(palette[ib % len(palette)], fallback=ROOT.kAzure + 1)
            bx = ROOT.TBox(x, y0, x + w, y1)
            bx.SetFillColor(col); bx.SetLineColor(col)
            bx.Draw("SAME")
            keep.append(bx)

            # NEW: draw percent label (skip if too narrow)
            _draw_pct(x, x + w, y0, y1, w)

            x += w

    # CMS + lumi
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(c.GetLeftMargin(), 0.985, "#bf{CMS} #it{Simulation}")
    keep.append(lat)

    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetTextSize(0.038)
        lat_lumi.DrawLatex(0.94, 0.985, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")
        keep.append(lat_lumi)

    leg.Draw()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = out_dir / out_name
    c.Modified(); c.Update()
    if ROOT.gPad:
        ROOT.gPad.RedrawAxis()
    c.SaveAs(str(out_pdf))
    try:
        c.Close()
    except Exception:
        pass
    _ = keep

def _plot_ma_migration_matrix_year(
    year: str,
    pred_to_true: Dict[int, Dict[int, float]],
    out_dir: Path,
    *,
    out_name: Optional[str] = None,
    zmin: float = 0.0,
    zmax: float = 100.0,
) -> None:
    """
    2D matrix: y=pred mA, x=true mA, content = percent (row-normalized per pred).
    Draw with COLZ+TEXT.
    """
    if not pred_to_true:
        return

    pred_mas = sorted(pred_to_true.keys())
    true_set = set()
    for p, tm in pred_to_true.items():
        true_set |= set((tm or {}).keys())
    true_mas = sorted(true_set)

    if not pred_mas or not true_mas:
        return

    nx = len(true_mas)
    ny = len(pred_mas)

    c = ROOT.TCanvas(f"c_migMat_{year}", "", 900, 820)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat("4.1f")  # one decimal
    c.SetMargin(0.12, 0.20, 0.14, 0.08)  # left, right, bottom, top
    c.SetTickx(); c.SetTicky()
    c.cd()

    h = ROOT.TH2F(f"h_migMat_{year}", "", nx, 0.0, float(nx), ny, 0.0, float(ny))
    h.SetDirectory(0)
    h.SetTitle("")
    h.GetXaxis().SetTitle("m_{a}^{True} (GeV)")
    h.GetYaxis().SetTitle("m_{a}^{Predict} (GeV)")
    h.GetZaxis().SetTitle("Fraction (%)")
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetZaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetTitleOffset(1.2)
    h.GetZaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetLabelSize(0.06)
    h.GetYaxis().SetLabelSize(0.06)
    h.GetZaxis().SetLabelSize(0.04)

    # axis bin labels
    for ix, tma in enumerate(true_mas, start=1):
        h.GetXaxis().SetBinLabel(ix, f"{tma}")
    for iy, pma in enumerate(pred_mas, start=1):
        h.GetYaxis().SetBinLabel(iy, f"{pma}")

    # fill percent: normalize per pred row
    for iy, pma in enumerate(pred_mas, start=1):
        tm = pred_to_true.get(pma) or {}
        tot = float(sum(tm.values()))
        if tot <= 0:
            continue
        for ix, tma in enumerate(true_mas, start=1):
            v = float(tm.get(tma, 0.0))
            h.SetBinContent(ix, iy, 100.0 * v / tot)

    h.SetMinimum(float(zmin))
    h.SetMaximum(float(zmax))

    try:
        ROOT.gStyle.SetPalette(ROOT.kBird)
    except Exception:
        pass

    h.SetMarkerSize(1.45)
    h.Draw("COLZ TEXT")
    if ROOT.gPad:
        ROOT.gPad.SetGrid(0, 0)

    # CMS + lumi
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.045)
    lat.DrawLatex(c.GetLeftMargin(), 0.9255, "#bf{CMS} #it{Simulation}")

    lumi_fb = _lumi_fb_for_year(year)
    if lumi_fb is not None:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetNDC()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetTextSize(0.04)
        lat_lumi.DrawLatex(1. - c.GetRightMargin(), 0.925, f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = out_dir / (out_name or f"mAmigration_matrix_{year}.pdf")
    c.Modified(); c.Update()
    c.SaveAs(str(out_pdf))
    try:
        c.Close()
    except Exception:
        pass

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

def main():
    parser = argparse.ArgumentParser(description="Plot mA migration (4-cat) and ΔR fraction bars per year.")
    parser.add_argument("--year", default="all", help="Year to plot (e.g. 2022preEE). Use 'all' to plot all years.")
    parser.add_argument("--pred-mas", default="same", help="Predicted mA list: 'same' or comma list like '1,2,3,...,30'.")
    parser.add_argument("--debug", action="store_true", help="Print diagnostics (why files/trees/branches are skipped).")  # NEW
    parser.add_argument("--max-files", type=int, default=None, help="Debug helper: only read first N ROOT files.")       # NEW
    parser.add_argument(
        "--score-prefix",
        default="auto",
        choices=["auto", "mva", "bdt"],
        help="Which score branch prefix to use for migration cuts: auto(prefer MVA), mva, bdt.",
    )  # NEW

    # NEW: eos path knobs (optional)
    parser.add_argument("--eos-user", default=None, help="EOS area letter for /eos/home-<letter>/... (default: $EOS_USER or 'p').")
    parser.add_argument("--eos-username", default=None, help="EOS username for /eos/home-*/<username>/... (default: $EOS_USERNAME or 'pelai').")

    args = parser.parse_args()

    # NEW: normalize input_root_dir to correct EOS format
    global input_root_dir
    if input_root_dir.startswith("/eos/home-") and "/HZa/root_P2Root/run3_mergedBDT" in input_root_dir:
        input_root_dir = _eos_path_default(
            "HZa/root_P2Root/run3_mergedBDT",
            eos_user=args.eos_user,
            eos_username=args.eos_username,
        )

    in_dir = Path(input_root_dir)
    out_dir = Path(output_plot_dir)

    # RENAME: clearer name, same payload/usage
    dr_frac_scan = _load_dr_fraction_scan_by_year_ma(
        Path(input_root_dir),
        years=YEAR_ORDER if args.year == "all" else [args.year],
        ma_dirs=signal_ma_dirs,
        tree_name="inclusive",
        dr_branch="var_dR_g1g2",
        dr_min=0.0,
        dr_max=4.0,
        dr_step=0.05,
        debug=args.debug,          # NEW
        max_files=args.max_files,  # NEW
    )

    if args.debug:
        print(f"[DEBUG] input_root_dir = {input_root_dir}")
        print(f"[DEBUG] optimized_BDT_Cut = {optimized_BDT_Cut}")


    # --- migration fractions by BDT_Score_mA_M{pred} passing pred cut ---
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

    pred_pass = _accumulate_pred_pass_sumw_by_year(
        Path(input_root_dir),
        years=years,
        true_ma_dirs=signal_ma_dirs,
        pred_mas=pred_mas,
        mva_cuts=mva_cuts,
        debug=args.debug,          # NEW
        max_files=args.max_files,  # NEW
        score_prefer=args.score_prefix,  # NEW
    )

    if args.debug and not pred_pass:
        print("[DEBUG] pred_pass is empty: likely原因=找不到tree/branch、mva cut缺、或全都沒通過cut(加權為0)。")

    for y in years:
        by_pred = pred_pass.get(y) or {}
        if by_pred:
            # _plot_ma_migration_4cat_bar_year(
            #     y,
            #     by_pred,
            #     out_dir,
            #     out_name=f"mAmigration_4cat_{y}.pdf",
            # )
            # # NEW: matrix plot (percent + numbers)
            _plot_ma_migration_matrix_year(
                y,
                by_pred,
                out_dir,
                out_name=f"mAmigration_matrix_{y}.pdf",
            )

    # NEW: Run3 combined migration plot (sumw across years)
    pred_run3 = _merge_pred_pass_years_to_run3(pred_pass, years, run3_label="Run3")
    by_pred_run3 = pred_run3.get("Run3") or {}
    if by_pred_run3:
        # _plot_ma_migration_4cat_bar_year(
        #     "Run3",
        #     by_pred_run3,
        #     out_dir,
        #     out_name="mAmigration_4cat_Run3.pdf",
        # )
        # # NEW: Run3 matrix plot
        _plot_ma_migration_matrix_year(
            "Run3",
            by_pred_run3,
            out_dir,
            out_name="mAmigration_matrix_Run3.pdf",
        )

if __name__ == "__main__":
    main()
