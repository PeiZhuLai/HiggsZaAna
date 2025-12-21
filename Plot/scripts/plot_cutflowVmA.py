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

baseDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/cutflow/cutflow_list"

outDir = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/cutflowVmA"

dataType = ["Data", "Bkg_MC", "Sig_MC"]

# Year of Signal
year_sig_2022 = ["2022preEE", "2022postEE"]
year_sig_2023 = ["2023preBPix", "2023postBPix"]
year_sig_2024 = ["2024"]
# Year of Bkg
year_DYG_2022 = ["2022preEE", "2022postEE"]
year_DYG_2023 = ["2023preBPix", "2023postBPix"]
year_DYG_2024 = ["2024"]
years_DYJet_2022 = ["2022preEE","2022postEE"]
year_DYJet_2023  = ["2023preBPix", "2023postBPix"]
years_DYJet_2024 = ["2024"]
# Year of Data
years_Data_2022 = ["2022preEE","2022postEE"]
year_Data_2023  = ["2023preBPix", "2023postBPix"]
years_Data_2024 = ["2024"]

# Name of Signal Sample
name_sig_2022 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
name_sig_2023 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
name_sig_2024 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
# Name of Bkg Sample
name_DYG_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
name_DYG_2023 = ["DYGto2LG_10to100"]
name_DYG_2024 = ["DYGto2LG_10to100"]
name_DYJet_2022 = ["DYJetsToLL"]
name_DYJet_2023 = ["DYJetsToLL"]
name_DYJet_2024 = ["DYJetsTo2E","DYJetsTo2Mu","DYJetsTo2Tau"]
# Name of Data Sample
name_Data_2022 = ["DYJetsToLL"]
name_Data_2023 = ["DYJetsToLL"]
name_Data_2024 = ["DYJetsToLL"]





# ===== Samples / 常數 (整合後唯一版本) =====
years_sig  = ["2022preEE"]          # 信號年份
ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
sig_samples = ["ALP_M5", "ALP_M15", "ALP_M30"]
INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_BDT"
optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

YEAR = "2022preEE"
XS_PB = 0.1
BR = 1.0
KEEP_FB = 1000.0
TEST_TO_ALL = 2.0
ALT_WEIGHT_CANDS = ["weight","genWeight","eventWeight"]
INPUT_BASE_TREE_NAME = "test"

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, "2024": 108.95,
            'combined_run3':170.84 }
LUMI_FB = float(lumiMap[YEAR])

# ==== Helpers ====
# 新增：ROOT 風格常數（對齊 plot_effisigma.py）
OFFSET = 0.01
MARKER_SIZE = 1.3
LINE_WIDTH = 3
# 新增：固定 y 軸範圍
Y_MIN = 0.8
Y_MAX = 3.5

# ====== 解析 BDT 門檻 ======
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

def _parse_mva_cuts(path: str) -> Dict[int,float]:
    """
    從 JSON 檔案讀取各 mA 的 MVAcut，容錯支援：
    - dict: { "results": [ {mA/ma/mass: ..., MVAcut/mvaCut/cut: ...}, ... ] }
    - list: [ { ... }, { ... } ]
    - dict: { "ALP_M5": 0.975, "ALP_M15": 0.98, ... }
    - dict: { "5": {"MVAcut": 0.975}, "15": {"cut": 0.98} }
    - 字串型數值亦可，例如 "5 GeV", "0.975 +/- 0.01"
    解析失敗時回傳 {}，並在主程式顯示錯誤訊息。
    """
    def _entry_to_pair(entry: dict) -> Tuple[Optional[int], Optional[float]]:
        # 優先直接取數值欄位
        mass_keys = ("mA", "ma", "mass", "massA")
        cut_keys  = ("MVAcut", "mvaCut", "mva_cut", "cut", "bdtCut", "bdt_cut", "best_MVAcut")
        ma = None
        for k in mass_keys:
            if k in entry:
                ma = _to_int(entry[k]); break
        # 從 sample/name 類欄位抽出，例如 "ALP_M5"
        if ma is None:
            for k in ("sample", "name", "label", "title"):
                if k in entry and isinstance(entry[k], str):
                    ma = _parse_ma(entry[k])
                    if ma is not None:
                        break
        # cut
        cut = None
        for k in cut_keys:
            if k in entry:
                cut = _to_float(entry[k]); break
        # 若 cut 在次層
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

    # 1) list 形式
    if isinstance(data, list):
        for ent in data:
            if not isinstance(ent, dict): continue
            ma, cut = _entry_to_pair(ent)
            if ma is not None and cut is not None:
                out[ma] = cut

    # 2) dict 且含 results/points/entries 清單
    if not out and isinstance(data, dict):
        for key in ("results", "points", "entries"):
            arr = data.get(key)
            if isinstance(arr, list):
                for ent in arr:
                    if not isinstance(ent, dict): continue
                    ma, cut = _entry_to_pair(ent)
                    if ma is not None and cut is not None:
                        out[ma] = cut
                if out:
                    break

    # 3) dict 的巢狀容器，例如 data["2022preEE"]["results"]
    if not out and isinstance(data, dict):
        for v in data.values():
            if isinstance(v, dict):
                for key in ("results", "points", "entries"):
                    arr = v.get(key)
                    if isinstance(arr, list):
                        for ent in arr:
                            if not isinstance(ent, dict): continue
                            ma, cut = _entry_to_pair(ent)
                            if ma is not None and cut is not None:
                                out[ma] = cut
                        if out:
                            break
            if out:
                break
    return out

# ====== 樣本與質量對應 ======
def _parse_ma(name: str) -> Optional[int]:
    import re
    m = re.search(r"[Mm](\d+)", name)
    return int(m.group(1)) if m else None

def _build_mass_map(samples: List[str]) -> Dict[int,str]:
    d = {}
    for s in samples:
        ma = _parse_ma(s)
        if ma is not None:
            d[ma] = s
    return d

# ====== 檔案列舉 ======
def _list_root_files(base: str, sample: str, years: List[str]) -> List[str]:
    files=[]
    sdir = os.path.join(base.rstrip("/"), sample)
    if not os.path.isdir(sdir):
        return files
    # 優先 {sample}/{year}.root
    for y in years:
        p = os.path.join(sdir, f"{y}.root")
        if os.path.isfile(p):
            files.append(p)
    if not files:  # fallback: search all
        for r,_,fnames in os.walk(sdir):
            for fn in fnames:
                if not fn.endswith(".root"): continue
                full = os.path.join(r,fn)
                if years and not any(y in fn for y in years): continue
                files.append(full)
    # 去重並排序
    out=[]
    seen=set()
    for f in files:
        if f not in seen:
            seen.add(f); out.append(f)
    return sorted(out)

def _pick_weight_branch(t) -> Optional[str]:
    keys = set(map(str, t.keys()))
    for k in ALT_WEIGHT_CANDS:
        if k in keys: return k
    return None

def _branch_exists(t, name: str) -> bool:
    return name in set(map(str, t.keys()))

# ====== 計算某個質量的 (w_pass_mu, w_pass_ele) ======
def _accumulate_pass_weights(ma: int,
                             sample_map: Dict[int,str],
                             mva_cuts: Dict[int,float]) -> Tuple[float,float]:
    if ma not in sample_map or ma not in mva_cuts:
        return 0.0,0.0
    sample = sample_map[ma]
    cut = mva_cuts[ma]
    files = _list_root_files(INPUT_BASE, sample, years_sig)
    w_pass_mu = 0.0
    w_pass_ele = 0.0
    for fp in files:
        try:
            with uproot.open(fp) as f:
                if INPUT_BASE_TREE_NAME not in f: continue
                t = f[INPUT_BASE_TREE_NAME]
                if not _branch_exists(t,"MVA_Score"): continue
                wname = _pick_weight_branch(t)
                # 動態要讀的分支
                branches = ["MVA_Score"]
                if wname: branches.append(wname)
                has_mu = _branch_exists(t,"z_mumu")
                has_ele = _branch_exists(t,"z_ee")
                if has_mu: branches.append("z_mumu")
                if has_ele: branches.append("z_ee")
                for arrs in t.iterate(branches, library="ak", step_size="200 MB"):
                    mva = arrs["MVA_Score"]
                    mask = mva >= cut
                    if wname:
                        w = ak.values_astype(arrs[wname], np.float64)
                    else:
                        w = ak.ones_like(mva, dtype=np.float64)
                    if has_mu:
                        mu_mask = mask & (ak.values_astype(arrs["z_mumu"], np.int8)==1)
                        w_pass_mu += float(ak.sum(w[mu_mask]))
                    if has_ele:
                        ele_mask = mask & (ak.values_astype(arrs["z_ee"], np.int8)==1)
                        w_pass_ele += float(ak.sum(w[ele_mask]))
        except Exception:
            continue
    return w_pass_mu, w_pass_ele

# ====== sumw -> Aε(%) 計算 ======
def _sumw_to_eff_percent(w_pass: float) -> float:
    denom = XS_PB * BR * LUMI_FB * KEEP_FB
    if denom <= 0: return 0.0
    return 100.0 * TEST_TO_ALL * w_pass / denom

def _build_sumw_series(channel: str,
                       ma_order: List[int],
                       sample_map: Dict[int,str],
                       mva_cuts: Dict[int,float]) -> Tuple[List[int],List[float]]:
    xs, ys = [], []
    for ma in ma_order:
        if ma not in sample_map or ma not in mva_cuts:  # 沒 sample 或沒 cut
            continue
        w_mu, w_ele = _accumulate_pass_weights(ma, sample_map, mva_cuts)
        w_use = w_mu if channel=="muon" else w_ele
        if w_use > 0:
            xs.append(ma)
            ys.append(_sumw_to_eff_percent(w_use))
    return xs, ys

def _quadratic_curve(xs: List[int], ys: List[float], n: int = 400) -> Tuple[np.ndarray, np.ndarray]:
    """
    產生二次插值的平滑曲線點集。
    優先使用 SciPy 的 make_interp_spline(k=2)；若不可用，則以滑動三點的局部二次擬合補上。
    """
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    # 確保 x 遞增
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if _HAVE_SCIPY and len(x) >= 3:
        spline = make_interp_spline(x, y, k=2)
        x_new = np.linspace(x.min(), x.max(), n)
        y_new = spline(x_new)
        return x_new, y_new

    # 後備方案：局部二次擬合，區段在相鄰點的中點間拼接
    if len(x) < 3:
        x_new = np.linspace(x.min(), x.max(), n)
        y_new = np.interp(x_new, x, y)
        return x_new, y_new

    x_new_list, y_new_list = [], []
    seg_pts = max(8, n // max(1, len(x) - 2))
    for i in range(len(x) - 2):
        xi = x[i:i+3]
        yi = y[i:i+3]
        coeffs = np.polyfit(xi, yi, deg=2)

        left = x[i] if i == 0 else 0.5 * (x[i] + x[i+1])
        right = x[i+2] if i == len(x) - 3 else 0.5 * (x[i+1] + x[i+2])

        xs_seg = np.linspace(left, right, seg_pts)
        ys_seg = np.polyval(coeffs, xs_seg)

        if x_new_list and xs_seg[0] <= x_new_list[-1][-1]:
            mask = xs_seg > x_new_list[-1][-1]
            xs_seg = xs_seg[mask]
            ys_seg = ys_seg[mask]

        x_new_list.append(xs_seg)
        y_new_list.append(ys_seg)

    x_new = np.concatenate(x_new_list) if x_new_list else x
    y_new = np.concatenate(y_new_list) if y_new_list else y
    return x_new, y_new

def _plot_line(x: List[int], y: List[float], channel_label: str, out_dir: Path) -> Path:
    from array import array as carray

    # 平滑二次曲線與原始點資料
    x_smooth, y_smooth = _quadratic_curve(x, y, n=400)
    xs = list(map(float, x))
    ys = list(map(float, y))
    sx = list(map(float, x_smooth))
    sy = list(map(float, y_smooth))

    # 建圖（線 + 點）
    g_line = ROOT.TGraph(len(sx), carray('d', sx), carray('d', sy))
    g_pts  = ROOT.TGraph(len(xs), carray('d', xs), carray('d', ys))

    # 使用十六進位色碼設定顏色
    color = ROOT.TColor.GetColor("#F564A9")
    marker = 20

    g_line.SetTitle("")
    g_line.SetLineColor(color)
    g_line.SetMarkerColor(color)
    g_line.SetLineWidth(LINE_WIDTH)
    g_line.SetMarkerStyle(marker)  # 供 legend 顯示

    g_pts.SetLineColor(0)
    g_pts.SetMarkerColor(color)
    g_pts.SetMarkerStyle(marker)
    g_pts.SetMarkerSize(MARKER_SIZE)

    # 畫布與邊界
    c = ROOT.TCanvas(f"c_{channel_label.lower()}", "", 800, 600)
    c.SetMargin(0.12 + OFFSET, 0.035, 0.14, 0.09)
    c.SetTickx()
    c.SetTicky()

    # 軸與樣式
    g_line.GetXaxis().SetTitle("m_{a} [GeV]")
    g_line.GetYaxis().SetTitle("Efficiency #times Acceptance (%)")
    g_line.Draw("AL")

    ax, ay = g_line.GetXaxis(), g_line.GetYaxis()
    for axis in (ax, ay):
        axis.CenterTitle(True)
        axis.SetTitleFont(42)
        axis.SetLabelFont(42)
        axis.SetTitleSize(0.055)
        axis.SetLabelSize(0.05)
    ax.SetTitleOffset(1.15)
    ay.SetTitleOffset(1.1 + 10*OFFSET)
    ax.SetLabelOffset(0.009)

    # 座標範圍：x 固定；y 改為固定範圍
    ax.SetLimits(0.0, 33.0)
    g_line.SetMinimum(Y_MIN)
    g_line.SetMaximum(Y_MAX)

    # 原始點
    g_pts.Draw("P SAME")

    # Legend
    leg = ROOT.TLegend(0.65, 0.70, 0.93, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.SetHeader(f"{channel_label} Channel", "L")
    leg.AddEntry(g_line, "2022preEE", "lp")
    leg.Draw()

    # CMS Preliminary 與 13.6 TeV
    lat = ROOT.TLatex()
    lat.SetTextFont(42)
    # 改：先留左側 CMS Preliminary 用左下對齊
    lat.SetTextAlign(11)
    lat.SetNDC()
    lat.SetTextSize(0.05)
    lat.DrawLatex(0.12 + OFFSET, 0.92, "#bf{CMS} #it{Preliminary}")
    # 新：右上角亮度改用右上對齊，避免被右邊界裁切
    lat.SetTextAlign(31)  # 3=Right, 1=Top
    x_lumi = 1.0 - c.GetRightMargin() - 0.005
    y_lumi = 0.92
    lat.DrawLatexNDC(x_lumi, y_lumi, "61.89 fb^{-1} (13.6 TeV)")

    c.Update()
    c.RedrawAxis()  # 確保標註與座標完整重繪

    # 輸出
    fname = f"sigEfficiencyVmA_{'muon' if channel_label.lower().startswith('muon') else 'ele'}"
    out_path = out_dir / f"{fname}.png"
    c.SaveAs(str(out_path))
    c.SaveAs(str(out_path.with_suffix(".pdf")))
    return out_path

def main():
    # 1. 解析各質量 BDT cut
    mva_cuts = _parse_mva_cuts(optimized_BDT_Cut)
    if not mva_cuts:
        print(f"[錯誤] 從 JSON 讀不到任何 MVA cut：{optimized_BDT_Cut}，請檢查格式/鍵名。")
        return
    else:
        print(f"[資訊] 成功載入 {len(mva_cuts)} 個 MVA cut：{sorted(mva_cuts.items())}")

    # 2. 建立 ma->sample 對應
    sample_map = _build_mass_map(sig_samples)
    out_dir = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/signal_eff_sumw")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 3. 計算並繪圖 (muon)
    x_mu, y_mu = _build_sumw_series("muon", ma_list, sample_map, mva_cuts)
    if x_mu:
        _plot_line(x_mu, y_mu, "Muon", out_dir)

    # 4. (electron)
    x_ele, y_ele = _build_sumw_series("ele", ma_list, sample_map, mva_cuts)
    if x_ele:
        _plot_line(x_ele, y_ele, "Electron", out_dir)

if __name__ == "__main__":
    main()
