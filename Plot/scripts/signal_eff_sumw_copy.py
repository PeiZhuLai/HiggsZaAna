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

# ===== Samples / 常數 (整合後唯一版本) =====
years_sig  = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]  # 信号
ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
# sig_samples = ["ALP_M5", "ALP_M15", "ALP_M30"]
sig_samples = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_BDT"
optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

YEAR = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
XS_PB = 0.1
BR = 1.0
KEEP_FB = 1000.0
TEST_TO_ALL = 2.0
ALT_WEIGHT_CANDS = ["weight","genWeight","eventWeight"]
SYS_Ele = ['weight_hlt_sf_up','weight_hlt_sf_down','weight_pu_reweight_sf_up','weight_pu_reweight_sf_down','weight_electron_wplid_sf_SelectedElectron_up','weight_electron_wplid_sf_SelectedElectron_down', 'weight_electron_reco_sf_SelectedElectron_up', 'weight_electron_reco_sf_SelectedElectron_down', 'weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_up', 'weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_down','weight_photon_id_sf_SelectedPhoton_up','weight_photon_id_sf_SelectedPhoton_down']
SYS_Ele_central = ['weight_hlt_sf_central','weight_pu_reweight_sf_central','weight_electron_wplid_sf_SelectedElectron_central', 'weight_electron_reco_sf_SelectedElectron_central', 'weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_central', 'weight_photon_id_sf_SelectedPhoton_central']
SYS_Mu = ['weight_hlt_sf_up','weight_hlt_sf_down','weight_pu_reweight_sf_up','weight_pu_reweight_sf_down','weight_muon_looseid_sf_SelectedMuon_up', 'weight_muon_looseid_sf_SelectedMuon_down', 'weight_muon_iso_sf_SelectedMuon_up', 'weight_muon_iso_sf_SelectedMuon_down', 'weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_up', 'weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_down','weight_photon_id_sf_SelectedPhoton_up','weight_photon_id_sf_SelectedPhoton_down']
SYS_Mu_central = ['weight_hlt_sf_central','weight_pu_reweight_sf_central','weight_muon_looseid_sf_SelectedMuon_central', 'weight_muon_iso_sf_SelectedMuon_central', 'weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_central','weight_photon_id_sf_SelectedPhoton_central']

INPUT_BASE_TREE_NAME = "test"

lumiMap = { '16':16.81,'16APV':19.52,'17':41.48,'18':59.83,'combined':137.65,
            '2022preEE':7.98,'2022postEE':26.70,'2023preBPix':17.79,'2023postBPix':9.45, '2024':108.95, 
            'Run3':170.84 }

# 改：LUMI_FB 應根據 years_sig 動態加總，而非固定 YEAR
def _sum_lumi_fb(years: List[str]) -> float:
    miss = [y for y in years if y not in lumiMap]
    if miss:
        raise KeyError(f"lumiMap 缺少年份鍵：{miss}")
    return float(sum(lumiMap[y] for y in years))

try:
    LUMI_FB = _sum_lumi_fb(years_sig)
except Exception:
    # fallback：若 years_sig 剛好是整包 run3 或 lumiMap 不完整，可退回 combined_run3
    LUMI_FB = float(lumiMap.get("Run3", 0.0))

# ==== Helpers ====
# 新增：ROOT 風格常數（對齊 plot_effisigma.py）
OFFSET = 0.01
MARKER_SIZE = 1.3
LINE_WIDTH = 3
# 新增：固定 y 軸範圍
Y_MIN = 0.0
Y_MAX_GROUPS = 2.5   # 3條線(2022/2023/2024)用
Y_MAX_5YEARS = 2.7   # 5條線(逐年)用：請改成你想要的值

# 新增：error bar 最小可見門檻（單位：y 軸的百分比）
ERR_ABS_FLOOR = 0.02   # 例如 0.02 (%)；太小會看不到
ERR_REL_FLOOR = 0.02   # 例如 2% 的相對誤差：err >= eff * 0.02

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
                             mva_cuts: Dict[int,float],
                             years: List[str]) -> Tuple[float,float,float,float]:
    if ma not in sample_map or ma not in mva_cuts:
        return 0.0, 0.0, 0.0, 0.0
    sample = sample_map[ma]
    cut = mva_cuts[ma]
    files = _list_root_files(INPUT_BASE, sample, years)
    w_pass_mu = 0.0
    w_pass_ele = 0.0
    w2_pass_mu = 0.0
    w2_pass_ele = 0.0
    for fp in files:
        try:
            with uproot.open(fp) as f:
                if INPUT_BASE_TREE_NAME not in f: continue
                t = f[INPUT_BASE_TREE_NAME]
                if not _branch_exists(t,"MVA_Score"): continue
                wname = _pick_weight_branch(t)
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

                    # 統計誤差：同時累積 sumw2
                    w2 = w * w

                    if has_mu:
                        mu_mask = mask & (ak.values_astype(arrs["z_mumu"], np.int8) == 1)
                        w_pass_mu += float(ak.sum(w[mu_mask]))
                        w2_pass_mu += float(ak.sum(w2[mu_mask]))
                    if has_ele:
                        ele_mask = mask & (ak.values_astype(arrs["z_ee"], np.int8) == 1)
                        w_pass_ele += float(ak.sum(w[ele_mask]))
                        w2_pass_ele += float(ak.sum(w2[ele_mask]))
        except Exception:
            continue
    return w_pass_mu, w_pass_ele, w2_pass_mu, w2_pass_ele

# --- 新增：系統通道的 sumw/sumw2（用 w_eff * weight） ---
def _accumulate_pass_sumw_systs(ma: int,
                                sample_map: Dict[int, str],
                                mva_cuts: Dict[int, float],
                                years: List[str],
                                sys_names: List[str],
                                sys_central_names: Optional[List[str]] = None) -> Dict[str, Tuple[float, float, float, float]]:
    """
    回傳 dict: sys_name -> (sumw_pass_mu, sumw_pass_ele, sumw2_pass_mu, sumw2_pass_ele)

    w_eff 的定義沿用原本邏輯，但同時累積 sumw2 = sum(w_eff^2)，
    讓誤差 bar 可以反映「sys 權重造成的統計波動」。
    """
    if ma not in sample_map or ma not in mva_cuts:
        return {}

    sample = sample_map[ma]
    cut = mva_cuts[ma]
    files = _list_root_files(INPUT_BASE, sample, years)

    out: Dict[str, Tuple[float, float, float, float]] = {s: (0.0, 0.0, 0.0, 0.0) for s in sys_names}

    base_to_central: Dict[str, str] = {}
    if sys_central_names:
        for c in sys_central_names:
            base = re.sub(r"_central$", "", c)
            base_to_central[base] = c

    for fp in files:
        try:
            with uproot.open(fp) as f:
                if INPUT_BASE_TREE_NAME not in f:
                    continue
                t = f[INPUT_BASE_TREE_NAME]
                if not _branch_exists(t, "MVA_Score"):
                    continue

                has_mu = _branch_exists(t, "z_mumu")
                has_ele = _branch_exists(t, "z_ee")

                sys_branches = [s for s in sys_names if _branch_exists(t, s)]
                if not sys_branches:
                    continue

                wname = _pick_weight_branch(t)

                central_branches = []
                if base_to_central:
                    for s in sys_branches:
                        base = re.sub(r"_(up|down)$", "", s)
                        c = base_to_central.get(base)
                        if c and _branch_exists(t, c):
                            central_branches.append(c)
                central_branches = sorted(set(central_branches))

                branches = ["MVA_Score"] + sys_branches + central_branches
                if wname:
                    branches.append(wname)
                if has_mu:
                    branches.append("z_mumu")
                if has_ele:
                    branches.append("z_ee")

                for arrs in t.iterate(branches, library="ak", step_size="200 MB"):
                    mva = arrs["MVA_Score"]
                    mask = mva >= cut

                    mu_mask = None
                    ele_mask = None
                    if has_mu:
                        mu_mask = mask & (ak.values_astype(arrs["z_mumu"], np.int8) == 1)
                    if has_ele:
                        ele_mask = mask & (ak.values_astype(arrs["z_ee"], np.int8) == 1)

                    if wname:
                        w_nom = ak.values_astype(arrs[wname], np.float64)
                    else:
                        w_nom = ak.ones_like(mva, dtype=np.float64)

                    for s in sys_branches:
                        w_sys_raw = ak.values_astype(arrs[s], np.float64)

                        base = re.sub(r"_(up|down)$", "", s)
                        c_name = base_to_central.get(base) if base_to_central else None

                        w_eff = None

                        # --- 優先 central ratio ---
                        if c_name and c_name in arrs.fields:
                            w_c = ak.values_astype(arrs[c_name], np.float64)
                            ok = ak.isfinite(w_c) & (w_c != 0) & ak.isfinite(w_sys_raw)
                            if ak.any(ok):
                                ratio = ak.ones_like(w_sys_raw, dtype=np.float64)
                                ratio = ak.where(ok, w_sys_raw / w_c, ratio)
                                w_eff = w_nom * ratio

                        # --- fallback：ratio/absolute 偵測 ---
                        if w_eff is None:
                            probe = w_sys_raw[mask]
                            is_ratio = False
                            try:
                                probe = probe[ak.isfinite(probe)]
                                if len(probe) > 0:
                                    med = float(ak.to_numpy(ak.median(probe)))
                                    if 0.2 < med < 5.0:
                                        q10 = float(np.quantile(ak.to_numpy(probe), 0.10))
                                        q90 = float(np.quantile(ak.to_numpy(probe), 0.90))
                                        if 0.0 < q10 and q90 < 10.0:
                                            is_ratio = True
                            except Exception:
                                is_ratio = False

                            w_eff = (w_nom * w_sys_raw) if is_ratio else w_sys_raw

                        w_eff2 = w_eff * w_eff

                        mu_sum, ele_sum, mu_sum2, ele_sum2 = out.get(s, (0.0, 0.0, 0.0, 0.0))
                        if mu_mask is not None:
                            mu_sum += float(ak.sum(w_eff[mu_mask]))
                            mu_sum2 += float(ak.sum(w_eff2[mu_mask]))
                        if ele_mask is not None:
                            ele_sum += float(ak.sum(w_eff[ele_mask]))
                            ele_sum2 += float(ak.sum(w_eff2[ele_mask]))
                        out[s] = (mu_sum, ele_sum, mu_sum2, ele_sum2)
        except Exception:
            continue

    out = {k: v for k, v in out.items() if (v[0] != 0.0 or v[1] != 0.0)}
    return out

# ====== sumw -> Aε(%) 計算 ======
def _sumw_to_eff_percent(w_pass: float) -> float:
    denom = XS_PB * BR * LUMI_FB * KEEP_FB
    if denom <= 0: return 0.0
    return 100.0 * TEST_TO_ALL * w_pass / denom

def _sumw2_to_efferr_percent(w2_pass: float) -> float:
    denom = XS_PB * BR * LUMI_FB * KEEP_FB
    if denom <= 0: return 0.0
    if w2_pass <= 0: return 0.0
    return 100.0 * TEST_TO_ALL * (float(np.sqrt(w2_pass)) / denom)

def _build_sumw_series(channel: str,
                       ma_order: List[int],
                       sample_map: Dict[int,str],
                       mva_cuts: Dict[int,float],
                       years: List[str]) -> Tuple[List[int],List[float],List[float]]:
    xs, ys, yerrs = [], [], []
    for ma in ma_order:
        if ma not in sample_map or ma not in mva_cuts:
            continue

        w_mu, w_ele, w2_mu, w2_ele = _accumulate_pass_weights(ma, sample_map, mva_cuts, years)
        if channel == "muon":
            w_use, w2_use = w_mu, w2_mu
        else:
            w_use, w2_use = w_ele, w2_ele

        if w_use <= 0:
            continue

        eff_nom = _sumw_to_eff_percent(w_use)
        stat_err = _sumw2_to_efferr_percent(w2_use)

        sys_names = SYS_Mu if channel == "muon" else SYS_Ele
        sys_central = SYS_Mu_central if channel == "muon" else SYS_Ele_central
        sumw_systs = _accumulate_pass_sumw_systs(ma, sample_map, mva_cuts, years, sys_names, sys_central)

        # envelope by base (up/down max), 同時蒐集 sys 統計誤差（用 w_eff 的 sumw2）
        deltas_by_base: Dict[str, float] = {}
        max_sys_stat_delta = 0.0

        for sys_name, (sw_mu, sw_ele, sw2_mu, sw2_ele) in sumw_systs.items():
            sw = sw_mu if channel == "muon" else sw_ele
            sw2 = sw2_mu if channel == "muon" else sw2_ele
            if sw == 0.0:
                continue

            eff_sys = _sumw_to_eff_percent(sw)
            delta = float(eff_sys - eff_nom)

            # sys 統計誤差：由 w_eff 的 sumw2 來
            stat_err_sys = _sumw2_to_efferr_percent(sw2)
            max_sys_stat_delta = max(max_sys_stat_delta, abs(stat_err_sys - stat_err))

            base = re.sub(r"_(up|down)$", "", sys_name)
            prev = deltas_by_base.get(base, 0.0)
            if abs(delta) > abs(prev):
                deltas_by_base[base] = delta

        syst2 = float(sum(d*d for d in deltas_by_base.values()))
        syst_err = float(np.sqrt(syst2)) if syst2 > 0 else 0.0

        # 新 error bar：把 sys 權重造成的統計波動也納入
        tot_err = float(np.sqrt(stat_err * stat_err + syst_err * syst_err + max_sys_stat_delta * max_sys_stat_delta))

        # 新增：避免 error bar 小到看不見（視覺化 floor，不改 central value）
        vis_floor = max(float(ERR_ABS_FLOOR), float(abs(eff_nom) * ERR_REL_FLOOR))
        tot_err = max(tot_err, vis_floor)

        xs.append(ma)
        ys.append(eff_nom)
        yerrs.append(tot_err)

    return xs, ys, yerrs

def _build_sumw_series_by_year(channel: str,
                               ma_order: List[int],
                               sample_map: Dict[int,str],
                               mva_cuts: Dict[int,float],
                               years_list: List[str]) -> Dict[str, Tuple[List[int], List[float], List[float]]]:
    out: Dict[str, Tuple[List[int], List[float], List[float]]] = {}
    for y in years_list:
        x, yy, ye = _build_sumw_series(channel, ma_order, sample_map, mva_cuts, years=[y])
        if x:
            out[y] = (x, yy, ye)
    return out

def _build_sumw_series_by_year_group(channel: str,
                                     ma_order: List[int],
                                     sample_map: Dict[int, str],
                                     mva_cuts: Dict[int, float],
                                     year_groups: Dict[str, List[str]]) -> Dict[str, Tuple[List[int], List[float], List[float]]]:
    out: Dict[str, Tuple[List[int], List[float], List[float]]] = {}
    for label, years in year_groups.items():
        x, yy, ye = _build_sumw_series(channel, ma_order, sample_map, mva_cuts, years=years)
        if x:
            out[label] = (x, yy, ye)
    return out

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

def _plot_lines_by_year(series_by_year: Dict[str, Tuple[List[int], List[float], List[float]]],
                        channel_label: str,
                        out_dir: Path,
                        name_suffix: str = "",
                        y_min: float = Y_MIN,
                        y_max: float = Y_MAX_GROUPS) -> Optional[Path]:
    from array import array as carray
    if not series_by_year:
        return None

    palette_hex_ele = [
        "#540D6E", "#EE4266", "#FFB640",
        "#3BCEAC", "#0EAD69",
    ]

    palette_hex_mu = [
        "#5C7AFF", "#242038", "#59D2FE",
        "#44E5E7", "#73FBD3",
    ]

    # 新增：依 channel 選擇調色盤（Electron/Muon 顏色不同）
    palette_hex = palette_hex_mu if channel_label.lower().startswith("muon") else palette_hex_ele

    line_syles = [1, 2, 3]
    marker_styles = [20,  21,  23,  33,  34,  47,   29]
    marker_size   = [1.4, 1.3, 1.6, 1.9, 1.7, 1.6, 1.9]

    years = list(series_by_year.keys())
    years.sort()  # 讓顏色/legend 順序穩定（可移除）

    c = ROOT.TCanvas(f"c_{channel_label.lower()}_by_year", "", 800, 600)
    c.SetMargin(0.12 + OFFSET, 0.035, 0.14, 0.09)
    c.SetTickx()
    c.SetTicky()

    # legend lines = 1(header) + 1 * number of actually drawable year-series (one combined entry per year)
    n_drawable = 0
    for y in years:
        x, yy, ye = series_by_year[y]
        if len(x) >= 2:
            n_drawable += 1
    n_lines = 1 + 1 * n_drawable

    # geometry: keep same top-right anchor, vary height with lines
    x1, y1 = 0.95, 0.85
    x0 = 0.63
    line_h = 0.052
    pad = 0.018
    h_leg = pad + n_lines * line_h
    y0 = max(0.14, y1 - h_leg)

    leg = ROOT.TLegend(x0, y0, x1, y1)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)

    if channel_label == "Muon":
        header = "H #rightarrow Za #rightarrow #mu^{+}#mu^{-} + 2#gamma"
    else:
        header = "H #rightarrow Za #rightarrow e^{+}e^{-} + 2#gamma"

    leg.SetHeader(header, "L")

    first = True
    g_store = []  # 避免被 GC

    # FIX: years 是 list[str]，需 enumerate 才能拿到 (index, year)
    for i, y in enumerate(years):
        x, yy, ye = series_by_year[y]
        if len(x) < 2:
            continue

        # styles (不足就循環)
        color = ROOT.TColor.GetColor(palette_hex[i % len(palette_hex)])
        ls = line_syles[i % len(line_syles)]
        ms = marker_styles[i % len(marker_styles)]
        msz = marker_size[i % len(marker_size)]

        # 1) smooth line (不加誤差)
        x_smooth, y_smooth = _quadratic_curve(x, yy, n=400)
        sx = list(map(float, x_smooth))
        sy = list(map(float, y_smooth))

        g_line = ROOT.TGraph(len(sx), carray('d', sx), carray('d', sy))
        g_line.SetTitle("")
        g_line.SetLineColor(color)
        g_line.SetLineWidth(LINE_WIDTH)
        g_line.SetLineStyle(ls)
        g_line.SetMarkerStyle(1)       # 線圖不需要 marker
        g_line.SetMarkerSize(0)

        draw_opt = "AL" if first else "L SAME"
        g_line.Draw(draw_opt)
        first = False
        g_store.append(g_line)

        # 2) original points + error bars
        px = list(map(float, x))
        py = list(map(float, yy))
        pye = list(map(float, ye))
        pxe = [0.0] * len(px)

        g_pts = ROOT.TGraphErrors(
            len(px),
            carray('d', px), carray('d', py),
            carray('d', pxe), carray('d', pye),
        )
        g_pts.SetTitle("")
        g_pts.SetLineColor(color)
        g_pts.SetLineWidth(LINE_WIDTH)   # 讓 legend 顯示線（用同樣粗細）
        g_pts.SetLineStyle(ls)           # 讓 legend 線型跟平滑線一致
        g_pts.SetMarkerColor(color)
        g_pts.SetMarkerStyle(ms)
        g_pts.SetMarkerSize(msz)
        # ROOT 版本相容：error bar 端點帽大小
        if hasattr(g_pts, "SetEndErrorSize"):
            g_pts.SetEndErrorSize(4)

        g_pts.Draw("P SAME")
        g_store.append(g_pts)

        # legend：單一 entry 同時顯示線+點（不要再分 line/pts）
        leg.AddEntry(g_pts, f"{y}", "lp")

    if not g_store:
        return None

    # 軸（用第一條線的 axis）
    g0 = g_store[0]
    g0.GetXaxis().SetTitle("m_{a} [GeV]")
    g0.GetYaxis().SetTitle("Efficiency #times Acceptance (%)")
    ax, ay = g0.GetXaxis(), g0.GetYaxis()
    for axis in (ax, ay):
        axis.CenterTitle(False)
        axis.SetTitleFont(42)
        axis.SetLabelFont(42)
        axis.SetTitleSize(0.055)
        axis.SetLabelSize(0.05)
    ax.SetTitleOffset(1.15)
    ay.SetTitleOffset(1.1 + 10*OFFSET)
    ax.SetLabelOffset(0.009)

    ax.SetLimits(0.0, 32.0)
    g0.SetMinimum(y_min)
    g0.SetMaximum(y_max)

    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetTextFont(42)
    lat.SetTextAlign(11)
    lat.SetNDC()
    lat.SetTextSize(0.05)
    lat.DrawLatex(0.12 + OFFSET, 0.92, "#bf{CMS} #it{Simulation}")
    lat.SetTextAlign(31)
    x_lumi = 1.0 - c.GetRightMargin() - 0.005
    y_lumi = 0.92
    lat.DrawLatexNDC(x_lumi, y_lumi, "13.6 TeV")

    c.Update()
    c.RedrawAxis()

    fname = f"sigEfficiencyVmA_{'muon' if channel_label.lower().startswith('muon') else 'ele'}_byYear{name_suffix}"
    out_path = out_dir / f"{fname}.pdf"
    c.SaveAs(str(out_path))
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

    # 3. 計算並繪圖：每個年份群組各一條線 (2022/2023/2024)
    year_groups = {
        "2022": ["2022preEE", "2022postEE"],
        "2023": ["2023preBPix", "2023postBPix"],
        "2024": ["2024"],
    }

    mu_series = _build_sumw_series_by_year_group("muon", ma_list, sample_map, mva_cuts, year_groups)
    _plot_lines_by_year(mu_series, "Muon", out_dir, y_max=Y_MAX_GROUPS)

    ele_series = _build_sumw_series_by_year_group("ele", ma_list, sample_map, mva_cuts, year_groups)
    _plot_lines_by_year(ele_series, "Electron", out_dir, y_max=Y_MAX_GROUPS)

    # 4. 新增：逐子年份（5 條線）各輸出一張圖
    mu_series_5 = _build_sumw_series_by_year("muon", ma_list, sample_map, mva_cuts, years_sig)
    _plot_lines_by_year(mu_series_5, "Muon", out_dir, name_suffix="_5years", y_max=Y_MAX_5YEARS)

    ele_series_5 = _build_sumw_series_by_year("ele", ma_list, sample_map, mva_cuts, years_sig)
    _plot_lines_by_year(ele_series_5, "Electron", out_dir, name_suffix="_5years", y_max=Y_MAX_5YEARS)

if __name__ == "__main__":
    main()
