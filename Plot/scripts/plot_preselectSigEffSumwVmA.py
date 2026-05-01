# ==== New imports/config for plotting ====
import os
from pathlib import Path
# 先安全匯入 PyROOT，避免 Cling 解析到重複 LLVM 參數
import ROOT
import numpy as np
import re
from typing import Dict, List, Tuple, Optional
import uproot
import awkward as ak
# 加回：SciPy 嘗試導入 (供 _quadratic_curve 使用)
try:
    from scipy.interpolate import make_interp_spline
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

# ===== Debug switches =====
DEBUG = False
DEBUG_MA = None          # e.g. 4; None -> all
DEBUG_MAX_FILES = 2      # per mA
DEBUG_MAX_CHUNKS = 2     # per file
DEBUG_PRINT_FIELDS_SAMPLE = 40  # print first N branch names

def _dbg(msg: str):
    if DEBUG:
        print(msg)

# ===== Samples / 常數 (整合後唯一版本) =====
years_sig  = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]  # 信号
ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
ma_interpolate = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
# sig_samples = ["ALP_M5", "ALP_M15", "ALP_M30"]
sig_samples = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal"
# optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
# 改：不再使用 optimized_BDT_Cut；改用固定 MVA cut（可自行調整或之後加 argparse）
DEFAULT_MVA_CUT = None  # type: Optional[float]

YEAR = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
XS_PB = 0.1
BR = 1.0
KEEP_FB = 1000.0
TEST_TO_ALL = 2.0
ALT_WEIGHT_CANDS = ["weight","genWeight","eventWeight"]
SYS_Ele = []  # will be filled by _discover_sys_branches() in main()
SYS_Ele_central = []  # will be filled by _discover_sys_branches() in main()
SYS_Mu = []  # will be filled by _discover_sys_branches() in main()
SYS_Mu_central = []  # will be filled by _discover_sys_branches() in main()

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

# 改：不要用全域 LUMI_FB 固定加總 years_sig；改成「每次用到時依 years 計算」
def _lumi_fb_for_years(years: List[str]) -> float:
    miss = [y for y in years if y not in lumiMap]
    if miss:
        raise KeyError(f"lumiMap 缺少年份鍵：{miss}")
    return float(sum(lumiMap[y] for y in years))

# （移除/停用舊的 LUMI_FB 固定值，避免逐年時計算錯分母）
# try:
#     LUMI_FB = _sum_lumi_fb(years_sig)
# except Exception:
#     LUMI_FB = float(lumiMap.get("Run3", 0.0))

# ==== Helpers ====
# 新增：ROOT 風格常數（對齊 plot_effisigma.py）
OFFSET = 0.01
MARKER_SIZE = 1.3
LINE_WIDTH = 3
# 新增：固定 y 軸範圍
Y_MIN = 0.0
Y_MAX_GROUPS = 12   # 3條線(2022/2023/2024)用
Y_MAX_5YEARS = 12   # 5條線(逐年)用：請改成你想要的值

# 新增：error bar 最小可見門檻（單位：y 軸的百分比）
ERR_ABS_FLOOR = 0.02   # 例如 0.02 (%)；太小會看不到
ERR_REL_FLOOR = 0.02   # 例如 2% 的相對誤差：err >= eff * 0.02
USE_VIS_FLOOR = False  # True 會強制 error bar 至少顯示到 floor（只為視覺化，不建議用於正式不確定度）

# 新增：是否列印系統誤差 breakdown（每個 mA 都會印，量可能很大）
PRINT_SYST_BREAKDOWN = False # False
# 新增：是否也列印每個 systematic base 的分量（更長）
PRINT_SYST_COMPONENTS = False # False

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

# 新增：允許 cut=None（代表不切）
def _build_fixed_mva_cuts(ma_values: List[int], cut: Optional[float]) -> Dict[int, Optional[float]]:
    return {int(m): (None if cut is None else float(cut)) for m in ma_values}

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
        if DEBUG:
            _dbg(f"[DBG] sample dir not found: {sdir}")
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
    out = sorted(out)
    if DEBUG:
        _dbg(f"[DBG] _list_root_files sample={sample} years={years} nfiles={len(out)}")
        for i, ff in enumerate(out[:DEBUG_MAX_FILES]):
            _dbg(f"      file[{i}] = {ff}")
        if len(out) > DEBUG_MAX_FILES:
            _dbg(f"      ... (showing first {DEBUG_MAX_FILES})")
    return out

def _discover_sys_branches(sample_map: Dict[int, str], years: List[str]) -> Tuple[List[str], List[str]]:
    """
    自動從一個可用的 ROOT 檔案中掃描 systematics weight branches。
    回傳：
      - sys_updown: 例如 weight_xxx_up / weight_xxx_down
      - sys_central: 例如 weight_xxx_central（若存在）
    好處：避免手動維護 SYS_* 清單，也避免因為字串被截斷導致 branch 名稱不匹配。
    """
    # 找到第一個能開的檔案當作「schema 代表」
    rep_file = None
    rep_sample = None
    for ma, sample in sample_map.items():
        files = _list_root_files(INPUT_BASE, sample, years)
        if files:
            rep_file = files[0]
            rep_sample = sample
            break
    if rep_file is None:
        print("[警告] 找不到任何輸入 ROOT 檔案，無法自動掃描 systematics branches。")
        return [], []

    try:
        with uproot.open(rep_file) as f:
            if INPUT_BASE_TREE_NAME not in f:
                print(f"[警告] {rep_file} 沒有 tree={INPUT_BASE_TREE_NAME}，無法掃描 systematics branches。")
                return [], []
            t = f[INPUT_BASE_TREE_NAME]
            fields = list(getattr(t, "keys", lambda: [])())
            # uproot v4: t.keys() 回傳 branch 名稱
            # 只抓 weight_* 且帶 _up/_down 的 branches；central 會另外嘗試從同一個 base 名稱推回來
            sys_updown = [b for b in fields if b.startswith("weight_") and re.search(r"_(up|down)$", b)]
            sys_updown = sorted(set(sys_updown))

            # central/nominal 的 branch 命名在不同 ntuple 裡常不一樣：
            #   - 有的叫 base（沒有 suffix）
            #   - 有的叫 base_central
            #   - 有的叫 base_nominal
            bases = sorted(set(re.sub(r"_(up|down)$", "", b) for b in sys_updown))
            sys_central = []
            for base in bases:
                for cand in (base, base + "_central", base + "_nominal"):
                    if cand in fields:
                        sys_central.append(cand)
            sys_central = sorted(set(sys_central))
            print(f"[資訊] 自動掃描 systematics branches：sample={rep_sample}, file={os.path.basename(rep_file)}")
            print(f"       up/down branches: {len(sys_updown)}  |  central branches: {len(sys_central)}")
            return sys_updown, sys_central
    except Exception as e:
        print(f"[警告] 掃描 systematics branches 失敗：{e}")
        return [], []


def _split_sys_branches_by_channel(sys_updown: List[str], sys_central: List[str]) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    依 branch 名稱將 systematics weight 分成：
      - Electron channel: (common + electron-only)
      - Muon channel:     (common + muon-only)
    使用者需求：沒有寫 electron 或 muon 的 → 兩個 channel 都要考慮。
    """
    def _tag(b: str) -> str:
        bl = b.lower()
        has_e = 'electron' in bl
        has_m = 'muon' in bl
        if has_e and has_m:
            return 'common'
        if has_e:
            return 'ele'
        if has_m:
            return 'mu'
        return 'common'

    common_ud, ele_ud, mu_ud = [], [], []
    for b in sys_updown:
        t = _tag(b)
        if t == 'ele':
            ele_ud.append(b)
        elif t == 'mu':
            mu_ud.append(b)
        else:
            common_ud.append(b)

    common_c, ele_c, mu_c = [], [], []
    for b in sys_central:
        t = _tag(b)
        if t == 'ele':
            ele_c.append(b)
        elif t == 'mu':
            mu_c.append(b)
        else:
            common_c.append(b)

    sys_ele = sorted(set(common_ud + ele_ud))
    sys_mu  = sorted(set(common_ud + mu_ud))
    cen_ele = sorted(set(common_c + ele_c))
    cen_mu  = sorted(set(common_c + mu_c))
    return sys_ele, sys_mu, cen_ele, cen_mu

def _pick_weight_branch(t) -> Optional[str]:
    keys = set(map(str, t.keys()))
    for k in ALT_WEIGHT_CANDS:
        if k in keys: return k
    return None

def _branch_exists(t, name: str) -> bool:
    return name in set(map(str, t.keys()))

def _pick_mva_score_branch(t, sample: str) -> Optional[str]:
    """
    優先讀 era 內分 mA 存的 branch：MVA_Score_{sample}（e.g. MVA_Score_mA_M4）
    若不存在則 fallback 到舊的 MVA_Score。
    """
    keys = set(map(str, t.keys()))
    cand = f"MVA_Score_{sample}"
    if cand in keys:
        return cand
    if "MVA_Score" in keys:
        return "MVA_Score"
    return None

# ====== 計算某個質量的 (w_pass_mu, w_pass_ele) ======
def _accumulate_pass_weights(ma: int,
                             sample_map: Dict[int,str],
                             mva_cuts: Dict[int,float],
                             years: List[str]) -> Tuple[float,float,float,float,float,float,float,float]:
    """
    回傳：
      (sumw_pass_mu, sumw_pass_ele, sumw2_pass_mu, sumw2_pass_ele,
       sumw_tot_mu,  sumw_tot_ele,  sumw2_tot_mu,  sumw2_tot_ele)
    """
    if ma not in sample_map or ma not in mva_cuts:
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    sample = sample_map[ma]
    cut = mva_cuts[ma]
    files = _list_root_files(INPUT_BASE, sample, years)

    if DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
        _dbg(f"[DBG] ==== _accumulate_pass_weights ma={ma} sample={sample} years={years} cut={cut} files={len(files)} ====")

    w_pass_mu = w_pass_ele = 0.0
    w2_pass_mu = w2_pass_ele = 0.0
    w_tot_mu = w_tot_ele = 0.0
    w2_tot_mu = w2_tot_ele = 0.0

    for ifile, fp in enumerate(files):
        if DEBUG_MAX_FILES is not None and ifile >= int(DEBUG_MAX_FILES) and DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
            _dbg(f"[DBG] stop after DEBUG_MAX_FILES={DEBUG_MAX_FILES}")
            break
        try:
            with uproot.open(fp) as f:
                if INPUT_BASE_TREE_NAME not in f:
                    if DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
                        _dbg(f"[DBG] skip (no tree): {fp}")
                    continue
                t = f[INPUT_BASE_TREE_NAME]

                if DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
                    try:
                        fields = list(t.keys())
                        _dbg(f"[DBG] opened file={os.path.basename(fp)} tree={INPUT_BASE_TREE_NAME} nbranches={len(fields)}")
                        _dbg(f"[DBG] branches(head)={fields[:min(len(fields), DEBUG_PRINT_FIELDS_SAMPLE)]}")
                    except Exception as e:
                        _dbg(f"[DBG] cannot list branches: {e}")

                mva_branch = None
                if cut is not None:
                    mva_branch = _pick_mva_score_branch(t, sample)
                    if not mva_branch:
                        if DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
                            _dbg(f"[DBG] skip (no MVA branch) sample={sample} file={os.path.basename(fp)}")
                        continue

                wname = _pick_weight_branch(t)

                has_mu = _branch_exists(t, "z_mumu")
                has_ele = _branch_exists(t, "z_ee")

                if DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
                    _dbg(f"[DBG] mva_branch={mva_branch} weight_branch={wname} has_mu={has_mu} has_ele={has_ele}")

                branches = []
                if mva_branch:
                    branches.append(mva_branch)
                if wname:
                    branches.append(wname)

                if has_mu:
                    branches.append("z_mumu")
                if has_ele:
                    branches.append("z_ee")

                for ichunk, arrs in enumerate(t.iterate(branches, library="ak", step_size="200 MB")):
                    if DEBUG_MAX_CHUNKS is not None and ichunk >= int(DEBUG_MAX_CHUNKS) and DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
                        _dbg(f"[DBG] stop after DEBUG_MAX_CHUNKS={DEBUG_MAX_CHUNKS}")
                        break

                    if mva_branch:
                        mva = arrs[mva_branch]
                        pass_mask = (mva >= float(cut))  # cut not None
                    else:
                        if wname and wname in arrs.fields:
                            base = arrs[wname]
                        elif has_mu and "z_mumu" in arrs.fields:
                            base = arrs["z_mumu"]
                        elif has_ele and "z_ee" in arrs.fields:
                            base = arrs["z_ee"]
                        else:
                            continue
                        pass_mask = ak.ones_like(base, dtype=np.bool_)

                    if wname and wname in arrs.fields:
                        w = ak.values_astype(arrs[wname], np.float64)
                    else:
                        w = ak.ones_like(pass_mask, dtype=np.float64)
                    w2 = w * w

                    if DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
                        n_evt = int(len(w)) if hasattr(w, "__len__") else -1
                        n_pass = int(ak.sum(pass_mask))
                        wsum = float(ak.sum(w))
                        wsum_pass = float(ak.sum(w[pass_mask]))
                        _dbg(f"[DBG] chunk={ichunk} nevt={n_evt} npass={n_pass} sumw={wsum:.6g} sumw_pass={wsum_pass:.6g}")

                    if has_mu and "z_mumu" in arrs.fields:
                        mu_sel = (ak.values_astype(arrs["z_mumu"], np.int8) == 1)
                        mu_pass = pass_mask & mu_sel
                        w_tot_mu += float(ak.sum(w[mu_sel]))
                        w2_tot_mu += float(ak.sum(w2[mu_sel]))
                        w_pass_mu += float(ak.sum(w[mu_pass]))
                        w2_pass_mu += float(ak.sum(w2[mu_pass]))

                    if has_ele and "z_ee" in arrs.fields:
                        ele_sel = (ak.values_astype(arrs["z_ee"], np.int8) == 1)
                        ele_pass = pass_mask & ele_sel
                        w_tot_ele += float(ak.sum(w[ele_sel]))
                        w2_tot_ele += float(ak.sum(w2[ele_sel]))
                        w_pass_ele += float(ak.sum(w[ele_pass]))
                        w2_pass_ele += float(ak.sum(w2[ele_pass]))
        except Exception as e:
            if DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
                _dbg(f"[DBG] exception reading {fp}: {e}")
            continue

    if DEBUG and (DEBUG_MA is None or ma == DEBUG_MA):
        _dbg(f"[DBG] RESULT ma={ma} w_pass_mu={w_pass_mu:.6g} w_tot_mu={w_tot_mu:.6g} | w_pass_ele={w_pass_ele:.6g} w_tot_ele={w_tot_ele:.6g}")

    return (w_pass_mu, w_pass_ele, w2_pass_mu, w2_pass_ele,
            w_tot_mu,  w_tot_ele,  w2_tot_mu,  w2_tot_ele)

# --- 新增：系統通道的 pass/total sumw/sumw2（用 w_eff） ---
def _accumulate_sumw_systs_pass_total(ma: int,
                                     sample_map: Dict[int, str],
                                     mva_cuts: Dict[int, float],
                                     years: List[str],
                                     sys_names: List[str],
                                     sys_central_names: Optional[List[str]] = None
                                     ) -> Dict[str, Tuple[float, float, float, float, float, float, float, float]]:
    """
    回傳 dict: sys_name ->
      (sumw_pass_mu, sumw_pass_ele, sumw2_pass_mu, sumw2_pass_ele,
       sumw_tot_mu,  sumw_tot_ele,  sumw2_tot_mu,  sumw2_tot_ele)
    """
    if ma not in sample_map or ma not in mva_cuts:
        return {}

    sample = sample_map[ma]
    cut = mva_cuts[ma]
    files = _list_root_files(INPUT_BASE, sample, years)

    out: Dict[str, Tuple[float, float, float, float, float, float, float, float]] = {
        s: (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) for s in sys_names
    }

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

                mva_branch = _pick_mva_score_branch(t, sample)
                if not mva_branch:
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

                branches = [mva_branch] + sys_branches + central_branches
                if wname:
                    branches.append(wname)
                if has_mu:
                    branches.append("z_mumu")
                if has_ele:
                    branches.append("z_ee")

                for arrs in t.iterate(branches, library="ak", step_size="200 MB"):
                    mva = arrs[mva_branch]
                    pass_mask = mva >= cut

                    mu_sel = None
                    ele_sel = None
                    if has_mu:
                        mu_sel = (ak.values_astype(arrs["z_mumu"], np.int8) == 1)
                    if has_ele:
                        ele_sel = (ak.values_astype(arrs["z_ee"], np.int8) == 1)

                    if wname:
                        w_nom = ak.values_astype(arrs[wname], np.float64)
                    else:
                        w_nom = ak.ones_like(mva, dtype=np.float64)

                    for s in sys_branches:
                        w_sys_raw = ak.values_astype(arrs[s], np.float64)
                        base = re.sub(r"_(up|down)$", "", s)
                        c_name = base_to_central.get(base) if base_to_central else None

                        w_eff = None

                        if c_name and c_name in arrs.fields:
                            w_c = ak.values_astype(arrs[c_name], np.float64)
                            ok = ak.isfinite(w_c) & (w_c != 0) & ak.isfinite(w_sys_raw)
                            if ak.any(ok):
                                ratio = ak.ones_like(w_sys_raw, dtype=np.float64)
                                ratio = ak.where(ok, w_sys_raw / w_c, ratio)
                                w_eff = w_nom * ratio

                        if w_eff is None:
                            probe = w_sys_raw[pass_mask]
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

                        (p_mu, p_ele, p2_mu, p2_ele,
                         t_mu, t_ele, t2_mu, t2_ele) = out.get(s, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

                        if mu_sel is not None:
                            mu_pass = pass_mask & mu_sel
                            t_mu += float(ak.sum(w_eff[mu_sel]))
                            t2_mu += float(ak.sum(w_eff2[mu_sel]))
                            p_mu += float(ak.sum(w_eff[mu_pass]))
                            p2_mu += float(ak.sum(w_eff2[mu_pass]))

                        if ele_sel is not None:
                            ele_pass = pass_mask & ele_sel
                            t_ele += float(ak.sum(w_eff[ele_sel]))
                            t2_ele += float(ak.sum(w_eff2[ele_sel]))
                            p_ele += float(ak.sum(w_eff[ele_pass]))
                            p2_ele += float(ak.sum(w_eff2[ele_pass]))

                        out[s] = (p_mu, p_ele, p2_mu, p2_ele, t_mu, t_ele, t2_mu, t2_ele)
        except Exception:
            continue

    out = {k: v for k, v in out.items() if (v[0] != 0.0 or v[1] != 0.0 or v[4] != 0.0 or v[5] != 0.0)}
    return out

def _build_sumw_series(channel: str,
                       ma_order: List[int],
                       sample_map: Dict[int,str],
                       mva_cuts: Dict[int,float],
                       years: List[str]) -> Tuple[List[int], List[float], List[float], List[float]]:
    """
    回傳 (xs, ys, yerr_lo, yerr_hi)

    改：y 軸改成「sum of weight 正規化」：
      y = sumw_pass / (XS_PB * BR * KEEP_FB * lumi(years))
    統計誤差用 sqrt(sumw2_pass) 做同樣正規化。
    """
    xs, ys, yerr_lo, yerr_hi = [], [], [], []

    lumi_fb = _lumi_fb_for_years(years)
    denom = float(XS_PB) * float(BR) * float(KEEP_FB) * float(lumi_fb)
    if DEBUG:
        _dbg(f"[DBG] _build_sumw_series channel={channel} years={years} lumi_fb={lumi_fb} denom={denom}")
    if denom <= 0:
        return xs, ys, yerr_lo, yerr_hi

    for ma in ma_order:
        if DEBUG_MA is not None and ma != DEBUG_MA:
            continue
        if ma not in sample_map or ma not in mva_cuts:
            if DEBUG:
                _dbg(f"[DBG] skip ma={ma} (not in sample_map or mva_cuts)")
            continue

        (w_pass_mu, w_pass_ele, w2_pass_mu, w2_pass_ele,
         w_tot_mu,  w_tot_ele,  w2_tot_mu,  w2_tot_ele) = _accumulate_pass_weights(ma, sample_map, mva_cuts, years)

        if channel == "muon":
            w_pass = w_pass_mu
            w2_pass = w2_pass_mu
        else:
            w_pass = w_pass_ele
            w2_pass = w2_pass_ele

        # 改：轉成 percent
        y = 100.0 * (2.0 * float(w_pass) / denom)
        stat = 100.0 * (float(np.sqrt(max(w2_pass, 0.0))) / denom)

        if DEBUG:
            _dbg(f"[DBG] ma={ma} channel={channel} w_pass={w_pass:.6g} w2_pass={w2_pass:.6g} -> y={y:.6g} stat={stat:.6g}")

        tot_up = stat
        tot_dn = stat

        if USE_VIS_FLOOR:
            vis_floor = max(float(ERR_ABS_FLOOR), float(abs(y) * ERR_REL_FLOOR))
            tot_up = max(tot_up, vis_floor)
            tot_dn = max(tot_dn, vis_floor)

        xs.append(ma)
        ys.append(y)
        yerr_lo.append(tot_dn)
        yerr_hi.append(tot_up)

    return xs, ys, yerr_lo, yerr_hi

def _build_sumw_series_by_year(channel: str,
                               ma_order: List[int],
                               sample_map: Dict[int,str],
                               mva_cuts: Dict[int,float],
                               years_list: List[str]) -> Dict[str, Tuple[List[int], List[float], List[float], List[float]]]:
    out: Dict[str, Tuple[List[int], List[float], List[float], List[float]]] = {}
    for y in years_list:
        x, yy, ylo, yhi = _build_sumw_series(channel, ma_order, sample_map, mva_cuts, years=[y])
        if x:
            out[y] = (x, yy, ylo, yhi)
    return out
def _build_sumw_series_by_year_group(channel: str,
                                     ma_order: List[int],
                                     sample_map: Dict[int, str],
                                     mva_cuts: Dict[int, float],
                                     year_groups: Dict[str, List[str]]) -> Dict[str, Tuple[List[int], List[float], List[float], List[float]]]:
    out: Dict[str, Tuple[List[int], List[float], List[float], List[float]]] = {}
    for label, years in year_groups.items():
        x, yy, ylo, yhi = _build_sumw_series(channel, ma_order, sample_map, mva_cuts, years=years)
        if x:
            out[label] = (x, yy, ylo, yhi)
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

def _plot_lines_by_year(series_by_year: Dict[str, Tuple[List[int], List[float], List[float], List[float]]],
                        channel_label: str,
                        out_dir: Path,
                        name_suffix: str = "",
                        y_min: float = Y_MIN,
                        y_max: float = Y_MAX_GROUPS) -> Optional[Path]:
    from array import array as carray
    if not series_by_year:
        return None

    palette_hex_mu = [
        "#540D6E", "#EE4266", "#FFB640",
        "#3BCEAC", "#086788",
    ]

    palette_hex_ele = [
        "#5C7AFF", "#242038", "#59D2FE",
        "#44E5E7", "#417B5A",
    ]

    # 新增：依 channel 選擇調色盤（Electron/Muon 顏色不同）
    palette_hex = palette_hex_mu if channel_label.lower().startswith("muon") else palette_hex_ele

    line_syles = [1, 2, 3]
    marker_styles = [20,  21,  23,  33,  34,  47,   29]
    marker_size   = [1.4, 1.3, 1.6, 1.9, 1.7, 1.6, 1.9]

    years = list(series_by_year.keys())
    # years.sort()  # 會導致 2022postEE 排到 2022preEE 前面，改用固定順序

    # 新增：固定 legend 顯示順序（只影響文字/entry 順序，不改點的樣式設定邏輯）
    _LEGEND_YEAR_ORDER = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
    years_set = set(years)
    years_in_order = [y for y in _LEGEND_YEAR_ORDER if y in years_set]
    years_rest = [y for y in years if y not in set(_LEGEND_YEAR_ORDER)]  # 保底：其他 key 維持原本插入順序
    years = years_in_order + years_rest

    # 新增：收集「平滑線在 ma_interpolate 上」的 y 值（只在 5years 模式輸出）
    interp_payload = {
        "channel": channel_label,
        "name_suffix": name_suffix,
        "ma_points": list(map(int, ma_interpolate)),
        "unit": "Efficiency x Acceptance (%)",
        "values": {},  # year -> { "11": y, ... }
    }

    c = ROOT.TCanvas(f"c_{channel_label.lower()}_by_year", "", 800, 600)
    c.SetMargin(0.12 + OFFSET, 0.035, 0.14, 0.09)
    c.SetTickx()
    c.SetTicky()

    # legend lines = 1(header) + 1 * number of actually drawable year-series (one combined entry per year)
    n_drawable = 0
    for y in years:
        x, yy, ylo, yhi = series_by_year[y]
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
        x, yy, ylo, yhi = series_by_year[y]
        if len(x) < 2:
            continue

        # styles (不足就循環)
        color = ROOT.TColor.GetColor(palette_hex[i % len(palette_hex)])
        ls = line_syles[i % len(line_syles)]
        ms = marker_styles[i % len(marker_styles)]
        msz = marker_size[i % len(marker_size)]

        # 1) smooth line (不加誤差)
        x_smooth, y_smooth = _quadratic_curve(x, yy, n=400)

        # 新增：輸出 ma_interpolate 的 y 值（用平滑線做線性內插）
        try:
            x_s = np.asarray(x_smooth, dtype=float)
            y_s = np.asarray(y_smooth, dtype=float)
            xmin, xmax = float(np.min(x_s)), float(np.max(x_s))

            per_year_items = []  # (ma_int, value) 先收集再排序，確保輸出順序 1,2,3,...
            for m in ma_interpolate:
                mf = float(m)
                if mf < xmin or mf > xmax:
                    continue
                per_year_items.append((int(m), float(np.interp(mf, x_s, y_s))))

            if per_year_items:
                per_year_items.sort(key=lambda t: t[0])
                interp_payload["values"][y] = {str(k): v for k, v in per_year_items}
        except Exception:
            pass

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
        pye_lo = list(map(float, ylo))
        pye_hi = list(map(float, yhi))
        pxe_lo = [0.0] * len(px)
        pxe_hi = [0.0] * len(px)

        g_pts = ROOT.TGraphAsymmErrors(
            len(px),
            carray('d', px), carray('d', py),
            carray('d', pxe_lo), carray('d', pxe_hi),
            carray('d', pye_lo), carray('d', pye_hi),
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
    # 改：y axis title 加上 (%)
    g0.GetYaxis().SetTitle("Preselection Signal Efficiency [%]")
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

    # 改：輸出檔名避免沿用 Eff
    fname = f"preselectSigEffSumwVmA_{'muon' if channel_label.lower().startswith('muon') else 'ele'}_byYear{name_suffix}"
    out_path = out_dir / f"{fname}.pdf"
    c.SaveAs(str(out_path))

    # 新增：在 5years 模式把 interpolated y 寫成 JSON（每個 channel 一份）
    # if name_suffix == "_5years" and interp_payload["values"]:
    #     json_name = f"{fname}_quadratic_interp_ma_points.json"
    #     jsonOutDir = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output")
    #     json_path = jsonOutDir / json_name
    #     meta = {
    #         "input_years": years,
    #         "mva_cut_json": optimized_BDT_Cut,
    #         "input_base": INPUT_BASE,
    #         "tree": INPUT_BASE_TREE_NAME,
    #         "xs_pb": XS_PB,
    #         "br": BR,
    #         "keep_fb": KEEP_FB,
    #         "test_to_all": TEST_TO_ALL,
    #     }
    #     try:
    #         with open(json_path, "w") as jf:
    #             # 重點：不要 sort_keys，不然 "10" 會又排到 "2" 前面
    #             json.dump({"meta": meta, **interp_payload}, jf, indent=2, sort_keys=False)
    #         print(f"[資訊] 已輸出 quadratic curve 內插點 JSON：{json_path}")
    #     except Exception as e:
    #         print(f"[警告] 寫出 JSON 失敗：{json_path} ({e})")

    return out_path

def main():
    # 改：preselect 不做 MVA cut（全通過）
    mva_cuts = _build_fixed_mva_cuts(ma_list, DEFAULT_MVA_CUT)
    print("[資訊] preselect 模式：不做 MVA cut（全通過）")

    if DEBUG:
        _dbg(f"[DBG] years_sig={years_sig}")
        _dbg(f"[DBG] XS_PB={XS_PB} BR={BR} KEEP_FB={KEEP_FB} INPUT_BASE={INPUT_BASE} tree={INPUT_BASE_TREE_NAME}")
        _dbg(f"[DBG] lumi(Run3 from map)={lumiMap.get('Run3')}  sum_lumi(years_sig)={_lumi_fb_for_years(years_sig)}")

    # 2. 建立 ma->sample 對應
    sample_map = _build_mass_map(sig_samples)

    # 2.1 自動掃描 systematics branches（避免手動 SYS_* 清單與 branch 名稱不一致）
    global SYS_Ele, SYS_Ele_central, SYS_Mu, SYS_Mu_central
    sys_updown, sys_central = _discover_sys_branches(sample_map, years_sig)
    # 依 branch 名稱把 e/μ 專屬系統分開；沒有寫 electron/muon 的視為 common（兩個 channel 都要考慮）
    SYS_Ele, SYS_Mu, SYS_Ele_central, SYS_Mu_central = _split_sys_branches_by_channel(sys_updown, sys_central)
    print(f"[資訊] syst branches split: common+ele={len(SYS_Ele)} | common+mu={len(SYS_Mu)} | central(ele)={len(SYS_Ele_central)} | central(mu)={len(SYS_Mu_central)}")
    out_dir = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/preselectSigEffSumwVmA")
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