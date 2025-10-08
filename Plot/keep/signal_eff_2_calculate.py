# ==== New imports/config for plotting ====
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot
import re
import ast
from typing import Dict, List, Tuple, Optional
import awkward as ak
import json  # 新增：讀入 TOTAL_JSON
from datetime import datetime, timezone  # 新增：輸出時間戳

# 你之前实际使用的组合：
years_sig  = ["2022preEE"]  # 信号

# 扫描的 ma 列表（背景也按你的原脚本扫）
ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]

# ===== Samples =====
total_sig_samples = ["HZaTo2l2g_M5_NanoAODv12", "HZaTo2l2g_M15_NanoAODv12", "HZaTo2l2g_M30_NanoAODv12"]
TOTAL_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sig_total_stats.json"

INPUT_BASE = "/eos/home-p/pelai/HZa/Root_Dataset/run3_BDT/"
sig_samples = ["ALP_M5", "ALP_M15", "ALP_M30"]

OUTPUT_BASE = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sig_eff.json"

optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/BDT_cut_all_ma_run3.txt"

# 如需其它：["ZGToLLG","WGToLNuG","ZG2JToG2L2J","EWKZ2J","TT","TTGJets","TGJets","ttWJets","ttZJets","WW","WZ","ZZ"]

MVA_CANDIDATES = ["MVA_Score"]
WEIGHT_CANDIDATES = ["weight"]

# ==== Helpers ====
# 设置全局字体大小和图形样式
plt.rcParams.update({
    'figure.figsize': (8, 6),
    'font.size': 14,
    'axes.titlesize': 17,
    'axes.labelsize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'figure.titlesize': 14,
    # 'patch.linewidth': 1.5,
    # 'patch.edgecolor': 'blue',
})

TOTAL_BASE_TREE_NAME = "Events"
INPUT_BASE_TREE_NAME = "test"
INPUT_BASE_Total_TREE_NAME = "inclusive"

ALT_WEIGHT_CANDS = ["weight", "genWeight", "eventWeight"]
UPROOT_STEP = "200 MB"  # iterate 每批大小

def parse_mva_cuts(txt_path: str) -> Dict[int, float]:
    with open(txt_path, "r") as f:
        content = f.read()
    m = re.search(r"\{.*\}", content, flags=re.S)
    if not m:
        return {}
    return {int(k): float(v) for k, v in ast.literal_eval(m.group(0)).items()}

def parse_ma_from_name(name: str) -> Optional[int]:
    m = re.search(r"[Mm](\d+)", name)
    return int(m.group(1)) if m else None

def build_ma_sample_maps(sig_samples: List[str], total_sig_samples: List[str]) -> Tuple[Dict[int,str], Dict[int,str]]:
    ma_to_input = {}
    for s in sig_samples:
        ma = parse_ma_from_name(s)
        if ma is not None:
            ma_to_input[ma] = s
    ma_to_total = {}
    for s in total_sig_samples:
        ma = parse_ma_from_name(s)
        if ma is not None:
            ma_to_total[ma] = s
    return ma_to_input, ma_to_total

def list_root_files(base: str, sample: str, years: List[str]) -> List[str]:
    # 僅在 {base}/{sample} 底下找；優先 {base}/{sample}/{year}.root；否則只收檔名含 year 的 .root
    files: List[str] = []
    base = base.rstrip("/")
    sample_dir = os.path.join(base, sample)

    if not os.path.exists(sample_dir):
        return files

    # 1) 優先使用單檔 {sample}/{year}.root
    for y in years:
        single = os.path.join(sample_dir, f"{y}.root")
        if os.path.isfile(single):
            files.append(single)

    # 2) 若還需要，掃描 sample 目錄，但僅保留檔名包含 year 的 .root
    if not files:
        for root, _, fnames in os.walk(sample_dir):
            for fn in fnames:
                if not fn.endswith(".root"):
                    continue
                full = os.path.join(root, fn)
                if years and not any((y in fn) or full.endswith(f"{y}.root") or (os.path.sep + y + os.path.sep in full) for y in years):
                    continue
                files.append(full)

    return sorted(list(dict.fromkeys(files)))

def pick_weight_branch(t) -> Optional[str]:
    keys = set(map(str, t.keys()))
    for k in ALT_WEIGHT_CANDS:
        if k in keys:
            return k
    return None

def get_branch_exists(t, name: str) -> bool:
    return name in set(map(str, t.keys()))

def analyze_input_files(files: List[str], mva_cut: float) -> Tuple[int, float, int, float, int, float, int, float]:
    # (n_input_total, w_input_total, n_pass_input, w_pass_input,
    #  n_pass_mu, w_pass_mu, n_pass_ele, w_pass_ele)
    # 分子來自 test；分母優先 inclusive，缺少則回退 test
    n_inclusive_tot = 0
    w_inclusive_tot = 0.0
    n_pass = 0
    w_pass = 0.0
    # 新增：channel-wise
    n_pass_mu = 0
    w_pass_mu = 0.0
    n_pass_ele = 0
    w_pass_ele = 0.0
    for fpath in files:
        try:
            with uproot.open(fpath) as f:
                incl_exists = INPUT_BASE_Total_TREE_NAME in f

                # 分母：inclusive（若存在）
                if incl_exists:
                    t_incl = f[INPUT_BASE_Total_TREE_NAME]
                    n_inclusive_tot += int(t_incl.num_entries)
                    wname_incl = pick_weight_branch(t_incl)
                    if wname_incl:
                        for arrs in t_incl.iterate([wname_incl], library="ak", step_size=UPROOT_STEP):
                            w_inclusive_tot += float(ak.sum(ak.values_astype(arrs[wname_incl], np.float64)))
                    else:
                        w_inclusive_tot += float(t_incl.num_entries)

                # 分子：test
                if INPUT_BASE_TREE_NAME not in f:
                    continue
                t = f[INPUT_BASE_TREE_NAME]
                has_mva = get_branch_exists(t, "MVA_Score")
                wname = pick_weight_branch(t)

                # 新增：確認 channel 分支存在與否
                has_z_mumu = get_branch_exists(t, "z_mumu")
                has_z_ee = get_branch_exists(t, "z_ee")

                if has_mva:
                    branches = ["MVA_Score"] + ([wname] if wname else [])
                    if has_z_mumu:
                        branches.append("z_mumu")
                    if has_z_ee:
                        branches.append("z_ee")

                    for arrs in t.iterate(branches, library="ak", step_size=UPROOT_STEP):
                        mva = arrs["MVA_Score"]
                        mask = mva >= mva_cut
                        if wname:
                            w = ak.values_astype(arrs[wname], np.float64)
                        else:
                            w = ak.values_astype(ak.ones_like(mva), np.float64)

                        n_pass += int(ak.sum(mask))
                        w_pass += float(ak.sum(w[mask]))

                        # 新增：channel-wise 計數
                        if has_z_mumu:
                            mu_mask = mask & (ak.values_astype(arrs["z_mumu"], np.int8) == 1)
                            n_pass_mu += int(ak.sum(mu_mask))
                            w_pass_mu += float(ak.sum(w[mu_mask]))
                        if has_z_ee:
                            ele_mask = mask & (ak.values_astype(arrs["z_ee"], np.int8) == 1)
                            n_pass_ele += int(ak.sum(ele_mask))
                            w_pass_ele += float(ak.sum(w[ele_mask]))

                        # 回退：若無 inclusive，用 test 總量當分母
                        if not incl_exists:
                            n_inclusive_tot += int(len(mva))
                            w_inclusive_tot += float(ak.sum(w))
                else:
                    # 沒有 MVA，仍需在缺 inclusive 時補上分母
                    if not incl_exists:
                        if wname:
                            total_w = 0.0
                            for arrs in t.iterate([wname], library="ak", step_size=UPROOT_STEP):
                                total_w += float(ak.sum(ak.values_astype(arrs[wname], np.float64)))
                            w_inclusive_tot += total_w
                            n_inclusive_tot += int(t.num_entries)
                        else:
                            n_inclusive_tot += int(t.num_entries)
                            w_inclusive_tot += float(t.num_entries)

        except Exception:
            continue
    return (
        n_inclusive_tot, w_inclusive_tot,
        n_pass, w_pass,
        n_pass_mu, w_pass_mu,
        n_pass_ele, w_pass_ele,
    )

def analyze_total_files(files: List[str]) -> Tuple[int, float]:
    # 不再使用 total 檔案讀取，保留函式但不呼叫
    n_tot = 0
    w_tot = 0.0
    return n_tot, w_tot

def safe_ratio(num: float, den: float) -> str:
    if den == 0:
        return "NA"
    return f"{num/den:.6g}"

def compute_eff(num: float, den: float, factor: float = 2.0) -> str:
    if den == 0:
        return "NA"
    return f"{factor * (num / den):.6g}"

def load_total_json(json_path: str) -> Dict[int, Dict[str, float]]:
    # 新增：從 TOTAL_JSON 取得每個 ma 的 n_events 與權重總和
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except Exception:
        return {}
    by_mass = data.get("by_mass", {})
    out: Dict[int, Dict[str, float]] = {}
    for k, v in by_mass.items():
        try:
            ma = int(k)
        except Exception:
            continue
        st = v.get("stats", {})
        n_events = int(st.get("n_events", 0))
        genw = st.get("genWeight", {}) or {}
        runs = st.get("runs", {}) or {}
        # 優先使用 genWeight.sum；缺少則回退 runs.genEventSumw_sum；再缺就用 n_events
        w_sum = float(genw.get("sum", runs.get("genEventSumw_sum", n_events)))
        out[ma] = {"n_events": float(n_events), "w_sum": float(w_sum)}
    return out

def main():
    os.makedirs(os.path.dirname(OUTPUT_BASE), exist_ok=True)

    # 1) 讀取 cut
    mva_cuts = parse_mva_cuts(optimized_BDT_Cut)

    # 2) 建立 ma 對應樣本
    ma_to_input, ma_to_total = build_ma_sample_maps(sig_samples, total_sig_samples)

    # 2.5) 從 TOTAL_JSON 讀取總事例數（pre-selection 前）
    total_json_stats = load_total_json(TOTAL_JSON)
    total_json_masses = set(total_json_stats.keys())

    # 3) 逐個 ma 計算
    # --- 取代原本的 lines/header 輸出成 JSON 結構 ---
    out = {
        "meta": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "input_base": INPUT_BASE,
            "total_json": TOTAL_JSON,
            "years": years_sig,
            "mva_branch": "MVA_Score",
            "weight_candidates": ALT_WEIGHT_CANDS,
            "optimized_BDT_Cut": optimized_BDT_Cut,
        },
        "by_mass": {}
    }

    candidate_mas = sorted(set(mva_cuts.keys()) & set(ma_to_input.keys()) & total_json_masses)
    for ma in candidate_mas:
        cut = mva_cuts[ma]
        input_sample = ma_to_input[ma]

        input_files = list_root_files(INPUT_BASE, input_sample, years_sig)

        n_input_tot, w_input_tot, n_pass_in, w_pass_in, n_pass_mu, w_pass_mu, n_pass_ele, w_pass_ele = analyze_input_files(input_files, cut)

        # 由 TOTAL_JSON 取得 pre-selection 總數（不再讀 EOS）
        jstats = total_json_stats.get(ma, {"n_events": 0.0, "w_sum": 0.0})
        n_total_tot = int(jstats["n_events"])
        w_total_tot = float(jstats["w_sum"])

        # 效率一律乘以 2（維持原邏輯，NA 以字串表示）
        ratio_input_unw = compute_eff(n_pass_in, n_input_tot)
        ratio_input_w   = compute_eff(w_pass_in, w_input_tot)
        ratio_total_unw = compute_eff(n_pass_in, n_total_tot)
        ratio_total_w   = compute_eff(w_pass_in, w_total_tot)

        # 新增：muon / ele channel 效率
        ratio_mu_input_unw = compute_eff(n_pass_mu, n_input_tot)
        ratio_mu_input_w   = compute_eff(w_pass_mu, w_input_tot)
        ratio_mu_total_unw = compute_eff(n_pass_mu, n_total_tot)
        ratio_mu_total_w   = compute_eff(w_pass_mu, w_total_tot)

        ratio_ele_input_unw = compute_eff(n_pass_ele, n_input_tot)
        ratio_ele_input_w   = compute_eff(w_pass_ele, w_input_tot)
        ratio_ele_total_unw = compute_eff(n_pass_ele, n_total_tot)
        ratio_ele_total_w   = compute_eff(w_pass_ele, w_total_tot)

        out["by_mass"][str(ma)] = {
            "cut": float(cut),
            "input": {
                "n_pass": int(n_pass_in),
                "w_pass": float(w_pass_in),
                "n_total": int(n_input_tot),
                "w_total": float(w_input_tot),
                "files": len(input_files),
                # 新增：channel-wise pass 數
                "n_pass_mu": int(n_pass_mu),
                "w_pass_mu": float(w_pass_mu),
                "n_pass_ele": int(n_pass_ele),
                "w_pass_ele": float(w_pass_ele),
            },
            "total": {
                "n_total": int(n_total_tot),
                "w_total": float(w_total_tot),
                "files": 0,
            },
            "efficiency": {
                "input_unweighted": ratio_input_unw,
                "input_weighted": ratio_input_w,
                "total_unweighted": ratio_total_unw,
                "total_weighted": ratio_total_w,
                # 新增：muon / ele 效率
                "muon": {
                    "input_unweighted": ratio_mu_input_unw,
                    "input_weighted": ratio_mu_input_w,
                    "total_unweighted": ratio_mu_total_unw,
                    "total_weighted": ratio_mu_total_w,
                },
                "ele": {
                    "input_unweighted": ratio_ele_input_unw,
                    "input_weighted": ratio_ele_input_w,
                    "total_unweighted": ratio_ele_total_unw,
                    "total_weighted": ratio_ele_total_w,
                },
            },
        }

    # --- 以 JSON 寫出 ---
    with open(OUTPUT_BASE, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

if __name__ == "__main__":
    main()