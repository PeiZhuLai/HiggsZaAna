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
import json
from datetime import datetime, timezone

# 你之前实际使用的组合：
years_sig  = ["2022preEE"]  # 信号

# 扫描的 ma 列表（背景也按你的原脚本扫）
ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]

# ===== Samples =====
TOTAL_BASE = "/eos/home-p/pelai/HZa/private_mc/signal/run3"
total_sig_samples = ["HZaTo2l2g_M5_NanoAODv12", "HZaTo2l2g_M15_NanoAODv12", "HZaTo2l2g_M30_NanoAODv12"]

OUTPUT_BASE = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sig_total_stats.json"

TOTAL_BASE_TREE_NAME = "Events"

def _parse_mass_from_name(name: str) -> Optional[int]:
    m = re.search(r"_M(\d+)_", name)
    return int(m.group(1)) if m else None

def _find_root_files(sample_dir: str, years_filter: List[str]) -> List[str]:
    files: List[str] = []
    if not os.path.isdir(sample_dir):
        return files
    for root, _, fnames in os.walk(sample_dir):
        for fn in fnames:
            if not fn.endswith(".root"):
                continue
            fp = os.path.join(root, fn)
            if years_filter and not any(y in fp for y in years_filter):
                continue
            files.append(fp)
    files.sort()
    return files

def _compute_stats_for_files(files: List[str], tree_name: str) -> Dict:
    total_events = 0
    sum_genw = 0.0
    sum_genw2 = 0.0
    sum_abs_genw = 0.0
    runs_genEventSumw = 0.0
    n_files_events_ok = 0
    n_files_total = len(files)
    bad_files: List[str] = []
    have_genw_any = False
    have_runs_sumw_any = False

    for f in files:
        try:
            with uproot.open(f) as rf:
                if tree_name not in rf:
                    bad_files.append(f"[no {tree_name}] " + f)
                    continue
                tree = rf[tree_name]
                total_events += int(tree.num_entries)
                n_files_events_ok += 1

                # genWeight 累計（若存在）
                tkeys = set(tree.keys())
                if "genWeight" in tkeys:
                    have_genw_any = True
                    for arrs in tree.iterate("genWeight", step_size="200 MB", library="np"):
                        gw = arrs["genWeight"]
                        if gw.size == 0:
                            continue
                        sum_genw += float(gw.sum(dtype=np.float64))
                        sum_genw2 += float((gw.astype(np.float64) ** 2).sum())
                        sum_abs_genw += float(np.abs(gw).sum(dtype=np.float64))

                # Runs 樹的 genEventSumw（若存在）
                if "Runs" in rf:
                    runs = rf["Runs"]
                    rkeys = set(runs.keys())
                    if "genEventSumw" in rkeys:
                        have_runs_sumw_any = True
                        vals = runs["genEventSumw"].array(library="np")
                        if vals.size:
                            runs_genEventSumw += float(vals.sum(dtype=np.float64))
        except Exception:
            bad_files.append(f"[error] " + f)
            continue

    out: Dict[str, object] = {
        "n_files_total": n_files_total,
        "n_files_ok": n_files_events_ok,
        "n_events": int(total_events),
        "genWeight": {
            "present": bool(have_genw_any),
            "sum": float(sum_genw) if have_genw_any else None,
            "sumw2": float(sum_genw2) if have_genw_any else None,
            "abs_sum": float(sum_abs_genw) if have_genw_any else None,
        },
        "runs": {
            "genEventSumw_present": bool(have_runs_sumw_any),
            "genEventSumw_sum": float(runs_genEventSumw) if have_runs_sumw_any else None,
        },
    }
    if bad_files:
        out["bad_files"] = bad_files
    return out

def _build_results() -> Dict:
    # 由樣本名解析質量，並逐一計算
    results_by_mass: Dict[str, Dict] = {}
    discovered_masses: List[int] = []
    years_filter: List[str] = []  # 不做年份過濾，包含所有 ROOT 檔
    for sample in total_sig_samples:
        mass = _parse_mass_from_name(sample)
        if mass is None:
            continue
        sample_dir = os.path.join(TOTAL_BASE, sample)
        files = _find_root_files(sample_dir, years_filter)  # 強制包含全部
        stats = _compute_stats_for_files(files, TOTAL_BASE_TREE_NAME)
        results_by_mass[str(mass)] = {
            "sample": sample,
            "stats": stats,
        }
        discovered_masses.append(mass)

    missing = [m for m in ma_list if str(m) not in results_by_mass]
    payload: Dict[str, object] = {
        "meta": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "total_base": TOTAL_BASE,
            "tree": TOTAL_BASE_TREE_NAME,
            "years_used": "ALL",  # 紀錄此次包含所有年份/子目錄
            "samples": total_sig_samples,
        },
        "by_mass": results_by_mass,
        "missing_masses": missing,
    }
    return payload

def main():
    payload = _build_results()
    outp = Path(OUTPUT_BASE)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("w") as f:
        json.dump(payload, f, indent=2)
    print(f"[OK] wrote stats to {outp} (masses: {', '.join(payload['by_mass'].keys())})")

if __name__ == "__main__":
    main()
