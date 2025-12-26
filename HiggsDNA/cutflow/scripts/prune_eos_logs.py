import os
import re
import argparse
import shutil
from collections import defaultdict

dataType = ["Bkg_MC","Sig_MC","Data"]

# Year of Signal
year_sig_2022 = ["2022preEE", "2022postEE"]
# year_sig_2023 = ["2023preBPix", "2023postBPix"]
# year_sig_2024 = ["2024"]
# Year of Bkg
# year_DYG_2022 = ["2022preEE", "2022postEE"]
# year_DYG_2023 = ["2023preBPix", "2023postBPix"]
# year_DYG_2024 = ["2024"]
# years_DYJet_2022 = ["2022preEE","2022postEE"]
# year_DYJet_2023  = ["2023preBPix", "2023postBPix"]
# years_DYJet_2024 = ["2024"]
# # Year of Data
# years_Data_2022 = ["2022preEE","2022postEE"]
# year_Data_2023  = ["2023preBPix", "2023postBPix"]
# years_Data_2024 = ["2024"]

# Name of Signal Sample
name_sig_2022 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
# name_sig_2023 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
# name_sig_2024 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
# Name of Bkg Sample
# name_DYG_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
# name_DYG_2023 = ["DYGto2LG_10to100"]
# name_DYG_2024 = ["DYGto2LG_10to100"]
# name_DYJet_2022 = ["DYJetsToLL"]
# name_DYJet_2023 = ["DYJetsToLL"]
# name_DYJet_2024 = ["DYJetsTo2E","DYJetsTo2Mu","DYJetsTo2Tau"]
# # Name of Data Sample
# name_Data_2022 = ["Data"]
# name_Data_2023 = ["Data"]
# name_Data_2024 = ["Data"]

# 当前目录（根目錄）
path = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs"

def iter_years_and_names(dtype: str):
    # 用你檔案上方定義的清單，統一在這裡映射
    if dtype == "Sig_MC":
        return [
            (year_sig_2022, name_sig_2022),
            (year_sig_2023, name_sig_2023),
            (year_sig_2024, name_sig_2024),
        ]
    if dtype == "Bkg_MC":
        return [
            (year_DYG_2022, name_DYG_2022),
            (year_DYG_2023, name_DYG_2023),
            (year_DYG_2024, name_DYG_2024),
            (years_DYJet_2022, name_DYJet_2022),
            (year_DYJet_2023, name_DYJet_2023),
            (years_DYJet_2024, name_DYJet_2024),
        ]
    if dtype == "Data":
        return [
            (years_Data_2022, name_Data_2022),
            (year_Data_2023, name_Data_2023),
            (years_Data_2024, name_Data_2024),
        ]
    return []

job_dir_pat = re.compile(r"^job_(\d+)$")
keep_ext = {".out", ".err", ".log"}  # files we consider pruning within each job dir

# NEW: match a trailing numeric "version" before extension.
# Supports: 10869318.0.log, 10869318.0.err, 10869318.0.out, etc.
file_ver_pat = re.compile(r"^(?P<ver>\d+(?:\.\d+)?)(?P<ext>\.out|\.err|\.log)$")

def prune_one_job_dir(job_dir: str, dry_run: bool, verbose: bool):
    if not os.path.isdir(job_dir):
        return

    ver_to_files = defaultdict(list)

    for fn in os.listdir(job_dir):
        m = file_ver_pat.match(fn)
        if not m:
            continue
        ext = m.group("ext")
        if ext not in keep_ext:
            continue

        # Compare versions numerically; accepts "10869318.0"
        ver = float(m.group("ver"))
        ver_to_files[ver].append(os.path.join(job_dir, fn))

    vers = sorted(ver_to_files.keys())
    if len(vers) <= 1:
        if verbose:
            print(f"[SKIP] {job_dir}: only {len(vers)} version group(s)")
        return

    keep_ver = max(vers)
    del_vers = [v for v in vers if v != keep_ver]
    print(f"[PRUNE] {job_dir} -> keep version: {keep_ver}; delete versions: {del_vers}")

    for v in del_vers:
        for fpath in sorted(ver_to_files[v]):
            if dry_run:
                print("  DRY-RUN delete:", fpath)
            else:
                try:
                    os.remove(fpath)
                    print("  deleted:", fpath)
                except FileNotFoundError:
                    print("  missing (skip):", fpath)

def prune_one_sample(sample_dir: str, dry_run: bool, verbose: bool):
    if not os.path.isdir(sample_dir):
        return

    # Now: for each job_<id> directory, prune log/err/out inside it by largest trailing number.
    for entry in os.listdir(sample_dir):
        if not job_dir_pat.match(entry):
            continue
        job_dir = os.path.join(sample_dir, entry)
        if not os.path.isdir(job_dir):
            continue
        prune_one_job_dir(job_dir, dry_run=dry_run, verbose=verbose)

def main():
    ap = argparse.ArgumentParser(
        description="Prune eos_logs: keep only the largest job_<id> group (.out/.err/.log) for each sample; delete all smaller job_<id> groups."
    )
    ap.add_argument("--base", default=path, help="Base eos_logs dir")
    ap.add_argument("--do-delete", action="store_true", help="Actually delete files (default: dry-run)")
    ap.add_argument("--verbose", action="store_true", help="More prints")
    args = ap.parse_args()

    base = args.base
    dry_run = (not args.do_delete)

    for dtype in dataType:
        dtype_dir = os.path.join(base, dtype)
        if not os.path.isdir(dtype_dir):
            if args.verbose:
                print(f"[SKIP] missing dtype dir: {dtype_dir}")
            continue

        for years_list, names_list in iter_years_and_names(dtype):
            for year in years_list:
                for name in names_list:
                    sample_dir = os.path.join(dtype_dir, f"{name}_{year}")
                    prune_one_sample(sample_dir, dry_run=dry_run, verbose=args.verbose)

if __name__ == "__main__":
    main()

# Usage:

# Dry-run (no deletion)
# python prune_eos_logs.py

# Action Delete
# python prune_eos_logs.py --do-delete