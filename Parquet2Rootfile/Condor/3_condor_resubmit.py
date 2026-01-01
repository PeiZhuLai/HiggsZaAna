#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter an HTCondor queue-from joblist (TSV) by checking whether the output ROOT file (2nd column) exists,
rewrite the joblist, then (optionally) resubmit via a submit file (e.g. 2_submit.sub).

Designed for your current workflow:
  - 1_make_joblist.py -> joblist.tsv with columns:
      ARG_IN \t ARG_OUT \t ARG_CORR \t ARG_SPLIT
  - 2_submit.sub contains:
      queue ARG_IN, ARG_OUT, ARG_CORR, ARG_SPLIT from joblist.tsv

Default behavior:
  - remove jobs whose output (ARG_OUT) already exists
  - overwrite joblist.tsv (with a timestamped backup)
  - condor_submit 2_submit.sub

Run:
  python3 3_condor_resubmit.py --joblist joblist.tsv --submit-sub 2_submit.sub
Optional:
  --check-root --tree Events --min-entries 1   (treat existing but broken/empty ROOT as "missing")
  --dry-run                                 (do not write/submit)
  --no-submit                               (write joblist, but do not condor_submit)
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
from datetime import datetime
from typing import Iterable, Tuple, List


def _is_root_file_good(path: str, *, tree_names: Tuple[str, ...] = ("Events", "tree", "events"),
                       min_entries: int = 1) -> bool:
    """
    Return True if the ROOT file exists and looks usable.
    Only used when --check-root is set.
    """
    if not os.path.exists(path):
        return False

    # quick sanity
    try:
        if os.path.getsize(path) < 1024:  # 1KB: likely broken
            return False
    except OSError:
        return False

    try:
        from ROOT import TFile  # lazy import so ROOT is only needed when requested
    except Exception as e:
        raise RuntimeError(
            "ROOT is required for --check-root but cannot be imported. "
            "Run this on lxplus/CMSSW env where PyROOT is available."
        ) from e

    tf = TFile.Open(path)
    if not tf or tf.IsZombie():
        return False

    # Prefer checking a TTree's entries if it exists
    for tn in tree_names:
        obj = tf.Get(tn)
        if obj:
            try:
                if int(obj.GetEntries()) >= int(min_entries):
                    return True
            except Exception:
                # Not a TTree or cannot read entries; ignore and keep trying
                pass

    # Fallback: file has some keys at least
    try:
        return int(tf.GetNkeys()) > 0
    except Exception:
        return False


def filter_joblist_tsv(joblist_in: str,
                       joblist_out: str,
                       *,
                       check_root: bool = False,
                       tree_names: Tuple[str, ...] = ("Events", "tree", "events"),
                       min_entries: int = 1,
                       keep_comments: bool = True,
                       dry_run: bool = False) -> Tuple[int, int, int]:
    """
    Filter joblist by removing rows whose output (2nd column) already exists (and is good, if check_root=True).

    Returns: (n_total_rows, n_kept_rows, n_removed_rows)
    """
    kept_lines: List[str] = []
    n_total = 0
    n_removed = 0

    with open(joblist_in, "r", encoding="utf-8", errors="replace") as fin:
        for raw in fin:
            line = raw.rstrip("\n")
            if not line.strip():
                if keep_comments:
                    kept_lines.append(raw)
                continue

            if keep_comments and line.lstrip().startswith("#"):
                kept_lines.append(raw)
                continue

            cols = line.split("\t")
            if len(cols) < 2:
                # malformed line -> keep (safer)
                kept_lines.append(raw)
                continue

            n_total += 1
            out_path = cols[1].strip()

            exists = os.path.exists(out_path)
            if exists and check_root:
                # treat broken ROOT as missing => keep it
                try:
                    exists = _is_root_file_good(out_path, tree_names=tree_names, min_entries=min_entries)
                except RuntimeError as e:
                    # If ROOT missing, be explicit and fail early
                    raise

            if exists:
                n_removed += 1
            else:
                kept_lines.append(raw)

    if dry_run:
        return n_total, (n_total - n_removed), n_removed

    # write
    with open(joblist_out, "w", encoding="utf-8") as fout:
        for l in kept_lines:
            fout.write(l)

    return n_total, (n_total - n_removed), n_removed


def backup_file(path: str) -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.bak_{ts}"
    shutil.copy2(path, bak)
    return bak


def make_submit_with_joblist(sub_in: str, joblist_path: str) -> str:
    """
    Create a temporary submit file where the queue-from line points to joblist_path.
    Returns the new submit filename.
    """
    with open(sub_in, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    q_re = re.compile(r"^\s*queue\s+.+\s+from\s+(.+?)\s*$", re.IGNORECASE)
    replaced = False
    new_lines = []
    for l in lines:
        m = q_re.match(l)
        if m:
            prefix = l[: m.start(1)]
            new_lines.append(prefix + joblist_path + "\n")
            replaced = True
        else:
            new_lines.append(l)

    if not replaced:
        raise RuntimeError(f"Cannot find a 'queue ... from <file>' line in {sub_in}")

    sub_out = sub_in.replace(".sub", f".joblist_{os.path.basename(joblist_path)}.sub")
    with open(sub_out, "w", encoding="utf-8") as f:
        f.writelines(new_lines)
    return sub_out


def count_jobs(joblist_path: str) -> int:
    n = 0
    with open(joblist_path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 2:
                n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description="Filter queue-from joblist by checking output existence, then resubmit.")
    ap.add_argument("--joblist", default="joblist.tsv", help="Input joblist TSV (queue-from).")
    ap.add_argument("--joblist-out", default="", help="Output joblist TSV. If empty, overwrite --joblist.")
    ap.add_argument("--submit-sub", default="2_submit.sub", help="Submit file to condor_submit.")
    ap.add_argument("--no-submit", action="store_true", help="Only rewrite joblist; do not condor_submit.")
    ap.add_argument("--dry-run", action="store_true", help="Do not write/submit; just print summary.")
    ap.add_argument("--check-root", action="store_true",
                    help="If output exists but ROOT file is broken/empty, keep it for resubmit.")
    ap.add_argument("--tree", default="Events",
                    help="Preferred TTree name to check entries for --check-root. (Will also try 'tree'/'events'.)")
    ap.add_argument("--min-entries", type=int, default=1, help="Minimum entries for ROOT tree when --check-root is set.")
    args = ap.parse_args()

    joblist_in = os.path.abspath(args.joblist)
    joblist_out = os.path.abspath(args.joblist_out) if args.joblist_out else joblist_in

    tree_names = (args.tree, "tree", "events", "Events")

    if not os.path.exists(joblist_in):
        raise FileNotFoundError(f"joblist not found: {joblist_in}")

    if not args.dry_run and joblist_out == joblist_in:
        bak = backup_file(joblist_in)
        print(f"[backup] {bak}")

    n_total, n_kept, n_removed = filter_joblist_tsv(
        joblist_in,
        joblist_out,
        check_root=args.check_root,
        tree_names=tree_names,
        min_entries=args.min_entries,
        dry_run=args.dry_run,
    )

    print("====================================")
    print(f"joblist_in : {joblist_in}")
    print(f"joblist_out: {joblist_out}")
    print(f"rows total : {n_total}")
    print(f"rows kept  : {n_kept}")
    print(f"rows removed (output exists): {n_removed}")
    print("====================================")

    if args.dry_run:
        print("[dry-run] Not writing or submitting.")
        return

    njobs = count_jobs(joblist_out)
    if njobs == 0:
        print("[info] No jobs to submit after filtering. Stop.")
        return

    if args.no_submit:
        print(f"[info] Wrote filtered joblist with {njobs} jobs. (no-submit)")
        return

    sub_to_submit = args.submit_sub
    # If submit file still points to "joblist.tsv" but we wrote a different joblist_out,
    # create a temporary submit file that points to the correct file.
    if os.path.basename(joblist_out) != "joblist.tsv":
        sub_to_submit = make_submit_with_joblist(args.submit_sub, joblist_out)
        print(f"[info] submit file rewritten: {sub_to_submit}")

    cmd = ["condor_submit", sub_to_submit]
    print("[submit]", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
