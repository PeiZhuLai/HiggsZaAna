#!/usr/bin/env python3
"""
Build a NanoAODv15 → parent MiniAOD UUID mapping for the merged-photon
friend-tree workflow.

Central NanoAODv15 production is *not* 1:1 with MiniAODv6 — a single
NanoAODv15 file may aggregate events from multiple MiniAOD parents, and a
MiniAOD's events may end up in multiple NanoAOD children. Because our
MLPhoton MLNanoAOD output is keyed by MiniAOD UUID, we need a per-NanoAOD-file
list of parent MiniAOD UUIDs in order to know which MLNanoAOD files to read
when injecting MLPhoton data into the analysis.

This script runs `dasgoclient parent file=<nano_lfn>` once per NanoAODv15 file
in the given dataset and writes a JSON file:

    {
        "/store/mc/.../<nano_uuid>.root": ["<mini_uuid_1>", "<mini_uuid_2>", ...],
        ...
    }

Usage:
    python build_parent_map.py --nano-dataset /Foo.../NANOAODSIM \\
                               --out /eos/.../parent_map_<tag>.json
    python build_parent_map.py --batch <yaml_or_json with tag→nano_dataset>
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

UUID_RE = re.compile(r"([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})")


def _check_dasgoclient() -> None:
    if shutil.which("dasgoclient") is None:
        raise SystemExit(
            "dasgoclient not on PATH — run inside a CMSSW env "
            "(`source /cvmfs/cms.cern.ch/cmsset_default.sh && cmsenv`)."
        )


def _das_query(query: str, timeout: int = 60, retries: int = 3) -> list[str]:
    """Run dasgoclient -query with retries on transient failures."""
    last_err = None
    for attempt in range(retries):
        try:
            r = subprocess.run(
                ["dasgoclient", "-query", query],
                check=True, capture_output=True, text=True, timeout=timeout,
            )
            return [ln.strip() for ln in r.stdout.splitlines() if ln.strip()]
        except subprocess.CalledProcessError as e:
            last_err = f"exit {e.returncode}: {e.stderr.strip()[:200]}"
        except subprocess.TimeoutExpired:
            last_err = "timeout"
        except Exception as e:
            last_err = repr(e)
        if attempt + 1 < retries:
            import time
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"dasgoclient failed after {retries} retries: {query!r} → {last_err}")


def _extract_uuid(path: str) -> str:
    """Return the UUID-like basename of a CMSSW file path, without .root."""
    base = os.path.basename(path).rsplit(".root", 1)[0]
    m = UUID_RE.search(base)
    return m.group(1) if m else base


def parent_uuids_for_one(nano_lfn: str) -> tuple[str, list[str]]:
    parents = _das_query(f"parent file={nano_lfn}")
    return nano_lfn, [_extract_uuid(p) for p in parents]


def build_map(nano_dataset: str, workers: int = 16) -> dict[str, list[str]]:
    nano_files = _das_query(f"file dataset={nano_dataset}")
    if not nano_files:
        raise SystemExit(f"no files in dataset {nano_dataset!r}")
    print(f"[parent_map] {nano_dataset}: {len(nano_files)} nano files")
    mapping: dict[str, list[str]] = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(parent_uuids_for_one, f): f for f in nano_files}
        done = 0
        for fut in as_completed(futures):
            try:
                lfn, parents = fut.result()
            except Exception as e:
                print(f"[parent_map] ERROR on {futures[fut]}: {e}", file=sys.stderr)
                continue
            mapping[lfn] = parents
            done += 1
            if done % 200 == 0 or done == len(nano_files):
                print(f"  progress: {done}/{len(nano_files)} files")
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nano-dataset", help="DAS path of the NanoAODv15 dataset")
    ap.add_argument("--out", required=True, help="output JSON path")
    ap.add_argument("--workers", type=int, default=16,
                    help="parallel dasgoclient calls (default 16)")
    args = ap.parse_args()

    _check_dasgoclient()
    if not args.nano_dataset:
        raise SystemExit("--nano-dataset is required")

    mapping = build_map(args.nano_dataset, workers=args.workers)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(mapping, f, indent=2, sort_keys=True)
    print(f"[parent_map] wrote {len(mapping)} entries to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
