#!/usr/bin/env python3
"""
Build a NanoAODv15 -> parent MiniAOD UUID map from *lumi overlap*, for datasets
where DAS has no file-level parentage.

`build_parent_map.py` asks DAS `parent file=<nano_lfn>` once per file. For some
datasets that query returns nothing at all -- e.g.

    /DYGto2LG-1Jets_Bin-MLL-50_.../RunIII2024Summer24NanoAODv15-...-v2/NANOAODSIM

returns rc=0 with empty output for every one of its 1495 files, which silently
produces a map whose values are all `[]`. A map like that is not an error
anywhere downstream: the MLPhoton friend join simply finds no files to read and
the analysis writes a parquet with the MLPhoton columns missing.

The dataset-level parentage *is* there, so the parent MiniAOD can be found and
the file-level correspondence rebuilt from the lumisections both sides report:
a MiniAOD file is a parent of a NanoAOD file iff they share at least one lumi.
Two `file,lumi dataset=...` queries are enough -- no per-file DAS calls.

MC events all sit in run 1 here, so a bare lumi number is a unique key. Use
--verify against a dataset that *does* have DAS parentage before trusting the
output of one that does not.

Output format is identical to build_parent_map.py:

    {"/store/mc/.../<nano_uuid>.root": ["<mini_uuid>", ...], ...}

Usage:
    # rebuild a broken map
    python build_parent_map_lumi.py \
        --nano-dataset /DYGto2LG-1Jets_.../NANOAODSIM \
        --out metadata/parent_maps/Bkg_DYGto2LG_10to100_2024.json

    # closure: rebuild a dataset that has real DAS parentage and compare
    python build_parent_map_lumi.py \
        --nano-dataset /DYto2E-2Jets_.../NANOAODSIM \
        --verify metadata/parent_maps/Bkg_DYJetsTo2E_2024.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict

UUID_RE = re.compile(r"([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})")
LUMI_LIST_RE = re.compile(r"\[([0-9,\s]*)\]")


def _check_dasgoclient() -> None:
    if shutil.which("dasgoclient") is None:
        raise SystemExit(
            "dasgoclient not on PATH -- run inside a CMSSW env "
            "(`source /cvmfs/cms.cern.ch/cmsset_default.sh`)."
        )


def _das(query: str, timeout: int = 900) -> list[str]:
    r = subprocess.run(["dasgoclient", "-query", query],
                       check=True, capture_output=True, text=True, timeout=timeout)
    return [ln.strip() for ln in r.stdout.splitlines() if ln.strip()]


def _uuid(path: str) -> str:
    base = os.path.basename(path).rsplit(".root", 1)[0]
    m = UUID_RE.search(base)
    return m.group(1) if m else base


def file_lumis(dataset: str) -> dict[str, set[int]]:
    """lfn -> set of lumisections, from one `file,lumi dataset=` query."""
    out: dict[str, set[int]] = {}
    for line in _das(f"file,lumi dataset={dataset}"):
        lfn, _, rest = line.partition(" ")
        m = LUMI_LIST_RE.search(rest)
        if not m:
            continue
        body = m.group(1).strip()
        lumis = {int(x) for x in body.split(",") if x.strip()} if body else set()
        # a file can appear on more than one line (multiple blocks)
        out.setdefault(lfn, set()).update(lumis)
    return out


def build_via_lumi(nano_dataset: str) -> tuple[dict[str, list[str]], dict]:
    parents = _das(f"parent dataset={nano_dataset}")
    if not parents:
        raise SystemExit(f"no dataset-level parent for {nano_dataset!r} -- nothing to match against")
    if len(parents) > 1:
        print(f"[lumi_map] WARNING: {len(parents)} parent datasets, using all of them")

    print(f"[lumi_map] nano  : {nano_dataset}")
    nano = file_lumis(nano_dataset)
    print(f"[lumi_map]         {len(nano)} files, {sum(len(v) for v in nano.values())} lumi entries")

    lumi_to_mini: dict[int, set[str]] = defaultdict(set)
    n_mini = 0
    for p in parents:
        print(f"[lumi_map] parent: {p}")
        mini = file_lumis(p)
        n_mini += len(mini)
        print(f"[lumi_map]         {len(mini)} files, {sum(len(v) for v in mini.values())} lumi entries")
        for lfn, lumis in mini.items():
            u = _uuid(lfn)
            for l in lumis:
                lumi_to_mini[l].add(u)

    mapping: dict[str, list[str]] = {}
    n_empty = 0
    for lfn, lumis in nano.items():
        hits: set[str] = set()
        for l in lumis:
            hits |= lumi_to_mini.get(l, set())
        mapping[lfn] = sorted(hits)
        if not hits:
            n_empty += 1

    stats = {
        "nano_files": len(nano),
        "mini_files": n_mini,
        "empty": n_empty,
        "total_links": sum(len(v) for v in mapping.values()),
    }
    print(f"[lumi_map] built {stats['nano_files']} entries, "
          f"{stats['total_links']} links, {stats['empty']} empty")
    if n_empty:
        print(f"[lumi_map] WARNING: {n_empty} nano files matched no parent -- "
              f"the friend join will silently skip those events")
    return mapping, stats


def verify(mapping: dict[str, list[str]], ref_path: str) -> int:
    """Compare against a map built from real DAS file-level parentage."""
    with open(ref_path) as f:
        ref = json.load(f)
    common = set(mapping) & set(ref)
    print(f"\n[verify] reference: {ref_path}")
    print(f"[verify] keys: built={len(mapping)} ref={len(ref)} common={len(common)}")
    if not common:
        print("[verify] FAIL: no overlapping keys")
        return 1

    exact = subset = superset = other = 0
    for k in common:
        a, b = set(mapping[k]), set(ref[k])
        if a == b:
            exact += 1
        elif a < b:
            subset += 1
        elif a > b:
            superset += 1
        else:
            other += 1
    n = len(common)
    print(f"[verify] exact match : {exact}/{n} ({100.0*exact/n:.2f}%)")
    print(f"[verify] built  subset of ref (missing parents): {subset}")
    print(f"[verify] built superset of ref (extra parents) : {superset}")
    print(f"[verify] disjoint/partial                      : {other}")
    ok = exact == n
    print(f"[verify] {'PASS' if ok else 'MISMATCH -- inspect before trusting'}")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nano-dataset", required=True)
    ap.add_argument("--out", help="write the map here (omit with --verify to only check)")
    ap.add_argument("--verify", help="existing parent map to compare against")
    args = ap.parse_args()

    _check_dasgoclient()
    mapping, _ = build_via_lumi(args.nano_dataset)

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            json.dump(mapping, f, indent=2, sort_keys=True)
        print(f"[lumi_map] wrote {len(mapping)} entries to {args.out}")

    if args.verify:
        return verify(mapping, args.verify)
    return 0


if __name__ == "__main__":
    sys.exit(main())
