#!/usr/bin/env python3
"""Build a NanoAODv15 -> MiniAODv6 UUID map for the DATA merged-photon friend
workflow, via the shared AOD grandparent.

Unlike private signal/bkg (NanoAOD's DAS parent IS the MiniAOD the MLNanoAOD was
produced from), central DATA reprocessing records the NanoAODv15's parent as the
*AOD* tier -- the MINIv6NANOv15 MiniAOD and NanoAOD are SIBLINGS, both children of
the same AOD. So `parent file=<nanoLFN>` returns AOD UUIDs that do NOT match the
MiniAOD UUIDs used to name the data MLNanoAOD files -> the friend join finds
nothing. (Confirmed: signal 33/33 overlap, data 0/6390.)

Fix = 2-hop via AOD:
  1. nano file  -> parent AOD LFNs            (parent file=<nanoLFN>)
  2. mini file  -> parent AOD LFNs            (parent file=<miniLFN>)  [invert -> AOD_UUID -> mini_UUID]
  3. nano -> [mini UUID] = union over the nano's AOD parents of AOD_UUID -> mini_UUID

Output JSON matches build_parent_map.py so mlphoton_loader consumes it unchanged:
    { "<nano LFN>": ["<mini UUID>", ...], ... }

Usage:
    python build_parent_map_data.py \\
        --nano-dataset /EGamma0/Run2024C-MINIv6NANOv15-v1/NANOAOD \\
        --out metadata/parent_maps/Data_2024_EGamma0_RunC_v1.json
    # --mini-dataset defaults to the nano dataset with /NANOAOD -> /MINIAOD
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

# Reuse the vetted dasgoclient helpers from the signal builder.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_parent_map import _check_dasgoclient, _das_query, _extract_uuid  # noqa: E402


def _parent_aods(lfn: str) -> tuple[str, list[str]]:
    """file -> its parent AOD UUIDs (as reported by DAS)."""
    return lfn, [_extract_uuid(p) for p in _das_query(f"parent file={lfn}")]


def _collect(files: list[str], workers: int, label: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_parent_aods, f): f for f in files}
        done = 0
        for fut in as_completed(futs):
            try:
                lfn, aods = fut.result()
                out[lfn] = aods
            except Exception as e:  # noqa: BLE001
                print(f"[{label}] ERROR {futs[fut]}: {e}", file=sys.stderr)
            done += 1
            if done % 500 == 0 or done == len(files):
                print(f"  [{label}] {done}/{len(files)}")
    return out


def build(nano_dataset: str, mini_dataset: str, workers: int) -> dict[str, list[str]]:
    nano_files = _das_query(f"file dataset={nano_dataset}")
    if not nano_files:
        raise SystemExit(f"no nano files in {nano_dataset!r}")

    # Detect the nano's file-level parent type. Some central reprocessings (e.g.
    # 2022/2023 NanoAODv15 made directly from the 22Sep2023 MiniAOD) record the
    # MiniAOD file as the nano's DIRECT parent, so `parent file=<nano>` returns
    # MINIAOD UUIDs we can use straight away (no AOD grandparent join). Others
    # (e.g. 2024 MINIv6NANOv15) record the AOD, needing the shared-AOD join below.
    sample_parents = _das_query(f"parent file={nano_files[0]}")
    if any("/MINIAOD" in p for p in sample_parents):
        print(f"[data parent_map] nano={len(nano_files)}; nano->MINI DIRECT "
              f"(file parent is MiniAOD)")
        nano_mini = _collect(nano_files, workers, "nano->MINI")
        mapping = {lfn: sorted(set(uuids)) for lfn, uuids in nano_mini.items()}
        n_empty = sum(1 for v in mapping.values() if not v)
        print(f"[data parent_map] {len(mapping)} nano files mapped; "
              f"{n_empty} with no mini match (direct)")
        return mapping

    mini_files = _das_query(f"file dataset={mini_dataset}")
    if not mini_files:
        raise SystemExit(f"no mini files in {mini_dataset!r}")
    print(f"[data parent_map] nano={len(nano_files)} mini={len(mini_files)}")

    # step 1: nano -> AOD uuids
    nano_aod = _collect(nano_files, workers, "nano->AOD")
    # step 2: mini -> AOD uuids, then invert to AOD_uuid -> [mini uuid]
    mini_aod = _collect(mini_files, workers, "mini->AOD")
    aod_to_mini: dict[str, list[str]] = {}
    for mini_lfn, aods in mini_aod.items():
        mu = _extract_uuid(mini_lfn)
        for a in aods:
            aod_to_mini.setdefault(a, []).append(mu)

    # step 3: nano -> mini uuids via shared AOD
    mapping: dict[str, list[str]] = {}
    n_empty = 0
    for nano_lfn, aods in nano_aod.items():
        minis: list[str] = []
        for a in aods:
            minis.extend(aod_to_mini.get(a, []))
        minis = sorted(set(minis))
        mapping[nano_lfn] = minis
        if not minis:
            n_empty += 1
    print(f"[data parent_map] {len(mapping)} nano files mapped; "
          f"{n_empty} with no mini match; "
          f"AOD keys with mini children={len(aod_to_mini)}")
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nano-dataset", required=True, help="DAS path of the NanoAODv15 dataset")
    ap.add_argument("--mini-dataset", default=None,
                    help="DAS path of the MiniAODv6 dataset (default: nano with /NANOAOD->/MINIAOD)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=24)
    args = ap.parse_args()
    _check_dasgoclient()
    mini = args.mini_dataset or args.nano_dataset.replace("/NANOAOD", "/MINIAOD")
    # Central DATA: the NanoAODv15 and its parent MiniAOD often carry different
    # processing strings (e.g. Run2022E-NanoAODv15-v1 vs Run2022E-22Sep2023-v1),
    # so the naive /NANOAOD->/MINIAOD substitution yields a non-existent dataset.
    # When the substituted mini has no files, resolve the real MiniAOD via the
    # DAS dataset-level parent of the nano.
    if not args.mini_dataset and not _das_query(f"file dataset={mini}"):
        parents = [p for p in _das_query(f"parent dataset={args.nano_dataset}")
                   if p.endswith("/MINIAOD")]
        if parents:
            mini = sorted(parents)[0]
            print(f"[data parent_map] resolved mini via DAS parent: {mini}")
    mapping = build(args.nano_dataset, mini, args.workers)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(mapping, f, indent=2, sort_keys=True)
    print(f"[data parent_map] wrote {len(mapping)} entries to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
