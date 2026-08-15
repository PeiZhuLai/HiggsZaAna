"""Closure test: offline python preprocessing must reproduce the online MLPhotons.

Reads a file produced by RecHitDumper/test/closure_cfg.py, which contains for
each event both the raw EB RecHits (rh_*) and the MLPhoton collection the C++
producer made from them (ml_*). This script re-clusters rh_* with
cluster_py.build_clusters(), runs the same two ONNX models, and compares.

This gate exists because everything downstream is built on the assumption that
the offline preprocessing IS the online preprocessing. If it does not hold,
a retrained model would be trained on images the detector never produces, and
nothing about the resulting mass scale would be trustworthy. So: run this first,
and do not proceed on a FAIL.

Input : closure ROOT file (TTree 'dumper/rechits')
Output: stdout report; optional --json <path> machine-readable summary

Usage (CERN, env `hzgdna` -- the only one with onnxruntime + uproot):
    /eos/home-p/pelai/App/Conda/.conda/envs/hzgdna/bin/python closure_test.py \\
        --input /eos/home-p/pelai/HZa/root_rechit/closure_M0p5.root \\
        --max-events 50
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_py import MLPhotonModels, build_clusters  # noqa: E402

DEFAULT_DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "RecoEgamma", "EgammaMLPhotonProducers", "data",
)

# float32 round-trips through ONNX and a python/C++ maths mismatch of order 1e-6
# are expected; anything looser than this means a real preprocessing difference.
RTOL = 1e-4
ATOL = 1e-6

COMPARE = [
    ("mass", "ml_mass"),
    ("moe", "ml_moe"),
    ("pt", "ml_pt"),
    ("eta", "ml_eta"),
    ("phi", "ml_phi"),
    ("diphotonScore", "ml_diphotonScore"),
    ("monophotonScore", "ml_monophotonScore"),
    ("hadronScore", "ml_hadronScore"),
    ("r1", "ml_r1"),
    ("r2", "ml_r2"),
    ("r3", "ml_r3"),
]


def rel_diff(a: float, b: float) -> float:
    denom = max(abs(a), abs(b), ATOL)
    return abs(a - b) / denom


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="closure ROOT file from closure_cfg.py")
    ap.add_argument("--tree", default="dumper/rechits")
    ap.add_argument("--classifier", default=os.path.join(DEFAULT_DATA_DIR, "classifier.onnx"))
    ap.add_argument("--regressor", default=os.path.join(DEFAULT_DATA_DIR, "regressor.onnx"))
    ap.add_argument("--max-events", type=int, default=20,
                    help="clustering is O(N^2) in crystals; start small")
    ap.add_argument("--match-dr", type=float, default=1e-3,
                    help="dR used to pair offline clusters with online MLPhotons "
                         "when the counts differ")
    ap.add_argument("--json", default=None, help="write a summary here")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    for p in (args.input, args.classifier, args.regressor):
        if not os.path.exists(p):
            print(f"[FATAL] missing input: {p}")
            return 2

    tree = uproot.open(args.input)[args.tree]
    branches = ["run", "lumi", "event", "rh_ieta", "rh_iphi", "rh_eta", "rh_phi",
                "rh_energy", "pv_z"] + [b for _, b in COMPARE]
    missing = [b for b in branches if b not in tree]
    if missing:
        print(f"[FATAL] tree {args.tree} lacks branches: {missing}")
        print("        -> the file was made with the pre-2026-08-09 dumper, or "
              "closure_cfg.py was run without mlPhotons/vertices. Re-dump.")
        return 2

    data = tree.arrays(branches, entry_stop=args.max_events, library="np")
    n_ev = len(data["event"])
    print(f"[closure] {args.input}")
    print(f"[closure] events: {n_ev}   models: {os.path.basename(args.classifier)}, "
          f"{os.path.basename(args.regressor)}")

    models = MLPhotonModels(args.classifier, args.regressor)

    n_cluster_match = 0
    n_cluster_mismatch = 0
    n_pairs = 0
    n_unmatched = 0
    worst = {name: (0.0, None) for name, _ in COMPARE}
    failures = []

    for iev in range(n_ev):
        rh_e = np.asarray(data["rh_energy"][iev], dtype=np.float64)
        n_hits = int(np.sum(rh_e > 0.0))
        pv_z = float(data["pv_z"][iev])
        if not np.isfinite(pv_z):
            failures.append((iev, "pv_z is NaN -- closure_cfg.py needs vertices=..."))
            continue

        clusters = build_clusters(data["rh_ieta"][iev], data["rh_iphi"][iev],
                                  data["rh_eta"][iev], data["rh_phi"][iev], rh_e)
        offline = models.predict(clusters, pv_z)
        online_n = len(np.asarray(data["ml_mass"][iev]))

        if len(offline) == online_n:
            n_cluster_match += 1
        else:
            n_cluster_mismatch += 1
            failures.append((iev, f"cluster count {len(offline)} offline vs {online_n} online "
                                  f"({n_hits} hits with E>0)"))

        # Online order follows the cluster loop, so index i should already line
        # up; fall back to a dR match so a single extra/missing cluster does not
        # smear every comparison downstream of it.
        on_eta = np.asarray(data["ml_eta"][iev], dtype=np.float64)
        on_phi = np.asarray(data["ml_phi"][iev], dtype=np.float64)

        for i, off in enumerate(offline):
            j = i if (i < online_n and len(offline) == online_n) else -1
            if j < 0 and online_n > 0:
                deta = on_eta - off["eta"]
                dphi = np.arctan2(np.sin(on_phi - off["phi"]), np.cos(on_phi - off["phi"]))
                dr = np.hypot(deta, dphi)
                k = int(np.argmin(dr))
                j = k if dr[k] < args.match_dr else -1
            if j < 0:
                n_unmatched += 1
                continue

            n_pairs += 1
            for name, br in COMPARE:
                a = float(off[name])
                b = float(np.asarray(data[br][iev])[j])
                d = rel_diff(a, b)
                if d > worst[name][0]:
                    worst[name] = (d, (iev, i, a, b))
                if d > RTOL and len(failures) < 200:
                    failures.append((iev, f"cluster {i}: {name} offline={a:.6g} "
                                         f"online={b:.6g} reldiff={d:.2e}"))

        if args.verbose:
            print(f"  ev {iev}: {n_hits} hits -> {len(offline)} clusters "
                  f"(online {online_n})")

    print()
    print(f"[closure] events with matching cluster count : {n_cluster_match}/{n_ev}")
    print(f"[closure] events with mismatching count      : {n_cluster_mismatch}")
    print(f"[closure] cluster pairs compared             : {n_pairs}")
    print(f"[closure] offline clusters with no online match: {n_unmatched}")
    print()
    print(f"{'variable':>18s} {'max reldiff':>12s}   worst case (ev, iclu, offline, online)")
    ok = True
    for name, _ in COMPARE:
        d, where = worst[name]
        flag = "" if d <= RTOL else "  <-- FAIL"
        if d > RTOL:
            ok = False
        loc = "" if where is None else f"ev{where[0]} c{where[1]} {where[2]:.6g} vs {where[3]:.6g}"
        print(f"{name:>18s} {d:12.3e}   {loc}{flag}")

    ok = ok and (n_cluster_mismatch == 0) and (n_unmatched == 0) and (n_pairs > 0)

    if failures:
        print(f"\n[closure] first {min(len(failures), 15)} of {len(failures)} issues:")
        for iev, msg in failures[:15]:
            print(f"   ev {iev}: {msg}")

    print()
    print(f"[closure] RESULT: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print("[closure] Do not train on this preprocessing until it passes.")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({
                "input": args.input,
                "n_events": n_ev,
                "n_cluster_count_match": n_cluster_match,
                "n_cluster_count_mismatch": n_cluster_mismatch,
                "n_pairs": n_pairs,
                "n_unmatched": n_unmatched,
                "max_reldiff": {k: worst[k][0] for k, _ in COMPARE},
                "pass": bool(ok),
            }, fh, indent=2)
        print(f"[closure] summary -> {args.json}")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
