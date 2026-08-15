"""Re-derive the per-mA ROI windows from a signal parquet production.

merged_p2root.py opens a window on MLPhoton_lead_mass for each m_a, defined as
the [q16, q84] of the signal reco distribution (the central 68%). Those numbers
are a property of the REGRESSOR: change the model and every window moves. The
committed ROI_WINDOWS were derived from the old model on 2026-07-01 and are not
valid for a re-trained one.

Prints a ready-to-paste ROI_WINDOWS dict, and -- when --compare is given -- the
old production's windows next to it so the change is visible.

A narrower window at the same containment means better mass resolution, which is
the whole point of the retraining; a window that moves without narrowing means
the mass scale shifted and the signal model has to be refit anyway.

Input : <base>/mA_MLNANO_<tag>_<era>/merged_nominal.parquet
Output: stdout

Usage (CERN, env `hza_ana`):
    python derive_roi_windows.py \\
        --base /eos/.../parquet_merged_DNA_v3/Sig_MC_MLNANO_all \\
        --compare /eos/.../parquet_merged_DNA_tmp/Sig_MC_MLNANO_all
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pyarrow.parquet as pq

DEFAULT_BASE = ("/eos/cms/store/group/phys_susy/pelai/HZa_merged/"
                "parquet_merged_DNA_v3/Sig_MC_MLNANO_all")
OLD_BASE = ("/eos/cms/store/group/phys_susy/pelai/HZa_merged/"
            "parquet_merged_DNA_tmp/Sig_MC_MLNANO_all")
# merged analysis is sub-GeV only (1_grand_merged.sh MASSES)
TAGS = ["M0p1", "M0p2", "M0p3", "M0p4", "M0p5", "M0p6", "M0p7", "M0p8", "M0p9"]
SEL = "pass_allcuts_merged_ML"
ROI_VAR = "MLPhoton_lead_mass"


def windows(base, era, qlo, qhi):
    out = {}
    for tag in TAGS:
        path = os.path.join(base, f"mA_MLNANO_{tag}_{era}", "merged_nominal.parquet")
        if not os.path.exists(path):
            out[tag] = None
            continue
        names = pq.read_schema(path).names
        if SEL not in names or ROI_VAR not in names:
            out[tag] = None
            continue
        t = pq.read_table(path, columns=[SEL, ROI_VAR])
        sel = t[SEL].to_numpy()
        v = t[ROI_VAR].to_numpy()
        m = sel & np.isfinite(v) & (v > -100)
        if m.sum() < 20:
            out[tag] = None
            continue
        lo, hi = np.quantile(v[m], [qlo, qhi])
        out[tag] = (float(lo), float(hi), int(m.sum()), float(np.median(v[m])))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=DEFAULT_BASE)
    ap.add_argument("--compare", default=None,
                    help=f"a second production to show alongside (e.g. {OLD_BASE})")
    ap.add_argument("--era", default="2024")
    ap.add_argument("--qlo", type=float, default=0.16)
    ap.add_argument("--qhi", type=float, default=0.84)
    args = ap.parse_args()

    new = windows(args.base, args.era, args.qlo, args.qhi)
    old = windows(args.compare, args.era, args.qlo, args.qhi) if args.compare else {}

    hdr = (f"{'tag':>6s} {'N':>7s} {'median':>8s} {'q16':>8s} {'q84':>8s} {'width':>8s}"
           + (f" | {'old width':>9s} {'ratio':>6s}" if old else ""))
    print(hdr)
    print("-" * len(hdr))
    for tag in TAGS:
        n = new.get(tag)
        if n is None:
            print(f"{tag:>6s}  MISSING")
            continue
        lo, hi, cnt, med = n
        w = hi - lo
        line = (f"{tag:>6s} {cnt:>7d} {med:>8.3f} {lo:>8.3f} {hi:>8.3f} {w:>8.3f}")
        o = old.get(tag)
        if old:
            if o is None:
                line += f" | {'--':>9s} {'--':>6s}"
            else:
                ow = o[1] - o[0]
                line += f" | {ow:>9.3f} {w/ow:>6.2f}"
        print(line)

    print()
    print("Paste into merged_p2root.py (ROI_WINDOWS):")
    print("ROI_WINDOWS = {")
    for i, tag in enumerate(TAGS):
        n = new.get(tag)
        if n is None:
            continue
        end = "," if i < len(TAGS) - 1 else ","
        print(f'    "{tag}": ({n[0]:.3f}, {n[1]:.3f}){end}')
    print("}")
    if old:
        print()
        print("ratio < 1 means the window narrowed at the same 68% containment, "
              "i.e. better mass resolution.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
