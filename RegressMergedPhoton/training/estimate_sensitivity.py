"""Signal-only sensitivity estimate for the retrained MLPhoton models.

WHAT THIS CAN AND CANNOT TELL YOU
---------------------------------
It compares the SIGNAL side only: yield through the selection, and the width of
the per-mA ROI window. The background side is unavailable -- both the data
friend parquet and the background MC were produced with the OLD model, and the
retrained regressor moves the whole MLPhoton_lead_mass scale, so those files say
nothing about what the background looks like under the new model.

So S/sqrt(B) cannot be computed. What is computed instead is the signal yield
ratio, and S/sqrt(B) under two EXPLICIT and opposite assumptions about B:

  scenario A  "background scales with the window"
      B_new/B_old = width_new/width_old
      i.e. the background density per unit mass is unchanged and only the
      narrower window matters. OPTIMISTIC: it ignores that the new classifier
      admits more clusters overall, background included.

  scenario B  "background scales with the selection, like the signal"
      B_new/B_old = N_new/N_old  (times the window ratio)
      i.e. the extra clusters the new classifier lets through are background in
      the same proportion as signal. PESSIMISTIC for the classifier: the ROC
      says it is more efficient at fixed contamination, so the true background
      growth should be smaller than the signal growth.

The honest answer lives between them. The ROC measurement (compare_classifiers.py:
at a fixed hadronic fake rate the new model keeps 1.5-1.7x more diphotons) is the
one number that does constrain the classifier's real gain -- quoted here for
reference but NOT folded in, because it is a cluster-level statement and the
analysis selection is not a pure fake-rate cut.

Input : old and new signal parquet productions + their ROI windows
Output: stdout table

Usage (CERN, env `hza_ana`):
    python estimate_sensitivity.py
"""

from __future__ import annotations

import argparse
import math
import os

import numpy as np
import pyarrow.parquet as pq

EOSB = "/eos/cms/store/group/phys_susy/pelai/HZa_merged"
OLD_BASE = f"{EOSB}/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all"
NEW_BASE = f"{EOSB}/parquet_merged_DNA_v4/Sig_MC_MLNANO_all"
TAGS = ["M0p1", "M0p2", "M0p3", "M0p4", "M0p5", "M0p6", "M0p7", "M0p8", "M0p9"]

# committed windows (old model) and the v4 re-derivation
ROI_OLD = {
    "M0p1": (0.154, 0.947), "M0p2": (0.176, 0.887), "M0p3": (0.226, 0.917),
    "M0p4": (0.326, 0.796), "M0p5": (0.413, 0.821), "M0p6": (0.489, 0.878),
    "M0p7": (0.583, 0.925), "M0p8": (0.666, 0.996), "M0p9": (0.747, 1.079),
}
ROI_NEW = {
    "M0p1": (0.132, 0.352), "M0p2": (0.168, 0.342), "M0p3": (0.222, 0.409),
    "M0p4": (0.263, 0.493), "M0p5": (0.287, 0.581), "M0p6": (0.304, 0.672),
    "M0p7": (0.311, 0.761), "M0p8": (0.325, 0.852), "M0p9": (0.330, 0.938),
}
SEL = "pass_allcuts_merged_ML"
ROI_VAR = "MLPhoton_lead_mass"
WEIGHT = "weight_central"


def yields(base, tag, window, era="2024"):
    path = os.path.join(base, f"mA_MLNANO_{tag}_{era}", "merged_nominal.parquet")
    if not os.path.exists(path):
        return None
    names = pq.read_schema(path).names
    cols = [c for c in (SEL, ROI_VAR, WEIGHT) if c in names]
    if SEL not in cols or ROI_VAR not in cols:
        return None
    t = pq.read_table(path, columns=cols)
    sel = t[SEL].to_numpy()
    v = t[ROI_VAR].to_numpy()
    w = t[WEIGHT].to_numpy() if WEIGHT in cols else np.ones(len(v))
    good = sel & np.isfinite(v) & (v > -100)
    inroi = good & (v >= window[0]) & (v < window[1])
    return {
        "n_sel": int(good.sum()),
        "n_roi": int(inroi.sum()),
        "s_sel": float(np.sum(w[good])),
        "s_roi": float(np.sum(w[inroi])),
        "width": window[1] - window[0],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", default=OLD_BASE)
    ap.add_argument("--new", default=NEW_BASE)
    ap.add_argument("--era", default="2024")
    args = ap.parse_args()

    hdr = (f"{'tag':>6s} | {'S_old':>10s} {'S_new':>10s} {'S ratio':>8s} | "
           f"{'w_old':>6s} {'w_new':>6s} {'w ratio':>8s} | "
           f"{'A: S/sqrtB':>11s} {'B: S/sqrtB':>11s}")
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for tag in TAGS:
        o = yields(args.old, tag, ROI_OLD[tag], args.era)
        n = yields(args.new, tag, ROI_NEW[tag], args.era)
        if o is None or n is None:
            print(f"{tag:>6s} | MISSING ({'old' if o is None else 'new'})")
            continue
        s_ratio = n["s_roi"] / o["s_roi"] if o["s_roi"] > 0 else float("nan")
        w_ratio = n["width"] / o["width"]
        # A: B scales with the window only
        a = s_ratio / math.sqrt(w_ratio) if w_ratio > 0 else float("nan")
        # B: B scales with selection growth AND the window
        n_ratio = n["n_sel"] / o["n_sel"] if o["n_sel"] > 0 else float("nan")
        b = s_ratio / math.sqrt(w_ratio * n_ratio) if w_ratio > 0 and n_ratio > 0 else float("nan")
        rows.append((tag, s_ratio, w_ratio, a, b))
        print(f"{tag:>6s} | {o['s_roi']:>10.3g} {n['s_roi']:>10.3g} {s_ratio:>8.1f} | "
              f"{o['width']:>6.3f} {n['width']:>6.3f} {w_ratio:>8.2f} | "
              f"{a:>11.2f} {b:>11.2f}")

    if rows:
        print()
        print("Scenario A (background density unchanged, only the window shrinks)"
              " -- OPTIMISTIC")
        print("Scenario B (background grows exactly like the selection does)"
              " -- PESSIMISTIC")
        print()
        print("Neither is right. A ignores that a looser classifier admits more")
        print("background; B ignores that the ROC shows the new classifier is")
        print("MORE efficient at fixed contamination (1.5-1.7x). The true number")
        print("needs the data friend parquet re-produced with the new models.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
