"""Measure the merged-photon mass bias of whatever regressor is currently in use.

Reads the HiggsDNA merged ML parquet (which already carries the truth-matched
MLPhoton, `MergedML_*`, and the gen ALP) and reports, per generated mass point:

    * the cluster energy the regressor saw
    * the TRUE m/E  (= m_gen / E_cluster)   vs the PREDICTED m/E
    * m_rec / m_gen                          -- the mass scale bias
    * IQR/median of m_rec                    -- the resolution
    * the fraction of events whose true m/E falls below 0.005, the lower edge of
      the training domain quoted in AN2020-159 Table 7

The last column is the diagnostic that matters: the current model was trained on
a 75-2500 GeV particle gun with m/E in [0.005, 0.07], while HZa merged clusters
sit at E ~ 40 GeV. Once the true m/E drops under the training floor the network
saturates and the reconstructed mass is pushed up, badly.

Run this before AND after retraining -- it is the number that says whether the
retraining bought anything.

Input : <base>/mA_MLNANO_<tag>_<era>/merged_nominal.parquet
Output: stdout table (optional --csv)

Usage (CERN, env `hza_ana`):
    /eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python bias_scan.py
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pyarrow.parquet as pq

DEFAULT_BASE = ("/eos/cms/store/group/phys_susy/pelai/HZa_merged/"
                "parquet_merged_DNA_tmp/Sig_MC_MLNANO_all")
DEFAULT_TAGS = ["M0p1", "M0p2", "M0p3", "M0p4", "M0p5", "M0p6", "M0p7", "M0p8", "M0p9"]

# AN2020-159 Table 7: the regression NN was trained over this m/E range only.
TRAIN_MOE_LO = 0.005
TRAIN_MOE_HI = 0.07

NEEDED = ["MergedML_mass", "MergedML_massEnergyRatio", "MergedML_dR_to_genA",
          "GenALP_mass", "pass_allcuts_merged_ML"]


def firsts(col) -> np.ndarray:
    """First element of each (jagged) row, NaN when empty."""
    return np.array([v[0] if (v is not None and len(v) > 0) else np.nan
                     for v in col.to_pylist()], dtype=float)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=DEFAULT_BASE)
    ap.add_argument("--era", default="2024")
    ap.add_argument("--tags", nargs="+", default=DEFAULT_TAGS)
    ap.add_argument("--match-dr", type=float, default=0.05,
                    help="truth-match dR between the MLPhoton and the gen ALP")
    ap.add_argument("--csv", default=None)
    args = ap.parse_args()

    hdr = (f"{'tag':5s} {'m_gen':>6s} {'N':>6s} {'E_clu':>6s} | {'moe_true':>8s} "
           f"{'moe_pred':>8s} | {'m_rec':>7s} {'ratio':>6s} {'IQR/med':>8s} | "
           f"{'below':>6s} {'above':>6s}")
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for tag in args.tags:
        path = os.path.join(args.base, f"mA_MLNANO_{tag}_{args.era}", "merged_nominal.parquet")
        if not os.path.exists(path):
            print(f"{tag:5s} MISSING {path}")
            continue

        names = pq.read_schema(path).names
        missing = [c for c in NEEDED if c not in names]
        if missing:
            print(f"{tag:5s} lacks columns {missing}")
            continue

        t = pq.read_table(path, columns=NEEDED)
        mrec = t["MergedML_mass"].to_numpy()
        moep = t["MergedML_massEnergyRatio"].to_numpy()
        dr = t["MergedML_dR_to_genA"].to_numpy()
        sel = t["pass_allcuts_merged_ML"].to_numpy()
        mgen = firsts(t["GenALP_mass"])

        keep = (sel & (dr < args.match_dr) & (mrec > -100) & (moep > 1e-9)
                & np.isfinite(mgen) & (mgen > 0))
        if keep.sum() < 5:
            print(f"{tag:5s} too few matched events ({keep.sum()})")
            continue

        # E_cluster is what the producer divided by; recover it exactly rather
        # than rebuilding it from pt/eta (which carry the vertex correction).
        e_clu = mrec[keep] / moep[keep]
        mg, mr, mp = mgen[keep], mrec[keep], moep[keep]
        moe_true = mg / e_clu

        q16, q50, q84 = np.quantile(mr, [0.16, 0.5, 0.84])
        row = dict(
            tag=tag, m_gen=float(np.median(mg)), n=int(keep.sum()),
            e_cluster=float(np.median(e_clu)),
            moe_true=float(np.median(moe_true)), moe_pred=float(np.median(mp)),
            m_rec=float(q50), ratio=float(q50 / np.median(mg)),
            iqr_over_med=float((q84 - q16) / q50),
            frac_below_train=float(np.mean(moe_true < TRAIN_MOE_LO)),
            frac_above_train=float(np.mean(moe_true > TRAIN_MOE_HI)),
        )
        rows.append(row)
        print(f"{row['tag']:5s} {row['m_gen']:6.3f} {row['n']:6d} {row['e_cluster']:6.1f} | "
              f"{row['moe_true']:8.5f} {row['moe_pred']:8.5f} | {row['m_rec']:7.3f} "
              f"{row['ratio']:6.2f} {row['iqr_over_med']:8.2f} | "
              f"{row['frac_below_train']*100:5.1f}% {row['frac_above_train']*100:5.1f}%")

    print()
    print(f"below/above = fraction of events with true m/E outside the "
          f"[{TRAIN_MOE_LO}, {TRAIN_MOE_HI}] training domain (AN2020-159 Table 7).")

    if args.csv and rows:
        import csv
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"[bias_scan] -> {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
