#!/usr/bin/env python3
"""Convert the sub-GeV merged-photon SIGNAL parquet into the flat ROOT trees the
merged dataVmc plotter reads.

Why a separate script.  Data and the backgrounds reach
`root_P2Root/run3_bdt_scored_mergedflag/` in two steps: `Parque2Root_za.py`
turns the *resolved* HiggsDNA parquet into `<sample>/<era>.root`, then
`Plot/scripts/add_merged_flag.py` joins the merged flags on (run, lumi, event)
from a friend parquet.  The signal cannot follow that route -- there is no
resolved sub-GeV production to join onto.  What there is instead is the merged
ML-NanoAOD production, and it already carries BOTH halves in one file: the
resolved `ALP_lead/sublead_photon_*` branches AND `pass_allcuts_merged_ML`.  So
one pass of the very same `decorate()` gives the plotter everything it needs.

Two things the resolved route supplies that have to be filled in here:

  has_merged_info   For data/bkg this separates "failed the merged selection"
                    from "was never in the friend parquet".  Every row here
                    comes out of the merged production itself, so it is 1
                    throughout -- there is no such thing as a missing row.
  weight            `decorate()` sets it from weight_central, same as the
                    resolved converter, so the plotter's weight branch resolves
                    identically for signal and background.

Only 2024 exists on the merged side (the ML-NanoAOD production is 2024-only),
which matches the data/bkg trees in the same directory.  The plotter's other
`years_sig` entries simply report MISS, exactly as they already do for the
backgrounds.

Trees written.  `Plot_Helper._run3_build_chain` reads SIGNAL out of the `test`
tree, not `inclusive` -- the resolved chain produces signal with
`Parque2Root_za.py --split`, keeps train+validation for the BDT and plots only
the held-out 30%.  The same `indices % 10` bucketing is reproduced here so the
sub-GeV points are split exactly like the resolved ones; writing `inclusive`
alone silently gives an empty chain (TChain reports the file as added and then
0 entries).

Input : <sig-dir>/mA_MLNANO_M<mass>_<era>/merged_nominal.parquet
Output: <out-dir>/mA_M<mass>/<era>.root   (trees `inclusive`, `train`,
        `validation`, `test`)

Usage (CERN, LCG_104 -- the same view `run_merged_datavmc.sh` sources):
    python Parque2Root_merged_signal.py
    python Parque2Root_merged_signal.py --masses 0p1,0p5 --dry-run
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Parque2Root_za import decorate, ensure_za_compatibility  # noqa: E402

SIG_DIR = ("/eos/cms/store/group/phys_susy/pelai/HZa_merged/"
           "parquet_merged_DNA_tmp/Sig_MC_MLNANO_all")
OUT_DIR = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_mergedflag"
MASSES = ["0p1", "0p2", "0p3", "0p4", "0p5", "0p6", "0p7", "0p8", "0p9"]
ERA = "2024"


def get_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sig-dir", default=SIG_DIR,
                        help="Merged signal parquet base directory.")
    parser.add_argument("--out-dir", default=OUT_DIR,
                        help="dataVmc input directory to write mA_M<mass>/<era>.root into.")
    parser.add_argument("--masses", default=",".join(MASSES),
                        help="Comma-separated flashgg mass tags, e.g. 0p1,0p5.")
    parser.add_argument("--era", default=ERA)
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would be written without writing it.")
    return parser.parse_args()


def writable_frame(data: pd.DataFrame) -> pd.DataFrame:
    """Drop what uproot cannot put in a flat TTree.

    `decorate()` coerces everything it can to numeric but deliberately leaves
    jagged/object columns alone rather than crashing on them, so they are still
    here and have to go before the write.
    """
    keep, dropped = [], []
    for column in data.columns:
        if pd.api.types.is_numeric_dtype(data[column]) or pd.api.types.is_bool_dtype(data[column]):
            keep.append(column)
        else:
            dropped.append(column)
    if dropped:
        print(f"    dropped {len(dropped)} non-numeric column(s): "
              f"{', '.join(dropped[:8])}{' ...' if len(dropped) > 8 else ''}")
    return data[keep]


def convert(mass: str, args) -> bool:
    src = os.path.join(args.sig_dir, f"mA_MLNANO_M{mass}_{args.era}", "merged_nominal.parquet")
    dst_dir = os.path.join(args.out_dir, f"mA_M{mass}")
    dst = os.path.join(dst_dir, f"{args.era}.root")

    print(f"\n=== m_a = {mass.replace('p', '.')} GeV")
    print(f"    in  : {src}")
    print(f"    out : {dst}")

    if not os.path.isfile(src):
        print("    [MISS] no merged_nominal.parquet -- skipped")
        return False

    data = pd.read_parquet(src)
    n_in = len(data)
    data = ensure_za_compatibility(data)
    data = decorate(data)

    # Every row is from the merged production, so it was always "processed".
    data["has_merged_info"] = np.uint8(1)
    if "pass_allcuts_merged_ML" not in data.columns:
        raise KeyError(f"{src} has no pass_allcuts_merged_ML -- wrong production?")
    data["pass_allcuts_merged_ML"] = data["pass_allcuts_merged_ML"].fillna(0).astype("uint8")

    data = writable_frame(data)

    passed = int(data["pass_allcuts_merged_ML"].sum())
    in_range = int(((data["pass_allcuts_merged_ML"] == 1)
                    & data["H_m"].between(95.0, 180.0)).sum())
    yield_ = float(data.loc[data["pass_allcuts_merged_ML"] == 1, "weight"].sum())
    print(f"    rows={n_in}  pass_merged_ML={passed}  in 95<m_llgg<180={in_range}  "
          f"sum(weight)={yield_:.3f}")

    if args.dry_run:
        print("    [dry-run] not written")
        return True

    # Same bucketing as Parque2Root_za.py main(): 0-4 train, 5-6 validation,
    # 7-9 test.  The plotter reads `test`, so this fixes the signal yield at 30%
    # of the sample -- matching how the resolved signal is already drawn.
    bucket = np.arange(len(data)) % 10
    splits = {
        "inclusive": data,
        "train": data[bucket < 5],
        "validation": data[(bucket >= 5) & (bucket < 7)],
        "test": data[bucket >= 7],
    }

    os.makedirs(dst_dir, exist_ok=True)
    with uproot.recreate(dst) as handle:
        for tree_name, frame in splits.items():
            handle[tree_name] = frame
    test_passed = int(splits["test"]["pass_allcuts_merged_ML"].sum())
    test_yield = float(splits["test"].loc[splits["test"]["pass_allcuts_merged_ML"] == 1,
                                          "weight"].sum())
    print(f"    wrote {len(data)} rows x {len(data.columns)} branches; "
          f"test tree {len(splits['test'])} rows, pass_merged_ML={test_passed}, "
          f"sum(weight)={test_yield:.3f}")
    return True


def main() -> None:
    args = get_args()
    masses = [m.strip() for m in args.masses.split(",") if m.strip()]
    print(f"[Config] masses={masses} era={args.era}")
    print(f"[Config] sig-dir={args.sig_dir}")
    print(f"[Config] out-dir={args.out_dir}")

    ok = [m for m in masses if convert(m, args)]
    print(f"\nDone: {len(ok)}/{len(masses)} mass point(s) converted.")
    if len(ok) != len(masses):
        missing = [m for m in masses if m not in ok]
        print(f"Missing input for: {missing}")
        sys.exit(1)


if __name__ == "__main__":
    main()
