#!/usr/bin/env python3
"""
Reconcile a merged-ML friend tag: chunk rows vs merged rows, merge/chunk mtimes,
and the MLPhoton columns that are the whole point of the reproduction.

Three ways this output can be quietly wrong, all guarded here:
  * counting merged_nominal.parquet among the chunks roughly doubles chunk_rows
    and makes a short merge look complete -- it is excluded by name;
  * a merged file whose mtime predates its newest chunk merged before that chunk
    was rewritten, so it is missing those events even when nothing errored;
  * a friend with no pass_allcuts_merged_ML column reads downstream as "nothing
    passed the merged selection" rather than "the join never ran".

⚠️ Run this with the `hza_ana` pyarrow (21.x), NOT LCG_104's:

    export PATH=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin:$PATH
    python reconcile_bkg_ml_friend.py

LCG_104's pyarrow cannot decode the encoding HiggsDNA writes and raises
ArrowNotImplementedError on perfectly good chunks -- on 2026-08-30 that read as
"all 460 chunks corrupt" when every one of them was fine. A genuinely damaged
file raises OSError ("Deserializing page header failed"), not
ArrowNotImplementedError; the check below keeps them apart so a stale reader is
never mistaken for data loss.

Usage:
    python reconcile_bkg_ml_friend.py [tag ...]      # default: all four 2024 backgrounds
    python reconcile_bkg_ml_friend.py --base <dir> Bkg_DYJetsTo2E_2024
"""
from __future__ import annotations

import argparse
import datetime as dt
import glob
import os
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

if tuple(int(x) for x in pa.__version__.split(".")[:2]) < (16, 0):
    sys.exit("pyarrow %s is too old to decode these files and will report good "
             "chunks as corrupt. Use hza_ana: export PATH=/eos/home-p/pelai/App/"
             "Conda/.conda/envs/hza_ana/bin:$PATH" % pa.__version__)

BASE = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend_ML"
TAGS = ["Bkg_DYGto2LG_10to100_2024", "Bkg_DYJetsTo2E_2024",
        "Bkg_DYJetsTo2Mu_2024", "Bkg_DYJetsTo2Tau_2024"]
NEEDED = ["run", "luminosityBlock", "event", "pass_allcuts_merged_ML",
          "MLPhoton_lead_mass", "MLPhoton_lead_diphotonScore", "n_MLPhoton_diphoton"]


def check(tag: str, base: str) -> bool:
    d = os.path.join(base, tag)
    merged = glob.glob(os.path.join(d, "*", "merged_nominal.parquet"))
    print("=" * 68)
    print(tag)
    if not merged:
        print("  no merged_nominal.parquet -- not merged yet")
        return False
    mg = merged[0]

    # find(), not ls glob: these dirs can hold enough files that a shell glob
    # overflows and silently reports zero.
    chunks = [f for f in glob.glob(os.path.join(d, "**", "*.parquet"), recursive=True)
              if os.path.basename(f) != "merged_nominal.parquet"]
    if not chunks:
        print("  no chunk parquet found -- cannot reconcile")
        return False

    crow = sum(pq.ParquetFile(f).metadata.num_rows for f in chunks)
    mrow = pq.ParquetFile(mg).metadata.num_rows

    # metadata alone cannot tell a healthy file from a corrupt one
    try:
        pq.read_table(mg, columns=["event"])
    except Exception as e:
        print("  READABLE     : FAIL -- %s: %s" % (type(e).__name__, str(e)[:90]))
        print("  => PROBLEM (merged file unreadable; check the chunks before re-merging)")
        return False
    rec_ok = crow == mrow
    print("  chunks       : %d files, %d rows" % (len(chunks), crow))
    print("  merged       : %d rows" % mrow)
    print("  RECONCILE    : %s" % ("PASS" if rec_ok else "FAIL  diff=%d" % (crow - mrow)))

    mt_m = os.path.getmtime(mg)
    mt_c = max(os.path.getmtime(f) for f in chunks)
    mt_ok = mt_m >= mt_c
    print("  merged mtime : %s" % dt.datetime.fromtimestamp(mt_m).strftime("%Y-%m-%d %H:%M:%S"))
    print("  newest chunk : %s" % dt.datetime.fromtimestamp(mt_c).strftime("%Y-%m-%d %H:%M:%S"))
    print("  MTIME        : %s" % ("PASS" if mt_ok else "FAIL -- merged predates a chunk"))

    names = pq.read_schema(mg).names
    miss = [c for c in NEEDED if c not in names]
    col_ok = not miss
    print("  ncols        : %d" % len(names))
    print("  COLUMNS      : %s" % ("PASS" if col_ok else "FAIL missing %s" % miss))

    pass_ok = False
    if "pass_allcuts_merged_ML" in names:
        t = pq.read_table(mg, columns=["pass_allcuts_merged_ML", "MLPhoton_lead_mass"])
        f = np.nan_to_num(t["pass_allcuts_merged_ML"].to_numpy().astype(float))
        m = t["MLPhoton_lead_mass"].to_numpy().astype(float)
        npass = int(f.sum())
        pass_ok = npass > 0
        print("  pass_merged  : %d / %d  (%.2f%%)  %s"
              % (npass, len(f), 100.0 * npass / max(len(f), 1),
                 "PASS" if pass_ok else "FAIL -- zero, check the join"))
        real = np.isfinite(m) & (m > -900)
        if real.sum():
            print("  lead_mass    : %d real, median %.4f, range [%.4f, %.4f]"
                  % (real.sum(), np.median(m[real]), m[real].min(), m[real].max()))

    ok = rec_ok and mt_ok and col_ok and pass_ok
    print("  => %s" % ("OK" if ok else "PROBLEM"))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("tags", nargs="*", default=None)
    ap.add_argument("--base", default=BASE)
    a = ap.parse_args()
    tags = a.tags if a.tags else TAGS
    results = [check(t, a.base) for t in tags]
    print("=" * 68)
    print("%d/%d tags OK" % (sum(results), len(results)))
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
