#!/usr/bin/env python3
"""Three-stage event-count reconciliation for the HZa signal re-production.

Rules this encodes (from CLAUDE.md, learned the hard way):
  * enumerate chunks with find, never `ls *.parquet` -- a glob over >10k files
    overflows and silently returns 0, which reads as "no chunks" rather than as
    an error;
  * EXCLUDE merged_nominal.parquet from the chunk sum, otherwise it is counted
    twice and the sample looks ~2x complete;
  * read only ParquetFile(f).metadata.num_rows, never the data;
  * merge-before-rechunk race: merged_nominal.parquet must be at least as new as
    every chunk it claims to contain. chunk_rows == merged_rows can still pass
    while merge silently missed a late chunk, so compare mtimes as well.

The ROOT stage is checked separately after p2root; this script covers
chunk -> merged and the mtime guard.

Input : <staging>/Sig_MC/<sample>_<era>/{job_*/output_job_*_nominal.parquet,
        merged_nominal.parquet}
Output: printed per-sample table + non-zero exit if anything fails

Usage : python reconcile_signal.py [--base <staging Sig_MC dir>] [--ref <old prod dir>]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pyarrow.parquet as pq

DEFAULT_BASE = ("/eos/project/h/htozg-dy-privatemc/pelai/HZa/"
                "parquet_DNA_tmp_fsrfix/Sig_MC")
DEFAULT_REF = "/eos/home-p/pelai/HZa/parquet_DNA/Sig_MC"
MERGED = "merged_nominal.parquet"


def chunk_files(sample_dir: Path):
    """All nominal chunks under the sample, excluding the merged file itself."""
    out = []
    for root, _dirs, files in os.walk(sample_dir):
        for name in files:
            if not name.endswith("_nominal.parquet"):
                continue
            if name == MERGED:
                continue
            out.append(Path(root) / name)
    return out


def num_rows(path: Path) -> int:
    return pq.ParquetFile(path).metadata.num_rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default=DEFAULT_BASE)
    parser.add_argument("--ref", default=DEFAULT_REF,
                        help="Previous production, for a size sanity check. "
                             "Counts are NOT expected to match exactly: the FSR "
                             "fix moves a few events across the selection.")
    parser.add_argument("--max-drift", type=float, default=1.0,
                        help="Percent change vs --ref that is still considered sane.")
    args = parser.parse_args()

    base = Path(args.base)
    ref = Path(args.ref)
    if not base.is_dir():
        print("[FAIL] base directory does not exist: %s" % base)
        return 2

    samples = sorted(d for d in base.iterdir() if d.is_dir())
    print("%-26s %9s %9s %-7s %-7s %10s" %
          ("sample", "chunks", "merged", "rows", "mtime", "vs ref"))
    print("%-26s %9s %9s %-7s %-7s %10s" %
          ("------", "------", "------", "-----", "-----", "------"))

    failures = []
    for sample in samples:
        chunks = chunk_files(sample)
        merged = sample / MERGED
        chunk_rows = sum(num_rows(f) for f in chunks) if chunks else 0

        if not merged.exists():
            print("%-26s %9d %9s %-7s %-7s %10s"
                  % (sample.name, chunk_rows, "-", "NO MERGE", "-", "-"))
            failures.append("%s: merged_nominal.parquet missing" % sample.name)
            continue

        merged_rows = num_rows(merged)
        rows_ok = chunk_rows == merged_rows

        # merge-before-rechunk race
        newest_chunk = max((f.stat().st_mtime for f in chunks), default=0.0)
        mtime_ok = merged.stat().st_mtime >= newest_chunk

        drift = "-"
        ref_merged = ref / sample.name / MERGED
        if ref_merged.exists():
            ref_rows = num_rows(ref_merged)
            if ref_rows:
                pct = 100.0 * (merged_rows - ref_rows) / ref_rows
                drift = "%+.2f%%" % pct
                if abs(pct) > args.max_drift:
                    failures.append("%s: %+.2f%% vs previous production (%d -> %d)"
                                    % (sample.name, pct, ref_rows, merged_rows))

        print("%-26s %9d %9d %-7s %-7s %10s"
              % (sample.name, chunk_rows, merged_rows,
                 "OK" if rows_ok else "MISMATCH",
                 "OK" if mtime_ok else "STALE", drift))

        if not rows_ok:
            failures.append("%s: chunk_rows %d != merged_rows %d"
                            % (sample.name, chunk_rows, merged_rows))
        if not mtime_ok:
            failures.append("%s: merged is OLDER than its newest chunk -- "
                            "re-merge before going downstream" % sample.name)

    print()
    if failures:
        print("[FAIL] %d problem(s):" % len(failures))
        for f in failures:
            print("   - %s" % f)
        return 1
    print("[OK] chunk == merged and no stale merge for all %d samples" % len(samples))
    return 0


if __name__ == "__main__":
    sys.exit(main())
