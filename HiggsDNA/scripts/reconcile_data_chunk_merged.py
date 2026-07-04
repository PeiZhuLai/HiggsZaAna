#!/usr/bin/env python3
"""Three-stage reconciliation (chunk vs merged stage) for the 2024 data ML
friend production.

For every Data_2024_<tag> under parquet_friend/, per CLAUDE.md rules:
  1. chunk_rows  = Sum of ParquetFile(f).metadata.num_rows over every
                   output_*_nominal.parquet, found with os.walk (NOT ls glob,
                   which overflows past ~10k files -> false 0), EXCLUDING
                   merged_nominal.parquet itself (else it double-counts ~2x).
  2. merged_rows = merged_nominal.parquet num_rows.
  3. assert chunk_rows == merged_rows, AND merged mtime >= max(chunk mtime)
     (the merge-before-rechunk race fingerprint: a late chunk newer than the
     merged file means merge missed it).

Only metadata is read (fast). The root stage (p2root) is validated later,
downstream of merged_p2root. Run in env `higgs-alp-ana` (pyarrow).

Exit 0 if all tags pass; exit 1 if any tag fails or is incomplete.
"""
import os
import sys
import glob
import pyarrow.parquet as pq

BASE = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_friend"
EXPECT_TAGS = 32


def find_chunks(inner):
    """All output_*_nominal.parquet under inner, excluding merged_nominal."""
    out = []
    for root, _dirs, files in os.walk(inner):
        for fn in files:
            if not fn.endswith(".parquet"):
                continue
            if fn == "merged_nominal.parquet":
                continue
            out.append(os.path.join(root, fn))
    return out


def rows_and_mtime(path):
    return pq.ParquetFile(path).metadata.num_rows, os.path.getmtime(path)


def main():
    tag_dirs = sorted(glob.glob(f"{BASE}/Data_2024_*"))
    if not tag_dirs:
        print(f"[FATAL] no tag dirs under {BASE}")
        return 1
    print(f"Found {len(tag_dirs)} tag dirs (expect {EXPECT_TAGS})\n")
    print(f"{'tag':32s} {'chunks':>7s} {'chunk_rows':>12s} {'merged_rows':>12s} "
          f"{'rows':>5s} {'mtime':>6s}")
    print("-" * 90)

    fails, oks, incomplete = [], [], []
    for td in tag_dirs:
        tag = os.path.basename(td)
        inner = os.path.join(td, f"{tag}_2024")
        merged = os.path.join(inner, "merged_nominal.parquet")
        if not os.path.isdir(inner) or not os.path.isfile(merged):
            print(f"{tag:32s} {'--':>7s} {'--':>12s} {'--':>12s} "
                  f"{'MISS':>5s} {'--':>6s}  (no merged_nominal yet)")
            incomplete.append(tag)
            continue
        chunks = find_chunks(inner)
        crows = 0
        cmax_mtime = 0.0
        for f in chunks:
            try:
                r, mt = rows_and_mtime(f)
            except Exception as e:
                print(f"{tag:32s}  [WARN] unreadable chunk {f}: {e}")
                continue
            crows += r
            cmax_mtime = max(cmax_mtime, mt)
        mrows, mmtime = rows_and_mtime(merged)
        rows_ok = (crows == mrows)
        mtime_ok = (mmtime >= cmax_mtime)
        verdict = "OK" if rows_ok else "FAIL"
        mtv = "OK" if mtime_ok else "LATE"
        print(f"{tag:32s} {len(chunks):7d} {crows:12d} {mrows:12d} "
              f"{verdict:>5s} {mtv:>6s}")
        if rows_ok and mtime_ok:
            oks.append(tag)
        else:
            fails.append((tag, crows, mrows, rows_ok, mtime_ok))

    print("-" * 90)
    print(f"\nSUMMARY: {len(oks)} OK, {len(fails)} FAIL, {len(incomplete)} "
          f"incomplete (of {len(tag_dirs)} dirs; expect {EXPECT_TAGS})")
    if incomplete:
        print("  incomplete (no merged_nominal):", ", ".join(incomplete))
    for tag, cr, mr, ro, mo in fails:
        if not ro:
            print(f"  [FAIL] {tag}: chunk_rows {cr} != merged_rows {mr} "
                  f"(diff {cr - mr}) -> merge missed chunk(s); re-merge this tag")
        if not mo:
            print(f"  [LATE] {tag}: a chunk is newer than merged_nominal "
                  f"-> merge-before-rechunk race; re-merge this tag")
    ok = (not fails and not incomplete and len(tag_dirs) >= EXPECT_TAGS)
    print("\nRESULT:", "ALL PASS — safe to proceed to p2root" if ok
          else "NOT CLEAN — fix flagged tags before downstream")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
