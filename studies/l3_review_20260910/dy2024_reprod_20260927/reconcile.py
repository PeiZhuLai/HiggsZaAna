"""Stage 1+2 reconciliation + completeness + merge-race check for the dyveto production.
Input : <base>/<sample>_2024/job_*/ (summary json = completion marker; chunk parquet), merged_nominal.parquet
Output: stdout; rc=0 only if every job ran, chunk_rows == merged_rows and no chunk is newer than merged.
Chunk list via os.walk (not a shell glob) and EXCLUDING merged_nominal.parquet (else ~2x).
"""
import os, sys, glob
import pyarrow.parquet as pq
base = sys.argv[1]; rc = 0
for s in sys.argv[2:]:
    d = os.path.join(base, s + "_2024")
    jobs = sorted(x for x in os.listdir(d) if x.startswith("job_"))
    notrun = [j for j in jobs if not glob.glob(os.path.join(d, j, "*_summary_%s.json" % j.replace("job_", "job")))]
    chunks = []
    for root, _, files in os.walk(d):
        for f in files:
            if f.endswith("_nominal.parquet") and f != "merged_nominal.parquet":
                chunks.append(os.path.join(root, f))
    crow = sum(pq.ParquetFile(f).metadata.num_rows for f in chunks)
    mp = os.path.join(d, "merged_nominal.parquet")
    mrow = pq.ParquetFile(mp).metadata.num_rows if os.path.exists(mp) else (0 if not chunks else -1)
    mm = os.path.getmtime(mp) if os.path.exists(mp) else 0
    late = [f for f in chunks if os.path.getmtime(f) > mm]
    ok = (not notrun) and crow == mrow and not late
    rc |= (not ok)
    print("%-14s jobs=%4d ran=%4d chunks=%4d chunk_rows=%8d merged_rows=%8d late_chunks=%d %s%s"
          % (s, len(jobs), len(jobs) - len(notrun), len(chunks), crow, mrow, len(late), "OK" if ok else "FAIL",
             (" notrun=" + ",".join(notrun[:20])) if notrun else ""))
sys.exit(rc)
