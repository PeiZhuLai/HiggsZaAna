"""Re-merge data friend chunks into merged_nominal.parquet.

The existing merged files hold only part of their chunks (22 of 32 2024
datasets; worst Muon0 RunG 0.5%), which is where the "data has only 47.5% ML
coverage" came from. Chunks are complete on disk, so a re-merge fixes it.

Streams one chunk per row group, writes <merged>.tmp and renames, then checks
merged_rows == sum(chunk_rows).

    python remerge_data_friend.py Data_2024_EGamma0_RunG_v2 [...]
"""
import os
import subprocess
import sys
import time

import pyarrow as pa
import pyarrow.parquet as pq

BASE = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend"


def read(f):
    for k in range(5):
        try:
            return pq.read_table(f)
        except OSError:
            time.sleep(10 * (k + 1))
    raise OSError(f"unreadable after retries: {f}")


def remerge(tag):
    sub = f"{BASE}/{tag}/{tag}_2024"
    out = subprocess.run(["find", sub, "-name", "output_job_*_nominal.parquet"],
                         capture_output=True, text=True, check=True).stdout.split()
    chunks = sorted(out)
    schema = pa.unify_schemas([pq.read_schema(f) for f in chunks])
    dest = f"{sub}/merged_nominal.parquet"
    tmp = dest + ".tmp"
    n_in = 0
    with pq.ParquetWriter(tmp, schema) as w:
        for f in chunks:
            t = read(f)
            n_in += t.num_rows
            w.write_table(t.select(schema.names).cast(schema) if t.schema != schema
                          else t)
    n_out = pq.ParquetFile(tmp).metadata.num_rows
    if n_out != n_in:
        raise RuntimeError(f"{tag}: wrote {n_out} rows but read {n_in}")
    os.replace(tmp, dest)
    print(f"{tag:32s} chunks={len(chunks):4d} rows={n_out}", flush=True)


if __name__ == "__main__":
    for tag in sys.argv[1:]:
        remerge(tag)
