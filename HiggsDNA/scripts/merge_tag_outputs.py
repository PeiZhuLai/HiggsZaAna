#!/usr/bin/env python3
"""
Standalone, HiggsDNA-faithful merge of one DATA task's per-job parquet chunks
into merged_nominal.parquet, bypassing the (occasionally-hanging) condor monitor.

It reuses the exact internal helpers of higgs_dna.job_management.task so the output
is byte-equivalent to what Task.merge_outputs() would have produced for a data sample:
  - unify per-chunk (normalized) schemas
  - add process_id (float64) and year (string) constant columns (values from job summary)
  - stream row-groups into a single ParquetWriter
No scale1fb / lumi scaling (that path is MC-only; these are data).

Usage:
  merge_tag_outputs.py <task_output_dir>
where <task_output_dir> is e.g.
  .../parquet_friend/Data_2024_EGamma0_RunD_v1/Data_2024_EGamma0_RunD_v1_2024
(the inner dir that holds the job_* subdirs). A tag's top dir is also accepted;
we descend into its single *_2024 subdir.

Exit 0 on success (merged rows == sum of chunk rows), non-zero otherwise.
"""
import glob
import json
import os
import sys

import pyarrow as pa
import pyarrow.parquet as pq

from higgs_dna.job_management.task import (
    _get_merged_output_schema,
    _replace_or_add_field,
    _iter_parquet_row_groups,
    _replace_or_add_column,
    _constant_column,
    _normalize_table_to_schema,
)


def resolve_task_dir(arg):
    arg = arg.rstrip("/")
    # if this dir already has job_* subdirs, use it
    if glob.glob(os.path.join(arg, "job_*")):
        return arg
    # else descend into the single *_2024 task subdir
    subs = [d for d in glob.glob(os.path.join(arg, "*_2024")) if os.path.isdir(d)]
    if len(subs) == 1:
        return subs[0]
    raise SystemExit("Cannot resolve task dir from %s (subs=%s)" % (arg, subs))


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    task_dir = resolve_task_dir(sys.argv[1])
    name = os.path.basename(task_dir)

    summaries = sorted(glob.glob(os.path.join(task_dir, "job_*", "*summary_job*.json")))
    outputs = []
    is_data = None
    process_id = None
    year = None
    for f in summaries:
        d = json.load(open(f))
        if not d.get("successful"):
            continue
        nom = d.get("outputs", {}).get("nominal")
        if not nom or not os.path.exists(nom):
            continue
        outputs.append(nom)
        if is_data is None:
            s = d.get("config", {}).get("sample", {})
            is_data = s.get("is_data")
            process_id = s.get("process_id")
            year = s.get("year")
    outputs = sorted(set(outputs))

    if not outputs:
        raise SystemExit("[merge_tag_outputs] no successful outputs in %s" % task_dir)
    if not is_data:
        raise SystemExit("[merge_tag_outputs] refusing: sample is not data (is_data=%s); MC needs scale1fb path." % is_data)

    print("[merge_tag_outputs] task=%s successful_chunks=%d process_id=%s year=%s"
          % (name, len(outputs), process_id, year))

    merged_output = os.path.join(task_dir, "merged_nominal.parquet")
    merged_schema = _get_merged_output_schema(outputs, True, name, "nominal")
    if process_id is not None:
        merged_schema = _replace_or_add_field(merged_schema, pa.field("process_id", pa.float64()))
    if year is not None:
        merged_schema = _replace_or_add_field(merged_schema, pa.field("year", pa.string()))

    tmp = merged_output + ".tmp"
    for p in (merged_output, tmp):
        if os.path.exists(p):
            os.remove(p)

    n_chunk_rows = 0
    writer = None
    try:
        for i, output in enumerate(outputs, 1):
            if i == 1 or i == len(outputs) or i % 25 == 0:
                print("[merge_tag_outputs] %d/%d" % (i, len(outputs)))
            for table in _iter_parquet_row_groups(output):
                n_chunk_rows += len(table)
                if process_id is not None:
                    table = _replace_or_add_column(
                        table, "process_id",
                        _constant_column(float(process_id), len(table), pa.float64()))
                if year is not None:
                    table = _replace_or_add_column(
                        table, "year",
                        _constant_column(str(year), len(table), pa.string()))
                table = _normalize_table_to_schema(table, merged_schema, name, "nominal", output)
                if writer is None:
                    writer = pq.ParquetWriter(tmp, merged_schema)
                writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()

    os.replace(tmp, merged_output)
    merged_rows = pq.ParquetFile(merged_output).metadata.num_rows
    print("[merge_tag_outputs] merged_rows=%d sum_chunk_rows=%d -> %s"
          % (merged_rows, n_chunk_rows, "OK" if merged_rows == n_chunk_rows else "MISMATCH"))
    if merged_rows != n_chunk_rows:
        raise SystemExit(2)
    print("[merge_tag_outputs] wrote %s" % merged_output)


if __name__ == "__main__":
    main()
