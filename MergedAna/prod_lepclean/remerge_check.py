"""Merge one HiggsDNA sample directory by scanning the disk, then reconcile.

Why not trust HiggsDNA's own --merge_outputs: it merges the job list stored in
analysis_manager.pkl, and run_analysis reports retired jobs as completed. The
2024 data friend ended up with 23/32 datasets partially merged that way (see
ref_hza_merged_data_friend_partial_merge). Here the inputs are whatever chunks
are on disk, and the verdict is fail-closed:

  PASS  iff  every job_* dir has its summary json (the completion marker)
        and  merged rows == sum of chunk rows
        and  (data only) no duplicated run:lumi:event

A job that selected zero events writes a summary but may write no chunk, so
the chunk count is reported but not required to equal the job count.

For MC the merge reproduces HiggsDNA's Task.merge_outputs exactly: weight_central is
multiplied by scale1fb * lumi (scale1fb = xs*bf*1000 / sum of the job sum_weights),
weight_central_no_lumi (scale1fb only), process_id and year are added, and every
systematic variation (output_job_*_<syst>.parquet -> merged_<syst>.parquet) is
written, not only nominal. The first version of this script merged raw chunks and
dropped the normalisation (signal weights ~0.85/event instead of ~0.02 -> limits off
by ~10^3); data has weight 1 and is merged nominal-only as before.

    python remerge_check.py <sample_dir> [--data] [--no-merge]

<sample_dir> is the directory holding job_*/ (e.g.
.../parquet_friend_lepclean/<tag>/<tag>_2024). Prints one line, exit 0 on PASS,
2 on FAIL.
"""
import argparse
import json
import re
import glob
import os
import subprocess
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def read(f, columns=None):
    for k in range(5):
        try:
            return pq.read_table(f, columns=columns)
        except OSError:
            time.sleep(10 * (k + 1))
    raise OSError(f"unreadable after retries: {f}")


def _set(t, name, col):
    i = t.schema.get_field_index(name)
    return t.set_column(i, name, col) if i >= 0 else t.append_column(name, col)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample_dir")
    ap.add_argument("--data", action="store_true", help="also check run:lumi:event duplicates")
    ap.add_argument("--no-merge", action="store_true", help="reconcile only, keep the merged file")
    a = ap.parse_args()

    d = a.sample_dir.rstrip("/")
    name = os.path.basename(d)
    jobdirs = [x for x in os.listdir(d) if x.startswith("job_")]
    n_jobs = len(jobdirs)
    # find, not a shell glob: >10k entries overflow a glob and come back as 0
    found = subprocess.run(["find", d, "-name", "*.json", "-path", "*/job_*/*"],
                           capture_output=True, text=True, check=True).stdout.split()
    done = {os.path.basename(os.path.dirname(f)) for f in found if "_summary_job" in f}
    missing = sorted(set(jobdirs) - done, key=lambda s: int(s.split("_")[1]))
    chunks = sorted(subprocess.run(["find", d, "-name", "output_job_*_nominal.parquet"],
                                   capture_output=True, text=True, check=True).stdout.split())
    chunk_rows = sum(pq.ParquetFile(f).metadata.num_rows for f in chunks)

    dest = f"{d}/merged_nominal.parquet"
    cfgs = sorted(glob.glob(f"{d}/job_*/*_config_job*.json"))
    sample = json.load(open(cfgs[0]))["sample"] if cfgs else {"is_data": True}
    is_mc = not sample.get("is_data", True)
    if is_mc:
        sumw = sum(json.load(open(f)).get("sum_weights", 0.0)
                   for f in found if "_summary_job" in f)
        scale1fb = sample["norm_factor"] * 1000.0 / sumw if sumw > 0 else 0.0
        lumi = float(sample["lumi"])
    if not a.no_merge and chunks:
        if is_mc:
            allc = subprocess.run(["find", d, "-name", "output_job_*.parquet"], capture_output=True,
                                  text=True, check=True).stdout.split()
            systs = sorted({os.path.basename(f).split("_", 3)[3][:-len(".parquet")] for f in allc})
        else:
            systs = ["nominal"]
        for syst in systs:
            files = sorted(f for f in subprocess.run(
                ["find", d, "-name", f"output_job_*_{syst}.parquet"], capture_output=True,
                text=True, check=True).stdout.split()
                if re.fullmatch(r"output_job_\d+_" + re.escape(syst) + r"\.parquet", os.path.basename(f)))
            if not files:
                continue
            schema = pa.unify_schemas([pq.read_schema(f) for f in files])
            if is_mc:
                for c in [n for n in schema.names if n.startswith("weight_")]:
                    schema = schema.set(schema.get_field_index(c), pa.field(c, pa.float64()))
                extra = [pa.field("weight_central_no_lumi", pa.float64()),
                         pa.field("process_id", pa.float64()), pa.field("year", pa.string())]
                for fld in extra:
                    schema = (schema.set(schema.get_field_index(fld.name), fld)
                              if fld.name in schema.names else schema.append(fld))
            out = f"{d}/merged_{syst}.parquet"
            with pq.ParquetWriter(out + ".tmp", schema) as w:
                for f in files:
                    t = read(f)
                    if is_mc:
                        cw = pc.cast(t.column("weight_central"), pa.float64())
                        t = _set(t, "weight_central", pc.multiply(cw, pa.scalar(scale1fb * lumi, pa.float64())))
                        t = _set(t, "weight_central_no_lumi", pc.multiply(cw, pa.scalar(scale1fb, pa.float64())))
                        t = _set(t, "process_id", pa.array([float(sample.get("process_id", 0))] * len(t), pa.float64()))
                        t = _set(t, "year", pa.array([str(sample.get("year", ""))] * len(t), pa.string()))
                    for n in schema.names:
                        if n not in t.column_names:
                            t = t.append_column(n, pa.nulls(len(t), type=schema.field(n).type))
                    w.write_table(t.select(schema.names).cast(schema))
            os.replace(out + ".tmp", out)
    merged_rows = pq.ParquetFile(dest).metadata.num_rows if os.path.exists(dest) else -1

    extra = ""
    dup = 0
    if os.path.exists(dest) and merged_rows > 0:
        cols = pq.read_schema(dest).names
        if "n_MLPhoton_preclean" in cols:
            pre = read(dest, ["n_MLPhoton_preclean"])["n_MLPhoton_preclean"].to_numpy()
            extra += f" nojoin={np.mean(pre == 0):.2%}"
        else:
            extra += " nojoin=NO_COLUMN"
        if a.data:
            t = read(dest, ["run", "luminosityBlock", "event"])
            k = np.stack([t[c].to_numpy().astype(np.int64) for c in t.column_names], 1)
            dup = len(k) - len(np.unique(k, axis=0))
            extra += f" dup={dup}"

    ok = (n_jobs > 0 and not missing and merged_rows == chunk_rows and dup == 0
          and "NO_COLUMN" not in extra)
    miss = ",".join(m.split("_")[1] for m in missing[:20]) + ("..." if len(missing) > 20 else "")
    if is_mc:
        extra += f" scale1fb={scale1fb:.6g} lumi={lumi}"
    print(f"{'PASS' if ok else 'FAIL'} {name} jobs={n_jobs} done={len(done & set(jobdirs))} "
          f"missing={len(missing)}[{miss}] chunks={len(chunks)} chunk_rows={chunk_rows} "
          f"merged_rows={merged_rows}{extra}", flush=True)
    sys.exit(0 if ok else 2)


if __name__ == "__main__":
    main()
