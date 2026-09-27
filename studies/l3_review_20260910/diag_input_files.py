#!/usr/bin/env python
"""Compare the INPUT FILE LISTS of the old and new Data productions, per era.

A yield drift can be physics (the FSR fix moving events across the selection)
or bookkeeping (a different set of input files). Only the second is a problem,
and the two are told apart by counting the files each production actually read.
2024 was still being collected when the old production ran, so a genuine
increase there is expected.
"""
import json, glob, os
NEW = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Data"
OLD = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA/Data"

def files_of(base, era):
    out = set()
    for cfg in glob.glob(os.path.join(base, era, "job_*", "*_config_job*.json")):
        try:
            out.update(json.load(open(cfg)).get("files") or [])
        except Exception:
            pass
    return out

print("%-20s %9s %9s %9s %9s" % ("era", "old_files", "new_files", "only_old", "only_new"))
for era in ("Data_2022preEE", "Data_2022postEE", "Data_2023preBPix", "Data_2023postBPix", "Data_2024"):
    o = files_of(OLD, era); n = files_of(NEW, era)
    print("%-20s %9d %9d %9d %9d" % (era, len(o), len(n), len(o - n), len(n - o)))
