#!/usr/bin/env python
"""List the 12 input files present in the old Data_2022preEE production but not the new one."""
import json, glob, os
NEW = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Data/Data_2022preEE"
OLD = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA/Data/Data_2022preEE"
def files_of(base):
    out=set()
    for cfg in glob.glob(os.path.join(base,"job_*","*_config_job*.json")):
        try: out.update(json.load(open(cfg)).get("files") or [])
        except Exception: pass
    return out
o, n = files_of(OLD), files_of(NEW)
miss = sorted(o - n)
print("missing from the new production: %d" % len(miss))
for f in miss: print(" ", f)
# which primary dataset do they belong to?
from collections import Counter
c = Counter("/".join(f.split("/store/")[-1].split("/")[:4]) for f in miss)
print("\nby dataset:")
for k, v in c.most_common(): print("  %-70s %d" % (k, v))
