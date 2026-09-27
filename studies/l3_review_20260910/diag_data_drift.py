#!/usr/bin/env python
"""Split the Data drift by lepton channel.

The FSR recovery is muon-only, so a drift caused by the FSR fix must appear in
the muon channel and be ~0 in the electron channel -- the same null test that
validated the sigma_eff improvement. A drift of the same size in BOTH channels
means the cause is something else (different input files, lumi/golden JSON, or
another code change between the two productions).
"""
import pyarrow.parquet as pq
import numpy as np, os

NEW = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Data"
OLD = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA/Data"

def counts(path):
    t = pq.read_table(path, columns=["z_ee", "z_mumu"])
    ee = np.asarray(t["z_ee"]).astype(bool)
    mm = np.asarray(t["z_mumu"]).astype(bool)
    return len(ee), int(ee.sum()), int(mm.sum())

print("%-20s %10s %10s %10s   %10s %10s" % ("era", "tot_old", "tot_new", "d_tot", "d_ee", "d_mumu"))
for era in ("Data_2022preEE", "Data_2022postEE", "Data_2023preBPix", "Data_2023postBPix", "Data_2024"):
    po = os.path.join(OLD, era, "merged_nominal.parquet")
    pn = os.path.join(NEW, era, "merged_nominal.parquet")
    if not (os.path.exists(po) and os.path.exists(pn)):
        print("%-20s  (missing: old=%s new=%s)" % (era, os.path.exists(po), os.path.exists(pn)));  continue
    to, eo, mo = counts(po)
    tn, en, mn = counts(pn)
    f = lambda a, b: "%+.2f%%" % (100.0 * (b - a) / a) if a else "n/a"
    print("%-20s %10d %10d %10s   %10s %10s" % (era, to, tn, f(to, tn), f(eo, en), f(mo, mn)))
