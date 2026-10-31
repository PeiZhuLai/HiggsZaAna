#!/usr/bin/env python3
"""Put the rows of the re-produced signal parquets back in the row order of the previous production.

Why: Parque2Root_BDT.py --split assigns train/validation/test by row index (row % 10). The BDT was
trained on the train split of the previous production (parquet_DNA_tmp_fsrfix). The re-production
selects exactly the same events but merges them in a different order, so the new test split would
contain events the BDT was trained on (mA_M10_2023preBPix: only 506 of 5466 test events in common)
and the signal efficiency would be biased. Reordering to the old (run, luminosityBlock, event) order
restores the identical splits.

For every <mA>_<era> and every merged_<syst>.parquet: the event sets must be identical (else abort
that file); the output is written to OUT/<mA>_<era>/merged_<syst>.parquet with the old row order.

usage: reorder_sig_sfunc_20261024.py NEW_SIG_MC_DIR OUT_SIG_MC_DIR [sample_year ...]
"""
import os
import sys

import numpy as np
import pyarrow.parquet as pq

OLD = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
KEY = ["run", "luminosityBlock", "event"]


def keys(table):
    cols = [table.column(k).to_numpy() for k in KEY]
    return np.rec.fromarrays(cols, names=KEY)


def main():
    new_dir, out_dir = sys.argv[1], sys.argv[2]
    samples = sys.argv[3:] or sorted(d for d in os.listdir(new_dir) if d.startswith("mA_M"))
    bad = 0
    for sy in samples:
        os.makedirs(os.path.join(out_dir, sy), exist_ok=True)
        for fname in sorted(f for f in os.listdir(os.path.join(new_dir, sy)) if f.startswith("merged_") and f.endswith(".parquet")):
            old_keys = keys(pq.read_table(os.path.join(OLD, sy, fname), columns=KEY))
            new_tab = pq.read_table(os.path.join(new_dir, sy, fname))
            new_keys = keys(new_tab)
            if len(old_keys) != len(new_keys) or len(np.unique(new_keys)) != len(new_keys):
                print("!! %s/%s: sizes %d vs %d or duplicate keys" % (sy, fname, len(old_keys), len(new_keys)))
                bad += 1
                continue
            order_new = np.argsort(new_keys, kind="stable")
            pos_old = np.searchsorted(new_keys[order_new], old_keys)
            pos_old = np.clip(pos_old, 0, len(new_keys) - 1)
            idx = order_new[pos_old]
            if not np.array_equal(new_keys[idx], old_keys):
                print("!! %s/%s: event sets differ" % (sy, fname))
                bad += 1
                continue
            out = new_tab.take(idx)
            assert np.array_equal(keys(out), old_keys)
            path = os.path.join(out_dir, sy, fname)
            pq.write_table(out, path + ".tmp")
            os.replace(path + ".tmp", path)
        print("ok %s" % sy, flush=True)
    print("SUMMARY bad=%d" % bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
