"""Stage 3: ROOT 'inclusive' entries == merged rows, ROOT newer than merged, and train+validation+test == inclusive.
Also prints old (no-veto) run3_bdt_inputs_fsrfix entries for reference.
Input : <parquet base> <root base> samples...   Output: stdout, rc
"""
import os, sys
import uproot, pyarrow.parquet as pq
IN, RO = sys.argv[1:3]; rc = 0
OLDR = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
for s in sys.argv[3:]:
    mp = os.path.join(IN, s + "_2024", "merged_nominal.parquet")
    mrow = pq.ParquetFile(mp).metadata.num_rows
    rp = os.path.join(RO, s, "2024.root")
    if not os.path.exists(rp):
        print("%-14s MISSING %s FAIL" % (s, rp)); rc = 1; continue
    f = uproot.open(rp)
    n = f["inclusive"].num_entries
    nsplit = sum(f[k].num_entries for k in ("train", "validation", "test") if k in f)
    br = f["inclusive"].keys()
    has_sb = all(b in br for b in ("sideband_rwgt", "weight_sideband_rwgt", "factor_sideband_rwgt"))
    fresh = os.path.getmtime(rp) >= os.path.getmtime(mp)
    old = uproot.open(os.path.join(OLDR, s, "2024.root"))["inclusive"].num_entries
    ok = n == mrow and fresh and nsplit == n and has_sb
    rc |= (not ok)
    print("%-14s merged=%8d root_inclusive=%8d train+val+test=%8d sideband_branches=%s fresh=%s | old(no veto)=%8d ratio=%.3f %s"
          % (s, mrow, n, nsplit, has_sb, fresh, old, n / old if old else float('nan'), "OK" if ok else "FAIL"))
sys.exit(rc)
