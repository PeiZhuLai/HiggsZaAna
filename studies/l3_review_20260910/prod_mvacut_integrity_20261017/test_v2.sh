#!/usr/bin/env bash
# Second test of apply_bdt_{data,sig}.py: read-once data step (#3) + score provenance check (#2) + local write (shared
# eos_safe_write in HiggsZaAna/Parquet2Rootfile). Data (30 mA) must be identical to the first test's outputs
# (_localwrite_test_20261017: original per-mA full reads), including the pandas index, branch order and dtypes;
# signal mA_M3 2022preEE must equal production (jagged branches by lengths + flat values). Production is only read.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/prod_mvacut_integrity_20261017
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T1=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/_localwrite_test_20261017
T2=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/_localwrite_test2_20261017
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd $FF && eval "$(scramv1 runtime -sh)"; set -u
s=$(date +%s); ( cd $FF && python3 MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py -o $T2/data ) > $D/logs/test2_data.log 2>&1; echo "data rc=$? $(( $(date +%s)-s )) s"
s=$(date +%s); ( cd $FF && python3 MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples mA_M3 --years 2022preEE -o $T2/sig ) > $D/logs/test2_sig.log 2>&1; echo "sig rc=$? $(( $(date +%s)-s )) s"
grep -h "Scores of" $D/logs/test2_data.log | head -2; grep -c "Scores of" $D/logs/test2_sig.log
env -i HOME=$HOME PATH=/usr/bin:/bin /eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3 - $T1 $T2 <<'PY'
import sys, uproot, numpy as np
T1, T2 = sys.argv[1:3]; P = "/eos/home-p/pelai/HZa/root_MVAcut"
def eq(x, y):
    if x.dtype == object:
        return [len(v) for v in x] == [len(v) for v in y] and np.array_equal(np.concatenate(list(x) + [np.zeros(0)]).astype(float), np.concatenate(list(y) + [np.zeros(0)]).astype(float), equal_nan=True)
    return x.dtype == y.dtype and np.array_equal(x, y, equal_nan=x.dtype.kind == "f")
def compare(a, b, skip=()):
    fa, fb = uproot.open(a), uproot.open(b)
    ta = [k for k, v in fa.items(cycle=False) if hasattr(v, "num_entries")]
    tb = [k for k, v in fb.items(cycle=False) if hasattr(v, "num_entries")]
    if sorted(ta) != sorted(tb): return "tree lists differ"
    for t in ta:
        if [k for k in fa[t].keys() if k not in skip] != [k for k in fb[t].keys() if k not in skip]: return t + ": branch list/order differs"
        for br in fa[t].keys():
            if br in skip: continue
            if not eq(fa[t][br].array(library="np"), fb[t][br].array(library="np")): return "%s:%s differs" % (t, br)
    return None
n = bad = 0
for m in range(1, 31):
    for f in ["output_2022preEE", "output_2022postEE", "output_2023preBPix", "output_2023postBPix", "output_2024", "run3"]:
        r = compare("%s/data/mA_M%d/%s.root" % (T2, m, f), "%s/data/mA_M%d/%s.root" % (T1, m, f)); n += 1
        if r: bad += 1; print("DATA DIFF mA%d %s: %s" % (m, f, r))
print("data: %d files compared with the original-logic outputs (all branches incl. index): %d differ" % (n, bad))
r = compare(T2 + "/sig/mA_M3/output_2022preEE.root", P + "/sig/mA_M3/output_2022preEE.root")
print("signal mA_M3 2022preEE vs production:", r or "identical")
PY
