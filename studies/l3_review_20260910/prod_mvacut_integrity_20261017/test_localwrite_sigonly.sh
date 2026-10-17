#!/usr/bin/env bash
# Signal-only rerun (data part already passed) of: Test of the local-write + xrdcp version of apply_bdt_{data,sig}.py: outputs go to a sandbox EOS dir (-o), then are
# compared tree by tree and branch by branch with the production outputs (checked good by check_prod_mvacut.sh).
# Data: all 30 mA; signal: mA_M3 2022preEE. Production outputs are only read. The pandas index branch is not compared
# (damaged in 3 production data files) but must read in every new file.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/prod_mvacut_integrity_20261017
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
T=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/_localwrite_test_20261017
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd $FF && eval "$(scramv1 runtime -sh)"; set -u
# data part not rerun: passed in test_localwrite.sh, logs/test_data.log kept
s=$(date +%s); ( cd $FF && python3 MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples mA_M3 --years 2022preEE -o $T/sig ) > $D/logs/test_sig.log 2>&1; echo "sig rc=$? $(( $(date +%s)-s )) s"
ls $FF/MVAcut/run3_ReReco_Sys/staging/ 2>/dev/null | head -3
env -i HOME=$HOME PATH=/usr/bin:/bin /eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3 - $T <<'PY'
import sys, uproot, numpy as np
T = sys.argv[1]; P = "/eos/home-p/pelai/HZa/root_MVAcut"; bad = 0; n = 0
def same(a, b):
    fa, fb = uproot.open(a), uproot.open(b)
    ta = {k for k, v in fa.items(cycle=False) if hasattr(v, "num_entries")}
    tb = {k for k, v in fb.items(cycle=False) if hasattr(v, "num_entries")}
    if ta != tb: return "tree lists differ"
    for t in ta:
        if set(fa[t].keys()) != set(fb[t].keys()): return t + ": branch lists differ"
        for br in fa[t].keys():
            x = fa[t][br].array(library="np")   # new output: every branch must read
            if br == "index": continue          # production: 3 data files have a damaged pandas index branch
            y = fb[t][br].array(library="np")
            if x.shape != y.shape or not np.array_equal(x, y, equal_nan=x.dtype.kind == "f"): return "%s:%s differs" % (t, br)
    return None
files = ["data/mA_M%d/%s.root" % (m, f) for m in range(1, 31) for f in ["output_2022preEE", "output_2022postEE", "output_2023preBPix", "output_2023postBPix", "output_2024", "run3"]]
files += ["sig/mA_M3/output_2022preEE.root"]
for f in files:
    r = same("%s/%s" % (T, f), "%s/%s" % (P, f)); n += 1
    if r: bad += 1; print("DIFF", f, r)
print("compared %d files with production: %d differ" % (n, bad))
PY
