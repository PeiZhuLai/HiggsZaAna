#!/usr/bin/env bash
# Production rerun of the data MVA cut with the local-write apply_bdt_data.py (user 2026-10-17), to replace the 3 data
# files whose pandas index branch was damaged by the EOS FUSE write. Physics content does not change, so the
# workspaces (ws/) and everything downstream stay valid. Steps:
#   1 server-side backup of root_MVAcut/data (ws/ included) to data.bak_preLocalWrite_20261017/data
#   2 apply_bdt_data.py with its production defaults (scores checked against the installed models)
#   3 every branch of the 180 files reads, and every file equals the test output _localwrite_test2_20261017/data
#     (identical to the original per-mA logic) bit for bit, index included
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/prod_mvacut_integrity_20261017
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
X=/eos/user/p/pelai/HZa/root_MVAcut
xrdfs eosuser.cern.ch stat $X/data.bak_preLocalWrite_20261017 >/dev/null 2>&1 && { echo "backup exists; stop"; exit 1; }
# file by file (server-side xrdcp -r copied the first file of every mA dir onto the dir path), adler32 compared
SRC=/eos/home-p/pelai/HZa/root_MVAcut/data
for f in $(cd $SRC && find . -type f -name '*.root' | sed 's#^\./##'); do
  xrdfs eosuser.cern.ch mkdir -p $X/data.bak_preLocalWrite_20261017/data/$(dirname $f)
  xrdcp -s root://eosuser.cern.ch/$X/data/$f root://eosuser.cern.ch/$X/data.bak_preLocalWrite_20261017/data/$f || { echo "backup of $f failed"; exit 1; }
  a=$(xrdfs eosuser.cern.ch query checksum $X/data/$f | awk '{print $2}'); b=$(xrdfs eosuser.cern.ch query checksum $X/data.bak_preLocalWrite_20261017/data/$f | awk '{print $2}')
  [ -n "$a" ] && [ "$a" = "$b" ] || { echo "backup checksum of $f differs ($a vs $b)"; exit 1; }
done
nb=$(find /eos/home-p/pelai/HZa/root_MVAcut/data.bak_preLocalWrite_20261017/data -type f -name '*.root' | wc -l); no=$(find /eos/home-p/pelai/HZa/root_MVAcut/data -type f -name '*.root' | wc -l)
echo "backup files $nb / original $no"; [ "$nb" -eq "$no" ] || { echo "backup incomplete"; exit 1; }
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd $FF && eval "$(scramv1 runtime -sh)"; set -u
s=$(date +%s); ( cd $FF && python3 MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py ) > $D/logs/rerun_prod_data.log 2>&1; rc=$?
echo "apply_bdt_data rc=$rc $(( $(date +%s)-s )) s"; [ $rc -eq 0 ] || exit 1
env -i HOME=$HOME PATH=/usr/bin:/bin /eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3 - <<'PY'
import uproot, numpy as np
P = "/eos/home-p/pelai/HZa/root_MVAcut/data"; T = "/eos/home-p/pelai/HZa/sandbox_varsearch_v1/_localwrite_test2_20261017/data"
n = bad = 0
for m in range(1, 31):
    for f in ["output_2022preEE", "output_2022postEE", "output_2023preBPix", "output_2023postBPix", "output_2024", "run3"]:
        a = uproot.open("%s/mA_M%d/%s.root" % (P, m, f))["DiphotonTree/Data_13p6TeV"]; b = uproot.open("%s/mA_M%d/%s.root" % (T, m, f))["DiphotonTree/Data_13p6TeV"]
        n += 1
        if a.keys() != b.keys(): bad += 1; print("BRANCHES differ mA%d %s" % (m, f)); continue
        for br in a.keys():
            x, y = a[br].array(library="np"), b[br].array(library="np")
            if len(x) != a.num_entries or x.dtype != y.dtype or not np.array_equal(x, y, equal_nan=x.dtype.kind == "f"):
                bad += 1; print("DIFF mA%d %s %s" % (m, f, br)); break
print("production data after rerun: %d files, every branch read; %d differ from the test output" % (n, bad))
PY
