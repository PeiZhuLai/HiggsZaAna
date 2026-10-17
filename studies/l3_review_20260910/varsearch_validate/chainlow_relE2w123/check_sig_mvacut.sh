#!/usr/bin/env bash
# After the damaged data file (redo_data_mvacut.sh): every m_a 1-3 signal MVA-cut file of chainlow_relE2w123 is read
# with ROOT through xrootd (every entry of every tree) and with uproot through FUSE; per-tree entries must agree and
# ROOT must print no error. Output: logs/check_sig_mvacut.log
set -uo pipefail
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chainlow_relE2w123
FU=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2w123/root_MVAcut/sig
XR=root://eosuser.cern.ch//eos/user/p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2w123/root_MVAcut/sig
env -i HOME=$HOME PATH=/usr/bin:/bin /eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3 - $FU > $CH/logs/check_sig_uproot.json <<'PY'
import sys, json, uproot
B = sys.argv[1]; out = {}
for m in (1, 2, 3):
    for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
        f = uproot.open("%s/mA_M%d/output_%s.root" % (B, m, y))
        out["%d_%s" % (m, y)] = {k: v.num_entries for k, v in f.items(cycle=False) if hasattr(v, "num_entries")}
print(json.dumps(out))
PY
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh
eval "$(cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && scramv1 runtime -sh)"; set -u
python3 - $XR $CH/logs/check_sig_uproot.json <<'PY' > $CH/logs/check_sig_mvacut.log 2>&1
import sys, json, ROOT
ROOT.gROOT.SetBatch(True)
X, U = sys.argv[1], json.load(open(sys.argv[2])); bad = 0; nt = 0
for key, trees in sorted(U.items()):
    m, y = key.split("_", 1)
    f = ROOT.TFile.Open("%s/mA_M%s/output_%s.root" % (X, m, y))
    for tn, nu in trees.items():
        t = f.Get(tn); nr = t.GetEntries() if t else -1
        for i in range(max(nr, 0)):
            if t.GetEntry(i) <= 0: nr = -2; break
        nt += 1
        if nr != nu: bad += 1; print("BAD mA%s %s %s: ROOT %d uproot %d" % (m, y, tn, nr, nu))
    f.Close()
print("trees checked %d, bad %d" % (nt, bad)); sys.exit(1 if bad else 0)
PY
rc=$?
grep -q "Error in <" $CH/logs/check_sig_mvacut.log && { echo "ROOT errors"; rc=1; }
tail -1 $CH/logs/check_sig_mvacut.log; exit $rc
