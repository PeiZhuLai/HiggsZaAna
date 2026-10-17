#!/usr/bin/env bash
# Redo the data MVA-cut step (m_a 1-3) of chainlow_relE2w123. In the 15:00 run apply_bdt_data.py wrote
# mA_M3/output_2022postEE.root through the EOS FUSE mount and the file came out damaged: uproot reads 15 events,
# ROOT reads 2 with "TBufferFile::ReadClassBuffer ... TTree read too few bytes"; hadd returned 0 anyway and run3.root
# held 156 = 169 - 15 + 2 events. A re-merge through xrootd gives the same error, so the era file itself is bad.
# Here apply_bdt_data.py writes to a local AFS directory (-o; its hadd runs there too); every era file and run3.root is
# read entry by entry with ROOT, the counts must equal the direct count of the scored data at the working points
# (167 / 267 / 169), then the files are xrdcp'd to EOS and read again through xrootd. The old EOS data dir is moved to
# root_MVAcut/data.bad_fuse_20261017.
set -euo pipefail
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chainlow_relE2w123
S=$CH/flashggFinalFit; L=$CH/logs; LOC=$CH/mvacut_data_local
XB=/eos/user/p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2w123/root_MVAcut; XR=root://eosuser.cern.ch/
ERAS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh
eval "$(cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && scramv1 runtime -sh)"
export PYTHONPATH="$S/tools:$S/Signal/tools"; set -u
rm -rf $LOC; mkdir -p $LOC
( cd $S && python3 $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py -o $LOC ) > $L/mva_data_redo_local.log 2>&1
grep -q "Error in <" $L/mva_data_redo_local.log && { echo "ROOT errors in mva_data_redo_local.log"; exit 1; }
check() {   # $1 = base (local dir or xrootd URL)
python3 - "$1" "$ERAS" <<'PY'
import sys, ROOT
ROOT.gROOT.SetBatch(True)
B, eras = sys.argv[1], sys.argv[2].split(); want = {1: 167, 2: 267, 3: 169}; bad = 0
def pairs(p):
    f = ROOT.TFile.Open(p); t = f.Get("DiphotonTree/Data_13p6TeV"); s = []
    for i in range(t.GetEntries()):
        assert t.GetEntry(i) > 0, (p, i)
        s.append((round(t.CMS_hza_mass, 5), int(t.index)))
    f.Close(); return s
for m in (1, 2, 3):
    src = []
    for e in eras: src += pairs("%s/mA_M%d/output_%s.root" % (B, m, e))
    mer = pairs("%s/mA_M%d/run3.root" % (B, m))
    ok = len(src) == len(mer) == want[m] and sorted(src) == sorted(mer)
    print("  mA%d era sum %d, run3 %d, expected %d, pairs equal %s -> %s" % (m, len(src), len(mer), want[m], sorted(src) == sorted(mer), "OK" if ok else "BAD"))
    bad |= not ok
sys.exit(bad)
PY
}
check $LOC 2>&1 | tee $L/mva_data_redo_check_local.log
grep -q "Error in <" $L/mva_data_redo_check_local.log && { echo "ROOT errors reading the local files"; exit 1; }
xrdfs eosuser.cern.ch stat $XB/data.bad_fuse_20261017 >/dev/null 2>&1 || xrdfs eosuser.cern.ch mv $XB/data $XB/data.bad_fuse_20261017
for m in 1 2 3; do
  xrdfs eosuser.cern.ch mkdir -p $XB/data/mA_M$m
  for f in $LOC/mA_M$m/*.root; do xrdcp -f -s $f $XR$XB/data/mA_M$m/$(basename $f); done
done
check $XR$XB/data 2>&1 | tee $L/mva_data_redo_check_eos.log
grep -q "Error in <" $L/mva_data_redo_check_eos.log && { echo "ROOT errors reading the EOS copies"; exit 1; }
echo "OK"
