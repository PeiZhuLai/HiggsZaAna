#!/usr/bin/env bash
# Regenerate the MVA-cut DATA trees for selected mass points only.
#   usage: bash regen_mvacut_data.sh 17 [18 ...]
#
# Why not just rerun apply_bdt_data.py: it processes all 30 mass points, and its final step
# hadds output_<era>.root -> run3.root directly on EOS, logging (not raising) on failure.
# On 2026-09-24 that produced a mA17 run3.root which opened fine, reported 90 entries, and
# threw an I/O error at entry 41 -- caught only when Tree2WS iterated it.
#
# Here: the module is loaded with its global mAs overridden (process_files() reads the
# global, so a pass_map restriction alone is not enough), the output goes to AFS scratch,
# every tree is read end to end there, then copied to EOS and read end to end again.
set -uo pipefail
[ $# -ge 1 ] || { echo "usage: $0 <mA> [<mA> ...]"; exit 2; }
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
APPLY=$CMSSW_SRC/flashggFinalFit/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py
EOSD=/eos/home-p/pelai/HZa/root_MVAcut/data
SCR=$D/scratch_regen_mvacut_$(date +%Y%m%d_%H%M%S)
mkdir -p "$SCR"
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"
set -u

MASSES="$*"
echo "[regen] masses: $MASSES   scratch: $SCR"
python3 - "$APPLY" "$SCR" $MASSES <<'PY' || { echo "[regen] apply failed"; exit 1; }
import importlib.util, sys
path, out, masses = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3:]]
spec = importlib.util.spec_from_file_location("apply_bdt_data", path)
m = importlib.util.module_from_spec(spec)
sys.argv = [path, "-o", out]
spec.loader.exec_module(m)
m.mAs = masses
args = m.get_args()
m.setup_logging(args.log_level)
cuts = m.parse_mva_cuts(m.optimized_BDT_Cut)
samples = ["mA_M%d" % x for x in masses]
pass_map, id_cols = m.build_pass_event_map(samples, m.years, args.inputFolder, cuts)
m.process_files(args.outputFolder, args.inputFolder, pass_map, id_cols, cuts)
m.hadd_outputs(args.outputFolder, samples, m.years)
PY

fullread() {  # $1 = file ; prints OK or the failure
python3 - "$1" <<'PY'
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
p = sys.argv[1]
try:
    f = ROOT.TFile.Open(p)
except OSError as e:
    print("FAIL open", e); sys.exit(1)
t = f.Get("DiphotonTree/Data_13p6TeV")
if not t: print("FAIL no tree"); sys.exit(1)
n = t.GetEntries()
for i in range(n):
    if t.GetEntry(i) <= 0: print("FAIL entry %d/%d" % (i, n)); sys.exit(1)
print("OK %d entries" % n)
PY
}

rc=0
for mA in $MASSES; do
  for f in "$SCR/mA_M$mA"/output_*.root "$SCR/mA_M$mA/run3.root"; do
    r=$(fullread "$f"); echo "  scratch $(basename $(dirname $f))/$(basename $f): $r"
    [[ "$r" == OK* ]] || rc=1
  done
done
[ $rc -eq 0 ] || { echo "[regen] scratch copy failed a full read -- NOT copying to EOS"; exit 1; }

for mA in $MASSES; do
  mkdir -p "$EOSD/mA_M$mA"
  for f in "$SCR/mA_M$mA"/output_*.root "$SCR/mA_M$mA/run3.root"; do
    cp -f "$f" "$EOSD/mA_M$mA/"
  done
  # the old workspace was built from the corrupt tree; remove it so nothing downstream can
  # mistake it for current (Tree2WS rebuilds it)
  rm -f "$EOSD/mA_M$mA/ws/run3.root"
  r=$(fullread "$EOSD/mA_M$mA/run3.root"); echo "  EOS mA_M$mA/run3.root: $r"
  [[ "$r" == OK* ]] || rc=1
done
[ $rc -eq 0 ] && echo "[regen] VERDICT: PASSED" || echo "[regen] VERDICT: FAILED"
exit $rc
