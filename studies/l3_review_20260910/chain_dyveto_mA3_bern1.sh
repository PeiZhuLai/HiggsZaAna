#!/usr/bin/env bash
# 2026-10-01 23:0x: with the S1 scheme-2 lnN in the mA3 datacard (chain_dyveto_mA3_s1fix.sh) the mA3
# high-stat bias (10000 toys) gave Bern2 -0.211 +- 0.011 (before: -0.181 with the per-event variation;
# every member moved by -0.013..-0.031, i.e. a shift from the asymmetric lnN 0.968/1.014, not toy noise).
# Per AN Sec 7.4 (same rule as mA11 Bern3->2, mA20 Bern6->5, mA25 Bern6->5; user decision 2026-10-01
# "Bern follow the AN and lower the order"): use the next lower order, Bern1 (GOF p = 0.17 > 0.01).
# The best-fit function is Exp1, so the Asimov median limit should not move.
#   1. backups; fTest bestBernOrder(3): 3 -> 1; rebuild
#   2. background refit mA3; check the envelope has Bern1 and no Bern2/Bern3
#   3. t2w + limits (all 30, as 1_makeLimits.sh wipes the output dir); report any median that moved
#   4. condor: bias + impact mA3; merge; verdict
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FR=$F/Background/ALP_BkgModel_ReReco/fit_results_run3
FZ=$C/root_t2w_dyveto_20260930
BK=$F/Background
SRCF=$BK/test/fTest_ALP_turnOn.cpp
TAG=before_mA3_bern1_$(date +%Y%m%d_%H%M)
HB=$L/dyveto_mA3_bern1_heartbeat.txt
hb(){ echo "$(date '+%F %T') $*" | tee -a "$HB"; }
die(){ hb "STOPPING: $*"; exit 1; }
fresh(){ [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv

hb "[1] backups ($TAG), fTest bestBernOrder(3) 3 -> 1, rebuild"
V=$C/output_combine_results_variants; mkdir -p $V
find $C/output_combine_results -maxdepth 1 -type f ! -regex '.*/higgsCombine[0-9]+\.AsymptoticLimits\.mH125\.38\.root' -exec mv {} $V/ \;
cp -rp $C/output_combine_results $V/output_combine_results_$TAG
# variants live OUTSIDE ALP_BkgModel_ReReco: sync_figures.sh copies that whole tree into the AN
mkdir -p $F/Background/fit_results_variants && cp -rp $FR/3 $F/Background/fit_results_variants/3_$TAG
cp -p $SRCF $SRCF.$TAG
cp -p $FZ/3_Datacard_leptons.root $FZ/3_Datacard_leptons.root.$TAG
mkdir -p $BN/bias_outputs_highstat_variants
[ -d $BN/bias_outputs_highstat/mA_3 ] && mv $BN/bias_outputs_highstat/mA_3 $BN/bias_outputs_highstat_variants/mA_3_$TAG
cp -p $C/output_impacts/3_impacts.json $C/output_impacts/3_impacts.json.$TAG 2>/dev/null
python3 - $SRCF <<'PYS' || die "fTest source edit failed"
import sys
p = sys.argv[1]; s = open(p).read()
a = "    case 3: return 3; case 9: return 6;"
assert s.count(a) == 1, "case 3 line not found"
s = s.replace(a, "    // [PZ 2026-10-01] m_a = 3: with the S1 scheme-2 lnN Bern2 = -0.211 (10000 toys), the only function\n"
                 "    // above threshold -> next lower order (AN Sec 7.4); Bern1 GOF p = 0.17.\n"
                 "    case 3: return 1; case 9: return 6;")
open(p, "w").write(s); print("  fTest: bestBernOrder(3)=1")
PYS
( cd $BK && make -f makefile -B > $L/mA3bern1_build.log 2>&1 ) || die "background build failed (mA3bern1_build.log)"
[ "$(stat -c %Y $BK/bin/fTest_ALP_turnOn)" -ge "$(stat -c %Y $SRCF)" ] || die "fTest binary older than its source"

hb "[2] background refit mA3"
T=$(date +%s)
bash $D/fit_bkg_single/fit_bkg_m3.sh > $L/mA3bern1_bkg.log 2>&1
MP=$FR/3/CMS-HGG_mva_13p6TeV_multipdf.root
fresh $MP $T && [ $(stat -c %s $MP) -ge 1024 ] || die "mA3 multipdf missing/stale/truncated"
fresh $FR/3/HZAmassInde_fTest/EnvelopeResults.txt $T || die "mA3 EnvelopeResults stale"
sed 's/^/    /' $FR/3/HZAmassInde_fTest/EnvelopeResults.txt | tee -a $HB
grep -q "pdf : Bern1 " $FR/3/HZAmassInde_fTest/EnvelopeResults.txt || die "mA3 envelope has no Bern1"
grep -qE "pdf : Bern[2-9] " $FR/3/HZAmassInde_fTest/EnvelopeResults.txt && die "mA3 envelope still has Bern>=2"

hb "[3] t2w + limits (all 30)"
T=$(date +%s)
( cd $C && sh 2_text2ws.sh > $L/mA3bern1_t2w.log 2>&1 && sh 1_makeLimits.sh > $L/mA3bern1_limits.log 2>&1 )
fresh $C/root_t2w/3_Datacard_leptons.root $T || die "mA3 workspace stale"
for m in $(seq 1 30); do fresh $C/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "limit mA$m stale"; done
python3 - $C/output_combine_results $V/output_combine_results_$TAG <<'PY' | tee -a $HB
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def q(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m))
    v = [e.limit for e in f.Get("limit")]; f.Close(); return v
for m in range(1, 31):
    n, o = q(sys.argv[1], m), q(sys.argv[2], m)
    if any(abs(a / b - 1) > 1e-4 for a, b in zip(n, o)):
        print("  limit moved mA%d %s -> %s" % (m, " ".join("%.5f" % x for x in o), " ".join("%.5f" % x for x in n)))
print("  limit comparison done")
PY

hb "[4] condor: bias + impact mA3"
cp -p $C/root_t2w/3_Datacard_leptons.root $FZ/
T0=$(date +%s)
bout=$(ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh 3 2>&1)
iout=$(bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh 3 2>&1)
CL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$bout
$iout" | sort -u | tr '\n' ' ')
[ -n "${CL// /}" ] || die "submission failed"
hb "  clusters $CL"
while true; do
  q=$(condor_q $CL -af JobStatus 2>/dev/null) || { sleep 300; continue; }
  idle=$(grep -c '^1$' <<<"$q"); run=$(grep -c '^2$' <<<"$q"); held=$(grep -c '^5$' <<<"$q")
  while read -r id mem; do
    [ -z "${id:-}" ] && continue
    if [ "$mem" -lt 8000 ]; then new=8000; elif [ "$mem" -lt 16000 ]; then new=16000; elif [ "$mem" -lt 20000 ]; then new=20000; else continue; fi
    hb "  held $id mem $mem -> $new"; condor_qedit $id RequestMemory $new >/dev/null; condor_release $id >/dev/null
  done < <(condor_q $CL -constraint 'JobStatus==5 && regexp("emory", HoldReason)' -af:j RequestMemory 2>/dev/null)
  hb "  idle=$idle running=$run held=$held"
  [ $((idle + run + held)) -eq 0 ] && break
  [ $((idle + run)) -eq 0 ] && [ "$held" -gt 0 ] && { hb "  only held left"; break; }
  sleep 1800
done
ROOT_DATACARD_PATH=$FZ bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh 3 > $L/mA3bern1_bias_merge.log 2>&1
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh > $L/mA3bern1_bias_collect.log 2>&1
python3 - $BN/bias_outputs_highstat $C/output_impacts $T0 <<'PY' | tee -a $HB
import json, os, sys
B, I, T0 = sys.argv[1], sys.argv[2], float(sys.argv[3]); bad = 0
p = "%s/mA_3/merged/BiasJson/3_gaussfit.json" % B
if os.path.exists(p) and os.path.getmtime(p) >= T0:
    fr = json.load(open(p))["fit_results"]
    items = sorted(((k, v["mean"], v["mean_err"]) for k, v in fr.items()), key=lambda x: -abs(x[1]))
    print("  mA3 bias " + " ".join("%s:%+.3f+-%.3f" % x for x in items)); bad += any(abs(v) > 0.2 for _, v, _ in items)
else:
    print("  mA3 bias MISSING/stale"); bad += 1
ip = "%s/3_impacts.json" % I
if os.path.exists(ip) and os.path.getmtime(ip) >= T0:
    d = json.load(open(ip)); r = [q for q in d["POIs"] if q["name"] == "r"][0]["fit"][1]
    print("  mA3 impact r=%.3f, %d parameters" % (r, len(d["params"])))
else:
    print("  mA3 impact MISSING/stale"); bad += 1
print("VERDICT: %s" % ("PASSED" if bad == 0 else "CHECK"))
PY
echo "$(date '+%F %H:%M') DONE mA3 Bern2 -> Bern1 (heartbeat dyveto_mA3_bern1_heartbeat.txt)" >> $L/PENDING_RESULTS.txt
