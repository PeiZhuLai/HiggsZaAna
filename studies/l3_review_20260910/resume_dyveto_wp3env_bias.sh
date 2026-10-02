#!/usr/bin/env bash
# 2024 DY+jets overlap-veto retrain, user decisions 2026-10-01:
#  * mA3 working point 0.992 -> 0.988: after the retrain the closure at 0.992 is Z = 1.59 (bootstrap
#    40 replicas: Z 1.62 +- 1.17, 42% above 2); 0.988 is the most stable (bootstrap Z 1.03 +- 1.07,
#    15% above 2; data-envelope GOF >= 0.92; R 1.16, lowest).
#  * envelope per AN Sec 7.4 after the high-stat bias (> 0.2): mA25 Bernstein capped 6 -> 5
#    (Bern6 +0.212), Laurent removed at mA24 (Lau1 -0.232) and mA25 (Lau1 -0.215) -- Lau1 is the
#    lowest order, as for mA8 Pow1 on 2026-09-25.
# Adapted from chain_wp_mA3_0992.sh: mA3 MVA cut -> ws -> bkg fit; bkg refit of mA24/25 with the
# rebuilt fTest; signal model mA3; all datacards + limits; bias + impact for every mass whose limit
# moved (must include 3, 24, 25); then wait, merge and summarize those biases.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FR=$F/Background/ALP_BkgModel_ReReco/fit_results_run3
FZ=$C/root_t2w_dyveto_20260930
BK=$F/Background
SRCF=$BK/test/fTest_ALP_turnOn.cpp
HB=$L/dyveto_wp3_env2425_heartbeat.txt
EOSMVA=/eos/home-p/pelai/HZa/root_MVAcut
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
NEWCUT=0.988
TAG=before_dyveto_wp0988_env2425_$(date +%Y%m%d_%H%M)
S1=$D/wp3_single; mkdir -p $S1
die() { echo "STOPPING: $*"; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

# 2026-10-01 15:0x resume: the chain stopped at its own check "mA24 limit did not change". That
# check was too strict: removing a non-best envelope member (Lau1) leaves the Asimov expected limit
# unchanged; the new t2w workspaces were verified to hold the new envelopes (mA24 Exp1/Pow1/Bern5,
# mA25 Pow1/Exp1/Bern5). Bias + impact for 3 24 25 from here.
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
TAG=before_dyveto_wp0988_env2425_20261001_1330
mkdir -p $BN/bias_outputs_highstat_variants
CHANGED="3 24 25"
echo "=============== [9] condor: bias + impact for changed masses: $CHANGED ==============="
for m in $CHANGED; do
  cp -p $FZ/${m}_Datacard_leptons.root $FZ/${m}_Datacard_leptons.root.$TAG 2>/dev/null
  cp -p $C/root_t2w/${m}_Datacard_leptons.root $FZ/
  [ -d $BN/bias_outputs_highstat/mA_$m ] && mv $BN/bias_outputs_highstat/mA_$m $BN/bias_outputs_highstat_variants/mA_${m}_$TAG
done
T0=$(date +%s)
bout=$(ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $CHANGED 2>&1)
iout=$(bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $CHANGED 2>&1)
CL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$bout
$iout" | sort -u | tr '\n' ' ')
[ -n "${CL// /}" ] || die "bias/impact submission failed"
echo "SUBMITTED $(date '+%F %T') clusters $CL" | tee -a $HB

echo "=============== [10] wait, merge, summarize bias for $CHANGED ==============="
while true; do
  q=$(condor_q $CL -af JobStatus 2>/dev/null) || { sleep 300; continue; }
  idle=$(grep -c '^1$' <<<"$q"); run=$(grep -c '^2$' <<<"$q"); held=$(grep -c '^5$' <<<"$q")
  while read -r id mem; do
    [ -z "${id:-}" ] && continue
    if [ "$mem" -lt 8000 ]; then new=8000; elif [ "$mem" -lt 16000 ]; then new=16000; elif [ "$mem" -lt 20000 ]; then new=20000; else continue; fi
    echo "$(date '+%F %T') held $id mem $mem -> $new" >> $HB; condor_qedit $id RequestMemory $new >/dev/null; condor_release $id >/dev/null
  done < <(condor_q $CL -constraint 'JobStatus==5 && regexp("emory", HoldReason)' -af:j RequestMemory 2>/dev/null)
  echo "$(date '+%F %T') idle=$idle running=$run held=$held" >> $HB
  [ $((idle + run + held)) -eq 0 ] && break
  [ $((idle + run)) -eq 0 ] && [ "$held" -gt 0 ] && { echo "only held left" >> $HB; break; }
  sleep 1800
done
ROOT_DATACARD_PATH=$FZ bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh $CHANGED > $L/wp3env_bias_merge.log 2>&1
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh > $L/wp3env_bias_collect.log 2>&1
python3 - $BN/bias_outputs_highstat $C/output_impacts $T0 $CHANGED <<'PYB' | tee -a $HB
import json, os, sys
B, I, T0, ms = sys.argv[1], sys.argv[2], float(sys.argv[3]), [int(x) for x in sys.argv[4:]]
bad = 0
for m in ms:
    p = "%s/mA_%d/merged/BiasJson/%d_gaussfit.json" % (B, m, m)
    ip = "%s/%d_impacts.json" % (I, m)
    if not os.path.exists(p): print("mA%d bias MISSING" % m); bad += 1; continue
    fr = json.load(open(p))["fit_results"]
    items = sorted(((k, v["mean"]) for k, v in fr.items()), key=lambda x: -abs(x[1]))
    r = None
    if os.path.exists(ip) and os.path.getmtime(ip) >= T0:
        r = [q for q in json.load(open(ip))["POIs"] if q["name"] == "r"][0]["fit"][1]
    over = [k for k, v in items if abs(v) > 0.2]
    bad += bool(over) + (r is None)
    print("mA%-2d bias %s | impact r=%s%s" % (m, " ".join("%s:%+.3f" % x for x in items),
          "%.3f" % r if r is not None else "MISSING", "  ABOVE 0.2: %s" % over if over else ""))
print("VERDICT: %s" % ("PASSED" if bad == 0 else "CHECK (%d issues)" % bad))
PYB
echo "$(date '+%F %H:%M') DONE dyveto mA3 0.988 + env mA24/25 (log dyveto_wp3_env2425.log, heartbeat dyveto_wp3_env2425_heartbeat.txt)" >> $L/PENDING_RESULTS.txt
