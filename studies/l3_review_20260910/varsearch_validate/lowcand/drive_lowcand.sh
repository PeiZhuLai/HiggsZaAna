#!/usr/bin/env bash
# Detached driver of the validation of the two minimal-change low-mass candidates (user request 2026-10-14), all heavy
# work on condor (the local slots are taken by chain4):
#   1  step 1 (condor cluster given as $1): float64 retrain + incl/test-DY evaluation -> <OUT>/<name>/{model.pkl,summary.txt}
#   2  install the model as models/model_<name>.{pkl,meta.json} (rescore_selected.py --low-model reads it there)
#   3  data + DY re-scoring with the candidate (6 sample dirs, 18 files) -> sandbox_varsearch_v1/run3_bdt_scored_<name>
#   4  cuts for 2, 1.5, 1, 0.7, 0.5 x the production data-sideband count (mA1-3) -> eval/closure_cuts_<name>.txt
#   5  closure scans (DROP=1, production signal model): 1 smoke job first, then the other 29 (both candidates)
#   6  mA1 bootstrap: 30 replicas at x1 and x0.5 per candidate (120 jobs), production mA1 datacard
#   7  summaries -> lowcand/summary_closure.txt
# Status: lowcand/logs/drive.status; DONE / STOP lines in ../logs/PENDING_RESULTS.txt (numbers only).
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
L=$V/lowcand; REG=$V/logs/PENDING_RESULTS.txt; ST=$L/logs/drive.status
OUT=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/lowcand; SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
DC1=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons.txt
NAMES="low_swapPt_relE2 low_etaW_relE1"
SAMPLES="Data DYGto2LG_10to100 DYJetsTo2E DYJetsTo2Mu DYJetsTo2Tau DYJetsToLL"
stat() { echo "$(date '+%F %T') $*" > $ST; echo "$(date '+%F %T') $*" >> $L/logs/drive.log; }
stop() { stat "STOP $*"; echo "$(date '+%F %H:%M') STOP  lowcand validation  $*" >> $REG; exit 1; }
waitc() { while [ "$(condor_q $1 -af ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
            h=$(condor_q $1 -af JobStatus 2>/dev/null | grep -c '^5$'); stat "$2: cluster $1 in queue $(condor_q $1 -af ProcId 2>/dev/null | wc -l) held $h"
            [ "$h" -gt 0 ] && [ "$(condor_q $1 -af JobStatus 2>/dev/null | grep -vc '^5$')" -eq 0 ] && stop "$2: all remaining jobs of cluster $1 held"
            sleep 120; done; }
xp() { echo "$1" | sed 's#^/eos/home-p/pelai/#/eos/user/p/pelai/#'; }
xcat() { xrdcp -f -s root://eosuser.cern.ch/$(xp $1) - 2>/dev/null; }   # read through xrootd (local FUSE can be stale)
sub() { (cd $L && condor_submit $1 LIST=$2 2>&1 | grep -oP 'submitted to cluster \K[0-9]+'); }

# 1
waitc $1 "step1"
# existence checks and copies through xrootd: the local EOS FUSE showed a freshly written dir as d????????? (10-14)
XOUT=/eos/user/p/pelai/HZa/sandbox_varsearch_v1/lowcand
for n in $NAMES; do for f in model.pkl summary.txt; do
  xrdfs eosuser.cern.ch stat $XOUT/$n/$f > /dev/null 2>&1 || stop "step1 output $f missing for $n ($L/logs/step1.*)"; done; done
# 2
for n in $NAMES; do
  xrdcp -f -s root://eosuser.cern.ch/$XOUT/$n/model.pkl $V/models/model_$n.pkl && xrdcp -f -s root://eosuser.cern.ch/$XOUT/$n/model.meta.json $V/models/model_$n.meta.json || stop "install $n"
  [ "$(xrdadler32 $V/models/model_$n.pkl | awk '{print $1}')" = "$(xrdfs eosuser.cern.ch query checksum $XOUT/$n/model.pkl | awk '{print $2}')" ] || stop "install $n: checksum differs"
done
# 3
for n in $NAMES; do for s in $SAMPLES; do echo "$s $SBE/run3_bdt_scored_$n $n"; done; done > $L/list_rescore.txt
id=$(sub rescore.sub list_rescore.txt); [ -n "$id" ] || stop "rescore submission failed"
waitc $id "rescore"
for n in $NAMES; do
  lst=$(xrdfs eosuser.cern.ch ls -R $(xp $SBE/run3_bdt_scored_$n) 2>/dev/null); k=$(echo "$lst" | grep -c '\.root$'); p=$(echo "$lst" | grep -c '\.partial$')
  [ "$k" -eq 18 ] && [ "$p" -eq 0 ] || stop "rescore $n: $k/18 files, $p partial ($L/logs/rescore.$id.*)"
done
# 4
for n in $NAMES; do
  env -i HOME=$HOME PATH=/usr/bin:/bin $PY $V/cuts_for_counts.py $SBE/run3_bdt_scored_$n 1,2,3 2,1.5,1,0.7,0.5 > $V/eval/closure_cuts_$n.txt || stop "cuts $n"
  [ $(wc -l < $V/eval/closure_cuts_$n.txt) -eq 3 ] || stop "cuts $n: not 3 lines"
done
# 5
for n in $NAMES; do while read m c1 c2 c3 c4 c5; do for c in $c1 $c2 $c3 $c4 $c5; do echo "scan cand_$n $m $c"; done; done < $V/eval/closure_cuts_$n.txt; done > $L/list_scan_all.txt
head -1 $L/list_scan_all.txt > $L/list_scan_smoke.txt; tail -n +2 $L/list_scan_all.txt > $L/list_scan_rest.txt
id=$(sub closure.sub list_scan_smoke.txt); [ -n "$id" ] || stop "smoke submission failed"
waitc $id "closure smoke"
read -r _ sm sma smc < $L/list_scan_smoke.txt
R=$SBE/closure_scan/${sm}_drop_mA$sma/cut_$(echo $smc | tr -d .)/result.txt
xcat $R | grep -q "dataN=[0-9].*closureZ=[0-9]" || stop "closure smoke job did not produce a result ($R, $L/logs/closure.$id.*)"
id=$(sub closure.sub list_scan_rest.txt); [ -n "$id" ] || stop "scan submission failed"
waitc $id "closure scans"
for n in $NAMES; do for m in 1 2 3; do
  D=$SBE/closure_scan/cand_${n}_drop_mA$m
  for c in $(awk -v m=$m '$1==m {$1=""; print}' $V/eval/closure_cuts_$n.txt); do x=$(xcat $D/cut_$(echo $c | tr -d .)/result.txt); [ -n "$x" ] && echo "$x" || echo "$c MISSING"; done > $L/scan_${n}_mA$m.txt
  grep -q MISSING $L/scan_${n}_mA$m.txt && stop "closure scan $n mA$m incomplete"
done; done
# 6
for n in $NAMES; do for c in $(awk '$1==1 {print $4, $6}' $V/eval/closure_cuts_$n.txt); do
  for s in $(seq 1 30); do echo "boot cand_$n 1 $c $s $DC1"; done; done; done > $L/list_boot.txt
id=$(sub closure.sub list_boot.txt); [ -n "$id" ] || stop "bootstrap submission failed"
waitc $id "bootstrap"
# 7
{ for n in $NAMES; do
    echo "===== $n"; cat $OUT/$n/summary.txt
    for m in 1 2 3; do echo "-- closure scan mA$m (DROP=1, production signal model)"; cut -c1-200 $L/scan_${n}_mA$m.txt; done
    for c in $(awk '$1==1 {print $4, $6}' $V/eval/closure_cuts_$n.txt); do
      echo "-- mA1 bootstrap cut $c"
      B=$SBE/closure_boot/cand_${n}_drop_mA1_cut$(echo $c | tr -d .)
      for s in $(seq 1 30); do xcat $B/rep_$s/result.txt; done > $L/boot_${n}_$(echo $c | tr -d .).txt
      $PY -c "
import numpy as np
rows=[l.split() for l in open('$L/boot_${n}_$(echo $c | tr -d .).txt') if 'rhat=' in l]
r=np.array([float(x[3][5:]) for x in rows if x[3][5:] not in ('NA','')]); z=np.array([float(x[2][2:]) for x in rows if x[2][2:] not in ('NA','')])
print('replicas %d/30  rhat median %+.3f mean %+.3f sd %.3f  frac(rhat>0) %.2f  |  Z mean %.2f' % (r.size, np.median(r), r.mean(), r.std(ddof=1), (r>0).mean(), z.mean()))"
    done
  done; } > $L/summary_closure.txt 2>&1
stat "DONE"
echo "$(date '+%F %H:%M') DONE  lowcand validation (step 1 + closure + mA1 bootstrap)  -> $L/summary_closure.txt" >> $REG
