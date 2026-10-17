#!/usr/bin/env bash
# Detached driver for the validation of a set of candidates (user 2026-10-15: high-mass candidates at mA25/30, low-mass
# subset candidates), generalization of drive_lowcand.sh. All heavy work on condor; existence checks through xrootd.
# Spec file: one candidate per line "<name> <low|high> <feature> ...".
#   1 step 1 (validate_cand.py: float64 retrain + incl/test-DY evaluation)       -> <OUT>/<name>/
#   2 install models/model_<name>.{pkl,meta.json} (adler32 checked)
#   3 data + DY re-scoring (6 sample dirs, 18 files; --low-model or --high-model) -> run3_bdt_scored_<name>
#   4 cuts for 2, 1.5, 1, 0.7, 0.5 x production data-sideband count             -> eval/closure_cuts_<name>.txt
#      low: mA1, 2, 3   high: mA25, 30
#   5 closure scans (DROP=1, production signal model and datacard), smoke job first
#   6 bootstrap, 30 replicas: low mA1 at x1 and x0.5; high mA25 and mA30 at x1 (production-like data count)
#   7 summary -> lowcand/summary_<tag>.txt
# Usage: drive_cands.sh <spec file> <tag>
SPEC=$1; TAG=$2
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
L=$V/lowcand; REG=$V/logs/PENDING_RESULTS.txt; ST=$L/logs/drive_$TAG.status
OUT=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/lowcand; SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
DCD=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons
SAMPLES="Data DYGto2LG_10to100 DYJetsTo2E DYJetsTo2Mu DYJetsTo2Tau DYJetsToLL"
NAMES=$(awk '{print $1}' $SPEC)
model_of() { awk -v n=$1 '$1==n {sub(/123$/, "", $2); print $2}' $SPEC; }   # low123 (10-15: window training at m_a 1-3) -> low
w23() { [ "$(awk -v n=$1 '$1==n {print $2}' $SPEC)" = low123 ]; }   # low123 also bootstraps m_a 2, 3 at x0.5
stat() { echo "$(date '+%F %T') $*" > $ST; echo "$(date '+%F %T') $*" >> $L/logs/drive_$TAG.log; }
stop() { stat "STOP $*"; echo "$(date '+%F %H:%M') STOP  candidates $TAG  $*" >> $REG; exit 1; }
xp() { echo "$1" | sed 's#^/eos/home-p/pelai/#/eos/user/p/pelai/#'; }
xcat() { xrdcp -f -s root://eosuser.cern.ch/$(xp $1) - 2>/dev/null; }
waitc() { while [ "$(condor_q $1 -af ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
            h=$(condor_q $1 -af JobStatus 2>/dev/null | grep -c '^5$'); stat "$2: cluster $1 in queue $(condor_q $1 -af ProcId 2>/dev/null | wc -l) held $h"
            [ "$h" -gt 0 ] && [ "$(condor_q $1 -af JobStatus 2>/dev/null | grep -vc '^5$')" -eq 0 ] && stop "$2: all remaining jobs of cluster $1 held"
            sleep 120; done; }
sub() { (cd $L && condor_submit $1 LIST=$2 2>&1 | grep -oP 'submitted to cluster \K[0-9]+'); }
masses_of() { [ "$(model_of $1)" = low ] && echo "1 2 3" || echo "25 30"; }

# 1
while read -r n mo feats; do echo "$L/validate_cand.py $mo $n $OUT $feats"; done < $SPEC > $L/args_step1_$TAG.txt
sed -e "s#queue ARGS from $L/args_step1.txt#queue ARGS from $L/args_step1_$TAG.txt#" -e "s#request_memory  = 8000MB#request_memory  = 12000MB#" $L/step1.sub > $L/step1_$TAG.sub
grep -q "args_step1_$TAG.txt" $L/step1_$TAG.sub || stop "step1 submit file not rewritten"
id=$(cd $L && condor_submit step1_$TAG.sub 2>&1 | grep -oP 'submitted to cluster \K[0-9]+')
[ -n "$id" ] || stop "step1 submission failed"
waitc $id "step1"
for n in $NAMES; do for f in model.pkl summary.txt; do
  xrdfs eosuser.cern.ch stat $(xp $OUT/$n/$f) > /dev/null 2>&1 || stop "step1 output $f missing for $n ($L/logs/step1.$id.*)"; done; done
# 2
for n in $NAMES; do
  xrdcp -f -s root://eosuser.cern.ch/$(xp $OUT/$n/model.pkl) $V/models/model_$n.pkl && xrdcp -f -s root://eosuser.cern.ch/$(xp $OUT/$n/model.meta.json) $V/models/model_$n.meta.json || stop "install $n"
  [ "$(xrdadler32 $V/models/model_$n.pkl | awk '{print $1}')" = "$(xrdfs eosuser.cern.ch query checksum $(xp $OUT/$n/model.pkl) | awk '{print $2}')" ] || stop "install $n: checksum differs"
done
# 3
for n in $NAMES; do for s in $SAMPLES; do
  if [ "$(model_of $n)" = low ]; then echo "$s $SBE/run3_bdt_scored_$n $n high"; else echo "$s $SBE/run3_bdt_scored_$n low $n"; fi
done; done > $L/list_rescore_$TAG.txt
id=$(sub rescore.sub list_rescore_$TAG.txt); [ -n "$id" ] || stop "rescore submission failed"
waitc $id "rescore"
for n in $NAMES; do
  lst=$(xrdfs eosuser.cern.ch ls -R $(xp $SBE/run3_bdt_scored_$n) 2>/dev/null); k=$(echo "$lst" | grep -c '\.root$'); p=$(echo "$lst" | grep -c '\.partial$')
  [ "$k" -eq 18 ] && [ "$p" -eq 0 ] || stop "rescore $n: $k/18 files, $p partial ($L/logs/rescore.$id.*)"
done
# 4
for n in $NAMES; do
  ms=$(masses_of $n | tr ' ' ',')
  env -i HOME=$HOME PATH=/usr/bin:/bin $PY $V/cuts_for_counts.py $SBE/run3_bdt_scored_$n $ms 2,1.5,1,0.7,0.5 > $V/eval/closure_cuts_$n.txt || stop "cuts $n"
  [ $(wc -l < $V/eval/closure_cuts_$n.txt) -eq $(masses_of $n | wc -w) ] || stop "cuts $n: wrong line count"
done
# 5
for n in $NAMES; do while read m c1 c2 c3 c4 c5; do for c in $c1 $c2 $c3 $c4 $c5; do echo "scan cand_$n $m $c"; done; done < $V/eval/closure_cuts_$n.txt; done > $L/list_scan_all_$TAG.txt
head -1 $L/list_scan_all_$TAG.txt > $L/list_scan_smoke_$TAG.txt; tail -n +2 $L/list_scan_all_$TAG.txt > $L/list_scan_rest_$TAG.txt
id=$(sub closure.sub list_scan_smoke_$TAG.txt); [ -n "$id" ] || stop "smoke submission failed"
waitc $id "closure smoke"
read -r _ sm sma smc < $L/list_scan_smoke_$TAG.txt
xcat $SBE/closure_scan/${sm}_drop_mA$sma/cut_$(echo $smc | tr -d .)/result.txt | grep -q "dataN=[0-9].*closureZ=[0-9]" || stop "closure smoke job did not produce a result"
id=$(sub closure.sub list_scan_rest_$TAG.txt); [ -n "$id" ] || stop "scan submission failed"
waitc $id "closure scans"
for n in $NAMES; do for m in $(masses_of $n); do
  D=$SBE/closure_scan/cand_${n}_drop_mA$m
  for c in $(awk -v m=$m '$1==m {$1=""; print}' $V/eval/closure_cuts_$n.txt); do x=$(xcat $D/cut_$(echo $c | tr -d .)/result.txt); [ -n "$x" ] && echo "$x" || echo "$c MISSING"; done > $L/scan_${n}_mA$m.txt
  grep -q MISSING $L/scan_${n}_mA$m.txt && stop "closure scan $n mA$m incomplete"
done; done
# 6
for n in $NAMES; do
  if w23 $n; then pts=$(awk '$1==1 {print "1:"$4, "1:"$6} $1==2 || $1==3 {print $1":"$6}' $V/eval/closure_cuts_$n.txt)
  elif [ "$(model_of $n)" = low ]; then pts=$(awk '$1==1 {print "1:"$4, "1:"$6}' $V/eval/closure_cuts_$n.txt)
  else pts=$(awk '{print $1":"$4}' $V/eval/closure_cuts_$n.txt); fi
  for pt in $pts; do m=${pt%%:*}; c=${pt##*:}; for s in $(seq 1 30); do echo "boot cand_$n $m $c $s $DCD/${m}_pruned_datacard_leptons.txt"; done; done
done > $L/list_boot_$TAG.txt
id=$(sub closure.sub list_boot_$TAG.txt); [ -n "$id" ] || stop "bootstrap submission failed"
waitc $id "bootstrap"
# 7
{ for n in $NAMES; do
    echo "===== $n ($(model_of $n))"; xcat $OUT/$n/summary.txt
    for m in $(masses_of $n); do echo "-- closure scan mA$m (DROP=1, production signal model)"; cut -c1-200 $L/scan_${n}_mA$m.txt; done
    awk '{print $1, ($1<=3 ? $4" "$6 : $4)}' $V/eval/closure_cuts_$n.txt | while read m cs; do for c in $cs; do
      if w23 $n; then [ "$m" != 1 ] && [ "$c" != "$(awk -v m=$m '$1==m {print $6}' $V/eval/closure_cuts_$n.txt)" ] && continue
      else [ "$(model_of $n)" = low ] && [ "$m" != 1 ] && continue; fi
      B=$SBE/closure_boot/cand_${n}_drop_mA${m}_cut$(echo $c | tr -d .)
      for s in $(seq 1 30); do xcat $B/rep_$s/result.txt; done > $L/boot_${n}_mA${m}_$(echo $c | tr -d .).txt
      echo "-- mA$m bootstrap cut $c: $($PY -c "
import numpy as np
rows=[l.split() for l in open('$L/boot_${n}_mA${m}_$(echo $c | tr -d .).txt') if 'rhat=' in l]
r=np.array([float(x[3][5:]) for x in rows if x[3][5:] not in ('NA','')]); z=np.array([float(x[2][2:]) for x in rows if x[2][2:] not in ('NA','')])
print('replicas %d/30  rhat median %+.4f sd %.4f  frac(rhat>0) %.2f  |  Z mean %.2f' % (r.size, np.median(r), r.std(ddof=1), (r>0).mean(), z.mean()))")"
    done; done
  done; } > $L/summary_$TAG.txt 2>&1
stat "DONE"
echo "$(date '+%F %H:%M') DONE  candidates $TAG (step 1 + closure + bootstrap)  -> $L/summary_$TAG.txt" >> $REG
