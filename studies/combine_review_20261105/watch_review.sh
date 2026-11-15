#!/usr/bin/env bash
# Watchdog for the Combine-review checks (bias rerun, HybridNew grid, local FitDiagnostics/r scan).
# Every 30 min: one status line to watch_status.txt. When a task finishes it is checked against its
# expected outputs and one line goes to PENDING_RESULTS.txt (DONE / DIED / INCOMPLETE). For the bias,
# a complete set of split files is merged (merge_bias_highstat.sh) and the injected value of every
# merged gaussfit.json is compared with the q50 recorded at submission. Numbers only, no interpretation.
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
BN=$FF/Combine/Checks/Bias_nominal
HN=$FF/Combine/Checks/review_20261105/hybridnew
REG=$W/PENDING_RESULTS.txt; STATUS=$W/watch_status.txt
BCL=$(cat $W/bias_cluster.txt); HCL=$(grep -oP 'cluster \K[0-9]+' $W/hybridnew_cluster.txt)
SUB=$(cat $W/bias_submitted_at.txt)
bias_done=0; hn_done=0; loc_done=0
reg() { echo "$(date '+%F %H:%M') $*" >> $REG; }

bias_finish() {
  local bad=0 m nexp ngot
  for m in $(seq 1 30); do
    nexp=$(( 10 * 2 * $(python3 -c "import json;print(len(json.load(open('$BN/bias_outputs_highstat_variants/mA_${m}_prePhidzee2_20261105/merged/BiasJson/${m}_gaussfit.json'))['fit_results']))") ))
    ngot=$(find $BN/bias_outputs_highstat/mA_$m/chunk*/BiasFits -name '*_split*_fits.root' -newermt "$SUB" -size +1k 2>/dev/null | wc -l)
    [ "$ngot" -eq "$nexp" ] || { echo "mA$m split files $ngot/$nexp" >> $W/logs/bias_check.txt; bad=1; }
  done
  if [ $bad -ne 0 ]; then reg "INCOMPLETE bias $BCL -> $W/logs/bias_check.txt"; return; fi
  for m in $(seq 1 30); do bash $FF/shellScripts/bias/Condor/merge_bias_highstat.sh $m >> $W/logs/bias_merge.log 2>&1 || echo "merge mA$m rc!=0" >> $W/logs/bias_check.txt; done
  ( BASE_DIR=$FF bash $FF/shellScripts/bias/Condor/collect_bias_results.sh ) > $W/logs/bias_collect.log 2>&1
  python3 - $BN $W > $W/bias_exp_vs_q50.txt <<'PY'
import json, sys
bn, w = sys.argv[1:3]
q = {int(l.split()[0]): float(l.split()[1]) for l in open(w + "/q50_at_bias_submission.txt")}
nbad = 0
for m in range(1, 31):
    d = json.load(open("%s/bias_outputs_highstat/mA_%d/merged/BiasJson/%d_gaussfit.json" % (bn, m, m)))
    ok = abs(d["exp"] - q[m]) < 1e-4
    nbad += not ok
    print("mA%d exp=%s q50=%s %s" % (m, d["exp"], q[m], "ok" if ok else "MISMATCH"))
print("mismatches", nbad)
PY
  reg "DONE bias $BCL -> merged in bias_outputs_highstat; $(tail -1 $W/bias_exp_vs_q50.txt) ($W/bias_exp_vs_q50.txt)"
}

for _ in $(seq 1 96); do   # up to 48 h
  b=$(condor_q $BCL -af JobStatus 2>/dev/null); bq=$(grep -c . <<<"$b"); bh=$(grep -c '^5$' <<<"$b")
  h=$(condor_q $HCL -af JobStatus 2>/dev/null); hq=$(grep -c . <<<"$h"); hh=$(grep -c '^5$' <<<"$h")
  nhn=$(ls $HN/mA*/higgsCombine_mA*.HybridNew.*.root 2>/dev/null | wc -l)
  nfd=$(ls $FF/Combine/Checks/review_20261105/fitdiag/*/fitDiagnostics_asimov.root 2>/dev/null | wc -l)
  nrs=$(ls $FF/Combine/Checks/review_20261105/rscan/*/higgsCombine_asimov.MultiDimFit.mH125.38.root 2>/dev/null | wc -l)
  alive=$(ps -u "$USER" -o cmd= --no-headers | grep -c '[r]un_fitdiag_rscan.sh' || true)
  echo "$(date '+%F %H:%M') bias q=$bq held=$bh | hybridnew q=$hq held=$hh files=$nhn/90 | fitdiag $nfd/60 rscan $nrs/30 alive=$alive" >> $STATUS
  [ "$bh" -gt 0 ] && reg "HELD bias $bh jobs in $BCL"
  [ "$hh" -gt 0 ] && reg "HELD hybridnew $hh jobs in $HCL"
  if [ $loc_done -eq 0 ]; then
    if [ -e $W/fitdiag_rscan.DONE ]; then reg "DONE fitdiag $nfd/60 rscan $nrs/30 -> $FF/Combine/Checks/review_20261105 (FAIL lines: $(cat $W/logs/fitdiag.log $W/logs/rscan.log | grep -c FAIL))"; loc_done=1
    elif [ "$alive" -eq 0 ]; then reg "DIED fitdiag/rscan after fitdiag $nfd/60 rscan $nrs/30"; loc_done=1; fi
  fi
  if [ $hn_done -eq 0 ] && [ "$hq" -eq 0 ]; then reg "DONE hybridnew $HCL files $nhn/90 -> $HN"; hn_done=1; fi
  if [ $bias_done -eq 0 ] && [ "$bq" -eq 0 ]; then bias_finish; bias_done=1; fi
  [ $bias_done -eq 1 ] && [ $hn_done -eq 1 ] && [ $loc_done -eq 1 ] && break
  sleep 1800
done
reg "WATCH END"
