#!/usr/bin/env bash
#
# V2 study monitor (systemd --user unit hza-v2trig-mon). Every 10 min writes logs/heartbeat.txt:
# driver unit state, condor jobs of THIS run (identified by the eos_logs path in Cmd, not by
# batch name -- the schedd is shared with other sessions), held jobs, and per-era completeness
# (job dirs vs chosen .out files that carry the nominal trigeff payloads).
# Held for memory -> qedit RequestMemory to 6000 once and release (logged).
# When no job of this run is left in the queue and the driver unit is gone, it runs the
# collector, the plots (OL and N2 denominators) and the plateau table, then writes
# logs/MONITOR_DONE. Gives up after 30 h (logs/MONITOR_TIMEOUT).
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v2_trigdenom_20260927
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
LOGS=$REPO/eos_logs/Sig_MC_trigdenomV2
HB=$D/logs/heartbeat.txt
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
PYROOT=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
CUT=$D/cutflow_list_trigdenomV2
PLOTS=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/trigEffCompareVlepPt_trigdenomV2
ERAS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
RELEASED=$D/logs/released_mem.txt
touch "$RELEASED"
t0=$(date +%s)
CONSTR='regexp("eos_logs/Sig_MC_trigdenomV2/", Cmd)'

while true; do
  now=$(date +%s)
  drv=$(systemctl --user is-active hza-v2trig-full 2>/dev/null)
  q=$(condor_q -constraint "$CONSTR" -af JobStatus 2>/dev/null)
  nq=$(printf '%s\n' "$q" | grep -c . ); nidle=$(printf '%s\n' "$q" | grep -c '^1$'); nrun=$(printf '%s\n' "$q" | grep -c '^2$'); nheld=$(printf '%s\n' "$q" | grep -c '^5$')
  {
    echo "=== $(date '+%F %T') elapsed=$(( (now-t0)/60 ))min driver=$drv queue=$nq idle=$nidle run=$nrun held=$nheld"
    for e in $ERAS; do
      sd=$LOGS/mA_M1_$e
      nj=$(ls -d $sd/job_* 2>/dev/null | wc -l)
      # log lines wrap, so the payload itself is verified by collect_trigdenom.py at the end;
      # here only count job dirs that already have a non-empty .out
      np=$(find $sd -mindepth 2 -maxdepth 2 -name '*.out' -size +0 2>/dev/null | xargs -r -n1 dirname | sort -u | wc -l)
      echo "  $e job_dirs=$nj jobs_with_nonempty_out=$np"
    done
  } > "$HB.tmp" && mv "$HB.tmp" "$HB"

  # The same user has ~2000 idle HZgamma jobs queued ahead; JobPrio only reorders our own
  # jobs, so bump the few jobs of this run (and any resubmission by the driver) to the front.
  condor_q -constraint "$CONSTR && JobStatus==1 && JobPrio < 10" -af ClusterId ProcId 2>/dev/null | \
    while read -r c p; do condor_prio -p 10 "$c.$p" >/dev/null 2>&1; done

  # held handling
  if [ "$nheld" -gt 0 ]; then
    condor_q -constraint "$CONSTR && JobStatus==5" -af ClusterId ProcId RequestMemory HoldReason 2>/dev/null | while read -r c p mem reason; do
      id="$c.$p"
      echo "$(date '+%F %T') HELD $id mem=$mem reason=$reason" >> $D/logs/held.log
      if echo "$reason" | grep -qi "memory" && ! grep -q "^$id$" "$RELEASED"; then
        echo "  condor_qedit $id RequestMemory 6000 ; condor_release $id" >> $D/logs/held.log
        condor_qedit "$id" RequestMemory 6000 >> $D/logs/held.log 2>&1
        condor_release "$id" >> $D/logs/held.log 2>&1
        echo "$id" >> "$RELEASED"
      fi
    done
  fi

  if [ "$nq" -eq 0 ] && [ "$drv" != "active" ] && [ "$drv" != "activating" ]; then
    echo "$(date '+%F %T') queue empty and driver $drv -> post-processing" >> "$HB"
    $PY $D/collect_trigdenom.py --logs $LOGS --out $CUT > $D/logs/collect_full.log 2>&1; crc=$?
    echo "collect rc=$crc" >> "$HB"
    if [ $crc -ne 0 ]; then
      echo "INCOMPLETE eras -> collecting with --allow-partial for inspection; NOT final" >> "$HB"
      $PY $D/collect_trigdenom.py --logs $LOGS --out ${CUT}_partial --allow-partial >> $D/logs/collect_full.log 2>&1
      echo "MONITOR_STOPPED_INCOMPLETE $(date '+%F %T')" > $D/logs/MONITOR_INCOMPLETE
      exit 1
    fi
    cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts
    for suf in OL N2; do
      env -i HOME=/afs/cern.ch/user/p/pelai PATH=/usr/bin:/bin PYTHONPATH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib \
        $PYROOT plot_trigEffVlepPt.py --ma 1 --in-dir $CUT --ord-suffix $suf --out $PLOTS/denom$suf > $D/logs/plot_$suf.log 2>&1
      echo "plot $suf rc=$?" >> "$HB"
    done
    env -i HOME=/afs/cern.ch/user/p/pelai PATH=/usr/bin:/bin PYTHONPATH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib \
      $PYROOT plot_trigEffVlepPt.py --ma 1 --in-dir $CUT --out $PLOTS/denomOld_sameRun > $D/logs/plot_old.log 2>&1
    echo "plot old-denom(same run) rc=$?" >> "$HB"
    for e in $ERAS; do
      env -i HOME=/afs/cern.ch/user/p/pelai PATH=/usr/bin:/bin $PYROOT $D/plot_denom_compare.py \
        --in-dir $CUT --era $e --out $PLOTS/denomCompare >> $D/logs/plot_compare.log 2>&1
      echo "plot compare $e rc=$?" >> "$HB"
    done
    $PY $D/plateau_table.py --dir $CUT --dir-old /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list \
      --json-out $D/plateau_table.json > $D/plateau_table.txt 2>&1
    echo "plateau rc=$?" >> "$HB"
    echo "MONITOR_DONE $(date '+%F %T')" > $D/logs/MONITOR_DONE
    exit 0
  fi
  if [ $(( now - t0 )) -gt 108000 ]; then echo "timeout" > $D/logs/MONITOR_TIMEOUT; exit 2; fi
  sleep 600
done
