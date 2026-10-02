#!/usr/bin/env bash
#
# Monitor for the full m_a >= 2 GeV trigdenom run (systemd --user unit hza-v2trig-mAge2-mon).
# Every 30 min writes logs/heartbeat.txt: driver unit state, condor jobs of THIS run (matched by
# the eos_logs path in Cmd; the schedd is shared with other sessions), held jobs, per-sample/era
# job dirs vs non-empty .out files.
# Held for memory -> RequestMemory 8000 -> 16000 -> 20000 MB and release (never condor_rm);
# still held at 20000 -> written to logs/held_giveup.txt for manual follow-up.
# Own idle jobs are bumped to JobPrio 10 so they run ahead of this user's other idle jobs.
# When the queue holds no job of this run and the driver unit is gone: postprocess_mAge2.sh
# (collect -> plots -> plateau tables), then logs/MONITOR_DONE. Gives up after 72 h.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v2_trigdenom_mAge2_20261002
LOGS=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs/Sig_MC_trigdenomV2mAge2
HB=$D/logs/heartbeat.txt
CONSTR='regexp("eos_logs/Sig_MC_trigdenomV2mAge2/", Cmd)'
t0=$(date +%s)
while true; do
  now=$(date +%s)
  drv=$(systemctl --user is-active hza-v2trig-mAge2-full 2>/dev/null)
  q=$(timeout 120 condor_q -constraint "$CONSTR" -af JobStatus 2>/dev/null); qrc=$?
  nq=$(printf '%s\n' "$q" | grep -c . ); nidle=$(printf '%s\n' "$q" | grep -c '^1$'); nrun=$(printf '%s\n' "$q" | grep -c '^2$'); nheld=$(printf '%s\n' "$q" | grep -c '^5$')
  {
    echo "=== $(date '+%F %T') elapsed=$(( (now-t0)/60 ))min driver=$drv condor_q_rc=$qrc queue=$nq idle=$nidle run=$nrun held=$nheld"
    for sd in $LOGS/mA_M*; do
      [ -d "$sd" ] || continue
      nj=$(ls -d $sd/job_* 2>/dev/null | wc -l)
      np=$(find $sd -mindepth 2 -maxdepth 2 -name '*.out' -size +0 2>/dev/null | xargs -r -n1 dirname | sort -u | wc -l)
      echo "  $(basename $sd) job_dirs=$nj with_out=$np"
    done
  } > "$HB.tmp" && mv "$HB.tmp" "$HB"
  timeout 120 condor_q -constraint "$CONSTR && JobStatus==1 && JobPrio < 10" -af ClusterId ProcId 2>/dev/null | \
    while read -r c p; do condor_prio -p 10 "$c.$p" >/dev/null 2>&1; done
  if [ "$nheld" -gt 0 ]; then
    timeout 120 condor_q -constraint "$CONSTR && JobStatus==5" -af ClusterId ProcId RequestMemory HoldReason 2>/dev/null | while read -r c p mem reason; do
      id="$c.$p"
      echo "$(date '+%F %T') HELD $id mem=$mem reason=$reason" >> $D/logs/held.log
      if echo "$reason" | grep -qi "memory"; then
        if   [ "$mem" -lt 16000 ]; then new=16000
        elif [ "$mem" -lt 20000 ]; then new=20000
        else echo "$id mem=$mem $reason" >> $D/logs/held_giveup.txt; continue; fi
        echo "  condor_qedit $id RequestMemory $new ; condor_release $id" >> $D/logs/held.log
        condor_qedit "$id" RequestMemory $new >> $D/logs/held.log 2>&1
        condor_release "$id" >> $D/logs/held.log 2>&1
      fi
    done
  fi
  if [ $qrc -eq 0 ] && [ "$nq" -eq 0 ] && [ "$drv" != "active" ] && [ "$drv" != "activating" ]; then
    echo "$(date '+%F %T') queue empty and driver $drv -> post-processing" >> "$HB"
    RUN=Sig_MC_trigdenomV2mAge2 bash $D/postprocess_mAge2.sh >> $D/logs/postprocess_full.log 2>&1
    echo "postprocess rc=$?" >> "$HB"
    echo "MONITOR_DONE $(date '+%F %T')" > $D/logs/MONITOR_DONE
    exit 0
  fi
  if [ $(( now - t0 )) -gt 259200 ]; then echo "timeout" > $D/logs/MONITOR_TIMEOUT; exit 2; fi
  sleep 1800
done
