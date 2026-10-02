#!/usr/bin/env bash
#
# Persistent orchestrator for the V3/V8 extra-background study (runs as the systemd --user
# unit hza-xbkg-orch; lxplus el9 kills setsid/nohup processes on logout).
#
#   A. wait for the small-batch TEST driver (unit hza-xbkg-test, 8 jobs = 2 per sample)
#   B. TEST post-chain: completeness, chunk==merged, p2root+scoring, root==merged, yields
#   C. GATE: every TEST job ran (summary present), chunk==merged, merged==root, scores finite
#      -> otherwise STOP here and report (nothing big is submitted)
#   D. FULL driver: 381 jobs (TTto2L2Nu 147 [10% subset], TTG 143, TTGG_Run3 62, ZGG 29)
#   E. FULL post-chain -> yields_extrabkg.{md,json}
# Throughout: held-job handler (memory holds: 3000->8000->16000->20000 + release; other
# holds: release at most 3 times per job, then leave held and report), heartbeat every 5 min.
#
# Heartbeat: $D/logs/heartbeat_orch.txt     Done marker: $D/logs/orch.done
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v3v8_extrabkg_20260927
HB=$D/logs/heartbeat_orch.txt
HOLDLOG=$D/logs/held_handler.log
RELCNT=$D/logs/held_release_counts.txt
touch "$RELCNT"
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1; set -u
hb() { echo "$(date '+%F %T') $*" >> "$HB"; }
# Identify my jobs by the executable path only (schedd is shared with other projects).
CONST='Owner=="pelai" && regexp("eos_logs/Bkg_MC_extraBkg",Cmd)'

next_rung() {
  if [ "$1" -lt 8000 ]; then echo 8000; elif [ "$1" -lt 16000 ]; then echo 16000
  elif [ "$1" -lt 20000 ]; then echo 20000; else echo stuck; fi; }

held_pass() {
  condor_q -constraint "$CONST && JobStatus==5" -af ClusterId ProcId HoldReasonCode RequestMemory 2>/dev/null |
  while read -r c p code mem; do
    id="$c.$p"
    reason=$(condor_q "$id" -af HoldReason 2>/dev/null | head -1)
    if [ "$code" = "34" ]; then
      nm=$(next_rung "${mem%.*}")
      if [ "$nm" = stuck ]; then echo "$(date '+%F %T') $id memory hold at ${mem} MB -> left held: $reason" >> "$HOLDLOG"; continue; fi
      echo "$(date '+%F %T') $id memory hold ${mem}->${nm} MB: condor_qedit $id RequestMemory $nm; condor_release $id | $reason" >> "$HOLDLOG"
      condor_qedit "$id" RequestMemory "$nm" >> "$HOLDLOG" 2>&1; condor_release "$id" >> "$HOLDLOG" 2>&1
    else
      n=$(grep -c "^$id\$" "$RELCNT")
      if [ "$n" -lt 3 ]; then
        echo "$id" >> "$RELCNT"
        echo "$(date '+%F %T') $id hold code=$code (release $((n+1))/3): condor_release $id | $reason" >> "$HOLDLOG"
        condor_release "$id" >> "$HOLDLOG" 2>&1
      else
        grep -q "LEFT_HELD $id" "$HOLDLOG" || echo "$(date '+%F %T') LEFT_HELD $id after 3 releases, code=$code | $reason" >> "$HOLDLOG"
      fi
    fi
  done; }

# The user's other sessions keep ~2000 idle jobs on this schedd with an older QDate and the
# user priority is deep in fair-share debt (~40 starts/h), so at JobPrio 0 these 8+381 short
# jobs would queue behind all of them. Raise only these jobs (own-user ordering, nobody
# else's jobs are affected).
prio_pass() {
  ids=$(condor_q -constraint "$CONST && JobStatus==1 && JobPrio<5" -af ClusterId 2>/dev/null | tr '\n' ' ')
  [ -n "$ids" ] && condor_prio +5 $ids >/dev/null 2>&1 && echo "$(date '+%F %T') condor_prio +5 on $(wc -w <<< "$ids") idle jobs" >> "$HOLDLOG"
  return 0; }

watch_until() {   # $1 = marker file to wait for; runs held handler + heartbeat
  local label=$1 marker=$2
  while [ ! -f "$marker" ]; do
    held_pass
    prio_pass
    q=$(condor_q -constraint "$CONST" -af JobStatus 2>/dev/null | sort | uniq -c | awk '{printf "%s:%s ",$2,$1}')
    hb "[$label] waiting for $(basename $marker); my jobs by JobStatus(1=idle,2=run,5=held): ${q:-none}"
    sleep 300
  done
  hb "[$label] $(basename $marker): $(cat $marker)"; }

hb "orchestrator start (pid $$, host $(hostname))"

# A + B
MODE=test bash $D/chain_post.sh > $D/logs/chain_test.log 2>&1 &
pid_ct=$!
watch_until test-chain "$D/logs/chain_test.done"
wait $pid_ct

# C. gate
gate_ok=1
grep -q "rc12=0 rc4=0 rc5=0" $D/logs/chain_test.done || gate_ok=0
grep -q "FAIL" $D/logs/reconcile_test.txt $D/logs/stage3_test.txt && gate_ok=0
nerr=$(grep -l "vector_read\|Traceback\|ENV_FETCH_FAIL\|ENV_EXTRACT_FAIL" /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs/Bkg_MC_extraBkgTEST/*/job_*/*.err 2>/dev/null | wc -l)
hb "[gate] test chain: $(cat $D/logs/chain_test.done); job .err files with errors: $nerr; gate_ok=$gate_ok"
if [ "$gate_ok" != 1 ]; then
  hb "STOP: small-batch test gate failed -- full production NOT submitted"
  echo "STOPPED_AT_TEST_GATE" > $D/logs/orch.done; exit 1
fi

# D + E
hb "[full] launching driver: 381 jobs expected (TTto2L2Nu 147, TTG 143, TTGG_Run3 62, ZGG 29)"
MODE=full bash $D/run_extrabkg_driver.sh > $D/logs/driver_full.log 2>&1 &
pid_df=$!
MODE=full bash $D/chain_post.sh > $D/logs/chain_full.log 2>&1 &
pid_cf=$!
watch_until full-chain "$D/logs/chain_full.done"
wait $pid_df; wait $pid_cf
echo "DONE $(cat $D/logs/chain_full.done)" > $D/logs/orch.done
hb "orchestrator done: $(cat $D/logs/orch.done)"
