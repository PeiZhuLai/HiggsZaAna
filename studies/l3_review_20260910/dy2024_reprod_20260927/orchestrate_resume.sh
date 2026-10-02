#!/usr/bin/env bash
# RECOVERY COPY (2026-09-29): use ONLY if unit hza-dyveto-orch died. Identical to orchestrate.sh except the first
# driver round uses RESUME=1 (orchestrate.sh round 0 uses CLEAN_ANALYSIS_STATE=1 = would delete the pickle and
# resubmit all 7802 jobs). Round numbers start at START_ROUND (default 10) so logs do not overwrite.
#
# Persistent orchestrator for the FULL 2024 DY+jets re-production with the MC overlap veto
# (systemd --user unit hza-dyveto-orch; lxplus el9 kills setsid/nohup processes on logout).
# Started only after the small-batch test passed (logs/accept_test.txt: ACCEPT_OK).
#
#   1. full driver (run_dyveto_driver.sh MODE=full): ~7802 fpo=1 jobs, driver resubmits removed
#      jobs (<=5 attempts, then "retired") and merges at the end (--merge_outputs)
#   2. throughout: held-job handler + own-job priority + heartbeat every 5 min
#        memory hold  : RequestMemory 3000->8000->12000->16000->20000 (condor_qedit + condor_release,
#                       and the job's batch_submit file is raised too, only up, so a driver resubmission
#                       keeps it); at 20000 -> left held, reported
#        other holds  : (OnExitHold: exit!=0, e.g. xrdcp failure) condor_release, at most 10 times/job,
#                       then left held and reported (no condor_rm)
#   3. completeness (summary json per job). Jobs without summary (retired / failed) -> relaunch the
#      driver with RESUME=1 (--unretire_jobs resubmits them; every config uses the global redirector
#      cms-xrd-global.cern.ch) and it re-merges; at most 4 extra rounds, then STOP (needs decision)
#   4. chunk_rows == merged_rows, merged mtime >= max(chunk mtime)       (reconcile.py)
#   5. p2root --split --sideband-reweight-mode always -> root_P2Root/run3_bdt_inputs_fsrfix_dyveto/<s>/2024.root
#      (one local process at a time)
#   6. root inclusive == merged, train+val+test == inclusive, sideband branches present (stage3.py)
#
# Heartbeat: logs/heartbeat_orch.txt     Done marker: logs/orch.done   (DONE ... or STOP ...)
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/dy2024_reprod_20260927
IN=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyveto2024
RO=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix_dyveto
SAMPLES="DYJetsTo2E DYJetsTo2Mu DYJetsTo2Tau"
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
HB=$D/logs/heartbeat_orch.txt
HOLDLOG=$D/logs/held_handler.log
RELCNT=$D/logs/held_release_counts.txt
touch "$RELCNT"
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1; set -u
export X509_USER_PROXY=/tmp/x509up_u175325
hb() { echo "$(date '+%F %T') $*" >> "$HB"; }
finish() { echo "$1 $(date '+%F %T')" > $D/logs/orch.done; hb "orchestrator end: $1"; exit "${2:-0}"; }
# Identify my jobs by executable path only (schedd shared with other projects/sessions).
CONST='Owner=="pelai" && regexp("eos_logs/Bkg_MC_dyveto2024/",Cmd)'

next_rung() {
  local m=$1
  for r in 8000 12000 16000 20000; do [ "$m" -lt "$r" ] && { echo $r; return; }; done
  echo stuck; }

raise_subfile() {  # $1 = Cmd (executable path), $2 = new memory; only raise
  local sub cur
  sub=$(ls "$(dirname "$1")"/*_batch_submit_job*.txt 2>/dev/null | head -1)
  [ -n "$sub" ] || return 0
  cur=$(awk '/^RequestMemory/{print $3}' "$sub")
  [ -n "$cur" ] && [ "$cur" -lt "$2" ] && sed -i "s/^RequestMemory = .*/RequestMemory = $2/" "$sub"
  return 0; }

held_pass() {
  local out rc
  out=$(condor_q -constraint "$CONST && JobStatus==5" -af:t ClusterId ProcId HoldReasonCode RequestMemory Cmd 2>/dev/null); rc=$?
  [ $rc -ne 0 ] && { hb "condor_q rc=$rc in held_pass -- skipped (NOT treated as empty)"; return 0; }
  [ -z "$out" ] && return 0
  while IFS=$'\t' read -r c p code mem cmd; do
    id="$c.$p"
    reason=$(condor_q "$id" -af HoldReason 2>/dev/null | head -1)
    if [ "$code" = "34" ] || grep -qi "memory" <<< "$reason"; then
      nm=$(next_rung "${mem%.*}")
      if [ "$nm" = stuck ]; then
        grep -q "LEFT_HELD_MEM $id" "$HOLDLOG" 2>/dev/null || echo "$(date '+%F %T') LEFT_HELD_MEM $id at ${mem} MB | $reason" >> "$HOLDLOG"
        continue
      fi
      echo "$(date '+%F %T') $id memory hold ${mem}->${nm} MB: condor_qedit $id RequestMemory $nm; condor_release $id | $reason" >> "$HOLDLOG"
      condor_qedit "$id" RequestMemory "$nm" >> "$HOLDLOG" 2>&1
      raise_subfile "$cmd" "$nm"
      condor_release "$id" >> "$HOLDLOG" 2>&1
    else
      n=$(grep -c "^$id\$" "$RELCNT")
      if [ "$n" -lt 10 ]; then
        echo "$id" >> "$RELCNT"
        echo "$(date '+%F %T') $id hold code=$code (release $((n+1))/10): condor_release $id | $reason" >> "$HOLDLOG"
        condor_release "$id" >> "$HOLDLOG" 2>&1
      else
        grep -q "LEFT_HELD $id" "$HOLDLOG" || echo "$(date '+%F %T') LEFT_HELD $id after 10 releases, code=$code | $reason" >> "$HOLDLOG"
      fi
    fi
  done <<< "$out"
  grep -q "LEFT_HELD" "$HOLDLOG" 2>/dev/null && grep "LEFT_HELD" "$HOLDLOG" > $D/logs/NEEDS_DECISION_left_held.txt
  return 0; }

# Own-user ordering only (does not affect anyone else): these short jobs ahead of my own older idle ones.
prio_pass() {
  local ids
  ids=$(condor_q -constraint "$CONST && JobStatus==1 && JobPrio<5" -af ClusterId 2>/dev/null | tr '\n' ' ')
  [ -n "${ids// /}" ] && condor_prio +5 $ids >/dev/null 2>&1
  return 0; }

status_line() {
  local q ch
  q=$(condor_q -constraint "$CONST" -af JobStatus 2>/dev/null | sort | uniq -c | awk '{printf "%s:%s ",$2,$1}')
  ch=$(for s in $SAMPLES; do n=$(find $IN/${s}_2024 -maxdepth 2 -name '*_summary_job*.json' 2>/dev/null | wc -l); printf "%s=%s " $s $n; done)
  echo "condor(1=idle,2=run,5=held): ${q:-none} | summaries: $ch | proxy_left=$(voms-proxy-info -file /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/x509up_u175325 -timeleft 2>/dev/null)s"; }

# 2026-09-29 05:30-05:45 the schedd bigbird26 answered every condor_q with
# "SECMAN:2007 Read failure during security negotiation". A driver started then dies at
# condor_submit and would burn the resume rounds, so each round first waits for the schedd.
wait_schedd() {
  local n=0
  until condor_q -totals >/dev/null 2>&1; do
    n=$((n+1)); [ $((n % 6)) -eq 1 ] && hb "schedd not answering condor_q -- waiting before launching the driver"
    sleep 300
  done; }

run_driver_round() {   # $1 = round, $2 = RESUME
  local r=$1 resume=$2 pid
  wait_schedd
  rm -f $D/logs/driver_full.done
  MODE=full RESUME=$resume bash $D/run_dyveto_driver.sh > $D/logs/driver_full_round$r.log 2>&1 &
  pid=$!
  hb "[round $r] driver launched (pid $pid, RESUME=$resume) log=logs/driver_full_round$r.log"
  while kill -0 $pid 2>/dev/null; do
    held_pass; prio_pass
    hb "[round $r] $(status_line)"
    for _ in $(seq 30); do kill -0 $pid 2>/dev/null || break; sleep 10; done
  done
  wait $pid; local rc=$?
  hb "[round $r] driver exited rc=$rc; marker: $(cat $D/logs/driver_full.done 2>/dev/null || echo MISSING)"; }

hb "orchestrator start (pid $$, host $(hostname))"
round=${START_ROUND:-10}
run_driver_round $round 1   # RESUME: keep analysis_manager.pkl, never CLEAN_ANALYSIS_STATE
while :; do
  $PY $D/reconcile.py "$IN" $SAMPLES > $D/logs/reconcile_full.txt 2>&1; rcr=$?
  hb "[reconcile round $round] rc=$rcr"; cat $D/logs/reconcile_full.txt >> "$HB"
  [ $rcr -eq 0 ] && break
  if grep -q "notrun=" $D/logs/reconcile_full.txt; then
    round=$((round+1))
    [ $round -gt $(( ${START_ROUND:-10} + 4 )) ] && finish "STOP: jobs still without summary after 4 resume rounds -> needs decision (see logs/reconcile_full.txt, logs/held_handler.log)" 3
    run_driver_round $round 1
    continue
  fi
  finish "STOP: every job ran but chunk_rows != merged_rows or late chunks (merge race) -> see logs/reconcile_full.txt" 4
done

# p2root, one local process at a time
for s in $SAMPLES; do
  src=$IN/${s}_2024/merged_nominal.parquet
  bash $D/p2root_local.sh "$src" "$RO/$s/2024.root" --split --sideband-reweight-mode always > $D/logs/p2root_full_$s.log 2>&1
  hb "[p2root] $s rc=$? -> $RO/$s/2024.root"
done
$PY $D/stage3.py "$IN" "$RO" $SAMPLES > $D/logs/stage3_full.txt 2>&1; rc3=$?
hb "[stage3] rc=$rc3"; cat $D/logs/stage3_full.txt >> "$HB"
[ $rc3 -ne 0 ] && finish "STOP: stage-3 (root vs merged) failed -> logs/stage3_full.txt" 5
finish "DONE: production + merge + p2root complete, three-stage reconciliation OK (logs/reconcile_full.txt, logs/stage3_full.txt)" 0
