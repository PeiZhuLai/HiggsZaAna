#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1
CON='regexp("MLPhoton.*Data_202[23]", JobBatchName)'
LOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/RegressMergedPhoton/condor/logs_monitor_resume_2223.txt
: > "$LOG"; log(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
zeros=0
while true; do
  mapfile -t st < <(condor_q pelai -constraint "$CON" -af JobStatus 2>/dev/null)
  tot=${#st[@]}; run=$(printf '%s\n' "${st[@]}"|grep -c '^2$'); idle=$(printf '%s\n' "${st[@]}"|grep -c '^1$'); held=$(printf '%s\n' "${st[@]}"|grep -c '^5$')
  log "queue=$tot (run=$run idle=$idle held=$held)"
  [ "$held" -gt 0 ] && log "  ⚠️ $held held"
  if [ "$tot" -eq 0 ]; then zeros=$((zeros+1)); else zeros=0; fi
  [ "$zeros" -ge 3 ] && { log "=== MLNanoAOD resume ALL DONE (3x confirm) ==="; break; }
  sleep 240
done
