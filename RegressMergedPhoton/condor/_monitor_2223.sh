#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1
CON='regexp("MLPhoton_Regressing_Data_202[23]", JobBatchName)'
OUTB=/eos/cms/store/group/phys_susy/pelai/HZa_merged/MLNanoAOD
LOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/RegressMergedPhoton/condor/logs_monitor_2223.txt
: > "$LOG"
log(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
while true; do
  mapfile -t st < <(condor_q pelai -constraint "$CON" -af JobStatus 2>/dev/null)
  tot=${#st[@]}
  idle=$(printf '%s\n' "${st[@]}" | grep -c '^1$')
  run=$(printf '%s\n' "${st[@]}"  | grep -c '^2$')
  held=$(printf '%s\n' "${st[@]}" | grep -c '^5$')
  nout=$(find "$OUTB" -path "*Data_202[23]*EE*.root" -o -path "*Data_202[23]*BPix*" -name "*.root" 2>/dev/null | wc -l)
  log "queue=$tot (idle=$idle run=$run held=$held)  outputs=$nout/50299"
  if [ "$held" -gt 0 ]; then log "  ⚠️ $held held jobs — check hold reason"; fi
  if [ "$tot" -eq 0 ]; then log "=== ALL DONE: queue empty, outputs=$nout ==="; break; fi
  sleep 300
done
