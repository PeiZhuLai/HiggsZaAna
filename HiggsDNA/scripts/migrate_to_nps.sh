#!/usr/bin/env bash
# One-off migration: HZa_merged outputs -> /eos/cms/store/group/phys_susy/NPS-25-014
# rsync (resumable) each source subdir to the NPS-25-014 workspace, then verify file counts.
# Sources are NOT deleted here -- deletion is a separate step after verification is confirmed.
set -uo pipefail

DST=/eos/cms/store/group/phys_susy/NPS-25-014/HZa_merged
PJ=/eos/project/h/htozg-dy-privatemc/pelai/HZa
HP=/eos/home-p/pelai/HZa
LOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/logs/migrate_nps_20260707.log

log(){ echo "[$(date +%m-%d\ %H:%M:%S)] $*" | tee -a "$LOG"; }

# label  src  dst
JOBS=(
  "MLNanoAOD           $PJ/MLNanoAOD            $DST/MLNanoAOD"
  "parquet_friend      $PJ/parquet_friend       $DST/parquet_friend"
  "parquet_friend_homep $HP/parquet_friend      $DST/parquet_friend_homep"
  "root_MVAcut         $HP/root_MVAcut          $DST/root_MVAcut"
)

log "=== migration start: ${#JOBS[@]} dirs -> $DST ==="
pids=()
for j in "${JOBS[@]}"; do
  set -- $j; label=$1; src=$2; dst=$3
  ( log "  [$label] rsync start ($src -> $dst)"
    rsync -a "$src/" "$dst/" \
      && log "  [$label] rsync DONE (exit 0)" \
      || log "  [$label] rsync EXIT=$? (re-run script to resume)"
  ) &
  pids+=($!)
done

# wait for all rsync jobs
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done

log "=== verification (src files == dst files) ==="
allok=1
for j in "${JOBS[@]}"; do
  set -- $j; label=$1; src=$2; dst=$3
  s=$(find "$src" -type f 2>/dev/null | wc -l)
  d=$(find "$dst" -type f 2>/dev/null | wc -l)
  if [ "$s" = "$d" ]; then log "  [$label] OK  src=$s dst=$d"
  else log "  [$label] MISMATCH  src=$s dst=$d  <-- re-run to resume"; allok=0; fi
done

if [ "$allok" = 1 ]; then log "=== ALL VERIFIED. Safe to delete sources (separate step). ==="
else log "=== INCOMPLETE. Re-run '$0' to resume rsync. ==="; fi
