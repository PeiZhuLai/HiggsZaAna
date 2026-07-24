#!/usr/bin/env bash
# Chunked deletion of the 3 big migrated sources that failed `eos rm -r` with
# EXIT=7 "directory tree exceeds the configured query limit". We delete each
# immediate child (per-tag) separately so every rm stays under the query limit,
# then remove the now-empty parent. All 3 are verified-backed in NPS-25-014.
# Run as the user (via `! bash ...`); auto-mode blocks agent-initiated mass deletes.
set -uo pipefail
DLOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/logs/delete_big_chunked_20260709.log
: > "$DLOG"
log(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$DLOG"; }

# ep  parent_path
BIG=(
  "eosproject-h /eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD"
  "eosproject-h /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_friend"
  "eoshome-p    /eos/user/p/pelai/HZa/parquet_friend"
)

for entry in "${BIG[@]}"; do
  set -- $entry; ep=$1; parent=$2
  log "=== chunked rm $parent (via $ep) ==="
  if ! eos root://$ep.cern.ch stat "$parent" >/dev/null 2>&1; then log "  parent gone, skip"; continue; fi
  # immediate children (per-tag dirs); delete each separately
  mapfile -t kids < <(eos root://$ep.cern.ch ls "$parent" 2>/dev/null)
  log "  ${#kids[@]} children to delete"
  for k in "${kids[@]}"; do
    [ -z "$k" ] && continue
    if eos root://$ep.cern.ch rm -r "$parent/$k" >>"$DLOG" 2>&1; then
      log "  OK  $k"
    else
      # child itself over the limit -> go one level deeper
      log "  child $k over limit? recursing one level"
      for gk in $(eos root://$ep.cern.ch ls "$parent/$k" 2>/dev/null); do
        eos root://$ep.cern.ch rm -r "$parent/$k/$gk" >>"$DLOG" 2>&1 \
          && log "    OK  $k/$gk" || log "    FAIL $k/$gk (re-run)"
      done
      eos root://$ep.cern.ch rm -r "$parent/$k" >>"$DLOG" 2>&1 && log "  OK  $k (after deep)" || log "  FAIL $k"
    fi
  done
  # remove now-empty parent
  eos root://$ep.cern.ch rm -r "$parent" >>"$DLOG" 2>&1 && log "  parent removed: $parent" || log "  parent NOT empty/failed: $parent (re-run)"
done

log "=== quota after ==="
eos root://eosproject-h.cern.ch quota /eos/project/h/htozg-dy-privatemc 2>/dev/null | grep -E "project .*%" | tee -a "$DLOG"
eos root://eoshome-p.cern.ch quota /eos/user/p/pelai 2>/dev/null | grep -E "%" | tee -a "$DLOG"
log "DONE (re-run this script if any FAIL / parent-not-empty remain)"
