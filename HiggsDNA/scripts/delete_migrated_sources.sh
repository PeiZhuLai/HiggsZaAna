#!/usr/bin/env bash
# Delete HZa_merged sources after verified migration to NPS-25-014 (server-side eos rm -r).
# Run as the user (via `! bash ...`) since auto-mode blocks agent-initiated mass EOS deletes.
# rm goes to the EOS recycle bin (restorable, frees quota immediately).
set -uo pipefail
DLOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/logs/delete_sources_20260708.log
PJ=/eos/project/h/htozg-dy-privatemc/pelai/HZa
HPU=/eos/user/p/pelai/HZa
: > "$DLOG"
del(){ local ep=$1 path=$2 note=$3
  echo "[$(date +%H:%M:%S)] rm $path  ($note)" | tee -a "$DLOG"
  eos root://$ep.cern.ch rm -r "$path" >>"$DLOG" 2>&1 && echo "  OK" | tee -a "$DLOG" \
    || echo "  EXIT=$? (re-run to retry)" | tee -a "$DLOG"
}

# --- verified-backed: src==dst confirmed, full copy in NPS-25-014 ---
del eosproject-h "$PJ/parquet_merged_DNA_tmp" "backed: NPS 1854==1854"
del eosproject-h "$PJ/parquet_friend"         "backed: NPS 203719==203719"
del eosproject-h "$PJ/MLNanoAOD"              "backed: NPS 178371==178371"
del eoshome-p    "$HPU/parquet_merged_DNA_tmp" "backed: dup of project (in NPS)"
del eoshome-p    "$HPU/root_MVAcut"            "backed: NPS 532==532"
del eoshome-p    "$HPU/parquet_friend"         "backed: NPS homep 184121==184121"

# --- user-authorized, NOT migrated (old / test residue) ---
del eoshome-p    "$HPU/MLNanoAOD"                 "OLD 266GB: superseded by project version, delete-not-move"
del eosproject-h "$PJ/parquet_friend_datamltest"  "test residue ~2MB"
del eosproject-h "$PJ/parquet_friend_smoke_mlfix" "test residue ~0.7MB"

echo "[$(date +%H:%M:%S)] ALL DELETES ISSUED" | tee -a "$DLOG"
echo "=== quota after ==="
eos root://eosproject-h.cern.ch quota /eos/project/h/htozg-dy-privatemc 2>/dev/null | grep -E "project|filled" | tee -a "$DLOG"
eos root://eoshome-p.cern.ch quota /eos/user/p/pelai 2>/dev/null | grep -iE "pelai|filled" | tee -a "$DLOG"
