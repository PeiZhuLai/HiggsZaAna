#!/usr/bin/env bash
# Pull two re-fitted 2024 electron TnP outputs into the HZa SF JSONs and re-merge.
#   hza_elid_gap_2024_sf               altSig bin03 (+eta gap, 7-35 GeV): new MC reference, data refit
#   hza_dielleg23trigger_nongap_2024_sf altBkg bin07 (eta 2.0-2.5, 7-20 GeV): passing template from the eta mirror bins
# Only these two raw files are copied. collect_sfunc_20261024.sh is not used because it would
# also re-copy the old photon ID raw files (photon ID now comes from the Z->ee measurement).
# Backups use cp -n so a re-run keeps the first (pre-change) copy.
# Default dry-run; --go to execute.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
T=/eos/home-p/pelai/HZa/root_TnP
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
BK=$D/backup_elidgap_d23_20261107
run() { if [ "$GO" -eq 1 ]; then "$@"; else echo "DRY-RUN $*"; fi; }
run mkdir -p $BK/2024_UL/custom_SF_raw
for f in hza_elid_gap_2024_sf hza_dielleg23trigger_nongap_2024_sf; do
  [ -s $T/$f/$f.json ] || { echo "!! missing $T/$f/$f.json"; exit 3; }
  run cp -pn $D/2024_UL/custom_SF_raw/$f.json $BK/2024_UL/custom_SF_raw/
  run cp -p $T/$f/$f.json $D/2024_UL/custom_SF_raw/
done
run cp -pn $D/2024_UL/*.json $BK/2024_UL/
run python3 $S/2_merge_custom_sf.py --electron-eras 2024
