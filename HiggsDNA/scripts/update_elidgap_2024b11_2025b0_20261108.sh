#!/usr/bin/env bash
# Pull two re-fitted electron ID gap TnP outputs into the HZa SF JSONs and re-merge.
#   hza_elid_gap_2024_sf  altSigBkg bin11 (+eta gap, 100-500 GeV): failing background fixed flat (alphaF_2 = 0)
#   hza_elid_gap_2025_sf  altSig bin00 (-eta gap, 7-35 GeV): failing background fixed to the nominal shape
# Backups use cp -n so a re-run keeps the first (pre-change) copy.
# Default dry-run; --go to execute.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
T=/eos/home-p/pelai/HZa/root_TnP
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
BK=$D/backup_elidgap_2024b11_2025b0_20261108
run() { if [ "$GO" -eq 1 ]; then "$@"; else echo "DRY-RUN $*"; fi; }
for y in 2024 2025; do
  f=hza_elid_gap_${y}_sf
  [ -s $T/$f/$f.json ] || { echo "!! missing $T/$f/$f.json"; exit 3; }
  run mkdir -p $BK/${y}_UL/custom_SF_raw
  run cp -pn $D/${y}_UL/custom_SF_raw/$f.json $BK/${y}_UL/custom_SF_raw/
  run cp -pn $D/${y}_UL/*.json $BK/${y}_UL/
  run cp -p $T/$f/$f.json $D/${y}_UL/custom_SF_raw/
done
run python3 $S/2_merge_custom_sf.py --electron-eras 2024 2025
