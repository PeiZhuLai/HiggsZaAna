#!/usr/bin/env bash
# Pull the 11-08 batch of user-flagged electron TnP refits into the HZa SF JSONs and re-merge (2024-2026).
#   2024: sielleg30trigger_gap (bin00 passing shoulder, bin07 failing shoulder), dielleg23trigger_gap (altSig bin08),
#         dielleg12trigger_gap (altSigBkg bin07, mirror-merged template), elminiIso0p15_nongap (altSig bin17), elminiIso0p1_nongap (altSig bin12)
#   2025: elid_gap (altSig bin08/bin11, flat failing background)
#   2026: elid_nongap (altSigBkg bin08), sielleg30trigger_nongap (altSigBkg bin35, altSig bin08), dielleg23trigger_gap (altSig bin08/bin11)
# Each HZa raw file was checked byte-identical to the pre-refit EOS copy, so only these refits come in.
# Backups use cp -n so a re-run keeps the first (pre-change) copy.
# Default dry-run; --go to execute.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
T=/eos/home-p/pelai/HZa/root_TnP
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
BK=$D/backup_ele_batch3_20261108
FILES="hza_sielleg30trigger_gap_2024_sf hza_dielleg23trigger_gap_2024_sf hza_dielleg12trigger_gap_2024_sf
       hza_elminiIso0p15_nongap_2024_sf hza_elminiIso0p1_nongap_2024_sf hza_elid_gap_2025_sf
       hza_elid_nongap_2026_sf hza_sielleg30trigger_nongap_2026_sf hza_dielleg23trigger_gap_2026_sf"
run() { if [ "$GO" -eq 1 ]; then "$@"; else echo "DRY-RUN $*"; fi; }
for f in $FILES; do [ -s $T/$f/$f.json ] || { echo "!! missing $T/$f/$f.json"; exit 3; }; done
for y in 2024 2025 2026; do run mkdir -p $BK/${y}_UL/custom_SF_raw; run cp -pn $D/${y}_UL/*.json $BK/${y}_UL/; done
for f in $FILES; do
  y=$(echo $f | grep -o '20[0-9][0-9]')
  run cp -pn $D/${y}_UL/custom_SF_raw/$f.json $BK/${y}_UL/custom_SF_raw/
  run cp -p $T/$f/$f.json $D/${y}_UL/custom_SF_raw/
done
run python3 $S/2_merge_custom_sf.py --electron-eras 2024 2025 2026
