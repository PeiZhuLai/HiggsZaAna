#!/usr/bin/env bash
# Pull the dielleg23trigger_gap 2026 altSig bin08/bin11 passing refit (shoulder Gaussian, monotonic background)
# into the HZa SF JSON and re-merge (2024-2026 electron).
# Precondition: the HZa raw file is byte-identical to the pre-refit EOS copy (checked below), so only this refit comes in.
# Backups use cp -n so a re-run keeps the first (pre-change) copy.
# Default dry-run; --go to execute.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
T=/eos/home-p/pelai/HZa/root_TnP
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
PRE=/afs/cern.ch/work/p/pelai/HZa/TnP/CMSSW_11_2_0/src/egm_tnp_analysis/diag_20261108/backup/hza_dielleg23trigger_gap_2026_sf_summary4
BK=$D/backup_dielleg23gap2026_passfix_20261108
f=hza_dielleg23trigger_gap_2026_sf; y=2026
run() { if [ "$GO" -eq 1 ]; then "$@"; else echo "DRY-RUN $*"; fi; }
[ -s $T/$f/$f.json ] || { echo "!! missing $T/$f/$f.json"; exit 3; }
cmp -s $PRE/$f.json $D/${y}_UL/custom_SF_raw/$f.json || { echo "!! HZa raw differs from pre-refit EOS copy"; exit 4; }
for yy in 2024 2025 2026; do run mkdir -p $BK/${yy}_UL/custom_SF_raw; run cp -pn $D/${yy}_UL/*.json $BK/${yy}_UL/; done
run cp -pn $D/${y}_UL/custom_SF_raw/$f.json $BK/${y}_UL/custom_SF_raw/
run cp -p $T/$f/$f.json $D/${y}_UL/custom_SF_raw/
run python3 $S/2_merge_custom_sf.py --electron-eras 2024 2025 2026
