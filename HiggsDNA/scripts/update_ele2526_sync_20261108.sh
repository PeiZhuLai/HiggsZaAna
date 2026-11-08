#!/usr/bin/env bash
# Sync all 2025 and 2026 electron TnP raw JSONs (16 measurements per year) from EOS into HZa and re-merge.
#   - The HZa copies were from 10-24. The EOS copies were regenerated on 10-27 with the |altMC - nominal MC|
#     fix (uncertainties only grow, central values bit-identical; elid also gains unc_pass_stat/unc_pass_syst),
#     the same fix that 2024 received on 10-26.
#   - hza_elid_gap_2026_sf also carries the 11-08 altSig refit of bin00/bin03 (failing background fixed to the
#     nominal shape), which changes sf_pass/sf_fail of those two cells.
# The muon raw files (hzg_mu*) in custom_SF_raw are not touched. collect_sfunc_20261024.sh is not used because it
# would also re-copy the old photon ID raw files (photon ID now comes from the Z->ee measurement).
# Backups use cp -n so a re-run keeps the first (pre-change) copy.
# Default dry-run; --go to execute.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
T=/eos/home-p/pelai/HZa/root_TnP
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
BK=$D/backup_ele2526_sync_20261108
MEAS="elid_gap elid_nongap elid_nongap_highpT elid_nongap_lowpT
      elminiIso0p1_gap elminiIso0p1_nongap elminiIso0p1_nongap_highpT
      elminiIso0p15_gap elminiIso0p15_nongap elminiIso0p15_nongap_highpT
      dielleg12trigger_gap dielleg12trigger_nongap dielleg23trigger_gap dielleg23trigger_nongap
      sielleg30trigger_gap sielleg30trigger_nongap"
run() { if [ "$GO" -eq 1 ]; then "$@"; else echo "DRY-RUN $*"; fi; }
n=0
for y in 2025 2026; do
  for m in $MEAS; do
    f=hza_${m}_${y}_sf
    [ -s $T/$f/$f.json ] || { echo "!! missing $T/$f/$f.json"; exit 3; }
    [ -s $D/${y}_UL/custom_SF_raw/$f.json ] || { echo "!! HZa has no $f.json (unexpected name)"; exit 4; }
    n=$((n+1))
  done
done
[ $n -eq 32 ] || { echo "!! expected 32 files, got $n"; exit 5; }
for y in 2025 2026; do
  run mkdir -p $BK/${y}_UL/custom_SF_raw
  run cp -pn $D/${y}_UL/*.json $BK/${y}_UL/
  for m in $MEAS; do
    f=hza_${m}_${y}_sf
    run cp -pn $D/${y}_UL/custom_SF_raw/$f.json $BK/${y}_UL/custom_SF_raw/
    run cp -p $T/$f/$f.json $D/${y}_UL/custom_SF_raw/
  done
done
run python3 $S/2_merge_custom_sf.py --electron-eras 2025 2026
