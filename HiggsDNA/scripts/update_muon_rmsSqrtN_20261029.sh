#!/usr/bin/env bash
# Update the HZa muon SF JSONs (2024-2026) after the muon fit-model uncertainty was changed to the EGM RMS method,
# RMS/sqrt(N) with RMS = sqrt(sum((c-m)^2)/(N-1)), N = 5 (spark_tnp/prepare.py, 2026-10-29; was sqrt(sum((c-m)^2))).
# Central values do not change, only the systematic uncertainties.
# Same mechanism as update_mutrig_altbkg_fix_20261016.sh: back up raw + final, copy the spark_tnp output into
# custom_SF_raw (iso, trigger) or straight into the era dir (muid), rebuild the final iso/trigger JSONs with
# 2_merge_custom_sf.py and validate them with correctionlib. 3_distribute_custom_sf.sh (copy into HZgamma) is NOT run here.
#
# Dry-run by default; add --go to write.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
H=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
MU=/eos/home-p/pelai/HZa/root_mTnP/efficiencies/muon/generalTracks/Z
BK=$H/backup_muon_rmsSqrtN_20261029
ERAS=(2024 2025 2026)

copies=()   # src|dst
finals=()   # files to back up
for era in "${ERAS[@]}"; do
  d=$H/${era}_UL
  copies+=("$MU/Run$era/NUM_HToZa_SignalMuons_DEN_TrackerMuons/hza_muid_${era}_scalefactors.json|$d/hza_muid_${era}_scalefactors.json")
  finals+=("$d/hza_muid_${era}_scalefactors.json")
  for leg in 8 17 24; do
    copies+=("$MU/Run$era/NUM_Mu${leg}leg_DEN_HToZa_SignalMuons/hzg_mu${leg}leg_${era}_efficiencies.json|$d/custom_SF_raw/hzg_mutrig${leg}_${era}_efficiencies.json")
    finals+=("$d/custom_SF_raw/hzg_mutrig${leg}_${era}_efficiencies.json" "$d/hza_mutrig${leg}_${era}_efficiencies.json")
  done
  for iso in 0p1 0p15; do
    copies+=("$MU/Run$era/NUM_MuIso${iso}_DEN_HToZa_SignalMuons_Trigger/hzg_muiso${iso}_${era}_efficiencies.json|$d/custom_SF_raw/hzg_muiso${iso}_${era}_efficiencies.json")
    finals+=("$d/custom_SF_raw/hzg_muiso${iso}_${era}_efficiencies.json" "$d/hza_muiso${iso}_${era}_efficiencies.json")
  done
done

for c in "${copies[@]}"; do
  src=${c%%|*}; dst=${c#*|}
  [ -f "$src" ] || { echo "missing $src"; exit 2; }
  [ -f "$dst" ] || { echo "missing target $dst"; exit 3; }
  echo "$(stat -c %y "$src" | cut -c1-16)  $src -> $dst"
done
for f in "${finals[@]}"; do [ -f "$f" ] || { echo "missing target $f"; exit 3; }; done
[ "$GO" -eq 1 ] || { echo "DRY-RUN (add --go to write): ${#copies[@]} copies, ${#finals[@]} files backed up"; exit 0; }

# never overwrite an existing backup: a rerun would otherwise back up the already-updated files
if [ ! -d "$BK" ]; then
  for f in "${finals[@]}"; do
    rel=${f#$H/}; mkdir -p "$BK/$(dirname "$rel")"; cp -p "$f" "$BK/$rel"
  done
fi
for c in "${copies[@]}"; do cp "${c%%|*}" "${c#*|}"; done

# the LCG setup script references unset variables: with -u it aborts this script silently
set +eu
source /cvmfs/sft.cern.ch/lcg/views/LCG_107_swan/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
set -eu
python3 - "$S/2_merge_custom_sf.py" "$H" "${ERAS[@]}" <<'EOF'
import importlib.util, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("merge", sys.argv[1])
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
base = Path(sys.argv[2])
outs = []
for era in sys.argv[3:]:
    for leg in ("8", "17", "24"):
        outs.append(m.muon_trigger_efficiency(base, era, leg))
    for iso in ("muiso0p1", "muiso0p15"):
        outs.append(m.muon_iso_efficiency(base, era, iso))
m.validate_with_correctionlib(outs + [base / f"{era}_UL" / f"hza_muid_{era}_scalefactors.json" for era in sys.argv[3:]])
for o in outs:
    print("wrote", o)
EOF
echo "backup: $BK"
