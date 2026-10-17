#!/usr/bin/env bash
# Update only the muon trigger efficiency JSONs touched by the 2026-10-16 single-bin AltBkg refits
# (FFTConv blow-up fixed via MUON_TRIGGER_ALTBKG_BLOWUP in spark_tnp/run_single_fit.py):
#   Mu24leg Run2025 abseta_4_pt_3, Mu17leg Run2024 abseta_2_pt_2, Mu17leg Run2025 abseta_4_pt_3.
# Same mechanism as update_mu{17,24}_sf_nodz_20261008.sh: back up raw + final, copy the spark_tnp
# output into custom_SF_raw, rebuild the final JSON with 2_merge_custom_sf.muon_trigger_efficiency
# and validate it with correctionlib. 3_distribute_custom_sf.sh (copy into HZgamma) is NOT run here.
#
# Dry-run by default; add --go to write.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
H=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
MU=/eos/home-p/pelai/HZa/root_mTnP/efficiencies/muon/generalTracks/Z
BK=$H/backup_mutrig_altbkg_fix_20261016
PAIRS=(24:2025 17:2024 17:2025)

for p in "${PAIRS[@]}"; do
  leg=${p%%:*}; era=${p#*:}
  src=$MU/Run$era/NUM_Mu${leg}leg_DEN_HToZa_SignalMuons/hzg_mu${leg}leg_${era}_efficiencies.json
  raw=$H/${era}_UL/custom_SF_raw/hzg_mutrig${leg}_${era}_efficiencies.json
  fin=$H/${era}_UL/hza_mutrig${leg}_${era}_efficiencies.json
  [ -f "$src" ] || { echo "missing $src"; exit 2; }
  [ -f "$raw" ] && [ -f "$fin" ] || { echo "missing target $raw or $fin"; exit 3; }
  echo "Mu$leg $era: $src ($(stat -c %y "$src" | cut -c1-16)) -> $raw -> $fin"
  [ "$GO" -eq 1 ] || continue
  # never overwrite an existing backup: a rerun would otherwise back up the already-updated files
  if [ ! -f "$BK/$era/$(basename "$fin")" ]; then
    mkdir -p "$BK/$era"
    cp -p "$raw" "$fin" "$BK/$era/"
  fi
  cp "$src" "$raw"
done
[ "$GO" -eq 1 ] || { echo "DRY-RUN (add --go to write)"; exit 0; }

# the LCG setup script references unset variables: with -u it aborts this script silently
set +eu
source /cvmfs/sft.cern.ch/lcg/views/LCG_107_swan/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
set -eu
python3 - "$S/2_merge_custom_sf.py" "$(dirname "$H")/data" "${PAIRS[@]}" <<'EOF'
import importlib.util, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("merge", sys.argv[1])
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
base = Path(sys.argv[2])
outs = [m.muon_trigger_efficiency(base, p.split(":")[1], p.split(":")[0]) for p in sys.argv[3:]]
m.validate_with_correctionlib(outs)
for o in outs:
    print("wrote", o)
EOF
echo "backup: $BK"
