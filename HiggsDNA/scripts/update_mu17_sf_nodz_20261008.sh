#!/usr/bin/env bash
# Update only the double-muon upper-leg (Mu17leg) trigger efficiency JSONs of 2024-2026 after the
# muon TnP re-measurement without pair_dz < 0.2 (the DZ filter is kept in Mu8leg only; spark_tnp
# run_mu17_nodz.sh). Every other SF is left untouched, which a full
# 1_collect + 2_merge would not guarantee (it rebuilds all photon/electron/muon outputs).
#
# Per era: back up custom_SF_raw/hzg_mutrig17_<era> and <era>_UL/hza_mutrig17_<era>, copy the new
# spark_tnp output into custom_SF_raw (as 1_collect_custom_sf.sh does), rebuild the final JSON with
# 2_merge_custom_sf.muon_trigger_efficiency and validate it with correctionlib.
# 3_distribute_custom_sf.sh (copy into HZgamma) is NOT run here.
#
# Dry-run by default; add --go to write.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
H=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
MU=/eos/home-p/pelai/HZa/root_mTnP/efficiencies/muon/generalTracks/Z
BK=$H/backup_mu17_nodz_20261008
ERAS=(2024 2025 2026)

for era in "${ERAS[@]}"; do
  src=$MU/Run$era/NUM_Mu17leg_DEN_HToZa_SignalMuons/hzg_mu17leg_${era}_efficiencies.json
  raw=$H/${era}_UL/custom_SF_raw/hzg_mutrig17_${era}_efficiencies.json
  fin=$H/${era}_UL/hza_mutrig17_${era}_efficiencies.json
  [ -f "$src" ] || { echo "missing $src"; exit 2; }
  echo "$era: $src ($(stat -c %y "$src" | cut -c1-16)) -> $raw -> $fin"
  [ "$GO" -eq 1 ] || continue
  # never overwrite an existing backup: a rerun would otherwise back up the already-updated raw file
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
python3 - "$S/2_merge_custom_sf.py" "$(dirname "$H")/data" "${ERAS[@]}" <<'EOF'
import importlib.util, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("merge", sys.argv[1])
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
base = Path(sys.argv[2])
outs = [m.muon_trigger_efficiency(base, era, "17") for era in sys.argv[3:]]
m.validate_with_correctionlib(outs)
for o in outs:
    print("wrote", o)
EOF
echo "backup: $BK"
