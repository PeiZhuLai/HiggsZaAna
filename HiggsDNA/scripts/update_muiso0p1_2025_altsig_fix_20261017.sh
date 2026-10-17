#!/usr/bin/env bash
# Update only the muon miniIso < 0.1 Run2025 efficiency JSON after the 2026-10-17 single-bin AltSig
# refit of |eta| < 0.9, pT 20-25 (spark_tnp: bin removed from MUISO_LEFT_PEAK_BINS, whose 35% nSigF
# ceiling had capped the failing signal). Same mechanism as update_mu{17,24}_sf_nodz_20261008.sh:
# back up raw + final, copy the spark_tnp output into custom_SF_raw, rebuild the final JSON with
# 2_merge_custom_sf.muon_iso_efficiency and validate it with correctionlib.
# 3_distribute_custom_sf.sh (copy into HZgamma) is NOT run here.
#
# Dry-run by default; add --go to write.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
H=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
MU=/eos/home-p/pelai/HZa/root_mTnP/efficiencies/muon/generalTracks/Z
BK=$H/backup_muiso0p1_2025_altsig_fix_20261017
era=2025

src=$MU/Run$era/NUM_MuIso0p1_DEN_HToZa_SignalMuons_Trigger/hzg_muiso0p1_${era}_efficiencies.json
raw=$H/${era}_UL/custom_SF_raw/hzg_muiso0p1_${era}_efficiencies.json
fin=$H/${era}_UL/hza_muiso0p1_${era}_efficiencies.json
[ -f "$src" ] || { echo "missing $src"; exit 2; }
[ -f "$raw" ] && [ -f "$fin" ] || { echo "missing target $raw or $fin"; exit 3; }
echo "$src ($(stat -c %y "$src" | cut -c1-16)) -> $raw -> $fin"
[ "$GO" -eq 1 ] || { echo "DRY-RUN (add --go to write)"; exit 0; }

# never overwrite an existing backup: a rerun would otherwise back up the already-updated files
if [ ! -f "$BK/$(basename "$fin")" ]; then
  mkdir -p "$BK"
  cp -p "$raw" "$fin" "$BK/"
fi
cp "$src" "$raw"

# the LCG setup script references unset variables: with -u it aborts this script silently
set +eu
source /cvmfs/sft.cern.ch/lcg/views/LCG_107_swan/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
set -eu
python3 - "$S/2_merge_custom_sf.py" "$(dirname "$H")/data" "$era" <<'EOF'
import importlib.util, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("merge", sys.argv[1])
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
out = m.muon_iso_efficiency(Path(sys.argv[2]), sys.argv[3], "muiso0p1")
m.validate_with_correctionlib([out])
print("wrote", out)
EOF
echo "backup: $BK"
