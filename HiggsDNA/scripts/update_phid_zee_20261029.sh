#!/usr/bin/env bash
# Replace only the photon ID SFs with the measurement on the dedicated Z->ee photon TnP ntuples
# (root_TnP_phidZee: the analysis photon ID, fits checked bin by bin by the user on 2026-10-29).
# Copies the 12 egm JSONs (highpt/lowpt x 5 eras, + 2023postBPix hole highpt/lowpt) into custom_SF_raw and
# rebuilds hza_phid_<tag>[Hole]_scalefactors.json with 2_merge_custom_sf.merge_photon, validated with correctionlib.
# Every other file under systematics/data must stay byte-identical: checked against an md5 snapshot taken before.
#
# ERAS="2022preEE 2022postEE" limits the eras (default: all five); the 2023postBPix hole maps go with 2023postBPix.
# Only the files of the selected eras are rewritten, so jobs of other eras that read their JSONs are not disturbed.
# Dry-run by default; add --go to write.
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
H=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts
R=/eos/home-p/pelai/HZa/root_TnP_phidZee
ERAS=(${ERAS:-2022preEE 2022postEE 2023preBPix 2023postBPix 2024})
BK=$H/backup_phid_zee_20261029_$(echo "${ERAS[*]}" | tr ' ' '_')
declare -A ERADIR=([2022preEE]=2022preEE_UL [2022postEE]=2022postEE_UL [2023preBPix]=2023preBPix_UL [2023postBPix]=2023postBPix_UL [2024]=2024_UL)

pairs=()
for era in "${ERAS[@]}"; do
  for m in hza_resolve_phid_${era}_sf hza_resolve_phid_lowpt_${era}_sf; do
    pairs+=("$R/$m/$m.json|$H/${ERADIR[$era]}/custom_SF_raw/$m.json")
  done
done
if [[ " ${ERAS[*]} " == *" 2023postBPix "* ]]; then
  for m in hza_resolve_phid_2023postBPixHole_sf hza_resolve_phid_lowpt_2023postBPixHole_sf; do
    pairs+=("$R/$m/$m.json|$H/2023postBPix_UL/custom_SF_raw/$m.json")
  done
fi

for p in "${pairs[@]}"; do src=${p%%|*}; dst=${p#*|}
  [ -f "$src" ] || { echo "missing $src"; exit 2; }
  [ -f "$dst" ] || { echo "missing target $dst"; exit 3; }
  echo "$(stat -c %y "$src" | cut -c1-16)  $src -> $dst"
done
[ "$GO" -eq 1 ] || { echo "DRY-RUN (add --go to write): ${#pairs[@]} raw files"; exit 0; }

[ -e "$BK" ] && { echo "backup exists: $BK, stop"; exit 4; }
mkdir -p "$BK"
# md5 snapshot of everything except the files this script is allowed to change
(cd "$H" && find . -name "*.json" -not -path "./backup_*" -type f | sort | xargs md5sum) > "$BK/md5_before.txt"
for p in "${pairs[@]}"; do dst=${p#*|}; rel=${dst#$H/}; mkdir -p "$BK/$(dirname "$rel")"; cp -p "$dst" "$BK/$rel"; done
for era in "${ERAS[@]}"; do for f in $H/${ERADIR[$era]}/hza_phid_*_scalefactors.json; do rel=${f#$H/}; mkdir -p "$BK/$(dirname "$rel")"; cp -p "$f" "$BK/$rel"; done; done
for p in "${pairs[@]}"; do cp "${p%%|*}" "${p#*|}"; done

set +eu
source /cvmfs/sft.cern.ch/lcg/views/LCG_107_swan/x86_64-el9-gcc13-opt/setup.sh >/dev/null 2>&1
set -eu
python3 - "$S/2_merge_custom_sf.py" "$H" "${ERAS[@]}" <<'EOF'
import importlib.util, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("merge", sys.argv[1])
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
base = Path(sys.argv[2])
eras = sys.argv[3:]
outs = [m.merge_photon(base, era) for era in eras] + ([m.merge_photon(base, "2023postBPix", hole=True)] if "2023postBPix" in eras else [])
m.validate_with_correctionlib(outs)
import correctionlib
for o in outs:
    correctionlib.CorrectionSet.from_file(str(o))
    print("wrote", o)
EOF
(cd "$H" && find . -name "*.json" -not -path "./backup_*" -type f | sort | xargs md5sum) > "$BK/md5_after.txt"
echo "changed files:"
{ diff <(sort -k2 "$BK/md5_before.txt") <(sort -k2 "$BK/md5_after.txt") || true; } | grep '^>' | awk '{print "  "$3}'
bad=$( { diff <(sort -k2 "$BK/md5_before.txt") <(sort -k2 "$BK/md5_after.txt") || true; } | grep '^>' | awk '{print $3}' | grep -vE 'custom_SF_raw/hza_resolve_phid_|/hza_phid_[0-9A-Za-z]+_scalefactors\.json$' || true)
[ -z "$bad" ] || { echo "!! unexpected changes:"; echo "$bad"; exit 5; }
echo "only photon ID files changed; backup: $BK"
