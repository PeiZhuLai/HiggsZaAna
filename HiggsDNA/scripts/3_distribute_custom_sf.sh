#!/usr/bin/env bash
# Distribute the electron and muon SF / efficiency JSONs produced by
# 2_merge_custom_sf.py (2024 and 2025 only) into the HZgamma framework, renaming
# the hza_ prefix to hzg_. For triggers only the *_efficiencies.json variant is
# copied.
#
# Source : $HiggsDNADir/<era>_UL/hza_<...>.json        (output of 2_merge_custom_sf.py)
# Dests  : $eosDir/<era>/hzg_<...>.json                 (EOS, year-only subdir)
#          $afsDir/<era>_UL/hzg_<...>.json              (HZgamma repo, _UL subdir)
set -euo pipefail

HiggsDNADir="${HiggsDNADir:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data}"
eosDir="${eosDir:-/eos/project/h/htozg-dy-privatemc/pelai/HZg/JSON_custom}"
afsDir="${afsDir:-/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/systematics/JSONs}"

eras=(2024 2025)

# Electron + muon outputs of 2_merge_custom_sf.py, grouped by output suffix.
sf_bases=(elid muid)  # *_scalefactors.json
eff_bases=(      # *_efficiencies.json
    dielleg12trigger dielleg23trigger sielleg30trigger  # electron triggers
    eliso0p1 eliso0p15                                  # electron iso
    muiso0p1 muiso0p15                                  # muon iso
    mutrig8 mutrig17 mutrig24                           # muon triggers
)

era_dir() { echo "${1}_UL"; }

# distribute <base> <era> <suffix>
distribute() {
    local base="$1" era="$2" suffix="$3"
    local fname="${base}_${era}_${suffix}.json"
    local src="$HiggsDNADir/$(era_dir "$era")/hza_${fname}"
    if [[ ! -f "$src" ]]; then
        echo "WARNING: missing source, skipped: $src" >&2
        return 0
    fi
    local eos_dst="$eosDir/$era"
    local afs_dst="$afsDir/$(era_dir "$era")"
    mkdir -p "$eos_dst" "$afs_dst"
    cp -f -- "$src" "$eos_dst/hzg_${fname}"
    cp -f -- "$src" "$afs_dst/hzg_${fname}"
    echo "copied $src"
    echo "    -> $eos_dst/hzg_${fname}"
    echo "    -> $afs_dst/hzg_${fname}"
}

for era in "${eras[@]}"; do
    for base in "${sf_bases[@]}"; do
        distribute "$base" "$era" "scalefactors"
    done
    for base in "${eff_bases[@]}"; do
        distribute "$base" "$era" "efficiencies"
    done
done

echo "done."
