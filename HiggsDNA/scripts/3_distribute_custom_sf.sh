#!/usr/bin/env bash
# Distribute the electron and muon SF / efficiency JSONs produced by
# 2_merge_custom_sf.py (electron 2024-2026, muon 2024-2025) into the HZgamma framework, renaming
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

# Electron and muon are measured for different sets of eras: 2026 has electron SFs
# (2026 data vs 2024 MC) but no muon TnP yet. Keeping the lists separate avoids
# emitting a WARNING for every muon file that was never supposed to exist.
electron_eras=(2024 2025 2026)
muon_eras=(2024 2025)

# Outputs of 2_merge_custom_sf.py, grouped by output suffix.
electron_sf_bases=(elid)                                       # *_scalefactors.json
electron_eff_bases=(                                           # *_efficiencies.json
    dielleg12trigger dielleg23trigger sielleg30trigger          # electron triggers
    eliso0p1 eliso0p15                                         # electron iso
)
muon_sf_bases=(muid)
muon_eff_bases=(
    muiso0p1 muiso0p15                                         # muon iso
    mutrig8 mutrig17 mutrig24                                  # muon triggers
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

for era in "${electron_eras[@]}"; do
    for base in "${electron_sf_bases[@]}"; do
        distribute "$base" "$era" "scalefactors"
    done
    for base in "${electron_eff_bases[@]}"; do
        distribute "$base" "$era" "efficiencies"
    done
done

for era in "${muon_eras[@]}"; do
    for base in "${muon_sf_bases[@]}"; do
        distribute "$base" "$era" "scalefactors"
    done
    for base in "${muon_eff_bases[@]}"; do
        distribute "$base" "$era" "efficiencies"
    done
done

echo "done."
