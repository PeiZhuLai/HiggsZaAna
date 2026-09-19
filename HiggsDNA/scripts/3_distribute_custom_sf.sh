#!/usr/bin/env bash
# Distribute the electron and muon SF / efficiency JSONs produced by
# 2_merge_custom_sf.py (electron and muon 2024-2026) into the HZgamma framework, renaming
# the hza_ prefix to hzg_. For triggers only the *_efficiencies.json variant is
# copied.
#
# Source : $HiggsDNADir/<era>_UL/hza_<...>.json        (output of 2_merge_custom_sf.py)
# Dests  : $eosDir/<era>/hzg_<...>.json                 (EOS, year-only subdir)
#          $afsDir/<era>_UL/hzg_<...>.json              (HZgamma repo, _UL subdir)
#
# Also publishes the in-house 2026 pileup weights to EOS (see the last section);
# those follow different rules and do not go through distribute().
set -euo pipefail

HiggsDNADir="${HiggsDNADir:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data}"
eosDir="${eosDir:-/eos/project/h/htozg-dy-privatemc/pelai/HZg/JSON_custom}"
afsDir="${afsDir:-/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/systematics/JSONs}"

# Electron and muon are measured for different sets of eras. 2026-08-29: muon TnP
# now covers 2026 as well (6/6 measurements, MC reused from 2025), so the two
# lists finally agree. Keeping them separate anyway -- photon still stops at 2024,
# and a future era will almost certainly land on one side before the other.
electron_eras=(2024 2025 2026)
muon_eras=(2024 2025 2026)

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

# ---------------------------------------------------------------------------
# In-house 2026 pileup weights
# ---------------------------------------------------------------------------
# These follow different rules from everything above, so they do not use
# distribute():
#   * the source is the HZgamma repo itself ($afsDir), not the hza_ output of
#     2_merge_custom_sf.py -- they were measured with
#     higgs_dna/scripts/pileup/ in the HZgamma repo, not by the TnP chain;
#   * there is no hza_ -> hzg_ rename. They are standard LUM-style correctionlib
#     payloads, so anyone can read them with correctionlib without knowing any
#     HZ naming convention; an hzg_ prefix would wrongly suggest a custom format;
#   * only the EOS copy is made -- the files already live in the HZgamma repo.
#
# Why they exist at all: there is no official 2026 LUM payload (neither in
# /cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/ nor as a Collisions26/PileUp
# directory), so HiggsDNA used to fall back to the 2025 one. Method and
# validation: doc/HZgamma/hzg_pileup_2026.md.
#
# 2026BD is the one to apply. The per-era B/D payloads are diagnostics and must
# NOT be applied while the MC is a single un-split Summer24 sample per year --
# each payload says so in its own "description" field too, so the warning
# travels with the file.
#
# The 2025 entry is the in-house RECOMPUTATION, published as a cross-check and
# as the reference for how the 2026 payload was produced. Production still
# applies the OFFICIAL Collisions25 payload: over the 267019 lumi sections both
# sides share, the two agree to 4e-5, and every visible difference traces to the
# LS sets (certification version) rather than to method or calibration. That is
# stated in the file's own description as well.
#
# Entries are "<era> <filename>": 2025 joining the list means this can no longer
# assume one era.
pileup_files=(
    "2026 puWeights_2026BD_Golden_Summer24_25ns_69200ub.json"
    "2026 puWeights_2026B_Golden_Summer24_25ns_69200ub.json"
    "2026 puWeights_2026D_Golden_Summer24_25ns_69200ub.json"
    "2025 puWeights_2025recalc_Golden_Summer24_25ns_69200ub.json"
)

distribute_pileup() {
    local era="$1" fname="$2"
    local src="$afsDir/$(era_dir "$era")/$fname"
    if [[ ! -f "$src" ]]; then
        echo "WARNING: missing source, skipped: $src" >&2
        return 0
    fi
    local eos_dst="$eosDir/$era"
    mkdir -p "$eos_dst"
    cp -f -- "$src" "$eos_dst/$fname"
    echo "copied $src"
    echo "    -> $eos_dst/$fname"
}

for entry in "${pileup_files[@]}"; do
    distribute_pileup $entry     # unquoted on purpose: splits "<era> <file>"
done

echo "done."
