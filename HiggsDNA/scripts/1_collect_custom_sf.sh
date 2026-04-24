#!/usr/bin/env bash
set -euo pipefail

HiggsDNADir="${HiggsDNADir:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data}"
egmSFDir="${egmSFDir:-/eos/home-p/pelai/HZa/root_TnP}"
muoSFDir="${muoSFDir:-/eos/home-p/pelai/HZa/root_mTnP/efficiencies/muon/generalTracks/Z/Run2024}"
scriptsDir="${scriptsDir:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts}"

eras=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)

era_dir() {
    case "$1" in
        2022preEE) echo "2022preEE_UL" ;;
        2022postEE) echo "2022postEE_UL" ;;
        2023preBPix) echo "2023preBPix_UL" ;;
        2023postBPix) echo "2023postBPix_UL" ;;
        2024) echo "2024_UL" ;;
        *)
            echo "Unknown era: $1" >&2
            return 1
            ;;
    esac
}

raw_dir() {
    echo "$HiggsDNADir/$(era_dir "$1")/custom_SF_raw"
}

rsync_json() {
    local src="$1"
    local dst="$2"

    mkdir -p "$dst"
    rsync -av -- "$src" "$dst/"
}

collect_egm_sf() {
    local era="$1"
    local sf_name="$2"

    rsync_json "$egmSFDir/$sf_name/$sf_name.json" "$(raw_dir "$era")"
}

collect_muo_json() {
    local era="$1"
    local relpath="$2"

    rsync_json "$muoSFDir/$relpath" "$HiggsDNADir/$(era_dir "$era")"
}

for era in "${eras[@]}"; do
    mkdir -p "$(raw_dir "$era")"
done

###-------------------
### ----- Photon -----
###-------------------
for era in "${eras[@]}"; do
    collect_egm_sf "$era" "hza_resolve_phid_${era}_sf"
    collect_egm_sf "$era" "hza_resolve_phid_lowpt_${era}_sf"
done

# 2023postBPix has a separate eta/phi hole region. The merge step will combine
# these with the nominal high-pT and low-pT photon SFs for the same era.
collect_egm_sf "2023postBPix" "hza_resolve_phid_2023postBPixHole_sf"
collect_egm_sf "2023postBPix" "hza_resolve_phid_lowpt_2023postBPixHole_sf"

###---------------------
### ----- Electron -----
###---------------------
# These four raw maps are merged into one hzg_elid_*_scalefactors.json for 2024.
# The high/low pT maps currently use the nongap_highpT/nongap_lowpT TnP names.
elid_components=(gap nongap nongap_highpT nongap_lowpT)
for component in "${elid_components[@]}"; do
    collect_egm_sf "2024" "hza_elid_${component}_2024_sf"
done

# Electron trigger custom SFs currently exist only for the 2024 campaign.
electron_trigger_sfs=(dielleg12trigger dielleg23trigger sielleg30trigger)
for trigger_sf in "${electron_trigger_sfs[@]}"; do
    collect_egm_sf "2024" "hza_${trigger_sf}_gap_2024_sf"
    collect_egm_sf "2024" "hza_${trigger_sf}_nongap_2024_sf"
done

###---------------------
### ----- Muon ---------
###---------------------
collect_muo_json "2024" "NUM_HToZa_SignalMuons_DEN_TrackerMuons/hza_muid_2024_scalefactors.json"
collect_muo_json "2024" "NUM_Mu8leg_DEN_HToZa_SignalMuons/hza_mutrig8_2024_efficiencies.json"
collect_muo_json "2024" "NUM_Mu17leg_DEN_HToZa_SignalMuons/hza_mutrig17_2024_efficiencies.json"
collect_muo_json "2024" "NUM_Mu24leg_DEN_HToZa_SignalMuons/hza_mutrig24_2024_efficiencies.json"

###------------------------
### ----- Convert ---------
###------------------------
python3 "$scriptsDir/2_merge_custom_sf.py"
