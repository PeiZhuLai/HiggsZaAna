#!/usr/bin/env bash
set -euo pipefail

HiggsDNADir="${HiggsDNADir:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data}"
egmSFDir="${egmSFDir:-/eos/home-p/pelai/HZa/root_TnP}"
muoSFDir="${muoSFDir:-/eos/home-p/pelai/HZa/root_mTnP/efficiencies/muon/generalTracks/Z}"
scriptsDir="${scriptsDir:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts}"

photon_eras=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)
electron_eras=(2024 2025)
muon_id_eras=(2024 2025)
muon_trigger_eras=(2024)
muon_iso_eras=(2024 2025)
all_eras=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024 2025)

era_dir() {
    case "$1" in
        2022preEE) echo "2022preEE_UL" ;;
        2022postEE) echo "2022postEE_UL" ;;
        2023preBPix) echo "2023preBPix_UL" ;;
        2023postBPix) echo "2023postBPix_UL" ;;
        2024) echo "2024_UL" ;;
        2025) echo "2025_UL" ;;
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

    rsync_json "$muoSFDir/Run${era}/$relpath" "$HiggsDNADir/$(era_dir "$era")"
}

for era in "${all_eras[@]}"; do
    mkdir -p "$(raw_dir "$era")"
done

###-------------------
### ----- Photon -----
###-------------------
for era in "${photon_eras[@]}"; do
    collect_egm_sf "$era" "hza_resolve_phid_${era}_sf"
    collect_egm_sf "$era" "hza_resolve_phid_lowpt_${era}_sf"
done

# 2023postBPix has a separate eta/phi hole region. The merge step will combine
# these with the nominal high-pT and low-pT photon SFs for the same era.
collect_egm_sf "2023postBPix" "hza_resolve_phid_2023postBPixHole_sf"
collect_egm_sf "2023postBPix" "hza_resolve_phid_lowpt_2023postBPixHole_sf"

###------------------------
### ----- Photon CSEV -----
###------------------------
for era in "${photon_eras[@]}"; do
    collect_egm_sf "$era" "hza_resolve_phcsev_lr9_summary_${era}_sf"
    collect_egm_sf "$era" "hza_resolve_phcsev_hr9_summary_${era}_sf"
done

# 2023postBPix also has separate CSEV maps for the eta/phi hole region.
collect_egm_sf "2023postBPix" "hza_resolve_phcsev_lr9_summary_2023postBPixHole_sf"
collect_egm_sf "2023postBPix" "hza_resolve_phcsev_hr9_summary_2023postBPixHole_sf"

###---------------------
### ----- Electron -----
###---------------------
# These four raw maps are merged into one hza_elid_*_scalefactors.json per era.
# The high/low pT maps currently use the nongap_highpT/nongap_lowpT TnP names.
elid_components=(gap nongap nongap_highpT nongap_lowpT)
for era in "${electron_eras[@]}"; do
    for component in "${elid_components[@]}"; do
        collect_egm_sf "$era" "hza_elid_${component}_${era}_sf"
    done
done

# The raw TnP directory names still end in _sf, but --exportJson now writes
# effdata/systdata/effmc/systmc correction names for trigger maps.
electron_trigger_effs=(dielleg12trigger dielleg23trigger sielleg30trigger)
for era in "${electron_eras[@]}"; do
    for trigger_eff in "${electron_trigger_effs[@]}"; do
        collect_egm_sf "$era" "hza_${trigger_eff}_gap_${era}_sf"
        collect_egm_sf "$era" "hza_${trigger_eff}_nongap_${era}_sf"
    done
done

electron_iso_effs=(elminiIso0p1 elminiIso0p15)
for era in "${electron_eras[@]}"; do
    for iso_eff in "${electron_iso_effs[@]}"; do
        collect_egm_sf "$era" "hza_${iso_eff}_gap_${era}_sf"
        collect_egm_sf "$era" "hza_${iso_eff}_nongap_${era}_sf"
    done
done

###---------------------
### ----- Muon ---------
###---------------------
for era in "${muon_id_eras[@]}"; do
    collect_muo_json "$era" "NUM_HToZa_SignalMuons_DEN_TrackerMuons/hza_muid_${era}_scalefactors.json"
done

for era in "${muon_trigger_eras[@]}"; do
    collect_muo_json "$era" "NUM_Mu8leg_DEN_HToZa_SignalMuons/hza_mutrig8_${era}_efficiencies.json"
    collect_muo_json "$era" "NUM_Mu17leg_DEN_HToZa_SignalMuons/hza_mutrig17_${era}_efficiencies.json"
    collect_muo_json "$era" "NUM_Mu24leg_DEN_HToZa_SignalMuons/hza_mutrig24_${era}_efficiencies.json"
done

for era in "${muon_iso_eras[@]}"; do
    collect_muo_json "$era" "NUM_MuIso0p1_DEN_HToZa_SignalMuons_Trigger/hzg_muiso0p1_${era}_efficiencies.json"
    collect_muo_json "$era" "NUM_MuIso0p15_DEN_HToZa_SignalMuons_Trigger/hzg_muiso0p15_${era}_efficiencies.json"
done

###------------------------
### ----- Convert ---------
###------------------------
python3 "$scriptsDir/2_merge_custom_sf.py"
