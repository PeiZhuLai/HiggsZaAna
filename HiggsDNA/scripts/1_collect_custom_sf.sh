

HiggsDNADir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data"
eras=(2022preEE_UL 2022postEE_UL 2023preBPix_UL 2023postBPix_UL 2024_UL)


for era in "${eras[@]}"; do
    mkdir -p "$HiggsDNADir/$era/custom_SF_raw"
done


###-------------------
### ----- Photon -----
###-------------------
egmSFDir="/eos/home-p/pelai/HZa/root_TnP"

rysnc -av $egmSFDir/hza_resolve_phid_2022preEE_sf/hza_resolve_phid_2022preEE_sf.json $HiggsDNADir/2022preEE_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_2022postEE_sf/hza_resolve_phid_2022postEE_sf.json $HiggsDNADir/2022postEE_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_2023preBPix_sf/hza_resolve_phid_2023preBPix_sf.json $HiggsDNADir/2023preBPix_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_2023postBPix_sf/hza_resolve_phid_2023postBPix_sf.json $HiggsDNADir/2023postBPix_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_2023postBPixHole_sf/hza_resolve_phid_2023postBPixHole_sf.json $HiggsDNADir/2023postBPix_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_2024_sf/hza_resolve_phid_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/

rysnc -av $egmSFDir/hza_resolve_phid_lowpt_2022preEE_sf/hza_resolve_phid_lowpt_2022preEE_sf.json $HiggsDNADir/2022preEE_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_lowpt_2022postEE_sf/hza_resolve_phid_lowpt_2022postEE_sf.json $HiggsDNADir/2022postEE_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_lowpt_2023preBPix_sf/hza_resolve_phid_lowpt_2023preBPix_sf.json $HiggsDNADir/2023preBPix_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_lowpt_2023postBPix_sf/hza_resolve_phid_lowpt_2023postBPix_sf.json $HiggsDNADir/2023postBPix_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_lowpt_2023postBPixHole_sf/hza_resolve_phid_lowpt_2023postBPixHole_sf.json $HiggsDNADir/2023postBPix_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_resolve_phid_lowpt_2024_sf/hza_resolve_phid_lowpt_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/
###---------------------
### ----- Electron ----- 
###---------------------
rysnc -av $egmSFDir/hza_elid_gap_2024_sf/hza_elid_gap_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_elid_nongap_2024_sf/hza_elid_nongap_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/

rysnc -av $egmSFDir/hza_elid_nongap_highpT_2024_sf/hza_elid_nongap_highpT_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_elid_nongap_lowpT_2024_sf/hza_elid_nongap_lowpT_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/

rysnc -av $egmSFDir/hza_dielleg12trigger_gap_2024_sf/hza_dielleg12trigger_gap_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_dielleg12trigger_nongap_2024_sf/hza_dielleg12trigger_nongap_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/

rysnc -av $egmSFDir/hza_dielleg23trigger_gap_2024_sf/hza_dielleg23trigger_gap_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_dielleg23trigger_nongap_2024_sf/hza_dielleg23trigger_nongap_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/

rysnc -av $egmSFDir/hza_sielleg30trigger_gap_2024_sf/hza_sielleg30trigger_gap_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/
rysnc -av $egmSFDir/hza_sielleg30trigger_nongap_2024_sf/hza_sielleg30trigger_nongap_2024_sf.json $HiggsDNADir/2024_UL/custom_SF_raw/
###---------------------
### ----- Muon ---------
###---------------------
muoSFDir="/eos/home-p/pelai/HZa/root_mTnP"

rysnc -av $muoSFDir/NUM_HToZa_SignalMuons_DEN_TrackerMuons/hza_muid_2024_scalefactors.json $HiggsDNADir/2024_UL/

rysnc -av $muoSFDir/NUM_Mu8leg_DEN_HToZa_SignalMuons/hza_mutrig8_2024_efficiencies.json $HiggsDNADir/2024_UL/

rysnc -av $muoSFDir/NUM_Mu17leg_DEN_HToZa_SignalMuons/hza_mutrig17_2024_efficiencies.json $HiggsDNADir/2024_UL/

rysnc -av $muoSFDir/NUM_Mu24leg_DEN_HToZa_SignalMuons/hza_mutrig24_2024_efficiencies.json $HiggsDNADir/2024_UL/

###------------------------
### ----- Convert ---------
###------------------------
scriptsDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts'

python3 scriptsDir/2_merge_custom_sf.py