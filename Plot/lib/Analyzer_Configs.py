import os
import sys
from ROOT import *

class Analyzer_Config:
    def __init__(self, channel, year, region, mva):
        self.channel            = channel
        self.year               = year
        self.version            = 'UL'
        self.region             = region
        self.mva                = mva
        self.sample_loc         = 'NONE'
        self.out_dir            = 'NONE'
        self.out_region_name    = 'NONE'
        self.root_output_name   = 'NONE'
        self.plot_output_path   = 'NONE'
        self.BDT_filename       = 'NONE'
        self.mvaCut             = {}
        self.sig_names          = []
        self.bkg_names          = []
        self.samp_names         = []
        self.sys_names          = []

        self.Config_Analyzer()

    def Config_Analyzer(self):
        if self.channel == 'inclusive' or self.channel == 'ggH' or self.channel == 'VBF' or self.channel == 'WH_3l' or self.channel == 'ZH_4l' or self.channel == 'ttH_lep':

            if self.region == 1:
                self.out_region_name = self.version + '_SR'
            elif self.region == 2:
                self.out_region_name = self.version + '_CR'
            elif self.mva:
                self.out_region_name = self.version + '_mva'
            else:
                self.out_region_name = self.version

            if self.year == '2016':
                self.sample_loc       = '/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_output/UL/16'
                self.out_dir          = 'plots_16UL'
                self.root_output_name = "ALP_plot_data16_{0}.root".format(self.out_region_name)
                # self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2016.pkl"
                self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/train_MVA/model_Za_BDT_passedEvents.pkl"
                #mvaCut = 0.8675
            elif self.year == '-2016':
                self.sample_loc       = '/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_output/UL/16APV'
                self.out_dir          = 'plots_16APVUL'
                self.root_output_name = "ALP_plot_data16APV_{0}}.root".format(self.out_region_name)
                # self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2016.pkl"
                self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/train_MVA/model_Za_BDT_passedEvents.pkl"
                #mvaCut = 0.8365
            elif self.year == '2017':
                self.sample_loc       = '/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_output/UL/17'
                self.out_dir          = 'plots_17UL'
                self.root_output_name = "ALP_plot_data17_{0}.root".format(self.out_region_name)
                # self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2017.pkl"
                self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/train_MVA/model_Za_BDT_passedEvents.pkl"                
                #mvaCut = 0.8365
            elif self.year == '2018':
                self.sample_loc       = '/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_output/UL/18'
                self.out_dir          = 'plots_18UL'
                self.root_output_name = "ALP_plot_data18_{0}.root".format(self.out_region_name)
                # self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/MVA/weight/nodR/model_ALP_BDT_param_2018.pkl"
                self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/train_MVA/model_Za_BDT_passedEvents.pkl"
                #mvaCut = 0.9766
            elif self.year == 'run2Rereco':
                self.sample_loc       = '/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_output/Rereco/run2'
            elif self.year == 'run2':
                self.sample_loc       = '/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_output/UL/run2'
                self.out_dir          = 'plots_run2UL'
                self.root_output_name = "ALP_plot_run2_{0}.root".format(self.out_region_name)
                # self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/MVA/weight/UL/model_ALP_BDT_param.pkl"
                self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/ALP/Analysis_code/train_MVA/model_Za_BDT_passedEvents.pkl"
                self.mvaCut           = {'M1':0.955, 'M2':0.98, 'M3':0.985, 'M4':0.98, 'M5':0.985, 'M6':0.99, 'M7':0.985, 'M8':0.99, 'M9':0.99, 'M10':0.99, 'M15':0.99, 'M20':0.99, 'M25':0.985, 'M30':0.98}
            elif self.year == 'run3':
                self.sample_loc       = '/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal'
                self.out_dir          = '/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/variables_dataVmc'
                self.root_output_name = "ALP_plot_run3_{0}.root".format(self.out_region_name)
                self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_run3.pkl"
                self.mvaCut           = {'M1':0.955, 'M2':0.98, 'M3':0.985, 'M4':0.98, 'M5':0.985, 'M6':0.99, 'M7':0.985, 'M8':0.99, 'M9':0.99, 'M10':0.99, 'M15':0.99, 'M20':0.99, 'M25':0.985, 'M30':0.98}
            elif self.year == 'run3_NFlow':
                self.sample_loc       = '/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_NFlow'
                self.out_dir          = '/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/variables_dataVmc_NFlow'
                self.root_output_name = "ALP_plot_run3_NFlow_{0}.root".format(self.out_region_name)
                self.BDT_filename     = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_run3_NFlow.pkl"
                self.mvaCut           = {'M1':0.955, 'M2':0.98, 'M3':0.985, 'M4':0.98, 'M5':0.985, 'M6':0.99, 'M7':0.985, 'M8':0.99, 'M9':0.99, 'M10':0.99, 'M15':0.99, 'M20':0.99, 'M25':0.985, 'M30':0.98}
            else:
                print('do not included at 2016/2017/2018!')
                exit(0)

            if self.year == 'run2Rereco':
                self.sig_names  = ['M1', 'M5', 'M15', 'M30']
            else:    
                self.sig_names  = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M15', 'M20', 'M25', 'M30']
                # self.sig_names  = ['M5', 'M15', 'M30']

            self.years_sig  = ["2022preEE","2022postEE","2023preBPix","2023postBPix", "2024"]
            self.years_22   = ["2022preEE", "2022postEE"]  # 背景（2022）
            self.years_23   = ["2023preBPix","2023postBPix", "2024"]  # 背景（2023）
            self.years_dyll = ["2022preEE","2022postEE","2023preBPix","2023postBPix", "2024"]
            self.bkg_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
            self.bkg_2023 = ["DYGto2LG_10to100"]
            self.bkg_dyll = ["DYJetsToLL"]
            # 2024 has no inclusive DYJetsToLL production; the DY+jets input is the
            # flavor-split DYJetsTo2E/2Mu/2Tau. Loaded under the DYJetsToLL sample key
            # for 2024 only (see Plot_Helper._run3_sources_for_sample). Keeps run3 DY+jets
            # consistent with select_lib.BKG_SAMPLES_BY_YEAR / the pseudo-data closure.
            self.bkg_dyll_2024 = ["DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]

            self.bkg_names  = ['DYJetsToLL', 'DYGto2LG']
            self.samp_names = self.bkg_names + self.sig_names + ['Data']
            self.plot_output_path = "{0}/plot_{1}".format(self.out_dir, self.out_region_name)
            # self.sys_names  = ['CMS_eff_g_up','CMS_eff_g_dn','CMS_pileup_up','CMS_pileup_dn','CMS_eff_lep_up','CMS_eff_lep_dn']
            self.sys_names  = ['weight_hlt_sf_up','weight_hlt_sf_down','weight_pu_reweight_sf_up','weight_pu_reweight_sf_down','weight_electron_wplid_sf_SelectedElectron_up','weight_electron_wplid_sf_SelectedElectron_down', 'weight_electron_reco_sf_SelectedElectron_up', 'weight_electron_reco_sf_SelectedElectron_down', 'weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_up', 'weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_down', 'weight_muon_looseid_sf_SelectedMuon_up', 'weight_muon_looseid_sf_SelectedMuon_down', 'weight_muon_iso_sf_SelectedMuon_up', 'weight_muon_iso_sf_SelectedMuon_down', 'weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_up', 'weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_down','weight_photon_id_sf_SelectedPhoton_up', 'weight_photon_id_sf_SelectedPhoton_down', 'weight_photon_csev_sf_SelectedPhoton_up', 'weight_photon_csev_sf_SelectedPhoton_down']
            # self.sys_names  = ['weight_pu_reweight_sf_up','weight_pu_reweight_sf_down'] # Too big
            # self.sys_names  = ['weight_electron_wplid_sf_SelectedElectron_up','weight_electron_wplid_sf_SelectedElectron_down']
            # self.sys_names  = ['weight_muon_looseid_sf_SelectedMuon_up','weight_muon_looseid_sf_SelectedMuon_down']
        else:
            print("channel is invalid: channel = %s" %self.channel)
            sys.exit()

    def Print_Config(self):
        print('Running analysis in channel: %s' %self.channel)
        print('getting ntuples from: %s' %self.sample_loc)
        print('using signals: ')
        print(self.sig_names)
        print('using backgrounds:')
        print(self.bkg_names)