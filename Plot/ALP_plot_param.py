####################################################
####################################################

import os
import sys
import numpy as np
import gc

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend, Total_Unc
from Analyzer_Helper import getMassSigma
import Analyzer_Configs as AC
import Plot_Configs     as PC

from Analyzer_ALP import PIso2D

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle
import copy 
import random

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-y", "--Year", dest="year", default="run2", help="which year's datasetes")
parser.add_argument("--region", dest="region", type=int, default=0, help="0 for full region, 1 for signal region, 2 for sideband region")
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument('--cut', dest='cut', action='store_true', default=False, help='apply mva')
parser.add_argument("--cutVal", dest="cutVal", type=float, default=0.0, help="mva cut value")
parser.add_argument("--mA", dest="mA", default="M5", help="ALP mass")
parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')
parser.add_argument('-ln', '--ln', dest='ln', action='store_true', default=False, help='log plot?')
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
args = parser.parse_args()



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


# load the model from disk
# model = pickle.load(open(BDT_filename, 'rb'))

def main():

    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year, args.region, args.mva)

    if not os.path.exists(analyzer_cfg.out_dir):
        os.makedirs(analyzer_cfg.out_dir)

    if not os.path.exists(analyzer_cfg.plot_output_path):
        os.makedirs(analyzer_cfg.plot_output_path)

    out_file = TFile( analyzer_cfg.out_dir + '/' + analyzer_cfg.root_output_name , "RECREATE")
    model = pickle.load(open(analyzer_cfg.BDT_filename, 'rb'))

    

    if args.cut: 
        analyzer_cfg.sig_names = [args.mA]
        analyzer_cfg.samp_names = analyzer_cfg.bkg_names + analyzer_cfg.sig_names + ['data']


    analyzer_cfg.Print_Config()
    ntuples = LoadNtuples(analyzer_cfg)
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)
    var_names = list(plot_cfg.var_title_map.keys())
    # print(var_names)

    if args.mva:
        var_mva = []
        for ALP_mass in analyzer_cfg.sig_names:
            var_names.append('mvaVal_'+ALP_mass)
            var_names.append('mvaVal_larger_'+ALP_mass)
            for r in ['1sigma', '1P5sigma', '2sigma', '3sigma']:
                var_mva.append('mvaVal_'+r+'_'+ALP_mass)
                var_mva.append('mvaVal_larger_'+r+'_'+ALP_mass)

        var_names = var_names + var_mva

    # get mass region
    sigma_low, sigma_hig = getMassSigma(analyzer_cfg)

    ### declare histograms

    histos = {}
    histos_sys = {}

    for var_name in var_names:
        histos[var_name] = {}
        
    for sample in analyzer_cfg.samp_names:
        histos['pho1Pt'][sample]    = TH1F('pho1Pt'    + '_' + sample, 'pho1Pt'    + '_' + sample, 25,  8., 50.)
        histos['pho1eta'][sample]    = TH1F('pho1eta'    + '_' + sample, 'pho1eta'    + '_' + sample, 20,  -3., 3.)
        histos['pho1phi'][sample]    = TH1F('pho1phi'    + '_' + sample, 'pho1phi'    + '_' + sample, 20,  -4., 4.)
        histos['pho1R9'][sample]    = TH1F('pho1R9'    + '_' + sample, 'pho1R9'    + '_' + sample, 25,  0.1, 1.)
        histos['pho1IetaIeta'][sample]    = TH1F('pho1IetaIeta'    + '_' + sample, 'pho1IetaIeta'    + '_' + sample, 15,  0., 0.06)
        histos['pho1IetaIeta55'][sample]    = TH1F('pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 15,  0., 0.06)
        histos['pho1PIso_noCorr'][sample]    = TH1F('pho1PIso_noCorr'    + '_' + sample, 'pho1PIso_noCorr'    + '_' + sample, 10, 0., 40.)
        histos['pho1CIso'][sample]    = TH1F('pho1CIso'    + '_' + sample, 'pho1CIso'    + '_' + sample, 10, 0., 0.7)
        histos['pho1NIso'][sample]    = TH1F('pho1NIso'    + '_' + sample, 'pho1NIso'    + '_' + sample, 10, 0., 3.0)
        histos['pho1HOE'][sample]    = TH1F('pho1HOE'    + '_' + sample, 'pho1HOE'    + '_' + sample, 10, 0., 0.032)
        histos['pho2Pt'][sample]    = TH1F('pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 12,  8., 30.)
        histos['pho2eta'][sample]    = TH1F('pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 20,  -3., 3.)
        histos['pho2phi'][sample]    = TH1F('pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 20,  -4., 4.)
        histos['pho2R9'][sample]    = TH1F('pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 25,  0.1, 1.)
        histos['pho2IetaIeta'][sample]    = TH1F('pho2IetaIeta'    + '_' + sample, 'pho2IetaIeta'    + '_' + sample, 15,  0., 0.06)
        histos['pho2IetaIeta55'][sample]    = TH1F('pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 15,  0., 0.06)
        histos['pho2PIso_noCorr'][sample]    = TH1F('pho2PIso_noCorr'    + '_' + sample, 'pho2PIso_noCorr'    + '_' + sample, 10, 0., 40.)
        histos['pho2CIso'][sample]    = TH1F('pho2CIso'    + '_' + sample, 'pho2CIso'    + '_' + sample, 10, 0., 0.7)
        histos['pho2NIso'][sample]    = TH1F('pho2NIso'    + '_' + sample, 'pho2NIso'    + '_' + sample, 10, 0., 3.)
        histos['pho2HOE'][sample]    = TH1F('pho2HOE'    + '_' + sample, 'pho2HOE'    + '_' + sample, 10, 0., 0.032)
        histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 20,  50., 130.)
        histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 25,  95., 180.)
        histos['H_pt'][sample]    = TH1F('H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 20,  0., 160.)
        histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 20, 0., 40.)
        histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 20, 0., 5)
        histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 20, 0., 100.)
        histos['var_dR_Za'][sample] = TH1F('var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 20, 0., 7.)
        histos['var_dR_g1Z'][sample] = TH1F('var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 20, 0., 7)
        histos['var_PtaOverMh'][sample] = TH1F('var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 25, 0., 0.75)
        histos['var_Pta'][sample] = TH1F('var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 20, 0., 60.)
        histos['var_MhMa'][sample] = TH1F('var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 20, 100., 200.)
        histos['var_MhMZ'][sample] = TH1F('var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 20, 145., 310.)
        histos['ALP_calculatedPhotonIso'][sample] = TH1F('ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 20, 0., 5.)
        # histos['param'][sample] = TH1F('param' + '_' + sample, 'param' + '_' + sample, 25, -0.3, 0.6)

        if args.mva:
            for ALP_mass in analyzer_cfg.sig_names:
                #mva_bin = np.hstack((np.arange(0.0, 0.95, 0.02),np.arange(0.95,1.0,0.005)))
                #histos['mvaVal_'+ALP_mass][sample]    = TH1F('mvaVal_'+ALP_mass    + '_' + sample, 'mvaVal_'+ALP_mass    + '_' + sample, mva_bin.shape[0]-1,  mva_bin)
                histos['mvaVal_'+ALP_mass][sample]           = TH1F('mvaVal_'+ALP_mass           + '_' + sample, 'mvaVal_'+ALP_mass    + '_' + sample, 240, -0.1, 1.1)
                histos['mvaVal_1sigma_'+ALP_mass][sample]    = TH1F('mvaVal_1sigma_'+ALP_mass    + '_' + sample, 'mvaVal_1sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_1P5sigma_'+ALP_mass][sample]  = TH1F('mvaVal_1P5sigma_'+ALP_mass  + '_' + sample, 'mvaVal_1P5sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_2sigma_'+ALP_mass][sample]    = TH1F('mvaVal_2sigma_'+ALP_mass    + '_' + sample, 'mvaVal_2sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_3sigma_'+ALP_mass][sample]    = TH1F('mvaVal_3sigma_'+ALP_mass    + '_' + sample, 'mvaVal_3sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)

                histos['mvaVal_larger_'+ALP_mass][sample]           = TH1F('mvaVal_larger_'+ALP_mass           + '_' + sample, 'mvaVal_larger_'+ALP_mass    + '_' + sample, 10, 0, 1.0)
                histos['mvaVal_larger_1sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_1sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_1sigma_'+ALP_mass    + '_' + sample, 10,  0., 1.0)
                histos['mvaVal_larger_1P5sigma_'+ALP_mass][sample]  = TH1F('mvaVal_larger_1P5sigma_'+ALP_mass  + '_' + sample, 'mvaVal_larger_1P5sigma_'+ALP_mass    + '_' + sample, 10,  0., 1.0)
                histos['mvaVal_larger_2sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_2sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_2sigma_'+ALP_mass    + '_' + sample, 10,  0., 1.0)
                histos['mvaVal_larger_3sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_3sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_3sigma_'+ALP_mass    + '_' + sample, 10,  0., 1.0)


   

    for var_name in var_names:
        histos_sys[var_name] = {}
        for sample in analyzer_cfg.samp_names:
            histos_sys[var_name][sample] = {}
            for sys in analyzer_cfg.sys_names:
                histos_sys[var_name][sample][sys] = copy.deepcopy(histos[var_name][sample])
                histos_sys[var_name][sample][sys].SetNameTitle(var_name+'_'+sample+'_'+sys, var_name+'_'+sample+'_'+sys)

    ### loop over samples and events
    mass_list = {'M1':1.0, 'M2':2.0, 'M3':3.0, 'M4':4.0, 'M5':5.0, 'M6':6.0, 'M7':7.0, 'M8':8.0, 'M9':9.0, 'M10':10.0, 'M15':15.0, 'M20':20.0, 'M25':25.0, 'M30':30.0}
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        print('\n\nOn sample: %s' %sample)
        print('total events: %d' %ntup.GetEntries())

        for iEvt in range( ntup.GetEntries() ):
    
            ntup.GetEvent(iEvt)
            # if (iEvt == 100): break


            if (iEvt % 100000 == 1):
                print("looking at event %d" %iEvt)

            
            if args.ele:
                if abs(ntup.l1_id) == 13: 
                    continue
            if args.mu:
                if abs(ntup.l1_id) == 11: 
                    continue
            

            weight = ntup.factor * ntup.pho1SFs * ntup.pho2SFs

            if (ntup.H_m > -90):
                #if ntup.dR_pho < 0.02: continue
                #if not ntup.passEleVeto: continue
                
                if not ntup.passChaHadIso: continue
                if not ntup.passNeuHadIso: continue
                if not ntup.passdR_gl: continue
                if not ntup.passHOverE: continue
                #if not ntup.pho1passCutBasedIDTight or ntup.pho2passCutBasedIDTight: continue
                #if (ntup.pho1IetaIeta > 0.00996) or (ntup.pho2IetaIeta > 0.00996): continue
                #if (ntup.pho1PIso_noCorr > (0.317 + 0.01512*ntup.pho1Pt + 0.00002259*ntup.pho1Pt*ntup.pho1Pt)) or (ntup.pho1PIso_noCorr > (0.317 + 0.01512*ntup.pho2Pt + 0.00002259*ntup.pho2Pt*ntup.pho2Pt)): continue
                if ntup.H_m>180. or ntup.H_m<95.: continue

                if  args.region == 1 and (ntup.H_m>135. or ntup.H_m<115.): continue
                if  args.region == 2 and (ntup.H_m<135. and ntup.H_m>115.): continue


                #if abs(ntup.l1_id)!=13: continue

                #print 'pho1: pt:' + str(ntup.pho1Pt) + ', eta: ' + str(ntup.pho1eta) + ', SFs: ' + str(pho1_SFs)
                MVA_value = {}
                if args.mva:
                    
                    for ALP_mass in analyzer_cfg.sig_names:

                        if args.cut and ALP_mass != args.mA: continue

                        param = (ntup.ALP_m - mass_list[ALP_mass])/ntup.H_m

                        # MVA_list = [ntup.pho1Pt, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt] # v_1
                        MVA_list = [ntup.pho1Pt, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr, ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt, param] # v_1
                        #MVA_list = [ntup.pho1Pt, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr,ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt, ntup.l1_pt, ntup.l2_pt, ntup.Z_m, param] # v_2
                        #MVA_list = [ntup.pho1Pt, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr,ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt, ntup.l1_pt, ntup.l2_pt, ntup.Z_m, param]
                        #MVA_list = [ntup.pho1Pt, ntup.pho1eta, ntup.pho1phi, ntup.pho1R9, ntup.pho1IetaIeta55, ntup.pho1PIso_noCorr ,ntup.pho2Pt, ntup.pho2eta, ntup.pho2phi, ntup.pho2R9, ntup.pho2IetaIeta55,ntup.pho2PIso_noCorr,ntup.ALP_calculatedPhotonIso, ntup.var_dR_Za, ntup.var_dR_g1g2, ntup.var_dR_g1Z, ntup.var_PtaOverMh, ntup.H_pt, param]
                        MVA_value[ALP_mass] = model.predict_proba([MVA_list])[:, 1]


                    if args.cut:
                        if MVA_value[args.mA] < analyzer_cfg.mvaCut[args.mA]: continue
                        #if MVA_value[args.mA] < args.cutVal: continue

                    for ALP_mass in analyzer_cfg.sig_names:
                        if args.cut and ALP_mass != args.mA: continue

                        if args.blind:
                            #if not (sample == 'data' and MVA_value[args.mA] > mvaCut[args.mA]):
                            if not (sample == 'data' and MVA_value[ALP_mass] > 0.95):
                                histos['mvaVal_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                                histos['mvaVal_larger_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                        else:
                            histos['mvaVal_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            
                        if ntup.H_m<(125.+sigma_hig[ALP_mass]) and ntup.H_m>(125.+sigma_low[ALP_mass]): 
                            histos['mvaVal_1sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_1sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )

                        if ntup.H_m<(125.+sigma_hig[ALP_mass]*1.5) and ntup.H_m>(125.+sigma_low[ALP_mass]*1.5): 
                            histos['mvaVal_1P5sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_1P5sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )

                        if ntup.H_m<(125.+sigma_hig[ALP_mass]*2.) and ntup.H_m>(125.+sigma_low[ALP_mass]*2.): 
                            histos['mvaVal_2sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_2sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )

                        if ntup.H_m<(125.+sigma_hig[ALP_mass]*3.) and ntup.H_m>(125.+sigma_low[ALP_mass]*3.): 
                            histos['mvaVal_3sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                            histos['mvaVal_larger_3sigma_'+ALP_mass][sample].Fill( MVA_value[ALP_mass], weight )
                        

                        # print(f"125.+sigma_hig[ALP_mass]: {125.+sigma_hig[ALP_mass]}")
                        # print(f"ntup.H_m: {ntup.H_m}")
                        # print(f"125.+sigma_low[ALP_mass]: {125.+sigma_low[ALP_mass]}")

                        #if sample == 'DYJetsToLL':
                            #histos['mvaVal_bkg'].Fill(mass_list[ALP_mass], MVA_value[ALP_mass], weight)
                
                var_map = {'Z_m':ntup.Z_m, 'H_m':ntup.H_m, 'ALP_m':ntup.ALP_m,'pho1Pt':ntup.pho1Pt, 'pho1eta':ntup.pho1eta, 'pho1phi':ntup.pho1phi, 'pho1R9':ntup.pho1R9, 'pho1IetaIeta':ntup.pho1IetaIeta, 'pho1IetaIeta55':ntup.pho1IetaIeta55,'pho1PIso_noCorr':ntup.pho1PIso_noCorr, 'pho1CIso':ntup.pho1CIso, 'pho1NIso':ntup.pho1NIso, 'pho1HOE':ntup.pho1HOE, 'pho2Pt':ntup.pho2Pt, 'pho2eta':ntup.pho2eta, 'pho2phi':ntup.pho2phi, 'pho2R9':ntup.pho2R9, 'pho2IetaIeta':ntup.pho2IetaIeta, 'pho2IetaIeta55':ntup.pho2IetaIeta55,'pho2PIso_noCorr':ntup.pho2PIso_noCorr, 'pho2CIso':ntup.pho2CIso, 'pho2NIso':ntup.pho2NIso, 'pho2HOE':ntup.pho2HOE,'ALP_calculatedPhotonIso':ntup.ALP_calculatedPhotonIso, 'var_dR_Za':ntup.var_dR_Za, 'var_dR_g1g2':ntup.var_dR_g1g2, 'var_dR_g1Z':ntup.var_dR_g1Z, 'var_PtaOverMh':ntup.var_PtaOverMh, 'var_Pta':ntup.var_Pta, 'var_MhMZ':ntup.var_MhMZ, 'H_pt':ntup.H_pt, 'var_PtaOverMa':ntup.var_PtaOverMa, 'var_MhMa':ntup.var_MhMa}
                
                if args.mva:
                    var_map_mva = {}
                    for ALP_mass in analyzer_cfg.sig_names:
                        var_map_mva['mvaVal_'+ALP_mass] = MVA_value[ALP_mass]
                        var_map_mva['mvaVal_larger_'+ALP_mass] = MVA_value[ALP_mass]
                        for r in ['1sigma', '1P5sigma', '2sigma', '3sigma']:
                            var_map_mva['mvaVal_'+r+'_'+ALP_mass] = MVA_value[ALP_mass]
                            var_map_mva['mvaVal_larger_'+r+'_'+ALP_mass] = MVA_value[ALP_mass]

                            
                    var_map.update(var_map_mva)

                histos['pho1Pt'][sample].Fill( ntup.pho1Pt, weight )
                histos['pho1eta'][sample].Fill( ntup.pho1eta, weight )
                histos['pho1phi'][sample].Fill( ntup.pho1phi, weight )
                histos['pho1R9'][sample].Fill( ntup.pho1R9, weight )
                histos['pho1IetaIeta'][sample].Fill( ntup.pho1IetaIeta, weight )
                histos['pho1IetaIeta55'][sample].Fill( ntup.pho1IetaIeta55, weight )
                histos['pho1PIso_noCorr'][sample].Fill( ntup.pho1PIso_noCorr, weight )
                histos['pho2Pt'][sample].Fill( ntup.pho2Pt, weight )
                histos['pho2eta'][sample].Fill( ntup.pho2eta, weight )
                histos['pho2phi'][sample].Fill( ntup.pho2phi, weight )
                histos['pho2R9'][sample].Fill( ntup.pho2R9, weight )
                histos['pho2IetaIeta'][sample].Fill( ntup.pho2IetaIeta, weight )
                histos['pho2IetaIeta55'][sample].Fill( ntup.pho2IetaIeta55, weight )
                histos['pho2PIso_noCorr'][sample].Fill( ntup.pho2PIso_noCorr, weight )

                histos['pho1CIso'][sample].Fill( ntup.pho1CIso, weight)
                histos['pho1NIso'][sample].Fill( ntup.pho1NIso, weight)
                histos['pho1HOE'][sample].Fill( ntup.pho1HOE, weight)
                histos['pho2CIso'][sample].Fill( ntup.pho2CIso, weight)
                histos['pho2NIso'][sample].Fill( ntup.pho2NIso, weight)
                histos['pho2HOE'][sample].Fill( ntup.pho2HOE, weight)

                if args.blind:
                    if not (sample == 'data' and (ntup.H_m<135. and ntup.H_m>115.)): 
                        histos['H_m'][sample].Fill( ntup.H_m, weight )
                else:        
                    histos['H_m'][sample].Fill( ntup.H_m, weight )

                histos['H_pt'][sample].Fill( ntup.H_pt, weight )
                histos['ALP_m'][sample].Fill( ntup.ALP_m, weight )
                histos['Z_m'][sample].Fill( ntup.Z_m, weight )

                histos['var_dR_Za'][sample].Fill( ntup.var_dR_Za, weight )
                histos['var_dR_g1g2'][sample].Fill( ntup.var_dR_g1g2, weight )
                histos['var_dR_g1Z'][sample].Fill( ntup.var_dR_g1Z, weight )
                histos['var_PtaOverMa'][sample].Fill( ntup.var_PtaOverMa, weight )
                histos['var_PtaOverMh'][sample].Fill( ntup.var_PtaOverMh, weight )
                histos['var_Pta'][sample].Fill( ntup.var_Pta, weight )
                histos['var_MhMa'][sample].Fill( ntup.var_MhMa, weight )
                histos['var_MhMZ'][sample].Fill( ntup.var_MhMZ, weight )
                histos['ALP_calculatedPhotonIso'][sample].Fill( ntup.ALP_calculatedPhotonIso, weight )

                param_val = {}
                if sample in analyzer_cfg.sig_names:
                    param_val['param'] = (ntup.ALP_m - mass_list[sample])/ntup.H_m
                else:
                    mass_random = random.choice(list(mass_list.values()))
                    param_val['param'] = (ntup.ALP_m - mass_random)/ntup.H_m
                
                var_map.update(param_val)

                # histos['param'][sample].Fill( param_val['param'], weight )
                    

                for sys_name in analyzer_cfg.sys_names:

                    if sys_name =='CMS_eff_g_up':
                        #weight = ntup.factor * (ntup.pho1SFs+ntup.pho1SFs_sys) * (ntup.pho2SFs+ntup.pho2SFs_sys)
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * (ntup.pho1SFs+ntup.pho1SFs_sys) * (ntup.pho2SFs+ntup.pho2SFs_sys)
                    elif sys_name =='CMS_eff_g_dn':
                        #weight = ntup.factor * (ntup.pho1SFs-ntup.pho1SFs_sys) * (ntup.pho2SFs-ntup.pho2SFs_sys)
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeight * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * (ntup.pho1SFs-ntup.pho1SFs_sys) * (ntup.pho2SFs-ntup.pho2SFs_sys)
                    elif sys_name =='CMS_pileup_up':
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeightUp * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                    elif sys_name =='CMS_pileup_dn':
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeightDn * ntup.l1_dataMC * ntup.l2_dataMC * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                    elif sys_name =='CMS_eff_lep_up':
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeight * (ntup.l1_dataMC+ntup.l1_dataMCErr) * (ntup.l2_dataMC+ntup.l2_dataMCErr) * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs
                    elif sys_name =='CMS_eff_lep_dn':
                        weight_sys = ntup.event_genWeight * ntup.event_pileupWeight * (ntup.l1_dataMC-ntup.l1_dataMCErr) * (ntup.l2_dataMC-ntup.l2_dataMCErr) * ntup.event_weight * ntup.pho1SFs * ntup.pho2SFs

                    for var in var_names:
                        histos_sys[var][sample][sys_name].Fill(var_map[var], weight_sys)

                



        ## End of for iEvt in range( ntup.GetEntries() )
    ## End of for sample in analyzer_cfg.samp_names

    ### save raw histograms
    raw_dir = out_file.mkdir('raw_plots')
    raw_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            plot_cfg.SetHistStyles(histos[var_name][sample], sample)
            histos[var_name][sample].Write()

    #### save sys histograms
    sys_dir = out_file.mkdir('sys_dir')
    sys_dir.cd()
    for var_name in var_names:
        for sample in analyzer_cfg.samp_names:
            for sys in analyzer_cfg.sys_names:
                plot_cfg.SetHistStyles(histos_sys[var_name][sample][sys], sample)
                histos_sys[var_name][sample][sys].Write()

    ### save stack plots and make ratio plots
    out_file.cd()
    lumi_label = MakeLumiLabel(plot_cfg.lumi)
    cms_label  = MakeCMSDASLabel()

    scaled_sig = {}
    for var_name in var_names:
        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        #scaled_sig = 0
        #ratio_plot = 0
        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['data'], stacks['all'].GetStack().Last(), var_name)
        legend = MakeLegend(plot_cfg, histos[var_name], scaled_sig)

        #### uncertainty graph
        total_unc = Total_Unc(stacks['bkg'], histos_sys[var_name], analyzer_cfg)

        if args.ln:
            canv_log = CreateCanvas(var_name+'_log')
            DrawOnCanv(canv_log, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, args.cut, args.mA, logY=True)

            canv_log.Write()
            SaveCanvPic(canv_log, analyzer_cfg.plot_output_path, var_name+'_log')
        else:
            canv = CreateCanvas(var_name)
            DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, args.cut, args.mA, logY=False)

            canv.Write()
            SaveCanvPic(canv, analyzer_cfg.plot_output_path, var_name)


    
    print('\n\n')
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()

    print('Done')


main()
