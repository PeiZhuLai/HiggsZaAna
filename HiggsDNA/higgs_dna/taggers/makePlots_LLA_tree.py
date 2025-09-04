#! /usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default=["Sync_1031_2018_ttH_v2.root"], nargs='*', help="List of input ggNtuplizer files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-m", "--maxevents", dest="maxevents", type=int, default=-1, help="Maximum number events to loop over")
parser.add_argument("-t", "--ttree", dest="ttree", default="Ana/passedEvents", help="TTree Name")
parser.add_argument("-xs", "--cross_section", dest="cross_section", default="1.0", help="the cross section of samples")
parser.add_argument("-L", "--Lumi", dest="Lumi", default="35.9", help="the luminosities to normalized")
parser.add_argument("-N", "--NEvent", dest="NEvent", default="0", help="number of events")
parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
args = parser.parse_args()

########################
import numpy as np
import ROOT
import os
import math
###########################
from deltaR import *
from array import array
import SystematicUnc_ALP as Sys
# from calculateMVA import bookMVA, calMVA

from sklearn.metrics import accuracy_score, roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import learning_curve
from sklearn.model_selection import ShuffleSplit

import seaborn as sns
from scipy.stats import ks_2samp
from scipy import stats

import xgboost as xgb
import random
import optuna

import pickle
# from sklearn.preprocessing import LabelEncoder

# Read in model saved from previous running of BDT


########## photon SFs
bins_eta = [-2.5,-2.0,-1.566,-1.4442, -0.8, 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5]
bins_pt = [10,15,20,35,50]

# photon ID
chCut_EB = 0.65
chCut_EE = 0.517
neuCut_EB = [0.317, 0.01512, 0.00002259]
neuCut_EE = [2.716, 0.0117, 0.000023]
hoeCut_EB = 0.02148
hoeCut_EE = 0.0321

sigieieCut_EB = 0.00996
sigieieCut_EE = 0.0271
phoCut_EB = [2.044, 0.004017]
phoCut_EE = [3.032, 0.0037]

if args.year == '2016':
    f_SFs = '/publicfs/cms/user/wangzebing/ALP/Plot/SFs_UL/egammaEffi_2016_postVFP.root'

elif args.year == '-2016':
    f_SFs = '/publicfs/cms/user/wangzebing/ALP/Plot/SFs_UL/egammaEffi_2016_preVFP.root'

elif args.year == '2017':
    f_SFs = '/publicfs/cms/user/wangzebing/ALP/Plot/SFs_UL/egammaEffi_2017.root'

elif args.year == '2018':
    f_SFs = '/publicfs/cms/user/wangzebing/ALP/Plot/SFs_UL/egammaEffi_2018.root'
    
else:
    print("do not include at 2016/2017/2018")
    exit(0)

# load the model from disk
BDT_filename = "/publicfs/cms/user/laipeizhu/ALP/Analysis_code/train_MVA/model_Za_BDT_passedEvents.pkl"
model = pickle.load(open(BDT_filename, 'rb'))


SFs = np.zeros((len(bins_eta)-1,len(bins_pt)-1))
SFs_statData = np.zeros((len(bins_eta)-1,len(bins_pt)-1))
SFs_statMC = np.zeros((len(bins_eta)-1,len(bins_pt)-1))
SFs_altBkgModel = np.zeros((len(bins_eta)-1,len(bins_pt)-1))
SFs_altSignalModel = np.zeros((len(bins_eta)-1,len(bins_pt)-1))
SFs_altMCEff = np.zeros((len(bins_eta)-1,len(bins_pt)-1))
SFs_altTagSelection = np.zeros((len(bins_eta)-1,len(bins_pt)-1))

files_SFs = ROOT.TFile(f_SFs)
his_SFs = files_SFs.Get('EGamma_SF2D')
his_statData = files_SFs.Get('statData')
his_statMC = files_SFs.Get('statMC')
his_altBkgModel = files_SFs.Get('altBkgModel')
his_altSignalModel = files_SFs.Get('altSignalModel')
his_altMCEff = files_SFs.Get('altMCEff')
his_altTagSelection = files_SFs.Get('altTagSelection')

for i in range(len(bins_eta)-1):
    for j in range(len(bins_pt)-1):
        #print str(i)+' '+str(j)+': '+str(his_SFs.GetBinContent(i+1,j+1))
        SFs[i][j] = his_SFs.GetBinContent(i+1,j+1)
        SFs_statData[i][j] = his_statData.GetBinContent(i+1,j+1)
        SFs_statMC[i][j] = his_statMC.GetBinContent(i+1,j+1)
        SFs_altBkgModel[i][j] = his_altBkgModel.GetBinContent(i+1,j+1)
        SFs_altSignalModel[i][j] = his_altSignalModel.GetBinContent(i+1,j+1)
        SFs_altMCEff[i][j] = his_altMCEff.GetBinContent(i+1,j+1)
        SFs_altTagSelection[i][j] = his_altTagSelection.GetBinContent(i+1,j+1)



###########################
if os.path.isfile('~/.rootlogon.C'): ROOT.gROOT.Macro(os.path.expanduser('~/.rootlogon.C'))
ROOT.gROOT.SetBatch()
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(000000)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.UseCurrentStyle()

sw = ROOT.TStopwatch()
sw.Start()

# Input ggNtuple
tchain = ROOT.TChain(args.ttree)
for filename in args.inputfiles: tchain.Add(filename)
print('Total number of events: ' + str(tchain.GetEntries()))

# Event weights
isMC = True
if 'Run2018' in filename:
    isMC = False
elif 'Run2017' in filename:
    isMC = False
elif 'Run2016' in filename:
    isMC = False
else:
    isMC = True

# get nEvents
nEvents = 0

for filename in args.inputfiles:

    files = ROOT.TFile(filename)
    n_his = files.Ana.Get('sumWeights')
    nEvents = nEvents + n_his.GetBinContent(1)


if isMC:
    cross_section = float(args.cross_section)
    lumi = float(args.Lumi)

    if int(args.NEvent):
        weight = cross_section * lumi * 1000.0 / float(args.NEvent)
    else:
        weight = cross_section * lumi * 1000.0 / nEvents
else:
    cross_section = 1.0
    weight = 1.0



print('events weight: '+str(weight))
print('events : '+str(tchain.GetEntries()))


# Output file and any histograms we want
file_out = ROOT.TFile(args.outputfile, 'recreate')

nEvents_total = ROOT.TH1D('nEvents_total', 'nEvents_total', 2, 0, 2)
nEvents_total.SetBinContent(1,nEvents)
h_weight = ROOT.TH1D('Events_weight', 'Events_weight', 2, 0, 2)
h_weight.SetBinContent(1,weight)
h_cross_section = ROOT.TH1D('cross_section', 'cross_section', 2, 0, 2)
h_cross_section.SetBinContent(1,cross_section)

# pass triger
h_n = ROOT.TH1D('nEvents_ntuple', 'nEvents_ntuple', 2, 0, 2)
h_n_weight = ROOT.TH1D('nEvents_weight_ntuple', 'nEvents_weight_ntuple', 2, 0, 2)
h_n_trig = ROOT.TH1D('nEvents_trig', 'nEvents_trig', 2, 0, 2)



npho = ROOT.TH1D('npho', 'npho', 10, 0., 10)
nlep = ROOT.TH1D('nlep', 'nlep', 10, 0., 10)

################################################################################################
Z_e_nocut = ROOT.TH1D('Z_e_nocut', 'Z_e_nocut', 100, 0, 500)
Z_mu_nocut = ROOT.TH1D('Z_mu_nocut', 'Z_mu_nocut', 100, 0, 500)

l1_id_lIso = ROOT.TH1D('l1_id_lIso', 'l1_id_lIso', 40, -20., 20.)
Z_e_lIso = ROOT.TH1D('Z_e_lIso', 'Z_e_lIso', 100, 0, 500)
Z_mu_lIso = ROOT.TH1D('Z_mu_lIso', 'Z_mu_lIso', 100, 0, 500)

Z_e_lIso_lTight = ROOT.TH1D('Z_e_lIso_lTight', 'Z_e_lIso_lTight', 100, 0, 500)
Z_mu_lIso_lTight = ROOT.TH1D('Z_mu_lIso_lTight', 'Z_mu_lIso_lTight', 100, 0, 500)

Z_50 = ROOT.TH1D('Z_50', 'Z_50', 100, 0, 500)

################################################################################################

################################################################################################




h_n.SetStats(1)
h_n_trig.SetStats(1)

Z_e_nocut.SetStats(1)
Z_mu_nocut.SetStats(1)
Z_e_lIso.SetStats(1)
Z_mu_lIso.SetStats(1)
Z_e_lIso_lTight.SetStats(1)
Z_mu_lIso_lTight.SetStats(1)

Z_50.SetStats(1)

################################################################################################

# photon cut tree
Run = array('l',[0])
LumiSect = array('l',[0])
Event = array('l',[0])
factor = array('f',[0.])
event_weight = array('f',[0.])
event_genWeight = array('f',[-1.])
event_pileupWeight = array('f',[-1.])
event_pileupWeightUp = array('f',[-1.])
event_pileupWeightDn = array('f',[-1.])
event_dataMCWeight = array('f',[-1.])
event_SelectEff_weight = array('f',[0.])
event_passedTrig = array('l',[0])
event_qcdweight = array('f',[0.])
event_nnloweight = array('f',[0.])
event_pdfweight_up = array('f',[0.])
event_pdfweight_dn = array('f',[0.])
event_qcdweight_up = array('f',[0.])
event_qcdweight_dn = array('f',[0.])

Z_m = array('f',[0.])
H_m = array('f',[0.])
ALP_m = array('f',[0.])
H_pt = array('f',[0.])

Z_noFSR_m = array('f',[0.])
H_noFSR_m = array('f',[0.])

l1_pt = array('f',[0.])
l1_eta = array('f',[0.])
l1_phi = array('f',[0.])
l1_mass = array('f',[0.])
l1_id = array('i',[0])
l1_dataMC = array('f',[0.])
l1_dataMCErr = array('f',[0.])
l1_scaleup = array('f',[0.])
l1_scaledn = array('f',[0.])
l1_smearup = array('f',[0.])
l1_smeardn = array('f',[0.])
l1_EPostCorr = array('f',[0.])
l1_EPreCorr = array('f',[0.])
l1_noFSR_pt = array('f',[0.])
l1_noFSR_eta = array('f',[0.])
l1_noFSR_phi = array('f',[0.])
l1_noFSR_mass = array('f',[0.])
l1_matchedR03_PdgId = array('f',[0.])
l1_matchedR03_MomId = array('f',[0.])
l1_matchedR03_MomMomId = array('f',[0.])



l2_pt = array('f',[0.])
l2_eta = array('f',[0.])
l2_phi = array('f',[0.])
l2_mass = array('f',[0.])
l2_id = array('i',[0])
l2_dataMCErr = array('f',[0.])
l2_scaleup = array('f',[0.])
l2_scaledn = array('f',[0.])
l2_smearup = array('f',[0.])
l2_smeardn = array('f',[0.])
l2_dataMC = array('f',[0.])
l2_EPostCorr = array('f',[0.])
l2_EPreCorr = array('f',[0.])
l2_noFSR_pt = array('f',[0.])
l2_noFSR_eta = array('f',[0.])
l2_noFSR_phi = array('f',[0.])
l2_noFSR_mass = array('f',[0.])
l2_matchedR03_PdgId = array('f',[0.])
l2_matchedR03_MomId = array('f',[0.])
l2_matchedR03_MomMomId = array('f',[0.])

# photon var
pho1eta = array('f',[0.])
pho1Pt = array('f',[0.])
pho1phi = array('f',[0.])
pho1R9 = array('f',[0.])
pho1IetaIeta = array('f',[0.])
pho1IetaIeta55 = array('f',[0.])
pho1HOE = array('f',[0.])
pho1CIso = array('f',[0.])
pho1NIso = array('f',[0.])
pho1PIso = array('f',[0.])
pho1PIso_noCorr = array('f',[0.])
pho1SFs = array('f',[0.])
pho1SFs_dR0P1 = array('f',[0.])
pho1SFs_dR0P15 = array('f',[0.])
pho1SFs_dR0P2 = array('f',[0.])
pho1SFs_dR0P25 = array('f',[0.])
pho1SFs_dR0P3 = array('f',[0.])
pho1SFs_sys = array('f',[0.])
pho1SFs_sys_replaced = array('f',[0.])
pho1scaleup = array('f',[0.])
pho1scaledn = array('f',[0.])
pho1smearup = array('f',[0.])
pho1smeardn = array('f',[0.])
pho1EPostCorr = array('f',[0.])
pho1ShowerShapeSys = array('f',[0.])
pho1trackIso = array('f',[0.])

pho2eta = array('f',[0.])
pho2Pt = array('f',[0.])
pho2phi = array('f',[0.])
pho2R9 = array('f',[0.])
pho2IetaIeta = array('f',[0.])
pho2IetaIeta55 = array('f',[0.])
pho2HOE = array('f',[0.])
pho2CIso = array('f',[0.])
pho2NIso = array('f',[0.])
pho2PIso = array('f',[0.])
pho2PIso_noCorr = array('f',[0.])
pho2SFs = array('f',[0.])
pho2SFs_dR0P1 = array('f',[0.])
pho2SFs_dR0P15 = array('f',[0.])
pho2SFs_dR0P2 = array('f',[0.])
pho2SFs_dR0P25 = array('f',[0.])
pho2SFs_dR0P3 = array('f',[0.])
pho2SFs_sys = array('f',[0.])
pho2SFs_sys_replaced = array('f',[0.])
pho2scaleup = array('f',[0.])
pho2scaledn = array('f',[0.])
pho2smearup = array('f',[0.])
pho2smeardn = array('f',[0.])
pho2EPostCorr = array('f',[0.])
pho2ShowerShapeSys = array('f',[0.])
pho2trackIso = array('f',[0.])

dR_g1l1 = array('f',[-1.])
dR_g1l2 = array('f',[-1.])
dR_g2l1 = array('f',[-1.])
dR_g2l2 = array('f',[-1.])

dR_pho = array('f',[0.])

var_dR_Za = array('f',[-1.])
var_dR_g1g2 = array('f',[-1.])
var_dR_g1Z = array('f',[-1.])
var_dEta_g1Z = array('f',[-1.])
var_dPhi_g1Z = array('f',[-1.])
var_PtaOverMa = array('f',[-1.])
var_PtaOverMh = array('f',[-1.])
var_Pta = array('f',[-1.])
var_MhMa = array('f',[-1.])
var_MhMZ = array('f',[-1.])

ALP_calculatedPhotonIso = array('f',[-1.])

# photon cut

passEleVeto = array('f',[0.])
passIeIe = array('f',[0.])
passHOverE = array('f',[0.])
passChaHadIso = array('f',[0.])
passNeuHadIso = array('f',[0.])
passPhoIso = array('f',[0.])
passdR_gl = array('f',[0.])
passdR_gg = array('f',[0.])
passH_m = array('f',[0.])
passBDT = array('f',[0.])
Val_BDT = array('f',[0.])

pho1passCutBasedIDTight = array('f',[-90.])
pho2passCutBasedIDTight = array('f',[-90.])
pho1photonCutBasedIDMedium = array('f',[-90.])
pho2photonCutBasedIDMedium = array('f',[-90.])

pho1_matchedR03_PdgId = array('f',[-90.])
pho1_matchedR03_MomId = array('f',[-90.])
pho1_matchedR03_MomMomId = array('f',[-90.])
pho1_matchedR = array('f',[-90.])
pho1_matchedPt = array('f',[-90.])

pho2_matchedR03_PdgId = array('f',[-90.])
pho2_matchedR03_MomId = array('f',[-90.])
pho2_matchedR03_MomMomId = array('f',[-90.])
pho2_matchedR = array('f',[-90.])
pho2_matchedPt = array('f',[-90.])
################################################################################################

passedEvents = ROOT.TTree("passedEvents","passedEvents")

################################################################################################
passedEvents.Branch("Run",Run,"Run/L")
passedEvents.Branch("LumiSect",LumiSect,"LumiSect/L")
passedEvents.Branch("Event",Event,"Event/L")
passedEvents.Branch("factor",factor,"factor/F")
passedEvents.Branch("event_weight",event_weight,"event_weight/F")
passedEvents.Branch("event_genWeight",event_genWeight,"event_genWeight/F")
passedEvents.Branch("event_pileupWeight",event_pileupWeight,"event_pileupWeight/F")
passedEvents.Branch("event_pileupWeightUp",event_pileupWeightUp,"event_pileupWeightUp/F")
passedEvents.Branch("event_pileupWeightDn",event_pileupWeightDn,"event_pileupWeightDn/F")
passedEvents.Branch("event_dataMCWeight",event_dataMCWeight,"event_dataMCWeight/F")
passedEvents.Branch("event_SelectEff_weight",event_SelectEff_weight,"event_SelectEff_weight/F")
passedEvents.Branch("event_passedTrig",event_passedTrig,"event_passedTrig/L")
passedEvents.Branch("event_qcdweight",event_qcdweight,"event_qcdweight/F")
passedEvents.Branch("event_nnloweight",event_nnloweight,"event_nnloweight/F")
passedEvents.Branch("event_pdfweight_up",event_pdfweight_up,"event_pdfweight_up/F")
passedEvents.Branch("event_pdfweight_dn",event_pdfweight_dn,"event_pdfweight_dn/F")
passedEvents.Branch("event_qcdweight_up",event_qcdweight_up,"event_qcdweight_up/F")
passedEvents.Branch("event_qcdweight_dn",event_qcdweight_dn,"event_qcdweight_dn/F")

passedEvents.Branch("Z_m",Z_m,"Z_m/F")
passedEvents.Branch("H_m",H_m,"H_m/F")
passedEvents.Branch("ALP_m",ALP_m,"ALP_m/F")
passedEvents.Branch("H_pt",H_pt,"H_pt/F")

passedEvents.Branch("Z_noFSR_m",Z_noFSR_m,"Z_noFSR_m/F")
passedEvents.Branch("H_noFSR_m",H_noFSR_m,"H_noFSR_m/F")


passedEvents.Branch("l1_pt",l1_pt,"l1_pt/F")
passedEvents.Branch("l1_eta",l1_eta,"l1_eta/F")
passedEvents.Branch("l1_phi",l1_phi,"l1_phi/F")
passedEvents.Branch("l1_mass",l1_mass,"l1_mass/F")
passedEvents.Branch("l1_id",l1_id,"l1_id/I")
passedEvents.Branch("l1_dataMC",l1_dataMC,"l1_dataMC/F")
passedEvents.Branch("l1_dataMCErr",l1_dataMCErr,"l1_dataMCErr/F")
passedEvents.Branch("l1_scaleup",l1_scaleup,"l1_scaleup/F")
passedEvents.Branch("l1_scaledn",l1_scaledn,"l1_scaledn/F")
passedEvents.Branch("l1_smearup",l1_smearup,"l1_smearup/F")
passedEvents.Branch("l1_smeardn",l1_smeardn,"l1_smeardn/F")
passedEvents.Branch("l1_EPostCorr",l1_EPostCorr,"l1_EPostCorr/F")
passedEvents.Branch("l1_EPreCorr",l1_EPreCorr,"l1_EPreCorr/F")
passedEvents.Branch("l1_noFSR_pt",l1_noFSR_pt,"l1_noFSR_pt/F")
passedEvents.Branch("l1_noFSR_eta",l1_noFSR_eta,"l1_noFSR_eta/F")
passedEvents.Branch("l1_noFSR_phi",l1_noFSR_phi,"l1_noFSR_phi/F")
passedEvents.Branch("l1_noFSR_mass",l1_noFSR_mass,"l1_noFSR_mass/F")
passedEvents.Branch("l1_matchedR03_PdgId",l1_matchedR03_PdgId,"l1_matchedR03_PdgId/F")
passedEvents.Branch("l1_matchedR03_MomId",l1_matchedR03_MomId,"l1_matchedR03_MomId/F")
passedEvents.Branch("l1_matchedR03_MomMomId",l1_matchedR03_MomMomId,"l1_matchedR03_MomMomId/F")


passedEvents.Branch("l2_pt",l2_pt,"l2_pt/F")
passedEvents.Branch("l2_eta",l2_eta,"l2_eta/F")
passedEvents.Branch("l2_phi",l2_phi,"l2_phi/F")
passedEvents.Branch("l2_mass",l2_mass,"l2_mass/F")
passedEvents.Branch("l2_id",l2_id,"l2_id/I")
passedEvents.Branch("l2_dataMC",l2_dataMC,"l2_dataMC/F")
passedEvents.Branch("l2_dataMCErr",l2_dataMCErr,"l2_dataMCErr/F")
passedEvents.Branch("l2_scaleup",l2_scaleup,"l2_scaleup/F")
passedEvents.Branch("l2_scaledn",l2_scaledn,"l2_scaledn/F")
passedEvents.Branch("l2_smearup",l2_smearup,"l2_smearup/F")
passedEvents.Branch("l2_smeardn",l2_smeardn,"l2_smeardn/F")
passedEvents.Branch("l2_EPostCorr",l2_EPostCorr,"l2_EPostCorr/F")
passedEvents.Branch("l2_EPreCorr",l2_EPreCorr,"l2_EPreCorr/F")
passedEvents.Branch("l2_matchedR03_PdgId",l2_matchedR03_PdgId,"l2_matchedR03_PdgId/F")
passedEvents.Branch("l2_matchedR03_MomId",l2_matchedR03_MomId,"l2_matchedR03_MomId/F")
passedEvents.Branch("l2_matchedR03_MomMomId",l2_matchedR03_MomMomId,"l2_matchedR03_MomMomId/F")


passedEvents.Branch("pho1eta",pho1eta,"pho1eta/F")
passedEvents.Branch("pho1Pt",pho1Pt,"pho1Pt/F")
passedEvents.Branch("pho1phi",pho1phi,"pho1phi/F")
passedEvents.Branch("pho1R9",pho1R9,"pho1R9/F")
passedEvents.Branch("pho1IetaIeta",pho1IetaIeta,"pho1IetaIeta/F")
passedEvents.Branch("pho1IetaIeta55",pho1IetaIeta55,"pho1IetaIeta55/F")
passedEvents.Branch("pho1HOE",pho1HOE,"pho1HOE/F")
passedEvents.Branch("pho1CIso",pho1CIso,"pho1CIso/F")
passedEvents.Branch("pho1NIso",pho1NIso,"pho1NIso/F")
passedEvents.Branch("pho1PIso",pho1PIso,"pho1PIso/F")
passedEvents.Branch("pho1PIso_noCorr",pho1PIso_noCorr,"pho1PIso_noCorr/F")
passedEvents.Branch("pho1SFs",pho1SFs,"pho1SFs/F")
passedEvents.Branch("pho1SFs_dR0P1",pho1SFs_dR0P1,"pho1SFs_dR0P1/F")
passedEvents.Branch("pho1SFs_dR0P15",pho1SFs_dR0P15,"pho1SFs_dR0P15/F")
passedEvents.Branch("pho1SFs_dR0P2",pho1SFs_dR0P2,"pho1SFs_dR0P2/F")
passedEvents.Branch("pho1SFs_dR0P25",pho1SFs_dR0P25,"pho1SFs_dR0P25/F")
passedEvents.Branch("pho1SFs_dR0P3",pho1SFs_dR0P3,"pho1SFs_dR0P3/F")
passedEvents.Branch("pho1SFs_sys",pho1SFs_sys,"pho1SFs_sys/F")
passedEvents.Branch("pho1SFs_sys_replaced",pho1SFs_sys_replaced,"pho1SFs_sys_replaced/F")
passedEvents.Branch("pho1scaleup",pho1scaleup,"pho1scaleup/F")
passedEvents.Branch("pho1scaledn",pho1scaledn,"pho1scaledn/F")
passedEvents.Branch("pho1smearup",pho1smearup,"pho1smearup/F")
passedEvents.Branch("pho1smeardn",pho1smeardn,"pho1smeardn/F")
passedEvents.Branch("pho1EPostCorr",pho1EPostCorr,"pho1EPostCorr/F")
passedEvents.Branch("pho1ShowerShapeSys",pho1ShowerShapeSys,"pho1ShowerShapeSys/F")
passedEvents.Branch("pho1trackIso",pho1trackIso,"pho1trackIso/F")

passedEvents.Branch("pho2eta",pho2eta,"pho2eta/F")
passedEvents.Branch("pho2Pt",pho2Pt,"pho2Pt/F")
passedEvents.Branch("pho2phi",pho2phi,"pho2phi/F")
passedEvents.Branch("pho2R9",pho2R9,"pho2R9/F")
passedEvents.Branch("pho2IetaIeta",pho2IetaIeta,"pho2IetaIeta/F")
passedEvents.Branch("pho2IetaIeta55",pho2IetaIeta55,"pho2IetaIeta55/F")
passedEvents.Branch("pho2HOE",pho2HOE,"pho2HOE/F")
passedEvents.Branch("pho2CIso",pho2CIso,"pho2CIso/F")
passedEvents.Branch("pho2NIso",pho2NIso,"pho2NIso/F")
passedEvents.Branch("pho2PIso",pho2PIso,"pho2PIso/F")
passedEvents.Branch("pho2PIso_noCorr",pho2PIso_noCorr,"pho2PIso_noCorr/F")
passedEvents.Branch("pho2SFs",pho2SFs,"pho2SFs/F")
passedEvents.Branch("pho2SFs_dR0P1",pho2SFs_dR0P1,"pho2SFs_dR0P1/F")
passedEvents.Branch("pho2SFs_dR0P15",pho2SFs_dR0P15,"pho2SFs_dR0P15/F")
passedEvents.Branch("pho2SFs_dR0P2",pho2SFs_dR0P2,"pho2SFs_dR0P2/F")
passedEvents.Branch("pho2SFs_dR0P25",pho2SFs_dR0P25,"pho2SFs_dR0P25/F")
passedEvents.Branch("pho2SFs_dR0P3",pho2SFs_dR0P3,"pho2SFs_dR0P3/F")
passedEvents.Branch("pho2SFs_sys",pho2SFs_sys,"pho2SFs_sys/F")
passedEvents.Branch("pho2SFs_sys_replaced",pho2SFs_sys_replaced,"pho2SFs_sys_replaced/F")
passedEvents.Branch("pho2scaleup",pho2scaleup,"pho2scaleup/F")
passedEvents.Branch("pho2scaledn",pho2scaledn,"pho2scaledn/F")
passedEvents.Branch("pho2smearup",pho2smearup,"pho2smearup/F")
passedEvents.Branch("pho2smeardn",pho2smeardn,"pho2smeardn/F")
passedEvents.Branch("pho2EPostCorr",pho2EPostCorr,"pho2EPostCorr/F")
passedEvents.Branch("pho2ShowerShapeSys",pho2ShowerShapeSys,"pho2ShowerShapeSys/F")
passedEvents.Branch("pho2trackIso",pho2trackIso,"pho2trackIso/F")

passedEvents.Branch("dR_g1l1",dR_g1l1,"dR_g1l1/F")
passedEvents.Branch("dR_g1l2",dR_g1l2,"dR_g1l2/F")
passedEvents.Branch("dR_g2l1",dR_g2l1,"dR_g2l1/F")
passedEvents.Branch("dR_g2l2",dR_g2l2,"dR_g2l2/F")

passedEvents.Branch("dR_pho",dR_pho,"dR_pho/F")

passedEvents.Branch("var_dR_Za",var_dR_Za,"var_dR_Za/F")
passedEvents.Branch("var_dR_g1g2",var_dR_g1g2,"var_dR_g1g2/F")
passedEvents.Branch("var_dR_g1Z",var_dR_g1Z,"var_dR_g1Z/F")
passedEvents.Branch("var_dEta_g1Z",var_dEta_g1Z,"var_dEta_g1Z/F")
passedEvents.Branch("var_dPhi_g1Z",var_dPhi_g1Z,"var_dPhi_g1Z/F")
passedEvents.Branch("var_PtaOverMa",var_PtaOverMa,"var_PtaOverMa/F")
passedEvents.Branch("var_PtaOverMh",var_PtaOverMh,"var_PtaOverMh/F")
passedEvents.Branch("var_Pta",var_Pta,"var_Pta/F")
passedEvents.Branch("var_MhMa",var_MhMa,"var_MhMa/F")
passedEvents.Branch("var_MhMZ",var_MhMZ,"var_MhMZ/F")

passedEvents.Branch("ALP_calculatedPhotonIso",ALP_calculatedPhotonIso,"ALP_calculatedPhotonIso/F")

passedEvents.Branch("passEleVeto",passEleVeto,"passEleVeto/F")
passedEvents.Branch("passIeIe",passIeIe,"passIeIe/F")
passedEvents.Branch("passHOverE",passHOverE,"passHOverE/F")
passedEvents.Branch("passChaHadIso",passChaHadIso,"passChaHadIso/F")
passedEvents.Branch("passNeuHadIso",passNeuHadIso,"passNeuHadIso/F")
passedEvents.Branch("passPhoIso",passPhoIso,"passPhoIso/F")
passedEvents.Branch("passdR_gl",passdR_gl,"passdR_gl/F")
passedEvents.Branch("passdR_gg",passdR_gg,"passdR_gg/F")
passedEvents.Branch("passH_m",passH_m,"passH_m/F")
passedEvents.Branch("passBDT",passBDT,"passBDT/F")
passedEvents.Branch("Val_BDT",Val_BDT,"Val_BDT/F")

passedEvents.Branch("pho1passCutBasedIDTight",pho1passCutBasedIDTight,"pho1passCutBasedIDTight/F")
passedEvents.Branch("pho2passCutBasedIDTight",pho2passCutBasedIDTight,"pho2passCutBasedIDTight/F")
passedEvents.Branch("pho1photonCutBasedIDMedium",pho1photonCutBasedIDMedium,"pho1photonCutBasedIDMedium/F")
passedEvents.Branch("pho2photonCutBasedIDMedium",pho2photonCutBasedIDMedium,"pho2photonCutBasedIDMedium/F")

passedEvents.Branch("pho1_matchedR03_PdgId",pho1_matchedR03_PdgId,"pho1_matchedR03_PdgId/F")
passedEvents.Branch("pho1_matchedR03_MomId",pho1_matchedR03_MomId,"pho1_matchedR03_MomId/F")
passedEvents.Branch("pho1_matchedR03_MomMomId",pho1_matchedR03_MomMomId,"pho1_matchedR03_MomMomId/F")
passedEvents.Branch("pho1_matchedR",pho1_matchedR,"pho1_matchedR/F")
passedEvents.Branch("pho1_matchedPt",pho1_matchedPt,"pho1_matchedPt/F")

passedEvents.Branch("pho2_matchedR03_PdgId",pho2_matchedR03_PdgId,"pho2_matchedR03_PdgId/F")
passedEvents.Branch("pho2_matchedR03_MomId",pho2_matchedR03_MomId,"pho2_matchedR03_MomId/F")
passedEvents.Branch("pho2_matchedR03_MomMomId",pho2_matchedR03_MomMomId,"pho2_matchedR03_MomMomId/F")
passedEvents.Branch("pho2_matchedR",pho2_matchedR,"pho2_matchedR/F")
passedEvents.Branch("pho2_matchedPt",pho2_matchedPt,"pho2_matchedPt/F")

mass_list = {'M1':1.0, 'M2':2.0, 'M3':3.0, 'M4':4.0, 'M5':5.0, 'M6':6.0, 'M7':7.0, 'M8':8.0, 'M9':9.0, 'M10':10.0, 'M15':15.0, 'M20':20.0, 'M25':25.0, 'M30':30.0}
mvaCut = {'M1':0.955, 'M2':0.98, 'M3':0.985, 'M4':0.98, 'M5':0.985, 'M6':0.99, 'M7':0.985, 'M8':0.99, 'M9':0.99, 'M10':0.99, 'M15':0.99, 'M20':0.99, 'M25':0.985, 'M30':0.98}
################################################################################################

#reader = bookMVA()
########################################################################################################################################
#Loop over all the events in the input ntuple
for ievent,event in enumerate(tchain):#, start=650000):
    if ievent > args.maxevents and args.maxevents != -1: break

    # if ievent == 500: break

    if ievent % 10000 == 0: print('Processing entry ' + str(ievent))

    #if ievent == 5076520: continue


    # Loop over all the electrons in an event

    # pho parameters
    foundpho1 = False
    foundpho2 = False
    pho1_maxPt = 0.0
    pho2_maxPt = 0.0
    pho1_index = 0
    pho2_index = 0
    pho_passEleVeto = True
    pho_passIeIe = True
    pho_passHOverE = True
    pho_passChaHadIso = True
    pho_passNeuHadIso = True
    passedPhoIso = True

    mva_value = -999.0

    # initialise Z parameters
    Nlep = 0
    lep=0

    Zmass = 91.1876
    dZmass = 9999.0
    n_Zs = 0
    Z_pt = []
    Z_eta = []
    Z_phi = []
    Z_mass = []
    Z_index = 0
    Z_lepindex1 = []
    Z_lepindex2 = []
    foundZ = False
    lep_leadindex = [] # lepindex[0] for leading, lepindex[1] for subleading

    # pass trigger
################################################################################################
    h_n.Fill(event.passedTrig)
    h_n_weight.Fill(event.passedTrig, event.pileupWeight*event.genWeight*weight)
    if (not event.passedTrig): continue # check trigger eff
    h_n_trig.Fill(event.passedTrig)

    # find all Z candidates
################################################################################################
    Nlep = event.lep_pt.size()
    nlep.Fill(Nlep)

    if Nlep!=2: continue

    for i in range(Nlep):

        for j in range(i+1,Nlep):

            if ((event.lep_id[i] + event.lep_id[j]) != 0): continue

            lifsr = ROOT.TLorentzVector()
            ljfsr = ROOT.TLorentzVector()
            lifsr.SetPtEtaPhiM(event.lepFSR_pt[i], event.lepFSR_eta[i], event.lepFSR_phi[i], event.lepFSR_mass[i])
            ljfsr.SetPtEtaPhiM(event.lepFSR_pt[j], event.lepFSR_eta[j], event.lepFSR_phi[j], event.lepFSR_mass[j])

            Z = ROOT.TLorentzVector()
            Z = (lifsr + ljfsr)

            if (Z.M()>0.0):
                n_Zs = n_Zs + 1
                Z_pt.append(Z.Pt())
                Z_eta.append(Z.Eta())
                Z_phi.append(Z.Phi())
                Z_mass.append(Z.M())
                Z_lepindex1.append(i)
                Z_lepindex2.append(j)

                foundZ = True
        # lep j
    # lep i

    if (not foundZ): continue

    # find Z
    for i in range(n_Zs):
        if (abs(Z_mass[i] - Zmass) <= dZmass):
            dZmass = abs(Z_mass[i] - Zmass)
            Z_index = i

    # find Z end

    if (event.lepFSR_pt[Z_lepindex1[Z_index]] >= event.lepFSR_pt[Z_lepindex2[Z_index]]):
        lep_leadindex.append(Z_lepindex1[Z_index])
        lep_leadindex.append(Z_lepindex2[Z_index])
    else:
        lep_leadindex.append(Z_lepindex2[Z_index])
        lep_leadindex.append(Z_lepindex1[Z_index])
################################################################################################



    l1_find = ROOT.TLorentzVector()
    l2_find = ROOT.TLorentzVector()
    Z_find = ROOT.TLorentzVector()
    
    l1_noFSR_find = ROOT.TLorentzVector()
    l2_noFSR_find = ROOT.TLorentzVector()
    Z_noFSR_find = ROOT.TLorentzVector()

    l1_find.SetPtEtaPhiM(event.lepFSR_pt[lep_leadindex[0]], event.lepFSR_eta[lep_leadindex[0]], event.lepFSR_phi[lep_leadindex[0]], event.lepFSR_mass[lep_leadindex[0]])
    l2_find.SetPtEtaPhiM(event.lepFSR_pt[lep_leadindex[1]], event.lepFSR_eta[lep_leadindex[1]], event.lepFSR_phi[lep_leadindex[1]], event.lepFSR_mass[lep_leadindex[1]])

    Z_find = (l1_find + l2_find)

    l1_noFSR_find.SetPtEtaPhiM(event.lep_pt[lep_leadindex[0]], event.lep_eta[lep_leadindex[0]], event.lep_phi[lep_leadindex[0]], event.lep_mass[lep_leadindex[0]])
    l2_noFSR_find.SetPtEtaPhiM(event.lep_pt[lep_leadindex[1]], event.lep_eta[lep_leadindex[1]], event.lep_phi[lep_leadindex[1]], event.lep_mass[lep_leadindex[1]])

    Z_noFSR_find = (l1_noFSR_find + l2_noFSR_find)

    if (abs(event.lep_id[lep_leadindex[0]]) == 11):
        Z_e_nocut.Fill(Z_find.M())
    if (abs(event.lep_id[lep_leadindex[0]]) == 13):
        Z_mu_nocut.Fill(Z_find.M())

        # Leptons Cuts
################################################################################################
    # pass lep isolation

    if (abs(event.lep_id[lep_leadindex[0]]) == 11):
        if (event.lep_RelIsoNoFSR[lep_leadindex[0]] > 0.35): continue
        if (event.lep_RelIsoNoFSR[lep_leadindex[1]] > 0.35): continue

        # pt Cut
        if (event.lepFSR_pt[lep_leadindex[0]] <= 25): continue
        if (event.lepFSR_pt[lep_leadindex[1]] <= 15): continue

    if (abs(event.lep_id[lep_leadindex[0]]) == 13):
        if (event.lep_RelIsoNoFSR[lep_leadindex[0]] > 0.35): continue
        if (event.lep_RelIsoNoFSR[lep_leadindex[1]] > 0.35): continue

        # pt Cut
        if (event.lepFSR_pt[lep_leadindex[0]] <= 20): continue
        if (event.lepFSR_pt[lep_leadindex[1]] <= 10): continue


    if (abs(event.lep_id[lep_leadindex[0]]) == 11):
        Z_e_lIso.Fill(Z_find.M())
    if (abs(event.lep_id[lep_leadindex[0]]) == 13):
        Z_mu_lIso.Fill(Z_find.M())

    # lep Tight ID Cut
    if (not (event.lep_tightId[lep_leadindex[0]])): continue
    if (not (event.lep_tightId[lep_leadindex[1]])): continue

    if (abs(event.lep_id[lep_leadindex[0]]) == 11):
        Z_e_lIso_lTight.Fill(Z_find.M())
    if (abs(event.lep_id[lep_leadindex[0]]) == 13):
        Z_mu_lIso_lTight.Fill(Z_find.M())
################################################################################################


    # m_Z > 50 GeV
################################################################################################
    if (Z_find.M() < 50): continue
    Z_50.Fill(Z_find.M())

################################################################################################


    # Find photon
############################################################

    N_pho = event.pho_pt.size()
    npho.Fill(N_pho)
    if (N_pho < 2): continue

    for i in range(N_pho):
        if (event.pho_hasPixelSeed[i] == 1): continue
        if (event.pho_pt[i]<10) : continue
        #dR_l1gi = deltaR(l1_find.Eta(), l1_find.Phi(), event.pho_eta[i], event.pho_phi[i])
        #dR_l2gi = deltaR(l2_find.Eta(), l2_find.Phi(), event.pho_eta[i], event.pho_phi[i])
        #if (dR_l1gi < 0.4) or (dR_l2gi < 0.4): continue
        if (event.pho_pt[i] > pho1_maxPt):
            pho1_maxPt = event.pho_pt[i]
            pho1_index = i
            foundpho1 = True

    for j in range(N_pho):
        if (event.pho_hasPixelSeed[j] == 1): continue
        if j == pho1_index: continue
        if (event.pho_pt[j]<10) : continue
        #dR_l1gj = deltaR(l1_find.Eta(), l1_find.Phi(), event.pho_eta[j], event.pho_phi[j])
        #dR_l2gj = deltaR(l2_find.Eta(), l2_find.Phi(), event.pho_eta[j], event.pho_phi[j])
        #if (dR_l1gj < 0.4) or (dR_l2gj < 0.4): continue
        if (event.pho_pt[j] > pho2_maxPt):
            pho2_maxPt = event.pho_pt[j]
            pho2_index = j
            foundpho2 = True

    if (foundpho1 and foundpho2):


    ################################################################################################

        pho1_find = ROOT.TLorentzVector()
        pho2_find = ROOT.TLorentzVector()

        #pho1_find.SetPtEtaPhiE(event.pho_pt[pho1_index], event.pho_eta[pho1_index], event.pho_phi[pho1_index], event.pho_pt[pho1_index] * np.cosh(event.pho_eta[pho1_index]))
        #pho2_find.SetPtEtaPhiE(event.pho_pt[pho2_index], event.pho_eta[pho2_index], event.pho_phi[pho2_index], event.pho_pt[pho2_index] * np.cosh(event.pho_eta[pho2_index]))

        pho1_find.SetPtEtaPhiM(event.pho_pt[pho1_index], event.pho_eta[pho1_index], event.pho_phi[pho1_index], 0.0)
        pho2_find.SetPtEtaPhiM(event.pho_pt[pho2_index], event.pho_eta[pho2_index], event.pho_phi[pho2_index], 0.0)

        ALP_find = ROOT.TLorentzVector()
        ALP_find = (pho1_find + pho2_find)
    #######################################################################################################

        # Higgs Candidates
    #######################################################################################################
        H_find = ROOT.TLorentzVector()
        H_find = (Z_find + ALP_find)

        H_noFSR_find = ROOT.TLorentzVector()
        H_noFSR_find = (Z_noFSR_find + ALP_find)
    #######################################################################################################

        # Photon Cuts
    #######################################################################################################

        if (((abs(event.pho_eta[pho1_index]) > 1.566) and (abs(event.pho_eta[pho1_index]) < 2.5)) or (abs(event.pho_eta[pho1_index]) < 1.4442)) and (((abs(event.pho_eta[pho2_index]) > 1.566) and (abs(event.pho_eta[pho2_index]) < 2.5)) or (abs(event.pho_eta[pho2_index]) < 1.4442)):
            #if (((abs(event.pho_eta[pho2_index]) >1.4442) and (abs(event.pho_eta[pho2_index]) < 1.566)) or (abs(event.pho_eta[pho2_index]) >2.5)): continue

            if (event.pho_EleVote[pho1_index] == 0 ): pho_passEleVeto = False
            if (event.pho_EleVote[pho2_index] == 0 ): pho_passEleVeto = False

            if deltaR(pho1_find.Eta(),pho1_find.Phi(),pho2_find.Eta(),pho2_find.Phi()) < 0.3:
                pho1_phoIso = event.pho_photonIso[pho1_index] - pho2_find.Pt()
                pho2_phoIso = event.pho_photonIso[pho2_index] - pho1_find.Pt()
            else:
                pho1_phoIso = event.pho_photonIso[pho1_index]
                pho2_phoIso = event.pho_photonIso[pho2_index]


            # photon 1
            # barrel
            if (abs(event.pho_eta[pho1_index]) < 1.4442):
                if (event.pho_full5x5_sigmaIetaIeta[pho1_index] > sigieieCut_EB): pho_passIeIe = False
                if (event.pho_hadronicOverEm[pho1_index] > hoeCut_EB): pho_passHOverE = False
                if (event.pho_chargedHadronIso[pho1_index] > chCut_EB ): pho_passChaHadIso = False
                if (event.pho_neutralHadronIso[pho1_index] > (neuCut_EB[0] + event.pho_pt[pho1_index]*neuCut_EB[1] + event.pho_pt[pho1_index]*event.pho_pt[pho1_index]*neuCut_EB[2])): pho_passNeuHadIso = False
                if (event.pho_photonIso[pho1_index] > (phoCut_EB[0] + event.pho_pt[pho1_index]*phoCut_EB[1])): passedPhoIso = False

            # endcap
            else:
                if (event.pho_full5x5_sigmaIetaIeta[pho1_index] > sigieieCut_EE): pho_passIeIe = False
                if (event.pho_hadronicOverEm[pho1_index] > hoeCut_EE): pho_passHOverE = False
                if (event.pho_chargedHadronIso[pho1_index] > chCut_EE ): pho_passChaHadIso = False
                if (event.pho_neutralHadronIso[pho1_index] > (neuCut_EE[0] + event.pho_pt[pho1_index]*neuCut_EE[1] + event.pho_pt[pho1_index]*event.pho_pt[pho1_index]*neuCut_EE[2])): pho_passNeuHadIso = False
                if (event.pho_photonIso[pho1_index] > (phoCut_EE[0] + event.pho_pt[pho1_index]*phoCut_EE[1])): passedPhoIso = False
            # photon 2
            # barrel
            if (abs(event.pho_eta[pho2_index]) < 1.4442):
                if (event.pho_full5x5_sigmaIetaIeta[pho2_index] > sigieieCut_EB): pho_passIeIe = False
                if (event.pho_hadronicOverEm[pho2_index] > hoeCut_EB): pho_passHOverE = False
                if (event.pho_chargedHadronIso[pho2_index] > chCut_EB ): pho_passChaHadIso = False
                if (event.pho_neutralHadronIso[pho2_index] > (neuCut_EB[0] + event.pho_pt[pho1_index]*neuCut_EB[1] + event.pho_pt[pho1_index]*event.pho_pt[pho1_index]*neuCut_EB[2])): pho_passNeuHadIso = False
                if (event.pho_photonIso[pho2_index] > (phoCut_EB[0] + event.pho_pt[pho1_index]*phoCut_EB[1])): passedPhoIso = False

            # endcap
            else:
                if (event.pho_full5x5_sigmaIetaIeta[pho2_index] > sigieieCut_EE): pho_passIeIe = False
                if (event.pho_hadronicOverEm[pho2_index] > hoeCut_EE): pho_passHOverE = False
                if (event.pho_chargedHadronIso[pho2_index] > chCut_EE ): pho_passChaHadIso = False
                if (event.pho_neutralHadronIso[pho2_index] > (neuCut_EE[0] + event.pho_pt[pho1_index]*neuCut_EE[1] + event.pho_pt[pho1_index]*event.pho_pt[pho1_index]*neuCut_EE[2])): pho_passNeuHadIso = False
                if (event.pho_photonIso[pho2_index] > (phoCut_EE[0] + event.pho_pt[pho1_index]*phoCut_EE[1])): passedPhoIso = False

            #H_Ceta.Fill(H_find.M())


            Run[0] = event.Run
            LumiSect[0] = event.LumiSect
            Event[0] = event.Event
            if isMC:
                lep_dataMC = event.lep_dataMC[lep_leadindex[0]] * event.lep_dataMC[lep_leadindex[1]]
            else:
                lep_dataMC = 1.0
            factor[0] = event.genWeight * event.pileupWeight * lep_dataMC * weight
            event_weight[0] = weight
            event_genWeight[0] = event.genWeight
            event_pileupWeight[0] = event.pileupWeight
            event_pileupWeightUp[0] = event.pileupWeightUp
            event_pileupWeightDn[0] = event.pileupWeightDn
            event_dataMCWeight[0] = lep_dataMC
            event_SelectEff_weight[0] = event.pileupWeight*event.genWeight*weight
            event_passedTrig[0] = event.passedTrig
            #event_qcdweight[0] = event.qcdWeights
            #event_nnloweight[0] = event.nnloWeights

            '''
            if isMC:
                ### pdf uncertainty
                pdf_weight_id = 48

                temp_PDFweight = 0.0
                for i in range(pdf_weight_id, pdf_weight_id+100):
                    temp_PDFweight = temp_PDFweight + event.pdfWeights[i] * event.pdfWeights[i]
                
                if temp_PDFweight/100.0 > 1.0:
                    event_pdfweight_up[0] = math.sqrt(temp_PDFweight/100.0)
                    event_pdfweight_dn[0] = 1.0/math.sqrt(temp_PDFweight/100.0)
                else:
                    event_pdfweight_up[0] = 1.0/math.sqrt(temp_PDFweight/100.0)
                    event_pdfweight_dn[0] = math.sqrt(temp_PDFweight/100.0)

                temp_qcdweight_up = 1.0
                temp_qcdweight_dn = 1.0

                for i in range(0,pdf_weight_id-5,5):
                    if event.pdfWeights[i] > temp_qcdweight_up:
                        temp_qcdweight_up = event.pdfWeights[i]
                    if event.pdfWeights[i] < temp_qcdweight_dn:
                        temp_qcdweight_dn = event.pdfWeights[i]

                event_qcdweight_up[0] = temp_qcdweight_up
                event_qcdweight_dn[0] = temp_qcdweight_dn


            else:
                event_pdfweight_up[0] = 1.0
                event_pdfweight_dn[0] = 1.0
                event_qcdweight_up[0] = 1.0
                event_qcdweight_dn[0] = 1.0
            '''

            Z_m[0] = Z_find.M()
            H_m[0] = H_find.M()
            ALP_m[0] = ALP_find.M()
            H_pt[0] = H_find.Pt()

            Z_noFSR_m[0] = Z_noFSR_find.M()
            H_noFSR_m[0] = H_noFSR_find.M()

            l1_id[0] = event.lep_id[lep_leadindex[0]]
            l2_id[0] = event.lep_id[lep_leadindex[1]]
            l1_pt[0] = event.lepFSR_pt[lep_leadindex[0]]
            l2_pt[0] = event.lepFSR_pt[lep_leadindex[1]]
            l1_eta[0] = event.lepFSR_eta[lep_leadindex[0]]
            l2_eta[0] = event.lepFSR_eta[lep_leadindex[1]]
            l1_phi[0] = event.lepFSR_phi[lep_leadindex[0]]
            l2_phi[0] = event.lepFSR_phi[lep_leadindex[1]]
            l1_mass[0] = event.lepFSR_mass[lep_leadindex[0]]
            l2_mass[0] = event.lepFSR_mass[lep_leadindex[1]]
            l1_dataMC[0] = event.lep_dataMC[lep_leadindex[0]]
            l2_dataMC[0] = event.lep_dataMC[lep_leadindex[1]]

            l1_noFSR_pt[0] = event.lep_pt[lep_leadindex[0]]
            l2_noFSR_pt[0] = event.lep_pt[lep_leadindex[1]]
            l1_noFSR_eta[0] = event.lep_eta[lep_leadindex[0]]
            l2_noFSR_eta[0] = event.lep_eta[lep_leadindex[1]]
            l1_noFSR_phi[0] = event.lep_phi[lep_leadindex[0]]
            l2_noFSR_phi[0] = event.lep_phi[lep_leadindex[1]]
            l1_noFSR_mass[0] = event.lep_mass[lep_leadindex[0]]
            l2_noFSR_mass[0] = event.lep_mass[lep_leadindex[1]]

            l1_dataMCErr[0] = event.lep_dataMCErr[lep_leadindex[0]]
            l1_scaleup[0] = event.lep_scaleup[lep_leadindex[0]]
            l1_scaledn[0] = event.lep_scaledn[lep_leadindex[0]]
            l1_smearup[0] = event.lep_smearup[lep_leadindex[0]]
            l1_smeardn[0] = event.lep_smeardn[lep_leadindex[0]]
            l1_EPostCorr[0] = event.lep_ecalTrkEnergyPostCorr[lep_leadindex[0]]
            l1_EPreCorr[0] = event.lep_ecalEnergyPreCorr[lep_leadindex[0]]
            l2_dataMCErr[0] = event.lep_dataMCErr[lep_leadindex[1]]
            l2_scaleup[0] = event.lep_scaleup[lep_leadindex[1]]
            l2_scaledn[0] = event.lep_scaledn[lep_leadindex[1]]
            l2_smearup[0] = event.lep_smearup[lep_leadindex[1]]
            l2_smeardn[0] = event.lep_smeardn[lep_leadindex[1]]
            l2_EPostCorr[0] = event.lep_ecalTrkEnergyPostCorr[lep_leadindex[1]]
            l2_EPreCorr[0] = event.lep_ecalEnergyPreCorr[lep_leadindex[1]]

            #if abs(event.lep_id[lep_leadindex[0]]) == 11:
            #    print "lep_smearup:",event.lep_smearup[lep_leadindex[0]], "l1_find.E():",l1_find.E(), "lep_ecalTrkEnergyPostCorr:",event.lep_ecalTrkEnergyPostCorr[lep_leadindex[0]], "lep_ecalEnergyPreCorr:",event.lep_ecalEnergyPreCorr[lep_leadindex[0]], "lep_ecalEnergy:",event.lep_ecalEnergy[lep_leadindex[0]], "lep_smearup/l1_find.E():",event.lep_smearup[lep_leadindex[0]]/l1_find.E()

            # print('len of lep_matchedR03_PdgId')
            # print(len(event.lep_matchedR03_PdgId))
            # print(f'lep_leadindex[0]= {lep_leadindex[0]}')

            if(len(event.lep_matchedR03_PdgId)!=0): 
                l1_matchedR03_PdgId[0] = event.lep_matchedR03_PdgId[lep_leadindex[0]]
                l1_matchedR03_MomId[0] = event.lep_matchedR03_MomId[lep_leadindex[0]]
                l1_matchedR03_MomMomId[0] = event.lep_matchedR03_MomMomId[lep_leadindex[0]]

                l2_matchedR03_PdgId[0] = event.lep_matchedR03_PdgId[lep_leadindex[1]]
                l2_matchedR03_MomId[0] = event.lep_matchedR03_MomId[lep_leadindex[1]]
                l2_matchedR03_MomMomId[0] = event.lep_matchedR03_MomMomId[lep_leadindex[1]]
            
            ################################## Photon variables ##############################################################

            dR_l1g1 = deltaR(l1_find.Eta(), l1_find.Phi(), pho1_find.Eta(), pho1_find.Phi())
            dR_l1g2 = deltaR(l1_find.Eta(), l1_find.Phi(), pho2_find.Eta(), pho2_find.Phi())
            dR_l2g1 = deltaR(l2_find.Eta(), l2_find.Phi(), pho1_find.Eta(), pho1_find.Phi())
            dR_l2g2 = deltaR(l2_find.Eta(), l2_find.Phi(), pho2_find.Eta(), pho2_find.Phi())

            dR_g1g2 = deltaR(pho1_find.Eta(), pho1_find.Phi(), pho2_find.Eta(), pho2_find.Phi())
            dR_g1Z = deltaR(pho1_find.Eta(), pho1_find.Phi(), Z_find.Eta(), Z_find.Phi())

            dR_g1l1[0] = dR_l1g1
            dR_g1l2[0] = dR_l1g2
            dR_g2l1[0] = dR_l2g1
            dR_g2l2[0] = dR_l2g2

            dR_pho[0] = dR_g1g2

            pho1Pt[0] = event.pho_pt[pho1_index]
            pho1eta[0] = event.pho_eta[pho1_index]
            pho1phi[0] = event.pho_phi[pho1_index]
            pho1R9[0] = event.pho_R9[pho1_index]
            pho1IetaIeta[0] = event.pho_sigmaIetaIeta[pho1_index]
            pho1IetaIeta55[0] = event.pho_full5x5_sigmaIetaIeta[pho1_index]
            pho1HOE[0] = event.pho_hadronicOverEm[pho1_index]
            pho1CIso[0] = event.pho_chargedHadronIso[pho1_index]
            pho1NIso[0] = event.pho_neutralHadronIso[pho1_index]
            pho1PIso[0] = pho1_phoIso
            pho1PIso_noCorr[0] = event.pho_photonIso[pho1_index]
            pho1scaleup[0] = event.pho_scaleup[pho1_index]
            pho1scaledn[0] = event.pho_scaledn[pho1_index]
            pho1smearup[0] = event.pho_smearup[pho1_index]
            pho1smeardn[0] = event.pho_smeardn[pho1_index]
            pho1EPostCorr[0] = event.pho_ecalEnergyPostCorr[pho1_index]
            pho1ShowerShapeSys[0] = Sys.showerShapeUncVal(event.pho_eta[pho1_index], event.pho_R9[pho1_index])[1]
            pho1trackIso[0] = event.pho_trackIso[pho1_index]

            pho2Pt[0] = event.pho_pt[pho2_index]
            pho2eta[0] = event.pho_eta[pho2_index]
            pho2phi[0] = event.pho_phi[pho2_index]
            pho2R9[0] = event.pho_R9[pho2_index]
            pho2IetaIeta[0] = event.pho_sigmaIetaIeta[pho2_index]
            pho2IetaIeta55[0] = event.pho_full5x5_sigmaIetaIeta[pho2_index]
            pho2HOE[0] = event.pho_hadronicOverEm[pho2_index]
            pho2CIso[0] = event.pho_chargedHadronIso[pho2_index]
            pho2NIso[0] = event.pho_neutralHadronIso[pho2_index]
            pho2PIso[0] = pho2_phoIso
            pho2PIso_noCorr[0] = event.pho_photonIso[pho2_index]
            pho2scaleup[0] = event.pho_scaleup[pho2_index]
            pho2scaledn[0] = event.pho_scaledn[pho2_index]
            pho2smearup[0] = event.pho_smearup[pho2_index]
            pho2smeardn[0] = event.pho_smeardn[pho2_index]
            pho2EPostCorr[0] = event.pho_ecalEnergyPostCorr[pho2_index]
            pho2ShowerShapeSys[0] = Sys.showerShapeUncVal(event.pho_eta[pho2_index], event.pho_R9[pho2_index])[1]
            pho2trackIso[0] = event.pho_trackIso[pho2_index]

            var_dR_Za[0] = deltaR(Z_find.Eta(), Z_find.Phi(), ALP_find.Eta(), ALP_find.Phi())
            var_dR_g1g2[0] = dR_g1g2
            var_dR_g1Z[0] = dR_g1Z
            var_dEta_g1Z[0] = pho1_find.Eta() - Z_find.Eta()
            var_dPhi_g1Z[0] = pho1_find.Phi() - Z_find.Phi()
            var_PtaOverMa[0] = ALP_find.Pt()/ALP_find.M()
            var_PtaOverMh[0] = ALP_find.Pt()/H_find.M()
            var_Pta[0] = ALP_find.Pt()
            var_MhMa[0] = ALP_find.M()+H_find.M()
            var_MhMZ[0] = Z_find.M()+H_find.M()

            ################# photon SFs

            pho1_SFs = 1.0
            pho1_SFs_sys = 0.0
            pho1_SFs_sys_replaced = 0.0
            pho2_SFs = 1.0
            pho2_SFs_sys = 0.0
            pho2_SFs_sys_replaced = 0.0
            if isMC:

                for i in range(len(bins_eta)-1):
                    for j in range(len(bins_pt)-1):
                        if pho1_find.Eta() > bins_eta[i] and pho1_find.Eta() < bins_eta[i+1]:
                            if pho1_find.Pt() > bins_pt[j] and pho1_find.Pt() < bins_pt[j+1]:
                                pho1_SFs = SFs[i][j]
                                pho1_SFs_sys = math.sqrt(SFs_statData[i][j]*SFs_statData[i][j]+SFs_statMC[i][j]*SFs_statMC[i][j]+SFs_altBkgModel[i][j]*SFs_altBkgModel[i][j]+SFs_altSignalModel[i][j]*SFs_altSignalModel[i][j]+SFs_altMCEff[i][j]*SFs_altMCEff[i][j]+SFs_altTagSelection[i][j]*SFs_altTagSelection[i][j])


                for i in range(len(bins_eta)-1):
                    for j in range(len(bins_pt)-1):
                        if pho2_find.Eta() > bins_eta[i] and pho2_find.Eta() < bins_eta[i+1]:
                            if pho2_find.Pt() > bins_pt[j] and pho2_find.Pt() < bins_pt[j+1]:
                                pho2_SFs = SFs[i][j]
                                pho2_SFs_sys = math.sqrt(SFs_statData[i][j]*SFs_statData[i][j]+SFs_statMC[i][j]*SFs_statMC[i][j]+SFs_altBkgModel[i][j]*SFs_altBkgModel[i][j]+SFs_altSignalModel[i][j]*SFs_altSignalModel[i][j]+SFs_altMCEff[i][j]*SFs_altMCEff[i][j]+SFs_altTagSelection[i][j]*SFs_altTagSelection[i][j])


            pho1SFs[0] = pho1_SFs
            pho1SFs_sys[0] = pho1_SFs_sys
            pho1SFs_sys_replaced[0] = pho1_SFs_sys
            pho2SFs[0] = pho2_SFs
            pho2SFs_sys[0] = pho2_SFs_sys
            pho2SFs_sys_replaced[0] = pho2_SFs_sys

            pt_mid = (event.pho_pt[pho1_index] + event.pho_pt[pho2_index])/2.0
            eta_mid = event.pho_eta[pho1_index]

            if pt_mid < bins_pt[-1]:
                for i in range(len(bins_eta)-1):
                        for j in range(len(bins_pt)-1):
                            if eta_mid > bins_eta[i] and eta_mid < bins_eta[i+1]:
                                if pt_mid > bins_pt[j] and pt_mid < bins_pt[j+1]:
                                    SFs_mid = SFs[i][j]
                                    pho_SFs_sys_mid = math.sqrt(SFs_statData[i][j]*SFs_statData[i][j]+SFs_statMC[i][j]*SFs_statMC[i][j]+SFs_altBkgModel[i][j]*SFs_altBkgModel[i][j]+SFs_altSignalModel[i][j]*SFs_altSignalModel[i][j]+SFs_altMCEff[i][j]*SFs_altMCEff[i][j]+SFs_altTagSelection[i][j]*SFs_altTagSelection[i][j])
            else:
                for i in range(len(bins_eta)-1):
                    if eta_mid > bins_eta[i] and eta_mid < bins_eta[i+1]:
                        SFs_mid = SFs[i][-1]
                        pho_SFs_sys_mid = math.sqrt(SFs_statData[i][-1]*SFs_statData[i][-1]+SFs_statMC[i][-1]*SFs_statMC[i][-1]+SFs_altBkgModel[i][-1]*SFs_altBkgModel[i][-1]+SFs_altSignalModel[i][-1]*SFs_altSignalModel[i][-1]+SFs_altMCEff[i][-1]*SFs_altMCEff[i][-1]+SFs_altTagSelection[i][-1]*SFs_altTagSelection[i][-1])
            #SFs_mid = his_SFs.GetBinContent(his_SFs.FindBin(pt_mid, eta_mid))

            if dR_g1g2 < 0.1:
                pho1SFs_dR0P1[0] = SFs_mid
                pho2SFs_dR0P1[0] = SFs_mid
            else:
                pho1SFs_dR0P1[0] = pho1_SFs
                pho2SFs_dR0P1[0] = pho2_SFs

            if dR_g1g2 < 0.15:
                pho1SFs_dR0P15[0] = SFs_mid
                pho2SFs_dR0P15[0] = SFs_mid
                pho1SFs_sys_replaced[0] = pho_SFs_sys_mid
                pho2SFs_sys_replaced[0] = pho_SFs_sys_mid
            else:
                pho1SFs_dR0P15[0] = pho1_SFs
                pho2SFs_dR0P15[0] = pho2_SFs

            if dR_g1g2 < 0.2:
                pho1SFs_dR0P2[0] = SFs_mid
                pho2SFs_dR0P2[0] = SFs_mid
            else:
                pho1SFs_dR0P2[0] = pho1_SFs
                pho2SFs_dR0P2[0] = pho2_SFs

            if dR_g1g2 < 0.25:
                pho1SFs_dR0P25[0] = SFs_mid
                pho2SFs_dR0P25[0] = SFs_mid
            else:
                pho1SFs_dR0P25[0] = pho1_SFs
                pho2SFs_dR0P25[0] = pho2_SFs

            if dR_g1g2 < 0.3:
                pho1SFs_dR0P3[0] = SFs_mid
                pho2SFs_dR0P3[0] = SFs_mid
            else:
                pho1SFs_dR0P3[0] = pho1_SFs
                pho2SFs_dR0P3[0] = pho2_SFs


            ################# calculate ALP's photon isolation #################
            checkDuplicate = ''
            ALP_calculatedPhotonIso_tmp = 0.0
            isReco = False

            for i in range(event.pho_PF_Id.size()):
                if str(event.pho_PF_Id[i]) in checkDuplicate: continue
                checkDuplicate = checkDuplicate + ' ' + str(event.pho_PF_Id[i])
                dR_AP = deltaR(ALP_find.Eta(), ALP_find.Phi(), event.pho_PF_eta[i], event.pho_PF_phi[i])

                # check if pfPhoton == recoPhoton
                for pho in range(N_pho):
                    dR_Pg = deltaR(event.pho_eta[pho], event.pho_phi[pho], event.pho_PF_eta[i], event.pho_PF_phi[i])
                    if dR_Pg < 0.08: isReco = True

                if isReco: continue

                if dR_AP <= 0.3:
                    ALP_calculatedPhotonIso_tmp += event.pho_PF_pt[i]

            #print checkDuplicate
            ALP_calculatedPhotonIso[0] = ALP_calculatedPhotonIso_tmp
            ################# END calculate ALP's photon isolation #################



            ################# define cut flow #################
            cutFlow = {}

            cutdR_gl = (dR_l1g1 > 0.4) and (dR_l1g2 > 0.4) and (dR_l2g1 > 0.4) and (dR_l2g2 > 0.4)

            cutdR_gg = dR_g1g2 > 0.02
            cutH_m = (H_find.M()>95) and (H_find.M()<180)

            cutFlow['cut1'] = pho_passChaHadIso
            cutFlow['cut2'] = pho_passNeuHadIso
            cutFlow['cut3'] = cutdR_gl
            cutFlow['cut4'] = cutdR_gg
            cutFlow['cut5'] = pho_passHOverE
            cutFlow['cut6'] = pho_passIeIe
            cutFlow['cut7'] = passedPhoIso
            cutFlow['cut8'] = cutH_m

            passEleVeto[0] = pho_passEleVeto
            passIeIe[0] = pho_passIeIe
            passHOverE[0] = pho_passHOverE
            passChaHadIso[0] = pho_passChaHadIso
            passNeuHadIso[0] = pho_passNeuHadIso
            passPhoIso[0] = passedPhoIso
            passdR_gl[0] = cutdR_gl
            passdR_gg[0] = cutdR_gg
            passH_m[0] = cutH_m

            pho1passCutBasedIDTight[0] = event.photonCutBasedIDTight[pho1_index]
            pho2passCutBasedIDTight[0] = event.photonCutBasedIDTight[pho2_index]
            pho1photonCutBasedIDMedium[0] = event.photonCutBasedIDMedium[pho1_index]
            pho2photonCutBasedIDMedium[0] = event.photonCutBasedIDMedium[pho2_index]

            if isMC:
                pho1_matchedR03_PdgId[0] = event.pho_matchedR03_PdgId[pho1_index]
                pho1_matchedR03_MomId[0] = event.pho_matchedR03_MomId[pho1_index]
                pho1_matchedR03_MomMomId[0] = event.pho_matchedR03_MomMomId[pho1_index]
                pho1_matchedR[0] = event.pho_matchedR[pho1_index]
                pho1_matchedPt[0] = event.pho_matchedPt[pho1_index]

                pho2_matchedR03_PdgId[0] = event.pho_matchedR03_PdgId[pho2_index]
                pho2_matchedR03_MomId[0] = event.pho_matchedR03_MomId[pho2_index]
                pho2_matchedR03_MomMomId[0] = event.pho_matchedR03_MomMomId[pho2_index]
                pho2_matchedR[0] = event.pho_matchedR[pho2_index]
                pho2_matchedPt[0] = event.pho_matchedPt[pho2_index]

            passBDT[0] = -99.
            Val_BDT[0] = -99.
            
            if pho_passChaHadIso and pho_passNeuHadIso and cutdR_gl and pho_passHOverE:
                # MVA_list = [event.pho_full5x5_sigmaIetaIeta[pho1_index], event.pho_photonIso[pho1_index], event.pho_full5x5_sigmaIetaIeta[pho2_index], event.pho_photonIso[pho2_index], ALP_calculatedPhotonIso_tmp, ALP_find.Pt()/ALP_find.M(), ALP_find.Pt()/H_find.M(), Z_find.M()+H_find.M()]
                param_mass = ''
                param_mass_value = 0.0
                for mass in mass_list.keys():
                    if mass == args.inputfiles[0].split('/')[-1].split('.')[0].split('_')[1]:
                        param_mass = mass
                        param_mass_value = mass_list[mass]
                        param = (ALP_find.M() - mass_list[mass])/H_find.M()

                        break
                    else:
                        continue
                    
                #print param_mass_value
                if param_mass_value:
                    MVA_list = [event.pho_pt[pho1_index], event.pho_R9[pho1_index], event.pho_full5x5_sigmaIetaIeta[pho1_index], event.pho_photonIso[pho1_index], event.pho_pt[pho2_index], event.pho_R9[pho2_index], event.pho_full5x5_sigmaIetaIeta[pho2_index], event.pho_photonIso[pho2_index], ALP_calculatedPhotonIso_tmp, var_dR_Za[0], dR_g1g2, dR_g1Z, var_PtaOverMh[0], H_find.Pt(), param]
                    MVA_value = model.predict_proba([MVA_list])[:, 1]
                    passBDT[0] = MVA_value>mvaCut[param_mass]
                    Val_BDT[0] = MVA_value[0]
                    # print (MVA_value[0])
            
            #######################################################################################################

            passedEvents.Fill()











file_out.Write()
file_out.Close()

sw.Stop()
print('Real time: ' + str(round(sw.RealTime() / 60.0,2)) + ' minutes')
print('CPU time:  ' + str(round(sw.CpuTime() / 60.0,2)) + ' minutes')
