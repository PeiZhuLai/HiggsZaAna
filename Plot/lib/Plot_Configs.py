import os
import sys
import Analyzer_Configs as AC
from ROOT import *

class Plot_Config:

    def __init__(self, ana_cfg, year):
        self.ana_cfg = ana_cfg
        self.colors  = {}
        self.var_title_map  = {}
        self.logY    = False
        self.year    = year
        if year == '2016':
            self.lumi    = '16.81'
        elif year == '-2016':
            self.lumi    = '19.52'
        elif year == '2017':
            self.lumi    = '41.48'
        elif year == '2018':
            self.lumi    = '59.83'
        elif year == 'run2Rereco':
            self.lumi    = '137.24'
        elif year == 'run2':
            self.lumi    = '138'
        elif year == '2022preEE':
            self.lumi    = '7.98'
        elif year == '2022postEE':
            self.lumi    = '26.67'
        elif year == '2023preBPix':
            self.lumi    = '17.79'
        elif year == '2023postBPix':
            self.lumi    = '9.45'
        elif year == 'run3':
            self.lumi    = '61.89'
        else:
            print('do not in 2016/2017/2018!')
            exit(0)
        self.sig_scale  = 5
        self.sig_weight = 0.4133

        self.LoadColors()

        self.LoadVariableNames()



    def LoadColors(self):
        self.colors["Data"] = kBlack
        self.colors["M1"] = TColor.GetColor("#3DFAFF")
        self.colors["M2"] = kYellow
        self.colors["M3"] = kGreen
        self.colors["M4"] = kOrange
        self.colors["M5"] = kPink
        self.colors["M6"] = kViolet
        self.colors["M7"] = kSpring
        self.colors["M8"] = kRed -7
        self.colors["M9"] = kYellow -7
        self.colors["M10"] = TColor.GetColor("#2C363F")
        self.colors["M15"] = kOrange -7
        self.colors["M20"] = TColor.GetColor("#E75A7C")
        self.colors["M25"] = kViolet -7
        self.colors["M30"] = TColor.GetColor("#FF934F")

        # self.colors["DYJetsToLL"]  =  kAzure + 7
        self.colors["DYJetsToLL"]  =  TColor.GetColor("#48BEFF") # blue
        self.colors["DYGto2LG"]    =  TColor.GetColor("#ffa90e") # orange


    def LoadVariableNames(self):
        self.var_title_map = {
            'Z_m': "m_{ll}",
            'H_m': "m_{ll#gamma#gamma}",
            'ALP_m': "m_{#gamma#gamma}",
            'pho1Pt': "p_{T,#gamma1}",
            'pho1eta': "#eta_{#gamma1}",
            'pho1phi': "#phi_{#gamma1}",
            'pho1R9': "R_{9,#gamma1}",
            'pho1IetaIeta55': "#sigma_{i#etai#eta,#gamma1}",
            'pho1ECALIso': "I_{ECAL,#gamma1}",
            'pho1CIso': "I_{ch,#gamma1}",
            'pho1HCALIso': "I_{HCAL,#gamma1}",
            'pho1HOE': "#gamma1 H/E",

            'pho2Pt': "p_{T,#gamma2}",
            'pho2eta': "#eta_{#gamma2}",
            'pho2phi': "#phi_{#gamma2}",
            'pho2R9': "R_{9,#gamma2}",
            'pho2IetaIeta55': "#sigma_{i#etai#eta,#gamma2}",
            'pho2ECALIso': "I_{ECAL, #gamma2}",
            'pho2CIso': "I_{ch, #gamma2}",
            'pho2HCALIso': "I_{HCAL, #gamma2}",
            'pho2HOE': "#gamma2 H/E",

            'ALP_calculatedPhotonIso': "I_{#gamma,ALPs}",
            'var_dR_Za': "#Delta R(Z,a)",
            'var_dR_g1g2': "#Delta R(#gamma1,#gamma2)",
            'var_dR_g1Z': "#Delta R(#gamma1, Z)",
            'var_PtaOverMh': "p_{T,a}/m_{ll#gamma#gamma}",
            'var_Pta': "p_{T,a}",
            'var_MhMZ': "m_{ll#gamma#gamma}+m_{ll}",
            'H_pt': "p_{T,H}",
            'var_PtaOverMa': "p_{T,a}/m_{#gamma#gamma}",
            'var_MhMa': "m_{ll#gamma#gamma}+m_{#gamma#gamma}",
            'param': "(m_{a}-m_{a,hyp})/m_{ll#gamma#gamma}"
        }

    def SetHistStyles(self, hist, sample):
        if sample == 'data' or sample == 'Data':
            hist.SetMarkerStyle(20)
        elif sample in self.ana_cfg.sig_names:
            hist.SetLineColor(self.colors[sample])
            hist.SetLineWidth(2)
            hist.SetFillColor(kGray)
        elif sample in self.ana_cfg.bkg_names:
            hist.SetFillColor(self.colors[sample])
            hist.SetLineColor(self.colors[sample])


    def SetHzaStyle():
        hzaStyle = TStyle("hzaPaperStyle","Hza Paper Style")

        hzaStyle.SetFrameFillColor(0)
        hzaStyle.SetStatColor(0)
        hzaStyle.SetOptStat(0)
        hzaStyle.SetTitleFillColor(0)
        hzaStyle.SetCanvasBorderMode(0)
        hzaStyle.SetPadBorderMode(0)
        hzaStyle.SetFrameBorderMode(0)
        hzaStyle.SetPadColor(kWhite)
        hzaStyle.SetCanvasColor(kWhite)

        hzaStyle.SetCanvasDefH(600) #Height of canvas
        hzaStyle.SetCanvasDefW(600) #Width of canvas
        hzaStyle.SetCanvasDefX(0)   #POsition on screen
        hzaStyle.SetCanvasDefY(0)

        hzaStyle.SetPadLeftMargin(0.13)
        hzaStyle.SetPadRightMargin(0.1)
        hzaStyle.SetPadTopMargin(0.085)
        hzaStyle.SetPadBottomMargin(0.12)

        # For hgg axis titles:
        hzaStyle.SetTitleColor(1, "XYZ")
        hzaStyle.SetTitleFont(42, "XYZ")
        hzaStyle.SetTitleSize(0.05, "XYZ")
        hzaStyle.SetTitleXOffset(0.95)#//0.9)
        hzaStyle.SetTitleYOffset(1.15)# // => 1.15 if exponents
        hzaStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
        hzaStyle.SetPadTickY(1)

        # For hgg axis labels:
        hzaStyle.SetLabelColor(1, "XYZ")
        hzaStyle.SetLabelFont(42, "XYZ")
        hzaStyle.SetLabelOffset(0.007, "XYZ")
        hzaStyle.SetLabelSize(0.04, "XYZ")

        # Legends
        hzaStyle.SetLegendBorderSize(0)
        hzaStyle.SetLegendFillColor(kWhite)
        hzaStyle.SetLegendFont(42)

        hzaStyle.SetFillColor(10)
        # Nothing for now
        hzaStyle.SetTextFont(42)
        hzaStyle.SetTextSize(0.03)
        hzaStyle.cd()
