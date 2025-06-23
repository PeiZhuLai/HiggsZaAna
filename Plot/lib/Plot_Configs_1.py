import os
import sys
import Analyzer_Configs as AC
from ROOT import *

class Plot_Config:

    def __init__(self, ana_cfg, year):
		self.ana_cfg = ana_cfg
		self.colors  = {}
		self.logY    = False
	self.year    = year
	if year == '2016':
		self.lumi    = '31.9'
	elif year == '2017':
		self.lumi    = '41.5'
	elif year == '2018':
    		self.lumi    = '56.96'
	elif year == 'run2':
		self.lumi    = '134.3'
	else:
		print 'do not in 2016/2017/2018!'
		exit(0)
    	self.sig_scale  = 5
	self.sig_weight = 0.4133

    	self.LoadColors()

    def LoadColors(self):
	self.colors["Data"] = kBlack

        #self.colors["ggH"] = kRed
        self.colors["M1"] = kRed
		self.colors["M2"] = kYellow
		self.colors["M3"] = kGreen
		self.colors["M4"] = kOrange
		self.colors["M5"] = kPink
		self.colors["M6"] = kViolet
		self.colors["M7"] = kSpring
		self.colors["M8"] = kRed -7
		self.colors["M9"] = kYellow -7
		self.colors["M10"] = kGreen -7
		self.colors["M15"] = kOrange -7
		self.colors["M20"] = kPink -7
        self.colors["M25"] = kViolet -7
		self.colors["M30"] = kSpring -7

        self.colors["VBF"] = kOrange + 7
        self.colors["ZH"] =  kBlue + 1
        self.colors["WH"] =  kGreen + 2
        self.colors["ttH"]  = kPink + 6

	# trying to use blue for Z, green for W, yellow for top, red for g/q

	self.colors["GGHToZG"] = kRed-7
	self.colors["VBFHToZG"] = kRed-7
	self.colors["QCD"] = kRed
	self.colors["DY1JetsToLL"]  =  kAzure + 7
        self.colors["DY2JetsToLL"]  =  kAzure + 7
        self.colors["DY3JetsToLL"]  =  kAzure + 7
        self.colors["DY4JetsToLL"]  =  kAzure + 7

        self.colors["DYJetsToLL"]  =  kAzure + 7

    	self.colors["TTTo2L2Nu"] = kYellow - 9
    	self.colors["ZGGToLLGG"]  = kViolet +6
        self.colors["TGJets"] =  kOrange + 6
        self.colors["VHToGG"] =  kRed - 7

	#self.colors["WGGJets"]  = kOrange + 6
    	#self.colors["WWToLNuQQ"] =  kSpring +6
        #self.colors["WWW"] =  kGreen - 9
    	self.colors["GJets"] =  kCyan - 7
    	self.colors["TTGG"] = kViolet - 9
        self.colors["TTJets"] = kOrange - 9
        self.colors["ZGToLLG"] = kViolet - 9
        self.colors["TTGJets"] = kPink + 6
        self.colors["ZZTo2L2Q"] = kOrange - 5
        self.colors["ZZTo4L"] = kSpring -1

    # def LoadColors(self):


    def SetHistStyles(self, hist, sample):
    	if sample == 'data':
    	    hist.SetMarkerStyle(20)
    	elif sample in self.ana_cfg.sig_names:
    	    hist.SetLineColor(self.colors[sample])
    	    hist.SetLineWidth(2)
    	    hist.SetFillColor(kGray)
    	elif sample in self.ana_cfg.bkg_names:
    	    hist.SetFillColor(self.colors[sample])
