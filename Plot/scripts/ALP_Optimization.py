####################################################
####################################################

import os
import sys
import numpy as np
from array import array
import math

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir))
PLOT_OUTPUT_DIR = os.path.join(PLOT_DIR, 'plots')

sys.path.insert(0, os.path.join(PLOT_DIR, 'lib'))
from ROOT import *
import Analyzer_Configs as AC
#import Plot_Configs     as PC

#from Analyzer_ALP import PIso2D, plot2D, plot2D_CONT

import CMS_lumi, tdrstyle

from xgboost import XGBClassifier
import pickle

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

#import argparse
from optparse import OptionParser
import json

#parser = argparse.ArgumentParser(description="A simple ttree plotter")
#parser.add_argument("-y", "--Year", dest="year", default="2017", help="which year's datasetes")
#parser.add_option("-c", "--nCats",   dest="nCats",    default=5,    type="int",    help="nCats"  )
#parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')

# Silent
ROOT.gErrorIgnoreLevel = 2000  # 等同於 ROOT.kError


parser = OptionParser()
parser.add_option("-y", "--Year", dest="year", default='run3', type="str", help="which year's datasetes")
parser.add_option('-o', "--outDir", dest='outDir', default="./optimize", type="string", help="outDir")
parser.add_option("--region", dest="region", type=int, default=0, help="0 for full region, 1 for signal region, 2 for sideband region")
parser.add_option('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_option('-p', '--plot', dest='plot', action='store_true', default=False, help='Plot?')
parser.add_option('-s', '--sigma', dest='sigma', action='store_true', default=False, help='mass region?')
parser.add_option('--sigVSscore', dest='sigVSscore', action='store_true', default=True, help='mass region?')
parser.add_option('--doOpt', dest='doOpt', action='store_true', default=False, help='whether optimize the category?')
parser.add_option('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_option('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
parser.add_option("-c", "--nCats",   dest="nCats",    default=5,    type="int",    help="nCats"  )

(options, args) = parser.parse_args()



gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


if options.year == '2016':
    file_out = 'plots_16'
    name_SR = "ALP_plot_data16_param_SR.root"
    name_CR = "ALP_plot_data16_param_CR.root"
elif options.year == '2017':
    file_out = 'plots_17'
    name_SR = "ALP_plot_data17_param_SR.root"
    name_CR = "ALP_plot_data17_param_CR.root"
elif options.year == '2018':
    file_out = 'plots_18'
    name_SR = "ALP_plot_data18_param_SR.root"
    name_CR = "ALP_plot_data18_param_CR.root"
elif options.year == 'run2':
    file_out = 'plots_run2UL'
    if options.ele:
        name_SR = "ALP_plot_run2_UL_SR_ele.root"
        name_CR = "ALP_plot_run2_UL_CR_ele.root"
    elif options.mu:
        name_SR = "ALP_plot_run2_UL_SR_mu.root"
        name_CR = "ALP_plot_run2_UL_CR_mu.root"
    else:
        name_SR = "ALP_plot_run2_UL_SR.root"
        name_CR = "ALP_plot_run2_UL_CR.root"

        #name_SR = "ALP_plot_run2_UL_onlyM10_SR.root"
        #name_CR = "ALP_plot_run2_UL_onlyM10_CR.root"
        
        #name_SR = "ALP_plot_run2_UL_Ns_10_SR.root"
        #name_CR = "ALP_plot_run2_UL_Ns_10_CR.root"
elif options.year == 'run3':
    file_out = os.path.join(PLOT_OUTPUT_DIR, 'plots_run3UL')
    if options.ele:
        name_SR = "ALP_plot_run3_UL_SR_ele.root"
        name_CR = "ALP_plot_run3_UL_CR_ele.root"
    elif options.mu:
        name_SR = "ALP_plot_run3_UL_SR_mu.root"
        name_CR = "ALP_plot_run3_UL_CR_mu.root"
    else:
        name_SR = "ALP_plot_run3_UL_SR.root"
        name_CR = "ALP_plot_run3_UL_CR.root"

        #name_SR = "ALP_plot_run3_UL_onlyM10_SR.root"
        #name_CR = "ALP_plot_run3_UL_onlyM10_CR.root"
        
        #name_SR = "ALP_plot_run3_UL_Ns_10_SR.root"
        #name_CR = "ALP_plot_run3_UL_Ns_10_CR.root"  
else:
    print ("do not include at 2016/2017/2018")
    exit(0)


def SetgStyle():

    gStyle.SetFrameFillColor(0)
    gStyle.SetStatColor(0)
    gStyle.SetOptStat(0)
    gStyle.SetTitleFillColor(0)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetPadColor(kWhite)
    gStyle.SetCanvasColor(kWhite)
    
    
    gStyle.SetCanvasDefH(600) #Height of canvas
    gStyle.SetCanvasDefW(600) #Width of canvas
    gStyle.SetCanvasDefX(0)   #POsition on screen
    gStyle.SetCanvasDefY(0)

    
    gStyle.SetPadLeftMargin(0.13)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadTopMargin(0.085)
    gStyle.SetPadBottomMargin(0.12)
    
    # For hgg axis titles:
    gStyle.SetTitleColor(1, "XYZ")
    gStyle.SetTitleFont(42, "XYZ")
    gStyle.SetTitleSize(0.04, "XYZ")
    gStyle.SetTitleXOffset(0.95)#//0.9)
    gStyle.SetTitleYOffset(1.15)# // => 1.15 if exponents
    
    # For hgg axis labels:
    gStyle.SetLabelColor(1, "XYZ")
    gStyle.SetLabelFont(42, "XYZ")
    gStyle.SetLabelOffset(0.007, "XYZ")
    gStyle.SetLabelSize(0.04, "XYZ")
    
    # Legends
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendFillColor(kWhite)
    gStyle.SetLegendFont(42)
    
    gStyle.SetFillColor(10)
    # Nothing for now
    gStyle.SetTextFont(42)
    gStyle.SetTextSize(0.03)


def open_file(file_name):
    in_file = TFile( file_name , "r")
    return in_file

def get_hist(file, mvaVal_name):
    hist = file.raw_plots.Get(mvaVal_name)
    return hist

def hist2graph(hist, mva_low = 0.1):
    
    Nbins = hist.GetNbinsX()
    bin_x_low = hist.FindBin(mva_low)
    xaxis = hist.GetXaxis()

    bin_x_Center=[]
    bin_y=[]
    for x in range(Nbins+1):
        # remove BDT score less than mva_low
        if x < bin_x_low: continue
        bin_x_Center.append(xaxis.GetBinCenter(x))
        bin_y.append(hist.GetBinContent(x))

    graph = TGraph(Nbins-bin_x_low+1, np.array(bin_x_Center), np.array(bin_y))

    return graph

def smooth(graph, hist, mva_low = 0.1):
    h_smooth = TH1F(hist.GetName()+"_smooth",hist.GetName()+"_smooth",hist.GetNbinsX(), -0.1, 1.1)
    h_smooth_up = TH1F(hist.GetName()+"_smooth_up",hist.GetName()+"_smooth_up",hist.GetNbinsX(), -0.1, 1.1)
    h_smooth_dn = TH1F(hist.GetName()+"_smooth_dn",hist.GetName()+"_smooth_dn",hist.GetNbinsX(), -0.1, 1.1)

    smoother = TGraphSmooth()
    g_smooth = smoother.SmoothSuper(graph)

    x = array('d', [0])
    y = array('d', [0])

    for i in range(1, hist.GetNbinsX()+1):

        if i < hist.FindBin(mva_low):
            y[0] = hist.GetBinContent(i)
        else:
            g_smooth.GetPoint(i - hist.FindBin(mva_low), x, y)

        if y[0] < 0:
            y[0] = 0.0

    
        if y[0] >= 0. and x[0] <=1.0 and x[0] >= 0.0:
            h_smooth.SetBinContent(i, y[0])

            h_smooth_up.SetBinContent(i, y + np.sqrt(y))
            if (y - np.sqrt(y)) > 0.: 
                h_smooth_dn.SetBinContent(i, y - np.sqrt(y))

            else: 
                h_smooth_dn.SetBinContent(i, 0.)

        else:
            h_smooth_up.SetBinContent(i, 0.)
            h_smooth_dn.SetBinContent(i, 0.)

 
    return [h_smooth, h_smooth_up, h_smooth_dn]


def compare(hist, hist_smooth, ma):
    if gROOT.FindObject("cc"):
        gROOT.FindObject("cc").Close()
    canv = TCanvas("cc", "cc", 800, 600)
    canv.cd()
    canv.SetLogy()
    SetgStyle()
    leftMargin=0.145
    rightMargin = 0.04
    topMargin = 0.08
    canv.SetMargin(leftMargin, rightMargin, 0.14, topMargin) # //left//right//bottom//top

    # canv.SetRightMargin(0.04)
    # canv.SetLeftMargin(0.145)
    # canv.SetTopMargin(0.05)
    # canv.SetBottomMargin(0.14)
    
    canv.SetTickx(1)
    canv.SetTicky(1)

    hist.SetMinimum(0.5)
    hist.SetMaximum( 10. * hist.GetMaximum() )

    hist.SetFillColor(0)
    hist.SetLineColor(1)
    hist.SetLineWidth(2)
    hist.GetXaxis().SetTitle('BDT Score')
    # hist.GetYaxis().SetTitle('Events')
    hist.GetYaxis().SetTitle(f'Events / {hist.GetBinWidth(1):.3f}')
    hist.GetXaxis().CenterTitle(True)
    hist.GetYaxis().CenterTitle(True)

    hist.GetXaxis().SetTitleColor(1)
    hist.GetXaxis().SetTitleFont(42)
    hist.GetXaxis().SetTitleSize(0.055)
    hist.GetXaxis().SetTitleOffset(1.15)

    hist.GetXaxis().SetLabelSize(0.05)
    hist.GetXaxis().SetLabelColor(1)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetXaxis().SetLabelOffset(0.007)

    hist.GetYaxis().SetTitleColor(1)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetYaxis().SetTitleSize(0.058)
    hist.GetYaxis().SetTitleOffset(1.2)

    hist.GetYaxis().SetLabelSize(0.05)
    hist.GetYaxis().SetLabelColor(1)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelOffset(0.007)

    hist.Draw("HIST")

    hist_smooth[1].SetLineColor(TColor.GetColor("#03CEA4"))
    hist_smooth[1].SetLineWidth(3)
    hist_smooth[1].SetLineStyle(1)
    hist_smooth[1].SetMarkerSize(0)
    hist_smooth[1].Draw("SAME")

    hist_smooth[0].SetLineColor(TColor.GetColor("#E40066"))
    hist_smooth[0].SetLineWidth(3)
    hist_smooth[0].SetLineStyle(1)
    hist_smooth[0].SetMarkerSize(0)
    hist_smooth[0].Draw("SAME")

    hist_smooth[2].SetLineColor(TColor.GetColor("#2F52E0"))
    hist_smooth[2].SetLineWidth(3)
    hist_smooth[2].SetLineStyle(1)
    hist_smooth[2].SetMarkerSize(0)
    hist_smooth[2].Draw("SAME")

    hist.GetXaxis().SetRangeUser(-0.1, 1.1)
    hist_smooth[0].GetXaxis().SetRangeUser(-0.1, 1.1)    
    hist_smooth[1].GetXaxis().SetRangeUser(-0.1, 1.1)    
    hist_smooth[2].GetXaxis().SetRangeUser(-0.1, 1.1)    
    
    #canv.cd()
    global legend
    legend = TLegend(0.65,0.62,0.94,0.88)
    #legend.AddEntry(hist, "Background")
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetFillColor(0)
    legend.SetTextSize(0.045)
    legend.AddEntry("", f"m_{{a}} = {ma.lstrip('M')} GeV", "")
    legend.AddEntry(hist, "Nominal", "l")
    legend.AddEntry(hist_smooth[1], "Smoothed + 1#sigma", "l")
    legend.AddEntry(hist_smooth[0], "Smoothed", "l")
    legend.AddEntry(hist_smooth[2], "Smoothed - 1#sigma", "l")
    legend.Draw("SAME")

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
 
    latex.SetTextAlign(13)
    latex.SetTextSize(0.045)  # 字體大小
    latex.SetTextFont(61)  # 粗體字 CMS 標籤
    latex.DrawLatex(leftMargin, 0.969, "CMS")

    latex.SetTextFont(52)  # 斜體字 Preliminary 標籤
    latex.DrawLatex(leftMargin + 0.09, 0.965, "Preliminary")

    latex.SetTextAlign(31)
    # latex.DrawLatex(0.73, 0.87, f"m_{{a}} = {ma.lstrip('M')} GeV")
    latex.DrawLatex(1. - rightMargin, 1. - topMargin + 0.01, ("170.84 fb^{-1} (13.6 TeV)"))
    
    # latex.DrawLatex(0.85, 0.88, r"m_{a} = %s GeV" % (ma.lstrip("M")))


    canv.Update()
    #print legend

    #cmsText     = "CMS";
    #cmsTextFont   = 61
    #latex = rt.TLatex()

    return canv


def GetHist(file, var_list):
    if isinstance(var_list, str):
        var_list = [var_list]
    
    hist_combined = None
    for i, var in enumerate(var_list):
        hist = get_hist(file, var)
        if i == 0:
            hist_combined = hist.Clone(f"combined_{var}")
        else:
            hist_combined.Add(hist)
    
    graph = hist2graph(hist_combined, 0.01)
    hist_smooth = smooth(graph, hist_combined, 0.01)

    return hist_combined, hist_smooth


def computeSignificance(s,b,d_noSmooth):
  significance = -999.
  if b>0 and s>0.: significance = (2*(s+b)*math.log(1+(s/b))) - 2*s
  if significance>0. and d_noSmooth>=0. and b>0.: return np.sqrt( (2*(s+b)*math.log(1+(s/b))) - 2*s )
  else: return -999.

def sumSignificance(partition, h_sig_SR, h_bkg_SR, h_data_SB_noSmooth):
  sum = 0.
  for pair in partition:
    s = 2*h_sig_SR.Integral(pair[0],pair[1]) # bing for reweight half of training signal to the corresponding luminosity 
    b = h_bkg_SR.Integral(pair[0],pair[1])
    d_noSmooth = h_data_SB_noSmooth.Integral(pair[0],pair[1])
    significance = computeSignificance(s,b,d_noSmooth)
    #print h_sig_SR.GetBinCenter(pair[0])-h_bdt_signal_SR.GetBinWidth(pair[0])/2.,significance,b
    if significance>0.: sum += significance*significance
    else: return -999.
  return np.sqrt(sum)

def getResults(h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB, nCats, nBins):

    significance_final = -999.
    partition_final = []
    sig_all = {}

    #1 categories
    if nCats == 1:

        for i in range(1,nBins+1):
            partition = [[i,nBins]]
            significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
            #print h_bdt_signal_SR.GetBinCenter(partition[0][0])-h_bdt_signal_SR.GetBinWidth(partition[0][0])/2,"1. --->",significance,h_bdt_data_SB.Integral(partition[0][0],nBins)
            #print h_bdt_signal_SR.GetBinCenter(partition[0][0])-h_bdt_signal_SR.GetBinWidth(partition[0][0])/2,"1. --->",significance
            sig_all[i] = significance
            
            if significance>significance_final:
                significance_final = significance
                partition_final = partition
        # output = nCats," - Best category: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2,"1. --->",significance_final,"signal total: ",h_bdt_signal_SR.Integral(1,nBins+1),"events:",h_bdt_signal_SR.GetEntries(),"signal cut",h_bdt_signal_SR.Integral(partition_final[0][0],nBins),"smoothed background: ",h_bdt_datamix_SR_weighted_smooth.Integral(partition_final[0][0],nBins),"data: ",h_bdt_data_SB.Integral(partition_final[0][0],nBins),"NEXT BIN: ","smoothed background: ",h_bdt_datamix_SR_weighted_smooth.Integral(partition_final[0][0]+1,nBins),"data: ",h_bdt_data_SB.Integral(partition_final[0][0]+1,nBins)
        cut_value    = h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2
        # print(f"significance_final: {significance_final:.3f}")
        # print(f"h_bdt_signal_SR.GetBinCenter(partition_final[0][0]) : {h_bdt_signal_SR.GetBinCenter(partition_final[0][0]):.3f}")
        # print(f"h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2 : {(h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2):.3f}")
        signal_total = h_bdt_signal_SR.Integral(1,nBins+1)
        events_total = h_bdt_signal_SR.GetEntries()
        signal_cut   = h_bdt_signal_SR.Integral(partition_final[0][0],nBins)
        data         = h_bdt_data_SB.Integral(partition_final[0][0],nBins)
        smoothed_background = h_bdt_datamix_SR_weighted_smooth.Integral(partition_final[0][0],nBins)
        
        next_bin_smoothed_background = h_bdt_datamix_SR_weighted_smooth.Integral(partition_final[0][0]+1,nBins)
        next_bin_data = h_bdt_data_SB.Integral(partition_final[0][0]+1,nBins)
        
        '''
        output = nCats," - Best category: ", cut_value, "1. --->", significance_final, \
            "signal total: ", signal_total, \
                "events:", events_total, \
                    "signal cut", signal_cut, \
                        "smoothed background: ", smoothed_background, \
                            "data: ", data, \
                                "NEXT BIN: ", \
                                    "smoothed background: ", next_bin_smoothed_background, \
                                        "data: ", next_bin_data"
        '''
        
        output = (
            f"{nCats} - Best category: {cut_value:.3f} 1. ---> {significance_final:.3f} "
            f"signal total: {signal_total:.3f}, "
            f"events: {events_total:.0f}, "
            f"signal cut: {signal_cut:.3f}, "
            f"smoothed_background: {smoothed_background:.3f}, "
            f"data: {data:.0f}, "
            f"NEXT BIN: "
            f"smoothed_background: {next_bin_smoothed_background:.3f}, "
            f"data: {next_bin_data:.0f}"
        )

        # print(output)
        # print (' '.join(map(str,output)))

    #2 categories
    elif nCats == 2:

        for i in range(1,nBins+1):
            for j in range(i+1,nBins+1):
                partition = [[1,i],[j,nBins]]
                if abs(i-j)==1:
                    #print partition
                    significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                    if significance>significance_final:
                        significance_final = significance
                        partition_final = partition

        output = nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2,h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2,"1. --->",significance_final
        # print (' '.join(map(str,output)))

    #3 categories
    elif nCats == 3:

        for i in range(1,nBins+1):
            for j in range(i+1,nBins+1):
                for k in range(j+1,nBins+1):
                    partition = [[1,i],[j,k-1],[k,nBins]]
                    if abs(i-j)==1:
                        #print partition
                        significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                        if significance>significance_final:
                            significance_final = significance
                            partition_final = partition

        output = nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, "1. --->",significance_final
        print (' '.join(map(str,output)))

    #4 categories
    elif nCats == 4:

        for i in range(1,nBins+1):
            for j in range(i+1,nBins+1):
                for k in range(j+1,nBins+1):
                    for d in range(k+1,nBins+1):
                        partition = [[1,i],[j,k-1],[k,d-1],[d,nBins]]
                        if abs(i-j)==1:
                            #print partition
                            significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                            if significance>significance_final:
                                significance_final = significance
                                partition_final = partition

        output = nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, "1. --->",significance_final
        print (' '.join(map(str,output)))

    #5 categories
    elif nCats == 5:
        for i in range(1,nBins+1):
            for j in range(i+1,nBins+1):
                for k in range(j+1,nBins+1):
                    for d in range(k+1,nBins+1):
                        for f in range(d+1,nBins+1):
                            partition = [[1,i],[j,k-1],[k,d-1],[d,f-1],[f,nBins]]
                            if abs(i-j)==1:
                            #print partition
                                significance = sumSignificance(partition, h_bdt_signal_SR, h_bdt_datamix_SR_weighted_smooth, h_bdt_data_SB)
                                if significance>significance_final:
                                    significance_final = significance
                                    partition_final = partition

        output = nCats," - Best categories: ",h_bdt_signal_SR.GetBinCenter(partition_final[0][0])-h_bdt_signal_SR.GetBinWidth(partition_final[0][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[1][0])-h_bdt_signal_SR.GetBinWidth(partition_final[1][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[2][0])-h_bdt_signal_SR.GetBinWidth(partition_final[2][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[3][0])-h_bdt_signal_SR.GetBinWidth(partition_final[3][0])/2, h_bdt_signal_SR.GetBinCenter(partition_final[4][0])-h_bdt_signal_SR.GetBinWidth(partition_final[4][0])/2, "1. --->",significance_final
        print (' '.join(map(str,output)))

    else:
        print ("Number of categories not supported, choose: 1, 2, 3, 4, 5, 6, 7, 8 or 9!")
        sys.exit()

    
    return partition_final, significance_final, output, sig_all



def main():

    #analyzer_cfg.sig_names
    # analyzer_cfg = AC.Analyzer_Config('inclusive',options.year)
    analyzer_cfg = AC.Analyzer_Config('inclusive', options.year, options.region, options.mva)

    file = {}
    file['SR'] = open_file(file_out + '/' + name_SR)
    file['CR'] = open_file(file_out + '/' + name_CR)


    signal_region = ['1sigma', '1P5sigma', '2sigma', '3sigma','all']
    hist_SR = {}
    hist_SR_smooth = {}

    hist_CR = {}
    hist_CR_smooth = {}

    hist_signal = {}

    for r in signal_region:
        hist_SR[r] = {}
        hist_SR_smooth[r] = {}

        hist_signal[r] = {}

    y = {}
    y_list = ['significance', 'boundary', 'Nsignal', 'Nbackground']
    for l in y_list:
        y[l] = {}

    for sig in analyzer_cfg.sig_names:

        #hist_signal[sig] = get_hist(file['SR'], "mvaVal_"+sig+"_"+sig)

        ###### smoothing his
        #hist_SR[sig], hist_SR_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        #hist_SR_1sigma[sig], hist_SR_1sigma_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        #hist_SR_1P5sigma[sig], hist_SR_1P5sigma_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        #hist_SR_2sigma[sig], hist_SR_2sigma_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        #hist_SR_3sigma[sig], hist_SR_3sigma_smooth[sig] = GetHist(file['SR'], "mvaVal_"+sig+"_DYJetsToLL")
        for r in signal_region:
            if r == 'all':
                hist_signal[r][sig] = get_hist(file['SR'], "mvaVal_"+sig+"_"+sig)
                hist_SR[r][sig], hist_SR_smooth[r][sig] = GetHist(file['SR'], ["mvaVal_"+sig+"_DYJetsToLL", "mvaVal_"+sig+"_DYGto2LG"])
            else:
                hist_signal[r][sig] = get_hist(file['SR'], "mvaVal_"+r+"_"+sig+"_"+sig)
                hist_SR[r][sig], hist_SR_smooth[r][sig] = GetHist(file['SR'], ["mvaVal_"+r+"_"+sig+"_DYJetsToLL", "mvaVal_"+r+"_"+sig+"_DYGto2LG"])
        hist_CR[sig], hist_CR_smooth[sig] = GetHist(file['CR'], "mvaVal_"+sig+"_Data")

        if options.plot:
            for r in signal_region:
                canv_SR = compare(hist_SR[r][sig], hist_SR_smooth[r][sig],sig)
                canv_SR.SaveAs(options.outDir+"/"+sig+"_"+r+"_DYJetsToLL_DYGto2LG_SR.pdf")
                # canv_SR.SaveAs(options.outDir+"/"+sig+"_"+r+"_DYJetsToLL_SR.png")
            
            canv_CR = compare(hist_CR[sig], hist_CR_smooth[sig],sig)
            canv_CR.SaveAs(options.outDir+"/"+sig+"_data_CR.pdf")
            # canv_CR.SaveAs(options.outDir+"/"+sig+"_data_CR.png")

        if options.doOpt:
            ###### category optimization
            nCats = options.nCats
            nBins = hist_signal['all'][sig].GetNbinsX()

            partition_final = {}
            significance_final = {}
            output = {}
            significance_all = {}

            partition_final_up = {}
            significance_final_up = {}
            output_up = {}
            significance_all_up = {}
            partition_final_dn = {}
            significance_final_dn = {}
            output_dn = {}
            significance_all_dn = {}

            for r in signal_region:

                # print(r)
                partition_final[r], significance_final[r], output[r], significance_all[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][0], hist_CR[sig], nCats, nBins)
                partition_final_up[r], significance_final_up[r], output_up[r], significance_all_up[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][1], hist_CR[sig], nCats, nBins)
                partition_final_dn[r], significance_final_dn[r], output_dn[r], significance_all_dn[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][2], hist_CR[sig], nCats, nBins)
                #partition_final[r], significance_final[r], output[r], significance_all[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][1], hist_CR[sig], nCats, nBins)
                #partition_final[r], significance_final[r], output[r], significance_all[r] = getResults(hist_signal[r][sig], hist_SR_smooth[r][sig][2], hist_CR[sig], nCats, nBins)

                #print "nSig: ", 

                if not os.path.exists(options.outDir):
                    os.makedirs(options.outDir)

                # 準備 JSON 內容
                def build_json_obj(partitions, sigf, smooth_idx):
                    # boundaries 以各分類的左邊界表示（與原本輸出一致）
                    boundaries = [
                        float(hist_signal[r][sig].GetBinCenter(p[0]) - hist_signal[r][sig].GetBinWidth(p[0]) / 2.0)
                        for p in partitions
                    ]
                    obj = {
                        "nCats": int(nCats),
                        "region": str(r),
                        "mass": str(sig),
                        "significance": round(float(sigf), 3) if sigf != -999. else -999.0,
                        "boundaries": [round(b, 3) for b in boundaries]
                    }
                    # 僅在 1 類別時，補充詳細欄位（對照舊 txt 內容）
                    if nCats == 1 and len(partitions) == 1 and partitions[0][0] >= 1:
                        low_bin = partitions[0][0]
                        signal_total = float(hist_signal[r][sig].Integral(1, nBins+1))
                        events_total = int(hist_signal[r][sig].GetEntries())
                        signal_cut = float(hist_signal[r][sig].Integral(low_bin, nBins))
                        smoothed_background = float(hist_SR_smooth[r][sig][smooth_idx].Integral(low_bin, nBins))
                        data_val = float(hist_CR[sig].Integral(low_bin, nBins))
                        next_bin_smoothed_background = float(hist_SR_smooth[r][sig][smooth_idx].Integral(low_bin+1, nBins))
                        next_bin_data = float(hist_CR[sig].Integral(low_bin+1, nBins))
                        obj["details"] = {
                            "signal_total": round(signal_total, 3),
                            "events": events_total,
                            "signal_cut": round(signal_cut, 3),
                            "smoothed_background": round(smoothed_background, 3),
                            "data": round(data_val, 3),
                            "next_bin_smoothed_background": round(next_bin_smoothed_background, 3),
                            "next_bin_data": round(next_bin_data, 3)
                        }
                    return obj

                # 正常、上偏、下偏三種 JSON
                obj_nom = build_json_obj(partition_final[r],     significance_final[r],     0)
                obj_up  = build_json_obj(partition_final_up[r],  significance_final_up[r],  1)
                obj_dn  = build_json_obj(partition_final_dn[r],  significance_final_dn[r],  2)

                # 寫出 json 檔案（取代原本 .txt）
                out_base = f"{options.outDir}/nCat{int(nCats)}_{r}_{sig}"
                with open(out_base + ".json", "w") as f_nom:
                    json.dump(obj_nom, f_nom, indent=2)
                with open(out_base + "_up.json", "w") as f_up:
                    json.dump(obj_up, f_up, indent=2)
                with open(out_base + "_dn.json", "w") as f_dn:
                    json.dump(obj_dn, f_dn, indent=2)

        if options.sigVSscore:
            
            if nCats == 1:
                # 不要硬編碼 'M5'，改用目前迴圈的 sig；若不存在則用任意可用的 mass 當 fallback
                if sig in hist_signal.get('all', {}):
                    w = hist_signal['all'][sig].GetBinWidth(1)
                else:
                    all_keys = list(hist_signal.get('all', {}).keys())
                    if len(all_keys) == 0:
                        continue
                    w = hist_signal['all'][all_keys[0]].GetBinWidth(1)

                '''
                plt.plot([(i-1)*w for i in significance_all['all'].keys()], significance_all['all'].values(), "o-", markersize=2., label='Significance')
                #plt.plot([(i-1)*w for i in significance_all['all'].keys()], significance_all_up['all'].values(), "o-", c='red')
                #plt.plot([(i-1)*w for i in significance_all['all'].keys()], significance_all_dn['all'].values(), "o-", c='green')
                plt.fill_between([(i-1)*w for i in significance_all['all'].keys()], significance_all_up['all'].values(), significance_all_dn['all'].values(), alpha=0.2, color='red',label=r'Significance$\pm 1\sigma$')
                plt.ylim([0.0, 150.0])
                plt.xlim([0.5, 1.0])
                plt.xlabel('boundary_'+sig)
                plt.ylabel('significance')
                plt.legend(loc='best')
                plt.plot([(partition_final['all'][0][0]-1)*w, (partition_final['all'][0][0]-1)*w], [0., 150], c='black', linestyle='--')
                plt.grid()
                plt.savefig(options.outDir+ '/sigVSscore_' + sig + '.png')
                plt.savefig(options.outDir+ '/sigVSscore_' + sig + '.pdf')
                plt.close('all')
                '''

                if gROOT.FindObject("cc1"):
                    gROOT.FindObject("cc1").Close()
                canv_sigVSscore = TCanvas("cc1", "cc1", 800, 600)
                leftMargin=0.12
                rightMargin=0.02
                topMargin=0.08
                canv_sigVSscore.SetMargin(leftMargin, rightMargin, 0.145, topMargin) #//left//right//bottom//top
                canv_sigVSscore.cd()
                # gStyle.SetPadTickX(1)
                # gStyle.SetPadTickY(1)
                canv_sigVSscore.SetTickx()
                canv_sigVSscore.SetTicky()

                x_low = 0.58
                x_high = 1.01
                num_arr = len(significance_all['all'].values())
                x_arr = array('d')
                y_arr = array('d')
                y_arr_err = array('d')
                count = 0
                burden = 0
                for i in significance_all['all'].keys():
                    if ((i-1)*w - 0.1 < x_low):
                        burden += 1
                    if ((i-1)*w - 0.1 >= x_low) and ((i-1)*w - 0.1 <= 1.):
                        x_arr.append((i-1)*w - 0.1)
                        y_arr.append(significance_all['all'][i])
                        # if ((i-1)*w - 0.1 > 0.94):
                            # print(f'x: {((i-1)*w - 0.1):.3f}')
                            # print(f"y: {significance_all['all'][i]:.3f}")
                        count += 1

                gr = TGraph( count, x_arr, y_arr)
                gr_err = TGraph( 2 * count)

                for i in range(count):
                    gr_err.SetPoint(i, x_arr[i], list(significance_all_up['all'].values())[i+burden])
                    gr_err.SetPoint(2*count-1-i, x_arr[i], list(significance_all_dn['all'].values())[i+burden])
                
                # 先算上沿的最大 y
                y_up = [list(significance_all_up['all'].values())[i+burden] for i in range(count)]
                ymax = max(y_up)

                # y_high = 65.
                y_low = 0.
                gr_err.SetLineColor(0)
                gr_err.SetFillColorAlpha(2,0.3)
                # gr_err.GetHistogram().SetMaximum(1.5 * ymax)
                # gr_err.GetHistogram().SetMinimum(y_low)
                # gr_err.GetHistogram().GetXaxis().SetRangeUser(x_low,x_high)

                gr_err.GetYaxis().SetTitle( f'Significance (AMS) / {w:.3f}' )
                gr_err.GetXaxis().SetTitle( 'BDT Score' )
                gr_err.GetXaxis().CenterTitle(True)
                gr_err.GetYaxis().CenterTitle(True)

                gr_err.GetYaxis().SetTitleColor(1)
                gr_err.GetYaxis().SetTitleFont(42)
                gr_err.GetYaxis().SetTitleSize(0.055)
                gr_err.GetYaxis().SetTitleOffset(0.98)

                gr_err.GetYaxis().SetLabelColor(1)
                gr_err.GetYaxis().SetLabelFont(42)
                gr_err.GetYaxis().SetLabelSize(0.05)
                gr_err.GetYaxis().SetLabelOffset(0.007)


                gr_err.GetXaxis().SetTitleColor(1)
                gr_err.GetXaxis().SetTitleFont(42)
                gr_err.GetXaxis().SetTitleSize(0.055)
                gr_err.GetXaxis().SetTitleOffset(1.2)

                gr_err.GetXaxis().SetLabelSize(0.05)
                gr_err.GetXaxis().SetLabelColor(1)
                gr_err.GetXaxis().SetLabelFont(42)
                gr_err.GetXaxis().SetLabelOffset(0.007)


                SetgStyle()
                gr.SetMarkerSize(0.8)

                gr_err.Draw('AF')

                y_up = [list(significance_all_up['all'].values())[i+burden] for i in range(count)]
                ymax = max(y_up)
                frame = gr_err.GetHistogram()  # 现在能拿到有效的 TH1 框
                frame.SetMaximum(2.0 * ymax)
                frame.SetMinimum(y_low)
                frame.GetXaxis().SetRangeUser(x_low, x_high)
                gPad.Modified()
                gPad.Update()

                Line = TLine((partition_final['all'][0][0]-1)*w - 0.1, y_low, (partition_final['all'][0][0]-1)*w - 0.1, 2.0 * ymax)
                Line.SetLineStyle(7)
                Line.SetLineWidth(3)
                Line.SetLineColor(4)
                Line.Draw("SAME")

                gr.Draw('LP')

                legend = TLegend(0.15,0.62,0.50,0.82)                     # Add a legend near the top right corner
                legend.SetFillColor(0)
                legend.SetLineColor(0)
                legend.SetFillStyle(0)
                legend.SetBorderSize(0)
                legend.SetTextFont(42)
                legend.SetTextSize(0.045)
                legend.SetLineWidth(0)                              # Remove the boundary on the legend
                legend.AddEntry(Line,"Working Point", "l")          # Add the data points, labelled as "Data"
                legend.AddEntry(gr,"Significance", "pl")                      # Add the MC histogram, labelled as "MC"
                legend.AddEntry(gr_err,"Significance #pm 1#sigma", "f") # Add the data points, labelled as "Data"
                legend.Draw("same")                                 # Draw the legend on the plot

                latex_cut = TLatex()
                latex_cut.SetNDC()
                latex_cut.SetTextColor(kBlack)
                latex_cut.SetTextFont(42)
                latex_cut.SetTextSize(0.045)
                latex_cut.SetTextAlign(13)
                latex_cut.DrawLatex(0.24, 0.87, f"m_{{a}} = {sig.lstrip('M')} GeV")

                latex_cut.SetTextSize(0.045)  # 字體大小
                latex_cut.SetTextFont(61)  # 粗體字 CMS 標籤
                latex_cut.DrawLatex(leftMargin, 0.969, "CMS")

                latex_cut.SetTextFont(52)  # 斜體字 Preliminary 標籤
                latex_cut.DrawLatex(leftMargin + 0.09, 0.965, "Preliminary") # + 0.09

                # latex_cut.DrawLatex(0.76,0.97,("138 fb^{-1} (13 TeV)"))
                latex_cut.SetTextAlign(31)
                latex_cut.DrawLatex(1. - rightMargin, 1. - topMargin + 0.01,("170.84 fb^{-1} (13.6 TeV)"))
                
                canv_sigVSscore.SaveAs(options.outDir+ '/sigVScore_' + sig + '.pdf')
                # canv_sigVSscore.SaveAs(options.outDir+ '/sigVSscore_' + sig + '.png')

        #make sigma plots
        if nCats == 1:
            for l in y_list:
                y[l][sig] = []
            for r in signal_region:
                
                s = 2*hist_signal['all'][sig].Integral(partition_final[r][0][0],nBins)
                b = hist_SR_smooth['all'][sig][0].Integral(partition_final[r][0][0],nBins)

                if b > 0.0:
                    y['boundary'][sig].append(hist_signal[r][sig].GetBinCenter(partition_final[r][0][0])-hist_signal[r][sig].GetBinWidth(partition_final[r][0][0])/2)
                    #y['significance'][sig].append(significance_final[r])

                    y['significance'][sig].append(np.sqrt((2*(s+b)*math.log(1+(s/b))) - 2*s))
                    y['Nsignal'][sig].append(hist_signal['all'][sig].Integral(partition_final[r][0][0],nBins))
                    y['Nbackground'][sig].append(hist_SR_smooth['all'][sig][0].Integral(partition_final[r][0][0],nBins))
                else:
                    b = hist_SR_smooth['all'][sig][0].Integral(partition_final['all'][0][0],nBins)

                    y['boundary'][sig].append(hist_signal['all'][sig].GetBinCenter(partition_final['all'][0][0])-hist_signal['all'][sig].GetBinWidth(partition_final['all'][0][0])/2)

                    y['significance'][sig].append(np.sqrt((2*(s+b)*math.log(1+(s/b))) - 2*s))
                    y['Nsignal'][sig].append(hist_signal['all'][sig].Integral(partition_final['all'][0][0],nBins))
                    y['Nbackground'][sig].append(hist_SR_smooth['all'][sig][0].Integral(partition_final['all'][0][0],nBins))



    if options.sigma and nCats == 1:
        for l in y_list:
            colors = ['red','orange','blue','purple','darkturquoise', 'yellow', 'pink', 'green', 'black', 'gray', 'violet', 'gold', 'cyan', 'brown']

            plt.xlabel('signal mass region')
            plt.ylabel(l)

            for sig in analyzer_cfg.sig_names:

                plt.plot(signal_region, y[l][sig], c=colors[analyzer_cfg.sig_names.index(sig)], label=sig)


            plt.grid()
            plt.legend()
            plt.savefig(options.outDir+ '/' + l + '.pdf')
            # plt.savefig(options.outDir+ '/' + l + '.png')
            plt.close('all')


main()
