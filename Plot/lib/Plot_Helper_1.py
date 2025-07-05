
import os
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *

#####################################################################

def LoadNtuples(ana_cfg):
    ntuples = {}
    for sample in ana_cfg.samp_names:
        ntuples[sample] = TChain("passedEvents","chain_" + sample)
	ntuples[sample]. Add(ana_cfg.sample_loc + '/ALP_%s.root' %sample)
    return ntuples


def MakeStack(histos, ana_cfg, var_name):
    stacks = {}
    stacks['all']  = THStack("h_stack_"+var_name, var_name)
    stacks['sig']  = THStack("h_stack_"+var_name, var_name)
    stacks['bkg']  = THStack("h_stack_"+var_name, var_name)
    for sample in ana_cfg.samp_names:
	stacks[sample]  = THStack("h_stack_"+var_name, var_name)

    for sample in ana_cfg.bkg_names:
        stacks['bkg'].Add(histos[sample])
        stacks['all'].Add(histos[sample])

    for sample in ana_cfg.sig_names:
        stacks['sig'].Add(histos[sample])
	stacks[sample].Add(histos[sample])
	#if var_name in ['Z_m', 'H_m', 'ALP_m']:
	    #stacks['all'].Add(histos[sample])
    return stacks


def CreateCanvas(canv_name):
    canv = TCanvas(canv_name, canv_name, 600,650)
    return canv


def MakeLumiLabel(lumi):
    tex = TLatex()
    tex.SetTextSize(0.035)
    tex.SetTextAlign(31)
    tex.DrawLatexNDC(0.90, 0.91, '%s fb^{-1} (13 TeV)' %lumi)
    return tex


def MakeCMSDASLabel():
    tex = TLatex()
    tex.SetTextSize(0.03)
    tex.DrawLatexNDC(0.15, 0.85, '#scale[1.5]{CMSDAS} H To Z + ALP')
    return tex


def ScaleSignal(plt_cfg, stack_sig, hist_data, var_name):
    #sig_hist = stack_sig.GetStack().Last()
    sig_hist = hist_data
    #sig_hist.Scale(plt_cfg.sig_scale)
   
    #print "all data: " + str(hist_data['data'].GetEntries())
    #print "all sig: " + str(sig_hist.GetEntries()*plt_cfg.sig_weight)

    #sig_hist.Scale(hist_data['data'].GetEntries()/(sig_hist.GetEntries()*plt_cfg.sig_weight))
    #sig_hist.SetLineColor(kRed)
    sig_hist.SetLineWidth(3)
    sig_hist.SetFillStyle(0)

    sig_hist.GetXaxis().SetTitle(var_name)
    sig_hist.GetXaxis().SetTitleSize(0.5)
    sig_hist.GetYaxis().SetTitle('Events / %.2f' %sig_hist.GetBinWidth(1))
    #sig_hist.GetYaxis().SetTitleOffset(0.20)
    return sig_hist


def MakeRatioPlot(h_data, h_MC, var_name):

    ratio_plot = TGraphAsymmErrors()
    ratio_plot.Divide(h_data, h_MC, "pois")
    ratio_plot.SetName("ratiograph_" + var_name)
    ratio_plot.SetMinimum(0.2)
    ratio_plot.SetMaximum(3.4)
    ratio_plot.SetMarkerStyle(20)

    #ratio_plot.GetXaxis().SetRangeUser( h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax() )
    ratio_plot.GetXaxis().SetLimits( h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax() )
    ratio_plot.GetXaxis().SetLabelSize(0.18)
    ratio_plot.GetXaxis().SetTitle(var_name)
    ratio_plot.GetXaxis().SetTitleSize(0.15)
    ratio_plot.GetXaxis().SetTitleOffset(0.7)

    ratio_plot.GetYaxis().SetNdivisions(505)
    ratio_plot.GetYaxis().SetLabelSize(0.1)
    ratio_plot.GetYaxis().SetTitle("data/MC")
    ratio_plot.GetYaxis().SetTitleSize(0.15)
    ratio_plot.GetYaxis().SetTitleOffset(0.2)

    return ratio_plot


def MakeLegend(plt_cfg, histos, scaled_signal):
    legend = TLegend(0.5,0.6,0.9,0.9)
    legend.SetNColumns(3)

    legend.AddEntry(histos["data"], "data")

    for sample in plt_cfg.ana_cfg.sig_names:
    	legend.AddEntry(scaled_signal[sample], sample )#+ " signal X%d" %plt_cfg.sig_scale )
    #legend.AddEntry(scaled_signal, "signal X%d" %plt_cfg.sig_scale)

    for sample in plt_cfg.ana_cfg.bkg_names:
        legend.AddEntry(histos[sample], sample )
    return legend


def DrawOnCanv(canv, var_name, plt_cfg, stacks, histos, scaled_sig, ratio_plot, legend, lumi_label, cms_label, logY):
    canv.SetBottomMargin(0.01)
    canv.cd()
    upper_pad = TPad("upperpad_"+var_name, "upperpad_"+var_name, 0,0.2, 1,1)
    upper_pad.SetBottomMargin(0.05)
    upper_pad.Draw()
    upper_pad.cd()
    if logY:
    	upper_pad.SetLogy()
    	stacks['all'].SetMinimum(1e-3)
        stacks['all'].SetMaximum(1e8)


    histos['data'].Draw('PE')
    histos['data'].GetXaxis().SetLabelSize(0)
    histos['data'].GetYaxis().SetTitle('Events / (%.2f GeV)' %histos['data'].GetBinWidth(1))
    histos['data'].GetYaxis().SetTitleSize(0.045)
    stacks['all'].Draw('HISTSAME')
    #stacks['all'].GetHistogram().GetXaxis().SetLabelSize(0)
    #stacks['all'].GetHistogram().GetYaxis().SetTitle('Events / (%.2f GeV)' %stacks['all'].GetHistogram().GetBinWidth(1))
    #stacks['all'].GetHistogram().GetYaxis().SetTitleSize(0.045)
    histos['data'].SetMarkerStyle(20)
    histos['data'].Draw('SAMEPE')
    for sample in plt_cfg.ana_cfg.sig_names:
	scaled_sig[sample].Draw('HISTSAME')

    legend.Draw()
    cms_label.DrawLatexNDC(0.15, 0.85, '#scale[1.5]{CMS} H To Z + ALP')
    cms_label.Draw('same')
    lumi_label.DrawLatexNDC(0.90, 0.91, '%s fb^{-1} (13 TeV)' %plt_cfg.lumi)
    lumi_label.Draw('same')

    canv.cd()
    lower_pad = TPad("lowerpad_"+var_name, "lowerpad_"+var_name, 0, 0.05, 1,0.22)
    lower_pad.SetTopMargin(0.)
    lower_pad.SetBottomMargin(0.25)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    ratio_plot.Draw("AP")


def SaveCanvPic(canv, save_dir, save_name):
    canv.cd()
    #canv.SaveAs(save_dir + '/' + save_name + '.pdf')
    canv.SaveAs(save_dir + '/' + save_name + '.png')
