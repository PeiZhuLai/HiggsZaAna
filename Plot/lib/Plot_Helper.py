import os
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *
import ROOT
import numpy as np

import CMS_lumi
from array import array

import math
#####################################################################

def ScaleBkgToData(histos, ana_cfg):
    """
    Scale background histograms to match the total yield of the data histogram.
    
    Parameters:
    - histos (dict): Dictionary of histograms for a given variable, keyed by sample name.
    - ana_cfg (Analyzer_Config): Configuration object containing sample names, including bkg_names and data.
    
    Returns:
    - None: Modifies the histograms in place.
    """
    # Get the data histogram
    data_hist = histos.get('Data')
    if not data_hist:
        print("[ScaleBkgToData] Warning: No 'Data' histogram found. Skipping scaling.")
        return
    
    # Calculate total data yield
    data_integral = data_hist.Integral()
    if data_integral <= 0:
        print("[ScaleBkgToData] Warning: Data integral is zero or negative. Skipping scaling.")
        return
    
    # Calculate total background yield
    bkg_integral = 0.0
    for sample in ana_cfg.bkg_names:
        if sample in histos:
            bkg_integral += histos[sample].Integral()
    
    if bkg_integral <= 0:
        print("[ScaleBkgToData] Warning: Background integral is zero or negative. Skipping scaling.")
        return
    
    # Compute scaling factor
    scale_factor = data_integral / bkg_integral
    print(f"[ScaleBkgToData] Scaling factor: {scale_factor:.4f} (Data integral: {data_integral:.2f}, Bkg integral: {bkg_integral:.2f})")
    
    # Apply scaling factor to each background histogram
    for sample in ana_cfg.bkg_names:
        if sample in histos:
            histos[sample].Scale(scale_factor)
            print(f"[ScaleBkgToData] Scaled histogram for {sample} by {scale_factor:.4f}")


def LoadNtuples(ana_cfg):
    ntuples = {}
    for sample in ana_cfg.samp_names:
        if "M" in sample: 
            ntuples[sample] = TChain("test","chain_" + sample)
            ntuples[sample].Add(ana_cfg.sample_loc + '/ALP_%s/run3.root' %sample)
        else: 
            ntuples[sample] = TChain("inclusive","chain_" + sample)
            ntuples[sample].Add(ana_cfg.sample_loc + '/%s/run3.root' %sample)
    return ntuples


def MakeStack(histos, ana_cfg, var_name):
    stacks = {}
    stacks['all']  = THStack("h_stack_all_"+var_name, "all_"+var_name)
    stacks['sig']  = THStack("h_stack_sig_"+var_name, "sig_"+var_name)
    stacks['bkg']  = THStack("h_stack_bkg_"+var_name, "bkg_"+var_name)

    for sample in ana_cfg.samp_names:
        stacks[sample] = THStack("h_stack_"+sample+"_"+var_name, sample+"_"+var_name)

    for sample in ana_cfg.bkg_names:
        stacks['bkg'].Add(histos[sample])
        stacks['all'].Add(histos[sample])

    for sample in ana_cfg.sig_names:
        stacks['sig'].Add(histos[sample])
        stacks[sample].Add(histos[sample])

    return stacks

def CreateCanvas(canv_name):
    canv = TCanvas(canv_name, canv_name, 800, 800)
    return canv

def MakeLumiLabel(lumi):
    tex = TLatex()
    tex.SetTextSize(0.035)
    tex.SetTextAlign(31)
    tex.DrawLatexNDC(0.9, 0.91, '%s fb^{-1} (13.6 TeV)' %lumi)
    return tex

def MakeCMSDASLabel():
    #tex = TLatex()
    #tex.SetTextSize(0.03)
    #tex.DrawLatexNDC(0.12, 0.85, '#scale[1.5]{CMSDAS} H To Z + ALP')
    #return tex

    onTop=False
    text='#bf{CMS} #scale[0.75]{#it{Simulation Preliminary}  H#rightarrow#gamma#gamma}'
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.1, 0.85 if not onTop else 0.93, text)
    return latex

def ScaleSignal(plt_cfg, stack_sig, hist_data, var_name):
    sig_hist = hist_data
    sig_hist.SetLineWidth(3)
    sig_hist.SetFillStyle(0)

    sig_hist.GetXaxis().SetTitle(var_name)
    sig_hist.GetXaxis().SetTitleSize(0.5)
    sig_hist.GetYaxis().SetTitle('Events / %.2f' %sig_hist.GetBinWidth(1))
    return sig_hist

def MakeRatioPlot(h_data, h_MC, var_name):
    ratio_plot = TGraphAsymmErrors()
    ratio_plot.Divide(h_data, h_MC, "pois")
    ratio_plot.SetName("ratiograph_" + var_name)
    ratio_plot.SetMinimum(0.4)
    ratio_plot.SetMaximum(1.6)
    ratio_plot.SetMarkerStyle(20)

    ratio_plot.GetXaxis().SetLimits( h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax() )
    ratio_plot.GetXaxis().SetLabelSize(0.14)
    ratio_plot.GetXaxis().SetTitle(var_name)
    ratio_plot.GetXaxis().SetTitleSize(0.16)
    ratio_plot.GetXaxis().SetTitleOffset(0.65)
    ratio_plot.GetXaxis().SetTickLength(0.08)

    ratio_plot.GetYaxis().SetNdivisions(505)
    ratio_plot.GetYaxis().SetLabelSize(0.14)
    ratio_plot.GetYaxis().SetTitle("Data / SM")
    ratio_plot.GetYaxis().CenterTitle(True)
    ratio_plot.GetYaxis().SetTitleSize(0.14)
    ratio_plot.GetYaxis().SetTitleOffset(0.535)

    return ratio_plot

def MakeLegend(plt_cfg, histos, scaled_signal):
    legend = TLegend(0.55,0.65,0.85,0.86)
    legend.SetNColumns(1)
    legend.AddEntry(histos["Data"], "Data", "PE")
    for sample in plt_cfg.ana_cfg.sig_names:
    #for sample in ["M1","M10","M20","M30"]:
        legend.AddEntry(scaled_signal[sample], sample, "l")

    for sample in plt_cfg.ana_cfg.bkg_names:
        legend.AddEntry(histos[sample], sample )

    return legend


def Get_StatUnc(hist):
    TH1.Sumw2

    BinTotal = hist.GetNbinsX()
    WidthBin = hist.GetBinWidth(1)

    XMean=[]
    YMean=[]
    YMeanNorm=[]
    X_ErrH=[]
    X_ErrL=[]
    Y_ErrH=[]
    Y_ErrL=[]
    YNorm_ErrH=[]
    YNorm_ErrL=[]

    xaxis = hist.GetXaxis()

    graph = TGraphAsymmErrors()

    for iBin in range(1, BinTotal+1):
        N = hist.GetBinContent(iBin)
        NErr = hist.GetBinError(iBin)
        XMean.append(xaxis.GetBinCenter(iBin))
        YMean.append(N)
        YMeanNorm.append(1.0)
        X_ErrH.append(0.5*WidthBin)
        X_ErrL.append(0.5*WidthBin)
        Y_ErrH.append(NErr)
        Y_ErrL.append(NErr)

        if N>0:
            errNorm = NErr/N
        else:
            errNorm = 0.
        YNorm_ErrH.append(errNorm)
        YNorm_ErrL.append(errNorm)


    #print YMean
    #print Y_ErrH

    #print 'stat:'
    #print 'high:'
    #print Y_ErrH
    #print Y_ErrL

    graph = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMean), np.array(X_ErrL), np.array(X_ErrH), np.array(Y_ErrL), np.array(Y_ErrH))
    graph_norm = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMeanNorm), np.array(X_ErrL), np.array(X_ErrH), np.array(YNorm_ErrH), np.array(YNorm_ErrL))

    return graph, graph_norm


def Get_SysUnc(hist, hist_up, hist_dn):
    TH1.Sumw2

    BinTotal = hist.GetNbinsX()
    WidthBin = hist.GetBinWidth(1)

    XMean=[]
    YMean=[]
    YMeanNorm=[]
    X_ErrH=[]
    X_ErrL=[]
    YSys_ErrH=[]
    YSys_ErrL=[]
    YSysNorm_ErrH=[]
    YSysNorm_ErrL=[]

    xaxis = hist.GetXaxis()

    graph_sys = TGraphAsymmErrors()

    for iBin in range(1, BinTotal+1):
        N = hist.GetBinContent(iBin)
        NErr = hist.GetBinError(iBin)

        dNp = hist_up.GetBinContent(iBin) - N
        dNm = hist_dn.GetBinContent(iBin) - N

        dplus = dNp if (dNp>dNm) else dNm
        if dplus < 0: dplus = -1.0*dplus
        dminus = dNp if (dNp<dNm) else dNm
        if dminus > 0: dminus = -1.0*dminus
        dminus = -1.0*dminus

        XMean.append(xaxis.GetBinCenter(iBin))
        YMean.append(N)
        YMeanNorm.append(1.0)
        X_ErrH.append(0.5*WidthBin)
        X_ErrL.append(0.5*WidthBin)
        YSys_ErrH.append(dplus)
        YSys_ErrL.append(dminus)

        if N>0:
            errNormH = dplus/N
            errNormL = dminus/N
        else:
            errNormH = 0.
            errNormL = 0.
        YSysNorm_ErrH.append(errNormH)
        YSysNorm_ErrL.append(errNormL)

    graph_sys = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMean), np.array(X_ErrL), np.array(X_ErrH), np.array(YSys_ErrH), np.array(YSys_ErrL))
    #graph_sys_norm = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMeanNorm), np.array(X_ErrL), np.array(X_ErrH), np.array(YSysNorm_ErrH), np.array(YSysNorm_ErrL))

    #print "high:"
    #print YSys_ErrH
    #print "low:"
    #print YSys_ErrL

    return YSys_ErrH, YSys_ErrL, YSysNorm_ErrH, YSysNorm_ErrL

def Total_Unc(hist_norm, hist_sys, analyzer_cfg):
    TH1.Sumw2
    hist_norm = hist_norm.GetStack().Last()

    stacks_sys = {}
    for sys in analyzer_cfg.sys_names:
        stacks_sys['bkgSys_'+sys]  = THStack("h_stack_bkgSys_"+sys, "bkgSys_"+sys)

    
    for sys in analyzer_cfg.sys_names:
        for sample in analyzer_cfg.bkg_names:
            stacks_sys['bkgSys_'+sys].Add(hist_sys[sample][sys])

    #print stacks_sys

    YSys_ErrH = {}
    YSys_ErrL = {}
    YSysNorm_ErrH = {}
    YSysNorm_ErrL = {}
    for i in range(len(analyzer_cfg.sys_names)//2):
        #print stacks_sys['bkgSys_'+analyzer_cfg.sys_names[2*i]].GetStack().Last()
        YSys_ErrH[i], YSys_ErrL[i], YSysNorm_ErrH[i], YSysNorm_ErrL[i] = Get_SysUnc(hist_norm, stacks_sys['bkgSys_'+analyzer_cfg.sys_names[2*i]].GetStack().Last(), stacks_sys['bkgSys_'+analyzer_cfg.sys_names[2*i+1]].GetStack().Last())


    BinTotal = hist_norm.GetNbinsX()
    WidthBin = hist_norm.GetBinWidth(1)

    XMean=[]
    YMean=[]
    YMeanNorm=[]
    X_ErrH=[]
    X_ErrL=[]
    Y_ErrH=[]
    Y_ErrL=[]
    YNorm_ErrH=[]
    YNorm_ErrL=[]

    xaxis_total = hist_norm.GetXaxis()

    graph_Total = TGraphAsymmErrors()
    graph_norm_Total = TGraphAsymmErrors()

    for iBin in range(1, BinTotal+1):
        N = hist_norm.GetBinContent(iBin)
        NErr = hist_norm.GetBinError(iBin)
        XMean.append(xaxis_total.GetBinCenter(iBin))
        YMean.append(N)
        YMeanNorm.append(1.0)
        X_ErrH.append(0.5*WidthBin)
        X_ErrL.append(0.5*WidthBin)
        
        NErrH = 0.0
        NErrL = 0.0
        for i in range(len(analyzer_cfg.sys_names)//2):
            NErrH = NErrH + YSys_ErrH[i][iBin-1] * YSys_ErrH[i][iBin-1] 
            NErrL = NErrL + YSys_ErrL[i][iBin-1] * YSys_ErrL[i][iBin-1]

        NErrH = np.sqrt(NErrH + NErr*NErr)
        NErrL = np.sqrt(NErrL + NErr*NErr)

        Y_ErrH.append(NErrH)
        Y_ErrL.append(NErrL)

        if N>0:
            errNormH = NErrH/N
            errNormL = NErrL/N
        else:
            errNormH = 0.
            errNormL = 0.
        YNorm_ErrH.append(errNormH)
        YNorm_ErrL.append(errNormL)

    graph_Total = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMean), np.array(X_ErrL), np.array(X_ErrH), np.array(Y_ErrL), np.array(Y_ErrH))
    graph_norm_Total = TGraphAsymmErrors(BinTotal, np.array(XMean), np.array(YMeanNorm), np.array(X_ErrL), np.array(X_ErrH), np.array(YNorm_ErrH), np.array(YNorm_ErrL))

    return [graph_Total, graph_norm_Total]

def Draw_unc(graph, color):
    graph.SetFillColor(color)
    graph.SetFillStyle(3001)
    graph.SetLineColor(color)
    #graph.Draw("SAME2")
    graph.Draw("SAME2")

# def Draw_unc(graph, color):
#     # gStyle.SetHatchesSpacing(1)
#     gStyle.SetHatchesLineWidth(2)
#     graph.SetFillColor(1)
#     graph.SetFillStyle(3004)
#     graph.SetLineColor(0)
#     graph.Draw("SAME2")

def DrawOnCanv(canv, var_name, plt_cfg, stacks, histos, scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, bdtCut, mA, logY):

    canv.SetBottomMargin(0.012)
    canv.cd()

    United_pad_LeftMargin = 0.16 
    United_pad_RightMargin = 0.05   
    United_pad_BottomMargin = 0.19
    United_pad_TopMargin = 0.085

    #upper_pad = TPad("upperpad_"+var_name, "upperpad_"+var_name, 0,0.2, 1,1)
    upper_pad = TPad("upperpad_"+var_name, "upperpad_"+var_name, 0,0.25, 1,1)
    upper_pad.SetLeftMargin(United_pad_LeftMargin)
    upper_pad.SetRightMargin(United_pad_RightMargin)
    upper_pad.SetBottomMargin(United_pad_BottomMargin)
    upper_pad.SetTopMargin(United_pad_TopMargin)
    upper_pad.Draw()
    upper_pad.cd()
    upper_pad.SetTickx(1)
    upper_pad.SetTicky(1)
    canv.SetTickx()
    canv.SetTicky()

    if logY:
        upper_pad.SetLogy()
        stacks['all'].SetMinimum(1e-2)
        stacks['all'].SetMaximum(9e9)

        histos['Data'].SetMinimum(1e-2)
        histos['Data'].SetMaximum(9e9)

    if histos['Data'].GetMaximum() > stacks['all'].GetMaximum():
        h_max = histos['Data'].GetMaximum()
    else:
        h_max = stacks['all'].GetMaximum()
    if h_max < stacks['sig'].GetMaximum():
        h_max = stacks['sig'].GetMaximum()

    histos['Data'].SetMaximum(h_max*1.4)
    stacks['all'].SetMaximum(h_max*1.4)

    histos['Data'].Draw('PE')
    histos['Data'].GetXaxis().SetLabelSize(0)
    histos['Data'].GetXaxis().SetTitleOffset(0.95)
    histos['Data'].GetYaxis().SetLabelSize(0.06)
    
    if var_name in ["H_m","ALP_m","Z_m"]:
        histos['Data'].GetYaxis().SetTitle('Events / (%.2f GeV)' %histos['Data'].GetBinWidth(1))
    else:
        histos['Data'].GetYaxis().SetTitle('Events')
    histos['Data'].GetYaxis().SetTitleSize(0.07)
    histos['Data'].GetYaxis().SetTitleFont(42)
    histos['Data'].GetYaxis().SetTitleOffset(1.15)
    stacks['all'].Draw('HISTSAME')

    if (var_name.split("_")[-1] in plt_cfg.ana_cfg.sig_names):
        scaled_sig[var_name.split("_")[-1]].Draw('HISTSAME')
    else:
        if bdtCut:
            scaled_sig[mA].Draw('HISTSAME')
        else:
            line_styles = [1, 5, 7, 9]
            # Run2 Paper
            # for sample, i in zip(["M1", "M10", "M20", "M30"], range(4)):
            for sample, i in zip(["M5", "M15", "M30"], range(3)):
                scaled_sig[sample].SetLineStyle(line_styles[i])
                scaled_sig[sample].SetLineWidth(4)
                scaled_sig[sample].Draw('HISTSAME')

    histos['Data'].SetMarkerStyle(20)
    histos['Data'].SetLineWidth(2)
    histos['Data'].GetXaxis().SetTickLength(0.04)
    
    ### Draw the uncertainties
    global stat_err, stat_err_norm
    stat_err,  stat_err_norm= Get_StatUnc(stacks['bkg'].GetStack().Last())

    # Draw_unc(total_unc[0], kGray+10)
    # Draw_unc(stat_err, kRed-10)
    # Draw_unc(total_unc[0], TColor.GetColor("#324376")) # Navy Blue 
    Draw_unc(stat_err, TColor.GetColor("#F76C5E")) # orange

    histos['Data'].Draw('SAMEPE')
    histos['Data'].Draw("AXIS SAME")

    if var_name.split("_")[-1] in plt_cfg.ana_cfg.sig_names:
        # legend.Clear()
        # legend.AddEntry(histos["Data"], "Data", "PE")
        # for sample_bkg in plt_cfg.ana_cfg.bkg_names:
        #     legend.AddEntry(histos[sample_bkg], "Z \rightarrow \ell^{+}\ell^{-}", "f")
        #     # legend.SetHeader("#splitline{Title on top line}{Title on second line}")

        # legend.AddEntry(stat_err,"Stat. uncertainty","f")
        # legend.AddEntry(total_unc[0],"Syst. uncertainty","f")
        # legend.AddEntry(scaled_sig[var_name.split("_")[-1]], r"m_{a} = %s GeV" % (var_name.split("_")[-1].lstrip("M")), "l" )

        # legend.SetBorderSize(0)
        # legend.SetTextFont(42)
        # legend.SetTextSize(0.045)
        # legend.SetFillColor(0)
        # legend.Draw()

        legend_1 = TLegend(0.18, 0.64, 0.48, 0.88)

        ROOT.SetOwnership(legend_1, False)
        legend_1.AddEntry(histos["Data"], "Data", "PE")

        legend_1.AddEntry(histos["DYJetsToLL"], r"Z + jets", "f")
        legend_1.AddEntry(histos["DYGto2LG"], r"Z + \gamma", "f")

        # for sample_bkg in plt_cfg.ana_cfg.bkg_names:
        #     legend_1.AddEntry(histos[sample_bkg], r"Z + jets", "f")
        
        # Statistic 
        legend_1.AddEntry(stat_err,"Stat. Uncer.","f")
        
        # Total Uncertainty
        # legend_1.AddEntry(total_unc[0],"Syst. Uncer.","f")

        legend_2 = TLegend(0.43, 0.80, 0.73, 0.88)
        ROOT.SetOwnership(legend_2, False)
        legend_2.AddEntry(scaled_sig[var_name.split("_")[-1]], r"m_{a} = %s GeV" % (var_name.split("_")[-1].lstrip("M")), "l" )

        legend_1.SetBorderSize(0)
        legend_1.SetFillStyle(0)
        legend_1.SetFillColor(0)
        legend_1.SetTextFont(42)
        legend_1.SetTextSize(0.05)
        legend_1.Draw("SAME")
        

        legend_2.SetBorderSize(0)
        legend_2.SetFillStyle(0)
        legend_2.SetFillColor(0)
        legend_2.SetTextFont(42)
        legend_2.SetTextSize(0.05)
        legend_2.Draw("SAME")


    else:
        legend_1 = TLegend(0.18, 0.64, 0.48, 0.88)
        ROOT.SetOwnership(legend_1, False)
        legend_1.AddEntry(histos["Data"], "Data", "PE")
        bkg_labels = {"DYJetsToLL": r"Z + jets",                    
                      "DYGto2LG": r"Z + \gamma"}    
            
        for sample_bkg in plt_cfg.ana_cfg.bkg_names:            
            legend_1.AddEntry(histos[sample_bkg], bkg_labels.get(sample_bkg, sample_bkg), "f")
            
        legend_1.AddEntry(stat_err,"Stat. Uncer.","f")
        # legend_1.AddEntry(total_unc[0],"Syst. Uncer.","f")
        # legend_1.AddEntry(total_unc[0],"Uncertainty","f")

        legend_2 = TLegend(0.53, 0.64, 0.83, 0.88)
        ROOT.SetOwnership(legend_2, False)
        if bdtCut:
            legend_2.AddEntry(scaled_sig[mA], r"m_{a} = %s GeV" % (mA.lstrip("M")), "l")
        else:
            # for s in ["M1","M10","M20","M30"]:
            for s in ["M5","M15","M30"]:
                legend_2.AddEntry(scaled_sig[s], r"m_{a} = %s GeV" % (s.lstrip("M")), "l")

        legend_1.SetBorderSize(0)
        legend_1.SetFillStyle(0)
        legend_1.SetFillColor(0)
        legend_1.SetTextFont(42)
        legend_1.SetTextSize(0.05)
        legend_1.Draw("SAME")
        

        legend_2.SetBorderSize(0)
        legend_2.SetFillStyle(0)
        legend_2.SetFillColor(0)
        legend_2.SetTextFont(42)
        legend_2.SetTextSize(0.05)
        legend_2.Draw("SAME")



    # CMS style
    CMS_lumi.cmsText = "CMS"
    #CMS_lumi.extraText = "Preliminary"
    #CMS_lumi.extraText = "Supplementary"
    #CMS_lumi.cmsText = ""
    CMS_lumi.extraText = ""
    #CMS_lumi.extraText_posX = 0.07
    #CMS_lumi.extraText = "Private"
    CMS_lumi.cmsTextSize = 0.95
    CMS_lumi.CMSText_posX = -0.03
    CMS_lumi.outOfFrame = True

    CMS_lumi.lumiText_posX = -0.010
    CMS_lumi.CMS_lumi(canv,5,0,plt_cfg.year)

    canv.cd()
    #lower_pad = TPad("lowerpad_"+var_name, "lowerpad_"+var_name, 0, 0.01, 1,0.22)
    lower_pad = TPad("lowerpad_"+var_name, "lowerpad_"+var_name, 0, 0, 1,0.35)
    lower_pad.SetTopMargin(0.00001)
    lower_pad.SetBottomMargin(0.36)
    lower_pad.SetLeftMargin(United_pad_LeftMargin)
    lower_pad.SetRightMargin(United_pad_RightMargin)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    lower_pad.SetTickx(1)
    lower_pad.SetTicky(1)

    if var_name.split("_")[-1] in plt_cfg.ana_cfg.sig_names:
        ratio_plot.GetXaxis().SetTitle('BDT Score')
    else:
        ratio_plot.GetXaxis().SetTitle(plt_cfg.var_title_map[var_name])
    ratio_plot.GetXaxis().SetTitleOffset(1.0)
    ratio_plot.Draw("APZ SAME")

    # Draw_unc(total_unc[1], kGray+10)
    # Draw_unc(stat_err_norm, kRed-10)
    # Draw_unc(total_unc[1], TColor.GetColor("#324376")) # Navy Blue
    Draw_unc(stat_err_norm, TColor.GetColor("#F76C5E")) # orange

    ratio_plot.Draw("SAMEPZ")

    
    

def SaveCanvPic(canv, save_dir, save_name):
    canv.cd()
    canv.SaveAs(save_dir + '/' + save_name + '.pdf')
    canv.SaveAs(save_dir + '/' + save_name + '.png')
    # print(f"[SaveCanvPic] Writing {save_name}")
    # canv.SaveAs(save_dir + '/' + save_name + '.eps')

    canv.Close()
