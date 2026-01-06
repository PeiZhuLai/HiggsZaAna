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

def ScaleBkgToData(histos, ana_cfg, syst_histos=None):
    """
    將背景總 yield scale 到 Data，並可選擇同步縮放系統變動直方圖。
    
    Parameters
    ----------
    histos : dict[str, TH1]
        單一變數的 (sample -> TH1) 直方圖。
    ana_cfg : Analyzer_Config
    syst_histos : dict[str, dict[str, TH1]] | None
        結構為 syst_histos[sample][sys]，若提供則一併套用相同 scale。
    
    Returns
    -------
    float : 使用的縮放因子；若未縮放則為 1.0
    """
    data_hist = histos.get('Data')
    if not data_hist:
        print("[ScaleBkgToData] Warning: No 'Data' histogram found. Skipping scaling.")
        return 1.0

    data_integral = data_hist.Integral()
    if data_integral <= 0:
        print("[ScaleBkgToData] Warning: Data integral <= 0. Skip.")
        return 1.0

    bkg_integral = 0.0
    for sample in ana_cfg.bkg_names:
        h = histos.get(sample)
        if h:
            bkg_integral += h.Integral()

    if bkg_integral <= 0:
        print("[ScaleBkgToData] Warning: Background integral <= 0. Skip.")
        return 1.0

    scale_factor = data_integral / bkg_integral
    if abs(scale_factor - 1.0) < 1e-12:
        print("[ScaleBkgToData] Scaling factor ~1. No action.")
        return 1.0

    print(f"[ScaleBkgToData] Scaling factor: {scale_factor:.4f} "
          f"(Data: {data_integral:.2f}  Bkg: {bkg_integral:.2f})")

    for sample in ana_cfg.bkg_names:
        if sample in histos:
            histos[sample].Scale(scale_factor)
            # 同步縮放 syst
            if syst_histos and sample in syst_histos:
                for sys_name, hsys in syst_histos[sample].items():
                    hsys.Scale(scale_factor)

    return scale_factor

def LoadNtuples_split_by_mass(ana_cfg):
    """
    为每个 mA in ana_cfg.sig_names 构建一组链：
      ntuples_by_mass[mA][sample] = TChain(...)
    规则（文件系统目录）：
      - 信号 Mx ： base/mA_Mx/<year>.root （年表用 ana_cfg.years_sig）
      - DYJetsToLL： base/DYJetsToLL/mA_Mx/<year>.root （年表用 ana_cfg.years_dyll）
      - DYGto2LG ： base/<子样本>/mA_Mx/<year>.root （2022/2023 各自子样本与年表）
      - Data     ： 先试 base/Data/mA_Mx/<year>.root；若不存在回退 base/Data/<year>.root
    """
    base = ana_cfg.sample_loc
    ntuples_by_mass = {}  # dict[mA][sample] -> TChain

    def _add_if_exists(ch, path):
        if os.path.exists(path):
            ch.Add(path)
            return 1
        else:
            print(f"[LoadNtuples@{ch.GetName()}][MISS] {path}")
            return 0

    for mA in ana_cfg.sig_names:
        mdir = f"mA_{mA}"
        group = {}

        # --- 信号（Mx） ---
        ch_sig = ROOT.TChain("test", f"chain_sig_{mA}")
        added = 0
        for y in getattr(ana_cfg, "years_sig", ["2022preEE"]):
            p = os.path.join(base, f"mA_{mA}", f"{y}.root")
            added += _add_if_exists(ch_sig, p)
        if added == 0:
            print(f"[WARN] signal {mA} empty.")
        group[mA] = ch_sig  # 用样本名=质量名保存

        # --- DYJetsToLL（背景） ---
        ch_dy = ROOT.TChain("inclusive", f"chain_DY_{mA}")
        added = 0
        for y in getattr(ana_cfg, "years_dyll", []):
            p = os.path.join(base, "DYJetsToLL", mdir, f"{y}.root")
            added += _add_if_exists(ch_dy, p)
        if added == 0:
            print(f"[WARN] DYJetsToLL@{mA} empty.")
        group["DYJetsToLL"] = ch_dy

        # --- DYGto2LG（背景，分 2022/2023 子样本） ---
        ch_dyg = ROOT.TChain("inclusive", f"chain_DYG_{mA}")
        added = 0
        for subs in getattr(ana_cfg, "bkg_2022", []):
            for y in getattr(ana_cfg, "years_22", []):
                p = os.path.join(base, subs, mdir, f"{y}.root")
                added += _add_if_exists(ch_dyg, p)
        for subs in getattr(ana_cfg, "bkg_2023", []):
            for y in getattr(ana_cfg, "years_23", []):
                p = os.path.join(base, subs, mdir, f"{y}.root")
                added += _add_if_exists(ch_dyg, p)
        if added == 0:
            print(f"[WARN] DYGto2LG@{mA} empty.")
        group["DYGto2LG"] = ch_dyg

        # --- Data ---
        ch_data = ROOT.TChain("inclusive", f"chain_Data_{mA}")
        added = 0
        for y in getattr(ana_cfg, "years_dyll", []):
            p1 = os.path.join(base, "Data", mdir, f"{y}.root")
            added += _add_if_exists(ch_data, p1)
            if added == 0:
                p2 = os.path.join(base, "Data", f"{y}.root")
                added += _add_if_exists(ch_data, p2)
        if added == 0:
            print(f"[WARN] Data@{mA} empty.")
        group["Data"] = ch_data

        ntuples_by_mass[mA] = group
        print(f"[LoadNtuples_split_by_mass] built group for {mA}: "
              f"sig={group[mA].GetEntries()}  "
              f"DY={group['DYJetsToLL'].GetEntries()}  "
              f"DYG={group['DYGto2LG'].GetEntries()}  "
              f"Data={group['Data'].GetEntries()}")
    return ntuples_by_mass

# ==== New helpers for run3 multi-year chaining ====
def _add_file_if_exists(chain, path):
    if os.path.exists(path):
        chain.Add(path)
        return True
    else:
        print(f"[Plot_Helper][MISS] {path}")
        return False

def _run3_build_chain(sample, ana_cfg):
    """
    Build a TChain for run3 according to updated rules:
      - 全部樣本（Data / DYJetsToLL / DYGto2LG / signal）都直接讀 era ROOT：
          base/<sample-or-subSample>/<year>.root
        不再使用 /mA_M*/ 子資料夾
      - MVA score 由 branch 決定（在 plot script 端讀 MVA_Score_mA_<mass>）
    """
    tree_name = "test" if (sample in getattr(ana_cfg, "sig_names", [])) else "inclusive"
    ch = TChain(tree_name, f"chain_{sample}")
    base = ana_cfg.sample_loc
    added = 0

    if sample in ana_cfg.sig_names:
        # signal: base/mA_<signalMass>/<year>.root (沿用你目前的信號資料夾結構)
        # 若你的 signal 其實是 base/<signalSample>/<year>.root，請在這裡把 "mA_" 拿掉即可
        for y in getattr(ana_cfg, "years_sig", ["2022preEE"]):
            path = os.path.join(base, f"mA_{sample}", f"{y}.root")
            if _add_file_if_exists(ch, path):
                added += 1

    elif sample == "DYJetsToLL":
        for y in getattr(ana_cfg, "years_dyll", []):
            path = os.path.join(base, "DYJetsToLL", f"{y}.root")
            if _add_file_if_exists(ch, path):
                added += 1

    elif sample == "DYGto2LG":
        # 2022 sub-samples
        for subs in getattr(ana_cfg, "bkg_2022", []):
            for y in getattr(ana_cfg, "years_22", []):
                path = os.path.join(base, subs, f"{y}.root")
                if _add_file_if_exists(ch, path):
                    added += 1
        # 2023 sub-samples
        for subs in getattr(ana_cfg, "bkg_2023", []):
            for y in getattr(ana_cfg, "years_23", []):
                path = os.path.join(base, subs, f"{y}.root")
                if _add_file_if_exists(ch, path):
                    added += 1

    elif sample == "Data":
        for y in getattr(ana_cfg, "years_dyll", []):
            path = os.path.join(base, "Data", f"{y}.root")
            if _add_file_if_exists(ch, path):
                added += 1

    else:
        # Fallback
        path = os.path.join(base, sample, "run3.root")
        if _add_file_if_exists(ch, path):
            added += 1

    if added == 0:
        print(f"[Plot_Helper][WARN] sample={sample} has no files added (chain empty).")
    else:
        print(f"[Plot_Helper] sample={sample} added files: {added}")
    return ch

def LoadNtuples(ana_cfg):
    ntuples = {}
    # run3 special aggregation
    if ana_cfg.year == 'run3':
        for sample in ana_cfg.samp_names:
            # Decide by sample type
            ntuples[sample] = _run3_build_chain(sample, ana_cfg)
        return ntuples

    # ===== Original behavior for non-run3 (保持舊邏輯) =====
    for sample in ana_cfg.samp_names:
        if "M" in sample: 
            ntuples[sample] = TChain("test","chain_" + sample)
            ntuples[sample].Add(ana_cfg.sample_loc + '/mA_%s/run3.root' %sample)
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
    tex.DrawLatexNDC(0.9, 0.91, f"{float(lumi):.2f} fb^{{-1}} (13.6 TeV)")
    return tex

def MakeCMSDASLabel():
    #tex = TLatex()
    #tex.SetTextSize(0.03)
    #tex.DrawLatexNDC(0.12, 0.85, '#scale[1.5]{CMSDAS} H To Z + mA')
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

    ratio_plot.GetYaxis().SetNdivisions(505)
    ratio_plot.GetYaxis().SetTitle("Data / SM")
    ratio_plot.GetYaxis().CenterTitle(True)
    ratio_plot.GetYaxis().SetLabelSize(0.14)
    ratio_plot.GetYaxis().SetTitleSize(0.14)
    ratio_plot.GetYaxis().SetTitleOffset(0.535)

    ratio_plot.GetXaxis().SetLimits( h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax() )
    ratio_plot.GetXaxis().SetTitle(var_name)
    ratio_plot.GetXaxis().SetLabelSize(0.14)
    ratio_plot.GetXaxis().SetTitleSize(0.16)
    ratio_plot.GetXaxis().SetTitleOffset(0.65)
    ratio_plot.GetXaxis().SetTickLength(0.08)

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

# ==== New helper: build systematics-only band (total^2 - stat^2) ====
def BuildSysOnlyBand(total_graph, stat_graph, name_suffix="_sysOnly"):
    """
    從 (stat⊕sys) 與 stat band 推得純系統 band：
      err_sys = sqrt( err_total^2 - err_stat^2 )
    PyROOT 相容：不再使用 ROOT.Double，而直接用 GetX()/GetY() 緩衝陣列。
    """
    if not total_graph or not stat_graph:
        print("[BuildSysOnlyBand] Missing graph(s).")
        return None
    n_tot  = total_graph.GetN()
    n_stat = stat_graph.GetN()
    if n_tot != n_stat:
        print(f"[BuildSysOnlyBand] Point size mismatch (total={n_tot}, stat={n_stat}). Will use min.")
    n = min(n_tot, n_stat)
    if n == 0:
        print("[BuildSysOnlyBand] No points.")
        return None

    xs_buf_tot = total_graph.GetX()
    ys_buf_tot = total_graph.GetY()
    xs_buf_stat = stat_graph.GetX()
    # ys_buf_stat = stat_graph.GetY()  # y 其實不需獨立

    xs   = []
    ys   = []
    exl  = []
    exh  = []
    eyl  = []
    eyh  = []

    for i in range(n):
        # 中心值 (採用 total, 兩者應一致)
        xs.append(float(xs_buf_tot[i]))
        ys.append(float(ys_buf_tot[i]))

        # X error
        exl.append(total_graph.GetErrorXlow(i))
        exh.append(total_graph.GetErrorXhigh(i))

        # Y errors
        t_up = total_graph.GetErrorYhigh(i)
        t_dn = total_graph.GetErrorYlow(i)
        s_up = stat_graph.GetErrorYhigh(i)
        s_dn = stat_graph.GetErrorYlow(i)

        add_up = math.sqrt(max(0.0, t_up * t_up - s_up * s_up))
        add_dn = math.sqrt(max(0.0, t_dn * t_dn - s_dn * s_dn))

        eyl.append(add_dn)
        eyh.append(add_up)

    g_sys = ROOT.TGraphAsymmErrors(
        n,
        np.array(xs), np.array(ys),
        np.array(exl), np.array(exh),
        np.array(eyl), np.array(eyh)
    )
    g_sys.SetName(total_graph.GetName() + name_suffix)
    return g_sys

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

def Draw_unc(graph, color, alpha=0.1, draw_option="2 SAME", on_top=True, fill_style=None):
    """
    繪製不確定度 band。
    修正：原實作在指定 fill_style 時不會呼叫 SetFillColorAlpha，導致 alpha 被忽略。
    新邏輯：
      1. 若物件支援 SetFillColorAlpha，先嘗試套用透明度。
      2. 若指定 fill_style：
         - 若為實心 (1001) 仍保留透明度。
         - 若為陰影/斜線樣式 (>=3000) ROOT 不支援透明度，回退為不透明。
      3. 若透明度設定失敗，回退到舊式 SetFillColor + 預設 / 指定 fill_style。
    """
    if not graph:
        print("[Draw_unc] Warning: graph is None")
        return

    graph.SetMarkerStyle(0)
    graph.SetLineColor(color)
    graph.SetLineWidth(1)

    alpha_applied = False

    # 嘗試透明度 (僅在有 SetFillColorAlpha 時)
    if hasattr(graph, "SetFillColorAlpha"):
        try:
            # 針對可見度：alpha 太小(例如 0.05) 會幾乎看不到，由使用者自行調整
            graph.SetFillColorAlpha(color, alpha)
            alpha_applied = True
        except Exception:
            alpha_applied = False

    # 如果指定了 fill_style
    if fill_style is not None:
        # 針對陰影/斜線樣式透明度其實沒作用，必要時覆蓋顏色
        if not alpha_applied or fill_style >= 3000:
            graph.SetFillColor(color)
        graph.SetFillStyle(fill_style)
    else:
        # 未指定 fill_style，若透明度失敗則回退預設樣式
        if not alpha_applied:
            graph.SetFillColor(color)
            graph.SetFillStyle(3004)  # 舊預設
        else:
            graph.SetFillStyle(1001)  # 透明實心

    opt = draw_option
    if "SAME" not in opt:
        opt += " SAME"
    graph.Draw(opt)
    if on_top:
        gPad.Modified()
        gPad.Update()

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
        
        if histos['Data'].GetMaximum() > stacks['all'].GetMaximum():
            h_max = histos['Data'].GetMaximum()
        else:
            h_max = stacks['all'].GetMaximum()
        if h_max < stacks['sig'].GetMaximum():
            h_max = stacks['sig'].GetMaximum()

        upper_pad.SetLogy()
        stacks['all'].SetMinimum(1e-2)
        histos['Data'].SetMinimum(1e-2)

        # stacks['all'].SetMaximum(1e10)
        # histos['Data'].SetMaximum(1e10)

        histos['Data'].SetMaximum(h_max*1.1e4)
        stacks['all'].SetMaximum(h_max*1.1e4)

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
    
    if var_name in ["H_m","ALP_m","Z_m"]:
        histos['Data'].GetYaxis().SetTitle('Events / %.2f GeV' %histos['Data'].GetBinWidth(1))
    else:
        histos['Data'].GetYaxis().SetTitle('Events / %.3f' %histos['Data'].GetBinWidth(1))
    histos['Data'].GetYaxis().SetTitleSize(0.07)
    histos['Data'].GetYaxis().SetLabelSize(0.055)
    histos['Data'].GetYaxis().SetTitleFont(42)
    histos['Data'].GetYaxis().SetTitleOffset(1.15)
    stacks['all'].Draw('HISTSAME')

    if (var_name.split("_")[-1] in plt_cfg.ana_cfg.sig_names):
        scaled_sig[var_name.split("_")[-1]].Draw('HISTSAME')
    else:
        if bdtCut:
            scaled_sig[mA].Draw('HISTSAME')
        else:
            line_styles = [1, 2, 5, 7]
            # Run2 Paper
            for sample, i in zip(["M1", "M10", "M20", "M30"], range(4)):
            # for sample, i in zip(["M5", "M15", "M30"], range(3)):
                scaled_sig[sample].SetLineStyle(line_styles[i])
                scaled_sig[sample].SetLineWidth(4)
                scaled_sig[sample].Draw('HISTSAME')

    histos['Data'].SetMarkerStyle(20)
    histos['Data'].SetLineWidth(2)
    histos['Data'].GetXaxis().SetTickLength(0.04)
    
    ### Draw the uncertainties
    global stat_err, stat_err_norm
    stat_err,  stat_err_norm = Get_StatUnc(stacks['bkg'].GetStack().Last())

    total_abs  = total_unc[0]   # (stat⊕syst)
    total_norm = total_unc[1]

    # 先畫總不確定度 (外層) 再畫統計 (內層)
    Draw_unc(total_abs,  TColor.GetColor("#A0A0A0"), alpha=0.3, fill_style=1001)   # Total (stat⊕syst)
    Draw_unc(stat_err,   TColor.GetColor("#404040"), alpha=0.90, fill_style=3354)   # Stat only

    histos['Data'].Draw('SAMEPE')
    histos['Data'].Draw("AXIS SAME")

    if var_name.split("_")[-1] in plt_cfg.ana_cfg.sig_names:
        # -------------------------------------------------------------------------------------------
        legend_1 = TLegend(0.66, 0.59, 0.97, 0.86)
        ROOT.SetOwnership(legend_1, False)
        legend_1.AddEntry(histos["Data"], "Data", "PE")
        bkg_labels = {"DYGto2LG": r"Z + \gamma",
                    "DYJetsToLL": r"Z + jets"}    
        
        for sample_bkg in plt_cfg.ana_cfg.bkg_names:            
            legend_1.AddEntry(histos[sample_bkg], bkg_labels.get(sample_bkg, sample_bkg), "f")

        legend_1.AddEntry(total_abs,"Total Unc.","f")
        legend_1.AddEntry(stat_err,"Stat. Unc.","f")
        # -------------------------------------------------------------------------------------------
        legend_2 = TLegend(0.40, 0.80, 0.66, 0.86)
        ROOT.SetOwnership(legend_2, False)
        legend_2.AddEntry(scaled_sig[var_name.split("_")[-1]], r"m_{a} = %s GeV" % (var_name.split("_")[-1].lstrip("M")), "l" )
        # -------------------------------------------------------------------------------------------
        legend_1.SetBorderSize(0)
        legend_1.SetFillStyle(0)
        legend_1.SetFillColor(0)
        legend_1.SetTextFont(42)
        legend_1.SetTextSize(0.05)
        legend_1.Draw("SAME")
        # -------------------------------------------------------------------------------------------
        legend_2.SetBorderSize(0)
        legend_2.SetFillStyle(0)
        legend_2.SetFillColor(0)
        legend_2.SetTextFont(42)
        legend_2.SetTextSize(0.05)
        legend_2.Draw("SAME")


    else:
        legend_1 = TLegend(0.29, 0.816, 0.50, 0.87)
        ROOT.SetOwnership(legend_1, False)
        legend_1.AddEntry(histos["Data"], "Data", "PE")
        # -------------------------------------------------------------------------------------------
        legend_2 = TLegend(0.44, 0.654, 0.73, 0.87)
        ROOT.SetOwnership(legend_2, False)
        bkg_labels = {"DYGto2LG": r"Z + \gamma",
                    "DYJetsToLL": r"Z + jets"}    
            
        for sample_bkg in plt_cfg.ana_cfg.bkg_names:            
            legend_2.AddEntry(histos[sample_bkg], bkg_labels.get(sample_bkg, sample_bkg), "f")
            
        legend_2.AddEntry(total_abs,"Total Unc.","f")
        legend_2.AddEntry(stat_err,"Stat. Unc.","f")
        # -------------------------------------------------------------------------------------------
        legend_3 = TLegend(0.67, 0.654, 0.97, 0.87)
        ROOT.SetOwnership(legend_3, False)
        if bdtCut:
            legend_3.AddEntry(scaled_sig[mA], r"m_{a} = %s GeV" % (mA.lstrip("M")), "l")
        else:
            for s in ["M1","M10","M20","M30"]:
                legend_3.AddEntry(scaled_sig[s], r"m_{a} = %s GeV" % (s.lstrip("M")), "l")
        # -------------------------------------------------------------------------------------------
        legend_1.SetBorderSize(0)
        legend_1.SetFillStyle(0)
        legend_1.SetFillColor(0)
        legend_1.SetTextFont(42)
        legend_1.SetTextSize(0.045)
        legend_1.Draw("SAME")
        # -------------------------------------------------------------------------------------------
        legend_2.SetBorderSize(0)
        legend_2.SetFillStyle(0)
        legend_2.SetFillColor(0)
        legend_2.SetTextFont(42)
        legend_2.SetTextSize(0.045)
        legend_2.Draw("SAME")
        # -------------------------------------------------------------------------------------------
        legend_3.SetBorderSize(0)
        legend_3.SetFillStyle(0)
        legend_3.SetFillColor(0)
        legend_3.SetTextFont(42)
        legend_3.SetTextSize(0.045)
        legend_3.Draw("SAME")

    # CMS style
    CMS_lumi.cmsText = "CMS"
    CMS_lumi.extraText = "Preliminary"
    #CMS_lumi.extraText = "Supplementary"
    #CMS_lumi.cmsText = ""
    # CMS_lumi.extraText = ""
    #CMS_lumi.extraText_posX = 0.07
    #CMS_lumi.extraText = "Private"
    CMS_lumi.cmsTextSize = 0.95
    CMS_lumi.CMSText_posX = -0.03
    CMS_lumi.outOfFrame = True

    CMS_lumi.lumiText_posX = -0.000 # 0.01 
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

    Draw_unc(total_norm, TColor.GetColor("#A0A0A0"), alpha=0.3, fill_style=1001)
    Draw_unc(stat_err_norm, TColor.GetColor("#404040"), alpha=0.90, fill_style=3354)
    ratio_plot.Draw("SAMEPZ")

    
    

def SaveCanvPic(canv, save_dir, save_name):
    canv.cd()
    canv.SaveAs(save_dir + '/' + save_name + '.pdf')
    canv.SaveAs(save_dir + '/' + save_name + '.png')
    # print(f"[SaveCanvPic] Writing {save_name}")
    # canv.SaveAs(save_dir + '/' + save_name + '.eps')

    canv.Close()
