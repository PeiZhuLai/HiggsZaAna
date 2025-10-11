####################################################
####################################################

import os
import sys
import numpy as np
import gc

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples_split_by_mass, LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend, Total_Unc, ScaleBkgToData
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
parser.add_argument("-y", "--Year", dest="year", default="run3", help="which year's datasetes")
parser.add_argument("--region", dest="region", type=int, default=0, help="0 for full region, 1 for signal region, 2 for sideband region")
parser.add_argument('-m', '--mva', dest='mva', action='store_true', default=False, help='use mva or not')
parser.add_argument('--cut', dest='cut', action='store_true', default=False, help='apply mva')
parser.add_argument("--cutVal", dest="cutVal", type=float, default=0.0, help="mva cut value")
parser.add_argument("--mA", dest="mA", default="M5", help="ALP mass")
parser.add_argument('-b', '--blind', dest='blind', action='store_true', default=False, help='Blind signal region')
parser.add_argument('-ln', '--ln', dest='ln', action='store_true', default=False, help='log plot?')
parser.add_argument('--ele', dest='ele', action='store_true', default=False, help='electron channel?')
parser.add_argument('--mu', dest='mu', action='store_true', default=False, help='muon channel?')
# 新增：MVA 偵錯輸出控制
parser.add_argument('--mvaDebug', dest='mva_debug', action='store_true', default=False, help='print per-mA fill info')
parser.add_argument('--mvaDebugN', dest='mva_debug_n', type=int, default=15, help='max debug prints per (sample,mA)')
args = parser.parse_args()


gROOT.SetBatch(True)
tdrstyle.setTDRStyle()


# load the model from disk
# model = pickle.load(open(BDT_filename, 'rb'))

def SideBandScaleBkgToData(histos, histos_sys, analyzer_cfg, signal_low=115., signal_high=135.):
    """
    只針對 histos['H_m']：
      將背景總和在 sideband (H_m < signal_low 或 H_m > signal_high) 的 yield
      正規化到資料在同一 sideband 的 yield。
    """
    if 'H_m' not in histos:
        print('[SideBandScale] H_m not found, skip.')
        return 1.0

    # 資料 sample key 可能大小寫不一致
    data_key = None
    for k in histos['H_m'].keys():
        if k.lower() == 'data':
            data_key = k
            break
    if data_key is None:
        print('[SideBandScale] Data hist not found, skip.')
        return 1.0

    axis = histos['H_m'][data_key].GetXaxis()
    nb = axis.GetNbins()
    # 找到 signal window 所對應的 bin
    sig_low_bin = axis.FindFixBin(signal_low)
    sig_high_bin = axis.FindFixBin(signal_high)

    def sideband_integral(h):
        # 左側
        left = 0.0
        if sig_low_bin > 1:
            left = h.Integral(1, sig_low_bin - 1)
        # 右側
        right = 0.0
        if sig_high_bin < nb:
            right = h.Integral(sig_high_bin + 1, nb)
        return left + right

    data_sb = sideband_integral(histos['H_m'][data_key])

    bkg_samples = [
        s for s in analyzer_cfg.samp_names
        if (s.lower() != 'data') and (s not in analyzer_cfg.sig_names)
    ]

    bkg_sb = 0.0
    for s in bkg_samples:
        bkg_sb += sideband_integral(histos['H_m'][s])

    if bkg_sb <= 0:
        print('[SideBandScale] Background sideband yield = 0, skip scaling.')
        return 1.0

    scale_factor = data_sb / bkg_sb
    # 套用縮放到背景 (包含 systematics)
    for s in bkg_samples:
        histos['H_m'][s].Scale(scale_factor)
        if histos_sys and 'H_m' in histos_sys:
            for sys_name in analyzer_cfg.sys_names:
                histos_sys['H_m'][s][sys_name].Scale(scale_factor)

    print(f'[SideBandScale] H_m sideband Data={data_sb:.3f}  Bkg(before)={bkg_sb:.3f}  Scale={scale_factor:.4f}')
    return scale_factor


# 新增：偵測指定 mA 的 MVA 分支名稱（回傳第一個存在的分支，否則預設為 "MVA_Score"）
def _resolve_mva_branch_for_mass(chain, mass_tag):
    """
    嘗試在 TChain 上尋找對應 mA 的 MVA 分支，常見候選：
      - MVA_Score_<M>
      - MVA_score_<M>
      - BDTScore_<M>
      - BDT_<M>
    找不到則回退為 "MVA_Score"。
    """
    if not chain or not chain.GetListOfBranches():
        return "MVA_Score"
    br_list = chain.GetListOfBranches()
    candidates = [
        f"MVA_Score_{mass_tag}",
        f"MVA_score_{mass_tag}",
        f"BDTScore_{mass_tag}",
        f"BDT_{mass_tag}",
        "MVA_Score",
    ]
    for name in candidates:
        if br_list.FindObject(name):
            return name
    return "MVA_Score"


def main():

    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year, args.region, args.mva)

    # 將欲使用的 ALP mass 傳給載檔函式（mva 模式下用於路徑 ALP_{mass}）
    analyzer_cfg.mva = bool(args.mva)
    analyzer_cfg.mva_alp_mass = str(args.mA) if args.mva else "M1"

    # 在 mva + run3 下，自動按 self.sig_names 跑每個 mA；其他情況維持單一目標質量
    if args.mva and analyzer_cfg.year == 'run3':
        target_masses = list(analyzer_cfg.sig_names)
    else:
        target_masses = [analyzer_cfg.mva_alp_mass] if args.mva else []

    if not os.path.exists(analyzer_cfg.out_dir):
        os.makedirs(analyzer_cfg.out_dir)

    if not os.path.exists(analyzer_cfg.plot_output_path):
        os.makedirs(analyzer_cfg.plot_output_path)

    out_file = TFile( analyzer_cfg.out_dir + '/' + analyzer_cfg.root_output_name , "RECREATE")
    # model = pickle.load(open(analyzer_cfg.BDT_filename, 'rb'))

    
    if args.cut: 
        analyzer_cfg.sig_names = [args.mA]
        # 修正大小寫，與其他地方一致使用 'Data'
        analyzer_cfg.samp_names = analyzer_cfg.bkg_names + analyzer_cfg.sig_names + ['Data']


    analyzer_cfg.Print_Config()
    if args.mva and analyzer_cfg.year == 'run3':
        ntuples_by_mass = LoadNtuples_split_by_mass(analyzer_cfg)
    else:
        ntuples = LoadNtuples(analyzer_cfg)
    plot_cfg = PC.Plot_Config(analyzer_cfg, args.year)
    var_names = list(plot_cfg.var_title_map.keys())

    if args.mva:
        # 僅針對選定質量新增 MVA 變數
        var_mva = []
        for ALP_mass in target_masses:
            var_names.append('mvaVal_'+ALP_mass)
            var_names.append('mvaVal_larger_'+ALP_mass)
            for r in ['1sigma', '1P5sigma', '2sigma', '3sigma']:
                var_mva.append('mvaVal_'+r+'_'+ALP_mass)
                var_mva.append('mvaVal_larger_'+r+'_'+ALP_mass)
        var_names = var_names + var_mva

    # get mass region
    sigma_low, sigma_hig = getMassSigma(analyzer_cfg)

    # 為每個 (sample, mA) 解析一次 MVA 分支
    mva_branch_map = {}  # dict[sample][mass] -> branch_name
    if args.mva:
        # 需要一份樣本清單；注意：信號樣本名稱就是 'M5'/'M15'/'M30'
        def _samples_for_mass(ALP_mass):
            return analyzer_cfg.bkg_names + [ALP_mass, 'Data']

        if analyzer_cfg.year == 'run3' and 'ntuples_by_mass' in locals():
            # 逐 mA、逐 sample 解析
            for ALP_mass in target_masses:
                for s in _samples_for_mass(ALP_mass):
                    mva_branch_map.setdefault(s, {})
                    chain = ntuples_by_mass[ALP_mass][s]
                    br = _resolve_mva_branch_for_mass(chain, ALP_mass)
                    mva_branch_map[s][ALP_mass] = br
                    print(f"[MVA] sample={s:>12s} mass={ALP_mass:>3s} uses branch '{br}'")
        else:
            # 回退：舊邏輯（無分質量目錄）
            for s in analyzer_cfg.samp_names:
                chain = ntuples.get(s, None)
                mva_branch_map[s] = {}
                for ALP_mass in target_masses:
                    br = _resolve_mva_branch_for_mass(chain, ALP_mass)
                    mva_branch_map[s][ALP_mass] = br
                    print(f"[MVA] sample={s:>12s} mass={ALP_mass:>3s} uses branch '{br}'")


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
        histos['pho1IetaIeta55'][sample]    = TH1F('pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 15,  0., 0.06)
        histos['pho1ECALIso'][sample]    = TH1F('pho1ECALIso'    + '_' + sample, 'pho1ECALIso'    + '_' + sample, 10, 0., 40.)
        histos['pho1CIso'][sample]    = TH1F('pho1CIso'    + '_' + sample, 'pho1CIso'    + '_' + sample, 10, 0., 0.7)
        histos['pho1HCALIso'][sample]  = TH1F('pho1HCALIso'    + '_' + sample, 'pho1HCALIso'    + '_' + sample, 10, 0., 3.0)
        histos['pho1HOE'][sample]    = TH1F('pho1HOE'    + '_' + sample, 'pho1HOE'    + '_' + sample, 10, 0., 0.032)
        histos['pho2Pt'][sample]    = TH1F('pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 12,  8., 30.)
        histos['pho2eta'][sample]    = TH1F('pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 20,  -3., 3.)
        histos['pho2phi'][sample]    = TH1F('pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 20,  -4., 4.)
        histos['pho2R9'][sample]    = TH1F('pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 25,  0.1, 1.)
        histos['pho2IetaIeta55'][sample]    = TH1F('pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 15,  0., 0.06)
        histos['pho2ECALIso'][sample]    = TH1F('pho2ECALIso'    + '_' + sample, 'pho2ECALIso'    + '_' + sample, 10, 0., 40.)
        histos['pho2CIso'][sample]    = TH1F('pho2CIso'    + '_' + sample, 'pho2CIso'    + '_' + sample, 10, 0., 0.7)
        histos['pho2HCALIso'][sample]    = TH1F('pho2HCALIso'    + '_' + sample, 'pho2HCALIso'    + '_' + sample, 10, 0., 3.)
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
        histos['ALP_calculatedPhotonIso'][sample] = TH1F('ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 20, 0., 125.)
        histos['param'][sample] = TH1F('param' + '_' + sample, 'param' + '_' + sample, 25, -0.3, 0.6)

        if args.mva:
            # 僅為選定質量建立 MVA 相關直方圖
            for ALP_mass in target_masses:
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
    # 新增：控制每個 (sample, mA) 的偵錯輸出次數
    debug_printed = {}

    # ========= RDataFrame 版本：替换逐事件循环 =========
    import ROOT
    ROOT.ROOT.EnableImplicitMT()   # 开启多线程

    def _apply_common_filters(df, args):
        # 质量窗与区域
        df = df.Filter("H_m>-90 && H_m>95 && H_m<180")
        if args.region == 1:
            df = df.Filter("H_m>115 && H_m<135")
        elif args.region == 2:
            df = df.Filter("H_m<115 || H_m>135")
        # 通道选择
        if args.ele:
            df = df.Filter("abs(z_mumu)!=1")
        if args.mu:
            df = df.Filter("abs(z_ee)!=1")
        return df

    # 统一把系统权重列定义出来（Data 不使用这些列也无妨）
    def _define_sys_weights(df, sample):
        """
        安全地定义中心权重 w 以及系统变体。
        - Data: 基本设为 w = 1.0（或若有 weight 分支则用它）
        - MC: 若某个 up/down/central 缺失 => 该系统权重退回 w
        """
        # 小工具：检查列是否存在
        def has_col(dframe, name: str) -> bool:
            return dframe.GetColumnType(name) != ""

        # 先定义中心权重 w
        if sample == "Data":
            # Data 通常不需要 event weight；如果你的 Data 里有 weight 分支，改成 "weight"
            if has_col(df, "weight"):
                df = df.Define("w", "weight")
            else:
                df = df.Define("w", "1.0")
        else:
            # MC：若没有 weight 分支，退回 1.0
            if has_col(df, "weight"):
                df = df.Define("w", "weight")
            else:
                df = df.Define("w", "1.0")

        # 辅助：有就按比值定义；缺任一列就退回 w
        def add_sys(df, up, down, cen, wname_up, wname_down):
            up_ok   = has_col(df, up)
            dn_ok   = has_col(df, down)
            cen_ok  = has_col(df, cen)
            if up_ok and cen_ok:
                df = df.Define(wname_up,   f"w * ({up} / {cen})")
            else:
                df = df.Define(wname_up,   "w")
            if dn_ok and cen_ok:
                df = df.Define(wname_down, f"w * ({down} / {cen})")
            else:
                df = df.Define(wname_down, "w")
            return df

        # 逐系统安全定义（名字来自你原代码）
        df = add_sys(df, "weight_hlt_sf_up", "weight_hlt_sf_down", "weight_hlt_sf_central",
                    "w_hlt_up", "w_hlt_down")

        df = add_sys(df, "weight_pu_reweight_sf_up", "weight_pu_reweight_sf_down",
                    "weight_pu_reweight_sf_central",
                    "w_pu_up", "w_pu_down")

        df = add_sys(df, "weight_electron_wplid_sf_SelectedElectron_up",
                        "weight_electron_wplid_sf_SelectedElectron_down",
                        "weight_electron_wplid_sf_SelectedElectron_central",
                    "w_el_id_up", "w_el_id_down")

        df = add_sys(df, "weight_electron_iso_sf_SelectedElectron_up",
                        "weight_electron_iso_sf_SelectedElectron_down",
                        "weight_electron_iso_sf_SelectedElectron_central",
                    "w_el_iso_up", "w_el_iso_down")

        df = add_sys(df, "weight_electron_reco_sf_SelectedElectron_up",
                        "weight_electron_reco_sf_SelectedElectron_down",
                        "weight_electron_reco_sf_SelectedElectron_central",
                    "w_el_reco_up", "w_el_reco_down")

        df = add_sys(df, "weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_up",
                        "weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_down",
                        "weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_central",
                    "w_el_idnm_up", "w_el_idnm_down")

        df = add_sys(df, "weight_muon_looseid_sf_SelectedMuon_up",
                        "weight_muon_looseid_sf_SelectedMuon_down",
                        "weight_muon_looseid_sf_SelectedMuon_central",
                    "w_mu_id_up", "w_mu_id_down")

        df = add_sys(df, "weight_muon_iso_sf_SelectedMuon_up",
                        "weight_muon_iso_sf_SelectedMuon_down",
                        "weight_muon_iso_sf_SelectedMuon_central",
                    "w_mu_iso_up", "w_mu_iso_down")

        df = add_sys(df, "weight_muon_reco_sf_SelectedMuon_up",
                        "weight_muon_reco_sf_SelectedMuon_down",
                        "weight_muon_reco_sf_SelectedMuon_central",
                    "w_mu_reco_up", "w_mu_reco_down")

        df = add_sys(df, "weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_up",
                        "weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_down",
                        "weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_central",
                    "w_mu_idnm_up", "w_mu_idnm_down")

        return df

    # 盲区处理：Data 在 SR 直接过滤掉
    def _apply_blind(df, sample, args):
        if sample == "Data" and args.blind:
            return df.Filter("H_m<=115. || H_m>=135.")
        return df

    # 统一直方图规格（与你原本一致）
    _hist_specs = {
    'pho1Pt':        (25,  8., 50., "ALP_lead_photon_pt"),
    'pho1eta':       (20, -3.,  3., "ALP_lead_photon_eta"),
    'pho1phi':       (20, -4.,  4., "ALP_lead_photon_phi"),
    'pho1R9':        (25,  0.1, 1., "ALP_lead_photon_r9"),
    'pho1IetaIeta55':(15,  0., .06, "ALP_lead_photon_sieie"),
    'pho1ECALIso':   (10,  0., 40., "ALP_lead_photon_ecalPFClusterIso"),
    'pho1CIso':      (10,  0., 0.7, "ALP_lead_photon_chiso"),
    'pho1HCALIso':   (10,  0., 3.0, "ALP_lead_photon_hcalPFClusterIso"),
    'pho1HOE':       (10,  0., 0.032,"ALP_lead_photon_hoe_PUcorr"),
    'pho2Pt':        (12,  8., 30., "ALP_sublead_photon_pt"),
    'pho2eta':       (20, -3.,  3., "ALP_sublead_photon_eta"),
    'pho2phi':       (20, -4.,  4., "ALP_sublead_photon_phi"),
    'pho2R9':        (25,  0.1, 1., "ALP_sublead_photon_r9"),
    'pho2IetaIeta55':(15,  0., .06, "ALP_sublead_photon_sieie"),
    'pho2ECALIso':   (10,  0., 40., "ALP_sublead_photon_ecalPFClusterIso"),
    'pho2CIso':      (10,  0., 0.7, "ALP_sublead_photon_chiso"),
    'pho2HCALIso':   (10,  0., 3.,  "ALP_sublead_photon_hcalPFClusterIso"),
    'pho2HOE':       (10,  0., 0.032,"ALP_sublead_photon_hoe_PUcorr"),
    'Z_m':           (20, 50., 130., "Z_mass"),
    'H_m':           (25, 95., 180., "H_m"),
    'H_pt':          (20,  0., 160., "H_pt"),
    'ALP_m':         (20,  0., 40., "ALP_m"),
    'var_dR_g1g2':   (20,  0.,  5.,  "var_dR_g1g2"),
    'var_PtaOverMa': (20,  0., 100., "var_PtaOverMa"),
    'var_dR_Za':     (20,  0.,  7.,  "var_dR_Za"),
    'var_dR_g1Z':    (20,  0.,  7.,  "var_dR_g1Z"),
    'var_PtaOverMh': (25,  0., 0.75, "var_PtaOverMh"),
    'var_Pta':       (20,  0.,  60., "var_Pta"),
    'var_MhMa':      (20, 100., 200., "var_MhMa"),
    'var_MhMZ':      (20, 145., 310., "var_MhMZ"),
    'ALP_calculatedPhotonIso': (20, 0., 125., "ALP_calculatedPhotonIso"),
    }

    # param 变量（依你原逻辑，用 mass_list 或随机挑一个）
    ROOT.gInterpreter.Declare(r"""
    static const double MASS_LIST[] = {1,2,3,4,5,6,7,8,9,10,15,20,25,30};
    static const int MASS_LIST_N = sizeof(MASS_LIST)/sizeof(double);
    inline double pick_mass(ULong64_t entry) {
    return MASS_LIST[ entry % MASS_LIST_N ];
    }
    """)

    def _define_param(df, sample):
        # 信号样本：用它自己的质量（如 M5 -> 5.0）
        if sample in ["M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M15","M20","M25","M30"]:
            mass_val = float(sample[1:])   # 去掉前面的 'M'
            return df.Define("param", f"(ALP_m - {mass_val})/H_m")
        # 背景样本：用 pick_mass(rdfentry_) 从固定表中按事件号取一个
        else:
            return df.Define("param", "(ALP_m - pick_mass(rdfentry_))/H_m")

    # 把一个 DataFrame 的变量批量投到 histos / histos_sys
    def _fill_hists_from_df(df, sample, var_names, sys_names, histos, histos_sys,
                            mva_fillers=None):
        """
        mva_fillers: 可选，形如 [{'name':'mvaVal_M5', 'expr':'MVA_Score_M5', 'bins':(240,-0.1,1.1), 'mask':None}, ...]
        """
        # 中心权重 + 所有系统
        weight_cols = [("w", None)]
        sys_map = {
            "weight_hlt_sf_up": "w_hlt_up",
            "weight_hlt_sf_down": "w_hlt_down",
            "weight_pu_reweight_sf_up": "w_pu_up",
            "weight_pu_reweight_sf_down": "w_pu_down",
            "weight_electron_wplid_sf_SelectedElectron_up": "w_el_id_up",
            "weight_electron_wplid_sf_SelectedElectron_down": "w_el_id_down",
            "weight_electron_iso_sf_SelectedElectron_up": "w_el_iso_up",
            "weight_electron_iso_sf_SelectedElectron_down": "w_el_iso_down",
            "weight_electron_reco_sf_SelectedElectron_up": "w_el_reco_up",
            "weight_electron_reco_sf_SelectedElectron_down": "w_el_reco_down",
            "weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_up": "w_el_idnm_up",
            "weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_down": "w_el_idnm_down",
            "weight_muon_looseid_sf_SelectedMuon_up": "w_mu_id_up",
            "weight_muon_looseid_sf_SelectedMuon_down": "w_mu_id_down",
            "weight_muon_iso_sf_SelectedMuon_up": "w_mu_iso_up",
            "weight_muon_iso_sf_SelectedMuon_down": "w_mu_iso_down",
            "weight_muon_reco_sf_SelectedMuon_up": "w_mu_reco_up",
            "weight_muon_reco_sf_SelectedMuon_down": "w_mu_reco_down",
            "weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_up": "w_mu_idnm_up",
            "weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_down": "w_mu_idnm_down",
        }
        for sname in sys_names:
            if sname in sys_map:
                weight_cols.append((sys_map[sname], sname))

        # 非 MVA 变量
        for v in var_names:
            if v.startswith("mvaVal"):  # MVA 变量单独处理
                continue
            if v not in _hist_specs:
                continue
            nb, lo, hi, expr = _hist_specs[v]
            # 中心
            h = df.Histo1D((f"{v}_{sample}", f"{v}_{sample}", nb, lo, hi), expr, "w")
            # 取回 TH1 并写入已有容器
            histos[v][sample] = h.GetValue().Clone()
            # 系统
            for wcol, sname in weight_cols[1:]:
                hsys = df.Histo1D((f"{v}_{sample}_{sname}", f"{v}_{sample}_{sname}", nb, lo, hi), expr, wcol)
                histos_sys[v][sample][sname] = hsys.GetValue().Clone()

        # MVA 专属变量（如有）
        if mva_fillers:
            for item in mva_fillers:
                name  = item["name"]
                expr  = item["expr"]
                nb, lo, hi = item["bins"]
                mask  = item.get("mask")  # 例如 1σ/1.5σ/2σ/3σ 用 H_m 区间
                df2 = df if mask is None else df.Filter(mask)
                h = df2.Histo1D((f"{name}_{sample}", f"{name}_{sample}", nb, lo, hi), expr, "w")
                histos[name][sample] = h.GetValue().Clone()
                for wcol, sname in weight_cols[1:]:
                    hsys = df2.Histo1D((f"{name}_{sample}_{sname}", f"{name}_{sample}_{sname}", nb, lo, hi), expr, wcol)
                    histos_sys[name][sample][sname] = hsys.GetValue().Clone()

    # ========= 主逻辑：两种输入结构分别处理 =========
    if args.mva and analyzer_cfg.year == 'run3' and 'ntuples_by_mass' in locals():
        # 按 mA 分目录
        for ALP_mass in target_masses:
            samples_this_mass = analyzer_cfg.bkg_names + [ALP_mass, 'Data']
            # 预构 MVA 名称与 σ 区间掩码
            mva_bins = (240, -0.1, 1.1)
            mva_bins_large = (10, 0., 1.0)
            # 每个 ALP_mass 对应的 MVA 分支名
            # 已在 mva_branch_map 里解析好
            sigmaL = 125. + sigma_low[ALP_mass]
            sigmaH = 125. + sigma_hig[ALP_mass]

            for sample in samples_this_mass:
                chain = ntuples_by_mass[ALP_mass][sample]
                df = ROOT.RDataFrame(chain)
                df = _apply_common_filters(df, args)
                df = _define_sys_weights(df, sample)
                df = _apply_blind(df, sample, args)
                df = _define_param(df, sample)

                # 定义 MVA 列
                if args.mva:
                    br = mva_branch_map.get(sample, {}).get(ALP_mass, "MVA_Score")
                    if df.GetColumnType(br) == "":  # 分支不存在时，给常数 -1
                        df = df.Define(f"__mva_{ALP_mass}", "-1.0")
                        mva_col = f"__mva_{ALP_mass}"
                    else:
                        mva_col = br

                # MVA 直方图任务清单
                mva_fillers = []
                if args.mva:
                    base = f"mvaVal_{ALP_mass}"
                    baseL = f"mvaVal_larger_{ALP_mass}"
                    mva_fillers += [
                        {"name": base,  "expr": mva_col, "bins": mva_bins},
                        {"name": baseL, "expr": mva_col, "bins": mva_bins_large},
                        {"name": f"mvaVal_1sigma_{ALP_mass}",        "expr": mva_col, "bins": mva_bins,       "mask": f"(H_m>{sigmaL}) && (H_m<{sigmaH})"},
                        {"name": f"mvaVal_larger_1sigma_{ALP_mass}", "expr": mva_col, "bins": mva_bins_large, "mask": f"(H_m>{sigmaL}) && (H_m<{sigmaH})"},
                        {"name": f"mvaVal_1P5sigma_{ALP_mass}",        "expr": mva_col, "bins": mva_bins,       "mask": f"(H_m>{125.+1.5*sigma_low[ALP_mass]}) && (H_m<{125.+1.5*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_larger_1P5sigma_{ALP_mass}", "expr": mva_col, "bins": mva_bins_large, "mask": f"(H_m>{125.+1.5*sigma_low[ALP_mass]}) && (H_m<{125.+1.5*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_2sigma_{ALP_mass}",        "expr": mva_col, "bins": mva_bins,       "mask": f"(H_m>{125.+2.0*sigma_low[ALP_mass]}) && (H_m<{125.+2.0*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_larger_2sigma_{ALP_mass}", "expr": mva_col, "bins": mva_bins_large, "mask": f"(H_m>{125.+2.0*sigma_low[ALP_mass]}) && (H_m<{125.+2.0*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_3sigma_{ALP_mass}",        "expr": mva_col, "bins": mva_bins,       "mask": f"(H_m>{125.+3.0*sigma_low[ALP_mass]}) && (H_m<{125.+3.0*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_larger_3sigma_{ALP_mass}", "expr": mva_col, "bins": mva_bins_large, "mask": f"(H_m>{125.+3.0*sigma_low[ALP_mass]}) && (H_m<{125.+3.0*sigma_hig[ALP_mass]})"},
                    ]

                _fill_hists_from_df(df, sample, var_names, analyzer_cfg.sys_names,
                                    histos, histos_sys, mva_fillers=mva_fillers)

    else:
        # 传统路径：不分 mA 的 ntuples
        for sample in analyzer_cfg.samp_names:
            chain = ntuples[sample]
            df = ROOT.RDataFrame(chain)
            df = _apply_common_filters(df, args)
            df = _define_sys_weights(df, sample)
            df = _apply_blind(df, sample, args)
            df = _define_param(df, sample)

            # 定义所有 target_masses 的 MVA 列（若使用）
            mva_fillers = []
            if args.mva:
                for ALP_mass in target_masses:
                    br = mva_branch_map.get(sample, {}).get(ALP_mass, "MVA_Score")
                    if df.GetColumnType(br) == "":
                        df = df.Define(f"__mva_{ALP_mass}", "-1.0")
                        mva_col = f"__mva_{ALP_mass}"
                    else:
                        mva_col = br
                    # 若 args.cut 为真，仅对当前 args.mA 应用 cut
                    if args.cut and (ALP_mass == args.mA):
                        # 直接 filter，在后续所有直方图上生效
                        cutval = analyzer_cfg.mvaCut[ALP_mass]
                        df = df.Filter(f"{mva_col} >= {cutval}")

                    # 为该 mass 准备直方图（同上）
                    mva_bins = (240, -0.1, 1.1)
                    mva_bins_large = (10, 0., 1.0)
                    sigmaL = 125. + sigma_low[ALP_mass]
                    sigmaH = 125. + sigma_hig[ALP_mass]
                    mva_fillers += [
                        {"name": f"mvaVal_{ALP_mass}",              "expr": mva_col, "bins": mva_bins},
                        {"name": f"mvaVal_larger_{ALP_mass}",       "expr": mva_col, "bins": mva_bins_large},
                        {"name": f"mvaVal_1sigma_{ALP_mass}",        "expr": mva_col, "bins": mva_bins,       "mask": f"(H_m>{sigmaL}) && (H_m<{sigmaH})"},
                        {"name": f"mvaVal_larger_1sigma_{ALP_mass}", "expr": mva_col, "bins": mva_bins_large, "mask": f"(H_m>{sigmaL}) && (H_m<{sigmaH})"},
                        {"name": f"mvaVal_1P5sigma_{ALP_mass}",        "expr": mva_col, "bins": mva_bins,       "mask": f"(H_m>{125.+1.5*sigma_low[ALP_mass]}) && (H_m<{125.+1.5*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_larger_1P5sigma_{ALP_mass}", "expr": mva_col, "bins": mva_bins_large, "mask": f"(H_m>{125.+1.5*sigma_low[ALP_mass]}) && (H_m<{125.+1.5*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_2sigma_{ALP_mass}",        "expr": mva_col, "bins": mva_bins,       "mask": f"(H_m>{125.+2.0*sigma_low[ALP_mass]}) && (H_m<{125.+2.0*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_larger_2sigma_{ALP_mass}", "expr": mva_col, "bins": mva_bins_large, "mask": f"(H_m>{125.+2.0*sigma_low[ALP_mass]}) && (H_m<{125.+2.0*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_3sigma_{ALP_mass}",        "expr": mva_col, "bins": mva_bins,       "mask": f"(H_m>{125.+3.0*sigma_low[ALP_mass]}) && (H_m<{125.+3.0*sigma_hig[ALP_mass]})"},
                        {"name": f"mvaVal_larger_3sigma_{ALP_mass}", "expr": mva_col, "bins": mva_bins_large, "mask": f"(H_m>{125.+3.0*sigma_low[ALP_mass]}) && (H_m<{125.+3.0*sigma_hig[ALP_mass]})"},
                    ]

            _fill_hists_from_df(df, sample, var_names, analyzer_cfg.sys_names,
                                histos, histos_sys, mva_fillers=mva_fillers)
    # ========= RDataFrame 版本结束 =========

    ## End of for sample in analyzer_cfg.samp_names

    # 新增：MVA 直方圖摘要（每個 mA x 每個樣本）
    if args.mva:
        print("\n[MVA] Per-mass histogram summary:")
        for ALP_mass in target_masses:
            hname = 'mvaVal_' + ALP_mass
            if hname not in histos:
                continue
            for s in analyzer_cfg.samp_names:
                h = histos[hname].get(s)
                if not h:
                    continue
                print(f"  mass={ALP_mass:>3s} sample={s:>10s} Entries={h.GetEntries():8.0f}  Integral={h.Integral():.3f}")

    for var in ['H_m', 'pho1Pt', 'Z_m']:
        entries = histos[var][sample].GetEntries()
        print(f"[{sample}] {var} histogram entries: {entries}")


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
        # 依變數選擇縮放策略：H_m 用 sideband，其他維持原本全區間的 ScaleBkgToData
        if var_name == 'H_m':
            scale_factor = SideBandScaleBkgToData(histos, histos_sys, analyzer_cfg, signal_low=115., signal_high=135.)
        else:
            scale_factor = ScaleBkgToData(histos[var_name], analyzer_cfg, histos_sys.get(var_name))
        if scale_factor != 1.0:
            print(f"[Info] Applied background scaling ({var_name}) = {scale_factor:.4f}")

        stacks = MakeStack(histos[var_name], analyzer_cfg, var_name)
        
        if stacks['all'].GetStack().GetEntries() == 0:
            stack_entry = stacks['all'].GetStack().GetEntries()
            print(f"stack_entry: {stack_entry}")
            print(f"[Warning] Stack for {var_name} is empty. Skipping drawing.")
            continue  # 空 stack 直接跳過後續

        for sample in analyzer_cfg.sig_names:
            scaled_sig[sample] = ScaleSignal(plot_cfg, stacks[sample], histos[var_name][sample], var_name)
        ratio_plot = MakeRatioPlot(histos[var_name]['Data'], stacks['all'].GetStack().Last(), var_name)
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


