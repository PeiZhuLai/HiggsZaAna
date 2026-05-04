####################################################
####################################################

import os
import sys
import numpy as np
import gc
import time  # NEW

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from Plot_Helper import LoadNtuples, MakeStack, CreateCanvas, DrawOnCanv, SaveCanvPic, MakeLumiLabel, MakeCMSDASLabel, ScaleSignal, MakeRatioPlot, MakeLegend, Total_Unc, ScaleBkgToData
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

pdfName_map = {'pho1Pt': '1_pho1Pt', 'pho1R9': '2_pho1R9', 'pho1IetaIeta55': '3_pho1IetaIeta55', 
'pho2Pt': '4_pho2Pt', 'pho2R9': '5_pho2R9', 'pho2IetaIeta55': '6_pho2IetaIeta55',
'pho1ECALIso': '7_pho1ECALIso', 'pho2ECALIso': '8_pho2ECALIso', 'ALP_calculatedPhotonIso': '9_ALP_calculatedPhotonIso',
'var_dR_Za': '10_var_dR_Za', 'var_dR_g1g2': '11_var_dR_g1g2', 'var_dR_g1Z': '12_var_dR_g1Z',
'var_PtaOverMh': '13_var_PtaOverMh', 'H_pt': '14_H_pt', 'param': '15_param', 'ALP_m': '16_ALP_m'}

# NEW: 輸出檔名依 pdfName_map 排序；沒對應就用原本 var_name
def _pdf_output_name(var_name: str) -> str:
    return pdfName_map.get(var_name, var_name)

MVA_LARGER_NBINS = 20
MVA_LARGER_XMIN = 0.0
MVA_LARGER_XMAX = 1.0
MVA_FULL_RANGE_BLIND_BINS = 2
MVA_FULL_RANGE_DATA_BLIND_THRESHOLD = (
    MVA_LARGER_XMAX
    - MVA_FULL_RANGE_BLIND_BINS * (MVA_LARGER_XMAX - MVA_LARGER_XMIN) / MVA_LARGER_NBINS
)

def _is_full_range_mva_larger_var(var_name, target_masses):
    return var_name in [f"mvaVal_larger_{mass}" for mass in target_masses]

def _visible_integral(hist):
    if not hist:
        return 0.0
    return float(hist.Integral(1, hist.GetNbinsX()))

def _cdf_edges_from_reference_hist(reference_hist):
    if not reference_hist or _visible_integral(reference_hist) <= 0.0:
        return None

    nbins = reference_hist.GetNbinsX()
    bin_weights = [
        max(0.0, float(reference_hist.GetBinContent(i_bin)))
        for i_bin in range(1, nbins + 1)
    ]
    total = sum(bin_weights)
    if total <= 0.0:
        return None

    cdf_edges = [0.0]
    running = 0.0
    for weight in bin_weights:
        running += weight
        cdf_edges.append(min(1.0, max(0.0, running / total)))
    cdf_edges[-1] = 1.0
    return cdf_edges

def _make_cdf_transformed_hist(hist, name, cdf_edges):
    transformed_hist = hist.Clone(name)
    transformed_hist.Reset("ICES")
    transformed_hist.SetDirectory(0)

    if not cdf_edges:
        return transformed_hist

    out_axis = transformed_hist.GetXaxis()
    out_nbins = transformed_hist.GetNbinsX()
    for i_bin in range(1, hist.GetNbinsX() + 1):
        content = float(hist.GetBinContent(i_bin))
        error = float(hist.GetBinError(i_bin))
        if content == 0.0 and error == 0.0:
            continue

        u_low = min(1.0, max(0.0, cdf_edges[i_bin - 1]))
        u_high = min(1.0, max(0.0, cdf_edges[i_bin]))
        if u_high <= u_low:
            out_bin = transformed_hist.FindFixBin(u_high)
            out_bin = min(out_nbins, max(1, out_bin))
            transformed_hist.SetBinContent(out_bin, transformed_hist.GetBinContent(out_bin) + content)
            transformed_hist.SetBinError(out_bin, np.sqrt(transformed_hist.GetBinError(out_bin) ** 2 + error ** 2))
            continue

        width = u_high - u_low
        for out_bin in range(1, out_nbins + 1):
            bin_low = out_axis.GetBinLowEdge(out_bin)
            bin_high = out_axis.GetBinUpEdge(out_bin)
            overlap = max(0.0, min(u_high, bin_high) - max(u_low, bin_low))
            if overlap <= 0.0:
                continue
            frac = overlap / width
            transformed_hist.SetBinContent(out_bin, transformed_hist.GetBinContent(out_bin) + content * frac)
            transformed_hist.SetBinError(out_bin, np.sqrt(transformed_hist.GetBinError(out_bin) ** 2 + (error * frac) ** 2))

    transformed_hist.SetEntries(hist.GetEntries())
    return transformed_hist

def _build_cdf_histo_maps(var_name, histos_var, histos_sys_var, analyzer_cfg):
    """
    Build normal histograms after a CDF-based score transform.
    The reference CDF is taken from the matching signal mass when available,
    matching the usual BDT score_t convention.
    """
    mass_tag = var_name.split("_")[-1]
    reference_hist = histos_var.get(mass_tag) or histos_var.get("Data")
    cdf_edges = _cdf_edges_from_reference_hist(reference_hist)

    histos_cdf = {}
    histos_sys_cdf = {}

    for sample in analyzer_cfg.samp_names:
        histos_cdf[sample] = _make_cdf_transformed_hist(histos_var[sample], f"{var_name}_{sample}_cdf", cdf_edges)
        histos_sys_cdf[sample] = {}

        for sys_name in analyzer_cfg.sys_names:
            histos_sys_cdf[sample][sys_name] = _make_cdf_transformed_hist(
                histos_sys_var[sample][sys_name],
                f"{var_name}_{sample}_{sys_name}_cdf",
                cdf_edges,
            )

    return histos_cdf, histos_sys_cdf

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
    run3 更新：同一個 era 檔內同時有多個 mA 的 MVA 分數：
      - 優先：MVA_Score_mA_<M>  (e.g. MVA_Score_mA_M4)
    其餘命名仍保留回退相容。
    """
    if not chain or not chain.GetListOfBranches():
        return "MVA_Score"
    br_list = chain.GetListOfBranches()
    candidates = [
        f"MVA_Score_mA_{mass_tag}",   # NEW primary
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
    start_time = time.time()  # NEW

    analyzer_cfg = AC.Analyzer_Config('inclusive', args.year, args.region, args.mva)

    analyzer_cfg.mva = bool(args.mva)
    analyzer_cfg.mva_alp_mass = str(args.mA) if args.mva else "M1"

    # 在 mva + run3/run3_NFlow 下，自動按 self.sig_names 跑每個 mA；其他情況維持單一目標質量
    if args.mva and analyzer_cfg.year in ['run3', 'run3_NFlow']:
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

    # NEW: 一律走 LoadNtuples（不再 split_by_mass；檔案本身就含多 mA 的 MVA branch）
    ntuples = LoadNtuples(analyzer_cfg)

    # NEW: 印出使用的樣本/背景清單，並檢查 ntuples 實際載入到哪些 sample
    print("\n[Samples] analyzer_cfg lists:")
    print(f"  bkg_names({len(analyzer_cfg.bkg_names)}): {analyzer_cfg.bkg_names}")
    print(f"  sig_names({len(analyzer_cfg.sig_names)}): {analyzer_cfg.sig_names}")
    print(f"  samp_names({len(analyzer_cfg.samp_names)}): {analyzer_cfg.samp_names}")
    print(f"  sys_names({len(analyzer_cfg.sys_names)}): {analyzer_cfg.sys_names}")

    loaded_keys = sorted(list(ntuples.keys())) if ntuples else []
    print("\n[Samples] ntuples loaded keys:")
    print(f"  loaded({len(loaded_keys)}): {loaded_keys}")

    cfg_set = set(analyzer_cfg.samp_names)
    loaded_set = set(loaded_keys)
    missing_in_ntuples = sorted(list(cfg_set - loaded_set))
    extra_in_ntuples = sorted(list(loaded_set - cfg_set))
    if missing_in_ntuples:
        print("\n[Warning][Samples] In cfg.samp_names but NOT loaded by LoadNtuples:")
        print(f"  missing({len(missing_in_ntuples)}): {missing_in_ntuples}")
    if extra_in_ntuples:
        print("\n[Info][Samples] Loaded by LoadNtuples but NOT in cfg.samp_names:")
        print(f"  extra({len(extra_in_ntuples)}): {extra_in_ntuples}")

    # 另外特別看 bkg 是否有缺
    bkg_expected = set([s for s in analyzer_cfg.bkg_names if s.lower() != "data"])
    bkg_loaded = set([s for s in loaded_keys if (s.lower() != "data") and (s in analyzer_cfg.bkg_names)])
    bkg_missing = sorted(list(bkg_expected - bkg_loaded))
    print("\n[Samples] bkg cross-check:")
    print(f"  expected_bkg({len(bkg_expected)}): {sorted(list(bkg_expected))}")
    print(f"  loaded_bkg({len(bkg_loaded)}): {sorted(list(bkg_loaded))}")
    if bkg_missing:
        print(f"  [Warning] missing_bkg({len(bkg_missing)}): {bkg_missing}")

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

    # 為每個 (sample, mA) 解析一次 MVA 分支（從同一條 chain 上找 MVA_Score_mA_<mass>）
    mva_branch_map = {}  # dict[sample][mass] -> branch_name
    if args.mva:
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
        histos['pho1Pt'][sample]    = TH1F('pho1Pt'    + '_' + sample, 'pho1Pt'    + '_' + sample, 42,  8., 50.)
        histos['pho1eta'][sample]    = TH1F('pho1eta'    + '_' + sample, 'pho1eta'    + '_' + sample, 30,  -3., 3.)
        histos['pho1phi'][sample]    = TH1F('pho1phi'    + '_' + sample, 'pho1phi'    + '_' + sample, 40,  -4., 4.)
        histos['pho1R9'][sample]    = TH1F('pho1R9'    + '_' + sample, 'pho1R9'    + '_' + sample, 25,  0.1, 1.)
        histos['pho1IetaIeta55'][sample]    = TH1F('pho1IetaIeta55'    + '_' + sample, 'pho1IetaIeta55'    + '_' + sample, 25,  0., 0.08)
        histos['pho1ECALIso'][sample]    = TH1F('pho1ECALIso'    + '_' + sample, 'pho1ECALIso'    + '_' + sample, 50, 0., 50.)
        histos['pho1CIso'][sample]    = TH1F('pho1CIso'    + '_' + sample, 'pho1CIso'    + '_' + sample, 35, 0., 0.35)
        histos['pho1HCALIso'][sample]  = TH1F('pho1HCALIso'    + '_' + sample, 'pho1HCALIso'    + '_' + sample, 25, 0., 4.0)
        histos['pho1HOE'][sample]    = TH1F('pho1HOE'    + '_' + sample, 'pho1HOE'    + '_' + sample, 25, 0., 0.01)
        histos['pho2Pt'][sample]    = TH1F('pho2Pt'    + '_' + sample, 'pho2Pt'    + '_' + sample, 22,  8., 30.)
        histos['pho2eta'][sample]    = TH1F('pho2eta'    + '_' + sample, 'pho2eta'    + '_' + sample, 30,  -3., 3.)
        histos['pho2phi'][sample]    = TH1F('pho2phi'    + '_' + sample, 'pho2phi'    + '_' + sample, 40,  -4., 4.)
        histos['pho2R9'][sample]    = TH1F('pho2R9'    + '_' + sample, 'pho2R9'    + '_' + sample, 25,  0.1, 1.)
        histos['pho2IetaIeta55'][sample]    = TH1F('pho2IetaIeta55'    + '_' + sample, 'pho2IetaIeta55'    + '_' + sample, 25,  0., 0.08)
        histos['pho2ECALIso'][sample]    = TH1F('pho2ECALIso'    + '_' + sample, 'pho2ECALIso'    + '_' + sample, 50, 0., 50.)
        histos['pho2CIso'][sample]    = TH1F('pho2CIso'    + '_' + sample, 'pho2CIso'    + '_' + sample, 35, 0., 0.35)
        histos['pho2HCALIso'][sample]    = TH1F('pho2HCALIso'    + '_' + sample, 'pho2HCALIso'    + '_' + sample, 25, 0., 4.0)
        histos['pho2HOE'][sample]    = TH1F('pho2HOE'    + '_' + sample, 'pho2HOE'    + '_' + sample, 25, 0., 0.1)
        histos['Z_m'][sample]    = TH1F('Z_m'    + '_' + sample, 'Z_m'    + '_' + sample, 80,  50., 130.)
        histos['H_m'][sample]    = TH1F('H_m'    + '_' + sample, 'H_m'    + '_' + sample, 85,  95., 180.)
        histos['H_pt'][sample]    = TH1F('H_pt'    + '_' + sample, 'H_pt'    + '_' + sample, 80,  0., 160.)
        histos['ALP_m'][sample] = TH1F('ALP_m' + '_' + sample, 'ALP_m' + '_' + sample, 40, 0., 40.)
        histos['var_dR_g1g2'][sample] = TH1F('var_dR_g1g2' + '_' + sample, 'var_dR_g1g2' + '_' + sample, 25, 0., 5)
        histos['var_PtaOverMa'][sample] = TH1F('var_PtaOverMa' + '_' + sample, 'var_PtaOverMa' + '_' + sample, 25, 0., 100.)
        histos['var_dR_Za'][sample] = TH1F('var_dR_Za' + '_' + sample, 'var_dR_Za' + '_' + sample, 35, 0., 7.)
        histos['var_dR_g1Z'][sample] = TH1F('var_dR_g1Z' + '_' + sample, 'var_dR_g1Z' + '_' + sample, 35, 0., 7)
        histos['var_PtaOverMh'][sample] = TH1F('var_PtaOverMh' + '_' + sample, 'var_PtaOverMh' + '_' + sample, 25, 0., 0.75)
        histos['var_Pta'][sample] = TH1F('var_Pta' + '_' + sample, 'var_Pta' + '_' + sample, 60, 0., 60.)
        histos['var_MhMa'][sample] = TH1F('var_MhMa' + '_' + sample, 'var_MhMa' + '_' + sample, 100, 100., 200.)
        histos['var_MhMZ'][sample] = TH1F('var_MhMZ' + '_' + sample, 'var_MhMZ' + '_' + sample, 130, 180., 310.)
        histos['ALP_calculatedPhotonIso'][sample] = TH1F('ALP_calculatedPhotonIso' + '_' + sample, 'ALP_calculatedPhotonIso' + '_' + sample, 25, 0., 125.)
        histos['param'][sample] = TH1F('param' + '_' + sample, 'param' + '_' + sample, 25, -0.3, 0.6)

        if args.mva:
            # 僅為選定質量建立 MVA 相關直方圖
            for ALP_mass in target_masses:
                histos['mvaVal_'+ALP_mass][sample]           = TH1F('mvaVal_'+ALP_mass           + '_' + sample, 'mvaVal_'+ALP_mass    + '_' + sample, 240, -0.1, 1.1)
                histos['mvaVal_1sigma_'+ALP_mass][sample]    = TH1F('mvaVal_1sigma_'+ALP_mass    + '_' + sample, 'mvaVal_1sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_1P5sigma_'+ALP_mass][sample]  = TH1F('mvaVal_1P5sigma_'+ALP_mass  + '_' + sample, 'mvaVal_1P5sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_2sigma_'+ALP_mass][sample]    = TH1F('mvaVal_2sigma_'+ALP_mass    + '_' + sample, 'mvaVal_2sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)
                histos['mvaVal_3sigma_'+ALP_mass][sample]    = TH1F('mvaVal_3sigma_'+ALP_mass    + '_' + sample, 'mvaVal_3sigma_'+ALP_mass    + '_' + sample, 240,  -0.1, 1.1)

                histos['mvaVal_larger_'+ALP_mass][sample]           = TH1F('mvaVal_larger_'+ALP_mass           + '_' + sample, 'mvaVal_larger_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)
                histos['mvaVal_larger_1sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_1sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_1sigma_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)
                histos['mvaVal_larger_1P5sigma_'+ALP_mass][sample]  = TH1F('mvaVal_larger_1P5sigma_'+ALP_mass  + '_' + sample, 'mvaVal_larger_1P5sigma_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)
                histos['mvaVal_larger_2sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_2sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_2sigma_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)
                histos['mvaVal_larger_3sigma_'+ALP_mass][sample]    = TH1F('mvaVal_larger_3sigma_'+ALP_mass    + '_' + sample, 'mvaVal_larger_3sigma_'+ALP_mass    + '_' + sample, MVA_LARGER_NBINS, MVA_LARGER_XMIN, MVA_LARGER_XMAX)


    for var_name in var_names:
        histos_sys[var_name] = {}
        for sample in analyzer_cfg.samp_names:
            histos_sys[var_name][sample] = {}
            for sys in analyzer_cfg.sys_names:
                histos_sys[var_name][sample][sys] = copy.deepcopy(histos[var_name][sample])
                histos_sys[var_name][sample][sys].SetNameTitle(var_name+'_'+sample+'_'+sys, var_name+'_'+sample+'_'+sys)

    ### loop over samples and events
    mass_list = {'M1':1.0, 'M2':2.0, 'M3':3.0, 'M4':4.0, 'M5':5.0, 'M6':6.0, 'M7':7.0, 'M8':8.0, 'M9':9.0, 'M10':10.0, 'M15':15.0, 'M20':20.0, 'M25':25.0, 'M30':30.0}
    search_mA_list = [float(m) for m in range(1, 31)] 
    # 新增：控制每個 (sample, mA) 的偵錯輸出次數
    debug_printed = {}

    # NEW: 移除/停用 split-by-mass 的事件迴圈，統一用原本單一鏈流程
    for sample in analyzer_cfg.samp_names:
        ntup = ntuples[sample] # just a short name
        print('\n\nOn sample: %s' %sample)
        print('total events: %d' %ntup.GetEntries())

        for iEvt in range( ntup.GetEntries() ):
    
            ntup.GetEvent(iEvt)
            # if (iEvt == 10): break

            if (iEvt % 100000 == 1):
                print("looking at event %d" %iEvt)

            
            if args.ele:
                if abs(ntup.z_mumu) == 1: 
                    continue
            if args.mu:
                if abs(ntup.z_ee) == 1: 
                    continue
            

            # weight = ntup.factor * ntup.pho1SFs * ntup.pho2SFs
            weight = ntup.weight

            if (ntup.H_m > -90):
                if ntup.H_m>180. or ntup.H_m<95.: continue

                if  args.region == 1 and (ntup.H_m>135. or ntup.H_m<115.): continue
                if  args.region == 2 and (ntup.H_m<135. and ntup.H_m>115.): continue

                MVA_value = {}
                if args.mva:
                    # 以各 mA 對應分支取得 MVA 分數（若分支不存在，回退 "MVA_Score"，再不行設為 -1）
                    for ALP_mass in target_masses:
                        br = mva_branch_map.get(sample, {}).get(ALP_mass, "MVA_Score")
                        try:
                            MVA_value[ALP_mass] = getattr(ntup, br)
                        except Exception:
                            # 兩層保險：先試指定分支，再試通用分支，最後給 -1.0
                            try:
                                MVA_value[ALP_mass] = getattr(ntup, "MVA_Score")
                            except Exception:
                                MVA_value[ALP_mass] = -1.0

                        # 新增：每個 (sample, mA) 前 args.mva_debug_n 筆偵錯印出
                        if args.mva_debug:
                            key = (sample, ALP_mass)
                            if debug_printed.get(key, 0) < args.mva_debug_n:
                                print(f"[MVA-FILL] sample={sample:>10s} mA={ALP_mass:>3s} branch={br:>20s} score={MVA_value[ALP_mass]:+.4f} H_m={ntup.H_m:6.2f} w={weight: .3e}")
                                debug_printed[key] = debug_printed.get(key, 0) + 1

                    if args.cut:
                        if MVA_value[args.mA] < analyzer_cfg.mvaCut[args.mA]: continue

                    for ALP_mass in target_masses:
                        if args.blind:
                            if not (sample == 'Data' and MVA_value[ALP_mass] >= MVA_FULL_RANGE_DATA_BLIND_THRESHOLD):
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

                var_map = {'Z_m':ntup.Z_mass, 'H_m':ntup.H_m, 'ALP_m':ntup.ALP_m,'pho1Pt':ntup.pho1Pt, 'pho1eta':ntup.ALP_lead_photon_eta, 'pho1phi':ntup.ALP_lead_photon_phi, 'pho1R9':ntup.ALP_lead_photon_r9, 'pho1IetaIeta':ntup.ALP_lead_photon_sieie, 'pho1IetaIeta55':ntup.ALP_lead_photon_sieie,'pho1ECALIso':ntup.ALP_lead_photon_ecalPFClusterIso, 'pho1CIso':ntup.ALP_lead_photon_chiso, 'pho1HCALIso':ntup.ALP_lead_photon_hcalPFClusterIso, 'pho1HOE':ntup.ALP_lead_photon_hoe_PUcorr, 'pho2Pt':ntup.ALP_sublead_photon_pt, 'pho2eta':ntup.ALP_sublead_photon_eta, 'pho2phi':ntup.ALP_sublead_photon_phi, 'pho2R9':ntup.ALP_sublead_photon_r9, 'pho2IetaIeta':ntup.ALP_sublead_photon_sieie, 'pho2IetaIeta55':ntup.ALP_sublead_photon_sieie,'pho2ECALIso':ntup.ALP_sublead_photon_ecalPFClusterIso, 'pho2CIso':ntup.ALP_sublead_photon_chiso, 'pho2HCALIso':ntup.ALP_sublead_photon_hcalPFClusterIso, 'pho2HOE':ntup.ALP_sublead_photon_hoe_PUcorr,'ALP_calculatedPhotonIso':ntup.ALP_calculatedPhotonIso, 'var_dR_Za':ntup.var_dR_Za, 'var_dR_g1g2':ntup.var_dR_g1g2, 'var_dR_g1Z':ntup.var_dR_g1Z, 'var_PtaOverMh':ntup.var_PtaOverMh, 'var_Pta':ntup.var_Pta, 'var_MhMZ':ntup.var_MhMZ, 'H_pt':ntup.H_pt, 'var_PtaOverMa':ntup.var_PtaOverMa, 'var_MhMa':ntup.var_MhMa}
                
                if args.mva:
                    var_map_mva = {}
                    for ALP_mass in target_masses:
                        var_map_mva['mvaVal_'+ALP_mass] = MVA_value[ALP_mass]
                        var_map_mva['mvaVal_larger_'+ALP_mass] = MVA_value[ALP_mass]
                        for r in ['1sigma', '1P5sigma', '2sigma', '3sigma']:
                            var_map_mva['mvaVal_'+r+'_'+ALP_mass] = MVA_value[ALP_mass]
                            var_map_mva['mvaVal_larger_'+r+'_'+ALP_mass] = MVA_value[ALP_mass]
                    var_map.update(var_map_mva)

                histos['pho1Pt'][sample].Fill( ntup.ALP_lead_photon_pt, weight )
                histos['pho1eta'][sample].Fill( ntup.ALP_lead_photon_eta, weight )
                histos['pho1phi'][sample].Fill( ntup.ALP_lead_photon_phi, weight )
                histos['pho1R9'][sample].Fill( ntup.ALP_lead_photon_r9, weight )
                histos['pho1IetaIeta55'][sample].Fill( ntup.ALP_lead_photon_sieie, weight )
                histos['pho1ECALIso'][sample].Fill( ntup.ALP_lead_photon_ecalPFClusterIso, weight )
                histos['pho2Pt'][sample].Fill( ntup.ALP_sublead_photon_pt, weight )
                histos['pho2eta'][sample].Fill( ntup.ALP_sublead_photon_eta, weight )
                histos['pho2phi'][sample].Fill( ntup.ALP_sublead_photon_phi, weight )
                histos['pho2R9'][sample].Fill( ntup.ALP_sublead_photon_r9, weight )
                histos['pho2IetaIeta55'][sample].Fill( ntup.ALP_sublead_photon_sieie, weight )
                histos['pho2ECALIso'][sample].Fill( ntup.ALP_sublead_photon_ecalPFClusterIso, weight )

                histos['pho1CIso'][sample].Fill( ntup.ALP_lead_photon_chiso, weight)
                histos['pho1HCALIso'][sample].Fill( ntup.ALP_lead_photon_hcalPFClusterIso, weight)
                histos['pho1HOE'][sample].Fill( ntup.ALP_lead_photon_hoe_PUcorr, weight)
                histos['pho2CIso'][sample].Fill( ntup.ALP_sublead_photon_chiso, weight)
                histos['pho2HCALIso'][sample].Fill( ntup.ALP_sublead_photon_hcalPFClusterIso, weight)
                histos['pho2HOE'][sample].Fill( ntup.ALP_sublead_photon_hoe_PUcorr, weight)

                if args.blind:
                    if not (sample == 'Data' and (ntup.H_m<135. and ntup.H_m>115.)): 
                        histos['H_m'][sample].Fill( ntup.H_m, weight )
                else:        
                    histos['H_m'][sample].Fill( ntup.H_m, weight )

                histos['H_pt'][sample].Fill( ntup.H_pt, weight )
                histos['ALP_m'][sample].Fill( ntup.ALP_m, weight )
                histos['Z_m'][sample].Fill( ntup.Z_mass, weight )

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
                    # mass_random = random.choice(search_mA_list)
                    mass_random = random.choice(list(mass_list.values()))
                    param_val['param'] = (ntup.ALP_m - mass_random)/ntup.H_m
                
                var_map.update(param_val)

                histos['param'][sample].Fill( param_val['param'], weight )
                    

                for sys_name in analyzer_cfg.sys_names:
                    if sample != "Data": 
                        if sys_name =='weight_hlt_sf_up':
                            weight_sys = weight * ntup.weight_hlt_sf_up / ntup.weight_hlt_sf_central
                        elif sys_name =='weight_hlt_sf_down':
                            weight_sys = weight * ntup.weight_hlt_sf_down / ntup.weight_hlt_sf_central
                        elif sys_name =='weight_pu_reweight_sf_up':
                            weight_sys = weight * ntup.weight_pu_reweight_sf_up / ntup.weight_pu_reweight_sf_central
                        elif sys_name =='weight_pu_reweight_sf_down':
                            weight_sys = weight * ntup.weight_pu_reweight_sf_down / ntup.weight_pu_reweight_sf_central
                        elif sys_name =='weight_electron_wplid_sf_SelectedElectron_up':
                            weight_sys = weight * ntup.weight_electron_wplid_sf_SelectedElectron_up / ntup.weight_electron_wplid_sf_SelectedElectron_central
                        elif sys_name =='weight_electron_wplid_sf_SelectedElectron_down':
                            weight_sys = weight * ntup.weight_electron_wplid_sf_SelectedElectron_down / ntup.weight_electron_wplid_sf_SelectedElectron_central
                        elif sys_name =='weight_electron_reco_sf_SelectedElectron_up':
                            weight_sys = weight * ntup.weight_electron_reco_sf_SelectedElectron_up / ntup.weight_electron_reco_sf_SelectedElectron_central                        
                        elif sys_name =='weight_electron_reco_sf_SelectedElectron_down':
                            weight_sys = weight * ntup.weight_electron_reco_sf_SelectedElectron_down / ntup.weight_electron_reco_sf_SelectedElectron_central
                        elif sys_name =='weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_up':
                            weight_sys = weight * ntup.weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_up / ntup.weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_central
                        elif sys_name =='weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_down':
                            weight_sys = weight * ntup.weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_down / ntup.weight_electron_wplid_sf_nomatch_SelectedGenNoRecoElectron_central

                        elif sys_name =='weight_muon_looseid_sf_SelectedMuon_up':
                            weight_sys = weight * ntup.weight_muon_looseid_sf_SelectedMuon_up / ntup.weight_muon_looseid_sf_SelectedMuon_central
                        elif sys_name =='weight_muon_looseid_sf_SelectedMuon_down':
                            weight_sys = weight * ntup.weight_muon_looseid_sf_SelectedMuon_down / ntup.weight_muon_looseid_sf_SelectedMuon_central
                        elif sys_name =='weight_muon_reco_sf_SelectedMuon_up':
                            weight_sys = weight * ntup.weight_muon_reco_sf_SelectedMuon_up / ntup.weight_muon_reco_sf_SelectedMuon_central
                        elif sys_name =='weight_muon_reco_sf_SelectedMuon_down':
                            weight_sys = weight * ntup.weight_muon_reco_sf_SelectedMuon_down / ntup.weight_muon_reco_sf_SelectedMuon_central
                        elif sys_name =='weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_up':
                            weight_sys = weight * ntup.weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_up / ntup.weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_central
                        elif sys_name =='weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_down':
                            weight_sys = weight * ntup.weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_down / ntup.weight_muon_looseid_sf_nomatch_SelectedGenNoRecoMuon_central

                        elif sys_name =='weight_photon_id_sf_SelectedPhoton_up':
                            weight_sys = weight * ntup.weight_photon_id_sf_SelectedPhoton_up / ntup.weight_photon_id_sf_SelectedPhoton_central
                        elif sys_name =='weight_photon_id_sf_SelectedPhoton_down':
                            weight_sys = weight * ntup.weight_photon_id_sf_SelectedPhoton_down / ntup.weight_photon_id_sf_SelectedPhoton_central
                        elif sys_name =='weight_photon_csev_sf_SelectedPhoton_up':
                            weight_sys = weight * ntup.weight_photon_csev_sf_SelectedPhoton_up / ntup.weight_photon_csev_sf_SelectedPhoton_central
                        elif sys_name =='weight_photon_csev_sf_SelectedPhoton_down':
                            weight_sys = weight * ntup.weight_photon_csev_sf_SelectedPhoton_down / ntup.weight_photon_csev_sf_SelectedPhoton_central

                        for var in var_names:
                            histos_sys[var][sample][sys_name].Fill(var_map[var], weight_sys)
                    else:
                        for var in var_names:
                            histos_sys[var][sample][sys_name].Fill(var_map[var], 1.)
                



        ## End of for iEvt in range( ntup.GetEntries() )
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
            SaveCanvPic(canv_log, analyzer_cfg.plot_output_path, _pdf_output_name(var_name) + '_log')
        else:
            canv = CreateCanvas(var_name)
            DrawOnCanv(canv, var_name, plot_cfg, stacks, histos[var_name], scaled_sig, ratio_plot, legend, lumi_label, cms_label, total_unc, args.cut, args.mA, logY=False)

            canv.Write()
            SaveCanvPic(canv, analyzer_cfg.plot_output_path, _pdf_output_name(var_name))

        if args.mva and _is_full_range_mva_larger_var(var_name, target_masses):
            histos_cdf, histos_sys_cdf = _build_cdf_histo_maps(
                var_name,
                histos[var_name],
                histos_sys[var_name],
                analyzer_cfg,
            )
            cdf_var_name = var_name + "_cdf"
            stacks_cdf = MakeStack(histos_cdf, analyzer_cfg, cdf_var_name)

            if stacks_cdf['all'].GetStack().GetEntries() == 0:
                stack_entry = stacks_cdf['all'].GetStack().GetEntries()
                print(f"stack_entry: {stack_entry}")
                print(f"[Warning] CDF stack for {var_name} is empty. Skipping drawing.")
                continue

            scaled_sig_cdf = {}
            for sample in analyzer_cfg.sig_names:
                scaled_sig_cdf[sample] = ScaleSignal(plot_cfg, stacks_cdf[sample], histos_cdf[sample], cdf_var_name)

            ratio_plot_cdf = MakeRatioPlot(histos_cdf['Data'], stacks_cdf['all'].GetStack().Last(), cdf_var_name)
            legend_cdf = MakeLegend(plot_cfg, histos_cdf, scaled_sig_cdf)
            total_unc_cdf = Total_Unc(stacks_cdf['bkg'], histos_sys_cdf, analyzer_cfg)
            cdf_y_axis_title = "Entries / %.3f" % histos_cdf['Data'].GetBinWidth(1)

            if args.ln:
                canv_cdf = CreateCanvas(cdf_var_name + "_log")
                DrawOnCanv(
                    canv_cdf,
                    cdf_var_name + "_log",
                    plot_cfg,
                    stacks_cdf,
                    histos_cdf,
                    scaled_sig_cdf,
                    ratio_plot_cdf,
                    legend_cdf,
                    lumi_label,
                    cms_label,
                    total_unc_cdf,
                    args.cut,
                    args.mA,
                    logY=True,
                    y_axis_title=cdf_y_axis_title,
                    axis_var_name=var_name,
                    x_axis_title="Transformed BDT Score",
                )
                canv_cdf.Write()
                SaveCanvPic(canv_cdf, analyzer_cfg.plot_output_path, _pdf_output_name(var_name) + "_cdf_log")
            else:
                canv_cdf = CreateCanvas(cdf_var_name)
                DrawOnCanv(
                    canv_cdf,
                    cdf_var_name,
                    plot_cfg,
                    stacks_cdf,
                    histos_cdf,
                    scaled_sig_cdf,
                    ratio_plot_cdf,
                    legend_cdf,
                    lumi_label,
                    cms_label,
                    total_unc_cdf,
                    args.cut,
                    args.mA,
                    logY=False,
                    y_axis_title=cdf_y_axis_title,
                    axis_var_name=var_name,
                    x_axis_title="Transformed BDT Score",
                )
                canv_cdf.Write()
                SaveCanvPic(canv_cdf, analyzer_cfg.plot_output_path, _pdf_output_name(var_name) + "_cdf")


    
    print('\n\n')
    #CountYield(analyzer_cfg, histos['ALP_m'])
    out_file.Close()

    elapsed = time.time() - start_time  # NEW
    print(f"Total runtime: {time.strftime('%H:%M:%S', time.gmtime(elapsed))}")  # NEW
    print('Done')


main()
