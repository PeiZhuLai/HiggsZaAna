import os
from Plot_Helper import MakeStack
import Plot_Configs as PC
import Analyzer_Configs as AC
from ROOT import *

#####################################################################
class MuPair:
    def __init__(self, vec, idx1, idx2):
        self.dimu_vec = vec
        self.mass     = vec.M()
        self.pt       = vec.Pt()
        self.mu1_idx  = idx1
        self.mu2_idx  = idx2

def getMassSigma(ana_cfg):
    """
    為每個訊號樣本計算 H 質量中心 (125 GeV) 週期內的 ±sigma_low / sigma_hig。
    run3: 嘗試 /ALP_{M}/<year>.root (years_sig)；舊樣式: /ALP_{M}/run3.root
    若檔案缺失或無法取得直方圖，使用預設 (-1.5, +1.5) GeV 視為 1σ。
    回傳:
        sigma_low[sample]  (負值)
        sigma_hig[sample]  (正值)
    """
    sigma_low = {}
    sigma_hig = {}

    # 預設 fallback 視窗
    default_low = -1.5
    default_high = 1.5

    for sample in ana_cfg.sig_names:
        # 收集可嘗試的檔案路徑
        candidate_paths = []
        if ana_cfg.year in ['run3', 'run3_NFlow']:
            for y in getattr(ana_cfg, "years_sig", []):
                candidate_paths.append(os.path.join(ana_cfg.sample_loc, f"mA_{sample}", f"{y}.root"))
        # 舊結構 (若仍存在)
        candidate_paths.append(os.path.join(ana_cfg.sample_loc, f"mA_{sample}", "run3.root"))

        opened_file = None
        for p in candidate_paths:
            if os.path.exists(p):
                opened_file = TFile.Open(p)
                if opened_file and not opened_file.IsZombie():
                    print(f"[getMassSigma] 使用檔案: {p}")
                    break
                else:
                    print(f"[getMassSigma][WARN] 檔案損毀或無法開啟: {p}")
                    opened_file = None
            else:
                print(f"[getMassSigma][MISS] {p}")

        if not opened_file:
            print(f"[getMassSigma][FALLBACK] {sample} 使用預設視窗 ({default_low}, {default_high}) GeV")
            sigma_low[sample] = default_low
            sigma_hig[sample] = default_high
            continue

        tree = opened_file.Get("test")
        if not tree or not isinstance(tree, TTree):
            print(f"[getMassSigma][WARN] {sample} 檔案中找不到 TTree 'test'，使用預設視窗")
            sigma_low[sample] = default_low
            sigma_hig[sample] = default_high
            opened_file.Close()
            continue

        # 建立暫時直方圖名稱避免覆寫
        hist_name = f"tmpHm_{sample}"
        draw_expr = f"H_m>>{hist_name}"
        cut_expr = "H_m>-90&&H_m>115&&H_m<135"
        # 使用 factor 權重若存在，否則權重=1
        if tree.GetListOfBranches().FindObject("factor"):
            weight = f"factor*({cut_expr})"
        else:
            weight = cut_expr

        nevt = tree.Draw(draw_expr, weight, "goff")
        if nevt <= 0:
            print(f"[getMassSigma][WARN] {sample} 無事件通過選擇，使用預設視窗")
            sigma_low[sample] = default_low
            sigma_hig[sample] = default_high
            opened_file.Close()
            continue

        hist = gDirectory.Get(hist_name)
        if not hist or hist.GetNbinsX() < 3:
            print(f"[getMassSigma][WARN] {sample} 直方圖無效，使用預設視窗")
            sigma_low[sample] = default_low
            sigma_hig[sample] = default_high
            opened_file.Close()
            continue

        # 與原程式一致：尋找最小對稱 bin-span 使得累積 > 68.3%
        total = hist.Integral()
        if total <= 0:
            print(f"[getMassSigma][WARN] {sample} 直方圖積分為 0，使用預設視窗")
            sigma_low[sample] = default_low
            sigma_hig[sample] = default_high
            opened_file.Close()
            continue

        center_bin = hist.FindBin(125.0)
        sigma_bin = 0
        max_span = min(center_bin - 1, hist.GetNbinsX() - center_bin)

        for i in range(max_span + 1):
            low_bin = center_bin - i
            high_bin = center_bin + i
            if low_bin < 1 or high_bin > hist.GetNbinsX():
                continue
            frac = hist.Integral(low_bin, high_bin) / total
            if frac >= 0.683:
                sigma_bin = i
                break

        # 轉換成 offset (注意保持原本 low 為負值, high 為正值慣例)
        low_edge = hist.GetBinLowEdge(center_bin - sigma_bin)
        high_edge = hist.GetBinLowEdge(center_bin + sigma_bin) + hist.GetBinWidth(center_bin + sigma_bin)
        # 例如 low_edge ~ 123.7 -> offset = 123.7 - 125 = -1.3
        sigma_low[sample] = low_edge - 125.0
        sigma_hig[sample] = high_edge - 125.0

        print(f"[getMassSigma] {sample}: low={sigma_low[sample]:.3f}  high={sigma_hig[sample]:.3f}")

        opened_file.Close()

    return sigma_low, sigma_hig
