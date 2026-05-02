# ==== New imports/config for plotting ====
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot
import json
from ROOT import TH1F, TGraph, TGraphSmooth, gROOT
gROOT.SetBatch(True)
from array import array

INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal/"

# 你之前实际使用的组合：
years_sig  = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]  # 信号
years_22   = ["2022preEE", "2022postEE"]  # 背景（2022）
years_23   = ["2023preBPix","2023postBPix", "2024"]
years_dyll = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]

# 扫描的 ma 列表（背景也按你的原脚本扫）
current_ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
interploate_ma_list = [11,12,13,14,16,17,18,19,21,22,23,24,26,27,28,29]
# interploate_ma_list = [11,12]
# interploate_ma_list = [5,15,30]
# ===== Samples =====
# sig_samples = ["mA_M5", "mA_M15", "mA_M30"]
sig_samples = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]

# 背景
bkg_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
bkg_2023 = ["DYGto2LG_10to100"]
bkg_dyll = ["DYJetsToLL"]
# 如需其它：["ZGToLLG","WGToLNuG","ZG2JToG2L2J","EWKZ2J","TT","TTGJets","TGJets","ttWJets","ttZJets","WW","WZ","ZZ"]

MVA_CANDIDATES = ["MVA_Score"]
WEIGHT_CANDIDATES = ["weight"]
optimized_BDT_Cut="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

# 新增：signal 的 BDT binning（用來把連續 cut 值對齊到確切邊界）
SIG_NBINS = 100
SIG_XMIN = -0.1
SIG_XMAX = 1.1

# 新增：是否輸出图档的开关（True: 绘图, False: 不绘图）
ENABLE_PLOTTING = False

# ==== Helpers ====
# 设置全局字体大小和图形样式
plt.rcParams.update({
    'figure.figsize': (8, 6),  # 图形大小

    'font.size': 14,
    'axes.titlesize': 17,  # 图标题，plt.title()
    'axes.labelsize': 18,  # 坐标轴标题/标签，plt.xlabel, plt.ylabel
    'xtick.labelsize': 14,  # x 轴刻度“数字或文字”大小
    'ytick.labelsize': 14,  # y 轴刻度“数字或文字”大小
    'legend.fontsize': 14,  # 图例文字大小 plt.legend()
    'figure.titlesize': 14,  # 整个 Figure 的总标题字体大小，用 plt.suptitle() 设定的顶层标题

    # 'patch.linewidth': 1.5 ,       # 直方图柱子的边框线宽
    # 'patch.edgecolor': 'blue',      # 直方图柱子的边框颜色

})

def mass_from_sig(sample: str) -> float:
    # ALP_M5 -> 5 (float)
    if "_M" in sample:
        try:
            return float(sample.split("_M", 1)[1])
        except Exception:
            pass
    return float("nan")

def first_tree(file_path: str):
    try:
        f = uproot.open(file_path)
    except Exception:
        return None
    # 找第一个 TTree
    for k, cls in f.classnames().items():
        if "TTree" in cls:
            return f[k]
    # 后备名称
    for name in ["Events", "tree", "t"]:
        if name in f:
            return f[name]
    return None

def read_arrays(file_path: str):
    """Return (mva, weight) numpy arrays, or (None, None) if not available."""
    tree = first_tree(file_path)
    if tree is None:
        return None, None
    # MVA
    mva = None
    for v in MVA_CANDIDATES:
        if v in tree.keys():
            try:
                mva = tree[v].array(library="np")
                break
            except Exception:
                continue
    if mva is None:
        return None, None
    # weights
    w = None
    for wv in WEIGHT_CANDIDATES:
        if wv in tree.keys():
            try:
                w = tree[wv].array(library="np")
                break
            except Exception:
                continue
    if w is None:
        w = np.ones_like(mva, dtype=np.float64)
    # 转为 1D
    try:
        mva = np.asarray(mva).astype(np.float64)
        w = np.asarray(w).astype(np.float64)
    except Exception:
        return None, None
    # 对齐长度
    n = min(len(mva), len(w))
    return mva[:n], w[:n]

def read_arrays_with_mass(file_path: str, mass_tag: str):
    """
    参照 plot_variable_dataVmc.py 的候选规则读取 MVA 分支（含 mass_tag），回传 (mva, weight)。
    若存在 H_m 分支，会在 115–135 GeV 视窗内进行过滤。
    """
    tree = first_tree(file_path)
    if tree is None:
        return None, None
    # 分支候选（依序尝试）
    candidates = [
        f"MVA_Score_{mass_tag}",
        f"MVA_score_{mass_tag}",
        f"BDTScore_{mass_tag}",
        f"BDT_{mass_tag}",
        "MVA_Score",
    ]
    mva = None
    for v in candidates:
        if v in tree.keys():
            try:
                mva = tree[v].array(library="np")
                break
            except Exception:
                continue
    if mva is None:
        return None, None

    # 权重
    w = None
    for wv in WEIGHT_CANDIDATES:
        if wv in tree.keys():
            try:
                w = tree[wv].array(library="np")
                break
            except Exception:
                continue
    if w is None:
        w = np.ones_like(mva, dtype=np.float64)

    # 读取 H_m 并在 115–135 GeV 视窗过滤（若分支存在）
    H = None
    try:
        if "H_m" in tree.keys():
            H = tree["H_m"].array(library="np")
    except Exception:
        H = None

    try:
        mva = np.asarray(mva).astype(np.float64)
        w   = np.asarray(w).astype(np.float64)
        if H is not None:
            H = np.asarray(H).astype(np.float64)
    except Exception:
        return None, None

    # 对齐长度（含 H）
    n = min(len(mva), len(w))
    if H is not None:
        n = min(n, len(H))
        H = H[:n]
    mva = mva[:n]
    w   = w[:n]

    # 套用 H_m 视窗遮罩
    if H is not None:
        mask = (H >= 115.0) & (H <= 135.0)
        if mask.size != 0:
            mva = mva[mask]
            w   = w[mask]

    return mva, w

def build_mass_edges(masses):
    """给定离散的质量中心，产生不等宽 y-edges。"""
    ms = np.array(sorted(set(float(m) for m in masses)))
    if len(ms) == 0:
        return None
    if len(ms) == 1:
        dm = 1.0
        return np.array([ms[0] - dm/2, ms[0] + dm/2])
    mids = (ms[1:] + ms[:-1]) / 2.0
    edges = np.concatenate(([ms[0] - (mids[0] - ms[0])], mids, [ms[-1] + (ms[-1] - mids[-1])]))
    return edges

# 新增：以 1 GeV 等宽分 bin，边界在整数 ±0.5，中心是整数
def build_uniform_mass_edges(masses, step=1.0):
    ms = [float(m) for m in masses]
    if len(ms) == 0:
        return None
    mmin = int(np.floor(min(ms)))
    mmax = int(np.ceil(max(ms)))
    start = mmin - 0.5
    stop = mmax + 0.5 + 1e-9  # 避免浮点边界遗漏
    return np.arange(start, stop, step)

def ensure_outdir():
    outdir = Path("../plots/interpolate_bkgYield")
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir

# === 新增：先画出 interploate_ma_list 的所有 bkg MVA_Score 分布 ===
def plot_bkg_mva_distributions():
    # 若关闭绘图，直接略过
    if not ENABLE_PLOTTING:
        return

    outdir = ensure_outdir() / "bkg_mva"
    outdir.mkdir(parents=True, exist_ok=True)

    # 三组背景与对应年份（与既有流程一致）
    groups = [
        (bkg_2022, years_22),
        (bkg_2023, years_23),
        (bkg_dyll, years_dyll),
    ]

    # 与后续 hist 设置一致的 binning
    bins = np.linspace(-0.1, 1.1, 241)

    for ma in interploate_ma_list:
        mass_tag = f"M{int(ma)}"
        plt.figure()
        plotted_any = False

        # 逐样本汇整跨年份的 MVA 与权重，并画在同一张图上
        for samples, years in groups:
            for sample in samples:
                mva_all = []
                w_all = []
                for y in years:
                    fpath = os.path.join(INPUT_BASE, sample, f"mA_M{int(ma)}", f"{y}.root")
                    if not os.path.exists(fpath):
                        continue
                    mva, w = read_arrays_with_mass(fpath, mass_tag)
                    if mva is None or len(mva) == 0:
                        continue
                    mva_all.append(mva)
                    w_all.append(w)
                if len(mva_all) == 0:
                    continue

                mva_all = np.concatenate(mva_all)
                w_all = np.concatenate(w_all)
                if mva_all.size == 0 or np.sum(w_all) <= 0:
                    continue

                plt.hist(
                    mva_all,
                    bins=bins,
                    weights=w_all,
                    histtype="step",
                    density=False,
                    label=f"{sample} (N={int(np.sum(w_all))})",
                    linewidth=1.5,
                )
                plotted_any = True

        if not plotted_any:
            plt.close()
            print(f"[plot] mA={ma}: no available inputs, skip.")
            continue

        plt.xlabel("MVA_Score")
        plt.ylabel("Normalized entries")
        plt.title(f"Background MVA distributions at mA={int(ma)}")
        plt.legend(ncol=2, fontsize=10)
        plt.grid(True, alpha=0.25)
        png = outdir / f"bkg_mva_ma{int(ma)}.png"
        pdf = outdir / f"bkg_mva_ma{int(ma)}.pdf"
        plt.savefig(png, bbox_inches="tight", dpi=150)
        plt.savefig(pdf, bbox_inches="tight")
        plt.close()
        print(f"[plot] saved: {png}")

# 以 ALP_Optimization.py 的做法：将 TH1F 转为 TGraph（只取 mva_low 以上）
def alp_hist2graph(hist, mva_low=0.01):
    nbins = hist.GetNbinsX()
    bin_low = hist.FindBin(mva_low)
    xs, ys = [], []
    xaxis = hist.GetXaxis()
    for ib in range(1, nbins + 1):
        if ib < bin_low:
            continue
        xs.append(xaxis.GetBinCenter(ib))
        ys.append(hist.GetBinContent(ib))
    # 若没有内容，给一个安全空图
    if len(xs) == 0:
        xs = [0.0]
        ys = [0.0]
    arrx = np.array(xs, dtype=np.float64)
    arry = np.array(ys, dtype=np.float64)
    return TGraph(len(arrx), arrx, arry)

# 以 ALP_Optimization.py 的 smooth 流程：SmoothSuper 并回填到新 TH1F
def alp_smooth(graph, hist, mva_low=0.01):
    h_smooth    = TH1F(hist.GetName()+"_smooth",    hist.GetTitle()+"_smooth",    hist.GetNbinsX(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    h_smooth_up = TH1F(hist.GetName()+"_smooth_up", hist.GetTitle()+"_smooth_up", hist.GetNbinsX(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    h_smooth_dn = TH1F(hist.GetName()+"_smooth_dn", hist.GetTitle()+"_smooth_dn", hist.GetNbinsX(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    h_smooth.Sumw2(); h_smooth_up.Sumw2(); h_smooth_dn.Sumw2()

    smoother = TGraphSmooth()
    g_smooth = smoother.SmoothSuper(graph)

    nbins = hist.GetNbinsX()
    bin_low = hist.FindBin(mva_low)
    x = array('d', [0.0])
    y = array('d', [0.0])

    for ib in range(1, nbins + 1):
        if ib < bin_low:
            # 低于 mva_low：沿用原始内容（对齐 ALP_Optimization 行为）
            y[0] = hist.GetBinContent(ib)
            in_range = True  # 不套用 x 范围门槛
        else:
            idx = ib - bin_low
            if idx < g_smooth.GetN():
                g_smooth.GetPoint(idx, x, y)
            else:
                y[0] = 0.0
                x[0] = 0.0
            if y[0] < 0.0:
                y[0] = 0.0
            # 加入与 ALP_Optimization 一致的范围门槛
            in_range = (y[0] >= 0.0 and x[0] <= 1.0 and x[0] >= 0.0)

        if in_range:
            h_smooth.SetBinContent(ib, y[0])
            yerr = np.sqrt(y[0]) if y[0] >= 0.0 else 0.0
            h_smooth_up.SetBinContent(ib, y[0] + yerr)
            h_smooth_dn.SetBinContent(ib, max(0.0, y[0] - yerr))
        else:
            # 超出范围：全部设 0（与 ALP_Optimization 一致）
            h_smooth.SetBinContent(ib, 0.0)
            h_smooth_up.SetBinContent(ib, 0.0)
            h_smooth_dn.SetBinContent(ib, 0.0)

    return h_smooth, h_smooth_up, h_smooth_dn

# 载入 cut 对照表，建立 mA->cut 的 dict，并提供线性内插
def load_cut_map(cut_json_path):
    with open(cut_json_path, "r") as f:
        data = json.load(f)
    pts = sorted([(float(e["mA"]), float(e["MVAcut"])) for e in data.get("results", [])], key=lambda x: x[0])
    return pts

def cut_for_mass(cut_pts, mass):
    m = float(mass)
    if not cut_pts:
        return 0.97  # 后备
    # 精确命中
    for mi, ci in cut_pts:
        if abs(mi - m) < 1e-9:
            return ci
    # 端点外：夹到最近端点
    if m <= cut_pts[0][0]:
        return cut_pts[0][1]
    if m >= cut_pts[-1][0]:
        return cut_pts[-1][1]
    # 区间线性内插
    for (m1, c1), (m2, c2) in zip(cut_pts[:-1], cut_pts[1:]):
        if m1 <= m <= m2:
            t = (m - m1) / (m2 - m1) if m2 > m1 else 0.0
            return c1 + t * (c2 - c1)
    return cut_pts[-1][1]

# 载入 cut 对照表为 dict，并提供「最近邻 current_ma_list」的 cut 查询
def load_cut_dict(cut_json_path):
    with open(cut_json_path, "r") as f:
        data = json.load(f)
    return {int(e["mA"]): float(e["MVAcut"]) for e in data.get("results", [])}

def nearest_cut_for_mass(cut_dict, mass, anchors):
    """
    对于内差点，取与 anchors(current_ma_list) 距离最近的 mA 的 cut 值；
    若 mass 本身就在 cut_dict，则直接回传该值。
    """
    m = int(round(float(mass)))
    if m in cut_dict:
        return cut_dict[m]
    # 最近邻
    nearest_anchor = min(anchors, key=lambda a: abs(a - m))
    return cut_dict.get(nearest_anchor, 0.97)

# 新增：严格小于且最近的 anchor；若无则回退到端点(anchors 最小值)
def lower_cut_for_mass(cut_dict, mass, anchors, default=0.97):
    m = float(mass)
    lowers = [a for a in anchors if a < m]
    anc = max(lowers) if lowers else min(anchors)
    return cut_dict.get(int(anc), default)

# 新增：严格大于且最近的 anchor；若无则回退到端点(anchors 最大值)
def upper_cut_for_mass(cut_dict, mass, anchors, default=0.97):
    m = float(mass)
    uppers = [a for a in anchors if a > m]
    anc = min(uppers) if uppers else max(anchors)
    return cut_dict.get(int(anc), default)

# 建立背景直方图：histos['mvaVal_M{ma}'][sample] = TH1F(...)
def build_interpolated_bkg_hists(mva_bins=240, mva_min=-0.1, mva_max=1.1, mva_branch_candidates=None, weight_candidates=None):
    # groups 定义：样本对应年份
    groups = [
        (bkg_2022, years_22),
        (bkg_2023, years_23),
        (bkg_dyll, years_dyll),
    ]
    histos = {}
    # 填充
    for ma in interploate_ma_list:
        mass_tag = f"M{int(ma)}"
        key = f"mvaVal_{mass_tag}"
        histos[key] = {}
        for samples, years in groups:
            for sample in samples:
                hname = f"{key}_{sample}"
                htitle = hname
                h = TH1F(hname, htitle, mva_bins, mva_min, mva_max)
                h.Sumw2()
                # 搜集各年资料
                filled = False
                for y in years:
                    fpath = os.path.join(INPUT_BASE, sample, f"mA_M{int(ma)}", f"{y}.root")
                    if not os.path.exists(fpath):
                        continue
                    # 以 mass_tag 解析对应的 MVA 分支（与 plot_variable_dataVmc.py 一致）
                    mva, w = read_arrays_with_mass(fpath, mass_tag)
                    if mva is None:
                        continue
                    # 逐项填入（含权重）
                    for val, wt in zip(mva, w):
                        try:
                            h.Fill(float(val), float(wt))
                            filled = True
                        except Exception:
                            continue
                if filled and (h.GetEntries() > 0 or h.Integral() > 0.0):
                    histos[key][sample] = h
    return histos

# 新增：将相同 binning 的 TH1F 字典加总为单一 TH1F
def sum_background_hists(hdict: dict, sum_name: str):
    """
    hdict: {sample: TH1F}
    回传: 加总后的 TH1F；若 hdict 为空则回传 None
    """
    if not hdict:
        return None
    # 取第一个样本的 binning 作为模板
    any_h = next(iter(hdict.values()))
    hsum = TH1F(sum_name, sum_name, any_h.GetNbinsX(), any_h.GetXaxis().GetXmin(), any_h.GetXaxis().GetXmax())
    hsum.Sumw2()
    for h in hdict.values():
        hsum.Add(h)
    return hsum

# 新增：将 cut 值先对齐到 signal bin 边界，再映射到背景直方图的确切低边 bin
def cut_value_to_low_bin(cut_val, h, sig_nbins=SIG_NBINS, sig_xmin=SIG_XMIN, sig_xmax=SIG_XMAX):
    """
    1) 以与 signal 完全相同的 binning (sig_nbins, sig_xmin, sig_xmax) 将 cut 对齐到最近的「边界」。
    2) 再用背景直方图 h 的 (nbins, xmin, xmax) 反解出对应的 bin 编号 low_bin。
       回传的是 Integral(low_bin, nbins) 的 low_bin。
    """
    # 对齐到 signal 的确切边界
    if cut_val <= sig_xmin:
        x_edge = sig_xmin
    elif cut_val >= sig_xmax:
        x_edge = sig_xmax
    else:
        sig_bw = (sig_xmax - sig_xmin) / float(sig_nbins)
        # 找到最近的 signal 边界（避免浮点误差）
        i_sig_edge = int(round((cut_val - sig_xmin) / sig_bw))
        x_edge = sig_xmin + i_sig_edge * sig_bw

    # 映射到背景直方图的 bin 编号
    nb = h.GetNbinsX()
    xmin = h.GetXaxis().GetXmin()
    xmax = h.GetXaxis().GetXmax()
    bw = (xmax - xmin) / float(nb)

    # 计算对应的「边界索引」后转成 ROOT bin 编号（1 起算）
    i_edge = int(round((x_edge - xmin) / bw))
    low_bin = i_edge + 1  # 边界右侧的第一个 bin

    # 边界情形处理
    if low_bin < 1:
        low_bin = 1
    return low_bin

# 根据平滑后的直方图，计算通过 cut 的加总数量
def integral_above_cut(h, cut_val):
    nb = h.GetNbinsX()
    # ROOT FindBin: 若 cut 剛好落在邊界，會回傳下一個 bin
    low_bin = h.FindBin(cut_val + 1e-6)
    if low_bin > nb:
        return 0.0
    return float(h.Integral(low_bin, nb))

# ==== Collectors ====
def collect_background_xyw():
    xs, ys, ws = [], [], []
    # 三类背景各自对应年份
    groups = [
        (bkg_2022, years_22),
        (bkg_2023, years_23),
        (bkg_dyll, years_dyll),
    ]
    missing = 0
    for samples, years in groups:
        for s in samples:
            for y in years:
                for ma in ma_list:
                    fpath = os.path.join(INPUT_BASE, s, f"mA_M{ma}", f"{y}.root")
                    if not os.path.exists(fpath):
                        missing += 1
                        continue
                    mva, w = read_arrays(fpath)
                    if mva is None:
                        continue
                    xs.append(mva)
                    ws.append(w)
                    ys.append(np.full_like(mva, float(ma), dtype=np.float64))
    if missing > 0:
        print(f"[bkg] missing files skipped: {missing}")
    if len(xs) == 0:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(ws)

def collect_signal_xyw():
    xs, ys, ws = [], [], []
    missing = 0
    for s in sig_samples:
        ma = mass_from_sig(s)
        for y in years_sig:
            fpath = os.path.join(INPUT_BASE, s, f"{y}.root")
            if not os.path.exists(fpath):
                missing += 1
                continue
            mva, w = read_arrays(fpath)
            if mva is None:
                continue
            xs.append(mva)
            ws.append(w)
            ys.append(np.full_like(mva, float(ma), dtype=np.float64))
    if missing > 0:
        print(f"[sig] missing files skipped: {missing}")
    if len(xs) == 0:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(ws)

# ==== Main ====
def main():
    outdir = ensure_outdir()

    # 先画出各背景在欲内插质量点的 MVA 分布（受开关控制）
    if ENABLE_PLOTTING:
        plot_bkg_mva_distributions()

    # y 轴（BDT）bin：0~1 共 100 bin；如需 -1~1 可改成 np.linspace(-1,1,101)
    bdt_edges = np.linspace(0.0, 1.0, 101)

    # ==== 新增：为 interploate_ma_list 产生背景直方图、平滑并依优化 cut 计算通过数量 ====
    print("[interp] building background histograms for interpolated masses...")
    histos = build_interpolated_bkg_hists()
    cut_dict = load_cut_dict(optimized_BDT_Cut)

    results = []
    for ma in interploate_ma_list:
        mass_tag = f"M{int(ma)}"
        key = f"mvaVal_{mass_tag}"
        # 三种 cut 选择策略
        cut_near = nearest_cut_for_mass(cut_dict, float(ma), current_ma_list)
        cut_low  = lower_cut_for_mass(cut_dict, float(ma), current_ma_list)
        cut_up   = upper_cut_for_mass(cut_dict, float(ma), current_ma_list)

        entry = {
            "mA": int(ma),
            "MVAcut_nearest": round(float(cut_near), 3),
            "MVAcut_lower": round(float(cut_low), 3),
            "MVAcut_upper": round(float(cut_up), 3)
            # "bkgYields": {}
        }

        # 取得此 mA 的所有样本直方图
        hdict = histos.get(key, {})

        # ---- 变更处：所有背景样本先加总 -> 平滑 -> 积分 ----
        h_sum_all = sum_background_hists(hdict, sum_name=f"{key}_ALLBKG_sum")
        if h_sum_all is not None and (h_sum_all.GetEntries() > 0 or h_sum_all.Integral() > 0.0):
            g_sum = alp_hist2graph(h_sum_all, mva_low=0.01)
            h_sum_smooth, _, _ = alp_smooth(g_sum, h_sum_all, mva_low=0.01)

            total_near = integral_above_cut(h_sum_smooth, cut_near)
            total_low  = integral_above_cut(h_sum_smooth, cut_low)
            total_up   = integral_above_cut(h_sum_smooth, cut_up)
        else:
            total_near = total_low = total_up = 0.0

        # ---- 参考：仍保留每个样本各自平滑后积分（使用最近邻 cut）----
        for sample, h in hdict.items():
            g = alp_hist2graph(h, mva_low=0.01)
            h_smooth, _, _ = alp_smooth(g, h, mva_low=0.01)
            val_near = integral_above_cut(h_smooth, cut_near)
            # entry["bkgYields"][sample] = round(val_near, 3)

        # 总量栏位：以「所有背景先加总后平滑」的结果为准
        # entry["bkgTotal"] = round(total_near, 3)
        entry["bkgTotal_nearest"] = round(total_near, 3)
        entry["bkgTotal_lower"]   = round(total_low, 3)
        entry["bkgTotal_upper"]   = round(total_up, 3)

        results.append(entry)

    out_json = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/bkg_interpolated_yields_run3.json"
    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    with open(out_json, "w") as f:
        json.dump({
            "optimized_BDT_Cut": optimized_BDT_Cut,
            "hist_def": {"nbins": 240, "xmin": -0.1, "xmax": 1.1, "mva_low": 0.01},
            "results": results
        }, f, indent=2)
    print(f"[interp] saved yields -> {out_json}")

if __name__ == "__main__":
    main()

