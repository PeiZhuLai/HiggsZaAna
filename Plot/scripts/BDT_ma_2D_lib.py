# ==== New imports/config for plotting ====
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot

INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/"
# 新增：merged BDT（同一 era 檔內含多個 mA 分支）
INPUT_BASE_MERGED = "/eos/home-p/pelai/HZa/root_P2Root/run3_mergedBDT/"

# Year of Signal
year_sig_2022 = ["2022preEE", "2022postEE"]
year_sig_2023 = ["2023preBPix", "2023postBPix"]
year_sig_2024 = ["2024"]
# Year of Bkg
year_DYG_2022 = ["2022preEE", "2022postEE"]
year_DYG_2023 = ["2023preBPix", "2023postBPix"]
year_DYG_2024 = ["2024"]
years_DYJet_2022 = ["2022preEE","2022postEE"]
year_DYJet_2023  = ["2023preBPix", "2023postBPix"]
years_DYJet_2024 = ["2024"]
# Year of Data
years_Data_2022 = ["2022preEE","2022postEE"]
year_Data_2023  = ["2023preBPix", "2023postBPix"]
years_Data_2024 = ["2024"]

# Name of Signal Sample
name_sig_2022 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
name_sig_2023 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
name_sig_2024 = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
# Name of Bkg Sample
name_DYG_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
name_DYG_2023 = ["DYGto2LG_10to100"]
name_DYG_2024 = ["DYGto2LG_10to100"]
name_DYJet_2022 = ["DYJetsToLL"]
name_DYJet_2023 = ["DYJetsToLL"]
name_DYJet_2024 = ["DYJetsTo2E","DYJetsTo2Mu","DYJetsTo2Tau"]
# Name of Data Sample
name_Data_2022 = ["Data"]
name_Data_2023 = ["Data"]
name_Data_2024 = ["Data"]

#-----------------------------------------------------
# 使用你提供的選項來組裝實際要跑的 samples/years
years_sig = year_sig_2022 + year_sig_2023 + year_sig_2024
sig_samples = name_sig_2022 + name_sig_2023 + name_sig_2024

# 改：背景要把 DYG + DYJets 全部跑完（各自對應年份）
bkg_samples_by_year = {
    "2022preEE":  name_DYG_2022 + name_DYJet_2022,
    "2022postEE": name_DYG_2022 + name_DYJet_2022,
    "2023preBPix":  name_DYG_2023 + name_DYJet_2023,
    "2023postBPix": name_DYG_2023 + name_DYJet_2023,
    "2024": name_DYG_2024 + name_DYJet_2024,
}
bkg_years_all = year_DYG_2022 + year_DYG_2023 + year_DYG_2024
#-----------------------------------------------------

#-----------------------------------------------------
# 扫描的 ma 列表（背景也按你的原脚本扫）
interpolate = True

if not interpolate:
    ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
else:
    ma_list = list(range(1, 31))

tuneStyle = False # False # True: 只產生樣式預覽圖，不讀資料

sig_ma_ticks =  [1,2,3,4,5,6,7,8,9,10,15,20,25,30]  # x 軸只顯示這些質量刻度
# sig_ma_ticks =  [5,15,30]  # x 軸只顯示這些質量刻度
bkg_ma_ticks = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
# bkg_ma_ticks = [5,15,30]

MVA_CANDIDATES = ["MVA_Score"]
WEIGHT_CANDIDATES = ["weight", "w"]
# 新增：背景 merged 檔內的 MVA 分支前綴
BKG_MVA_BRANCH_PREFIXES = ["MVA_Score_mA_M", "MVA_Score_ma_M"]
#-----------------------------------------------------

# ==== Helpers ====
# 设置全局字体大小和图形样式
plt.rcParams.update({
    'figure.figsize': (10, 6),  # 图形大小
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
    # 支援：ALP_M5 / mA_M5 -> 5
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
    # 找第一個 TTree
    for k, cls in f.classnames().items():
        if "TTree" in cls:
            return f[k]
    # 後備名稱
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
    # 轉為 1D
    try:
        mva = np.asarray(mva).astype(np.float64)
        w = np.asarray(w).astype(np.float64)
    except Exception:
        return None, None
    # 對齊長度
    n = min(len(mva), len(w))
    return mva[:n], w[:n]

# 新增：從同一個 tree 讀指定 ma 的 MVA 分支（背景 merged 檔用）
def read_mva_for_ma_from_tree(tree, ma: int):
    """
    Return mva numpy array for given ma from a merged-background tree.
    Tries exact name and common uproot-suffixed variants like '.0'.
    """
    if tree is None:
        return None

    keys = set(tree.keys())
    # 候選分支名（含常見的 uproot 版本尾碼）
    cand = []
    for pfx in BKG_MVA_BRANCH_PREFIXES:
        base = f"{pfx}{int(ma)}"
        cand.extend([base, base + ".0", base + ".1", base + "_0", base + "_1"])

    for br in cand:
        if br in keys:
            try:
                arr = tree[br].array(library="np")
                return np.asarray(arr, dtype=np.float64)
            except Exception:
                return None
    return None

def build_mass_edges(masses):
    """給定離散的質量中心，產生不等寬 y-edges。"""
    ms = np.array(sorted(set(float(m) for m in masses)))
    if len(ms) == 0:
        return None
    if len(ms) == 1:
        dm = 1.0
        return np.array([ms[0] - dm/2, ms[0] + dm/2])
    mids = (ms[1:] + ms[:-1]) / 2.0
    edges = np.concatenate(([ms[0] - (mids[0] - ms[0])], mids, [ms[-1] + (ms[-1] - mids[-1])]))
    return edges

# 新增：以 1 GeV 等寬分 bin，邊界在整數 ±0.5，中心是整數
def build_uniform_mass_edges(masses, step=1.0):
    ms = [float(m) for m in masses]
    if len(ms) == 0:
        return None
    mmin = int(np.floor(min(ms)))
    mmax = int(np.ceil(max(ms)))
    start = mmin - 0.5
    stop = mmax + 0.5 + 1e-9  # 避免浮點邊界遺漏
    return np.arange(start, stop, step)

# 新增：加權皮爾森相關係數
def weighted_corr(x, y, w):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    sw = w.sum()
    if sw <= 0:
        return np.nan
    mx = np.sum(w * x) / sw
    my = np.sum(w * y) / sw
    dx = x - mx
    dy = y - my
    cov = np.sum(w * dx * dy) / sw
    vx = np.sum(w * dx * dx) / sw
    vy = np.sum(w * dy * dy) / sw
    if vx <= 0 or vy <= 0:
        return np.nan
    return cov / np.sqrt(vx * vy)

def ensure_outdir():
    outdir = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/BDT_ma_2D")
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir

def draw_hist2d(
    x, y, w, x_edges, y_edges, title, out_png,
    y_label="BDT Score", x_tick_masses=None,
    corr_text=None, corr_loc="upper right", corr_pos=None,
    xtick_labelsize=None, xlabel_size=None, style_only=False
):
    # 新增：style_only=True 時不根據輸入 (x,y,w) 計數，直接用假資料生成 2D 畫布
    if not style_only:
        H, xe, ye = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=w)
    else:
        xe, ye = x_edges, y_edges
        nx = len(xe) - 1
        ny = len(ye) - 1
        # 產生正值漸層 (避免 LogNorm 問題)
        gx = np.linspace(0.3, 1.0, max(nx, 1))
        gy = np.linspace(0.3, 1.0, max(ny, 1))
        H = np.outer(gx[:nx], gy[:ny])
        if H.size == 0:
            # 邊界不足時用最小 1x1
            H = np.array([[1.0]])

    # 更健壯的 LogNorm 邊界處理，避免 vmin == vmax 或沒有正數時崩潰
    if np.any(H > 0):
        vmin = max(np.min(H[H > 0]), 1e-1)
        vmax = float(np.max(H))
    else:
        vmin, vmax = 1e-1, 1.0
    if not (vmax > vmin):
        vmax = vmin * 10.0

    plt.figure(figsize=(8, 6))
    pcm = plt.pcolormesh(
        xe, ye, H.T,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap="viridis", shading="auto"
    )
    plt.colorbar(pcm, label="Events", pad=0.01)
    plt.xlabel(r"$m_{a}$ [GeV]")
    if xlabel_size is not None:
        plt.gca().xaxis.label.set_size(xlabel_size)
    plt.ylabel(y_label)

    # 只顯示指定的質量刻度（自動過濾超出邊界者）
    if x_tick_masses is not None and len(x_tick_masses) > 0:
        xmin_c = x_edges[0] + 0.5
        xmax_c = x_edges[-1] - 0.5
        ticks = [int(m) for m in x_tick_masses if (xmin_c <= float(m) <= xmax_c)]
        plt.xticks(ticks, [str(m) for m in ticks], rotation=0)
    else:
        xcenters = 0.5 * (x_edges[:-1] + x_edges[1:])
        plt.xticks(xcenters, [f"{int(round(c))}" for c in xcenters], rotation=0)

    if xtick_labelsize is not None:
        plt.gca().tick_params(axis="x", labelsize=xtick_labelsize)

    if title:
        plt.title(title)
    plt.tight_layout()

    # 圖框內 correlation 標示
    if corr_text:
        ax = plt.gca()
        if corr_pos is not None:
            # 新增：corr_pos 使用 axes fraction (0~1) 直接指定文字位置，會覆蓋 corr_loc
            x_pos, y_pos = corr_pos
            ha, va = "left", "bottom"
        else:
            loc = (corr_loc or "upper right").lower()
            x_pos, y_pos = 0.98, 0.98
            ha, va = "right", "top"
            if "lower" in loc:
                y_pos = 0.02
                va = "bottom"
            if "left" in loc:
                x_pos = 0.02
                ha = "left"
            elif "center" in loc:
                x_pos = 0.5
                ha = "center"
        ax.text(
            x_pos, y_pos, corr_text,
            ha=ha, va=va, transform=ax.transAxes,
            fontsize=16,
            fontweight='bold',
            bbox=dict(facecolor="white", edgecolor="white", alpha=0.75)
        )

    # 左上角 CMS + Preliminary，右上角亮度與能量
    fig = plt.gcf()
    x0, y0 = 0.13, 0.97
    t_cms = fig.text(x0, y0, "CMS", ha="left", va="top", fontsize=19, fontweight="bold")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    cms_bb = t_cms.get_window_extent(renderer=renderer)
    dx = cms_bb.width / fig.bbox.width + 0.006
    fig.text(x0 + dx, y0, "Preliminary", ha="left", va="top", fontsize=19)
    fig.text(0.84, 0.965, r"$170.84\,\mathrm{fb}^{-1}\ (13.6\ \mathrm{TeV})$", ha="right", va="top", fontsize=16)
    plt.subplots_adjust(left=0.13, right=0.98, bottom=0.14, top=0.92)

    # 只輸出 PDF（不輸出 PNG）
    out_pdf = out_png.with_suffix(".pdf")
    plt.savefig(out_pdf, bbox_inches="tight", pad_inches=0.05)
    plt.close()

# 新增：樣式預覽（不讀取資料）
def make_style_previews(outdir):
    # y 軸（BDT）bin：0~1 共 100 bin
    bdt_edges = np.linspace(0.0, 1.0, 101)

    # 背景樣式預覽
    bkg_x_edges = build_uniform_mass_edges(ma_list, step=1.0)
    draw_hist2d(
        x=None, y=None, w=None,
        x_edges=bkg_x_edges, y_edges=bdt_edges,
        title=None,
        out_png=outdir / "BDT_vs_ma_background_style.png",
        y_label="Background BDT Score",
        x_tick_masses=ma_list,
        corr_text="Correlation = 0.000",
        corr_loc="upper right",
        xtick_labelsize=(10 if interpolate else None),
        style_only=True
    )
    print(f"[style] saved -> {outdir/'BDT_vs_ma_background_style.pdf'}")

    # 訊號樣式預覽
    sig_masses = [mass_from_sig(s) for s in sig_samples]
    sig_x_edges = build_uniform_mass_edges(sig_masses, step=1.0)
    draw_hist2d(
        x=None, y=None, w=None,
        x_edges=sig_x_edges, y_edges=bdt_edges,
        title=None,
        out_png=outdir / "BDT_vs_ma_signal_style.pdf",
        y_label="Signal BDT Score",
        x_tick_masses=ma_list,
        corr_text="Correlation = 0.000",
        corr_loc="lower right",
        xtick_labelsize=(10 if interpolate else None),
        style_only=True
    )
    print(f"[style] saved -> {outdir/'BDT_vs_ma_signal_style.png'}")

# ==== Collectors ====
def collect_background_xyw():
    xs, ys, ws = [], [], []
    missing = 0

    # 改：跑完所有背景（DYG + DYJets），每個檔內含 MVA_Score_mA_M*
    for y in bkg_years_all:
        for s in bkg_samples_by_year.get(y, []):
            fpath = os.path.join(INPUT_BASE_MERGED, s, f"{y}.root")
            if not os.path.exists(fpath):
                missing += 1
                continue

            tree = first_tree(fpath)
            if tree is None:
                continue

            # weights（同檔共用）
            w_all = None
            for wv in WEIGHT_CANDIDATES:
                if wv in tree.keys():
                    try:
                        w_all = np.asarray(tree[wv].array(library="np"), dtype=np.float64)
                        break
                    except Exception:
                        w_all = None

            for ma in ma_list:
                mva = read_mva_for_ma_from_tree(tree, ma)
                if mva is None:
                    continue

                if w_all is None:
                    w = np.ones_like(mva, dtype=np.float64)
                else:
                    n = min(len(mva), len(w_all))
                    mva = mva[:n]
                    w = w_all[:n]

                xs.append(mva)  # x = BDT
                ws.append(w)
                ys.append(np.full_like(mva, float(ma), dtype=np.float64))  # y = mA

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

    # 新增：樣式調整模式，完全不讀資料，直接產出預覽圖後結束
    if tuneStyle:
        make_style_previews(outdir)
        return

    # y 軸（BDT）bin：0~1 共 100 bin；如需 -1~1 可改成 np.linspace(-1,1,101)
    bdt_edges = np.linspace(0.0, 1.0, 101)

    # ---- 背景：x=ma, y=BDT ----
    bx, by, bw = collect_background_xyw()  # 目前 bx=mva, by=ma
    if bx.size > 0:
        bkg_x_edges = build_uniform_mass_edges(ma_list, step=1.0)
        corr_bkg = weighted_corr(by, bx, bw)
        draw_hist2d(
            by, bx, bw, bkg_x_edges, bdt_edges,
            title=None,
            out_png=outdir / "BDT_vs_ma_background.png",
            y_label="Background BDT Score",
            x_tick_masses=bkg_ma_ticks,
            corr_text=f"Correlation = {corr_bkg:.3f}",
            corr_pos=(0.53, 0.5),  
            xtick_labelsize=(14 if interpolate else None)  # interpolate=True 時縮小 x 刻度字體
            # 如需同時調整 x 標籤大小可加：xlabel_size=(16 if interpolate else None)
        )
        print(f"[bkg] saved -> {outdir/'BDT_vs_ma_background.png'} (corr={corr_bkg:.3f})")
    else:
        print("[bkg] no entries found.")

    # ---- 訊號：x=ma, y=BDT ----
    sx, sy, sw = collect_signal_xyw()  # 目前 sx=mva, sy=ma
    if sx.size > 0:
        # extract ma from ma_list
        sig_x_edges = build_uniform_mass_edges(ma_list, step=1.0)
        corr_sig = weighted_corr(sy, sx, sw)
        draw_hist2d(
            sy, sx, sw, sig_x_edges, bdt_edges,
            title=None,
            out_png=outdir / "BDT_vs_ma_signal.png",
            y_label="Signal BDT Score",
            x_tick_masses=sig_ma_ticks,  # 只顯示你提供的質量清單，超出邊界會自動過濾
            corr_text=f"Correlation = {corr_sig:.3f}",
            corr_pos=(0.53, 0.5),
            xtick_labelsize=(14 if interpolate else None)
            # 如需同時調整 x 標籤大小可加：xlabel_size=(16 if interpolate else None)
        )
        print(f"[sig] saved -> {outdir/'BDT_vs_ma_signal.png'} (corr={corr_sig:.3f})")
    else:
        print("[sig] no entries found.")

if __name__ == "__main__":
    main()

