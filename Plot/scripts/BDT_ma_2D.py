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

# 你之前实际使用的组合：
years_sig  = ["2022preEE"]  # 信号
years_22   = ["2022preEE", "2022postEE"]  # 背景（2022）
years_23   = ["2023preBPix","2023postBPix"]
years_dyll = ["2022preEE","2022postEE","2023preBPix","2023postBPix"]

# 扫描的 ma 列表（背景也按你的原脚本扫）
interpolate = True

if not interpolate:
    ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
else:
    ma_list = list(range(1, 31))

tuneStyle = False # False # True: 只產生樣式預覽圖，不讀資料

# sig_ma_ticks =  [1,2,3,4,5,6,7,8,9,10,15,20,25,30]  # x 軸只顯示這些質量刻度
sig_ma_ticks =  [5,15,30]  # x 軸只顯示這些質量刻度
# bkg_ma_ticks = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
bkg_ma_ticks = [5,15,30]
# ===== Samples =====
sig_samples = ["ALP_M5", "ALP_M15", "ALP_M30"]

# 背景
bkg_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
bkg_2023 = ["DYGto2LG_10to100"]
bkg_dyll = ["DYJetsToLL"]
# 如需其它：["ZGToLLG","WGToLNuG","ZG2JToG2L2J","EWKZ2J","TT","TTGJets","TGJets","ttWJets","ttZJets","WW","WZ","ZZ"]

MVA_CANDIDATES = ["MVA_Score"]
WEIGHT_CANDIDATES = ["weight", "event_weight", "wgt_nominal", "genWeight", "genweight", "totalWeight", "w"]

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
    outdir = Path("../plots/BDT_ma_2D")
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
    fig.text(0.84, 0.965, r"$61.89\,\mathrm{fb}^{-1}\ (13.6\ \mathrm{TeV})$", ha="right", va="top", fontsize=16)
    plt.subplots_adjust(left=0.13, right=0.98, bottom=0.14, top=0.92)

    plt.savefig(out_png, dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.05)
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
    print(f"[style] saved -> {outdir/'BDT_vs_ma_background_style.png'}")

    # 訊號樣式預覽
    sig_masses = [mass_from_sig(s) for s in sig_samples]
    sig_x_edges = build_uniform_mass_edges(sig_masses, step=1.0)
    draw_hist2d(
        x=None, y=None, w=None,
        x_edges=sig_x_edges, y_edges=bdt_edges,
        title=None,
        out_png=outdir / "BDT_vs_ma_signal_style.png",
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
    # 三類背景各自對應年份
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
                    fpath = os.path.join(INPUT_BASE, s, f"ALP_M{ma}", f"{y}.root")
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
        # extract ma from sig_samples
        # sig_masses = [mass_from_sig(s) for s in sig_samples]
        # sig_x_edges = build_uniform_mass_edges(sig_masses, step=1.0)
        # -------------------------------
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

