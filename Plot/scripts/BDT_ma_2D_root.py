# ==== New imports/config for plotting ====
import os
from pathlib import Path
import numpy as np
import uproot

# NEW: PyROOT plotting (replace matplotlib)
import ROOT

INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/"

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

bkg_2022 = name_DYG_2022
bkg_2023 = name_DYG_2023
years_22 = year_DYG_2022
years_23 = year_DYG_2023

bkg_dyll = name_DYJet_2022 + name_DYJet_2023 + name_DYJet_2024
years_dyll = years_DYJet_2022 + year_DYJet_2023 + years_DYJet_2024
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
#-----------------------------------------------------

# ==== Helpers ====
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

# NEW: ROOT global style (roughly match your matplotlib look)
def setup_root_style():
    ROOT.gROOT.SetBatch(True)
    st = ROOT.gStyle
    st.SetOptStat(0)
    st.SetTitle(0)

    st.SetPadLeftMargin(0.13)
    st.SetPadRightMargin(0.15)   # leave room for colorbar like matplotlib
    st.SetPadBottomMargin(0.14)
    st.SetPadTopMargin(0.08)

    st.SetLabelFont(42, "XYZ")
    st.SetTitleFont(42, "XYZ")
    st.SetTitleSize(0.055, "X")
    st.SetTitleSize(0.055, "Y")
    st.SetLabelSize(0.045, "X")
    st.SetLabelSize(0.045, "Y")
    st.SetLabelSize(0.040, "Z")

    st.SetNumberContours(255)
    # similar to viridis-ish palette (fallback if unavailable)
    try:
        ROOT.gStyle.SetPalette(ROOT.kViridis)
    except Exception:
        st.SetPalette(ROOT.kBird)

def _edges_to_roots(edges):
    arr = np.asarray(edges, dtype=np.float64)
    return (len(arr) - 1, arr)

def _safe_logz_range(h2):
    # emulate your vmin/vmax logic
    zmin = None
    zmax = h2.GetMaximum()
    if zmax <= 0:
        return 1e-1, 1.0
    # find smallest positive
    zmin = None
    nbx, nby = h2.GetNbinsX(), h2.GetNbinsY()
    for ix in range(1, nbx + 1):
        for iy in range(1, nby + 1):
            v = h2.GetBinContent(ix, iy)
            if v > 0 and (zmin is None or v < zmin):
                zmin = v
    if zmin is None:
        zmin = 1e-1
    zmin = max(float(zmin), 1e-1)
    zmax = float(zmax)
    if not (zmax > zmin):
        zmax = zmin * 10.0
    return zmin, zmax

def draw_hist2d_root(
    x, y, w, x_edges, y_edges, title, out_pdf,
    y_label="BDT Score", x_tick_masses=None,
    corr_text=None, corr_loc="upper right", corr_pos=None,
    xtick_labelsize=None, xlabel_size=None, style_only=False
):
    setup_root_style()

    # canvas
    c = ROOT.TCanvas("c", "c", 800, 600)
    c.SetLogz(True)

    nbx, xarr = _edges_to_roots(x_edges)
    nby, yarr = _edges_to_roots(y_edges)

    h = ROOT.TH2D("h2", "", nbx, xarr, nby, yarr)
    h.Sumw2()

    if style_only:
        # fill a positive gradient to make logz stable
        for ix in range(1, nbx + 1):
            for iy in range(1, nby + 1):
                gx = 0.3 + 0.7 * (ix / max(nbx, 1))
                gy = 0.3 + 0.7 * (iy / max(nby, 1))
                h.SetBinContent(ix, iy, gx * gy)
    else:
        # fill from arrays
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if w is None:
            w = np.ones_like(x, dtype=np.float64)
        else:
            w = np.asarray(w, dtype=np.float64)
        n = min(len(x), len(y), len(w))
        for i in range(n):
            h.Fill(float(x[i]), float(y[i]), float(w[i]))

    # axis titles
    h.GetXaxis().SetTitle("m_{a} [GeV]")
    h.GetYaxis().SetTitle(y_label)
    if xlabel_size is not None:
        # ROOT uses relative sizes; map your "pt-size" loosely
        h.GetXaxis().SetTitleSize(0.055 if xlabel_size is None else 0.060)

    # x ticks: show only selected masses
    if x_tick_masses is not None and len(x_tick_masses) > 0:
        # Use labels by bin center with custom labels (TGaxis-ish workaround):
        # keep numeric axis, but hide most labels by setting divisions and then relabel bins.
        h.GetXaxis().SetNdivisions(510)
        # Clear all bin labels then set selected bins
        for ix in range(1, nbx + 1):
            h.GetXaxis().SetBinLabel(ix, "")
        xmin_c = x_edges[0] + 0.5
        xmax_c = x_edges[-1] - 0.5
        for m in x_tick_masses:
            fm = float(m)
            if not (xmin_c <= fm <= xmax_c):
                continue
            ib = h.GetXaxis().FindBin(fm)
            h.GetXaxis().SetBinLabel(ib, str(int(fm)))
        h.LabelsOption("h", "X")

    if xtick_labelsize is not None:
        # ROOT label size is fraction of pad; approximate
        h.GetXaxis().SetLabelSize(0.040 if xtick_labelsize >= 14 else 0.032)

    # z-range for log
    zmin, zmax = _safe_logz_range(h)
    h.SetMinimum(zmin)
    h.SetMaximum(zmax)

    # draw
    h.Draw("COLZ")

    # correlation box
    if corr_text:
        latex = ROOT.TLatex()
        latex.SetNDC(True)
        latex.SetTextFont(42)
        latex.SetTextSize(0.045)
        latex.SetTextColor(ROOT.kBlack)

        # background box (white)
        x_pos, y_pos = 0.98, 0.98
        align = 33  # right-top
        if corr_pos is not None:
            x_pos, y_pos = float(corr_pos[0]), float(corr_pos[1])
            align = 11  # left-bottom (match your matplotlib default for manual pos)
        else:
            loc = (corr_loc or "upper right").lower()
            if "lower" in loc:
                y_pos = 0.18
                align = 13  # right-bottom
            if "left" in loc:
                x_pos = 0.20
                align = 13 if "lower" in loc else 11
            elif "center" in loc:
                x_pos = 0.50
                align = 23

        # Draw "bbox" via TPaveText
        pave = ROOT.TPaveText(x_pos - 0.22, y_pos - 0.06, x_pos, y_pos, "NDC")
        pave.SetFillColor(ROOT.kWhite)
        pave.SetFillStyle(1001)
        pave.SetLineColor(ROOT.kWhite)
        pave.SetTextFont(42)
        pave.SetTextSize(0.045)
        pave.SetTextAlign(12)
        pave.AddText(corr_text)
        pave.Draw("same")

    # CMS + Preliminary (top-left), lumi (top-right)
    cms = ROOT.TLatex()
    cms.SetNDC(True)
    cms.SetTextFont(62)
    cms.SetTextSize(0.055)
    cms.DrawLatex(0.13, 0.965, "CMS")

    prelim = ROOT.TLatex()
    prelim.SetNDC(True)
    prelim.SetTextFont(42)
    prelim.SetTextSize(0.050)
    prelim.DrawLatex(0.205, 0.965, "Preliminary")

    lumi = ROOT.TLatex()
    lumi.SetNDC(True)
    lumi.SetTextFont(42)
    lumi.SetTextSize(0.045)
    lumi.SetTextAlign(31)
    lumi.DrawLatex(0.98, 0.965, "61.89 fb^{-1} (13.6 TeV)")

    c.Update()
    c.SaveAs(str(out_pdf))
    c.Close()

# 新增：樣式預覽（不讀取資料）-> ROOT 版
def make_style_previews(outdir):
    bdt_edges = np.linspace(0.0, 1.0, 101)

    bkg_x_edges = build_uniform_mass_edges(ma_list, step=1.0)
    draw_hist2d_root(
        x=None, y=None, w=None,
        x_edges=bkg_x_edges, y_edges=bdt_edges,
        title=None,
        out_pdf=outdir / "BDT_vs_ma_background_style.pdf",
        y_label="Background BDT Score",
        x_tick_masses=ma_list,
        corr_text="Correlation = 0.000",
        corr_loc="upper right",
        xtick_labelsize=(10 if interpolate else None),
        style_only=True
    )

    sig_masses = [mass_from_sig(s) for s in sig_samples]
    sig_x_edges = build_uniform_mass_edges(sig_masses, step=1.0)
    draw_hist2d_root(
        x=None, y=None, w=None,
        x_edges=sig_x_edges, y_edges=bdt_edges,
        title=None,
        out_pdf=outdir / "BDT_vs_ma_signal_style.pdf",
        y_label="Signal BDT Score",
        x_tick_masses=ma_list,
        corr_text="Correlation = 0.000",
        corr_loc="lower right",
        xtick_labelsize=(10 if interpolate else None),
        style_only=True
    )

# ==== Collectors ====
def collect_background_xyw():
    xs, ys, ws = [], [], []
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
        draw_hist2d_root(
            by, bx, bw, bkg_x_edges, bdt_edges,
            title=None,
            out_pdf=outdir / "BDT_vs_ma_background.pdf",
            y_label="Background BDT Score",
            x_tick_masses=bkg_ma_ticks,
            corr_text=f"Correlation = {corr_bkg:.3f}",
            corr_pos=(0.53, 0.5),
            xtick_labelsize=(14 if interpolate else None)
        )
        print(f"[bkg] saved -> {outdir/'BDT_vs_ma_background.pdf'} (corr={corr_bkg:.3f})")
    else:
        print("[bkg] no entries found.")

    # ---- 訊號：x=ma, y=BDT ----
    sx, sy, sw = collect_signal_xyw()  # 目前 sx=mva, sy=ma
    if sx.size > 0:
        sig_x_edges = build_uniform_mass_edges(ma_list, step=1.0)
        corr_sig = weighted_corr(sy, sx, sw)
        draw_hist2d_root(
            sy, sx, sw, sig_x_edges, bdt_edges,
            title=None,
            out_pdf=outdir / "BDT_vs_ma_signal.pdf",
            y_label="Signal BDT Score",
            x_tick_masses=sig_ma_ticks,  # 只顯示你提供的質量清單，超出邊界會自動過濾
            corr_text=f"Correlation = {corr_sig:.3f}",
            corr_pos=(0.53, 0.5),
            xtick_labelsize=(14 if interpolate else None)
        )
        print(f"[sig] saved -> {outdir/'BDT_vs_ma_signal.pdf'} (corr={corr_sig:.3f})")
    else:
        print("[sig] no entries found.")

if __name__ == "__main__":
    main()

