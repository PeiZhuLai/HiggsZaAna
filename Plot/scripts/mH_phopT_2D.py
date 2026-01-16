# ==== Imports/config for mH vs mZ plotting only ====
import os
from pathlib import Path
from typing import Optional

import numpy as np
import uproot

# --- PyROOT ---
import ROOT

INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_sumStudy/"
treeName = "inclusive"

# Years
years_sig = ["2023preBPix"]
bkg_years_all = ["2023preBPix"]

# Samples
sig_samples = ["mA_M1","mA_M2","mA_M3","mA_M4","mA_M5","mA_M6","mA_M7","mA_M8","mA_M9","mA_M10", "mA_M15", "mA_M20", "mA_M25", "mA_M30"]
bkg_samples_by_year = {
    "2023preBPix": ["DYGto2LG_10to100", "DYJetsToLL"],
}

WEIGHT_CANDIDATES = ["weight", "w"]
H_MASS_CANDIDATES = ["H_mass", "Hmass", "mH", "Higgs_mass", "massH"]
Z_MASS_CANDIDATES = ["Z_mass", "Zmass", "mZ", "Zcand_mass", "massZ"]

# NEW: photon pT candidates (lead/sublead)
PHO_LEAD_PT_CANDIDATES = ["ALP_lead_photon_pt", "lead_photon_pt", "photon1_pt", "pho1_pt"]
PHO_SUBLEAD_PT_CANDIDATES = ["ALP_sublead_photon_pt", "sublead_photon_pt", "photon2_pt", "pho2_pt"]

def mass_from_sig(sample: str) -> float:
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
    for k, cls in f.classnames().items():
        if "TTree" in cls:
            return f[k]
    for name in ["Events", "tree", "t", treeName]:
        if name in f:
            return f[name]
    return None

def read_hz_arrays(file_path: str):
    """Return (x_pt, h_mass, weight) numpy arrays, or (None, None, None) if not available.
    x_pt is concatenation of (lead_photon_pt, sublead_photon_pt); h_mass/weight are duplicated to match.
    """
    tree = first_tree(file_path)
    if tree is None:
        return None, None, None

    h = None
    for br in H_MASS_CANDIDATES:
        if br in tree.keys():
            try:
                h = tree[br].array(library="np")
                break
            except Exception:
                continue

    # NOTE: y-axis remains mH; x-axis becomes photon pT (lead+sublead)
    pt1 = None
    for br in PHO_LEAD_PT_CANDIDATES:
        if br in tree.keys():
            try:
                pt1 = tree[br].array(library="np")
                break
            except Exception:
                continue

    pt2 = None
    for br in PHO_SUBLEAD_PT_CANDIDATES:
        if br in tree.keys():
            try:
                pt2 = tree[br].array(library="np")
                break
            except Exception:
                continue

    if h is None or pt1 is None or pt2 is None:
        return None, None, None

    w = None
    for wv in WEIGHT_CANDIDATES:
        if wv in tree.keys():
            try:
                w = tree[wv].array(library="np")
                break
            except Exception:
                continue
    if w is None:
        w = np.ones_like(h, dtype=np.float64)

    try:
        h = np.asarray(h, dtype=np.float64)
        pt1 = np.asarray(pt1, dtype=np.float64)
        pt2 = np.asarray(pt2, dtype=np.float64)
        w = np.asarray(w, dtype=np.float64)
    except Exception:
        return None, None, None

    n = min(len(h), len(pt1), len(pt2), len(w))
    h = h[:n]
    pt1 = pt1[:n]
    pt2 = pt2[:n]
    w = w[:n]

    # concatenate lead/sublead into one x-array; duplicate h and w to match
    x_pt = np.concatenate([pt1, pt2])
    y_h = np.concatenate([h, h])
    ww = np.concatenate([w, w])
    return x_pt, y_h, ww

def ensure_outdir_mH_mZ():
    outdir = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/mH_phopT_2D")
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir

def _root_init_style():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetNumberContours(255)

    # 優先嘗試使用 kBird（不同 ROOT 版本可能叫法不同；這裡做最大相容）
    # 1) 若 ROOT 有定義 kBird 常數：直接用
    # 2) 否則嘗試用名稱找 palette id
    # 3) 都失敗就 fallback 到 Viridis
    bird_id = None
    try:
        bird_id = int(getattr(ROOT, "kBird"))
    except Exception:
        bird_id = None

    if bird_id is None:
        try:
            pal = ROOT.TColor.GetColorPalette(0)  # 觸發 palette 系統初始化（某些環境需要）
            # ROOT 6.26+ 可能支援以名稱取 palette（若無此 API 會丟例外）
            bird_id = int(ROOT.TColor.GetPalette("bird"))  # noqa: F821
        except Exception:
            bird_id = None

    if bird_id is not None:
        ROOT.gStyle.SetPalette(bird_id)
    else:
        ROOT.gStyle.SetPalette(ROOT.kViridis)

    ROOT.gStyle.SetTitleFontSize(0.045)
    ROOT.gStyle.SetLabelSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleSize(0.055, "XYZ")

    # Title offset（全域）
    ROOT.gStyle.SetTitleOffset(1.05, "X")
    ROOT.gStyle.SetTitleOffset(1.15, "Y")
    ROOT.gStyle.SetTitleOffset(1.35, "Z")

    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.18)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.08)

def _draw_cms_labels(
    canvas=None,
    lumi_fb: Optional[float] = None,
    lumi_text: Optional[str] = "17.8 fb^{-1} (13.6 TeV)",
):
    """Draw CMS + lumi text safely without relying on outer-scope variables."""
    latex = ROOT.TLatex()
    latex.SetNDC(True)

    left = 0.13
    if canvas is not None:
        try:
            left = float(canvas.GetLeftMargin())
        except Exception:
            pass

    # "CMS" + "Simulation"
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(13)  # left, top
    latex.DrawLatex(left, 0.965, "#bf{CMS} #it{Simulation}")

    # lumi top-right
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(31)  # right, top
    latex.SetNDC(True)
    if lumi_fb is not None:
        lumi_str = f"{lumi_fb:.2f} fb^{{-1}} (13.6 TeV)"
    else:
        lumi_str = lumi_text or ""
    latex.DrawLatex(0.82, 0.93, lumi_str)

def _draw_ma_label(canvas=None, ma_label: Optional[str] = None):
    """Draw mA label at top-right (below lumi). ma_label should be like 'All' or '10'."""
    if not ma_label:
        return
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(31)  # right, top
    latex.DrawLatex(0.80, 0.87, f"m_{{a}} = {ma_label} GeV" if ma_label != "All" else 'm_{a} = All')

def draw_hist2d_mH_mZ(
    z, h, w, z_edges, h_edges, out_pdf: Path, title=None, label_kind: Optional[str] = None, ma_label: Optional[str] = None
):
    """
    PyROOT version: fill TH2D with (x_pt, mH) and weight.
    """
    if z is None or h is None or w is None or len(z) == 0:
        return

    _root_init_style()

    nz = len(z_edges) - 1
    nh = len(h_edges) - 1
    z_arr = np.asarray(z_edges, dtype=np.float64)
    h_arr = np.asarray(h_edges, dtype=np.float64)

    h2 = ROOT.TH2D(
        "h2_mH_phopT",
        "" if title is None else title,
        nz, z_arr,
        nh, h_arr,
    )
    h2.Sumw2()

    # Fill: x=photon pT, y=mH
    for xi, hi, wi in zip(z, h, w):
        if not (np.isfinite(xi) and np.isfinite(hi) and np.isfinite(wi)):
            continue
        h2.Fill(float(xi), float(hi), float(wi))

    if label_kind == "signal":
        ztitle = "Signal MC Events / 2 GeV"
    elif label_kind == "background":
        ztitle = "Background MC Events / 2 GeV"
    else:
        ztitle = "Events"

    h2.GetXaxis().SetTitle("p_{T}^{#gamma} [GeV]")
    h2.GetYaxis().SetTitle("m_{ll#gamma#gamma} [GeV]")
    h2.GetZaxis().SetTitle(ztitle)

    h2.GetXaxis().SetTitleOffset(1.2)
    h2.GetYaxis().SetTitleOffset(1.1)
    h2.GetZaxis().SetTitleOffset(1.2)

    c = ROOT.TCanvas("c_mH_phopT", "c_mH_phopT", 800, 600)
    c.SetLogz(True)

    maxv = h2.GetMaximum()
    if maxv > 0:
        h2.SetMinimum(1e-1)

    h2.Draw("COLZ")

    # Reference lines: keep mH=95,180 horizontal; drop mZ=50 vertical; drop diagonal (was defined in mZ-mH plane)
    x_min = float(h2.GetXaxis().GetXmin())
    x_max = float(h2.GetXaxis().GetXmax())
    y_min = float(h2.GetYaxis().GetXmin())
    y_max = float(h2.GetYaxis().GetXmax())

    line_mh_95 = ROOT.TLine(x_min, 95.0, x_max, 95.0)
    line_mh_95.SetLineStyle(2)
    line_mh_95.SetLineWidth(3)
    line_mh_95.Draw("same")

    line_mh_180 = ROOT.TLine(x_min, 180.0, x_max, 180.0)
    line_mh_180.SetLineStyle(2)
    line_mh_180.SetLineWidth(3)
    line_mh_180.Draw("same")

    # NEW: photon_pT / mH = 10/110  <=>  mH = 11 * photon_pT
    slope = 125.0 / 10.0  # 11
    # find visible segment of y = slope * x within [x_min,x_max]x[y_min,y_max]
    x1 = max(x_min, y_min / slope)
    x2 = min(x_max, y_max / slope)
    if x2 > x1:
        ratio_line = ROOT.TLine(x1, slope * x1, x2, slope * x2)
        ratio_line.SetLineColor(ROOT.kRed + 1)
        ratio_line.SetLineWidth(3)
        ratio_line.SetLineStyle(9)
        ratio_line.Draw("same")

        # 文字標籤：放在線段中點附近，沿法向量偏移避免壓線
        y1 = slope * x1
        y2 = slope * x2
        mx, my = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        dx, dy = (x2 - x1), (y2 - y1)
        nx, ny = -dy, dx  # 法向量
        n = (nx * nx + ny * ny) ** 0.5
        if n > 0:
            nx, ny = nx / n, ny / n

        # 偏移量（用座標軸範圍的比例，對不同 range 較穩）
        off = 0.05 * min((x_max - x_min), (y_max - y_min))
        # 原本是 +off 會在你現在看到的「左邊」；改成 -off 放到線的另一側（右邊）
        tx, ty = mx - off * nx, my - off * ny

        # clamp 到框內（留一點邊界）
        margin = 0.02
        tx = max(x_min + margin * (x_max - x_min), min(tx, x_max - margin * (x_max - x_min)))
        ty = max(y_min + margin * (y_max - y_min), min(ty, y_max - margin * (y_max - y_min)))

        lab = ROOT.TLatex()
        lab.SetNDC(False)  # 使用資料座標
        lab.SetTextFont(42)
        lab.SetTextSize(0.05)
        lab.SetTextAlign(22)  # center, middle
        # 角度跟著線段方向（度）
        try:
            import math
            lab.SetTextAngle(math.degrees(math.atan2(dy, dx)))
        except Exception:
            lab.SetTextAngle(-31)

        lab.SetTextColor(ROOT.kRed + 1)
        lab.DrawLatex(float(tx), float(ty), "p_{T}^{#gamma} / m_{ll#gamma#gamma} > 10/125")

    _draw_cms_labels(canvas=c)
    _draw_ma_label(canvas=c, ma_label=ma_label)

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    c.Print(str(out_pdf))

    c.Close()
    ROOT.gDirectory.Delete("h2_mH_phopT;*")

def collect_background_hzw():
    hs, xs, ws = [], [], []
    missing = 0

    for y in bkg_years_all:
        for s in bkg_samples_by_year.get(y, []):
            fpath = os.path.join(INPUT_BASE, s, f"{y}.root")
            if not os.path.exists(fpath):
                missing += 1
                continue
            x, h, w = read_hz_arrays(fpath)
            if h is None:
                continue
            hs.append(h); xs.append(x); ws.append(w)

    if missing > 0:
        print(f"[bkg] missing files skipped: {missing}")
    if len(hs) == 0:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(hs), np.concatenate(xs), np.concatenate(ws)

def collect_signal_hzw_for_sample(sample: str):
    hs, xs, ws = [], [], []
    missing = 0
    for y in years_sig:
        fpath = os.path.join(INPUT_BASE, sample, f"{y}.root")
        if not os.path.exists(fpath):
            missing += 1
            continue
        x, h, w = read_hz_arrays(fpath)
        if h is None:
            continue
        hs.append(h); xs.append(x); ws.append(w)

    if missing > 0:
        print(f"[sig:{sample}] missing files skipped: {missing}")
    if len(hs) == 0:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(hs), np.concatenate(xs), np.concatenate(ws)

def collect_all_signal_hzw():
    """Collect (mH, x_pt, w) for all signal mass points combined."""
    hs, xs, ws = [], [], []
    missing = 0
    for s in sig_samples:
        for y in years_sig:
            fpath = os.path.join(INPUT_BASE, s, f"{y}.root")
            if not os.path.exists(fpath):
                missing += 1
                continue
            x, h, w = read_hz_arrays(fpath)
            if h is None:
                continue
            hs.append(h); xs.append(x); ws.append(w)

    if missing > 0:
        print(f"[sig:all] missing files skipped: {missing}")
    if len(hs) == 0:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(hs), np.concatenate(xs), np.concatenate(ws)

def main():
    outdir_hz = ensure_outdir_mH_mZ()

    # x-axis: photon pT (adjust range/step as needed)
    x_edges = np.arange(0.0, 100.0 + 2.0 + 1e-9, 2.0)
    # y-axis: mH
    h_edges = np.arange(70.0, 190.0 + 2.0 + 1e-9, 2.0)

    # background：1 張
    bh, bx, bw = collect_background_hzw()
    if bh.size > 0:
        out_pdf = outdir_hz / "mH_vs_phopT_background.pdf"
        draw_hist2d_mH_mZ(bx, bh, bw, x_edges, h_edges, out_pdf, title=None, label_kind="background")
        print(f"[bkg] saved -> {out_pdf}")
    else:
        print("[bkg] no entries found for mH/photon pT.")

    # signal：所有質量點加總：1 張
    ah, ax, aw = collect_all_signal_hzw()
    if ah.size > 0:
        out_pdf = outdir_hz / "mH_vs_phopT_signal_allmA.pdf"
        draw_hist2d_mH_mZ(
            ax, ah, aw, x_edges, h_edges, out_pdf,
            title=None, label_kind="signal", ma_label="All"
        )
        print(f"[sig:all] saved -> {out_pdf}")
    else:
        print("[sig:all] no entries found for combined mH/photon pT.")

    # signal：每個質量點 1 張
    for s in sig_samples:
        ma = mass_from_sig(s)
        sh, sx, sw = collect_signal_hzw_for_sample(s)
        if sh.size == 0:
            continue
        out_pdf = outdir_hz / f"mH_vs_phopT_signal_ma{int(ma)}.pdf"
        draw_hist2d_mH_mZ(
            sx, sh, sw, x_edges, h_edges, out_pdf,
            title=None, label_kind="signal", ma_label=str(int(ma))
        )
        print(f"[sig] saved -> {out_pdf}")

if __name__ == "__main__":
    main()

