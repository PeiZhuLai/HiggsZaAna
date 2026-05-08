# ==== New imports/config for plotting ====
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot

INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal/"
INPUT_BASE_MERGED = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal/"
INPUT_BASE_FALLBACKS = [
    "/eos/home-p/pelai/HZa/root_P2Root/run3_mergedBDT/",
]

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
years_sig = year_sig_2022 + year_sig_2023 + year_sig_2024
sig_samples = list(dict.fromkeys(name_sig_2022 + name_sig_2023 + name_sig_2024))

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
interpolate = True

if not interpolate:
    ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
else:
    ma_list = list(range(1, 31))

tuneStyle = False

sig_ma_ticks =  [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
bkg_ma_ticks = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]
# bkg_ma_ticks = [5,15,30]

MVA_CANDIDATES = ["MVA_Score"]
WEIGHT_CANDIDATES = ["weight", "w"]
MVA_BRANCH_PREFIXES = [
    "MVA_Score_mA_M",
    "MVA_Score_ma_M",
    "BDT_Score_mA_M",
    "BDT_Score_ma_M",
]
#-----------------------------------------------------

# ==== Helpers ====
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 14,
    'axes.titlesize': 17,
    'axes.labelsize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'figure.titlesize': 14,
})

def mass_from_sig(sample: str) -> float:
    if "_M" in sample:
        try:
            return float(sample.split("_M", 1)[1])
        except Exception:
            pass
    return float("nan")

def _unique_paths(paths):
    out = []
    seen = set()
    for p in paths:
        if p is None:
            continue
        norm = str(p).rstrip("/")
        if norm and norm not in seen:
            out.append(norm)
            seen.add(norm)
    return out

def _input_bases(primary):
    return _unique_paths([primary, INPUT_BASE, INPUT_BASE_MERGED] + INPUT_BASE_FALLBACKS)

def resolve_root_path(sample: str, year: str, primary_base: str):
    first_candidate = None
    for base in _input_bases(primary_base):
        path = os.path.join(base, sample, f"{year}.root")
        if first_candidate is None:
            first_candidate = path
        if os.path.exists(path):
            return path, True
    return first_candidate, False

def first_tree(file_path: str, preferred_names=None):
    try:
        f = uproot.open(file_path)
    except Exception:
        return None
    names = list(preferred_names or []) + ["test", "inclusive", "Events", "tree", "t"]
    seen = set()
    for name in names:
        if name in seen:
            continue
        seen.add(name)
        if name in f:
            try:
                if "TTree" in f[name].classname:
                    return f[name]
            except Exception:
                continue
    for k, cls in f.classnames().items():
        if "TTree" in cls:
            return f[k]
    return None

def read_weight_from_tree(tree, size):
    if tree is None:
        return None
    for wv in WEIGHT_CANDIDATES:
        if wv in tree.keys():
            try:
                return np.asarray(tree[wv].array(library="np"), dtype=np.float64)
            except Exception:
                return None
    return np.ones(int(size), dtype=np.float64)

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
    try:
        mva = np.asarray(mva).astype(np.float64)
        w = np.asarray(w).astype(np.float64)
    except Exception:
        return None, None
    n = min(len(mva), len(w))
    return mva[:n], w[:n]

def read_mva_for_ma_from_tree(tree, ma: int):
    """
    Return mva numpy array for given ma from a merged-background tree.
    Tries exact name and common uproot-suffixed variants like '.0'.
    """
    if tree is None:
        return None

    keys = set(tree.keys())
    cand = []
    for pfx in MVA_BRANCH_PREFIXES:
        base = f"{pfx}{int(ma)}"
        cand.extend([base, base + ".0", base + ".1", base + "_0", base + "_1"])
    cand.extend(["MVA_Score", "MVA_Score.0", "MVA_Score.1"])

    for br in cand:
        if br in keys:
            try:
                arr = tree[br].array(library="np")
                return np.asarray(arr, dtype=np.float64)
            except Exception:
                return None
    return None

def build_mass_edges(masses):
    ms = np.array(sorted(set(float(m) for m in masses)))
    if len(ms) == 0:
        return None
    if len(ms) == 1:
        dm = 1.0
        return np.array([ms[0] - dm/2, ms[0] + dm/2])
    mids = (ms[1:] + ms[:-1]) / 2.0
    edges = np.concatenate(([ms[0] - (mids[0] - ms[0])], mids, [ms[-1] + (ms[-1] - mids[-1])]))
    return edges

def build_uniform_mass_edges(masses, step=1.0):
    ms = [float(m) for m in masses]
    if len(ms) == 0:
        return None
    mmin = int(np.floor(min(ms)))
    mmax = int(np.ceil(max(ms)))
    start = mmin - 0.5
    stop = mmax + 0.5 + 1e-9
    return np.arange(start, stop, step)

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
    if not style_only:
        H, xe, ye = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=w)
    else:
        xe, ye = x_edges, y_edges
        nx = len(xe) - 1
        ny = len(ye) - 1
        gx = np.linspace(0.3, 1.0, max(nx, 1))
        gy = np.linspace(0.3, 1.0, max(ny, 1))
        H = np.outer(gx[:nx], gy[:ny])
        if H.size == 0:
            H = np.array([[1.0]])

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

    if corr_text:
        ax = plt.gca()
        if corr_pos is not None:
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

    fig = plt.gcf()
    x0, y0 = 0.13, 0.97
    t_cms = fig.text(x0, y0, "CMS", ha="left", va="top", fontsize=19, fontweight="bold")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    cms_bb = t_cms.get_window_extent(renderer=renderer)
    dx = cms_bb.width / fig.bbox.width + 0.006
    fig.text(x0 + dx, y0, "Preliminary", ha="left", va="top", fontsize=19)
    fig.text(0.84, 0.965, r"$172.13\,\mathrm{fb}^{-1}\ (13.6\ \mathrm{TeV})$", ha="right", va="top", fontsize=16)
    plt.subplots_adjust(left=0.13, right=0.98, bottom=0.14, top=0.92)

    out_pdf = out_png.with_suffix(".pdf")
    plt.savefig(out_pdf, bbox_inches="tight", pad_inches=0.05)
    plt.close()

def make_style_previews(outdir):
    bdt_edges = np.linspace(0.0, 1.0, 101)

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
    no_tree = 0
    no_branch = 0

    for y in bkg_years_all:
        for s in bkg_samples_by_year.get(y, []):
            fpath, found = resolve_root_path(s, y, INPUT_BASE_MERGED)
            if not found:
                missing += 1
                continue

            tree = first_tree(fpath, preferred_names=["inclusive"])
            if tree is None:
                no_tree += 1
                continue

            w_all = None

            for ma in ma_list:
                mva = read_mva_for_ma_from_tree(tree, ma)
                if mva is None:
                    no_branch += 1
                    continue

                if w_all is None:
                    w_all = read_weight_from_tree(tree, len(mva))
                if w_all is None:
                    w = np.ones(len(mva), dtype=np.float64)
                else:
                    n = min(len(mva), len(w_all))
                    mva = mva[:n]
                    w = w_all[:n]

                mask = np.isfinite(mva) & np.isfinite(w)
                if mask.size != 0:
                    mva = mva[mask]
                    w = w[mask]
                if len(mva) == 0:
                    continue

                xs.append(mva)
                ws.append(w)
                ys.append(np.full_like(mva, float(ma), dtype=np.float64))

    if missing > 0:
        print(f"[bkg] missing files skipped: {missing}")
    if no_tree > 0 or no_branch > 0:
        print(f"[bkg] skipped trees={no_tree}, missing MVA branches={no_branch}")
    if len(xs) == 0:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(ws)

def collect_signal_xyw():
    xs, ys, ws = [], [], []
    missing = 0
    no_tree = 0
    no_branch = 0
    for s in sig_samples:
        ma = mass_from_sig(s)
        for y in years_sig:
            fpath, found = resolve_root_path(s, y, INPUT_BASE)
            if not found:
                missing += 1
                continue

            tree = first_tree(fpath, preferred_names=["test", "inclusive"])
            if tree is None:
                no_tree += 1
                continue

            mva = read_mva_for_ma_from_tree(tree, int(ma))
            if mva is None:
                no_branch += 1
                continue

            w = read_weight_from_tree(tree, len(mva))
            if w is None:
                w = np.ones(len(mva), dtype=np.float64)
            n = min(len(mva), len(w))
            mva = mva[:n]
            w = w[:n]
            mask = np.isfinite(mva) & np.isfinite(w)
            if mask.size != 0:
                mva = mva[mask]
                w = w[mask]
            if len(mva) == 0:
                continue

            xs.append(mva)
            ws.append(w)
            ys.append(np.full_like(mva, float(ma), dtype=np.float64))
    if missing > 0:
        print(f"[sig] missing files skipped: {missing}")
    if no_tree > 0 or no_branch > 0:
        print(f"[sig] skipped trees={no_tree}, missing MVA branches={no_branch}")
    if len(xs) == 0:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(ws)

# ==== Main ====
def main():
    outdir = ensure_outdir()
    print(f"[input] ROOT base candidates: {', '.join(_input_bases(INPUT_BASE))}")

    if tuneStyle:
        make_style_previews(outdir)
        return

    bdt_edges = np.linspace(0.0, 1.0, 101)

    bx, by, bw = collect_background_xyw()
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
            xtick_labelsize=(14 if interpolate else None)
        )
        print(f"[bkg] saved -> {outdir/'BDT_vs_ma_background.png'} (corr={corr_bkg:.3f})")
    else:
        print("[bkg] no entries found.")

    sx, sy, sw = collect_signal_xyw()
    if sx.size > 0:
        # extract ma from ma_list
        sig_x_edges = build_uniform_mass_edges(ma_list, step=1.0)
        corr_sig = weighted_corr(sy, sx, sw)
        draw_hist2d(
            sy, sx, sw, sig_x_edges, bdt_edges,
            title=None,
            out_png=outdir / "BDT_vs_ma_signal.png",
            y_label="Signal BDT Score",
            x_tick_masses=sig_ma_ticks,
            corr_text=f"Correlation = {corr_sig:.3f}",
            corr_pos=(0.53, 0.5),
            xtick_labelsize=(14 if interpolate else None)
        )
        print(f"[sig] saved -> {outdir/'BDT_vs_ma_signal.png'} (corr={corr_sig:.3f})")
    else:
        raise RuntimeError(
            "[sig] no entries found. Check signal ROOT path candidates and "
            "MVA_Score_mA_M* branches."
        )

if __name__ == "__main__":
    main()
