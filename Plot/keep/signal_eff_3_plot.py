# ==== New imports/config for plotting ====
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot
import re
import ast
from typing import Dict, List, Tuple, Optional
import awkward as ak
import json  # 新增：讀入 TOTAL_JSON
from datetime import datetime, timezone  # 新增：輸出時間戳

# 嘗試使用 SciPy 的二次樣條插值；若無則會後備到本地實作
try:
    from scipy.interpolate import make_interp_spline
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

# 你之前实际使用的组合：
years_sig  = ["2022preEE"]  # 信号

# 扫描的 ma 列表（背景也按你的原脚本扫）
ma_list = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]

# ===== Samples =====

INPUT_BASE = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sig_eff.json"


# 如需其它：["ZGToLLG","WGToLNuG","ZG2JToG2L2J","EWKZ2J","TT","TTGJets","TGJets","ttWJets","ttZJets","WW","WZ","ZZ"]

# ==== Helpers ====
# 设置全局字体大小和图形样式
plt.rcParams.update({
    'figure.figsize': (8, 6),
    'font.size': 14,
    'axes.labelsize': 20,
    'xtick.labelsize': 17,
    'ytick.labelsize': 17,
    'legend.fontsize': 16,
    'legend.title_fontsize': 16,  # 修正鍵名（原本 legend.titlesize 無效）
    'figure.titlesize': 14,
    # 'patch.linewidth': 1.5,
    # 'patch.edgecolor': 'blue',
})
# 新增：全域線寬設定，方便之後調整
LINE_WIDTH = 2.5

def _load_sig_eff(json_path: str) -> Dict:
    with open(json_path, "r") as f:
        return json.load(f)

def _available_masses(data: Dict, prefer_order: List[int]) -> List[int]:
    by_mass = data.get("by_mass", {})
    present = sorted([int(k) for k in by_mass.keys() if str(k).isdigit()])
    if prefer_order:
        # 依 ma_list 順序保留 JSON 內存在的點
        return [m for m in prefer_order if m in present]
    return present

def _get_total_unweighted(data: Dict, ma: int, channel: str) -> Optional[float]:
    try:
        entry = data["by_mass"][str(ma)]
        if channel == "muon":
            return float(entry["efficiency"]["muon"]["total_unweighted"]) * 100.0
        elif channel == "ele":
            return float(entry["efficiency"]["ele"]["total_unweighted"]) * 100.0
        else:
            return None
    except Exception:
        return None

def _build_series(data: Dict, channel: str, prefer_order: List[int]) -> Tuple[List[int], List[float]]:
    masses = _available_masses(data, prefer_order)
    xs, ys = [], []
    for m in masses:
        v = _get_total_unweighted(data, m, channel)
        if v is not None:
            xs.append(m)
            ys.append(v)
    return xs, ys

def _quadratic_curve(xs: List[int], ys: List[float], n: int = 400) -> Tuple[np.ndarray, np.ndarray]:
    """
    產生二次插值的平滑曲線點集。
    優先使用 SciPy 的 make_interp_spline(k=2)；若不可用，則以滑動三點的局部二次擬合補上。
    """
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    # 確保 x 遞增
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if _HAVE_SCIPY and len(x) >= 3:
        spline = make_interp_spline(x, y, k=2)
        x_new = np.linspace(x.min(), x.max(), n)
        y_new = spline(x_new)
        return x_new, y_new

    # 後備方案：局部二次擬合，區段在相鄰點的中點間拼接
    if len(x) < 3:
        x_new = np.linspace(x.min(), x.max(), n)
        y_new = np.interp(x_new, x, y)
        return x_new, y_new

    x_new_list, y_new_list = [], []
    seg_pts = max(8, n // max(1, len(x) - 2))
    for i in range(len(x) - 2):
        xi = x[i:i+3]
        yi = y[i:i+3]
        coeffs = np.polyfit(xi, yi, deg=2)

        left = x[i] if i == 0 else 0.5 * (x[i] + x[i+1])
        right = x[i+2] if i == len(x) - 3 else 0.5 * (x[i+1] + x[i+2])

        xs_seg = np.linspace(left, right, seg_pts)
        ys_seg = np.polyval(coeffs, xs_seg)

        if x_new_list and xs_seg[0] <= x_new_list[-1][-1]:
            mask = xs_seg > x_new_list[-1][-1]
            xs_seg = xs_seg[mask]
            ys_seg = ys_seg[mask]

        x_new_list.append(xs_seg)
        y_new_list.append(ys_seg)

    x_new = np.concatenate(x_new_list) if x_new_list else x
    y_new = np.concatenate(y_new_list) if y_new_list else y
    return x_new, y_new

def _plot_line(x: List[int], y: List[float], channel_label: str, out_dir: Path) -> Path:
    plt.figure()
    # 以二次曲線連線，並保留原始資料點為 marker
    x_smooth, y_smooth = _quadratic_curve(x, y, n=400)
    color = '#FF0000'
    plt.plot(
        x_smooth, y_smooth,
        linestyle='-',
        color=color,
        linewidth=LINE_WIDTH,  # 新增：設定線寬
        label="2022preEE"
    )
    plt.plot(
        x, y,
        linestyle='None',
        marker='o',
        color=color,
        markersize=8,
        label="_nolegend_"
    )
    plt.xlabel(r"$m_{a}$ [GeV]", labelpad=13)
    plt.ylabel(r"$Efficiency \ \times \ Acceptance \ (\%)$", labelpad=15)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True, length=10)  # 主刻度長度
    plt.minorticks_on()
    plt.tick_params(axis='both', which='minor', direction='in', top=True, right=True, length=6)  # 次刻度長度
    fig = plt.gcf()
    # 將左上角文字拆成粗體 "CMS" 與一般字重 "Preliminary"
    x0, y0 = 0.13, 0.97
    t_cms = fig.text(x0, y0, "CMS", ha="left", va="top", fontsize=19, fontweight="bold")
    fig.canvas.draw()  # 取得 renderer 以量測文字寬度
    renderer = fig.canvas.get_renderer()
    cms_bb = t_cms.get_window_extent(renderer=renderer)
    dx = cms_bb.width / fig.bbox.width + 0.006  # 加一點間距
    fig.text(x0 + dx, y0, "Preliminary", ha="left", va="top", fontsize=19)
    # 右上角亮度與能量維持不變
    fig.text(0.98, 0.965, r"$62.5\,\mathrm{fb}^{-1}\ (13.6\ \mathrm{TeV})$", ha="right", va="top", fontsize=16)
    plt.xlim(0, 32)  # 固定 x 軸範圍到 [0, 32]
    plt.ylim(0.8, 3.5)
    plt.subplots_adjust(left=0.13, right=0.98, bottom=0.14, top=0.92)
    leg = plt.legend(loc="best", title=f"{channel_label} Channel", title_fontsize=18, frameon=False)
    fname = f"sig_eff_{'muon' if channel_label.lower().startswith('muon') else 'ele'}.png"
    out_path = out_dir / fname
    plt.savefig(out_path, dpi=150)
    plt.savefig(out_path.with_suffix(".pdf"))
    plt.close()
    return out_path

def main():
    data = _load_sig_eff(INPUT_BASE)
    out_dir = Path("../plots/signal_eff")
    # 新增：確保輸出資料夾存在
    out_dir.mkdir(parents=True, exist_ok=True)

    # muon
    x_mu, y_mu = _build_series(data, channel="muon", prefer_order=ma_list)
    if x_mu and y_mu:
        _plot_line(x_mu, y_mu, "Muon", out_dir)

    # electron
    x_ele, y_ele = _build_series(data, channel="ele", prefer_order=ma_list)
    if x_ele and y_ele:
        _plot_line(x_ele, y_ele, "Electron", out_dir)

if __name__ == "__main__":
    main()
