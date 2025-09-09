import argparse
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # 避免無顯示環境問題
import matplotlib.pyplot as plt
from typing import Optional  # 新增

# 原始路徑
# PARQUET_PATH = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/Parquet/Sig_MC/ALP_M5_2022preEE/merged_nominal.parquet")
# PARQUET_PATH = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/Parquet/Bkg_MC/DYJetsToLL_2022postEE/merged_nominal.parquet")
PARQUET_PATH = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/Parquet/Bkg_MC/DYGto2LG_50to100_2022postEE/merged_nominal.parquet")

def load_column(parquet_path: Path, column: str) -> pd.Series:
    if not parquet_path.is_file():
        raise FileNotFoundError(f"找不到檔案: {parquet_path}")
    # 只讀需要的欄位避免占記憶體
    try:
        df = pd.read_parquet(parquet_path, columns=[column])
    except Exception as e:
        raise RuntimeError(f"讀取 parquet 失敗: {e}")
    if column not in df.columns:
        raise KeyError(f"欄位 {column} 不存在於檔案 {parquet_path}")
    return df[column]

def plot_distribution(series: pd.Series, out_path: Path, bins: int = 80,
                      max_val: Optional[float] = None,
                      logy: bool = False, show: bool = False,
                      xmin: Optional[float] = None, xmax: Optional[float] = None,
                      normalize: bool = True):  # 新增參數
    data = series.dropna()
    if max_val is not None:
        data = data[data <= max_val]
    if xmin is not None:
        data = data[data >= xmin]
    if xmax is not None:
        data = data[data <= xmax]
    plt.figure(figsize=(7,5))
    plt.hist(
        data,
        bins=bins,
        histtype="stepfilled",
        alpha=0.7,
        color="#1f77b4",
        density=normalize  # 歸一化
    )
    plt.xlabel("ALP_calculatedPhotonIso")
    plt.ylabel("Normalized entries" if normalize else "Entries")  # 動態標籤
    plt.title("Distribution of ALP_calculatedPhotonIso" + (" (normalized)" if normalize else ""))
    if logy:
        plt.yscale("log")
    if xmin is not None or xmax is not None:
        cur_min, cur_max = plt.xlim()
        plt.xlim(xmin if xmin is not None else cur_min,
                 xmax if xmax is not None else cur_max)
    plt.grid(alpha=0.25, linestyle="--")
    plt.tight_layout()
    plt.savefig(out_path)
    if show:
        plt.show()
    plt.close()

def parse_args():
    ap = argparse.ArgumentParser(description="繪製 ALP_calculatedPhotonIso 分布")
    ap.add_argument("--parquet", type=Path, default=PARQUET_PATH, help="Parquet 檔案路徑")
    ap.add_argument("--column", default="ALP_PhotonIso", help="目標欄位名稱")
    ap.add_argument("--bins", type=int, default=10, help="直方圖 bins 數量")
    ap.add_argument("--max", dest="max_val", type=float, default=None, help="裁切最大值 (排除更大的)")
    ap.add_argument("--logy", action="store_true", help="Y 軸使用對數尺度")
    ap.add_argument("--out", type=Path, default=Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/check_area/Plot/ALP_PhotonIso_hist.png"), help="輸出圖檔名稱")
    ap.add_argument("--show", action="store_true", help="顯示圖 (需互動環境)")
    ap.add_argument("--xmin", type=float, default=-0.5, help="x 軸最小值 (同時裁切資料)")
    ap.add_argument("--xmax", type=float, default=None, help="x 軸最大值 (同時裁切資料)")  # 改成預設 None
    # 新增：控制是否歸一化（預設開啟），可用 --no-normalize 關閉
    ap.add_argument("--normalize", dest="normalize", action="store_true", default=True, help="將直方圖歸一化 (預設)")
    ap.add_argument("--no-normalize", dest="normalize", action="store_false", help="關閉直方圖歸一化")
    ap.add_argument("--nonzero-only", action="store_true", help="只使用非 0 事件作圖")
    ap.add_argument("--topk", type=int, default=0, help="列印前 K 大非 0 值 (0 表示不列印)")
    ap.add_argument("--auto-range", action="store_true", help="用分位數自動設定 x 軸範圍")
    ap.add_argument("--qmin", type=float, default=0.00, help="自動範圍的下分位數（含），例如 0.00")
    ap.add_argument("--qmax", type=float, default=0.99, help="自動範圍的上分位數（含），例如 0.99")
    return ap.parse_args()

def main():
    args = parse_args()
    try:
        series = load_column(args.parquet, args.column)
    except Exception as e:
        print(f"[錯誤] {e}", file=sys.stderr)
        sys.exit(1)

    # 基本統計
    nz_mask = series.fillna(0) > 0
    n_all = series.shape[0]
    n_nz = int(nz_mask.sum())
    print(f"讀取 {args.parquet} 欄位 {args.column}，有效筆數: {series.notna().sum()}，非 0 筆數: {n_nz}/{n_all}")

    # 只使用非 0
    if args.nonzero_only:
        series = series[nz_mask]

    # 依需求列印前 K 大非 0 值
    if args.topk > 0:
        top_vals = series[series > 0].nlargest(args.topk)
        print(f"[非 0 前 {args.topk} 大]：")
        for i, v in enumerate(top_vals):
            print(f"  #{i+1}: {float(v):.6g}")

    # 自動 x 範圍（避免裁掉尾端）
    xmin, xmax = args.xmin, args.xmax
    if args.auto_range:
        s = series.dropna()
        s = s[s > 0] if args.nonzero_only else s
        if len(s) > 0:
            q_lo = float(s.quantile(args.qmin))
            q_hi = float(s.quantile(args.qmax))
            xmin = q_lo if xmin is None else xmin
            xmax = q_hi if xmax is None else xmax
            print(f"[auto-range] x in [{xmin:.6g}, {xmax:.6g}]（由分位數 qmin={args.qmin}, qmax={args.qmax} 推得）")

    if args.max_val is not None:
        kept = (series <= args.max_val).sum()
        print(f"套用最大值裁切 {args.max_val} 後剩 {kept} 筆 (原始 {len(series)})")

    out_path = args.out
    if out_path.is_dir():
        out_path = out_path / f"{args.column}_hist.png"

    plot_distribution(series, out_path, bins=args.bins, max_val=args.max_val,
                      logy=args.logy, show=args.show,
                      xmin=xmin, xmax=xmax,
                      normalize=args.normalize)  # 傳入是否歸一化

    print(f"已輸出圖檔: {out_path.resolve()}")

if __name__ == "__main__":
    main()