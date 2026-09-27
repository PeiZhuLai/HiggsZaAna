#!/usr/bin/env python3
"""Direct sigma_eff comparison of the signal samples, before vs after the FSR fix.

This is the fast answer to "how much did the resolution improve": it reads the
merged signal parquets of the old and the new production and computes the
effective resolution of m(llgg) and m(ll) straight from the events, without
going through p2root / BDT / the flashgg signal fit. The signal-model effSigma
that ends up in the AN still has to come from the full chain afterwards -- this
measures the input to it.

Built-in control: the FSR recovery is applied to muons only, so the ELECTRON
channel must not move. If it does, something other than the FSR fix changed.

sigma_eff = half-width of the narrowest interval containing 68.3% of the events.

Input : <old>/<sample>_<era>/merged_nominal.parquet
        <new>/<sample>_<era>/merged_nominal.parquet
Output: printed table + <out>/sigma_eff_comparison.txt

Usage : python compare_sigma_eff.py [--new DIR] [--old DIR] [--out DIR]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

NEW = ("/eos/project/h/htozg-dy-privatemc/pelai/HZa/"
       "parquet_DNA_tmp_fsrfix/Sig_MC")
OLD = "/eos/home-p/pelai/HZa/parquet_DNA/Sig_MC"
MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
COLUMNS = ["z_mumu", "z_ee", "H_mass", "Z_mass", "Z_noFSR_mass"]


def sigma_eff(values: np.ndarray) -> float:
    values = np.sort(np.asarray(values, dtype=float))
    n = len(values)
    if n < 20:
        return float("nan")
    k = int(round(0.683 * n))
    widths = values[k:] - values[: n - k]
    return float(np.min(widths) / 2.0)


def read(path: Path):
    if not path.exists():
        return None
    table = pq.read_table(path, columns=COLUMNS)
    return {c: table[c].combine_chunks().to_numpy(zero_copy_only=False)
            for c in table.column_names}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--new", default=NEW)
    parser.add_argument("--old", default=OLD)
    parser.add_argument("--out", default="/afs/cern.ch/work/p/pelai/HZa/"
                                         "HiggsZaAna/studies/l3_review_20260910/logs_fsrfix")
    args = parser.parse_args()

    lines = []
    header = ("%-6s %-14s %7s %9s %9s %8s   %9s %9s %8s"
              % ("ma", "era", "Nmu", "seff_mu0", "seff_mu1", "d[%]",
                 "seff_ee0", "seff_ee1", "d[%]"))
    lines.append("sigma_eff of m(llgg): 0 = before the FSR fix, 1 = after")
    lines.append("the electron column is the control -- it must not move")
    lines.append(header)
    lines.append("-" * len(header))

    deltas_mu, deltas_ee = [], []
    for ma in MASSES:
        for era in ERAS:
            name = "mA_M%d_%s" % (ma, era)
            old = read(Path(args.old) / name / "merged_nominal.parquet")
            new = read(Path(args.new) / name / "merged_nominal.parquet")
            if old is None or new is None:
                lines.append("%-6s %-14s %7s  (missing %s)"
                             % (ma, era, "-", "old" if old is None else "new"))
                continue

            row = [ma, era]
            cells = []
            for flag in ("z_mumu", "z_ee"):
                s0 = sigma_eff(old["H_mass"][old[flag] == 1])
                s1 = sigma_eff(new["H_mass"][new[flag] == 1])
                pct = 100.0 * (s1 - s0) / s0 if s0 and np.isfinite(s0) else float("nan")
                cells.append((s0, s1, pct))
                (deltas_mu if flag == "z_mumu" else deltas_ee).append(pct)
            n_mu = int((new["z_mumu"] == 1).sum())
            lines.append("%-6s %-14s %7d %9.4f %9.4f %+8.2f   %9.4f %9.4f %+8.2f"
                         % (row[0], row[1], n_mu,
                            cells[0][0], cells[0][1], cells[0][2],
                            cells[1][0], cells[1][1], cells[1][2]))

    def summarize(label, values):
        v = np.asarray([x for x in values if np.isfinite(x)])
        if not len(v):
            return "%s: no finite entries" % label
        return ("%s: mean %+.2f%%  median %+.2f%%  best %+.2f%%  worst %+.2f%%  (n=%d)"
                % (label, v.mean(), np.median(v), v.min(), v.max(), len(v)))

    lines.append("")
    lines.append(summarize("muon channel  ", deltas_mu))
    lines.append(summarize("electron ctrl ", deltas_ee))
    ee = np.asarray([x for x in deltas_ee if np.isfinite(x)])
    if len(ee) and np.max(np.abs(ee)) > 0.5:
        lines.append("WARNING: the electron control moved by up to %.2f%% -- "
                     "the FSR fix is muon-only, so something else changed."
                     % np.max(np.abs(ee)))

    text = "\n".join(lines)
    print(text)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "sigma_eff_comparison.txt").write_text(text + "\n")
    print("\nwrote", out / "sigma_eff_comparison.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
