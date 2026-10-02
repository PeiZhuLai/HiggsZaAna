#!/usr/bin/env python3
"""FSR recovery on the FINAL samples: without vs with the FSR photon, same events.

The final production stores both the dressed and the undressed reconstruction in
every event (Z_mass / Z_noFSR_mass, H_mass / H_noFSR_mass), so the "before" of the
before/after comparison does not need the deleted pre-fix products: it is the
same event reconstructed without adding the FSR photon to the muon.

Input : signal  parquet_DNA_tmp_fsrfix/Sig_MC/mA_M<ma>_<era>/merged_nominal.parquet
        data    parquet_DNA_tmp_fsrfix_fpo1/Data/Data_<era>/merged_nominal.parquet
        DY      parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC/{DYJetsToLL,DYGto2LG_10to100}_<era>
                (2024 DY+jets from parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyveto2024,
                 i.e. the re-production with the overlap veto)
Output: <out>/fsr_numbers_signal.txt, <out>/fsr_numbers_datamc.txt, <out>/fsr_numbers.json

sigma_eff = half-width of the narrowest interval holding 68.3% of the (unweighted) events.
"dressed" = an FSR photon was added, |Z_mass - Z_noFSR_mass| > 0.01 GeV.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

P = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"
SIG = f"{P}/parquet_DNA_tmp_fsrfix/Sig_MC"
DATA = f"{P}/parquet_DNA_tmp_fsrfix_fpo1/Data"
BKG = f"{P}/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC"
BKG24 = f"{P}/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyveto2024"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
COLS = ["z_mumu", "z_ee", "Z_mass", "Z_noFSR_mass", "H_mass", "H_noFSR_mass",
        "gamma_fsr_pt", "weight_central"]


def bkg_files(era):
    if era == "2024":
        dy = [f"{BKG24}/DYJetsTo2{f}_2024/merged_nominal.parquet" for f in ("E", "Mu", "Tau")]
    else:
        dy = [f"{BKG}/DYJetsToLL_{era}/merged_nominal.parquet"]
    return dy + [f"{BKG}/DYGto2LG_10to100_{era}/merged_nominal.parquet"]


def read(paths):
    if isinstance(paths, str):
        paths = [paths]
    parts = []
    for p in paths:
        t = pq.read_table(p, columns=COLS)
        parts.append({c: t[c].combine_chunks().to_numpy(zero_copy_only=False).astype(float)
                      for c in COLS})
    return {c: np.concatenate([d[c] for d in parts]) for c in COLS}


def sigma_eff(v):
    v = np.sort(np.asarray(v, float))
    n = len(v)
    if n < 20:
        return float("nan")
    k = int(round(0.683 * n))
    return float(np.min(v[k:] - v[:n - k]) / 2.0)


def summary(d):
    mu = d["z_mumu"] == 1
    ee = d["z_ee"] == 1
    dressed = mu & (np.abs(d["Z_mass"] - d["Z_noFSR_mass"]) > 0.01)
    out = dict(n_mu=int(mu.sum()), n_ee=int(ee.sum()),
               dressed_frac=float(dressed.sum() / max(mu.sum(), 1)),
               ee_dressed=int((ee & (np.abs(d["Z_mass"] - d["Z_noFSR_mass"]) > 0.01)).sum()),
               fsr_pt_median=float(np.median(d["gamma_fsr_pt"][dressed])) if dressed.any() else float("nan"))
    for tag, sel in (("mu", mu), ("ee", ee), ("dressed", dressed)):
        for var, b, a in (("mll", "Z_noFSR_mass", "Z_mass"), ("mllgg", "H_noFSR_mass", "H_mass")):
            s0, s1 = sigma_eff(d[b][sel]), sigma_eff(d[a][sel])
            out[f"seff_{var}_{tag}"] = [s0, s1, 100.0 * (s1 - s0) / s0 if s0 == s0 and s0 > 0 else float("nan")]
    # Z-peak core: sigma_eff of m_ll inside 80-100 GeV of the undressed mass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "output"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {"signal": {}, "data": {}, "dy": {}}

    lines = ["sigma_eff [GeV] without -> with FSR recovery (change in %); same events",
             "%-3s %-13s %6s %7s %6s | %-24s | %-24s | %-24s | %-24s | %s" % (
                 "ma", "era", "Nmu", "dress%", "pTmed", "mll mu", "mllgg mu",
                 "mll mu dressed", "mllgg mu dressed", "mllgg ee (control)")]
    for ma in MASSES:
        allmu = []
        for era in ERAS:
            d = read(f"{SIG}/mA_M{ma}_{era}/merged_nominal.parquet")
            s = summary(d)
            res["signal"][f"{ma}_{era}"] = s
            f = lambda k: "%6.3f->%6.3f (%+5.2f)" % tuple(s[k])
            lines.append("%-3d %-13s %6d %7.2f %6.2f | %s | %s | %s | %s | %s" % (
                ma, era, s["n_mu"], 100 * s["dressed_frac"], s["fsr_pt_median"],
                f("seff_mll_mu"), f("seff_mllgg_mu"), f("seff_mll_dressed"),
                f("seff_mllgg_dressed"), f("seff_mllgg_ee")))
    dm = np.array([v["seff_mllgg_mu"][2] for v in res["signal"].values()])
    dl = np.array([v["seff_mll_mu"][2] for v in res["signal"].values()])
    de = np.array([v["seff_mllgg_ee"][2] for v in res["signal"].values()])
    nee = sum(v["ee_dressed"] for v in res["signal"].values())
    lines += ["",
              "mllgg mu : mean %+.2f%% median %+.2f%% best %+.2f%% worst %+.2f%% improved %d/%d"
              % (dm.mean(), np.median(dm), dm.min(), dm.max(), (dm < 0).sum(), len(dm)),
              "mll   mu : mean %+.2f%% median %+.2f%% best %+.2f%% worst %+.2f%% improved %d/%d"
              % (dl.mean(), np.median(dl), dl.min(), dl.max(), (dl < 0).sum(), len(dl)),
              "mllgg ee : max |change| %.4f%%, dressed electron-channel events: %d" % (np.abs(de).max(), nee)]
    for ma in MASSES:
        v = [res["signal"][f"{ma}_{e}"] for e in ERAS]
        lines.append("ma=%-2d  dressed %.2f-%.2f%%  mllgg mu %+.2f..%+.2f%%  mll mu %+.2f..%+.2f%%" % (
            ma, 100 * min(x["dressed_frac"] for x in v), 100 * max(x["dressed_frac"] for x in v),
            min(x["seff_mllgg_mu"][2] for x in v), max(x["seff_mllgg_mu"][2] for x in v),
            min(x["seff_mll_mu"][2] for x in v), max(x["seff_mll_mu"][2] for x in v)))
    (out / "fsr_numbers_signal.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-20:]))

    lines = ["%-6s %-13s %8s %7s %6s | %-24s" % ("sample", "era", "Nmu", "dress%", "pTmed", "mll mu")]
    for era in ERAS:
        for key, paths in (("data", f"{DATA}/Data_{era}/merged_nominal.parquet"), ("dy", bkg_files(era))):
            s = summary(read(paths))
            res[key][era] = s
            lines.append("%-6s %-13s %8d %7.2f %6.2f | %6.3f->%6.3f (%+5.2f)" % (
                key, era, s["n_mu"], 100 * s["dressed_frac"], s["fsr_pt_median"], *s["seff_mll_mu"]))
    (out / "fsr_numbers_datamc.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    (out / "fsr_numbers.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
