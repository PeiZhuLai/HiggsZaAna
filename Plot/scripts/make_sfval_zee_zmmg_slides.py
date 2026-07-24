#!/usr/bin/env python3
"""Beamer slides for the zee_zmmg SF-validation plots.

每 (era x SF 群組) 一張投影片，該群組的 2-4 個變數圖排成 grid，依 era 分 section。
圖來源: Plot/plots/SF_validation_zee_zmmg/<era>/<sf>/<var>.pdf（用絕對路徑，故 tex 可放任意目錄）。
輸出: <out>/sfval_zee_zmmg_slides.tex（+ pdflatex 編成 .pdf）。
"""
import os, argparse

PLOT_ROOT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/SF_validation_zee_zmmg"

# era 顯示名 + lumi 標籤（fb^-1, 13.6 TeV）
ERAS = [
    ("2022preEE",    "2022 preEE",    "7.98"),
    ("2022postEE",   "2022 postEE",   "26.67"),
    ("2023preBPix",  "2023 preBPix",  "18.06"),
    ("2023postBPix", "2023 postBPix", "9.69"),
    ("2024",         "2024",          "108.83"),
    ("run3_all",     "Run 3 (all)",   "171.24"),
]

EE  = ["lead_electron_pt", "sublead_electron_pt", "lead_electron_eta", "mee"]
MU  = ["lead_muon_pt", "sublead_muon_pt", "lead_muon_eta", "mmmg"]
PHO = ["photon_pt", "photon_eta"]

# (dir, 顯示標題, 變數清單)
GROUPS = [
    ("electron_reco_sf",    "Electron reco SF",    EE),
    ("electron_id_sf",      "Electron ID SF",      EE),
    ("electron_trigger_sf", "Electron trigger SF", EE),
    ("muon_reco_sf",        "Muon reco SF",        MU),
    ("muon_id_sf",          "Muon ID SF",          MU),
    ("muon_trigger_sf",     "Muon trigger SF",     MU),
    ("photon_id_sf",        "Photon ID SF",        PHO),
    ("photon_csev_sf",      "Photon CSEV SF",      PHO),
]

HEADER = r"""\documentclass[aspectratio=169,10pt]{beamer}
\usetheme{default}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{footline}[frame number]
\usepackage{graphicx}
\usepackage[T1]{fontenc}
\setbeamerfont{frametitle}{size=\normalsize}
\title{Lepton/Photon SF validation}
\subtitle{$Z\to ee$ and $Z\to\mu\mu\gamma$ control sample (zee\_zmmg)}
\author{Pei-Zhu Lai}
\date{\today}
\begin{document}
\frame{\titlepage}
"""

FOOTER = r"\end{document}" + "\n"


def esc(s):
    return s.replace("_", r"\_")


def frame(era_dir, era_disp, lumi, g_title, gdir, varlist):
    paths = [os.path.join(PLOT_ROOT, era_dir, gdir, f"{v}.pdf") for v in varlist]
    paths = [(p if os.path.exists(p) else None) for p in paths]
    title = f"{esc(era_disp)} --- {esc(g_title)}"
    subtitle = rf"{lumi}\,fb$^{{-1}}$ (13.6 TeV)"
    out = [rf"\begin{{frame}}{{{title}}}", r"\centering"]
    present = [p for p in paths if p]
    if len(varlist) >= 4:  # 2x2
        h = "0.36\\textheight"
        rows = [present[i:i + 2] for i in range(0, len(present), 2)]
    else:                  # 1xN（photon 2 圖）
        h = "0.62\\textheight"
        rows = [present]
    body = []
    for row in rows:
        imgs = [rf"\includegraphics[height={h},keepaspectratio]{{{p}}}" for p in row]
        body.append(r"\hfill ".join(imgs))
    out.append((r"\\[2pt]" + "\n").join(body))
    out.append(rf"\\[1pt]{{\scriptsize {subtitle}}}")
    out.append(r"\end{frame}")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(PLOT_ROOT, "slides"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    tex = [HEADER]
    for era_dir, era_disp, lumi in ERAS:
        tex.append(rf"\section{{{esc(era_disp)}}}")
        tex.append(rf"\begin{{frame}}\vfill\centering\Large {esc(era_disp)}\\[4pt]"
                   rf"\normalsize {lumi}\,fb$^{{-1}}$ (13.6 TeV)\vfill\end{{frame}}")
        for gdir, g_title, varlist in GROUPS:
            tex.append(frame(era_dir, era_disp, lumi, g_title, gdir, varlist))
    tex.append(FOOTER)

    texpath = os.path.join(args.out, "sfval_zee_zmmg_slides.tex")
    with open(texpath, "w") as f:
        f.write("\n".join(tex))
    print(f"[tex] {texpath}  frames={len(ERAS)*(len(GROUPS)+1)+1}")


if __name__ == "__main__":
    main()
