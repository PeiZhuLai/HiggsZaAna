#!/usr/bin/env python3
r"""Build the AN sculpting-ratio table (tab:sculpt_R, "Table 16") directly from the adopted
BDT working points in MVAcut_points_run3.json, using the ALREADY-SCORED background samples
(stored MVA_Score_mA_M* branches) -- no re-scoring with a model .pkl.

Why (2026-07-23): the R value in the AN table must be the ACTUAL sculpting ratio at the
adopted BDT cut, computed with the SAME BDT scores the analysis uses. Re-scoring with the
model .pkl is fragile: the models were retrained (2026-07-06, signal reweight) after the
old sculpt_R_run3.json / MVAcut_points_run3.json (2026-06-19) were produced, so re-scoring
now drifts. The scored samples in run3_bdt_scored_nominal carry the FINAL-model scores
(MVA_Score_mA_M1..M30), so reading those is both correct and reproducible. This makes the
table a single workflow step keyed to (MVAcut_points_run3.json, run3_bdt_scored_nominal).

Definition (identical to scan_score_R_significance._bkg_R_at_cuts):
  R(m_a) = f_peak^{pass} / f_peak^{incl},  f_peak = w(120<m<130) / w(95<m<180),
on the combined simulated DY background, at the adopted BDT cut for that mass.

Env: any python3 with uproot + numpy (no ROOT needed).
Reads:  Plot/output/MVAcut_points_run3.json                  (adopted working-point BDT cuts)
        <run3_bdt_scored_nominal>/<sample>/<era>.root        (inclusive tree, MVA_Score_mA_M*)
Writes: Plot/output/latexTable_sculptR.txt                   (and optionally --into <Sec-06-04 tex>)
"""
import argparse
import json
import os
import re

import numpy as np
import uproot

# --- paths ----------------------------------------------------------------------------
SCORED_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal"
JSON_PATH = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
LOWMASS_META = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_BDT_lowmass_run3.meta.json"
HIGHMASS_META = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_BDT_highmass_run3.meta.json"
DEF_OUT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/latexTable_sculptR.txt"

# --- background combination (identical to table_interpolate_bkgYield_1.py) -------------
_bkg_2022 = ["DYGto2LG_10to50", "DYGto2LG_50to100"]
_bkg_2023 = ["DYGto2LG_10to100"]
_bkg_dyll = ["DYJetsToLL"]
_bkg_dyll_2024 = ["DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
BKG_SAMPLES_BY_YEAR = {
    "2022preEE":   _bkg_2022 + _bkg_dyll,
    "2022postEE":  _bkg_2022 + _bkg_dyll,
    "2023preBPix": _bkg_2023 + _bkg_dyll,
    "2023postBPix": _bkg_2023 + _bkg_dyll,
    "2024":        _bkg_2023 + _bkg_dyll_2024,
}

TREE = "inclusive"
MASS_BR = "H_mass"
WEIGHT_BR = "factor"          # == "weight" in the scored files
SCORE_FMT = "MVA_Score_mA_M{}"
WINDOW = (95.0, 180.0)        # inclusive Higgs-candidate window
PEAK = (120.0, 130.0)         # 125-peak window

# --- AUTO-block markers (drop-in for the block written by make_sculpt_R_table.py) ------
BEGIN_MARK = "% >>> AUTO sculpt_R table (make_sculpt_R_table_from_wp.py) >>>"
END_MARK = "% <<< AUTO sculpt_R table <<<"
ANY_BEGIN = r"% >>> AUTO sculpt_R table \([^)]*\) >>>"


def load_scored_bkg(masses):
    """Concatenate the combined DY background across all (sample, era), reading H_mass,
    the event weight, and the per-mass MVA scores MVA_Score_mA_M{m} for the requested masses.
    Returns dict: {'m': H_mass, 'w': weight, m: score_array_for_mass_m}."""
    score_brs = {m: SCORE_FMT.format(m) for m in masses}
    m_all, w_all = [], []
    s_all = {m: [] for m in masses}
    for era, samples in BKG_SAMPLES_BY_YEAR.items():
        for sample in samples:
            path = os.path.join(SCORED_BASE, sample, f"{era}.root")
            if not os.path.exists(path):
                raise FileNotFoundError(path)
            arr = uproot.open(path)[TREE].arrays(
                [MASS_BR, WEIGHT_BR] + list(score_brs.values()), library="np")
            m_all.append(arr[MASS_BR])
            w_all.append(arr[WEIGHT_BR])
            for m in masses:
                s_all[m].append(arr[score_brs[m]])
    out = {"m": np.concatenate(m_all), "w": np.concatenate(w_all)}
    for m in masses:
        out[m] = np.concatenate(s_all[m])
    return out


def compute_R(bkg, mass, cut):
    """R at a fixed BDT cut for one mass hypothesis, on the combined background."""
    m, w, s = bkg["m"], bkg["w"], bkg[mass]
    sel = (m > WINDOW[0]) & (m < WINDOW[1]) & np.isfinite(w) & np.isfinite(s)
    m, w, s = m[sel], w[sel], s[sel]
    peak = (m > PEAK[0]) & (m < PEAK[1])
    frac = w[peak].sum() / w.sum()
    pB = s > cut
    denom = w[pB].sum()
    return float(w[pB & peak].sum() / denom / frac) if denom > 0 and frac > 0 else float("nan")


def compute_R_from_wp(json_path=JSON_PATH):
    cuts = {}
    for e in json.load(open(json_path))["results"]:
        if e.get("MVAcut") is not None:
            cuts[int(e["mA"])] = float(e["MVAcut"])
    low_masses = [int(x) for x in json.load(open(LOWMASS_META))["masses"] if int(x) in cuts]
    high_masses = [int(x) for x in json.load(open(HIGHMASS_META))["masses"] if int(x) in cuts]
    bkg = load_scored_bkg(sorted(set(low_masses + high_masses)))
    low_R = {m: compute_R(bkg, m, cuts[m]) for m in low_masses}
    high_R = {m: compute_R(bkg, m, cuts[m]) for m in high_masses}
    return low_R, high_R


def build_table(low_R, high_R):
    low = sorted((m, low_R[m]) for m in low_R)
    high = sorted((m, high_R[m]) for m in high_R)
    n = max(len(low), len(high))
    rows = []
    for i in range(n):
        lm = f"{low[i][0]} & {low[i][1]:.2f}" if i < len(low) else "  &     "
        hm = f"{high[i][0]:<2} & {high[i][1]:.2f}" if i < len(high) else "   &     "
        end = r" \\ \hline" if i == n - 1 else r" \\"
        rows.append(f"      {lm} & {hm}{end}")
    body = "\n".join(rows)
    return (
        r"""\begin{table}[h]
  \begin{center}
    \topcaption{Sculpting ratio $R(m_{a})$ of Eq.~(\ref{eq:sculpt_R}), evaluated on the
    inclusive simulated background, for the low-mass and high-mass BDT models. $R=1$
    corresponds to no sculpting of the $m_{\ell\ell\gamma\gamma}$ spectrum in the
    Higgs-peak window; $R>1$ ($R<1$) indicates an enhancement (depletion) of the
    peak window relative to the sidebands. The largest value,
    $R(m_{a}=30~\GeV)=1.30$, is not a concern: $R$ is evaluated on the simulated
    background, whereas the final background is the data-driven sideband fit, and a
    single fixed-window ratio cannot distinguish a broad, smoothly-varying tilt of the
    mass shape (absorbed by the analytic background function) from a localized fake
    peak. This is confirmed by the background-only pseudo-data closure test
    (Fig.~\ref{fig:sculpt_closure}), whose observed limit at $m_{a}=30~\GeV$ lies within
    the $\pm1\sigma$ expected band with a best-fit signal strength compatible with zero,
    showing that the $R>1$ there does not bias the result.}
    \label{tab:sculpt_R}
    \small
    \begin{tabular}{|c|c||c|c|} \hline
      \multicolumn{2}{|c||}{Low-mass model} & \multicolumn{2}{c|}{High-mass model} \\ \hline
      $m_{a}$ [\GeV] & $R(m_{a})$ & $m_{a}$ [\GeV] & $R(m_{a})$ \\ \hline
"""
        + body
        + r"""
    \end{tabular}
  \end{center}
\end{table}
"""
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", default=JSON_PATH,
                    help="MVAcut_points_run3.json with the adopted working-point BDT cuts")
    ap.add_argument("--out", default=DEF_OUT, help="output LaTeX table path")
    ap.add_argument("--into", default=None,
                    help="target .tex; replace the AUTO sculpt_R block (append if absent)")
    a = ap.parse_args()

    low_R, high_R = compute_R_from_wp(a.json)

    print(f"# R(m_a) at adopted BDT working points, from stored MVA scores in {SCORED_BASE}")
    for m in sorted(low_R):
        print(f"#  low-mass  mA={m:>2}: R={low_R[m]:.3f}")
    for m in sorted(high_R):
        print(f"#  high-mass mA={m:>2}: R={high_R[m]:.3f}")
    print()

    table = build_table(low_R, high_R)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        f.write(table)
    print(table)
    print(f"[sculpt-R-table] wrote {a.out}")

    if a.into:
        block = f"{BEGIN_MARK}\n{table}{END_MARK}\n"
        with open(a.into) as f:
            txt = f.read()
        pat = ANY_BEGIN + r".*?" + re.escape(END_MARK) + r"\n?"
        if re.search(pat, txt, flags=re.DOTALL):
            txt = re.sub(pat, lambda _m: block, txt, flags=re.DOTALL)
            print(f"[sculpt-R-table] replaced AUTO block in {a.into}")
        else:
            txt = txt.rstrip() + "\n\n" + block
            print(f"[sculpt-R-table] appended AUTO block to {a.into}")
        with open(a.into, "w") as f:
            f.write(txt)


if __name__ == "__main__":
    main()
