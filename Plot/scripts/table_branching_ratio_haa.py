#!/usr/bin/env python3
r"""Compute the exotic Higgs branching fractions Br(H->Za) and Br(H->aa) as a function of
the ALP mass, following the ALP effective-field-theory results of
Bauer, Neubert, Thamm, JHEP 12 (2017) 044 [arXiv:1708.00443] (AN ref. \cite{Bauer:2017ris}).

Motivation (L3 convener comment, AN L92): the ALP model does *not* predict H->Za
exclusively. H->Za and H->aa are two independent exotic Higgs decays governed by
*different* Wilson coefficients, so there is no model-independent relation fixing the
ratio of their branching fractions. This script tabulates both, side by side, over the
Run 3 signal mass grid, at representative couplings, and emits the AN LaTeX table
(latexTable_BrHaa.txt) that is \input in Sec-02-00-ALPModel.tex.

Partial widths used:
  H->Za  (Eq. 5.6 of the paper; already in the AN):
    Gamma(H->Za) = m_H^3 / (16 pi) * (cZh/Lambda)^2
                   * lambda^{3/2}(m_Z^2/m_H^2, m_a^2/m_H^2),
    with lambda(x,y) = (1 - x - y)^2 - 4 x y.

  H->aa  (Eq. 5.12 of the paper; dimension-6 Higgs-portal coupling C_ah^eff):
    Gamma(H->aa) = v^2 m_H^3 / (32 pi) * (C_ah^eff/Lambda^2)^2
                   * (1 - 2 m_a^2/m_H^2)^2 * sqrt(1 - 4 m_a^2/m_H^2).

Branching fraction: Br(H->X) = Gamma(H->X) / (Gamma(H->X) + Gamma_SM),
with Gamma_SM = 4.1e-3 GeV the SM Higgs total width.

Reference couplings (defaults):
  cZh/Lambda        = 2e-5   GeV^-1   (as in the existing Br(H->Za) table / AN example)
  C_ah^eff/Lambda^2 = 0.62e-6 GeV^-2  = 0.62 TeV^-2
                      (chosen so Br(H->aa) = 10% for a light ALP, per the paper).

Usage:
  python3 table_branching_ratio_haa.py                # print + write latexTable_BrHaa.txt
  python3 table_branching_ratio_haa.py --cahL2 1.34e-6 --czhL 2e-5
  python3 table_branching_ratio_haa.py --into <Sec-02-00-ALPModel.tex>   # optional in-place
"""
import argparse
import math
import os
import re

# ---- physical constants (GeV) --------------------------------------------------------
V_EW = 246.0        # Higgs field vacuum expectation value
M_H = 125.38        # Higgs boson mass (analysis value)
M_Z = 91.1876       # Z boson mass
GAMMA_SM = 4.1e-3   # SM Higgs total width

# ---- reference couplings -------------------------------------------------------------
CZH_OVER_LAMBDA = 2.0e-5     # c^eff_Zh / Lambda        [GeV^-1]
CAH_OVER_LAMBDA2 = 0.62e-6   # C^eff_ah / Lambda^2      [GeV^-2]  (= 0.62 TeV^-2)

# ALP mass grid used for the Run 3 signal (GeV)
MASS_GRID = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]

DEF_OUT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/latexTable_BrHaa.txt"
BEGIN_MARK = "% >>> AUTO BrHaa table (table_branching_ratio_haa.py) >>>"
END_MARK = "% <<< AUTO BrHaa table <<<"


def _lam(x, y):
    """Kaellen-type function lambda(x, y) = (1 - x - y)^2 - 4 x y."""
    return (1.0 - x - y) ** 2 - 4.0 * x * y


def gamma_hZa(m_a, czh_over_lambda=CZH_OVER_LAMBDA):
    """Partial width Gamma(H->Za) in GeV. Requires m_a < m_H - m_Z (~34.2 GeV)."""
    if m_a >= M_H - M_Z:
        return 0.0
    lam = _lam(M_Z ** 2 / M_H ** 2, m_a ** 2 / M_H ** 2)
    if lam <= 0.0:
        return 0.0
    return M_H ** 3 / (16.0 * math.pi) * czh_over_lambda ** 2 * lam ** 1.5


def gamma_haa(m_a, cah_over_lambda2=CAH_OVER_LAMBDA2):
    """Partial width Gamma(H->aa) in GeV. Requires m_a < m_H / 2 (~62.7 GeV)."""
    beta2 = 1.0 - 4.0 * m_a ** 2 / M_H ** 2
    if beta2 <= 0.0:
        return 0.0
    phase = (1.0 - 2.0 * m_a ** 2 / M_H ** 2) ** 2 * math.sqrt(beta2)
    return V_EW ** 2 * M_H ** 3 / (32.0 * math.pi) * cah_over_lambda2 ** 2 * phase


def branching_fraction(gamma):
    """Br(H->X) = Gamma(H->X) / (Gamma(H->X) + Gamma_SM)."""
    return gamma / (gamma + GAMMA_SM)


def _fmt(x):
    """Format a positive number as 'm.mm \\times 10^{e}' for LaTeX."""
    if x <= 0.0:
        return "0"
    exp = math.floor(math.log10(x))
    mant = x / 10.0 ** exp
    # guard against mantissa rounding up to 10.00 (e.g. 9.997e-2 -> 1.00e-1)
    if round(mant, 2) >= 10.0:
        mant /= 10.0
        exp += 1
    return rf"{mant:.2f} \times 10^{{{exp}}}"


def build_table(czh_over_lambda=CZH_OVER_LAMBDA, cah_over_lambda2=CAH_OVER_LAMBDA2):
    rows = []
    for m_a in MASS_GRID:
        br_za = branching_fraction(gamma_hZa(m_a, czh_over_lambda))
        br_aa = branching_fraction(gamma_haa(m_a, cah_over_lambda2))
        rows.append(rf"{m_a} & ${_fmt(br_za)}$ & ${_fmt(br_aa)}$ \\")
    body = "\n".join(rows)
    return (
        r"""\begin{table}[htbp]
  \centering
  \caption{Comparison of the exotic Higgs branching fractions \Br(\hToZa) and \Br(\hToaa),
  computed with the ALP effective field theory of Ref.~\cite{Bauer:2017ris}, as a function
  of the ALP mass over the Run~3 signal mass grid. The two channels are governed by
  independent Wilson coefficients (\cZh and \cah, respectively); the representative values
  $\cZhLambda = 2 \times 10^{-5}\GeV^{-1}$ and $\cahLambda = 6.2 \times 10^{-7}\GeV^{-2}$
  ($= 0.62\,\text{TeV}^{-2}$, chosen to give a $10\%$ \hToaa branching fraction for a light
  ALP~\cite{Bauer:2017ris}) are used. The two rates carry no fixed ratio: it depends on the
  underlying ultraviolet completion.}
  \label{tab:BrHaa}
  \begin{tabular}{|c|c|c|}
    \hline
    $m_{a}$ [GeV] & $\Br(\hToZa)$ & $\Br(\hToaa)$ \\
    \hline\hline
"""
        + body
        + r"""
    \hline
  \end{tabular}
\end{table}
"""
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--czhL", type=float, default=CZH_OVER_LAMBDA,
                    help="c^eff_Zh / Lambda [GeV^-1] (default %(default)g)")
    ap.add_argument("--cahL2", type=float, default=CAH_OVER_LAMBDA2,
                    help="C^eff_ah / Lambda^2 [GeV^-2] (default %(default)g)")
    ap.add_argument("--out", default=DEF_OUT, help="output LaTeX table path")
    ap.add_argument("--into", default=None,
                    help="target .tex; replace block between AUTO markers (append if absent)")
    a = ap.parse_args()

    table = build_table(a.czhL, a.cahL2)

    # human-readable console summary
    print(f"# Gamma_SM = {GAMMA_SM:.3e} GeV, m_H = {M_H} GeV, v = {V_EW} GeV")
    print(f"# cZh/Lambda = {a.czhL:.3e} GeV^-1, C_ah^eff/Lambda^2 = {a.cahL2:.3e} GeV^-2")
    print(f"# {'m_a[GeV]':>8}  {'Br(H->Za)':>12}  {'Br(H->aa)':>12}")
    for m_a in MASS_GRID:
        br_za = branching_fraction(gamma_hZa(m_a, a.czhL))
        br_aa = branching_fraction(gamma_haa(m_a, a.cahL2))
        print(f"# {m_a:>8}  {br_za:>12.3e}  {br_aa:>12.3e}")
    print()

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        f.write(table)
    print(table)
    print(f"[BrHaa-table] wrote {a.out}")

    if a.into:
        block = f"{BEGIN_MARK}\n{table}{END_MARK}\n"
        with open(a.into) as f:
            txt = f.read()
        if BEGIN_MARK in txt and END_MARK in txt:
            txt = re.sub(re.escape(BEGIN_MARK) + r".*?" + re.escape(END_MARK) + r"\n?",
                         lambda _m: block, txt, flags=re.DOTALL)
            print(f"[BrHaa-table] replaced AUTO block in {a.into}")
        else:
            txt = txt.rstrip() + "\n\n" + block
            print(f"[BrHaa-table] appended AUTO block to {a.into}")
        with open(a.into, "w") as f:
            f.write(txt)


if __name__ == "__main__":
    main()
