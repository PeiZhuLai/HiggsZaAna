#!/usr/bin/env python3
"""Read Plot/output/sculpt_R_run3.json (written by scan_score_R_significance.py) and emit the
AN sculpting-ratio R(m_a) LaTeX table for the low-mass + high-mass BDT models.

Prints the table to stdout and writes it to Plot/output/sculpt_R_table.tex. Optionally
(--into <tex>) replaces the table block (between the BEGIN/END markers) inside a target .tex.
"""
import json, argparse, os, re

DEF_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sculpt_R_run3.json"
DEF_OUT  = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sculpt_R_table.tex"
BEGIN_MARK = "% >>> AUTO sculpt_R table (make_sculpt_R_table.py) >>>"
END_MARK   = "% <<< AUTO sculpt_R table <<<"


def build_table(jpath):
    d = json.load(open(jpath))
    low  = sorted((int(k), float(v)) for k, v in d["low_mass"].items())
    high = sorted((int(k), float(v)) for k, v in d["high_mass"].items())
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
    peak window relative to the sidebands.}
    \label{tab:sculpt_R}
    \small
    \begin{tabular}{|c|c||c|c|} \hline
      \multicolumn{2}{|c||}{Low-mass model} & \multicolumn{2}{c|}{High-mass model} \\ \hline
      $m_{a}$ [\GeV] & $R(m_{a})$ & $m_{a}$ [\GeV] & $R(m_{a})$ \\ \hline
""" + body + r"""
    \end{tabular}
  \end{center}
\end{table}
"""
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", default=DEF_JSON)
    ap.add_argument("--out",  default=DEF_OUT)
    ap.add_argument("--into", default=None, help="target .tex; replace block between AUTO markers (append if absent)")
    a = ap.parse_args()

    table = build_table(a.json)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        f.write(table)
    print(table)
    print(f"[sculpt-R-table] wrote {a.out}")

    if a.into:
        block = f"{BEGIN_MARK}\n{table}{END_MARK}\n"
        with open(a.into) as f:
            txt = f.read()
        if BEGIN_MARK in txt and END_MARK in txt:
            # NB: replacement passed as a function so LaTeX backslashes in `block` are not
            # interpreted as regex escapes (re.sub with a string replacement would choke on \e etc.).
            txt = re.sub(re.escape(BEGIN_MARK) + r".*?" + re.escape(END_MARK) + r"\n?",
                         lambda _m: block, txt, flags=re.DOTALL)
            print(f"[sculpt-R-table] replaced AUTO block in {a.into}")
        else:
            txt = txt.rstrip() + "\n\n" + block
            print(f"[sculpt-R-table] appended AUTO block to {a.into}")
        with open(a.into, "w") as f:
            f.write(txt)


if __name__ == "__main__":
    main()
