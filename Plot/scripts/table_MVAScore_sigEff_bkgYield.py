import os
import json

mAs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
# mAs = [5, 15, 30]
# 新增：輸出 LaTeX 表格
output_latex = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/latexTable_MVAcut_eff_bkg.txt"

rows = []
for mA in mAs:
    input_file_name = f"/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/optimize_run3UL/nCat1_all_M{mA}.json"
    if not os.path.exists(input_file_name):
        continue
    try:
        with open(input_file_name, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError:
        continue

    boundaries = data.get("boundaries", [])
    if not boundaries:
        continue
    mvacut = float(boundaries[0])

    details = data.get("details", {}) or {}
    signal_cut = details.get("signal_cut")
    signal_total = details.get("signal_total")
    # 優先用 details.smoothed_background，否則回退到頂層
    bkg = details.get("smoothed_background", data.get("smoothed_background"))

    if signal_cut is None or signal_total is None or bkg is None:
        continue
    if signal_total == 0:
        continue

    eff = signal_cut / signal_total
    rows.append(f"{mA} & {mvacut:.3f} & {eff:.3f} & {bkg:.3f} \\\\")

# 輸出 LaTeX 表格到 txt
os.makedirs(os.path.dirname(output_latex), exist_ok=True)
latex_table = (
    "\\begin{table}[h]\n"
    "\\begin{center}\n"
    "\\topcaption{Summary of the minimum BDT output value, signal efficiency, and smoothed MC background yields with respect to a selection on the BDT output for each nominal signal hypothesis.}\n"
    "\\label{tab:boundary}\n"
    "\\small\n"
    "\\begin{tabular}{|c|c|c|c|}\n"
    "\\hline\n"
    "$m_a$ [GeV] & \\shortstack[c]{Min. BDT \\\\ output value} & \\shortstack[c]{Signal efficiency \\\\ w.r.t. BDT selection} & \\shortstack[c]{Smoothed MC background yields \\\\ w.r.t. BDT selection}\\\\\n"
    "\\hline\n"
    + ("\n".join(rows) if rows else "")
    + ("\n" if rows else "")
    + "\\hline\n"
    "\\end{tabular}\n"
    "\\end{center}\n"
    "\\end{table}\n"
)

with open(output_latex, 'w') as f:
    f.write(latex_table)
