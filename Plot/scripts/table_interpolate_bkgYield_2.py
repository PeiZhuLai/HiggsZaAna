import os
import json

output_latex = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/latexTable_bkg_interpolated_yields.txt"

# 讀入單一 JSON，解析 results
input_file_name = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/bkg_interpolated_yields_run3.json"
rows = []
if os.path.exists(input_file_name):
    try:
        with open(input_file_name, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError:
        data = {}
    results = data.get("results", []) or []
    # 依 mA 排序並輸出四欄
    for item in sorted(results, key=lambda x: x.get("mA", 0)):
        mA = item.get("mA")
        lower = item.get("bkgTotal_lower")
        nearest = item.get("bkgTotal_nearest")
        upper = item.get("bkgTotal_upper")
        if mA is None or lower is None or nearest is None or upper is None:
            continue
        rows.append(f"{mA} & {lower:.3f} & {nearest:.3f} & {upper:.3f} \\\\")

os.makedirs(os.path.dirname(output_latex), exist_ok=True)
latex_table = (
    "\\begin{table}[h]\n"
    "\\centering\n"
    "\\begin{tabular}{cccc}\n"
    "\\hline\n"
    "$m_a$ [GeV] & Left & Nominal & Right \\\\\n"
    "\\hline\n"
    + ("\n".join(rows) if rows else "")
    + ("\n" if rows else "")
    + "\\hline\n"
    "\\end{tabular}\n"
    "\\caption{Summary of the smoothed MC background yields with respect to different BDT threshold for each intermediate signal hypothesis. Nominal means the threshold is selected at the same point as the closest nominal mass point. Left and right mean the threshold is selected from the smaller ALP mass threshold and larger ALP mass threshold sides of the nominal mass point, respectively.}\n"
    "\\end{table}\n"
)

with open(output_latex, 'w') as f:
    f.write(latex_table)
