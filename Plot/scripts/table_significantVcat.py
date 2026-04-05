import os
import json

mAs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
# mAs = [5, 15, 30]
cats = [1, 2]
# 準備 JSON 輸出結構
results = []
# 新增：準備 LaTeX 表格輸出與彙整 cat1/2 significance 的對照表
output_latex = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/latexTable_significanceVcat.txt"
sig_map = {}

# 以 JSON 解析輸入檔，從 boundaries 與 significance 取值
for cat in cats:
    for mA in mAs:
        input_file_name = f"/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/optimize_run3UL/nCat{cat}_all_M{mA}.json"
        if not os.path.exists(input_file_name):
            continue
        try:
            with open(input_file_name, 'r') as f:
                data = json.load(f)
        except json.JSONDecodeError:
            continue

        boundaries = data.get("boundaries", [])
        if not boundaries or data.get("significance") is None:
            continue

        mvacut = float(boundaries[0])
        significance = float(data["significance"])

        results.append({
            "mA": mA,
            "MVAcut": round(mvacut, 3),
            "Significance": round(significance, 3)
        })
        # 新增：彙整 cat1/2 significance
        sig_map.setdefault(mA, {})[cat] = significance

# 新增：輸出 LaTeX 表格到 txt
os.makedirs(os.path.dirname(output_latex), exist_ok=True)
rows = []
for mA in sorted(sig_map.keys()):
    s1 = sig_map[mA].get(1)
    s2 = sig_map[mA].get(2)
    if s1 is not None and s2 is not None:
        improvement = (s2 - s1)/s1 * 100 if s1 != 0 else 0
        rows.append(f"{mA} & {s1:.3f} & {s2:.3f} & +{improvement:.1f}\% \\\\")
    elif s1 is not None:
        rows.append(f"{mA} & {s1:.3f} & -- & -- \\\\")
    elif s2 is not None:
        rows.append(f"{mA} & -- & {s2:.3f} & -- \\\\")
    # 若都沒有則略過

latex_table = (
    "\\begin{table}[h]\n"
    "\\begin{center}\n"
    "\\topcaption{Approximate mean significance for different numbers of categories, and its improvement as a function of the number of categories.}\n"
    "\\label{tab:significanceVcat}\n"
    "\\begin{tabular}{cccc}\n"
    "\\hline\n"
    "$m_a$ [GeV] & \multicolumn{2}{c}{Significance} & Improvement \\\ \n"
    "\\cline{2-3}\n"
    "& 1 cat. & 2 cats. &  \\\\\n"
    "\\hline\n"
    + ("\n".join(rows) if rows else "")
    + ("\n" if rows else "")
    + "\\hline\n"
    "\\end{tabular}\n"
    "\\end{center}\n"
    "\\end{table}"
)

with open(output_latex, 'w') as f:
    f.write(latex_table)
