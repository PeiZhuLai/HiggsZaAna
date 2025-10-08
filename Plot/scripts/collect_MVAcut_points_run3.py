import os
import json

# mAs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
mAs = [5, 15, 30]

# 準備 JSON 輸出結構
results = []
output_json = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

# Loop through each component in the list and extract the value
for mA in mAs:
    input_file_name = f"/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/optimize_run3UL/categorize_all_M{mA}.txt"
    with open(input_file_name, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 8:
                continue
            # 擷取 MVAcut 與 Significance（沿用原本欄位位置）
            BDT_value = float(parts[4])
            significance = float(parts[7])
            results.append({
                "mA": mA,
                "MVAcut": round(BDT_value, 3),
                "Significance": round(significance, 3)
            })
            break  # 只取第一筆有效資料

# 寫入單一 JSON 檔案
with open(output_json, 'w') as out:
    json.dump({"results": results}, out, indent=2)
