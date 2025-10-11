import os
import json

# mAs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
mAs = [5, 15, 30]

# 準備 JSON 輸出結構
results = []
output_json = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"

# 以 JSON 解析輸入檔，從 boundaries 與 significance 取值
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
    if not boundaries or data.get("significance") is None:
        continue

    mvacut = float(boundaries[0])
    significance = float(data["significance"])

    results.append({
        "mA": mA,
        "MVAcut": round(mvacut, 3),
        "Significance": round(significance, 3)
    })

# 寫入單一 JSON 檔案
os.makedirs(os.path.dirname(output_json), exist_ok=True)
with open(output_json, 'w') as out:
    json.dump({"results": results}, out, indent=2)
