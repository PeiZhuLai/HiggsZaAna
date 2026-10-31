# smoke 驗收（run_phidzee2_20261031.sh 的 SHORT=1 產物）：
#   每個 task 目錄都有 summary（完成標記；fpo=1 時 0 事件通過不寫 parquet 是正常）；
#   有 merged_nominal.parquet 的 MC：photon ID 權重欄（含 stat/syst 拆分）存在、weight_central 無 NaN；
#   control MC 另分 z_ee / z_mumu 兩區看 NaN（6 月那版 z_mumu 區 99.7% NaN，electron SF 對空集合回 NaN）。
import glob, os, sys
import numpy as np
import pyarrow.parquet as pq

BASES = ["/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_phidzee2_smoke"]  # phidzee2 只重產 sig
bad = 0
for base in BASES:
    for unit in sorted(glob.glob(base + "/*")):
        for task in sorted(d for d in glob.glob(unit + "/*") if os.path.isdir(d)):
            name = task.replace(base + "/", "")
            summ = [f for f in glob.glob(task + "/**/*summary*.json", recursive=True)]
            merged = task + "/merged_nominal.parquet"
            line = "%-60s summary=%d" % (name, len(summ))
            if not summ:
                bad += 1; print("NOSUMMARY " + line); continue
            if not os.path.exists(merged):
                print("ok(0 evt?) " + line + " no merged parquet"); continue
            sch = pq.read_schema(merged).names
            cols = [c for c in ["weight_central", "z_ee", "z_mumu"] if c in sch]
            phid = sorted(c for c in sch if "zaid" in c or "photon_id_sf" in c)  # control MC: weight key photon_id_sf -> photon_zaid_sf
            t = pq.read_table(merged, columns=cols + [c for c in phid if c.endswith("central") or "_central" in c][:3])
            n = t.num_rows
            msg = line + " rows=%d" % n
            is_data = "/data" in name or name.split("/")[1].startswith("Data")
            if "weight_central" in cols and not is_data:
                w = t.column("weight_central").to_numpy(zero_copy_only=False)
                nan = np.isnan(w).sum(); msg += " NaN(w)=%d" % nan
                if nan: bad += 1
                for reg in ("z_ee", "z_mumu"):
                    if reg in cols:
                        m = t.column(reg).to_numpy(zero_copy_only=False).astype(bool)
                        msg += " %s:%d/NaN %d" % (reg, m.sum(), np.isnan(w[m]).sum())
            if not is_data:
                has_split = any("zaid_sf_syst" in c for c in phid)
                msg += " phid_cols=%d split=%s" % (len(phid), has_split)
                for c in t.column_names:
                    if "zaid" in c or "photon_id_sf" in c:
                        v = t.column(c).to_numpy(zero_copy_only=False)
                        msg += " <%s>=%.3f" % (c.replace("weight_", ""), np.nanmean(v))
                if not phid: bad += 1; msg = "NOPHID " + msg
            print(msg)
print("problems", bad); sys.exit(1 if bad else 0)
