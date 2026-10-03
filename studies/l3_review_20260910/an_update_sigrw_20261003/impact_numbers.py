"""Expected-impact summary (Sec 9.3): best-fit r on Asimov and the rank of CMS_hza_mva_reweight / QCD scale.
usage: impact_numbers.py <impacts_dir>   Input: <dir>/<m>_impacts.json. Output: stdout"""
import json, sys
D = sys.argv[1]; rs = []; top2 = 0; cnt = {}
for m in range(1, 31):
    d = json.load(open("%s/%d_impacts.json" % (D, m)))
    r = d["POIs"][0]["fit"][1]; rs.append(r)
    ps = sorted(d["params"], key=lambda p: -abs(p["impact_r"]))
    names = [p["name"] for p in ps]
    rk = names.index("CMS_hza_mva_reweight") + 1 if "CMS_hza_mva_reweight" in names else None
    top2 += rk is not None and rk <= 2
    for n in names[:2]: cnt[n] = cnt.get(n, 0) + 1
    print("%2d r=%.4f  top3: %s  | mva_reweight rank %s impact %.4f" % (m, r, ", ".join("%s(%.4f)" % (p["name"], p["impact_r"]) for p in ps[:3]), rk, [p["impact_r"] for p in ps if p["name"] == "CMS_hza_mva_reweight"][0]))
print("best-fit r range: %.4f - %.4f" % (min(rs), max(rs)))
print("mva_reweight in top 2 at %d/30 points" % top2)
print("top-2 counts:", sorted(cnt.items(), key=lambda x: -x[1]))
