#!/usr/bin/env python3
"""Summary of the minimal-change low-mass scan (minscan/config.jsonl, results in
/eos/home-p/pelai/HZa/bdt_varsearch_20261013/results/minscan_low_20261014): gate violation and objective with the
variable-search definitions (driver_condor.gate_violation / z_fit / objective, reference = the stored production
scores), plus per-mass R, R_lo, closure-proxy pulls, data GOF p, and for mA1 the double-peak numbers.
Env: higgs-alp-ana. Output: minscan/summary.txt"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "bdt_varsearch"))
import driver_condor as dc  # noqa: E402
R = "/eos/home-p/pelai/HZa/bdt_varsearch_20261013/results/minscan_low_20261014"
cfgs = [json.loads(l) for l in open(os.path.join(HERE, "config.jsonl")) if l.strip()]
res = {c["id"]: json.load(open(os.path.join(R, c["id"] + ".json"))) for c in cfgs if os.path.exists(os.path.join(R, c["id"] + ".json"))}
ref = res["ms_ref_stored"]
zref = {m: dc.z_fit(r) for m, r in ref["per_mass"].items()}
L = ["%d/%d results" % (len(res), len(cfgs)), ""]
L.append("%-45s %6s %6s | %-58s | mA1 win/ovl (presel %.2f/%.2f)" % ("candidate", "viol", "obj", "mA1/2/3: Z/Zp  R  R_lo  maxApull  Bp",
         ref["per_mass"]["1"]["presel_win_frac"], ref["per_mass"]["1"]["presel_overlap"]))
for c in cfgs:
    r = res.get(c["id"])
    if r is None:
        L.append("%-45s missing" % c["id"]); continue
    v, o = dc.gate_violation(r), dc.objective(r, zref)
    pm = []
    for m in ("1", "2", "3"):
        x = r["per_mass"][m]
        pm.append("%.2f %.2f %.2f %+.1f %.3f" % (dc.z_fit(x) / zref[m], x["R"], x["R_lo"], max(x["A_pulls"], key=abs), x["B_p"]))
    x1 = r["per_mass"]["1"]
    L.append("%-45s %6.2f %6.2f | %s | %.2f/%.2f" % (c["id"], v, o, " ; ".join(pm), x1["win_frac"], x1["overlap"]))
txt = "\n".join(L)
open(os.path.join(HERE, "summary.txt"), "w").write(txt + "\n")
print(txt)
