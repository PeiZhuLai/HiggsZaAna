#!/usr/bin/env python3
"""Summary of the high-mass minimal-change scan (minscan_high_win/config.jsonl): gate violation and objective with the
variable-search definitions (reference = stored production high-mass scores), Z/Z_prod range over the 11 anchors, R
range, and the closure-proxy pulls at mA25 and mA30 (where the production closure shows r_hat ~ +0.03).
Env: higgs-alp-ana. Output: minscan_high_win/summary.txt"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "bdt_varsearch"))
import driver_condor as dc  # noqa: E402
R = "/eos/home-p/pelai/HZa/bdt_varsearch_20261013/results/minscan_high_win_20261017"
cfgs = [json.loads(l) for l in open(os.path.join(HERE, "config.jsonl")) if l.strip()]
res = {c["id"]: json.load(open(os.path.join(R, c["id"] + ".json"))) for c in cfgs if os.path.exists(os.path.join(R, c["id"] + ".json"))}
zref = {m: dc.z_fit(r) for m, r in res["mh_ref_stored"]["per_mass"].items()}
rows = []
for c in cfgs:
    r = res.get(c["id"])
    if r is None:
        rows.append((9e9, 0, "%-55s missing" % c["id"])); continue
    v, o = dc.gate_violation(r), dc.objective(r, zref)
    pm = r["per_mass"]
    zs = [dc.z_fit(pm[m]) / zref[m] for m in pm]
    Rs = [pm[m]["R"] for m in pm]
    pulls = lambda m: " ".join("%+.1f" % p for p in pm[m]["A_pulls"])
    worst = max(pm, key=lambda m: dc.gate_violation({"per_mass": {m: pm[m]}}))
    rows.append((v if v > dc.TOL else -1, -o, "%-55s %6.2f %6.3f | Z/Zp %.2f-%.2f | R %.2f-%.2f | pulls mA25 %s mA30 %s | worst mA%s"
                 % (c["id"], v, o, min(zs), max(zs), min(Rs), max(Rs), pulls("25"), pulls("30"), worst)))
L = ["%d/%d results" % (len(res), len(cfgs)), "%-55s %6s %6s" % ("candidate (sorted: passing by objective, then by violation)", "viol", "obj")]
L += [x[2] for x in sorted(rows)]
txt = "\n".join(L)
open(os.path.join(HERE, "summary.txt"), "w").write(txt + "\n")
print(txt)
