"""Summary numbers of the background-only pseudo-data closure (Sec 6.4 sec:sculpt_closure).
usage: closure_numbers.py <closure_results.json>   Output: stdout"""
import json, sys
d = json.load(open(sys.argv[1]))
in1 = []; in2 = []; down = []; up = []; rows = []
for m in range(1, 31):
    v = d[str(m)]; L = v["limit"]; o = L["obs"]
    i1 = L["exp_m1"] <= o <= L["exp_p1"]; i2 = L["exp_m2"] <= o <= L["exp_p2"]
    in1 += [m] if i1 else []; in2 += [m] if i2 else []
    if not i1: (down if o < L["exp_m1"] else up).append(m)
    rows.append((v["signif"], m, v["rhat"]["r"], v["rhat"]["rLoErr"], v["rhat"]["rHiErr"]))
    print("%2d obs=%.4f exp=%.4f [%.4f,%.4f] Z=%.2f rhat=%+.3f -%.3f +%.3f %s" % (m, o, L["exp_med"], L["exp_m1"], L["exp_p1"], v["signif"], v["rhat"]["r"], v["rhat"]["rLoErr"], v["rhat"]["rHiErr"], "" if i1 else ("DOWN" if o < L["exp_m1"] else "UP")))
print("within 1sigma: %d, within 2sigma: %d" % (len(in1), len(in2)))
print("outside 1 sigma, stronger (down):", down, " weaker (up):", up)
rows.sort(reverse=True)
print("largest Z:", ", ".join("mA%d Z=%.2f rhat=%+.3f" % (m, z, r) for z, m, r, _, _ in rows[:5]))
print("max |rhat| for mA>=3:", max((abs(r), m) for z, m, r, _, _ in rows if m >= 3), " excluding 5:", max((abs(r), m) for z, m, r, _, _ in rows if m >= 3 and m != 5))
print("mA1,2,3 Z:", [round(d[str(m)]["signif"], 2) for m in (1, 2, 3)], "rhat mA1:", d["1"]["rhat"])
