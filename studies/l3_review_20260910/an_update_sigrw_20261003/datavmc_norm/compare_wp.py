#!/usr/bin/env python3
"""Working-point comparison after the signal-normalization rerun of ALP_Optimization.

adopted   : the analysis configuration Plot/output/MVAcut_points_run3.json (not modified)
candidate : what chain_dyveto_flashgg.sh [2/9] would write = collect (nCat1_all_M*.json boundaries)
            followed by the R scan's chosen_wp for mA1-3 (FIXED_WP for mA2/3, JSON for mA1)
pure opt  : significance maximum without the low-mass override (ALP_Optimization copy with the
            override disabled), old merged histograms vs new merged histograms
Also: nCat1 / nCat2 significance old vs new (old = backup of optimize_run3UL)."""
import json, os, re, sys
W = os.path.dirname(os.path.abspath(__file__))
PLOT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot"
OPT = PLOT + "/plots/optimize_run3UL"
OPTBAK = PLOT + "/plots_variants/optimize_run3UL_stale_preSigRwNorm_20261003"
MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]

adopted = {int(r["mA"]): float(r["MVAcut"]) for r in json.load(open(PLOT + "/output/MVAcut_points_run3.json"))["results"]}
cand = {int(r["mA"]): float(r["MVAcut"]) for r in json.load(open(W + "/MVAcut_points_run3_candidate.json"))["results"]}
# R scan (copy, no --write-json): chosen working point for mA1-3 (columns: mA R1cut Z@R1 Rmin-cut Rmin Z@Rmin WPcut R@WP Z@WP)
scan_wp = {}
for line in open(W + "/wp_scan.log"):
    p = line.split()
    if len(p) >= 9 and p[0] in ("1", "2", "3"):
        scan_wp[int(p[0])] = float(p[6])
for m, c in scan_wp.items():
    cand[m] = c


def jl(path):
    return json.load(open(path)) if os.path.exists(path) else None


def pure(d, m):
    j = jl("%s/nCat1_all_M%d.json" % (d, m))
    return (j["boundaries"][0], j["significance"]) if j else (None, None)


print("%4s %9s %10s %7s | %11s %11s | %9s %9s %8s | %9s %9s" % (
    "mA", "adopted", "candidate", "change", "pureOpt_old", "pureOpt_new", "Z1_old", "Z1_new", "Z1 n/o", "Z2_old", "Z2_new"))
changed = []
for m in MASSES:
    a, c = adopted.get(m), cand.get(m)
    ch = "SAME" if (a is not None and c is not None and abs(a - c) < 1e-9) else "CHANGED"
    if ch != "SAME":
        changed.append(m)
    po, pn = pure(W + "/opt_nooverride_old", m), pure(W + "/opt_nooverride_new", m)
    o1, n1 = jl("%s/nCat1_all_M%d.json" % (OPTBAK, m)), jl("%s/nCat1_all_M%d.json" % (OPT, m))
    o2, n2 = jl("%s/nCat2_all_M%d.json" % (OPTBAK, m)), jl("%s/nCat2_all_M%d.json" % (OPT, m))
    z1o, z1n = o1["significance"], n1["significance"]
    print("%4d %9.3f %10.3f %7s | %11s %11s | %9.3f %9.3f %8.4f | %9.3f %9.3f" % (
        m, a, c, ch, po[0], pn[0], z1o, z1n, z1n / z1o, o2["significance"], n2["significance"]))
print()
print("CHANGED working points: %s" % (changed if changed else "none"))
sys.exit(1 if changed else 0)
