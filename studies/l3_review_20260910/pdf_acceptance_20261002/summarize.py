#!/usr/bin/env python
"""Collect results/<sample>.json into a per-(m_a, era) table and a summary.

Columns (all symmetric Hessian, NNPDF31_nnlo_as_0118_nf_4_mc_hessian, 100 eigenvectors):
  acc_wp      acceptance uncertainty for the datacard events (BDT test split, above the WP)  <- main
  stat        bootstrap standard deviation of acc_wp
  acc_presel  same at preselection (selection used by the QCD-scale study)
  acc_wpincl  same cut applied to all BDT splits (3.3x statistics)
  yield_wp    selected-yield variation (acceptance x cross section)
  xs          generated cross-section variation (denominator only)
Output: results/summary.json, results/summary.txt
"""
import glob, json, os
import numpy as np

W = os.path.dirname(os.path.abspath(__file__))
R = sorted(glob.glob(os.path.join(W, "results", "mA_M*.json")), key=lambda p: (int(os.path.basename(p).split("_")[1][1:]), p))
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
rows = []
for p in R:
    r = json.load(open(p))
    c = r["checks"]
    rows.append(dict(m_a=r["m_a"], era=r["era"],
                     acc_wp=r["wp"]["acc_unc"], stat=r["wp"]["acc_unc_boot_std"],
                     acc_presel=r["presel"]["acc_unc"], acc_wpincl=r["wp_incl"]["acc_unc"],
                     yield_wp=r["wp"]["yield_unc"], xs=r["wp"]["xs_unc"], n_wp=r["wp"]["n"],
                     wp_weight_diff=abs(r["wp_weight"]["acc_unc"] - r["wp"]["acc_unc"]),
                     unmatched=c["presel_unmatched"] + c["wp_unmatched"] + c["wp_incl_unmatched"],
                     m0dev=max(abs(c["presel_member0_minmax"][0] - 1), abs(c["presel_member0_minmax"][1] - 1)),
                     rmscheck=c["presel_parquetUp_minus_unit_vs_std_maxabs"],
                     wp_cut=c["wp_json_cut"], wp_min=c["wp_min_score_in_tree"],
                     wpdiff=c["wp_tree_vs_test_pass_cut_diff"]))
lines = ["%4s %-13s %8s %7s %10s %10s %9s %7s %6s" % ("m_a", "era", "acc_wp%", "stat%", "acc_presel", "acc_wpincl", "yield_wp%", "xs%", "n_wp")]
for x in rows:
    lines.append("%4d %-13s %8.3f %7.3f %10.3f %10.3f %9.3f %7.3f %6d" % (
        x["m_a"], x["era"], 100 * x["acc_wp"], 100 * x["stat"], 100 * x["acc_presel"], 100 * x["acc_wpincl"],
        100 * x["yield_wp"], 100 * x["xs"], x["n_wp"]))


def stats(k):
    v = np.array([x[k] for x in rows])
    return dict(mean=float(v.mean()), median=float(np.median(v)), max=float(v.max()), min=float(v.min()),
                argmax="m_a=%d %s" % (rows[int(v.argmax())]["m_a"], rows[int(v.argmax())]["era"]))


summ = {k: stats(k) for k in ["acc_wp", "acc_presel", "acc_wpincl", "yield_wp", "xs", "stat"]}
lines.append("")
lines.append("over %d (m_a, era) points" % len(rows))
for k, s in summ.items():
    lines.append("%-11s mean %.3f%%  median %.3f%%  min %.3f%%  max %.3f%% (%s)" % (k, 100 * s["mean"], 100 * s["median"], 100 * s["min"], 100 * s["max"], s["argmax"]))
# per mass, era-averaged
lines.append("")
lines.append("per m_a (mean over eras): acc_wp%")
for m in sorted(set(x["m_a"] for x in rows)):
    v = [x["acc_wp"] for x in rows if x["m_a"] == m]
    lines.append("  m_a=%2d  %.3f  (n_era=%d)" % (m, 100 * np.mean(v), len(v)))
lines.append("")
lines.append("validation: unmatched events total %d; max |member0-1| %.2e; max |parquet RMS - recomputed std| %.2e; "
             "max |acc(weight)-acc(weight_central)| %.2e; WP tree minus (test & score>cut) in [%d, %d]; "
             "min score in WP tree >= JSON cut everywhere: %s"
             % (sum(x["unmatched"] for x in rows), max(x["m0dev"] for x in rows), max(x["rmscheck"] for x in rows),
                max(x["wp_weight_diff"] for x in rows), min(x["wpdiff"] for x in rows), max(x["wpdiff"] for x in rows),
                all(x["wp_min"] > x["wp_cut"] for x in rows)))
txt = "\n".join(lines)
print(txt)
open(os.path.join(W, "results", "summary.txt"), "w").write(txt + "\n")
json.dump({"rows": rows, "summary": summ}, open(os.path.join(W, "results", "summary.json"), "w"), indent=1)
