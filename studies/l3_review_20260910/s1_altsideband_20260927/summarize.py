#!/usr/bin/env python3
"""Summarize s1_altsideband_yields.csv: Run-3 combination per mass and the proposed envelope.

Run-3 value per mass = nominal-yield-weighted mean over the five eras (channel 'all'); the
reweighting is derived once on the Run-3 sample, so one value per mass correlated across eras is
the natural granularity. Stat error = same weighting of the per-era bootstrap errors (in quadrature).
Proposed envelope (sidebands that exclude 115<m_llgg<135): up = max(0, dy_zmhsb, dy_zmnarrowhsb),
down = min(0, dy_zmhsb, dy_zmnarrowhsb).
"""
import json, sys, math
import numpy as np, pandas as pd

S = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/s1_altsideband_20260927"
d = pd.read_csv(f"{S}/results/s1_altsideband_yields.csv")
alts = ["zm", "zmnarrow", "zmhsb", "zmnarrowhsb"]
d["yield_nom"] = d["yield_none"] * d["eff_nominal"] / d["eff_none"]
for k in ["none", "nominal"] + alts:
    pass
out = []
for (m, ch), g in d.groupby(["mA", "ch"]):
    w = g["yield_nom"].to_numpy()
    w = w / w.sum()
    r = {"mA": int(m), "ch": ch}
    for k in alts + ["nominal"]:
        key = f"dy_{k}" if k != "nominal" else "dy_nominal_vs_none"
        ekey = f"dyerr_{k}" if k != "nominal" else "dyerr_nominal_vs_none"
        r[key] = float(np.sum(w * g[key]))
        r[ekey] = float(math.sqrt(np.sum((w * g[ekey]) ** 2))) if False else float(np.sqrt(np.sum((w * g[ekey]) ** 2)))
        r[key + "_eramin"] = float(g[key].min()); r[key + "_eramax"] = float(g[key].max())
    r["cur_up"] = float(np.sum(w * g["cur_up"])); r["cur_dn"] = float(np.sum(w * g["cur_dn"]))
    for k in ["nominal"] + alts:
        r[f"dmean_{k}"] = float(np.sum(w * (g[f"mean_{k}"] - g["mean_none" if k == "nominal" else "mean_nominal"])))
        r[f"rseff_{k}"] = float(np.sum(w * (g[f"seff_{k}"] / g["seff_none" if k == "nominal" else "seff_nominal"] - 1)))
        r[f"rrms_{k}"] = float(np.sum(w * (g[f"rms_{k}"] / g["rms_none" if k == "nominal" else "rms_nominal"] - 1)))
    hs = [r["dy_zmhsb"], r["dy_zmnarrowhsb"]]
    r["env_up"] = max(0.0, *hs); r["env_dn"] = min(0.0, *hs)
    r["env_sym"] = max(abs(x) for x in hs)
    r["env_all_sym"] = max(abs(r[f"dy_{k}"]) for k in alts)
    out.append(r)
o = pd.DataFrame(out).sort_values(["ch", "mA"])
o.to_csv(f"{S}/results/s1_altsideband_run3_summary.csv", index=False)

a = o[o.ch == "all"]
pct = lambda x: "%+.1f" % (100 * x)
print("Run-3 (lumi/yield-weighted over 5 eras), channel ele+mu, relative yield after BDT WP [%]")
print("mA | zm | zmnarrow | zm&hsb | zmnarrow&hsb | stat(1sig) | nom-vs-none | current(up/dn) | proposed(dn/up)")
for _, r in a.iterrows():
    print(f"{r.mA:>2} | {pct(r.dy_zm)} | {pct(r.dy_zmnarrow)} | {pct(r.dy_zmhsb)} | {pct(r.dy_zmnarrowhsb)} | "
          f"{100*max(r.dyerr_zmhsb, r.dyerr_zmnarrowhsb):.1f} | {pct(r.dy_nominal_vs_none)} | "
          f"{pct(r.cur_up)}/{pct(r.cur_dn)} | {pct(r.env_dn)}/{pct(r.env_up)}")
print()
print("Shape (Run-3, ele+mu): d(mean) [GeV] and relative sigma_eff change [%] vs nominal")
for _, r in a.iterrows():
    print(f"{r.mA:>2} | nom-vs-none dmean={r.dmean_nominal:+.3f} dseff={100*r.rseff_nominal:+.2f}% | " + " | ".join(
        f"{k}: {r['dmean_'+k]:+.3f}/{100*r['rseff_'+k]:+.2f}%" for k in alts))
print()
print("era spread (min..max) of dy_zmnarrowhsb and dy_zmhsb per mass [%]:")
for _, r in a.iterrows():
    print(f"{r.mA:>2}: zm&hsb {pct(r.dy_zmhsb_eramin)}..{pct(r.dy_zmhsb_eramax)}   zmnarrow&hsb {pct(r.dy_zmnarrowhsb_eramin)}..{pct(r.dy_zmnarrowhsb_eramax)}")
print()
print("ele vs mu (Run-3) proposed envelope dn/up [%]:")
for m in sorted(o.mA.unique()):
    e = o[(o.mA == m) & (o.ch == "ele")].iloc[0]; u = o[(o.mA == m) & (o.ch == "mu")].iloc[0]
    print(f"{m:>2}: ele {pct(e.env_dn)}/{pct(e.env_up)}  mu {pct(u.env_dn)}/{pct(u.env_up)}")
