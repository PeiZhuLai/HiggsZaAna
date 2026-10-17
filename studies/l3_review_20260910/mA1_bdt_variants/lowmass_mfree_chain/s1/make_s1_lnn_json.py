#!/usr/bin/env python3
"""S1 scheme 2 -> per-mA asymmetric lnN for CMS_hza_mva_reweight.

Two references (--reference):
  none     (user decision 2026-09-27): the datacard signal is NOT reweighted; envelope of
           {nominal h_m, z_m&hsb, z_m_narrow&hsb} relative to no reweight.
  nominal  (user decision 2026-10-02, option B): the signal is reweighted with the nominal h_m
           sideband reweight (preselected yield preserved per channel); envelope of
           {z_m&hsb, z_m_narrow&hsb, none} relative to the nominal reweight.
    d_k = eff_k / eff_ref - 1,  kappa_up = 1 + max(0, d_k),  kappa_down = 1 + min(0, d_k)
eff_k = BDT-WP efficiency with reweight k normalized to the preselected yield per channel
(eval_alt_reweight.py). One value per mA, channel 'all', eras combined with the nominal-
reweighted yield as weight -> correlated across eras and channels.

Usage: make_s1_lnn_json.py [--reference none|nominal] <out.json> <yields.csv> [<yields.csv> ...]
       later CSVs override the mass points they contain (e.g. a re-evaluation at a new WP).
"""
import argparse, json, time
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--reference", choices=("none", "nominal"), default="none")
ap.add_argument("out")
ap.add_argument("csvs", nargs="+")
args = ap.parse_args()
REF = args.reference
VARS = ("nominal", "zmhsb", "zmnarrowhsb") if REF == "none" else ("zmhsb", "zmnarrowhsb", "none")

d = None
for src in args.csvs:
    x = pd.read_csv(src)
    x = x[x.ch == "all"].copy()
    d = x if d is None else pd.concat([d[~d.mA.isin(set(x.mA))], x], ignore_index=True)
missing = [f"eff_{k}" for k in set(VARS) | {REF, "none", "nominal"} if f"eff_{k}" not in d.columns]
if missing:
    raise SystemExit(f"CSV lacks {missing}")
d["yield_nom"] = d["yield_none"] * d["eff_nominal"] / d["eff_none"]
vals, table = {}, []
for m, g in d.groupby("mA"):
    w = g["yield_nom"].to_numpy(float); w = w / w.sum()
    dk = {k: float(np.sum(w * (g[f"eff_{k}"] / g[f"eff_{REF}"] - 1.0))) for k in VARS}
    up = 1.0 + max(0.0, *dk.values()); dn = 1.0 + min(0.0, *dk.values())
    cuts = sorted(set(round(float(c), 6) for c in g["wp"]))
    if len(cuts) != 1:
        raise SystemExit(f"mA{m}: rows were evaluated at different working points {cuts}")
    # mvacut = working point the envelope was evaluated at; apply_bdt_sig.py refuses a mismatch
    vals[str(int(m))] = {"up": round(up, 5), "down": round(dn, 5), "mvacut": cuts[0]}
    table.append((int(m), dk, up, dn))
json.dump({"scheme": "S1 scheme 2: envelope of {%s} vs %s reweight" % (", ".join(VARS), REF),
           "reference": REF, "source": args.csvs, "created": time.strftime("%F %T"), "values": vals},
          open(args.out, "w"), indent=1)
print("mA | " + " | ".join(VARS) + f"  (vs {REF}, %) | kappa_down / kappa_up")
for m, dk, up, dn in table:
    print(f"{m:>2} | " + " | ".join(f"{100*dk[k]:+.1f}" for k in VARS) + f" | {dn:.3f} / {up:.3f}")
print("wrote", args.out)
