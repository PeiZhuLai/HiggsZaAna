#!/usr/bin/env python3
"""S1 scheme 2 (user decision 2026-09-27) -> per-mA asymmetric lnN for CMS_hza_mva_reweight.

The datacard central value is the un-reweighted signal yield. The nuisance is the envelope of
the signal efficiency after the BDT WP under the {nominal h_m, z_m&hsb, z_m_narrow&hsb}
reweights, each relative to NO reweight:
    d_k = eff_k / eff_none - 1,  k in {nominal, zmhsb, zmnarrowhsb}
    kappa_up = 1 + max(0, d_k),  kappa_down = 1 + min(0, d_k)
One value per mA, channel 'all' (ele+mu together), combined over the five eras with the same
yield weighting as summarize.py -> correlated across eras and channels.

Usage: make_s1_lnn_json.py <yields.csv from eval_alt_reweight.py> <out.json>
"""
import json, sys, time
import numpy as np, pandas as pd

VARS = ("nominal", "zmhsb", "zmnarrowhsb")
src, out = sys.argv[1], sys.argv[2]
d = pd.read_csv(src)
d = d[d.ch == "all"].copy()
missing = [f"eff_{k}" for k in VARS + ("none",) if f"eff_{k}" not in d.columns]
if missing:
    sys.exit(f"{src} lacks {missing}")
d["yield_nom"] = d["yield_none"] * d["eff_nominal"] / d["eff_none"]
vals, table = {}, []
for m, g in d.groupby("mA"):
    w = g["yield_nom"].to_numpy(float); w = w / w.sum()
    dk = {k: float(np.sum(w * (g[f"eff_{k}"] / g["eff_none"] - 1.0))) for k in VARS}
    up = 1.0 + max(0.0, *dk.values()); dn = 1.0 + min(0.0, *dk.values())
    cuts = sorted(set(round(float(c), 6) for c in g["wp"]))
    if len(cuts) != 1:
        sys.exit(f"mA{m}: rows were evaluated at different working points {cuts}")
    # mvacut = working point the envelope was evaluated at; apply_bdt_sig.py refuses a mismatch
    vals[str(int(m))] = {"up": round(up, 5), "down": round(dn, 5), "mvacut": cuts[0]}
    table.append((int(m), dk, up, dn))
json.dump({"scheme": "S1 scheme 2: envelope of {nominal, z_m&hsb, z_m_narrow&hsb} vs no reweight",
           "source": src, "created": time.strftime("%F %T"), "values": vals}, open(out, "w"), indent=1)
print("mA | nominal | zm&hsb | zmnarrow&hsb  (vs none, %) | kappa_down / kappa_up")
for m, dk, up, dn in table:
    print(f"{m:>2} | {100*dk['nominal']:+.1f} | {100*dk['zmhsb']:+.1f} | {100*dk['zmnarrowhsb']:+.1f} | {dn:.3f} / {up:.3f}")
print("wrote", out)
