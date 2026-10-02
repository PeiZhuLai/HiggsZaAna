#!/usr/bin/env python3
"""L3 review S1: which reweight steps drive the alt-vs-nominal signal-efficiency difference.

For each variant, the per-event weight factorizes into per-variable pieces
R = prod_v R_v (R_v = product over iterations of that variable's step factor x norm).
Reports, on the signal test tree after the BDT WP:
  swap_v : eff(R_nom * R_alt,v / R_nom,v) / eff(R_nom) - 1   (only variable v taken from alt)
  dy_excl: alt vs nominal with the control-region-defining variables {Z_m, H_m} and param
           removed from BOTH (object-level steps only)
Also the mean per-variable factor on the signal (R_v mean) to see which steps act on the signal.
"""
import sys, json, argparse
import numpy as np, pandas as pd, uproot
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/s1_altsideband_20260927")
from eval_alt_reweight import VARIANTS, SCORED, READ, wp_cuts, FIT_WINDOW, DEFAULT_WINDOW  # noqa
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from sideband_reweight import SidebandReweighter, _lookup, _finite_float  # noqa


def per_var(rw, frame):
    out = {}
    for it in rw.iterations:
        for st in it["steps"]:
            v = st["var"]
            f = _lookup(rw._dataframe_values(frame, v), st["edges"], st["clipped_factors"]) * _finite_float(st.get("normalization_scale", 1.0), 1.0)
            out[v] = out.get(v, np.ones(len(frame))) * f
    return out


ap = argparse.ArgumentParser()
ap.add_argument("--masses", default="1,3,5,10,20,30")
ap.add_argument("--eras", default="2024")
ap.add_argument("--out", required=True)
a_ = ap.parse_args()
cuts = wp_cuts()
rws = {k: SidebandReweighter.from_json(v) for k, v in VARIANTS.items()}
rows = []
for m in map(int, a_.masses.split(",")):
    lo, hi = FIT_WINDOW.get(m, DEFAULT_WINDOW)
    for era in a_.eras.split(","):
        a = uproot.open(f"{SCORED}/mA_M{m}/{era}.root")["test"].arrays(READ + [f"MVA_Score_mA_M{m}"], library="pd")
        a["param"] = (a["ALP_m"] - m) / a["H_m"]
        w = a["weight"].to_numpy(float)
        sel = (a[f"MVA_Score_mA_M{m}"].to_numpy() > cuts[m]) & (a["H_mass"].between(lo, hi).to_numpy())
        ch = ((a["n_electrons"] == 2) | (a["n_muons"] == 2)).to_numpy()
        sel &= ch
        pv = {k: per_var(rw, a) for k, rw in rws.items()}

        def eff(r):
            return np.sum((w * r)[sel]) / np.sum((w * r)[ch])
        Rn = np.prod(list(pv["nominal"].values()), axis=0)
        en = eff(Rn)
        vars_ = list(pv["nominal"].keys())
        for k in rws:
            row = {"mA": m, "era": era, "variant": k}
            Rk = np.prod(list(pv[k].values()), axis=0)
            row["dy_total"] = eff(Rk) / en - 1
            drop = {"Z_m", "H_m", "param"}
            Rn_x = np.prod([pv["nominal"][v] for v in vars_ if v not in drop], axis=0)
            Rk_x = np.prod([pv[k][v] for v in vars_ if v not in drop], axis=0)
            row["dy_excl_Zm_Hm_param"] = eff(Rk_x) / eff(Rn_x) - 1
            Rn_z = np.prod([pv["nominal"][v] for v in vars_ if v not in {"Z_m", "H_m"}], axis=0)
            Rk_z = np.prod([pv[k][v] for v in vars_ if v not in {"Z_m", "H_m"}], axis=0)
            row["dy_excl_Zm_Hm"] = eff(Rk_z) / eff(Rn_z) - 1
            row["nom_vs_none_excl_Zm_Hm_param"] = eff(Rn_x) / eff(np.ones_like(Rn_x)) - 1
            for v in vars_:
                r = Rn * pv[k][v] / pv["nominal"][v]
                row[f"swap_{v}"] = eff(r) / en - 1
                row[f"meanR_{v}"] = float(np.sum(w[ch] * pv[k][v][ch]) / np.sum(w[ch]))
            rows.append(row)
        print(m, era, "done", flush=True)
d = pd.DataFrame(rows)
d.to_csv(a_.out, index=False)
pd.set_option("display.width", 250)
cols = ["mA", "era", "variant", "dy_total", "dy_excl_Zm_Hm", "dy_excl_Zm_Hm_param", "nom_vs_none_excl_Zm_Hm_param"]
print(d[cols].to_string(float_format=lambda x: "%+.4f" % x))
sw = [c for c in d.columns if c.startswith("swap_")]
print(d[["mA", "variant"] + sw].to_string(float_format=lambda x: "%+.3f" % x))
mr = [c for c in d.columns if c.startswith("meanR_")]
print(d[["mA", "variant"] + mr].to_string(float_format=lambda x: "%.3f" % x))
