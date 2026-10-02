#!/usr/bin/env python3
"""
V2 study: plateau efficiencies of the trigger-efficiency curves, old vs new denominators.

Plateau = sum(pass_trigger)/sum(in_bin) over pT bins with low edge >= PLATEAU[leg]
(lead: 40 GeV, sublead: 30 GeV, both well above the single- and double-lepton offline
thresholds 35/25 (e) and 25/20 (mu) for single, 25/15 and 20/10 for double), plus the
open overflow bin. Uncertainty: Clopper-Pearson 68% (unweighted counts).

Usage: python plateau_table.py --dir <cutflow json dir> [--dir-old <prod dir>] [--eras ...]
"""
import argparse
import json
import os
import re

from scipy.stats import beta

PLATEAU = {"lead": 40, "sublead": 30}
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]


def cp(k, n, cl=0.682689492137086):
    a = (1 - cl) / 2
    lo = beta.ppf(a, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a, k + 1, n - k) if k < n else 1.0
    return lo, hi


def plateau(block, prefix, leg):
    b = (block.get(prefix) or {}).get("bins") or {}
    n = k = 0.0
    for name, rec in b.items():
        m = re.match(r"pt(\d+)to", name)
        if m and int(m.group(1)) >= PLATEAU[leg]:
            n += rec["in_bin"]
            k += rec["pass_trigger"]
    if n <= 0:
        return None
    lo, hi = cp(int(k), int(n))
    return k / n, lo, hi, int(n)


def fmt(r):
    if r is None:
        return "      --        "
    e, lo, hi, n = r
    return f"{100*e:6.2f}+{100*(hi-e):.2f}-{100*(e-lo):.2f} (n={n})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--dir-old", default=None)
    ap.add_argument("--eras", default=",".join(ERAS))
    ap.add_argument("--json-out", default=None)
    a = ap.parse_args()
    rows = []
    for era in a.eras.split(","):
        fn = os.path.join(a.dir, f"cutflow_Sig_MC_mA_M1_{era}.json")
        if not os.path.exists(fn):
            print(f"[missing] {fn}")
            continue
        blk = json.load(open(fn))["trigeff_nominal"]
        old = None
        if a.dir_old:
            fo = os.path.join(a.dir_old, f"cutflow_Sig_MC_mA_M1_{era}.json")
            if os.path.exists(fo):
                old = json.load(open(fo)).get("trigeff_nominal", {})
        print(f"\n=== {era}  (plateau: lead pT>={PLATEAU['lead']}, sublead pT>={PLATEAU['sublead']} GeV) ===")
        for lep, trigs in (("ele", ("double_ele", "OR_ele")), ("mu", ("double_mu", "OR_mu"))):
            for leg in ("lead", "sublead"):
                for trig in trigs:
                    r_prod = plateau(old, f"trigeff_{lep}_{leg}_{trig}", leg) if old is not None else None
                    r_old = plateau(blk, f"trigeff_{lep}_{leg}_{trig}", leg)
                    r_n2 = plateau(blk, f"trigeff_{lep}_{leg}N2_{trig}", leg)
                    r_ol = plateau(blk, f"trigeff_{lep}_{leg}OL_{trig}", leg)
                    print(f"{lep:3s} {leg:7s} {trig:10s} prodJSON {fmt(r_prod)} | old-denom {fmt(r_old)} | N2 {fmt(r_n2)} | OL {fmt(r_ol)}")
                    rows.append(dict(era=era, lep=lep, leg=leg, trig=trig,
                                     prod=r_prod, old=r_old, n2=r_n2, ol=r_ol))
    if a.json_out:
        json.dump(rows, open(a.json_out, "w"), indent=1)


if __name__ == "__main__":
    main()
