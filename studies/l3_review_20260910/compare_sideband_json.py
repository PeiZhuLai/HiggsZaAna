#!/usr/bin/env python
"""
Compare two sideband reweight JSONs (pre-FSR-fix vs FSR-fix).

Why this matters: Parque2Root_BDT.py writes sideband-reweight branches into every ROOT
file using this JSON, but the JSON is itself derived from Data/run3.root and
All_Bkg/run3.root -- i.e. from p2root output. The FSR-fix ROOT files were converted
using the OLD JSON (derived 2026-06-03 from the pre-FSR-fix production).

That circularity is only a problem if the JSON actually moved. If the reweight factors
are stable to well within their own statistical noise, the ROOT files stand. If they
moved, p2root has to be rerun with the regenerated JSON.

Usage: python compare_sideband_json.py OLD.json NEW.json
"""
import json, sys
import numpy as np


def main(old_p, new_p):
    o = json.load(open(old_p))
    n = json.load(open(new_p))
    print("old: %s   (created %s)" % (old_p, o.get("created")))
    print("new: %s   (created %s)" % (new_p, n.get("created")))
    print()

    on, nn = o.get("normalization", {}), n.get("normalization", {})
    print("--- normalization ---")
    for k in ("data_sideband_yield", "bkg_initial_sideband_yield",
              "bkg_final_sideband_yield", "data_over_bkg_initial"):
        a, b = on.get(k), nn.get(k)
        if isinstance(a, (int, float)) and isinstance(b, (int, float)) and a:
            print("  %-28s %14.4f -> %14.4f   %+.3f%%" % (k, a, b, 100.0 * (b - a) / a))
        else:
            print("  %-28s %s -> %s" % (k, a, b))

    print()
    print("--- per-variable reweight factors, all iterations pooled ---")
    print("  %-26s %6s %10s %10s %10s" % ("variable", "bins", "max|d|%", "rms|d|%", "mean|d|%"))

    # The JSON stores iterations[i]["steps"] = [{var, edges, clipped_factors, ...}].
    # An earlier version of this script looked for iterations[i][var]["weights"], found
    # nothing, and printed "worst change 0.000%" -- a PASS produced by an empty input.
    # Hence the explicit emptiness guard at the end.
    def collect(doc):
        out = {}
        for it in doc.get("iterations", []):
            for st in it.get("steps", []):
                v = st.get("var")
                f = st.get("clipped_factors")
                if v is None or not f:
                    continue
                out.setdefault(v, []).append(np.asarray(f, dtype=float))
        return out

    co, cn = collect(o), collect(n)
    shared = sorted(set(co) & set(cn))
    worst, nvar, ncmp = 0.0, 0, 0
    for var in shared:
        ds = []
        for a, b in zip(co[var], cn[var]):
            if a.size != b.size:
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                d = np.abs(np.where(a != 0, (b - a) / a, 0.0)) * 100.0
            d = d[np.isfinite(d)]
            if d.size:
                ds.append(d)
        if not ds:
            continue
        d = np.concatenate(ds)
        worst = max(worst, float(d.max()))
        nvar += 1
        ncmp += d.size
        print("  %-26s %6d %10.3f %10.3f %10.3f"
              % (var, d.size, d.max(), np.sqrt((d ** 2).mean()), d.mean()))

    print()
    only_o, only_n = sorted(set(co) - set(cn)), sorted(set(cn) - set(co))
    if only_o or only_n:
        print("variables only in old: %s" % only_o)
        print("variables only in new: %s" % only_n)
    print("variables compared: %d   bin comparisons: %d" % (nvar, ncmp))
    if nvar == 0 or ncmp == 0:
        print("ERROR: nothing was actually compared -- this is not a PASS.")
        return 2
    print("worst single-bin reweight change: %.3f%%" % worst)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
