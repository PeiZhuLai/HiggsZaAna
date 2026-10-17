#!/usr/bin/env python3
"""Train and evaluate a variant of a selected model with some features removed (same recipe, events and float64
inputs as train_selected.py; same metrics, gates and DY versions as eval_selected.py). selected.json is not changed:
the variant is added to the in-memory configuration of both modules.
Env: higgs-alp-ana. Output: models/model_<name>.{pkl,meta.json}, eval/<name>_new_<incl|test>.json,
eval/summary_<name>.txt
Usage: train_eval_variant.py <base model low|high> <name> <feature to drop> [...]
"""
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import train_selected as ts  # noqa: E402
import eval_selected as es  # noqa: E402
dc = es.dc


def main():
    base, name, drop = sys.argv[1], sys.argv[2], sys.argv[3:]
    cfg = dict(ts.SEL[base])
    missing = [f for f in drop if f not in cfg["feats"]]
    if missing:
        sys.exit("not in the %s model: %s" % (base, missing))
    cfg["feats"] = [f for f in cfg["feats"] if f not in drop]
    ts.SEL[name] = es.SEL[name] = cfg
    pkl = os.path.join(HERE, "models", "model_%s.pkl" % name)
    if not os.path.exists(pkl):
        sys.argv = [sys.argv[0], name]
        ts.main()
    clf = pickle.load(open(pkl, "rb"))
    S = es.load_cache("/eos/home-p/pelai/HZa/sandbox_varsearch_v1/eval_samples64")
    n_target = json.load(open(os.path.join(es.SLIM, "meta.json")))["prod_data_sb_count"]
    out = {}
    for variant in ("incl", "test"):
        out[variant] = es.run(S, name, clf, variant, n_target)
        json.dump(out[variant], open(os.path.join(HERE, "eval", "%s_new_%s.json" % (name, variant)), "w"), indent=1)
    ref = {(w, v): json.load(open(os.path.join(HERE, "eval", "%s_%s_%s.json" % (base, w, v))))
           for w in ("new", "prod") for v in ("incl", "test")}
    zref = {m: dc.z_fit(r) for m, r in ref[("prod", "incl")]["per_mass"].items()}
    lines = ["%s = %s minus %s: %s" % (name, base, ", ".join(drop), ", ".join(cfg["feats"]))]
    for tag, r in (("%s incl" % name, out["incl"]), ("%s test" % name, out["test"]), ("%s new incl" % base, ref[("new", "incl")]),
                   ("%s new test" % base, ref[("new", "test")]), ("prod incl", ref[("prod", "incl")]), ("prod test", ref[("prod", "test")])):
        lines.append("%-17s violation %6.2f  objective %.3f" % (tag, dc.gate_violation(r), dc.objective(r, zref)))
    lines.append(" mA | Z/Zprod variant base | R test variant   base        prod | A pulls test (variant)  | B p    | win/ovl (presel)")
    for m in out["test"]["per_mass"]:
        t, b, p = out["test"]["per_mass"][m], ref[("new", "test")]["per_mass"][m], ref[("prod", "test")]["per_mass"][m]
        lines.append("%3s | %7.2f %5.2f | %5.2f+-%.2f %5.2f+-%.2f %5.2f | %5.1f %5.1f %5.1f | %.4f | %.2f/%.2f (%.2f/%.2f)"
                     % (m, dc.z_fit(out["incl"]["per_mass"][m]) / zref[m], dc.z_fit(ref[("new", "incl")]["per_mass"][m]) / zref[m],
                        t["R"], t["R_err"], b["R"], b["R_err"], p["R"], *t["A_pulls"], t["B_p"], t["win_frac"], t["overlap"],
                        t["presel_win_frac"], t["presel_overlap"]))
    txt = "\n".join(lines)
    open(os.path.join(HERE, "eval", "summary_%s.txt" % name), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
