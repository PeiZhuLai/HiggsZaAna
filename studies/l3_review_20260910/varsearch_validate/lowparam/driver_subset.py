#!/usr/bin/env python3
"""Low-mass BDT variables as a SUBSET of the production high-mass training variables (user rule 2026-10-14): the
high-mass model is unchanged; the low-mass model may only drop variables of the high-mass set (same definitions),
param = (m_gg - m_a)/m_llgg always kept; m_a = 1 signal trained on 115-133 GeV only (needed against the double peak).
No fixed number of variables.

--mode backward: start from the full high-mass set (16 incl. param); each round drops each remaining variable once
                 (one candidate per variable) and keeps the best candidate; runs the whole path down to [param].
--mode forward:  start from [param]; each round adds each remaining high-mass variable once; whole path up to 16.
Best candidate of a round: among total gate violation <= TOL the highest objective, else the smallest violation
(gates, objective, Z reference = stored production low-mass scores: the variable-search definitions, imported from
driver_lowparam.py = bdt_varsearch/driver_condor.py). The final choice is made afterwards from the recorded path.
Jobs: bdt_varsearch/run_candidate.sh + search_job.py, unchanged.
State (AFS): <state>/round_<k>/..., state.json, DONE; results (EOS): <results>/round_<k>/<id>.json
Usage: python3 driver_subset.py --mode backward|forward --state <afs dir> --results <eos dir> --samples <eos dir> --code <afs dir>
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import driver_lowparam as d  # noqa: E402

HIGH16 = ["pho1Pt_oHm", "pho1R9", "pho1IetaIeta55", "pho1PIso_noCorr", "pho2Pt_oHm", "pho2R9", "pho2IetaIeta55",
          "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "var_dR_Za", "var_dR_g1g2", "var_dR_g1Z", "var_PtaOverMh",
          "H_pt_oHm", "pho_pt_asym_oHm", "param"]
REFS = [{"id": "ref_low_stored", "model": "low", "feats": [], "stored_prod": True},
        {"id": "ref_low_prodvars", "model": "low", "feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]},
        {"id": "ref_low_prodvars_win", "model": "low", "win_mA1": True,
         "feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]},
        {"id": "ref_low_high16_win", "model": "low", "win_mA1": True, "feats": list(HIGH16)}]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=("backward", "forward"))
    ap.add_argument("--state", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--code", required=True)
    ap.add_argument("--results", required=True)
    a = ap.parse_args()
    os.makedirs(a.state, exist_ok=True)
    sp = os.path.join(a.state, "state.json")
    start = list(HIGH16) if a.mode == "backward" else ["param"]
    st = json.load(open(sp)) if os.path.exists(sp) else {"round": 0, "current": start, "history": []}
    while True:
        k = st["round"]
        cur = st["current"]
        if a.mode == "backward":
            if cur == ["param"]:
                break
            cands = [[f for f in cur if f != x] for x in cur if x != "param"]
            tag = "drop"
            moved = [x for x in cur if x != "param"]
        else:
            if len(cur) == len(HIGH16):
                break
            moved = [x for x in HIGH16 if x not in cur]
            cands = [cur + [x] for x in moved]
            tag = "add"
        cfgs = [dict(r) for r in REFS] if k == 0 else []
        for x, feats in zip(moved, cands):
            cfgs.append({"id": "low_%s_r%02d_%s_%s" % (a.mode, k, tag, x), "model": "low", "feats": feats, "win_mA1": True})
        rdir = os.path.join(a.state, "round_%02d" % k)
        cl = None
        if not os.path.exists(os.path.join(rdir, "config.jsonl")):
            cl = d.submit(rdir, cfgs, a)
            if cl:
                open(os.path.join(rdir, "cluster.txt"), "w").write(cl + "\n")
        elif os.path.exists(os.path.join(rdir, "cluster.txt")):
            cl = open(os.path.join(rdir, "cluster.txt")).read().strip()
        rres = d.wait_round(rdir, cfgs, a, cl)
        res = {}
        for c in cfgs:
            p = os.path.join(rres, c["id"] + ".json")
            if os.path.exists(p):
                res[c["id"]] = json.load(open(p))
        if k == 0:
            st["zref"] = {m: d.z_fit(r) for m, r in res["ref_low_stored"]["per_mass"].items()}
            st["refs"] = {i: {"violation": d.gate_violation(res[i]), "objective": d.objective(res[i], st["zref"])}
                          for i in res if i.startswith("ref_")}
            d.log("references: " + json.dumps(st["refs"]))
        cand = [(i, r) for i, r in res.items() if i.startswith("low_%s_r%02d_" % (a.mode, k))]
        if not cand:
            d.log("no results this round -> stop"); break
        scored = [(d.gate_violation(r), d.objective(r, st["zref"]), i, r) for i, r in cand]
        scored = [s for s in scored if d.finite(s[1])] or scored
        passing = [s for s in scored if s[0] <= d.TOL]
        best = max(passing, key=lambda s: s[1]) if passing else min(scored, key=lambda s: (s[0], -s[1]))
        v, obj, i, r = best
        st["current"] = r["config"]["feats"]
        st["history"].append({"round": k, "id": i, "feats": st["current"], "n": len(st["current"]), "violation": v,
                              "objective": obj, "n_passing": len(passing), "n_cand": len(scored)})
        d.log("%s round %d: %s -> %d vars, obj %.3f viol %.2f (%d/%d pass)" % (a.mode, k, i, len(st["current"]), obj, v, len(passing), len(scored)))
        st["round"] = k + 1
        json.dump(st, open(sp + ".partial", "w"), indent=1)
        os.replace(sp + ".partial", sp)
    d.log("FINISHED")
    open(os.path.join(a.state, "DONE"), "w").write(time.strftime("%F %T") + "\n")


if __name__ == "__main__":
    main()
