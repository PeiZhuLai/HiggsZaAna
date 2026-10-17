#!/usr/bin/env python3
"""Forward-selection driver for the HZa BDT input-variable search on IHEP (hep_sub arrays).

Runs on an IHEP login node (hzgdna env), one process for both models. Each round, for each model that
is still growing, every pool variable not yet selected is added to the current set as one candidate;
the candidates of a round go out as ONE hep_sub array (search_job.py). When all results are in, the
best candidate is kept:
  gates (every anchor mass of the model): |R-1|, |R_lo-1|, |R_hi-1| <= 2 sigma (sigma floor 0.05),
          data-sideband shape p >= 0.01; m_a = 1 also: overlap share <= preselection and window fraction
          >= preselection (no BDT-enriched second population)
  among candidates passing every gate: highest objective = mean over masses of Z / Z_production
          (Z_production from the stored production scores at the same data-sideband count)
  if none passes: smallest total gate violation (sum of excess sigmas / log10 p deficits)
A model stops when the best objective improves by < 1% for two rounds, or at 12 variables.
Round 0 also runs the references: stored production scores, production variables retrained, and the
sandbox win_mfree (low) / shared-definition (high) models.

State: <state>/round_<k>/{config.jsonl, results/, submit.log}, <state>/state.json, heartbeat.txt, driver.log
Restartable: completed rounds are read back from state.json; a round whose results all exist is not resubmitted.
Usage (IHEP): nohup python driver.py --state <dir> --samples <dir> --code <dir> > driver.log 2>&1 &
"""
import argparse
import json
import math
import os
import re
import subprocess
import time

HEP = "/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin"
POOL = ["pho1Pt", "pho2Pt", "pho_pt_asym", "pho1R9", "pho2R9", "pho1IetaIeta55", "pho2IetaIeta55",
        "pho1PIso_noCorr", "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "ALP_PhotonIso", "ALP_pt", "ALP_eta",
        "var_PtaOverMa", "var_dR_g1g2", "var_dR_Za", "var_dR_g1Z", "var_dEta_g1Z", "var_dPhi_g1Z",
        "min_lg_deltaR", "max_lg_deltaR", "H_pt", "pho1_relEerr", "pho2_relEerr"] + \
       ["ALP_%s_photon_%s" % (w, v) for w in ("lead", "sublead") for v in
        ("mvaID", "eta", "hoe", "s4", "sieip", "etaWidth", "phiWidth", "esEffSigmaRR", "chiso",
         "ecalPFClusterIso", "hcalPFClusterIso", "pfChargedIsoWorstVtx", "trkSumPtHollowConeDR03", "lep_near_dR")]
REFS = [
    {"id": "ref_low_stored", "model": "low", "feats": [], "stored_prod": True},
    {"id": "ref_high_stored", "model": "high", "feats": [], "stored_prod": True},
    {"id": "ref_low_prodvars", "model": "low", "feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]},
    {"id": "ref_low_winmfree", "model": "low", "win_mA1": True, "feats": ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt"]},
    {"id": "ref_high_prodvars", "model": "high", "feats": ["pho1Pt_oHm", "pho1R9", "pho1IetaIeta55", "pho1PIso_noCorr",
        "pho2Pt_oHm", "pho2R9", "pho2IetaIeta55", "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "var_dR_Za", "var_dR_g1g2",
        "var_dR_g1Z", "var_PtaOverMh", "H_pt_oHm", "pho_pt_asym_oHm", "param"]},
    {"id": "ref_high_shared", "model": "high", "feats": ["pho1Pt", "pho1R9", "pho1IetaIeta55", "pho1PIso_noCorr",
        "pho2Pt_oHm", "pho2R9", "pho2IetaIeta55", "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "var_dR_Za", "var_dR_g1g2",
        "var_dR_g1Z", "var_PtaOverMh", "H_pt_oHm", "pho_pt_asym_oHm"]},
]
MAXVARS = 12


def log(msg):
    print(time.strftime("%F %T"), msg, flush=True)


def gate_violation(res):
    v = 0.0
    for m, r in res["per_mass"].items():
        for k in ("R", "R_lo", "R_hi"):
            x, e = r[k], max(r[k + "_err"], 0.05)
            if not (x == x):
                v += 10; continue
            v += max(0.0, abs(x - 1) / e - 2)
        p = r["data_p"]
        v += 10 if not (p == p) else max(0.0, math.log10(0.01 / max(p, 1e-300)))
        if m == "1":
            v += max(0.0, (r["overlap"] - r["presel_overlap"]) / 0.02)
            v += max(0.0, (r["presel_win_frac"] - r["win_frac"]) / 0.02)
    return v


def objective(res, zref):
    zs = [res["per_mass"][m]["Z"] / zref[m] for m in res["per_mass"] if zref.get(m)]
    return sum(zs) / len(zs) if zs else float("nan")


def submit(rdir, cfgs, a):
    """Write config + job script for one array and submit it; returns the cluster id (or None)."""
    os.makedirs(os.path.join(rdir, "results"), exist_ok=True)
    os.makedirs(os.path.join(rdir, "logs"), exist_ok=True)
    with open(os.path.join(rdir, "config.jsonl"), "w") as f:
        for c in cfgs:
            f.write(json.dumps(c) + "\n")
    job = os.path.join(rdir, "job.sh")
    with open(job, "w") as f:
        f.write("#!/bin/bash\n# one search candidate per ProcId; per-job log next to the results\n"
                "exec > %s/logs/job_$1.log 2>&1\nexport OMP_NUM_THREADS=1 NTHREADS=1\n"
                "%s/bin/python %s/search_job.py --samples %s --config %s/config.jsonl --index $1 --out %s/results\n"
                % (rdir, a.env, a.code, a.samples, rdir, rdir))
    os.chmod(job, 0o755)
    cmd = [HEP + "/hep_sub", job, "-g", "cms", "-m", "3500", "-n", str(len(cfgs)), "-argu", "%{ProcId}",
           "-o", os.path.join(rdir, "logs", "condor.out"), "-e", os.path.join(rdir, "logs", "condor.err")]
    log("submit: " + " ".join(cmd))
    out = subprocess.run(cmd, capture_output=True, text=True)
    open(os.path.join(rdir, "submit.log"), "a").write(out.stdout + out.stderr)
    m = re.search(r"submitted to cluster (\d+)", out.stdout)
    log("cluster %s (%d jobs)" % (m.group(1) if m else "NONE", len(cfgs)))
    return m.group(1) if m else None


def queue_status(clusters):
    """Status letters (hep_q ST column) of the jobs of the given clusters; None if hep_q failed."""
    out = subprocess.run([HEP + "/hep_q", "-u", "laipeizhu"], capture_output=True, text=True)
    if out.returncode != 0 or "jobs;" not in out.stdout:
        return None
    st = {}
    for line in out.stdout.splitlines():
        f = line.split()
        if len(f) >= 6 and "." in f[0] and f[0].split(".")[0] in clusters:
            st[f[0]] = f[5]
    return st


def handle_held(st):
    """Held-job SOP: raise the memory request 8000 -> 16000 -> 20000 MB and release; stop at 20000."""
    held = [j for j, s in st.items() if s == "H"]
    for j in held:
        n = HELD_TRIES.get(j, 0)
        if n >= 3:
            continue
        mem = (8000, 16000, 20000)[n]
        HELD_TRIES[j] = n + 1
        e = subprocess.run([HEP + "/hep_edit", "-m", str(mem), j], capture_output=True, text=True)
        r = subprocess.run([HEP + "/hep_release", j], capture_output=True, text=True)
        log("held %s -> hep_edit -m %d (%s) + hep_release (%s)" % (j, mem, e.stdout.strip()[-80:], r.stdout.strip()[-80:]))


def wait_round(rdir, cfgs, a, cluster):
    want = {c["id"] for c in cfgs}
    clusters = {cluster} if cluster else set()
    tries = 0
    t_last_submit = time.time()   # hep_q may not list a fresh array yet: no "queue empty" verdict for 10 min
    while True:
        have = {f[:-5] for f in os.listdir(os.path.join(rdir, "results")) if f.endswith(".json")}
        st = queue_status(clusters)
        n_run = sum(1 for v in (st or {}).values() if v == "R")
        n_idle = sum(1 for v in (st or {}).values() if v == "I")
        n_held = sum(1 for v in (st or {}).values() if v == "H")
        open(os.path.join(a.state, "heartbeat.txt"), "w").write(
            "%s %s results %d/%d  queue: %d running %d idle %d held\n"
            % (time.strftime("%F %T"), os.path.basename(rdir), len(have & want), len(want), n_run, n_idle, n_held))
        if want <= have:
            return
        if st is None:
            time.sleep(120); continue
        if n_held:
            handle_held(st)
        if not st and time.time() - t_last_submit > 600:   # nothing left in the queue but results missing
            missing = [c for c in cfgs if c["id"] not in have]
            if tries >= 3:
                log("giving up on %d candidates after %d resubmissions: %s" % (len(missing), tries, [c["id"] for c in missing][:5]))
                return
            tries += 1
            log("round %s: queue empty, %d results missing -> resubmission %d" % (os.path.basename(rdir), len(missing), tries))
            rr = rdir + "_retry%d" % tries
            cl = submit(rr, missing, a)
            t_last_submit = time.time()
            if cl:
                clusters.add(cl)
            # retry jobs write into rr/results; move finished ones back as they appear
            RETRY_DIRS.append(rr)
        for rr in RETRY_DIRS:
            if rr.startswith(rdir + "_retry"):
                for f in os.listdir(os.path.join(rr, "results")):
                    dst = os.path.join(rdir, "results", f)
                    if f.endswith(".json") and not os.path.exists(dst):
                        os.replace(os.path.join(rr, "results", f), dst)
        time.sleep(60)


HELD_TRIES = {}
RETRY_DIRS = []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--code", required=True)
    ap.add_argument("--env", default="/publicfs/cms/user/laipeizhu/App/Conda/envs/hzgdna")
    a = ap.parse_args()
    os.makedirs(a.state, exist_ok=True)
    sp = os.path.join(a.state, "state.json")
    st = json.load(open(sp)) if os.path.exists(sp) else {"round": 0, "models": {
        "low": {"current": [], "history": [], "stall": 0, "done": False, "win_mA1": True},
        "high": {"current": [], "history": [], "stall": 0, "done": False, "win_mA1": False}}}
    while True:
        k = st["round"]
        rdir = os.path.join(a.state, "round_%02d" % k)
        cfgs = list(REFS) if k == 0 else []
        for model, ms in st["models"].items():
            if ms["done"]:
                continue
            for f in POOL:
                if f in ms["current"]:
                    continue
                feats = ms["current"] + [f]
                cfgs.append({"id": "%s_r%02d_%s" % (model, k, f), "model": model, "feats": feats, "win_mA1": ms["win_mA1"]})
        if not cfgs:
            log("all models done"); break
        cl = None
        if not os.path.exists(os.path.join(rdir, "config.jsonl")):
            cl = submit(rdir, cfgs, a)
            if cl:
                open(os.path.join(rdir, "cluster.txt"), "w").write(cl + "\n")
        elif os.path.exists(os.path.join(rdir, "cluster.txt")):
            cl = open(os.path.join(rdir, "cluster.txt")).read().strip()
        wait_round(rdir, cfgs, a, cl)
        res = {}
        for c in cfgs:
            p = os.path.join(rdir, "results", c["id"] + ".json")
            if os.path.exists(p):
                res[c["id"]] = json.load(open(p))
        if k == 0:
            st["zref"] = {model: {m: r["Z"] for m, r in res["ref_%s_stored" % model]["per_mass"].items()}
                          for model in ("low", "high") if "ref_%s_stored" % model in res}
            st["refs"] = {i: {"violation": gate_violation(res[i]), "objective": objective(res[i], st["zref"]["low" if "_low_" in i else "high"])}
                          for i in res if i.startswith("ref_")}
            log("references: " + json.dumps(st["refs"]))
        for model, ms in st["models"].items():
            if ms["done"]:
                continue
            cand = [(i, r) for i, r in res.items() if i.startswith(model + "_r")]
            if not cand:
                log("%s: no results this round -> stop" % model); ms["done"] = True; continue
            scored = [(gate_violation(r), objective(r, st["zref"][model]), i, r) for i, r in cand]
            passing = [s for s in scored if s[0] == 0]
            best = max(passing, key=lambda s: s[1]) if passing else min(scored, key=lambda s: (s[0], -s[1]))
            v, obj, i, r = best
            prev = ms["history"][-1]["objective"] if ms["history"] else None
            ms["current"] = r["config"]["feats"]
            ms["history"].append({"round": k, "added": ms["current"][-1], "feats": ms["current"], "violation": v,
                                  "objective": obj, "n_passing": len(passing), "id": i})
            log("%s round %d: + %s -> obj %.3f viol %.2f (%d/%d candidates pass)" % (model, k, ms["current"][-1], obj, v, len(passing), len(scored)))
            improved = prev is None or (obj == obj and obj > prev * 1.01 and v <= ms["history"][-2]["violation"] + 1e-9 if len(ms["history"]) > 1 else True)
            ms["stall"] = 0 if improved else ms["stall"] + 1
            if ms["stall"] >= 2 or len(ms["current"]) >= MAXVARS:
                ms["done"] = True; log("%s: stop (stall %d, %d vars)" % (model, ms["stall"], len(ms["current"])))
        st["round"] = k + 1
        json.dump(st, open(sp + ".partial", "w"), indent=1); os.replace(sp + ".partial", sp)
    log("FINISHED")
    open(os.path.join(a.state, "DONE"), "w").write(time.strftime("%F %T") + "\n")


if __name__ == "__main__":
    main()
