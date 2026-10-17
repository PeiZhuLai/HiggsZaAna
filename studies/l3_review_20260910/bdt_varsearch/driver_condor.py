#!/usr/bin/env python3
"""Forward-selection driver for the HZa BDT input-variable search on CERN condor (one cluster per round).

Runs on lxplus under systemd-run (python3 of the system is enough), one process for both models. Each round, for each model that
is still growing, every pool variable not yet selected is added to the current set as one candidate;
the candidates of a round go out as ONE condor cluster (run_candidate.sh -> search_job.py). When all results are in, the
best candidate is kept:
  gates (every anchor mass of the model, see gate_violation): R, cheap MC closure pulls, data-sideband
          smooth-fit GOF; m_a = 1 also: overlap share <= preselection and window fraction >= preselection
          (no BDT-enriched second population)
  among candidates with total violation <= TOL (0.5): highest objective = mean over masses of Z / Z_production
          (Z_production from the stored production scores at the same data-sideband count)
  if none passes: smallest total gate violation (sum of excess sigmas / log10 p deficits)
A model stops when for two rounds neither the gate violation drops nor the objective improves by >= 1%, or at 12 variables.
Round 0 also runs the references: stored production scores, production variables retrained, and the
sandbox win_mfree (low) / shared-definition (high) models.

State (AFS): <state>/round_<k>/{config.jsonl, args.txt, search.sub, logs/}, state.json, heartbeat.txt;
results (EOS): <results>/round_<k>/<id>.json
Restartable: completed rounds are read back from state.json; a round whose results all exist is not resubmitted.
Usage: systemd-run --user ... python3 driver_condor.py --state <afs dir> --results <eos dir> --samples <eos dir> --code <afs dir>
"""
import argparse
import json
import math
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import search_job  # noqa: E402  (gate-B fit, shared with the jobs)

POOL = ["pho1Pt", "pho2Pt", "pho_pt_asym", "pho1R9", "pho2R9", "pho1IetaIeta55", "pho2IetaIeta55",
        "pho1PIso_noCorr", "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "ALP_pt", "ALP_eta",
        "var_PtaOverMa", "var_dR_g1g2", "var_dR_Za", "var_dR_g1Z", "var_dEta_g1Z", "var_dPhi_g1Z",
        "min_lg_deltaR", "max_lg_deltaR", "H_pt", "pho1_relEerr", "pho2_relEerr"] + \
       ["ALP_%s_photon_%s" % (w, v) for w in ("lead", "sublead") for v in
        ("eta", "hoe", "s4", "sieip", "etaWidth", "phiWidth", "esEffSigmaRR", "chiso",
         "hcalPFClusterIso", "pfChargedIsoWorstVtx", "trkSumPtHollowConeDR03")] + ["ALP_sublead_photon_lep_near_dR"]
# Left out by choice: photon mvaID (the analysis uses cut-based photon ID).
# Left out as exact duplicates of pool variables (checked event by event): ALP_PhotonIso = ALP_calculatedPhotonIso,
# ALP_{lead,sublead}_photon_ecalPFClusterIso = pho{1,2}PIso_noCorr, ALP_lead_photon_lep_near_dR = min_lg_deltaR.
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
TOL = 0.5   # a candidate passes if its total gate violation (sigmas / log10 p beyond the thresholds) is <= TOL


def log(msg):
    print(time.strftime("%F %T"), msg, flush=True)


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def gate_violation(res):
    """Sum of gate excesses over the anchor masses (0 = every gate passed):
      C  |R - 1| <= 2 sigma (sigma floor 0.05)
      L  one-sided: R_lo >= 1 - 2 sigma. A 105-115 depletion moves the turn-on toward the blinded window, where the
         erf-step families cannot follow it; the sandbox win_mfree (mA2 R_lo 0.46 +- 0.16) failed bias with -0.29
         while passing A and B. R_lo > 1 (production high, 2-3) is not penalized.
      A  cheap MC closure: |pull| <= 2 in 115-125, 120-130, 125-135 (passing DY at data-level statistics,
         log-polynomial fitted to the sidebands and extrapolated; search_job.closure_proxy)
      B  passing data sidebands described by a smooth log-polynomial: GOF p >= 0.01 / N_masses
         (1% for the whole model; production high mA15 has p = 0.005 from one 100-105 GeV bin)
      D  m_a = 1: overlap share <= preselection, window fraction >= preselection
    R_lo, R_hi and the data shape relative to preselection are recorded only: the production high-mass BDT moves the
    background smoothly to low m_llgg (R_lo 2-3) and still passes closure and bias at all 30 points."""
    v = 0.0
    pmin = 0.01 / len(res["per_mass"])
    for m, r in res["per_mass"].items():
        x, e = r["R"], max(r["R_err"], 0.05)
        v += max(0.0, abs(x - 1) / e - 2) if finite(x) and finite(e) else 10
        x, e = r["R_lo"], max(r["R_lo_err"], 0.05)
        v += max(0.0, (1 - x) / e - 2) if finite(x) and finite(e) else 10
        for pu in r.get("A_pulls", [float("nan")]):
            v += max(0.0, abs(pu) - 2) if finite(pu) else 10
        p = r.get("B_p", float("nan"))
        v += max(0.0, math.log10(pmin / max(p, 1e-300))) if finite(p) else 10
        if m == "1":
            v += max(0.0, (r["overlap"] - r["presel_overlap"]) / 0.02)
            v += max(0.0, (r["presel_win_frac"] - r["win_frac"]) / 0.02)
    return v


def z_fit(r):
    """Proxy Z = S_win / sqrt(B_win) with B_win from the smooth fit of the passing data sidebands (gate B fit, refitted
    from the stored B_obs) integrated over 115-133. The DY-ratio estimate stored as r["Z"] breaks down when few DY
    events pass (mA4: B_win 0.008 from 0.4 effective events -> Z 24x production)."""
    import numpy as np
    y = np.asarray(r["B_obs"], float)
    xc = search_job.SB_DATA_LO + 2.5
    order, (b, cov, dev) = search_job.ftest_order(xc, y)
    if b is None:
        return float("nan")
    xs = np.arange(115.5, 133.0, 1.0)
    B = float(np.exp(search_job.design(xs, order) @ b).sum() / 5.0)   # fit is counts per 5 GeV
    return r["S_win"] / math.sqrt(B) if B > 0 else float("nan")


def objective(res, zref):
    """Mean over masses of Z / Z_production (z_fit); -inf if any mass has no finite Z."""
    zs = [z_fit(res["per_mass"][m]) / zref[m] for m in res["per_mass"] if zref.get(m)]
    if not zs or not all(finite(z) for z in zs):
        return float("-inf")
    return sum(zs) / len(zs)


def submit(rdir, cfgs, a):
    """Write config + condor submit file for one round (one cluster); returns the ClusterId or None."""
    os.makedirs(os.path.join(rdir, "logs"), exist_ok=True)
    res = os.path.join(a.results, os.path.basename(rdir))
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(rdir, "config.jsonl"), "w") as f:
        for c in cfgs:
            f.write(json.dumps(c) + "\n")
    with open(os.path.join(rdir, "args.txt"), "w") as f:
        for i, c in enumerate(cfgs):
            f.write("%d %s\n" % (i, c["model"]))
    sub = os.path.join(rdir, "search.sub")
    bn = "hza_varsearch_" + os.path.basename(rdir)
    with open(sub, "w") as f:
        f.write("""universe              = vanilla
executable            = {code}/run_candidate.sh
arguments             = $(IDX) {res} $(MODEL)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files  = {code}/search_job.py, {rdir}/config.jsonl
transfer_output_files = ""
output                = {rdir}/logs/job_$(IDX).out
error                 = {rdir}/logs/job_$(IDX).err
log                   = {rdir}/logs/cluster.log
request_cpus          = 1
request_memory        = 2000
request_disk          = 4000000
+JobFlavour           = "longlunch"
batch_name            = {bn}
on_exit_hold          = (ExitBySignal == True) || (ExitCode != 0)
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)
queue IDX, MODEL from {rdir}/args.txt
""".format(code=a.code, res=res, rdir=rdir, bn=bn))
    out = subprocess.run(["condor_submit", sub], capture_output=True, text=True, cwd=rdir)
    open(os.path.join(rdir, "submit.log"), "a").write(out.stdout + out.stderr)
    m = re.search(r"submitted to cluster (\d+)", out.stdout)
    log("submit %s: %d jobs -> cluster %s" % (bn, len(cfgs), m.group(1) if m else "NONE (" + out.stderr[-200:] + ")"))
    return m.group(1) if m else None


def queue_status(clusters):
    """{ 'cluster.proc': (JobStatus, NumJobStarts, MemoryUsage) } for the given clusters; None if condor_q failed."""
    if not clusters:
        return {}
    out = subprocess.run(["condor_q"] + sorted(clusters) + ["-af", "ClusterId", "ProcId", "JobStatus", "NumJobStarts", "MemoryUsage"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        return None
    st = {}
    for line in out.stdout.splitlines():
        f = line.split()
        if len(f) >= 4:
            st["%s.%s" % (f[0], f[1])] = (int(f[2]), int(f[3]) if f[3].isdigit() else 0, f[4] if len(f) > 4 else "")
    return st


def handle_held(st):
    """Held after 3 starts = the job itself fails; remove it (the driver resubmits missing results later)."""
    dead = [j for j, (s, n, _) in st.items() if s == 5 and n >= 3]
    if dead:
        subprocess.run(["condor_rm"] + dead, capture_output=True, text=True)
        log("removed %d jobs held after 3 starts: %s" % (len(dead), dead[:10]))


def wait_round(rdir, cfgs, a, cluster):
    want = {c["id"] for c in cfgs}
    res = os.path.join(a.results, os.path.basename(rdir))
    clusters = {cluster} if cluster else set()
    tries, maxmem = 0, 0
    t_last_submit = time.time()
    while True:
        have = {f[:-5] for f in os.listdir(res) if f.endswith(".json") and not f.startswith(".")}
        st = queue_status(clusters)
        cnt = {}
        for v in (st or {}).values():
            cnt[v[0]] = cnt.get(v[0], 0) + 1
            if v[2].isdigit():
                maxmem = max(maxmem, int(v[2]))
        open(os.path.join(a.state, "heartbeat.txt"), "w").write(
            "%s %s results %d/%d  queue: %d running %d idle %d held  max MemoryUsage %d MB\n"
            % (time.strftime("%F %T"), os.path.basename(rdir), len(have & want), len(want),
               cnt.get(2, 0), cnt.get(1, 0), cnt.get(5, 0), maxmem))
        if want <= have:
            log("%s complete (max MemoryUsage %d MB)" % (os.path.basename(rdir), maxmem))
            return res
        if st is None:
            time.sleep(120); continue
        handle_held(st)
        live = [j for j, v in st.items() if not (v[0] == 5 and v[1] >= 3)]
        if not live and time.time() - t_last_submit > 600:
            missing = [c for c in cfgs if c["id"] not in have]
            if tries >= 3:
                log("giving up on %d candidates after %d resubmissions: %s" % (len(missing), tries, [c["id"] for c in missing][:5]))
                return res
            tries += 1
            rr = rdir + "_retry%d" % tries
            log("%s: queue empty, %d results missing -> resubmission %d" % (os.path.basename(rdir), len(missing), tries))
            cl = submit(rr, missing, a)
            t_last_submit = time.time()
            if cl:
                clusters.add(cl)
        # results of retries land in their own EOS dir; move them into the round's results
        for k in range(1, tries + 1):
            rres = os.path.join(a.results, os.path.basename(rdir) + "_retry%d" % k)
            if os.path.isdir(rres):
                for f in os.listdir(rres):
                    if f.endswith(".json") and not f.startswith(".") and not os.path.exists(os.path.join(res, f)):
                        os.replace(os.path.join(rres, f), os.path.join(res, f))
        time.sleep(90)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--code", required=True)
    ap.add_argument("--results", required=True, help="EOS dir for the result JSONs (one subdir per round)")
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
        cpath = os.path.join(rdir, "config.jsonl")
        if os.path.exists(cpath):
            # a round run for another selection path (e.g. after a change of the selection rule): keep it only if it
            # contains every candidate of this path, else move it and all later rounds aside and submit afresh
            have = {(c["id"], tuple(c["feats"])) for c in (json.loads(l) for l in open(cpath) if l.strip())}
            if not all((c["id"], tuple(c["feats"])) in have for c in cfgs):
                tag = "superseded_" + time.strftime("%Y%m%d_%H%M%S")
                for base in (a.state, a.results):
                    os.makedirs(os.path.join(base, tag), exist_ok=True)
                    for d in sorted(os.listdir(base)):
                        mm = re.match(r"round_(\d+)", d)
                        if mm and int(mm.group(1)) >= k:
                            os.replace(os.path.join(base, d), os.path.join(base, tag, d))
                log("round_%02d was run for another selection path -> moved rounds >= %d to %s" % (k, k, tag))
        if not os.path.exists(cpath):
            cl = submit(rdir, cfgs, a)
            if cl:
                open(os.path.join(rdir, "cluster.txt"), "w").write(cl + "\n")
        elif os.path.exists(os.path.join(rdir, "cluster.txt")):
            cl = open(os.path.join(rdir, "cluster.txt")).read().strip()
        rres = wait_round(rdir, cfgs, a, cl)
        res = {}
        for c in cfgs:
            p = os.path.join(rres, c["id"] + ".json")
            if os.path.exists(p):
                res[c["id"]] = json.load(open(p))
        if k == 0:
            st["zref"] = {model: {m: z_fit(r) for m, r in res["ref_%s_stored" % model]["per_mass"].items()}
                          for model in ("low", "high") if "ref_%s_stored" % model in res}
            st["refs"] = {i: {"violation": gate_violation(res[i]), "objective": objective(res[i], st["zref"]["low" if "_low_" in i else "high"])}
                          for i in res if i.startswith("ref_")}
            log("references: " + json.dumps(st["refs"]))
        for model, ms in st["models"].items():
            if ms["done"]:
                continue
            cand = [(i, r) for i, r in res.items() if i.startswith(model + "_r")]
            cand = [(i, r) for i, r in cand if r["config"]["feats"][:-1] == ms["current"]]   # result of this path
            if not cand:
                log("%s: no results this round -> stop" % model); ms["done"] = True; continue
            scored = [(gate_violation(r), objective(r, st["zref"][model]), i, r) for i, r in cand]
            scored = [s for s in scored if finite(s[1])] or scored   # a Z that cannot be computed is never selected
            passing = [s for s in scored if s[0] <= TOL]
            best = max(passing, key=lambda s: s[1]) if passing else min(scored, key=lambda s: (s[0], -s[1]))
            v, obj, i, r = best
            prev = ms["history"][-1]["objective"] if ms["history"] else None
            ms["current"] = r["config"]["feats"]
            ms["history"].append({"round": k, "added": ms["current"][-1], "feats": ms["current"], "violation": v,
                                  "objective": obj, "n_passing": len(passing), "id": i})
            log("%s round %d: + %s -> obj %.3f viol %.2f (%d/%d candidates pass)" % (model, k, ms["current"][-1], obj, v, len(passing), len(scored)))
            # progress: within tolerance, a >= 1% higher objective; above it, fewer gate violations
            pv = ms["history"][-2]["violation"] if len(ms["history"]) > 1 else None
            if prev is None:
                improved = True
            elif v <= TOL:
                improved = pv > TOL or (finite(obj) and obj > prev * 1.01)
            else:
                improved = v < pv - 1e-6
            ms["stall"] = 0 if improved else ms["stall"] + 1
            if ms["stall"] >= 2 or len(ms["current"]) >= MAXVARS:
                ms["done"] = True; log("%s: stop (stall %d, %d vars)" % (model, ms["stall"], len(ms["current"])))
        st["round"] = k + 1
        json.dump(st, open(sp + ".partial", "w"), indent=1); os.replace(sp + ".partial", sp)
    log("FINISHED")
    open(os.path.join(a.state, "DONE"), "w").write(time.strftime("%F %T") + "\n")


if __name__ == "__main__":
    main()
