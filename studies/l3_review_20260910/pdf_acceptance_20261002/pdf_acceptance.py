#!/usr/bin/env python
"""PDF uncertainty on the SIGNAL ACCEPTANCE, per (m_a, era), separated from the
cross-section part (companion of ../qcd_scale_acceptance.py).

Why the input NanoAOD is needed
    The analysis parquet keeps only the per-event central value and the
    central +- standard deviation over the PDF members (LHEPdfWeight_Unit/Up/Down,
    HiggsDNA/higgs_dna/analysis.py attach_lhe_weights). The per-member ratio
    A_i/A_0 cannot be formed from that, so the individual members are read back
    from the input NanoAOD for the selected events (matched on run, lumi, event).

PDF set
    LHEPdfWeight doc string: "LHA IDs 325500 - 325600" = NNPDF31_nnlo_as_0118_nf_4_mc_hessian,
    101 members, ErrorType symmhessian (LHAPDF .info). Member 0 is the central set,
    members 1..100 are symmetric Hessian eigenvectors; no alpha_s members are stored.
    Prescription: delta = sqrt( sum_{i=1..100} (A_i/A_0 - 1)^2 ), symmetric.

Method (same as the QCD-scale study)
    A_i / A_0 = (sum_sel w * r_i / sum_sel w * r_0) / (sum_gen w_i / sum_gen w_0)
    r_i = LHEPdfWeight[i] (w_var / w_nominal), w = weight_central of the analysis.
    Denominator from the Runs tree: LHEPdfSumw is normalized to genEventSumw, so
    sum_gen w_i = sum_files genEventSumw * LHEPdfSumw[i].

Selections ("variants")
    presel   : every event of merged_nominal.parquet (what the QCD-scale study used)
    wp       : the events entering the datacard = root_MVAcut/sig/mA_M<m>/output_<era>.root,
               nominal ele + mu trees (BDT test split, MVA score above the working point)
    wp_incl  : scored inclusive tree (all splits) with MVA_Score_mA_M<m> > working point
               (3.3x the statistics of "wp"; cross-check of the split choice)

Usage
    python pdf_acceptance.py fetch   <sample> [--nproc 6]   # read NanoAOD, cache per file
    python pdf_acceptance.py compute <sample>               # per-sample result JSON
"""
import argparse, glob, json, os, sys, time, hashlib
import numpy as np

SIG = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
SCORED = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
MVACUT = "/eos/home-p/pelai/HZa/root_MVAcut/sig"
WPJSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
WORK = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(WORK, "cache")
RES = os.path.join(WORK, "results")
NMEM = 101
TREES = ["DiphotonTree/ggh_125_Za_ele_13p6TeV_cat0", "DiphotonTree/ggh_125_Za_mu_13p6TeV_cat0"]


def input_files(sample):
    fs = []
    for cfg in glob.glob(os.path.join(SIG, sample, "job_*", "*_config_job*.json")):
        fs += json.load(open(cfg)).get("files") or []
    return sorted(set(fs))


def cache_path(sample, url):
    h = hashlib.md5(url.encode()).hexdigest()[:8]
    return os.path.join(CACHE, sample, "%s_%s.npz" % (os.path.basename(url)[:-5], h))


def presel_ids(sample):
    import pyarrow.parquet as pq
    t = pq.read_table(os.path.join(SIG, sample, "merged_nominal.parquet"),
                      columns=["run", "luminosityBlock", "event", "weight_central",
                               "LHEPdfWeight_Unit", "LHEPdfWeight_Up", "LHEPdfWeight_Down"])
    return {k: np.asarray(t[k]) for k in t.column_names}


def keyof(run, lumi, evt):
    run = np.asarray(run, dtype=np.uint64); lumi = np.asarray(lumi, dtype=np.uint64)
    evt = np.asarray(evt, dtype=np.uint64)
    # packed key: run < 2^11, lumi < 2^21, event < 2^32 (asserted, never silently truncated)
    assert (run < (1 << 11)).all() and (lumi < (1 << 21)).all() and (evt < (1 << 32)).all(), "key overflow"
    return (run << np.uint64(53)) | (lumi << np.uint64(32)) | evt


# --------------------------------------------------------------------------- fetch
def fetch_one(args):
    sample, url, keys = args
    out = cache_path(sample, url)
    if os.path.exists(out):
        return url, "cached", 0.0
    import uproot, awkward as ak
    t0 = time.time()
    last = None
    lfn = "/store" + url.split("/store", 1)[1]
    redirs = ["root://xrootd-cms.infn.it/", "root://cms-xrd-global.cern.ch/", "root://xrootd-cms.infn.it/",
              "root://cmsxrootd.fnal.gov/"]
    for attempt in range(4):
        try:
            f = uproot.open(redirs[attempt] + lfn, timeout=300)
            R = f["Runs"]
            gsw = np.asarray(R["genEventSumw"].array(library="np"), dtype=float)
            psw = R["LHEPdfSumw"].array(library="ak")
            npsw = ak.to_numpy(ak.num(psw))
            if not (npsw == NMEM).all():
                raise RuntimeError("Runs LHEPdfSumw length %s" % np.unique(npsw))
            psw = ak.to_numpy(psw).astype(float)
            E = f["Events"]
            a = E.arrays(["run", "luminosityBlock", "event", "genWeight", "LHEPdfWeight"])
            k = keyof(a["run"], a["luminosityBlock"], a["event"])
            m = np.isin(k, keys)
            nm = ak.to_numpy(ak.num(a["LHEPdfWeight"]))
            bad = int((nm[m] != NMEM).sum())
            pdf = ak.to_numpy(a["LHEPdfWeight"][m][nm[m] == NMEM]).astype(np.float32) if bad == 0 else None
            if bad:
                raise RuntimeError("%d selected events with nLHEPdfWeight != %d" % (bad, NMEM))
            tmp = out + ".part%d" % os.getpid()
            np.savez_compressed(tmp, key=k[m], run=ak.to_numpy(a["run"][m]), pdf=pdf,
                                genWeight=ak.to_numpy(a["genWeight"][m]),
                                genEventSumw=gsw, LHEPdfSumw=psw, nevents=len(k),
                                sum_genWeight_events=float(ak.sum(a["genWeight"])),
                                n_nmem_bad_all=int((nm != NMEM).sum()))
            os.replace(tmp + ".npz" if not tmp.endswith(".npz") else tmp, out)
            return url, "ok", time.time() - t0
        except Exception as e:  # noqa: BLE001
            last = "%s: %s" % (type(e).__name__, str(e)[:200])
            if not isinstance(e, (OSError, TimeoutError, ImportError, ValueError)):
                break  # programming / content error: do not retry
            time.sleep(10 * (attempt + 1))
    return url, "FAIL " + last, time.time() - t0


def cmd_fetch(sample, nproc):
    from multiprocessing import Pool
    # import the xrootd stack once before forking: a first import inside 6 forked
    # workers at once failed with "Unable to load filesystem ... fsspec_xrootd"
    import uproot, awkward, fsspec_xrootd, XRootD.client  # noqa: F401
    os.makedirs(os.path.join(CACHE, sample), exist_ok=True)
    p = presel_ids(sample)
    keys = np.unique(keyof(p["run"], p["luminosityBlock"], p["event"]))
    files = input_files(sample)
    print("[fetch] %s: %d files, %d preselected events" % (sample, len(files), len(keys)), flush=True)
    nfail = 0
    # hard wall-clock limit: an xrootd read once hung > 15 min despite uproot's timeout;
    # a timed-out file counts as failed, the pool is terminated, run_all.sh retries the sample
    deadline = time.time() + 600 + 120 * len(files) / max(nproc, 1)
    with Pool(nproc) as pool:
        jobs = [(u, pool.apply_async(fetch_one, ((sample, u, keys),))) for u in files]
        for u, j in jobs:
            try:
                url, st, dt = j.get(timeout=max(1.0, deadline - time.time()))
            except Exception as e:  # noqa: BLE001  (multiprocessing.TimeoutError)
                url, st, dt = u, "FAIL timeout %s" % type(e).__name__, -1.0
            if st.startswith("FAIL"):
                nfail += 1
            print("[fetch] %s %-10s %6.1fs %s" % (sample, st[:200], dt, os.path.basename(url)), flush=True)
        pool.terminate()
    print("[fetch] %s done, %d failed" % (sample, nfail), flush=True)
    return nfail


# --------------------------------------------------------------------------- compute
def hessian(acc):
    d = acc[1:] - 1.0
    return float(np.sqrt((d ** 2).sum()))


def wp_cut(ma):
    for r in json.load(open(WPJSON))["results"]:
        if int(r["mA"]) == int(ma):
            return float(r["MVAcut"])
    raise KeyError(ma)


def cmd_compute(sample, nboot=200):
    import uproot
    ma = int(sample.split("_")[1][1:]); era = sample.split("_")[2]
    files = input_files(sample)
    miss = [u for u in files if not os.path.exists(cache_path(sample, u))]
    if miss:
        sys.exit("[compute] %s: %d/%d files not fetched" % (sample, len(miss), len(files)))
    K, P, GW = [], [], []
    G = np.zeros(NMEM); G0 = 0.0; nev = 0; runs = set()
    for u in files:
        z = np.load(cache_path(sample, u))
        K.append(z["key"]); P.append(z["pdf"]); GW.append(z["genWeight"]); runs |= set(np.unique(z["run"]).tolist())
        G += (z["genEventSumw"][:, None] * z["LHEPdfSumw"]).sum(axis=0); G0 += z["genEventSumw"].sum()
        nev += int(z["nevents"])
    K = np.concatenate(K); P = np.concatenate(P).astype(float); GW = np.concatenate(GW)
    if len(np.unique(K)) != len(K):
        sys.exit("[compute] %s: duplicate (lumi,event) keys among matched NanoAOD events" % sample)
    order = np.argsort(K); K = K[order]; P = P[order]; GW = GW[order]
    gen = G / G0
    gen_ratio = gen / gen[0]
    checks = {"n_files": len(files), "n_nano_events": nev, "runs": sorted(runs),
              "gen_member0_over_genEventSumw": float(gen[0])}

    def lookup(run, lumi, evt, label):
        k = keyof(run, lumi, evt)
        idx = np.searchsorted(K, k); idx[idx >= len(K)] = 0
        ok = K[idx] == k
        checks["%s_n" % label] = int(len(k)); checks["%s_unmatched" % label] = int((~ok).sum())
        checks["%s_dupkeys" % label] = int(len(k) - len(np.unique(k)))
        if (~ok).any():
            sys.exit("[compute] %s: %d %s events not found in NanoAOD" % (sample, (~ok).sum(), label))
        return idx

    def measure(idx, w, label, boot=True):
        r = P[idx]
        sel = (w[:, None] * r).sum(axis=0)
        sel_ratio = sel / sel[0]
        acc = sel_ratio / gen_ratio
        out = {"n": int(len(idx)), "sumw": float(w.sum()),
               "acc_ratio": acc.tolist(), "sel_ratio": sel_ratio.tolist(),
               "acc_unc": hessian(acc), "yield_unc": hessian(sel_ratio), "xs_unc": hessian(gen_ratio),
               "acc_max_member_dev": float(np.abs(acc[1:] - 1).max())}
        if boot and nboot:
            rng = np.random.default_rng(12345)
            b = []
            for _ in range(nboot):
                pw = w * rng.poisson(1.0, len(w))
                s = (pw[:, None] * r).sum(axis=0)
                b.append(hessian((s / s[0]) / gen_ratio))
            out["acc_unc_boot_std"] = float(np.std(b)); out["acc_unc_boot_mean"] = float(np.mean(b))
        return out

    res = {"sample": sample, "m_a": ma, "era": era, "pdf_set": "NNPDF31_nnlo_as_0118_nf_4_mc_hessian (325500-325600)",
           "prescription": "symmetric Hessian, sqrt(sum_{i=1}^{100}(A_i/A_0-1)^2)", "gen_ratio": gen_ratio.tolist()}

    # presel
    p = presel_ids(sample)
    idx = lookup(p["run"], p["luminosityBlock"], p["event"], "presel")
    r = P[idx]
    checks["presel_unit_vs_member0_maxabs"] = float(np.abs(p["LHEPdfWeight_Unit"] - r[:, 0]).max())
    sd = r[:, 1:].std(axis=1)
    checks["presel_parquetUp_minus_unit_vs_std_maxabs"] = float(np.abs((p["LHEPdfWeight_Up"] - p["LHEPdfWeight_Unit"]) - sd).max())
    checks["presel_member0_minmax"] = [float(r[:, 0].min()), float(r[:, 0].max())]
    checks["parquet_rms_band_mean"] = float(np.mean(sd))
    res["presel"] = measure(idx, np.asarray(p["weight_central"], dtype=float), "presel")
    # same-sample analogue of what the parquet Up/Down would give (RMS, not Hessian)
    wc = np.asarray(p["weight_central"], dtype=float)
    res["presel"]["parquet_rms_yield_up"] = float((wc * p["LHEPdfWeight_Up"]).sum() / (wc * p["LHEPdfWeight_Unit"]).sum() - 1)

    # wp: datacard events
    fn = os.path.join(MVACUT, "mA_M%d" % ma, "output_%s.root" % era)
    cols = ["run", "luminosityBlock", "event", "weight_central", "weight", "MVA_Score_mA_M%d" % ma]
    with uproot.open(fn) as f:
        arrs = [f[t].arrays(cols, library="np") for t in TREES]
    W = {c: np.concatenate([a[c] for a in arrs]) for c in cols}
    cut = wp_cut(ma)
    checks["wp_file"] = fn; checks["wp_file_mtime"] = time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(fn)))
    checks["wp_json_cut"] = cut; checks["wp_min_score_in_tree"] = float(W[cols[-1]].min())
    checks["wp_n_ele_mu"] = [int(len(a["run"])) for a in arrs]
    idx = lookup(W["run"], W["luminosityBlock"], W["event"], "wp")
    res["wp"] = measure(idx, W["weight_central"].astype(float), "wp")
    res["wp_weight"] = measure(idx, W["weight"].astype(float), "wp_weight", boot=False)

    # wp_incl: all splits with the same cut
    fn2 = os.path.join(SCORED, "mA_M%d" % ma, "%s.root" % era)
    with uproot.open(fn2) as f:
        sc = cols[-1]
        a = f["inclusive"].arrays(["run", "luminosityBlock", "event", "weight_central", sc], library="np")
        tt = f["test"].arrays([sc], library="np")
    checks["test_pass_cut"] = int((tt[sc] > cut).sum())
    checks["wp_tree_vs_test_pass_cut_diff"] = int(sum(checks["wp_n_ele_mu"]) - checks["test_pass_cut"])
    m = a[sc] > cut
    idx = lookup(a["run"][m], a["luminosityBlock"][m], a["event"][m], "wp_incl")
    res["wp_incl"] = measure(idx, a["weight_central"][m].astype(float), "wp_incl")
    res["checks"] = checks
    os.makedirs(RES, exist_ok=True)
    out = os.path.join(RES, "%s.json" % sample)
    json.dump(res, open(out, "w"), indent=1)
    print("[compute] %s  acc(presel)=%.3f%%  acc(wp)=%.3f%% +- %.3f(stat)  acc(wp_incl)=%.3f%%  yield(wp)=%.3f%%  xs=%.3f%%"
          % (sample, 100 * res["presel"]["acc_unc"], 100 * res["wp"]["acc_unc"], 100 * res["wp"].get("acc_unc_boot_std", 0),
             100 * res["wp_incl"]["acc_unc"], 100 * res["wp"]["yield_unc"], 100 * res["wp"]["xs_unc"]), flush=True)
    print("[compute] checks:", json.dumps(checks), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["fetch", "compute", "both"])
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--nproc", type=int, default=6)
    ap.add_argument("--nboot", type=int, default=200)
    o = ap.parse_args()
    rc = 0
    for s in o.samples:
        if o.cmd in ("fetch", "both"):
            if cmd_fetch(s, o.nproc):
                rc = 1
                continue
        if o.cmd in ("compute", "both"):
            cmd_compute(s, o.nboot)
    sys.exit(rc)
