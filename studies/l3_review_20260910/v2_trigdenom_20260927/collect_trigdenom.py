#!/usr/bin/env python3
"""
V2 study: collect the trigeff payloads (nominal syst) of the dedicated m_a = 1 GeV signal run
into cutflow_Sig_MC_mA_M1_<era>.json files that plot_trigEffVlepPt.py can read (--in-dir).

Why not scripts/5_collect_cutflow.py main(): its base dir is hard-coded to eos_logs/{Sig_MC,...}
(CUTFLOW_BASEDIR is not used by main()), and it writes into Plot/output/cutflow_list, which is
the input of the current AN figures. Here we import its parsers and only change the paths.

Completeness gate: every job_* directory must have a chosen .out that contains the nominal
trigeff payloads; otherwise the era is reported as INCOMPLETE and not written (unless --allow-partial).

Usage:
  python collect_trigdenom.py --logs <HiggsDNA>/eos_logs/Sig_MC_trigdenomV2 --out <dir> [--eras ...]
"""
import argparse
import importlib.util
import json
import os
import sys
from collections import OrderedDict

REPO = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]


def load_collector():
    spec = importlib.util.spec_from_file_location("cc5", os.path.join(REPO, "scripts", "5_collect_cutflow.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def collect_dedup(cc, files):
    """Like 5_collect_cutflow.collect_trigeff, but each distinct nominal payload is counted once
    per file. The HiggsDNA job log prints every CutFlow JSON line twice (the whole block of
    syst payloads appears twice in each .out, verified on job 12757239.0: n_events=15000 and
    two identical 'nominal' payloads per prefix), so collect_trigeff doubles in_bin and
    pass_trigger. Efficiencies are unaffected, but the binomial errors computed from the
    doubled counts would be too small by sqrt(2). Identical duplicates are dropped; payloads
    that differ (e.g. separate chunks of one file) are all kept."""
    acc = OrderedDict()
    n_dups = 0
    for fp in files:
        seen = set()
        txt = open(fp, errors="ignore").read()
        for p in cc.iter_cutflow_payloads_from_text(txt, syst="nominal"):
            prefix = p.get("cut_type_prefix")
            if not prefix:
                continue
            key = json.dumps(p, sort_keys=True)
            if key in seen:
                n_dups += 1
                continue
            seen.add(key)
            d = acc.setdefault(prefix, {"bins": OrderedDict(), "overall": {"in_bin": 0.0, "pass_trigger": 0.0}})
            d["overall"]["in_bin"] += float(p["overall"]["in_bin"])
            d["overall"]["pass_trigger"] += float(p["overall"]["pass_trigger"])
            for b in p.get("bins", []):
                bb = d["bins"].setdefault(b["ptbin"], {"in_bin": 0.0, "pass_trigger": 0.0})
                bb["in_bin"] += float(b["in_bin"])
                bb["pass_trigger"] += float(b["pass_trigger"])
    out = OrderedDict()
    for prefix, d in sorted(acc.items()):
        o = d["overall"]
        out[prefix] = {
            "overall": {**o, "eff": (o["pass_trigger"] / o["in_bin"]) if o["in_bin"] > 0 else 0.0},
            "bins": OrderedDict((k, {**v, "eff": (v["pass_trigger"] / v["in_bin"]) if v["in_bin"] > 0 else 0.0})
                                for k, v in d["bins"].items()),
        }
    print(f"    dropped {n_dups} duplicated payload lines")
    return cc.TrigEffResult(by_prefix=out, n_jobs_used=len(files), out_files_used=list(files))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="mA_M1")
    ap.add_argument("--eras", default=",".join(ERAS))
    ap.add_argument("--allow-partial", action="store_true")
    args = ap.parse_args()

    cc = load_collector()
    os.makedirs(args.out, exist_ok=True)
    bad = 0
    for era in args.eras.split(","):
        sdir = os.path.join(args.logs, f"{args.sample}_{era}")
        if not os.path.isdir(sdir):
            print(f"[MISSING] {sdir}")
            bad += 1
            continue
        job_dirs = sorted(d for d in os.listdir(sdir) if d.startswith("job_"))
        outs = list(cc.iter_job_out_files(sdir))
        good, empty = [], []
        for fp in outs:
            txt = open(fp, errors="ignore").read()
            has = any(p.get("cut_type_prefix") for p in cc.iter_cutflow_payloads_from_text(txt, syst="nominal"))
            (good if has else empty).append(fp)
        status = "OK" if (len(good) == len(job_dirs) and not empty) else "INCOMPLETE"
        print(f"[{status}] {era}: job_dirs={len(job_dirs)} out_chosen={len(outs)} with_payload={len(good)} without={len(empty)}")
        for fp in empty[:5]:
            print("    no payload:", fp)
        if status != "OK":
            bad += 1
            if not args.allow_partial:
                continue
        res = collect_dedup(cc, good)
        payload = OrderedDict(
            meta=OrderedDict(
                source="V2 trigdenom study 2026-09-27",
                logs=sdir,
                n_job_dirs=len(job_dirs),
                n_jobs_used=res.n_jobs_used,
                status=status,
                out_files_used=res.out_files_used,
            ),
            trigeff_nominal=cc.render_trigeff_json_dict(res),
        )
        fn = os.path.join(args.out, f"cutflow_Sig_MC_{args.sample}_{era}.json")
        with open(fn, "w") as f:
            json.dump(payload, f, indent=1)
        print(f"    -> {fn}  prefixes={len(res.by_prefix)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
