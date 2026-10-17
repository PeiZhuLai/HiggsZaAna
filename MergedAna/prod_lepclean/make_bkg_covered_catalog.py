"""Catalog of the 2024 DY backgrounds restricted to NanoAOD files whose MLNanoAOD
friends are ALL present, plus a matching base config.

Why: the bkg MLNanoAOD is incomplete (88-99.7% of the parent MiniAOD files have a
friend). With the strict loader a job fails if any of its files has a parent
without a friend; without it those events silently look like "no MLPhoton" and
the MC yield is biased low. Running only fully covered files keeps the selection
unbiased; the dataVmc normalisation then needs the coverage factor
sumw(all files) / sumw(processed files) per sample (add_merged_flag.py).

Writes (HiggsDNA repo, never touching the shared zgamma_tutorial.json):
    metadata/samples/za_merged_bkgmc_2024_covered.json
    metadata/za_merged_bkgmc_run3_covered.json
and prints the per-sample coverage.

    python make_bkg_covered_catalog.py [--redirector root://cms-xrd-global.cern.ch/]
"""
import argparse
import json
import os
import re

REPO = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
MLNANO = "/eos/cms/store/group/phys_susy/pelai/HZa_merged/MLNanoAOD"
SAMPLES = {"DYJetsTo2E": "Bkg_DYJetsTo2E_2024", "DYJetsTo2Mu": "Bkg_DYJetsTo2Mu_2024",
           "DYJetsTo2Tau": "Bkg_DYJetsTo2Tau_2024", "DYGto2LG_10to100": "Bkg_DYGto2LG_10to100_2024"}
UUID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--redirector", default="root://cms-xrd-global.cern.ch/")
    a = ap.parse_args()

    base_cat = json.load(open(f"{REPO}/metadata/samples/zgamma_tutorial.json"))
    base_cat = base_cat.get("samples", base_cat)
    out_cat = {}
    for sample, tag in SAMPLES.items():
        pm = json.load(open(f"{REPO}/metadata/parent_maps/{tag}.json"))
        have = {m.group(0) for f in os.listdir(f"{MLNANO}/{tag}") for m in [UUID.search(f)] if m}
        covered = sorted(lfn for lfn, parents in pm.items()
                         if parents and all(UUID.search(p) and UUID.search(p).group(0) in have
                                            for p in (parents if isinstance(parents, list) else [parents])))
        entry = {k: v for k, v in base_cat[sample].items() if k != "files"}
        entry["files"] = {"2024": [a.redirector + lfn for lfn in covered]}
        out_cat[sample] = entry
        print(f"{sample:18s} nano files {len(pm):5d}  fully covered {len(covered):5d} ({len(covered)/len(pm):.1%})")

    cat_rel = "metadata/samples/za_merged_bkgmc_2024_covered.json"
    json.dump(out_cat, open(f"{REPO}/{cat_rel}", "w"), indent=2)
    cfg = json.load(open(f"{REPO}/metadata/za_merged_bkgmc_run3.json"))
    cfg["samples"]["catalog"] = cat_rel
    json.dump(cfg, open(f"{REPO}/metadata/za_merged_bkgmc_run3_covered.json", "w"), indent=2)
    print("wrote", cat_rel, "and metadata/za_merged_bkgmc_run3_covered.json")


if __name__ == "__main__":
    main()
