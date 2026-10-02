"""Write the catalogs + configs for the 2024 DY+jets re-production with the MC overlap veto.

Input : HiggsDNA/metadata/samples/zgamma_tutorial.json (DYJetsTo2E/2Mu/2Tau entries, unchanged)
        HiggsDNA/metadata/za_bkgmc_run3.json (tagger/systematics config of the FSR-fix Bkg_MC production)
        old production job configs (parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC/<s>_2024/job_N) for the test files
Output: HiggsDNA/metadata/samples/za_dyveto_2024.json        (full: DAS dataset names, same as zgamma_tutorial)
        HiggsDNA/metadata/samples/za_dyveto_2024_test.json   (3 explicit files/sample, global redirector)
        HiggsDNA/metadata/za_bkgmc_dyveto2024.json / _test.json (za_bkgmc_run3.json with only the sample block changed)
        test_files.json (sample -> [(old job, old rows, file)]) for the acceptance script

Separate catalog files: SampleManager rewrites <catalog>_sample_manager_full.json, so running on
zgamma_tutorial.json would overwrite the validated production's resolved file list.
Global redirector: INFN failed for some 2024 files on 2026-09-27 (v3v8 study, TTGG job_23).
"""
import json, copy
R = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/metadata/"
S = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/dy2024_reprod_20260927/"
OLD = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC/"
GLOBAL = "root://cms-xrd-global.cern.ch/"
SAMPLES = ["DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
TEST_JOBS = {"DYJetsTo2E": [3, 5, 6], "DYJetsTo2Mu": [1, 2, 3], "DYJetsTo2Tau": [557, 130, 356]}

cat = json.load(open(R + "samples/zgamma_tutorial.json"))
full = {s: {"xs": {"2024": cat[s]["xs"]["2024"]}, "files": {"2024": cat[s]["files"]["2024"]}} for s in SAMPLES}
for s in SAMPLES:
    for k in cat[s]:
        if k not in ("xs", "files", "fpo"):
            full[s][k] = cat[s][k]
json.dump(full, open(R + "samples/za_dyveto_2024.json", "w"), indent=4)

import pyarrow.parquet as pq, os
test, tf = {}, {}
for s in SAMPLES:
    files = []
    tf[s] = []
    for j in TEST_JOBS[s]:
        c = json.load(open(f"{OLD}{s}_2024/job_{j}/{s}_2024_config_job{j}.json"))
        assert len(c["files"]) == 1
        f = c["files"][0].replace("root://xrootd-cms.infn.it/", GLOBAL)
        p = f"{OLD}{s}_2024/job_{j}/output_job_{j}_nominal.parquet"
        n = pq.ParquetFile(p).metadata.num_rows if os.path.exists(p) else 0
        files.append(f); tf[s].append([j, n, f])
    test[s] = {"_comment": "SMALL-BATCH TEST (3 files) of %s 2024 with the MC overlap veto" % s,
               "xs": {"2024": cat[s]["xs"]["2024"]}, "files": {"2024": files}}
json.dump(test, open(R + "samples/za_dyveto_2024_test.json", "w"), indent=4)
json.dump(tf, open(S + "test_files.json", "w"), indent=2)

base = json.load(open(R + "za_bkgmc_run3.json"))
for suf, catname in (("", "za_dyveto_2024.json"), ("_test", "za_dyveto_2024_test.json")):
    c = copy.deepcopy(base)
    c["samples"] = {"catalog": "metadata/samples/" + catname, "sample_list": SAMPLES, "years": ["2024"]}
    json.dump(c, open(R + "za_bkgmc_dyveto2024%s.json" % suf, "w"), indent=4)
print(json.dumps(tf, indent=1))
