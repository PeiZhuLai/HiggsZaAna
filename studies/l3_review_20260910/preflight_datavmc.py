#!/usr/bin/env python
"""
Preflight for dataVmc: resolve every job's inputs with the REAL resolver and require
that each resolves to a non-empty list of files that all exist.

Why this has to run before submission: Plot_Helper loads inputs with
_add_file_if_exists(), which SKIPS a missing file silently. A job pointing at a file that
is not there does not fail -- its component just comes out smaller (or zero) in the
plots. On 2026-09-24 the job generator still listed the retired 2022 DYGto2LG_10to50 /
_50to100 slices, which do not exist in the FSR-fix production; nothing downstream would
have said so.

The resolver functions are pulled out of 1_prepare_dataVmc.py by AST (that module runs
argparse at import time) so that this checks the code the jobs actually run, not a
re-implementation of it.
"""
import ast, os, sys
PLOT = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot"
JOBS = sys.argv[1] if len(sys.argv) > 1 else PLOT + "/Condor/dataVmc_jobs.txt"
sys.path.insert(0, PLOT + "/lib")
import Analyzer_Configs as AC
from Plot_Helper import _run3_sources_for_sample

src = open(PLOT + "/scripts/1_prepare_dataVmc.py").read()
wanted = {"_append_unique", "_parse_source_selectors", "_format_source_selectors", "_parse_sample_filter"}
mod = ast.Module(body=[n for n in ast.parse(src).body
                       if isinstance(n, ast.FunctionDef) and n.name in wanted], type_ignores=[])
ns = {}
exec(compile(mod, "1_prepare_dataVmc.py[extract]", "exec"), ns)
missing_fn = wanted - set(ns)
if missing_fn:
    sys.exit("could not extract %s" % missing_fn)

cfg = AC.Analyzer_Config("inclusive", "run3", 1, False)
print("sample_loc:", cfg.sample_loc)
bad, nfiles, rows = [], 0, 0
for line in open(JOBS):
    parts = line.split()
    if len(parts) != 4 or parts[0].startswith("#"):
        continue
    rows += 1
    region, tag, stag, samples = parts
    try:
        selected, filters = ns["_parse_sample_filter"](samples, cfg)
    except Exception as e:
        bad.append((region, stag, "PARSE_ERROR", str(e)[:90])); continue
    cfg.run3_source_filters = filters
    for sample in (selected or []):
        srcs = _run3_sources_for_sample(sample, cfg)
        if not srcs:
            bad.append((region, stag, "EMPTY_SOURCE_LIST", sample)); continue
        for directory, year in srcs:
            p = os.path.join(cfg.sample_loc, directory, "%s.root" % year)
            nfiles += 1
            if not os.path.exists(p):
                bad.append((region, stag, "MISSING_FILE", p))
print("job rows: %d   resolved files: %d   problems: %d" % (rows, nfiles, len(bad)))
for b in bad[:30]:
    print("  %-4s %-26s %-18s %s" % b)
print("VERDICT: %s" % ("PASSED" if not bad and rows > 0 else "FAILED"))
sys.exit(0 if (not bad and rows > 0) else 1)
