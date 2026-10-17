#!/usr/bin/env bash
# Test of Parque2Root_BDT.py with local write + verify + xrdcp + model fingerprint: score mA_M3 2022preEE (the
# production job's input and --split) into a sandbox EOS dir with the scoring env (higgs-alp-ana, as on condor), then
# (1) compare every branch common with the production scored file, (2) the provenance check must pass through the
# fingerprint, (3) a copy with an altered fingerprint must be refused. Production files are only read.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/prod_mvacut_integrity_20261017
PQ=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile
IN=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC/mA_M3_2022preEE/merged_nominal.parquet
OUT=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/_scorer_test_20261017/mA_M3/2022preEE.root
PROD=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M3/2022preEE.root
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
s=$(date +%s)
( cd $PQ && env -i HOME=$HOME PATH=/usr/bin:/bin $PY Parque2Root_BDT.py -i $IN -o $OUT --split ) > $D/logs/test_scorer.log 2>&1
echo "scorer rc=$? $(( $(date +%s)-s )) s"
cd $PQ && env -i HOME=$HOME PATH=/usr/bin:/bin $PY - $OUT $PROD $D <<'PYC'
import sys, json, os, uproot, numpy as np
import score_provenance as sp
out, prod, D = sys.argv[1:4]
fo, fp = uproot.open(out), uproot.open(prod)
to = {k for k, v in fo.items(cycle=False) if hasattr(v, "num_entries")}
tp = {k for k, v in fp.items(cycle=False) if hasattr(v, "num_entries")}
print("trees new %s prod %s" % (sorted(to), sorted(tp)))
nd = nc = 0
for t in sorted(to & tp):
    ko, kp = set(fo[t].keys()), set(fp[t].keys())
    if t == "inclusive": print("branches only new:", sorted(ko - kp)[:10], " only prod:", sorted(kp - ko)[:10])
    if fo[t].num_entries != fp[t].num_entries: print("ENTRIES differ", t); nd += 1; continue
    for b in sorted(ko & kp):
        x = fo[t][b].array(library="np"); y = fp[t][b].array(library="np"); nc += 1
        if x.dtype == object:
            same = [len(v) for v in x] == [len(v) for v in y] and np.array_equal(np.concatenate(list(x) + [np.zeros(0)]).astype(float), np.concatenate(list(y) + [np.zeros(0)]).astype(float), equal_nan=True)
        else:
            same = x.dtype == y.dtype and np.array_equal(x, y, equal_nan=x.dtype.kind == "f")
        if not same: nd += 1; print("DIFF", t, b)
print("common branches compared %d, differ %d" % (nc, nd))
print("fingerprint key present:", sp.FINGERPRINT_KEY in fo)
print("provenance:", sp.check_scored_file(out))
fake = D + "/fake_fp.root"
fp_ = json.loads(str(fo[sp.FINGERPRINT_KEY])); fp_["low"]["md5"] = "0" * 32
with uproot.recreate(fake) as f:
    f["inclusive"] = fo["inclusive"].arrays(["H_mass"], library="pd", entry_stop=10); sp.write_fingerprint(f, fp_)
try: sp.check_scored_file(fake); print("altered fingerprint NOT refused")
except sp.ScoreProvenanceError as e: print("altered fingerprint refused:", str(e)[-80:])
os.remove(fake)
PYC
