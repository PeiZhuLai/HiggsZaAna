#!/usr/bin/env bash
# Read-only: every production scored file (run3_bdt_scored_fsrfix, 1208) must carry the score_provenance fingerprint
# equal to the installed low/high models; no *.w123partial left over; the backup tree holds 1208 files.
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017
X=/eos/user/p/pelai/HZa/root_P2Root
L=$(xrdfs eosuser.cern.ch ls -R $X/run3_bdt_scored_fsrfix 2>/dev/null)
echo "production .root: $(echo "$L" | grep -c '\.root$'), leftover .w123partial: $(echo "$L" | grep -c 'w123partial')"
echo "backup .root: $(xrdfs eosuser.cern.ch ls -R $X/run3_bdt_scored_fsrfix_before_w123_20261017 2>/dev/null | grep -c '\.root$')"
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile && env -i HOME=$HOME PATH=/usr/bin:/bin /eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3 - <<'PY'
import glob, score_provenance as sp
files = sorted(glob.glob("/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/*/*.root"))
bad = []
for p in files:
    try:
        if sp._fingerprint_status(p) != "fingerprint": bad.append(p + " (no fingerprint)")
    except sp.ScoreProvenanceError as e: bad.append(str(e)[:150])
print("fingerprint checked %d files, bad %d" % (len(files), len(bad)))
for b in bad[:10]: print("  ", b)
PY
