#!/usr/bin/env bash
#
# Ask a running condor job's own worker node whether it can open a CMS file.
#
# WHY this exists as a script
#   Three earlier attempts at this probe were invalid and each one produced a
#   confident-looking but meaningless answer:
#     1. condor_ssh_to_job's shell inherits the SUBMIT environment (getenv=True),
#        not the wrapper's, so X509_USER_PROXY pointed at the lxplus /tmp path
#        that does not exist on the worker -- the probe ran with no proxy at all.
#     2. same again.
#     3. the probe looked for the job config in $_CONDOR_SCRATCH_DIR, but configs
#        live in the /eos output dir, so it silently fell back to a placeholder
#        path and "proved" that a NON-EXISTENT file fails to open.
#   So: the proxy is set explicitly, the file is one verified to exist (it opens
#   in 2.9 s from lxplus, 989393 entries), and NOTHING is filtered out of the
#   output -- a failed ssh and a hung open have to be distinguishable.
#
# USAGE
#   bash probe_worker_xrootd.sh            # dry run: prints what it would do
#   DRY_RUN=0 bash probe_worker_xrootd.sh  # probe
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
OUT="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/probe_worker.log"
F="/store/data/Run2023C/EGamma1/NANOAOD/22Sep2023_v4-v1/2540000/8ba7b401-41ef-4f68-b26d-34a3f633f7d6.root"

read -r -d '' PY <<PYEOF
import time, uproot, sys, os
print("PROBE proxy=%s exists=%s" % (os.environ.get("X509_USER_PROXY"), os.path.exists(os.environ.get("X509_USER_PROXY",""))))
sys.stdout.flush()
for red in ["root://xrootd-cms.infn.it/","root://cms-xrd-global.cern.ch/"]:
    t=time.time()
    try:
        f=uproot.open(red+"$F", timeout=40)
        print("PROBE %s OK %.1fs n=%d"%(red,time.time()-t,f["Events"].num_entries))
    except Exception as e:
        print("PROBE %s FAIL %.1fs %s %s"%(red,time.time()-t,type(e).__name__,str(e)[:80]))
    sys.stdout.flush()
print("PROBE DONE")
PYEOF
B64=$(printf '%s' "$PY" | base64 -w0)

# Not every running job accepts condor_ssh_to_job ("does not support remote
# access" -- it depends on how the starter was configured on that node), so try
# candidates until one connects. That refusal is NOT a probe result and must not
# be mistaken for one.
cands=$(condor_q -run -af ClusterId ProcId JobBatchName 2>/dev/null | awk '$3 ~ /^(DY|Data_20)/ {print $1"."$2}' | shuf | head -8)
echo "candidates: $(wc -w <<< "$cands")"
[ -z "${cands// }" ] && { echo "no running production job to probe"; exit 1; }
if [ "$DRY_RUN" = "1" ]; then echo "would probe one of them with the verified file $F"; exit 0; fi

for id in $cands; do
  echo "=== $(date '+%F %T') probing $id ===" | tee -a "$OUT"
  out=$(timeout 200 condor_ssh_to_job "$id" \
    "export X509_USER_PROXY=\$_CONDOR_SCRATCH_DIR/x509up_u175325; echo $B64 | base64 -d > \$_CONDOR_SCRATCH_DIR/probe.py; timeout 130 \$_CONDOR_SCRATCH_DIR/hza_ana_env/bin/python -s \$_CONDOR_SCRATCH_DIR/probe.py; echo PROBE_EXIT=\$?" 2>&1)
  echo "$out" | tee -a "$OUT"
  case "$out" in
    *"does not support remote access"*|*"Failed to connect"*) echo "  (ssh refused, trying another)" ;;
    *) break ;;
  esac
done
