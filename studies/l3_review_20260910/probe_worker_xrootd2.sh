#!/usr/bin/env bash
#
# Probe v2: find out WHERE the worker-side xrootd open stalls.
#
# Established by v1 (logs_fsrfix/probe_worker.log):
#   PROBE proxy=/pool/condor/dir_113164/x509up_u175325 exists=True
#   PROBE_EXIT=124
# i.e. with a valid proxy and the B-mode env, uproot.open on a file that lxplus
# opens in 2.9 s (989393 entries) hung for the full 130 s and never returned --
# no exception, no error string. v2 asks the XRootD client itself where it stops.
#
# Each step prints a STEP line BEFORE it runs, so a hang is attributable to a
# specific step rather than showing up as silence. Nothing is filtered.
#
# USAGE
#   bash probe_worker_xrootd2.sh             # dry run
#   DRY_RUN=0 nohup bash probe_worker_xrootd2.sh >/dev/null 2>&1 &
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
OUT="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/probe_worker2.log"
F="/store/data/Run2023C/EGamma1/NANOAOD/22Sep2023_v4-v1/2540000/8ba7b401-41ef-4f68-b26d-34a3f633f7d6.root"

read -r -d '' REMOTE <<REOF
set -u
E=\$_CONDOR_SCRATCH_DIR/hza_ana_env
export X509_USER_PROXY=\$_CONDOR_SCRATCH_DIR/x509up_u175325
echo "STEP0 host=\$(hostname) proxy_exists=\$([ -f \$X509_USER_PROXY ] && echo yes || echo NO)"
echo "STEP1 CA dir /etc/grid-security/certificates: \$(ls /etc/grid-security/certificates 2>/dev/null | wc -l) files"
echo "STEP1b X509_CERT_DIR=\${X509_CERT_DIR:-unset}  XRD_ env: \$(env | grep -c '^XRD_')"
echo "STEP2 proxy validity:"
timeout 30 \$E/bin/openssl x509 -noout -dates -in \$X509_USER_PROXY 2>&1 | head -2
echo "STEP3 xrdfs stat via infn redirector (60s budget):"
timeout 60 \$E/bin/xrdfs xrootd-cms.infn.it stat "$F" 2>&1 | head -6
echo "STEP3 rc=\$?"
echo "STEP4 xrdfs locate (60s budget):"
timeout 60 \$E/bin/xrdfs xrootd-cms.infn.it locate -d "$F" 2>&1 | head -6
echo "STEP4 rc=\$?"
echo "STEP5 xrdcp first 10 MB to scratch (90s budget):"
timeout 90 \$E/bin/xrdcp -f --nopbar "root://xrootd-cms.infn.it/$F" \$_CONDOR_SCRATCH_DIR/probe.root 2>&1 | tail -4
echo "STEP5 rc=\$? size=\$(stat -c%s \$_CONDOR_SCRATCH_DIR/probe.root 2>/dev/null || echo none)"
echo "STEP6 uproot.open with XRD_LOGLEVEL=Debug (100s budget), last 25 log lines:"
cat > \$_CONDOR_SCRATCH_DIR/p2.py <<'PYEOF'
import time, uproot
t=time.time()
try:
    f=uproot.open("root://xrootd-cms.infn.it/$F", timeout=40)
    print("UPROOT OK %.1fs n=%d"%(time.time()-t, f["Events"].num_entries))
except Exception as e:
    print("UPROOT FAIL %.1fs %s %s"%(time.time()-t, type(e).__name__, str(e)[:120]))
PYEOF
XRD_LOGLEVEL=Debug timeout 100 \$E/bin/python -s \$_CONDOR_SCRATCH_DIR/p2.py 2>&1 | tail -25
echo "STEP6 rc=\$?"
echo "PROBE2 DONE"
REOF
B64=$(printf '%s' "$REMOTE" | base64 -w0)

cands=$(condor_q -run -af ClusterId ProcId JobBatchName 2>/dev/null | awk '$3 ~ /^(DY|Data_20)/ {print $1"."$2}' | shuf | head -8)
echo "candidates: $(wc -w <<< "$cands")"
[ -z "${cands// }" ] && { echo "no running production job to probe"; exit 1; }
if [ "$DRY_RUN" = "1" ]; then echo "would probe one of: $cands"; exit 0; fi

for id in $cands; do
    echo "=== $(date '+%F %T') probe2 on $id ===" | tee -a "$OUT"
    out=$(timeout 480 condor_ssh_to_job "$id" "echo $B64 | base64 -d > \$_CONDOR_SCRATCH_DIR/p2.sh; bash \$_CONDOR_SCRATCH_DIR/p2.sh" 2>&1)
    echo "$out" | tee -a "$OUT"
    case "$out" in
        *"does not support remote access"*|*"Failed to connect"*) echo "  (ssh refused, next)" | tee -a "$OUT" ;;
        *) break ;;
    esac
done
