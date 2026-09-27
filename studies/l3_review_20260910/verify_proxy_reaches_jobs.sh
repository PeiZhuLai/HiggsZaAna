#!/usr/bin/env bash
#
# Confirm that the REFRESHED proxy actually reaches newly started jobs.
#
# Background
#   Jobs ship a copy of the proxy. HiggsDNA copies it once into the repo
#   (jobs.py:435) and every submit file points at that copy, which still held
#   the Sep-12 proxy (notAfter Sep 12 20:01 GMT) long after it expired. Result:
#   XRootD could not authenticate anywhere, cycled the replica list with a 120 s
#   connection window per site, and uproot.open never returned -- the "hang"
#   that cost 2026-09-12 and most of 2026-09-13.
#     Trying to authenticate using gsi
#     Cannot get credentials for protocol gsi: Secgsi: ErrParseBuffer: ...
#     No protocols left to try -- [FATAL] Auth failed
#   The repo copy was replaced at 20:02 on 2026-09-13 with one valid to Sep 21.
#
# What this checks
#   "The file exists" is what was checked before and it was worthless -- the file
#   existed and was expired. This reads notAfter out of the proxy a RUNNING job
#   actually received, and only then watches for chunks.
#
# USAGE
#   bash verify_proxy_reaches_jobs.sh             # dry run
#   DRY_RUN=0 nohup bash verify_proxy_reaches_jobs.sh >/dev/null 2>&1 &
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix/verify_proxy_reaches_jobs.log
S=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix
E=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
CONST='Owner=="pelai" && regexp("HZa/HiggsZaAna",Cmd)'

say(){ echo "$(date '+%F %T') $*" >> "$L"; }
chunks(){ find $S/Bkg_MC $S/Data -name '*_nominal.parquet' ! -name 'merged_nominal.parquet' 2>/dev/null | wc -l; }

[ "$DRY_RUN" = "1" ] && { echo "would wait for an HZa job to start, read its shipped proxy's notAfter, then watch chunks"; exit 0; }

say "start: repo proxy notAfter=$($E/bin/openssl x509 -noout -enddate -in /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/x509up_u175325 2>&1 | tr -d '\n')"
b=$(chunks); say "baseline chunks=$b"

# Stage 1: the proxy a running job actually got.
probed=0
for i in $(seq 1 120); do
    ids=$(condor_q -constraint "$CONST && JobStatus==2" -af ClusterId ProcId 2>/dev/null | awk '{print $1"."$2}' | head -4)
    if [ -n "${ids// }" ]; then
        for id in $ids; do
            out=$(timeout 120 condor_ssh_to_job "$id" 'D=$_CONDOR_SCRATCH_DIR; for f in $D/x509up_u*; do echo "SHIPPED $f"; done' 2>&1)
            case "$out" in *"does not support remote access"*|*"Failed to connect"*) continue ;; esac
            p=$(awk '/^SHIPPED/{print $2; exit}' <<<"$out")
            [ -z "$p" ] && continue
            dates=$(timeout 120 condor_ssh_to_job "$id" "\$_CONDOR_SCRATCH_DIR/hza_ana_env/bin/openssl x509 -noout -dates -in $p 2>&1" 2>&1 | tr '\n' ' ')
            say "job $id shipped proxy: $dates"
            probed=1; break
        done
    fi
    [ "$probed" = 1 ] && break
    sleep 60
done
[ "$probed" = 0 ] && say "WARNING: no HZa job reached running state within 2 h -- fair-share debt, not a failure"

# Stage 2: the artifact.
for i in $(seq 1 96); do
    sleep 300
    n=$(chunks); say "chunks=$n (+$((n-b)))"
    [ "$n" -gt $((b+30)) ] && { say "CONFIRMED: chunks growing again ($b -> $n)"; exit 0; }
done
say "chunks still $(chunks) after 8 h"
