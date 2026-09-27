#!/usr/bin/env bash
#
# Re-run the bkg + data production at fpo=1, into a NEW staging directory.
#
# WHY fpo=1
#   With the copy-then-read change the jobs need a lot of memory at fpo=4:
#       MemoryUsage  median 7325 MB  p90 9766 MB   (20-job validation, 2026-09-15)
#   so RequestMemory had to go to 12000, and the CERN schedd silently widens that
#   to RequestCpus=4 (~3 GB per core). Measured slot supply at 21:55:
#       fits 4 cpu / 12 GB :   165 unclaimed slots
#       fits 2 cpu /  6 GB :   256
#       fits 1 cpu /  3 GB : 1392
#   We were getting 14 concurrent and ~35 chunks/h, i.e. ~50 h for the remaining
#   2000 chunks. At fpo=1 each job holds a quarter of the data (~1.8 GB), fits a
#   1-core slot, and the pool has 8.4x more of those.
#   Core-hours are a wash: 28k x 1 cpu x ~7 min vs 2k x 4 cpu x 27 min.
#
# WHY A NEW DIRECTORY, AND WHY THIS IS A FULL RE-RUN
#   fpo is pinned in analysis_manager.pkl and the CLI is ignored in its favour,
#   so changing it requires CLEAN_ANALYSIS_STATE=1, which renumbers every job.
#   Task.merge_outputs() merges self.outputs -- the pickle's job list -- not a
#   filesystem glob, so the 5044 chunks already produced under the fpo=4
#   numbering cannot be adopted by the new job structure. They are left in place
#   untouched as a fallback; nothing reads from both directories at once.
#
# USAGE
#   bash submit_fpo1_production.sh             # dry run
#   DRY_RUN=0 bash submit_fpo1_production.sh   # run
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
STAGES="${STAGES:-bkg data}"
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
NEW=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana

export CONDA_PREFIX="$ENVDIR"
export PATH="$ENVDIR/bin:$PATH"
export PYTHONPATH="$REPO"
export X509_USER_PROXY="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"
export HZA_BMODE_PACK="${HZA_BMODE_PACK:-root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz}"
export HZA_STAGE_INPUTS=1
# 3000 keeps RequestCpus at 1 -- the whole point of the change. Anything above
# ~3 GB gets widened by the schedd and lands back in the small slot pool.
export HIGGSDNA_CONDOR_REQ_MEMORY=3000

left=$(voms-proxy-info -timeleft 2>/dev/null || echo 0)
echo "proxy left: ${left}s"
[ "${left:-0}" -lt 86400 ] && { echo "[ERROR] proxy under 24 h"; exit 1; }
# The repo copy is what ships with each job; /tmp being fresh is not enough.
cp -p "$X509_USER_PROXY" "$REPO/x509up_u$(id -u)" 2>/dev/null
"$ENVDIR/bin/openssl" x509 -noout -enddate -in "$REPO/x509up_u$(id -u)"

for st in $STAGES; do
    case "$st" in
        bkg)  script=run_ana_bkgmc.sh; out="$NEW/Bkg_MC" ;;
        data) script=run_ana_data.sh;  out="$NEW/Data" ;;
    esac
    echo "--- $st -> $out (fpo=1, mem 3000, workday)"
    if [ "$DRY_RUN" = "1" ]; then
        echo "    would run: OUTDIR=$out FPO=1 CLEAN_ANALYSIS_STATE=1 CONDOR_REQ_MEMORY=3000 bash $REPO/scripts/$script"
        continue
    fi
    mkdir -p "$out"
    ( cd "$REPO" && setsid nohup env \
        "CONDA_PREFIX=$ENVDIR" "PATH=$ENVDIR/bin:$PATH" "PYTHONPATH=$REPO" \
        "HZA_BMODE_PACK=$HZA_BMODE_PACK" "HZA_STAGE_INPUTS=1" \
        "CONDOR_REQ_MEMORY=3000" "OUTDIR=$out" "FPO=1" \
        "CLEAN_ANALYSIS_STATE=1" "UNRETIRE_JOBS=1" "DRY_RUN=0" \
        bash "scripts/$script" > "$D/logs_fsrfix/${st}_fpo1.log" 2>&1 < /dev/null & )
    echo "    launched"
done
