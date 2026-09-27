#!/usr/bin/env bash
# dataVmc, second attempt. The first attempt (2026-09-24 00:42) reported RC=0 and
# "DATAVMC CHAIN DONE" while producing nothing new:
#  * make_dataVmc_condor_jobs.py REGENERATES dataVmc.submit from a template, which put
#    back LOCALIZE/SETUP=auto and a nonexistent AFS CONDA_TARBALL -> 81/81 jobs exit 2
#    ("CONDA_TARBALL set but missing"). The condor .out/.err were 0 bytes because the
#    wrapper tees everything into Plot/plots/logs_split/.
#  * the template had no JobBatchName, so the batch-name wait loop matched nothing,
#    returned at once, and the merge ran immediately -- on 192 partials dated
#    2026-05-28, i.e. the pre-FSR-fix production. The merged files it wrote carry
#    today's mtime and May's content.
#  * the generator also still listed the retired 2022 DYGto2LG_10to50/_50to100 slices.
# Fixes live in the generator now. This chain: parks the stale partials, preflights every
# input with the real resolver, smoketests 3 jobs on PRODUCT and MECHANISM, submits the
# rest, and waits on ClusterId -- which a regenerated submit file cannot break.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/Condor
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/variables_dataVmc
LS=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/logs_split

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u

wait_cluster() {   # $1 cluster, $2 label
  local n held last=-1
  while true; do
    n=$(condor_q "$1" -format "%d\n" ProcId 2>/dev/null | wc -l)
    held=$(condor_q "$1" -constraint 'JobStatus==5' -format "%d\n" ProcId 2>/dev/null | wc -l)
    if [ "$n" != "$last" ]; then echo "  [$(date '+%H:%M:%S')] $2 (cluster $1): $n in queue ($held held)"; last=$n; fi
    [ "$n" -eq 0 ] && return 0
    if [ "$held" -gt 0 ] && [ "$held" -eq "$n" ]; then
      condor_q "$1" -constraint 'JobStatus==5' -af HoldReason 2>/dev/null | sort | uniq -c | head -3
      return 1
    fi
    sleep 120
  done
}

# $1 = jobs file, $2 = epoch; verifies each job's partial + log are fresh and clean
check_jobs() {
  local jobs="$1" since="$2" bad=0 n=0 region tag stag samples suffix part log
  while read -r region tag stag samples; do
    [ -z "${region:-}" ] && continue
    n=$((n+1))
    case "$region" in SR) suffix=UL_SR;; CR) suffix=UL_CR;; mva) suffix=UL_mva;; esac
    part="$V/ALP_plot_run3_${suffix}_${tag}_part_${stag}.root"
    log="$LS/${tag}_${region}_${stag}.log"
    if [ ! -s "$part" ] || [ "$(stat -c %Y "$part")" -lt "$since" ]; then echo "  STALE/MISSING partial: $part"; bad=$((bad+1)); continue; fi
    if [ ! -s "$log" ] || [ "$(stat -c %Y "$log")" -lt "$since" ]; then echo "  STALE/MISSING log: $log"; bad=$((bad+1)); continue; fi
    grep -q "\[FAILED\]" "$log" && { echo "  FAILED in log: $log"; bad=$((bad+1)); continue; }
    grep -q "Fetch conda tarball" "$log" || { echo "  NOT LOCALIZED: $log"; bad=$((bad+1)); continue; }
    grep -q "Activate conda env" "$log" && { echo "  ACTIVATED FROM EOS: $log"; bad=$((bad+1)); continue; }
    grep -q "\[DONE\]" "$log" || { echo "  NO [DONE]: $log"; bad=$((bad+1)); continue; }
  done < "$jobs"
  echo "  checked $n jobs, problems $bad"
  [ "$n" -gt 0 ] && [ "$bad" -eq 0 ]
}

echo "=============== [1/6] park stale partials and contaminated merges ==============="
PARK="$V/stale_preFSRfix_20260924"
mkdir -p "$PARK"
nmv=$(ls "$V"/*_part_*.root 2>/dev/null | wc -l)
mv "$V"/*_part_*.root "$PARK"/ 2>/dev/null
for f in "$V"/ALP_plot_run3_UL_{SR,CR,mva}_sideband_rwgt.root; do [ -f "$f" ] && mv "$f" "$PARK"/; done
echo "  moved $nmv partials + merged files to $PARK"
echo "  partials left in place: $(ls "$V"/*_part_*.root 2>/dev/null | wc -l)"

echo
echo "=============== [2/6] preflight: every input resolves and exists ==============="
python "$D/preflight_datavmc.py" || { echo "PREFLIGHT FAILED -- not submitting"; exit 1; }

echo
echo "=============== [3/6] smoketest: 3 jobs ==============="
cd "$P"
{
  grep -E "^CR sideband_rwgt data "                  dataVmc_jobs.txt | head -1
  grep -E "^SR sideband_rwgt sig_m1 "                dataVmc_jobs.txt | head -1
  grep -E "^SR sideband_rwgt dyg_10to100_2022preEE " dataVmc_jobs.txt | head -1
} > dataVmc_smoke_jobs.txt
sed -e 's|dataVmc_jobs.txt|dataVmc_smoke_jobs.txt|' -e 's|"HZa_dataVmc"|"HZa_dataVmc_smoke"|' dataVmc.submit > dataVmc_smoke.submit
cat dataVmc_smoke_jobs.txt
T0=$(date +%s)
out=$(condor_submit dataVmc_smoke.submit 2>&1); echo "$out" | tail -1
CL=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out")
[ -z "$CL" ] && { echo "smoketest submission failed"; exit 1; }
wait_cluster "$CL" "smoke" || { echo "smoketest held -- stopping"; exit 1; }
check_jobs dataVmc_smoke_jobs.txt "$T0" || { echo "SMOKETEST FAILED (product or mechanism) -- not submitting the full set"; exit 1; }

echo
echo "=============== [4/6] full submission ==============="
T1=$(date +%s)
out=$(condor_submit dataVmc.submit 2>&1); echo "$out" | tail -1
CL=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out")
[ -z "$CL" ] && { echo "full submission failed"; exit 1; }
wait_cluster "$CL" "dataVmc" || echo "  (held jobs -- the per-job check below decides)"

echo
echo "=============== [5/6] every job: fresh partial, clean log, localized env ==============="
check_jobs dataVmc_jobs.txt "$T1" || { echo "JOB CHECK FAILED -- not merging"; exit 1; }

echo
echo "=============== [6/6] merge + draw ==============="
T2=$(date +%s)
bash 3_merge_dataVmc_condor.sh
echo "raw exit: $? (not a verdict)"
ok=1
for r in SR CR mva; do
  f="$V/ALP_plot_run3_UL_${r}_sideband_rwgt.root"
  if [ -s "$f" ] && [ "$(stat -c %Y "$f")" -ge "$T2" ]; then echo "  merged OK  $f"; else echo "  merged STALE/MISSING  $f"; ok=0; fi
done
for r in SR CR mva; do
  d=$V/plot_UL_${r}/sideband_rwgt
  n=$(find "$d" -name '*.pdf' -newermt "@$T2" 2>/dev/null | wc -l)
  echo "  fresh pdfs in plot_UL_${r}/sideband_rwgt: $n"
  [ "$n" -eq 0 ] && ok=0
done
echo "VERDICT: $([ $ok -eq 1 ] && echo PASSED || echo FAILED)"
echo "=============== DATAVMC2 DONE $(date '+%F %T') ==============="
[ $ok -eq 1 ]
