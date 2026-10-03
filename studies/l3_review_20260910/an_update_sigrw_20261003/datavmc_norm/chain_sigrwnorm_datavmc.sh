#!/usr/bin/env bash
# Regenerate ONLY the signal partials of dataVmc final tag sideband_rwgt with the signal
# normalization N (and the --optimizeBranches fix) in Plot/scripts/1_prepare_dataVmc.py, then
# merge + draw, then rerun the working-point optimization the way chain_dyveto_flashgg.sh [2/9]
# does -- WITHOUT writing Plot/output/MVAcut_points_run3.json (collect and the R scan are run on
# copies that write into this directory).
#
# Mechanism copied from chain_dyveto_datavmc.sh: generate jobs, preflight, 1-job smoketest checked
# on product AND mechanism, wait on ClusterId, per-job freshness/log check, merge.
# Background and data partials of sideband_rwgt (2026-09-29/30) are kept as they are.
#
#   systemd-run --user --unit=hza-sigrwnorm-dvmc --collect bash -c "bash <this> > <dir>/chain.log 2>&1"
#   START_STAGE=N skips stages < N.
set -uo pipefail
FT=sideband_rwgt
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/an_update_sigrw_20261003/datavmc_norm
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
P=$PLOT/Condor
V=$PLOT/plots/variables_dataVmc
LS=$PLOT/plots/logs_split
OPT=$PLOT/plots/optimize_run3UL
OPTBAK=$PLOT/plots_variants/optimize_run3UL_stale_preSigRwNorm_20261003   # outside Plot/plots (sync_figures.sh)
PARK=$V/stale_preSigRwNorm_${FT}_20261003
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
MASS14="1 2 3 4 5 6 7 8 9 10 15 20 25 30"
SIGJOBS=$P/dataVmc_sigrwnorm_jobs.txt
SIGSUB=$P/dataVmc_sigrwnorm.submit
START=${START_STAGE:-0}
run() { [ "$START" -le "$1" ]; }
die() { echo "STOPPING: $*"; echo "VERDICT: FAILED"; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
echo "[env] python=$(which python) host=$(hostname) start=$(date '+%F %T')"

wait_cluster() {   # $1 cluster, $2 label
  local n held last=-1
  while true; do
    n=$(condor_q "$1" -format "%d\n" ProcId 2>/dev/null | wc -l)
    held=$(condor_q "$1" -constraint 'JobStatus==5' -format "%d\n" ProcId 2>/dev/null | wc -l)
    if [ "$n" != "$last" ]; then echo "  [$(date '+%H:%M:%S')] $2 (cluster $1): $n in queue ($held held)"; last=$n; fi
    [ "$n" -eq 0 ] && return 0
    if [ "$held" -gt 0 ]; then
      # memory holds: 8000 -> 16000 -> 20000 (never condor_rm)
      condor_q "$1" -constraint 'JobStatus==5' -af ClusterId ProcId RequestMemory HoldReason 2>/dev/null | while read -r c p mem reason; do
        echo "  HELD $c.$p mem=$mem reason: $reason"
        if echo "$reason" | grep -qi "memory"; then
          if [ "$mem" -lt 8000 ]; then new=8000; elif [ "$mem" -lt 16000 ]; then new=16000; elif [ "$mem" -lt 20000 ]; then new=20000; else new=0; fi
          if [ "$new" -gt 0 ]; then
            echo "    condor_qedit $c.$p RequestMemory $new && condor_release $c.$p"
            condor_qedit "$c.$p" RequestMemory "$new" && condor_release "$c.$p"
          else
            echo "    already at 20000 MB -- leaving it held"
          fi
        fi
      done
      held=$(condor_q "$1" -constraint 'JobStatus==5' -format "%d\n" ProcId 2>/dev/null | wc -l)
      [ "$held" -gt 0 ] && [ "$held" -eq "$n" ] && { condor_q "$1" -constraint 'JobStatus==5' -af HoldReason | sort | uniq -c | head -3; return 1; }
    fi
    sleep 120
  done
}

# $1 = jobs file, $2 = epoch; each job's partial + log must be fresh and clean, and the log must
# show the conda localization AND the new signal-normalization code path
check_jobs() {
  local jobs="$1" since="$2" bad=0 n=0 region tag stag samples suffix part log
  while read -r region tag stag samples; do
    [ -z "${region:-}" ] && continue
    n=$((n+1))
    case "$region" in SR) suffix=UL_SR;; CR) suffix=UL_CR;; mva) suffix=UL_mva;; esac
    part="$V/ALP_plot_run3_${suffix}_${tag}_part_${stag}.root"
    log="$LS/${tag}_${region}_${stag}.log"
    if ! fresh "$part" "$since"; then echo "  STALE/MISSING partial: $part"; bad=$((bad+1)); continue; fi
    if ! fresh "$log" "$since"; then echo "  STALE/MISSING log: $log"; bad=$((bad+1)); continue; fi
    grep -q "\[FAILED\]" "$log" && { echo "  FAILED in log: $log"; bad=$((bad+1)); continue; }
    grep -q "Fetch conda tarball" "$log" || { echo "  NOT LOCALIZED: $log"; bad=$((bad+1)); continue; }
    grep -q "Activate conda env" "$log" && { echo "  ACTIVATED FROM EOS: $log"; bad=$((bad+1)); continue; }
    grep -q "\[DONE\]" "$log" || { echo "  NO [DONE]: $log"; bad=$((bad+1)); continue; }
    grep -q "\[SignalRwNorm\] Loaded 140 " "$log" || { echo "  NORMS NOT LOADED: $log"; bad=$((bad+1)); continue; }
    grep -q "sideband reweight variables read from enabled branches" "$log" || { echo "  BRANCH CHECK MISSING: $log"; bad=$((bad+1)); continue; }
    nera=$(grep -c "^\[SignalRwNorm\] sample=$samples era=" "$log")
    [ "$nera" -eq 5 ] || { echo "  only $nera/5 eras normalized: $log"; bad=$((bad+1)); continue; }
  done < "$jobs"
  echo "  checked $n jobs, problems $bad"
  [ "$n" -gt 0 ] && [ "$bad" -eq 0 ]
}

if run 0; then
echo "=============== [0/7] generate ${FT} jobs, keep the signal lines ==============="
python "$P/make_dataVmc_condor_jobs.py" --condor-dir "$P" --final-tags ${FT} || die "job generation failed"
[ "$(wc -l < $P/dataVmc_jobs.txt)" -eq 75 ] || die "expected 75 ${FT} jobs"
awk -v ft="$FT" '$2!=ft' "$P/dataVmc_jobs.txt" | grep -q . && die "jobs file has non-${FT} lines"
grep -q 'JobBatchName = "HZa_dataVmc"' "$P/dataVmc.submit" || die "submit file lacks JobBatchName -- generator regressed"
awk '$3 ~ /^sig_/' "$P/dataVmc_jobs.txt" > "$SIGJOBS"
[ "$(wc -l < $SIGJOBS)" -eq 42 ] || die "expected 42 signal jobs (14 masses x SR/CR/mva)"
sed -e "s|queue region_key, final_tag, sample_tag, samples from .*|queue region_key, final_tag, sample_tag, samples from $SIGJOBS|" \
    -e 's|"HZa_dataVmc"|"HZa_dataVmc_sigrwnorm"|' "$P/dataVmc.submit" > "$SIGSUB"
grep -q "from $SIGJOBS" "$SIGSUB" || die "signal submit file not rewritten"
echo "  signal jobs: $(wc -l < $SIGJOBS) -> $SIGJOBS"
fi
export DATA_VMC_MAKE_JOBS_ARGS="--final-tags ${FT}"

if run 1; then
echo; echo "=============== [1/7] park the old ${FT} signal partials and the merged files ==============="
mkdir -p "$PARK"
nmv=0
for f in "$V"/ALP_plot_run3_UL_{SR,CR,mva}_${FT}_part_sig_m*.root; do [ -f "$f" ] && { mv "$f" "$PARK"/; nmv=$((nmv+1)); }; done
for f in "$V"/ALP_plot_run3_UL_{SR,CR,mva}_${FT}.root "$V"/ALP_plot_run3_UL_{SR,CR,mva}_${FT}_plots.root; do [ -f "$f" ] && cp -p "$f" "$PARK"/; done
echo "  moved $nmv signal partials, copied the merged files, into $PARK"
nbd=$(ls "$V"/ALP_plot_run3_UL_*_${FT}_part_*.root 2>/dev/null | wc -l)
echo "  ${FT} background/data partials left in place: $nbd (expect 33)"
[ "$nbd" -eq 33 ] || die "background/data partials are not all there"
fi

if run 2; then
echo; echo "=============== [2/7] preflight: every signal input resolves and exists ==============="
python "$D/preflight_datavmc.py" "$SIGJOBS" || die "PREFLIGHT FAILED"
fi

if run 3; then
echo; echo "=============== [3/7] smoketest: 1 job ==============="
cd "$P"
grep -E "^SR ${FT} sig_m5 " "$SIGJOBS" > dataVmc_sigrwnorm_smoke_jobs.txt
sed -e "s|from $SIGJOBS|from $P/dataVmc_sigrwnorm_smoke_jobs.txt|" -e 's|"HZa_dataVmc_sigrwnorm"|"HZa_dataVmc_sigrwnorm_smoke"|' "$SIGSUB" > dataVmc_sigrwnorm_smoke.submit
cat dataVmc_sigrwnorm_smoke_jobs.txt
T0=$(date +%s); echo "$T0" > $W/T0_smoke
out=$(condor_submit dataVmc_sigrwnorm_smoke.submit 2>&1); echo "$out" | tail -1
CL=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out")
[ -z "$CL" ] && die "smoketest submission failed"
echo "SMOKE_CLUSTER=$CL"
wait_cluster "$CL" "smoke" || die "smoketest held"
check_jobs dataVmc_sigrwnorm_smoke_jobs.txt "$T0" || die "SMOKETEST FAILED (product or mechanism)"
grep -A8 "filled events per" "$LS/${FT}_SR_sig_m5.log" | head -14
fi

if run 4; then
echo; echo "=============== [4/7] remaining signal jobs ==============="
cd "$P"
grep -v -E "^SR ${FT} sig_m5 " "$SIGJOBS" > dataVmc_sigrwnorm_rest_jobs.txt
[ "$(wc -l < dataVmc_sigrwnorm_rest_jobs.txt)" -eq 41 ] || die "expected 41 remaining jobs"
sed -e "s|from $SIGJOBS|from $P/dataVmc_sigrwnorm_rest_jobs.txt|" "$SIGSUB" > dataVmc_sigrwnorm_rest.submit
T1=$(date +%s); echo "$T1" > $W/T1_rest
out=$(condor_submit dataVmc_sigrwnorm_rest.submit 2>&1); echo "$out" | tail -1
CL=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out")
[ -z "$CL" ] && die "submission failed"
echo "REST_CLUSTER=$CL"
wait_cluster "$CL" "signal" || echo "  (held jobs -- the per-job check below decides)"
fi

if run 5; then
echo; echo "=============== [5/7] every signal job: fresh partial, clean log, localized env, N applied ==============="
T0=$(cat $W/T0_smoke)
check_jobs "$SIGJOBS" "$T0" || die "JOB CHECK FAILED -- not merging"
fi

if run 6; then
echo; echo "=============== [6/7] merge + draw ==============="
T2=$(date +%s); echo "$T2" > $W/T2_merge
cd "$P"
bash 3_merge_dataVmc_condor.sh > $W/merge_draw.log 2>&1
echo "  raw exit: $? (not a verdict; log $W/merge_draw.log)"
ok=1
for r in SR CR mva; do
  f="$V/ALP_plot_run3_UL_${r}_${FT}.root"
  if fresh "$f" "$T2"; then echo "  merged OK  $f"; else echo "  merged STALE/MISSING  $f"; ok=0; fi
done
for r in SR CR mva; do
  n=$(find "$V/plot_UL_${r}/${FT}" -name '*.pdf' -newermt "@$T2" 2>/dev/null | wc -l)
  echo "  fresh pdfs in plot_UL_${r}/${FT}: $n"; [ "$n" -eq 0 ] && ok=0
done
[ $ok -eq 1 ] || die "merge/draw incomplete"
fi

if run 7; then
echo; echo "=============== [7/7] working points (ALP_Optimization -c 2, -c 1; collect + R scan into $W) ==============="
T=$(date +%s)
if [ ! -d "$OPTBAK" ]; then cp -a "$OPT" "$OPTBAK" || die "backup of optimize_run3UL failed"; fi
echo "  backup: $OPTBAK ($(ls $OPTBAK | wc -l) files)"
mv $OPT/nCat*_all_M*.json "$OPTBAK"/ 2>/dev/null   # chain_dyveto_flashgg.sh [2/9]: stale copies out of the way
cd $PLOT
export PYTHONPATH="${PYTHONPATH:-}:$PLOT/lib"
$ALP_PY scripts/ALP_Optimization.py -y run3 -o $OPT --region 2 -p --sigVSscore -s --doOpt -c 2 --inputTag ${FT} > $W/wp_opt_c2.log 2>&1
$ALP_PY scripts/ALP_Optimization.py -y run3 -o $OPT --region 2 -p --sigVSscore -s --doOpt -c 1 --inputTag ${FT} > $W/wp_opt_c1.log 2>&1
bad=0; for m in $MASS14; do fresh $OPT/nCat1_all_M$m.json $T || { echo "  missing/stale nCat1_all_M$m.json"; bad=1; }; fresh $OPT/nCat2_all_M$m.json $T || { echo "  missing/stale nCat2_all_M$m.json"; bad=1; }; done
[ $bad -eq 0 ] || die "ALP_Optimization did not produce every mass point (logs: $W/wp_opt_c*.log)"
# collect_MVAcut_points_run3.py and scan_score_R_significance.py write the analysis MVA cut
# configuration; run copies that write here instead
sed "s|^output_json = .*|output_json = \"$W/MVAcut_points_run3_candidate.json\"|" scripts/collect_MVAcut_points_run3.py > $W/collect_MVAcut_points_run3_copy.py
grep -q "MVAcut_points_run3_candidate.json" $W/collect_MVAcut_points_run3_copy.py || die "collect copy not redirected"
$ALP_PY $W/collect_MVAcut_points_run3_copy.py > $W/wp_collect.log 2>&1
fresh $W/MVAcut_points_run3_candidate.json $T || die "collect did not write the candidate JSON"
sed "s|^OUTDIR   = .*|OUTDIR   = \"$W/scan_score_R\"|" scripts/scan_score_R_significance.py > $W/scan_score_R_significance_copy.py
grep -q "^OUTDIR   = \"$W/scan_score_R\"" $W/scan_score_R_significance_copy.py || die "scan copy not redirected"
$ALP_PY $W/scan_score_R_significance_copy.py > $W/wp_scan.log 2>&1 || die "R scan failed ($W/wp_scan.log)"
echo "  MVAcut_points_run3.json untouched: $(stat -c '%y' $PLOT/output/MVAcut_points_run3.json)"
fi
echo
echo "VERDICT: PASSED"
echo "=============== SIGRWNORM DATAVMC CHAIN DONE $(date '+%F %T') ==============="
