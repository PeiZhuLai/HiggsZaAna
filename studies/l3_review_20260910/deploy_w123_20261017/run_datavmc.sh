#!/usr/bin/env bash
# Step E: data/simulation (dataVmc) and BDT-input plots for the low_relE2_w123 low-mass BDT.
#   1  wait until the scored files dataVmc reads are patched: the 273 of the first batch (data, DY, mA_M1-3 + variations)
#      and the 55 nominal signal files of m_a 4-30 (their MVA_Score_mA_M1..3 are filled in the m_a 1-3 histograms)
#   2  backup: Plot/plots/variables_dataVmc (tar on EOS), optimize_run3UL*, train_vars_sigVbkg, plots_MVA/run3_lowmass
#   3  condor: 1_prepare_dataVmc.py for both tags (nominal, sideband_rwgt) x SR/CR/mva = 150 jobs (Plot/Condor,
#      logs Plot/Condor/logs); jobs whose output was not rewritten are resubmitted (at most 3 rounds)
#   4  after the production chain has finished (local process limit): 3_merge_dataVmc_condor.sh (merge, draw,
#      optimization region 1 -> optimize_run3UL_<tag>) and ALP_Optimization region 2, -c 2 and -c 1 ->
#      optimize_run3UL (the AN directory, as 1_runPlot.sh). collect_MVAcut_points_run3.py and
#      scan_score_R_significance.py --write-json are NOT run: they would overwrite the m_a 1-3 working points.
#      The working-point JSON must be byte-identical before and after.
#   5  BDT-input plots: plot_training_variables.py pho2_relEerr; low-mass diagnostics (plot_mva_diagnostics.py +
#      plot_overtrain_lowmass.py) into plots_MVA/run3_lowmass
#   6  checks: every dataVmc PDF of the CR/mva regions has the ratio pad ("Data / SM", pdftotext), the new
#      19_pho2_relEerr plots exist for both tags
# Markers: logs/datavmc.DONE | logs/datavmc.FAILED; heartbeat logs/datavmc_heartbeat.txt
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
CD=$P/Condor
MVA=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA
WPJ=$P/output/MVAcut_points_run3.json
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
EB=root://eosuser.cern.ch//eos/user/p/pelai/HZa/backup_w123_20261017
L=$D/logs; BK=$D/backup_prod
mkdir -p $L $BK
hb() { echo "$(date '+%F %T') $*" | tee -a $L/datavmc_heartbeat.txt; }
die() { hb "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/datavmc.FAILED; exit 1; }
rm -f $L/datavmc.DONE $L/datavmc.FAILED
START=${START_STAGE:-1}
run() { [ "$START" -le "$1" ]; }

if run 1; then
hb "[1] wait for the scored files dataVmc reads"
while true; do
  n=0
  while read -r s e; do [ -e $D/state/${s}__${e}.ok ] && n=$((n + 1)); done < <(cat $L/patch_first.txt; grep -E '^mA_M([4-9]|10|15|20|25|30) ' $L/patch_all.txt)
  [ $n -eq 328 ] && break
  grep -q " FAIL " $L/patch_progress.txt && die "score patch has failed files"
  hb "    $n/328 patched"; sleep 600
done
hb "    328/328 patched"
fi

if run 2; then
hb "[2] backups"
xrdfs eosuser.cern.ch mkdir -p /eos/user/p/pelai/HZa/backup_w123_20261017
if ! xrdfs eosuser.cern.ch stat /eos/user/p/pelai/HZa/backup_w123_20261017/variables_dataVmc_before.tar >/dev/null 2>&1; then
  tar -C $P/plots -cf - variables_dataVmc | xrdcp -s - $EB/variables_dataVmc_before.tar || die "variables_dataVmc backup"
fi
[ -s $BK/plots_before_datavmc.tar ] || tar -C $P/plots -cf $BK/plots_before_datavmc.tar optimize_run3UL optimize_run3UL_nominal optimize_run3UL_sideband_rwgt train_vars_sigVbkg || die "plot backup"
[ -s $BK/plots_MVA_lowmass_before.tar ] || tar -C $MVA/plots_MVA -cf $BK/plots_MVA_lowmass_before.tar run3_lowmass || die "plots_MVA backup"
until grep -q "\[1\] working points" $L/chain_heartbeat.txt && grep -q "\[2\] S1" $L/chain_heartbeat.txt; do
  [ -e $L/chain.FAILED ] && die "production chain failed before the working points were set"; sleep 120; done
$ALP_PY -c "
import json, sys
r = {int(e['mA']): float(e['MVAcut']) for e in json.load(open('$WPJ'))['results']}
sys.exit(0 if [r[1], r[2], r[3]] == [0.98115, 0.98837, 0.99226] else 1)" || die "working-point JSON does not hold the w123 cuts"
md5sum $WPJ > $L/wpjson.md5
fi

if run 3; then
hb "[3] condor dataVmc jobs (nominal + sideband_rwgt)"
( cd $CD && NO_SUBMIT=1 SKIP_ENV_PACK=1 bash 1_submit_dataVmc_condor.sh --final-tags nominal,sideband_rwgt ) > $L/datavmc_make.log 2>&1 || die "job generation (datavmc_make.log)"
nj=$(wc -l < $CD/dataVmc_jobs.txt); [ $nj -eq 150 ] || die "dataVmc_jobs.txt has $nj jobs, expected 150"
T0=$(date +%s); echo $T0 > $L/datavmc_T0
round=0
sub=$CD/dataVmc.submit
while true; do
  out=$(cd $CD && condor_submit $sub 2>&1); echo "$out" >> $L/datavmc_submit.log
  CL=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out"); [ -n "$CL" ] || die "condor_submit failed (datavmc_submit.log)"
  hb "    round $round: cluster $CL ($(grep -oP '^\d+(?= job)' <<< "$out") jobs)"
  while true; do
    q=$(condor_q $CL -af JobStatus 2>/dev/null) || { sleep 300; continue; }
    idle=$(grep -c '^1$' <<< "$q"); runn=$(grep -c '^2$' <<< "$q"); held=$(grep -c '^5$' <<< "$q")
    while read -r id mem; do
      [ -z "${id:-}" ] && continue
      if [ "$mem" -lt 8000 ]; then new=8000; elif [ "$mem" -lt 16000 ]; then new=16000; elif [ "$mem" -lt 20000 ]; then new=20000; else continue; fi
      hb "    held $id (memory) $mem -> $new"; condor_qedit $id RequestMemory $new > /dev/null; condor_release $id > /dev/null
    done < <(condor_q $CL -constraint 'JobStatus==5 && regexp("emory", HoldReason)' -af:j RequestMemory 2>/dev/null)
    hb "    cluster $CL idle=$idle running=$runn held=$held"
    [ $((idle + runn + held)) -eq 0 ] && break
    [ $((idle + runn)) -eq 0 ] && [ $held -gt 0 ] && { condor_q $CL -hold -af HoldReason | sort | uniq -c > $L/datavmc_held_$CL.txt; condor_rm $CL > /dev/null; break; }
    sleep 1200
  done
  # a job output counts only if it was rewritten after T0 (the previous production's partial files are still on disk,
  # so 2_resubmit_dataVmc_condor.sh, which checks existence, cannot see a failed job)
  ( cd $CD && awk '{print $1, $2, $3}' dataVmc_jobs.txt | while read -r r t s; do
      case $r in SR) sfx=UL_SR;; CR) sfx=UL_CR;; mva) sfx=UL_mva;; esac
      f=$P/plots/variables_dataVmc/ALP_plot_run3_${sfx}_${t}_part_${s}.root
      [ -s $f ] && [ $(stat -c %Y $f) -ge $T0 ] || grep -E "^$r $t $s " dataVmc_jobs.txt; done ) > $CD/dataVmc_w123_redo_jobs.txt
  stale=$(grep -c . $CD/dataVmc_w123_redo_jobs.txt)
  hb "    after round $round: $stale job outputs not refreshed"
  [ "$stale" -eq 0 ] && break
  round=$((round + 1)); [ $round -gt 3 ] && die "dataVmc outputs still missing after 3 resubmissions"
  sed "s#/dataVmc_jobs.txt#/dataVmc_w123_redo_jobs.txt#" $CD/dataVmc.submit > $CD/dataVmc_w123_redo.submit
  grep -q dataVmc_w123_redo_jobs.txt $CD/dataVmc_w123_redo.submit || die "redo submit file"
  sub=$CD/dataVmc_w123_redo.submit
done
hb "    150/150 partial outputs refreshed"
fi

if run 4; then
hb "[4] wait for the production chain, then merge / draw / optimization"
while [ ! -e $L/chain.DONE ] && [ ! -e $L/chain.FAILED ]; do sleep 600; done
( cd $CD && RUN_MERGE_PLOTS=1 RUN_DATAVMC_PLOTS=1 RUN_OPTIMIZATION=1 bash 3_merge_dataVmc_condor.sh ) > $L/datavmc_merge.log 2>&1 || die "3_merge_dataVmc_condor.sh failed (datavmc_merge.log)"
export PYTHONPATH="$P/lib:$MVA/scripts"
for c in 2 1; do
  ( cd $P && $ALP_PY scripts/ALP_Optimization.py -y run3 -o $P/plots/optimize_run3UL --region 2 -p --sigVSscore -s --doOpt -c $c --inputTag sideband_rwgt ) > $L/datavmc_opt_c$c.log 2>&1
  rc=$?; [ $rc -eq 0 ] || [ $rc -eq 139 ] || die "ALP_Optimization -c $c rc=$rc"
done
md5sum -c $L/wpjson.md5 > /dev/null || die "the working-point JSON was modified by the dataVmc / optimization step"
hb "    merge, draw and optimization done; working-point JSON unchanged"
fi

if run 5; then
hb "[5] BDT-input plots"
( cd $P/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY plot_training_variables.py pho2_relEerr ) > $L/train_vars_w123.log 2>&1 || die "plot_training_variables.py"
( cd $MVA/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY plot_mva_diagnostics.py --model $MVA/using/model_Za_BDT_lowmass_run3.pkl \
    --meta $MVA/using/model_Za_BDT_lowmass_run3.meta.json --outdir $MVA/plots_MVA/run3_lowmass && \
  env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY plot_overtrain_lowmass.py ) > $L/mva_diag_lowmass.log 2>&1 || die "low-mass diagnostics"
fi

if run 6; then
hb "[6] checks"
T0=$(cat $L/datavmc_T0)
bad=0; nchk=0
for t in nominal sideband_rwgt; do
  for f in $P/plots/variables_dataVmc/plot_UL_CR/$t/*_log.pdf $P/plots/variables_dataVmc/plot_UL_mva/$t/mvaVal_larger_M*_log.pdf; do
    [ -e "$f" ] || continue
    nchk=$((nchk + 1))
    [ $(stat -c %Y "$f") -ge $T0 ] || { echo "stale $f"; bad=$((bad + 1)); continue; }
    case $f in *sigma_*) continue;; esac
    pdftotext "$f" - 2>/dev/null | grep -q "Data / SM" || { echo "no ratio pad $f"; bad=$((bad + 1)); }
  done
  [ -s $P/plots/variables_dataVmc/plot_UL_CR/$t/19_pho2_relEerr_log.pdf ] || { echo "missing 19_pho2_relEerr_log.pdf ($t)"; bad=$((bad + 1)); }
done > $L/datavmc_pdf_check.txt
bad=$(grep -c . $L/datavmc_pdf_check.txt)
hb "    PDF check: $nchk files, $bad problems (datavmc_pdf_check.txt)"
[ $bad -eq 0 ] || die "dataVmc PDF check failed"
ls $P/plots/train_vars_sigVbkg/19_pho2_relEerr_* > /dev/null 2>&1 || die "training-variable plot of pho2_relEerr missing"
fi

date '+%F %T' > $L/datavmc.DONE
hb "DATAVMC DONE"
echo "$(date '+%F %H:%M') DONE  deploy_w123 dataVmc + BDT-input plots -> $L/datavmc_heartbeat.txt" >> /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/logs/PENDING_RESULTS.txt
