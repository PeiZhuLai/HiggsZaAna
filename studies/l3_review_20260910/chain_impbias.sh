#!/usr/bin/env bash
# Watch the FSR-fix impacts + high-stat bias condor clusters, then merge/collect the bias.
#   impacts : 12747408 (mA5) + 12747410 (29 mA)
#   bias    : 12747409 (mA5 x 10 chunks) + 12747411 (29 mA x 10 chunks)
# Done-markers are products, never exit codes or queue emptiness:
#   impacts  Combine/output_impacts/<mA>_impacts.json newer than impbias.t0
#   bias     bias_outputs_highstat/mA_<m>/chunk<c>/BiasFits/*_split*_fits.root newer than impbias.t0
# The merge reads the same frozen workspace the bias jobs used (root_t2w_fsrfix_20260925).
# collect_bias_results.sh copies from BOTH bias_outputs (the 2026-08-02 1000-toy study) and
# bias_outputs_highstat into plots_bias, so the pre-FSR-fix plots are moved aside first.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FZ=$C/root_t2w_fsrfix_20260925
CLUSTERS="12747408 12747409 12747410 12747411"
T0=$(cat $L/impbias.t0)
HB=$L/impbias_heartbeat.txt

count_products() {
  imp=0; for m in $(seq 1 30); do f=$C/output_impacts/${m}_impacts.json; [ -s $f ] && [ $(stat -c %Y $f) -ge $T0 ] && imp=$((imp+1)); done
  chunks=0; for m in $(seq 1 30); do for c in $(seq 0 9); do
    compgen -G "$BN/bias_outputs_highstat/mA_${m}/chunk${c}/BiasFits/*_split*_fits.root" >/dev/null && chunks=$((chunks+1)); done; done
}

while true; do
  q=$(condor_q $CLUSTERS -af JobStatus 2>/dev/null)
  if [ $? -ne 0 ]; then sleep 300; continue; fi   # schedd busy: do not mistake it for "empty"
  idle=$(grep -c '^1$' <<<"$q"); run=$(grep -c '^2$' <<<"$q"); held=$(grep -c '^5$' <<<"$q")
  count_products
  {
    echo "$(date '+%F %T')  idle=$idle running=$run held=$held"
    echo "  impacts json fresh : $imp/30"
    echo "  bias chunks done   : $chunks/300"
    [ "$held" -gt 0 ] && condor_q $CLUSTERS -hold -af ClusterId ProcId NumJobStarts HoldReason 2>/dev/null | head -20 | sed 's/^/  HELD /'
  } > $HB
  [ $((idle + run + held)) -eq 0 ] && break
  # held with NumJobStarts>=3 will never be released by periodic_release; stop watching only when
  # nothing else is left, so the report below lists them
  [ $((idle + run)) -eq 0 ] && [ "$held" -gt 0 ] && { echo "  only held jobs left" >> $HB; break; }
  sleep 1800
done

echo "=============== impacts/bias condor finished $(date '+%F %T') ==============="
cat $HB
count_products
missing_imp=""; for m in $(seq 1 30); do f=$C/output_impacts/${m}_impacts.json; { [ -s $f ] && [ $(stat -c %Y $f) -ge $T0 ]; } || missing_imp="$missing_imp $m"; done
missing_bias=""; for m in $(seq 1 30); do n=0; for c in $(seq 0 9); do
  compgen -G "$BN/bias_outputs_highstat/mA_${m}/chunk${c}/BiasFits/*_split*_fits.root" >/dev/null && n=$((n+1)); done
  [ $n -eq 10 ] || missing_bias="$missing_bias mA$m:$n/10"; done
echo "impacts missing :${missing_imp:- none}"
echo "bias incomplete :${missing_bias:- none}"

echo "=============== merge high-stat bias (frozen workspace) ==============="
ROOT_DATACARD_PATH=$FZ bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh $(seq 1 30) > $L/bias_merge.log 2>&1
echo "  merge rc=$?  (log bias_merge.log)"
nj=0; bad=""; for m in $(seq 1 30); do
  j=$BN/bias_outputs_highstat/mA_${m}/merged/BiasJson/${m}_gaussfit.json
  p=$BN/bias_outputs_highstat/mA_${m}/merged/plots_bias/${m}_r_vs_bias.pdf
  if [ -s $j ] && [ -s $p ]; then nj=$((nj+1)); else bad="$bad $m"; fi; done
echo "  merged mass points with gaussfit json + r_vs_bias pdf: $nj/30${bad:+  missing:$bad}"

echo "=============== collect bias plots ==============="
[ -d $BN/bias_outputs ] && mv $BN/bias_outputs $BN/bias_outputs_preFSRfix_20260925
[ -d $BN/plots_bias ]   && mv $BN/plots_bias   $BN/plots_bias_preFSRfix_20260925
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh

if [ -z "$missing_imp" ] && [ -z "$missing_bias" ] && [ $nj -eq 30 ]; then
  echo "VERDICT: PASSED"
else
  echo "VERDICT: FAILED (see lists above)"
fi
echo "Not run: Update AN (needs the user's go-ahead to push)."
