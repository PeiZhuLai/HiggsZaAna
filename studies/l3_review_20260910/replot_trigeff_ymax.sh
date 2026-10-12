#!/usr/bin/env bash
# Re-draw every trigger-efficiency figure after raising the left y-axis maximum
# (plot_trigEffVlepPt.py default --ymax-left 160 -> 200, combined overlay follows it;
# plot_denom_compare.py frame 125 -> 170), so the ~100% plateau no longer touches the legend
# or the denominator note. Inputs and options are the ones used to make the existing figures:
#   production      : Plot/output/cutflow_list (script defaults)
#   m_a = 1 V2 study: v2_trigdenom_20260927/cutflow_list_trigdenomV2   (OL, N2, old denominator, compare)
#   m_a >= 2 V2     : v2_trigdenom_mAge2_20261002/cutflow_list_Sig_MC_trigdenomV2mAge2 (OL, old denominator)
# Runs serially (one process). Usage: bash replot_trigeff_ymax.sh
set -uo pipefail
ANA=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna
STUDY=$ANA/studies/l3_review_20260910
PYROOT=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
CUT1=$STUDY/v2_trigdenom_20260927/cutflow_list_trigdenomV2
CUT2=$STUDY/v2_trigdenom_mAge2_20261002/cutflow_list_Sig_MC_trigdenomV2mAge2
PLOTS=$ANA/Plot/plots/trigEffCompareVlepPt_trigdenomV2
LOG=$STUDY/logs_replot_trigeff_ymax
mkdir -p $LOG
: > $LOG/status.txt

run() {  # run <tag> <args...>
  local tag=$1; shift
  ( cd $ANA/Plot/scripts && env -i HOME=/afs/cern.ch/user/p/pelai PATH=/usr/bin:/bin \
      PYTHONPATH=$ANA/Plot/lib $PYROOT "$@" > $LOG/$tag.log 2>&1 )
  echo "$tag rc=$?" | tee -a $LOG/status.txt
}

run production    plot_trigEffVlepPt.py
run mA1_OL        plot_trigEffVlepPt.py --ma 1 --in-dir $CUT1 --ord-suffix OL --out $PLOTS/denomOL
run mA1_N2        plot_trigEffVlepPt.py --ma 1 --in-dir $CUT1 --ord-suffix N2 --out $PLOTS/denomN2
run mA1_old       plot_trigEffVlepPt.py --ma 1 --in-dir $CUT1 --out $PLOTS/denomOld_sameRun
run mAge2_OL      plot_trigEffVlepPt.py --in-dir $CUT2 --ord-suffix OL --out $PLOTS/denomOL
run mAge2_old     plot_trigEffVlepPt.py --in-dir $CUT2 --out $PLOTS/denomOld_sameRun
for e in 2022preEE 2022postEE 2023preBPix 2023postBPix 2024; do
  run compare_$e  $STUDY/v2_trigdenom_20260927/plot_denom_compare.py --in-dir $CUT1 --era $e --out $PLOTS/denomCompare
done
echo DONE >> $LOG/status.txt
