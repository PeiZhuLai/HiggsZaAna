#!/usr/bin/env bash
# Post-processing of the m_a >= 2 GeV trigdenom run (run only after every condor job is done):
#   1. collect the nominal trigeff payloads from HiggsDNA/eos_logs/<run> into cutflow JSONs
#      (../v2_trigdenom_20260927/collect_trigdenom.py, same dedup/completeness gate as m_a = 1)
#   2. plot the OL-denominator curves (same plotter/options as the m_a = 1 AN figure) and the
#      old denominator from the same events, for every mass point
#   3. plateau table per mass point
# Usage: RUN=Sig_MC_trigdenomV2mAge2 [SAMPLES=mA_M2,...] bash postprocess_mAge2.sh
set -uo pipefail
RUN="${RUN:-Sig_MC_trigdenomV2mAge2}"
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v2_trigdenom_mAge2_20261002
O=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v2_trigdenom_20260927
LOGS=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs/$RUN
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
PYROOT=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
CUT=$D/cutflow_list_$RUN
PLOTS=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/trigEffCompareVlepPt_trigdenomV2
[ "$RUN" != "Sig_MC_trigdenomV2mAge2" ] && PLOTS=$D/plots_$RUN
SAMPLES="${SAMPLES:-mA_M2,mA_M3,mA_M4,mA_M5,mA_M6,mA_M7,mA_M8,mA_M9,mA_M10,mA_M15,mA_M20,mA_M25,mA_M30}"
ERAS="${ERAS:-2022preEE,2022postEE,2023preBPix,2023postBPix,2024}"
mkdir -p $D/logs
: > $D/logs/post_${RUN}.log
for s in ${SAMPLES//,/ }; do
  m=${s#mA_M}
  $PY $O/collect_trigdenom.py --logs $LOGS --out $CUT --sample $s --eras $ERAS >> $D/logs/post_${RUN}.log 2>&1
  echo "$s collect rc=$?" | tee -a $D/logs/post_${RUN}.log
  cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts
  for yr in ${ERAS//,/ }; do
    [ -f $CUT/cutflow_Sig_MC_${s}_${yr}.json ] || continue
    env -i HOME=/afs/cern.ch/user/p/pelai PATH=/usr/bin:/bin PYTHONPATH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib \
      $PYROOT plot_trigEffVlepPt.py --ma $m --year $yr --in-dir $CUT --ord-suffix OL --out $PLOTS/denomOL >> $D/logs/plot_OL_${RUN}.log 2>&1
    echo "$s $yr plot OL rc=$?" | tee -a $D/logs/post_${RUN}.log
    env -i HOME=/afs/cern.ch/user/p/pelai PATH=/usr/bin:/bin PYTHONPATH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib \
      $PYROOT plot_trigEffVlepPt.py --ma $m --year $yr --in-dir $CUT --out $PLOTS/denomOld_sameRun >> $D/logs/plot_old_${RUN}.log 2>&1
  done
  $PY $D/plateau_table.py --sample $s --eras $ERAS --dir $CUT \
     --dir-old /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list \
     --json-out $CUT/plateau_${s}.json > $CUT/plateau_${s}.txt 2>&1
  echo "$s plateau rc=$?" | tee -a $D/logs/post_${RUN}.log
done
