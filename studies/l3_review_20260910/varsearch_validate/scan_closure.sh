#!/usr/bin/env bash
# Closure + data-GOF scan for the variable-search BDTs (copy of lowmass_mfree_chain/scan_closure.sh, any m_a):
#   data envelope at the cut   build_data_ws_cut.py + fTest_ALP_turnOn (100 GOF toys) -> GOF per member
#   closure                    background-only DY pseudo-data (build_pseudodata_ws.py) + envelope refit (30 toys)
#                              + the production datacard with bkg/data swapped to it -> combine Significance
#                              (closure Z) and MultiDimFit (r_hat)
# DROP=1 builds the pseudo-data without the DY events of the BDT training pool (select_lib CLOSURE_DROP_TRAIN).
# The signal model is the production one. Production code and outputs are only read.
# Usage: [DROP=1] [NPAR=3] scan_closure.sh <prod|new|nodr> <mA> "<cut> <cut> ..."
# Output: /eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_scan/<model>[_drop]_mA<m>/{cut_<tag>/, summary.txt}
set -u
MODEL=$1; MA=$2; CUTS=$3; NPAR=${NPAR:-3}; DROP=${DROP:-0}
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
case $MODEL in
  prod) SB=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix ;;
  new)  SB=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored ;;
  nodr) SB=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored_nodR ;;
  cand_*) SB=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored_${MODEL#cand_}; [ -d $SB ] || { echo "no scored copy $SB"; exit 1; } ;;
  *) echo "model must be prod, new, nodr or cand_<name>"; exit 1 ;;
esac
TAGD=""; [ "$DROP" = 1 ] && TAGD=_drop
OUT=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_scan/${MODEL}${TAGD}_mA${MA}; mkdir -p $OUT
realDC=$FFIT/Datacard/output_Datacard_leptons/${MA}_pruned_datacard_leptons.txt
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $FFIT/.. && eval `scramv1 runtime -sh`
source $FFIT/setup.sh >/dev/null 2>&1
set -u

one_cut() {
  local cut=$1 tag=$(echo $1 | tr -d '.') OB=$OUT/cut_$(echo $1 | tr -d '.')
  rm -rf $OB; mkdir -p $OB
  ( cd $V/closure_scripts && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 CLOSURE_SCORED_BASE=$SB CLOSURE_DROP_TRAIN=$DROP \
      $CPY build_data_ws_cut.py --mA $MA --cut $cut --outbase $OB/data ) > $OB/data.build.log 2>&1
  mkdir -p $OB/data/fit/ft
  ( cd $FFIT/Background && ./bin/fTest_ALP_turnOn -i $OB/data/mA_M$MA/ws/run3.root \
      --saveMultiPdf $OB/data/fit/CMS-HGG_mva_13p6TeV_multipdf.root -D $OB/data/fit/ft \
      --mass_ALP $MA -c 1 --isFlashgg 0 --isData 0 -f data, --datfile $OB/$(basename $OB)_fTest.dat \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $OB/data.ftest.log 2>&1
  ( cd $V/closure_scripts && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 CLOSURE_SCORED_BASE=$SB CLOSURE_DROP_TRAIN=$DROP \
      $CPY build_pseudodata_ws.py --mA $MA --cut $cut --outbase $OB/ps ) > $OB/ps.build.log 2>&1
  mkdir -p $OB/ps/fit/ft
  ( cd $FFIT/Background && ./bin/fTest_ALP_turnOn -i $OB/ps/mA_M$MA/ws/run3.root \
      --saveMultiPdf $OB/ps/fit/CMS-HGG_mva_13p6TeV_multipdf.root -D $OB/ps/fit/ft \
      --mass_ALP $MA -c 1 --isFlashgg 0 --isData 0 -f data, --datfile $OB/$(basename $OB)_fTest.dat \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 30 ) > $OB/ps.ftest.log 2>&1
  sed "s#$FFIT/Background/ALP_BkgModel_ReReco/fit_results_run3/$MA/#$OB/ps/fit/#g" "$realDC" > $OB/ps/datacard.txt
  grep -q "$OB/ps/fit/CMS-HGG_mva_13p6TeV_multipdf.root" $OB/ps/datacard.txt || { echo "$cut FAIL datacard" > $OB/result.txt; return; }
  ( cd $OB/ps && text2workspace.py datacard.txt -o ws_t2w.root > t2w.log 2>&1
    combine -M Significance ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n R$tag > sig.log 2>&1
    combine -M MultiDimFit ws_t2w.root --algo singles --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n R$tag > rhat.log 2>&1
    rm -f ws_t2w.root )
  local nd=$(grep -oE 'data N=[0-9]+' $OB/data.build.log | grep -oE '[0-9]+')
  local np=$(grep -oE 'pseudo-data N=[0-9]+' $OB/ps.build.log | grep -oE '[0-9]+')
  local Z=$(grep -E 'Significance:' $OB/ps/sig.log | sed -E 's/.*Significance: //')
  local rh=$(grep -E '^\s+r :' $OB/ps/rhat.log | sed -E 's/^\s+r :\s+//; s/ \(68%\)//')
  local dg=$(grep -oE 'pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+' $OB/data/fit/ft/EnvelopeResults.txt 2>/dev/null | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' ')
  local pg=$(grep -oE 'pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+' $OB/ps/fit/ft/EnvelopeResults.txt 2>/dev/null | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' ')
  local fb=$(grep -c 'fallback' $OB/data/fit/ft/EnvelopeResults.txt 2>/dev/null)
  echo "$cut dataN=${nd:-NA} psN=${np:-NA} closureZ=${Z:-NA} rhat=${rh:-NA} data_fallback=${fb:-NA} | data GOF: $dg | ps GOF: $pg" > $OB/result.txt
}

running=0
for cut in $CUTS; do
  one_cut $cut & running=$((running+1))
  if [ $running -ge $NPAR ]; then wait -n; running=$((running-1)); fi
done
wait
for cut in $CUTS; do cat $OUT/cut_$(echo $cut | tr -d '.')/result.txt; done > $OUT/summary.txt
echo "SCAN DONE $MODEL$TAGD mA$MA $(date '+%F %T')" >> $OUT/summary.txt
