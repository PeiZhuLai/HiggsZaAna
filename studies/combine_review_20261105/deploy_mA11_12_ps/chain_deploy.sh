#!/usr/bin/env bash
# Production deployment of the mA11/mA12 decision of the Combine review: the multi-start refit of the data fit also
# starts the exponent of the one-term shapes (Pow1, Lau1, Exp1) from 15/50/85% of its range at these two mass points.
# Reason: the stored Pow1 was a local minimum (steep exponent near -9.2); the Lau1 solution lies inside the Pow1 ranges
# with a lower NLL (mA11 478.975 vs 479.81, mA12 525.625 vs 527.117) and the multi-start only scanned turn-on x width.
# The steep Pow1 truth failed the bias requirement at both points. With the exponent scan Pow1 converges to the Lau1
# solution, Lau1/Exp1/Bern5 are unchanged, every truth passes the bias requirement (range_mA11_12_pscan/bias/summary_ps.txt)
# and the sideband goodness of fit passes (range_mA11_12_pscan/gof/).
# Local stages only; condor checks (bias, impacts) are submitted afterwards by a separate script.
#  1  build: fTest rebuilt from src/fTest_ALP_turnOn.deploy.cpp (side-study source with the mass list 11, 12 written into
#     the code, src/make_source.py); library and PdfModelBuilder.cc unchanged. Regression: mA1, mA2 and mA19 (not
#     redeployed) must reproduce production exactly.
#  2  staging fTest of mA11 and mA12 -> work/fit_results/<m> (production untouched); must be identical to the side study
#     range_mA11_12_pscan/fit_results/<m>.
#     The chain stops here (STOP_AFTER=2 by default) for the inspection.
#  3  install: production fit_results_run3/<m> and its AllFitResults files -> fit_results_run3_bak_deploy_mA11_12_ps/,
#     staged directories copied in; makeBkgPlots blind only.
#  4  text2workspace + blind expected limits + check_params for mA11 and mA12, both cards; old files backed up.
# Every overwritten product is backed up first (suffix bak_deploy_mA11_12_ps); an existing backup target stops the chain.
# Markers: logs/stage<N>.DONE, logs/chain.DONE or logs/chain.FAILED; heartbeat logs/heartbeat.txt.
# START_STAGE=<N> resumes, STOP_AFTER=<N> stops after stage N. At most 4 local processes.
set -uo pipefail
ST=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105
D=$ST/deploy_mA11_12_ps; L=$D/logs; STG=$ST/open_ranges_5pts; PS=$ST/range_mA11_12_pscan
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
FF=$CMSSW_SRC/flashggFinalFit; BG=$FF/Background; FR=$BG/ALP_BkgModel_ReReco/fit_results_run3
FRB=$BG/ALP_BkgModel_ReReco/fit_results_run3_bak_deploy_mA11_12_ps
C=$FF/Combine; MV=/eos/home-p/pelai/HZa/root_MVAcut/data
SRCP=$PS/src/fTest_ALP_turnOn.pscan.cpp; SRC=$D/src/fTest_ALP_turnOn.deploy.cpp
SIDE=$PS/fit_results
SFX=bak_deploy_mA11_12_ps; REG=$ST/PENDING_RESULTS.txt; M=CMS-HGG_mva_13p6TeV_multipdf.root
MASSES="11 12"; REGM="1 2 19"
START=${START_STAGE:-1}; STOP=${STOP_AFTER:-2}
mkdir -p $L $D/src $D/regress $D/work/fit_results
run() { [ "$START" -le "$1" ] && [ "$1" -le "$STOP" ]; }
hb() { echo "$(date '+%F %T') $*" | tee -a $L/heartbeat.txt; }
die() { hb "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/chain.FAILED; echo "$(date '+%F %H:%M') DIED deploy_mA11_12_ps: $*  -> $L/heartbeat.txt" >> $REG; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }
newdir() { [ -e "$1" ] && die "backup/target $1 already exists"; true; }
envelope() { grep -oE "pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+" $1 | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' '; }
scan() { grep -m1 'exponent scan in the multi-start' $1 | awk '{print $NF}'; }
rm -f $L/chain.DONE $L/chain.FAILED
unset EXPSCAN_MASSES ORDER2_MASSES

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_SRC && eval $(scramv1 runtime -sh)
set -u

ftest_one() {  # <m> <out dir>; production command line of shellScripts/bkg/fit_bkg.sh (mhLow 100 for mA1, 95 otherwise, 100 GOF toys)
  local m=$1 P=$2 lo=95
  [ "$m" = 1 ] && lo=100
  mkdir -p $P/HZAmassInde_fTest
  ( cd $BG && ./bin/fTest_ALP_turnOn -i $MV/mA_M$m/ws/run3.root --saveMultiPdf $P/CMS-HGG_mva_13p6TeV_multipdf.root \
      -D $P/HZAmassInde_fTest --mass_ALP $m -c 1 --isFlashgg 0 --isData 0 -f data, \
      --mhLow $lo --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $P/ftest.stdout.log 2>&1
  echo "rc=$?" > $P/.ftest_rc
}
export -f ftest_one; export BG MV

if run 1; then
hb "[1] build"
cd $BG
cmp -s src/PdfModelBuilder.cc $PS/Background/src/PdfModelBuilder.cc || die "production PdfModelBuilder.cc differs from the side-study copy"
[ -e bin/fTest_ALP_turnOn.$SFX ] || cp -p bin/fTest_ALP_turnOn bin/fTest_ALP_turnOn.$SFX || die "backup fTest binary"
python3 $D/src/make_source.py $SRCP $SRC > $L/make_source.log 2>&1 || die "make_source.py (see $L/make_source.log)"
diff $SRCP $SRC > $L/fTest_source_vs_pscan.diff
diff $ST/deploy_mA11/src/fTest_ALP_turnOn.deploy.cpp $SRC > $L/fTest_source_vs_previous_production.diff
( cd $BG && make -f makefile vars ) > $L/make_vars.txt
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $L/make_vars.txt); LIBS=$(sed -n 's/^LIBS: *//p' $L/make_vars.txt)
T=$(date +%s)
g++ -o $BG/bin/fTest_ALP_turnOn $CXXF $SRC -L$BG/lib -lBackgroundProfileFitting $LIBS > $L/build_ftest.log 2>&1 || die "g++ fTest (see $L/build_ftest.log)"
fresh $BG/bin/fTest_ALP_turnOn $T || die "fTest binary not rebuilt"
readelf -d $BG/bin/fTest_ALP_turnOn | grep -q "rpath: \[$BG/lib" || die "fTest rpath is not the production library"
hb "    fTest rebuilt; regression mA$REGM against production"
for m in $REGM; do rm -rf $D/regress/$m; echo "$m $D/regress/$m"; done | xargs -P3 -L1 bash -c 'ftest_one $0 $1'
for m in $REGM; do
  grep -q 'rc=0' $D/regress/$m/.ftest_rc || die "regression fTest mA$m failed"
  [ "$(scan $D/regress/$m/ftest.stdout.log)" = off ] || die "exponent scan not off at mA$m"
  python3 $STG/compare_ws.py $FR/$m/$M $D/regress/$m/$M > $L/regress_$m.txt 2>&1
  tail -1 $L/regress_$m.txt | grep -qx IDENTICAL || die "mA$m differs from production ($L/regress_$m.txt)"
done
hb "    regression: mA$REGM IDENTICAL to production, exponent scan off"
date '+%F %T' > $L/stage1.DONE
fi

if run 2; then
hb "[2] staging fTest mA $MASSES -> $D/work/fit_results"
[ -f $L/stage1.DONE ] || die "stage 1 not done"
T=$(date +%s)
for m in $MASSES; do rm -rf $D/work/fit_results/$m; echo "$m $D/work/fit_results/$m"; done | xargs -P4 -L1 bash -c 'ftest_one $0 $1'
for m in $MASSES; do
  grep -q 'rc=0' $D/work/fit_results/$m/.ftest_rc && fresh $D/work/fit_results/$m/$M $T || die "staging fTest mA$m failed"
  [ "$(scan $D/work/fit_results/$m/ftest.stdout.log)" = on ] || die "exponent scan not on at mA$m"
  hb "    mA$m $(envelope $D/work/fit_results/$m/HZAmassInde_fTest/EnvelopeResults.txt)"
  python3 $STG/compare_ws.py $SIDE/$m/$M $D/work/fit_results/$m/$M > $L/compare_side_$m.txt 2>&1
  tail -1 $L/compare_side_$m.txt | grep -qx IDENTICAL || die "staging mA$m differs from the side study ($L/compare_side_$m.txt)"
done
date '+%F %T' > $L/stage2.DONE
fi

if run 3; then
hb "[3] install mA $MASSES"
[ -f $L/stage2.DONE ] || die "stage 2 not done"
newdir $FRB; mkdir -p $FRB/AllFitResults_blind $FRB/AllFitResults_unblind
for m in $MASSES; do
  [ -d $FR/$m ] || die "no production directory $FR/$m"
  [ -s $D/work/fit_results/$m/$M ] || die "no staged multipdf mA$m"
  mv $FR/$m $FRB/$m || die "move $m to backup"
  for b in blind unblind; do
    find $FR/AllFitResults_$b -maxdepth 1 -type f -name "*_${m}.*" -exec mv -t $FRB/AllFitResults_$b/ {} + || die "move AllFitResults_$b $m"
  done
  cp -a $D/work/fit_results/$m $FR/$m || die "install $m"
  cmp -s $D/work/fit_results/$m/$M $FR/$m/$M || die "installed multipdf mA$m differs from staging"
done
plots_one() {
  local m=$1 P=$FR/$1 lo=95
  [ "$m" = 1 ] && lo=100
  ( cd $BG && ./bin/makeBkgPlots_ALP -b $P/CMS-HGG_mva_13p6TeV_multipdf.root -d $P/BkgPlots_blind --total_OutDir $FR/AllFitResults_blind \
      -o $P/BkgPlots_blind.root --sqrts 13p6TeV --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --maVal $m \
      --mhLow $lo --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi 172.13 -c 0 --isFlashgg 0 --doBands ) > $P/makeBkgPlots.stdout.log 2>&1
  echo "rc=$?" > $P/.plots_rc
}
export -f plots_one; export FR
T=$(date +%s)
printf '%s\n' $MASSES | xargs -P4 -I{} bash -c 'plots_one {}'
for m in $MASSES; do
  grep -q 'rc=0' $FR/$m/.plots_rc && fresh $FR/$m/BkgPlots_blind.root $T || die "makeBkgPlots mA$m failed"
  fresh $FR/AllFitResults_blind/bkgplot_$m.pdf $T || die "no fresh AllFitResults_blind/bkgplot_$m.pdf"
done
hb "    installed; makeBkgPlots blind done (old unblind plots of these masses only in $FRB/AllFitResults_unblind)"
date '+%F %T' > $L/stage3.DONE
fi

if run 4; then
hb "[4] text2workspace + blind limits mA $MASSES"
[ -f $L/stage3.DONE ] || die "stage 3 not done"
for k in xsec wilson; do
  newdir $C/root_t2w_${k}_$SFX; newdir $C/output_combine_results_${k}_$SFX
  mkdir -p $C/root_t2w_${k}_$SFX $C/output_combine_results_${k}_$SFX
  for m in $MASSES; do
    mv $C/root_t2w_$k/${m}_Datacard_leptons.root $C/root_t2w_${k}_$SFX/ || die "backup t2w $k $m"
    mv $C/output_combine_results_$k/higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_${k}_$SFX/ || die "backup limit $k $m"
  done
done
lim_one() {
  local m=$1 k=$2 sfx="" W
  [ $k = xsec ] && sfx=_xsec
  W=$D/work/lim_${k}_$m; rm -rf $W; mkdir -p $W
  ( cd $C && python3 RunText2Workspace.py --batch local --queue workday --mA $m --card $k ) > $L/t2w_${k}_$m.log 2>&1
  ( cd $W && combine -M AsymptoticLimits $FF/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons$sfx.txt --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/limit_${k}_$m.log 2>&1 \
    && mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_$k/ )
  python3 $FF/Datacard/check_params.py $FF/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons$sfx.txt $C/root_t2w_$k/${m}_Datacard_leptons.root > $L/params_${k}_$m.txt 2>&1
}
export -f lim_one; export C D L FF
T=$(date +%s)
for m in $MASSES; do printf '%s xsec\n%s wilson\n' $m $m; done | xargs -P4 -L1 bash -c 'lim_one $0 $1'
for m in $MASSES; do for k in xsec wilson; do
  fresh $C/root_t2w_$k/${m}_Datacard_leptons.root $T || die "t2w $k mA$m"
  fresh $C/output_combine_results_$k/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "limit $k mA$m"
  tail -1 $L/params_${k}_$m.txt | grep -qx 'ALL OK' || die "check_params $k mA$m ($L/params_${k}_$m.txt)"
done; done
hb "    4 workspaces, 4 blind expected limits, check_params ALL OK"
date '+%F %T' > $L/stage4.DONE
fi
hb "chain stages $START-$STOP done"
[ "$STOP" -ge 4 ] && date '+%F %T' > $L/chain.DONE
echo "$(date '+%F %H:%M') DONE deploy_mA11_12_ps stages $START-$STOP -> $L/heartbeat.txt" >> $REG
