#!/usr/bin/env bash
# Production deployment of the Combine-review changes (user-approved plan: ~/doc/HZa/review_comments/Combine_review/deploy_plan_20261111.md).
# Local stages only; the condor checks (bias 110 jobs, impacts 60 jobs) are submitted separately after this chain.
#  1  build: PdfModelBuilder.cc opens the step/smearing ranges at mA10, 15, 21-23 (24-29 already); library rebuilt
#     (make lib: only PdfModelBuilder.o changes, the interface is unchanged since 10-01); fTest_ALP_turnOn rebuilt from
#     the 10-23 production source + patch_databinned + patch_multistart + patch_extraenv (EXTRA_ENV_PDFS unset: no effect).
#     Regression: mA19 (not redeployed) must reproduce production, mA10 must reproduce staging fit_results_ms.
#  2  background: fTest of mA10, 15, 21-29 (old directories -> fit_results_run3_bak_deploy_20261111/), each multipdf
#     identical to staging fit_results_ms; makeBkgPlots blind only (no --unblind)
#  3  signal workspaces renamed -> Signal/outdir_{ele,mu}/signalFit/output_named/ (originals untouched)
#  4  datacards: makeDatacard (lumi 3-nuisance scheme, signal MC stat) + xsec cards -> output_Datacard_leptons_prenaming/,
#     apply_naming.py (rename, trigger split, unwired params dropped, signal shapes -> output_named) -> output_Datacard_leptons/
#  REDO_CARDS=1 REDO_LIMITS=1 START_STAGE=4: regenerate after a makeDatacard change, keeping the pre-deployment backups.
#  5  text2workspace (60) + blind expected limits (60) + check_params (every param acts on model_s); compare with backups
# Every overwritten product is backed up first (suffix bak_deploy_20261111); a backup target that exists stops the chain.
# Markers: logs/stage<N>.DONE, logs/chain.DONE or logs/chain.FAILED; heartbeat logs/heartbeat.txt. START_STAGE=<N> resumes.
# At most 5 local processes.
set -uo pipefail
ST=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105
D=$ST/deploy_20261111; L=$D/logs; STG=$ST/open_ranges_5pts
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
FF=$CMSSW_SRC/flashggFinalFit; BG=$FF/Background; FR=$BG/ALP_BkgModel_ReReco/fit_results_run3
FRB=$BG/ALP_BkgModel_ReReco/fit_results_run3_bak_deploy_20261111
C=$FF/Combine; DC=$FF/Datacard/output_Datacard_leptons
MV=/eos/home-p/pelai/HZa/root_MVAcut/data
REV=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_envrevert_20261023/src/fTest_ALP_turnOn.revert.cpp
SFX=bak_deploy_20261111; REG=$ST/PENDING_RESULTS.txt
MASSES="10 15 21 22 23 24 25 26 27 28 29"; ALL=$(seq 1 30)
START=${START_STAGE:-1}
mkdir -p $L $D/src $D/regress $D/work
run() { [ "$START" -le "$1" ]; }
hb() { echo "$(date '+%F %T') $*" | tee -a $L/heartbeat.txt; }
die() { hb "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/chain.FAILED; echo "$(date '+%F %H:%M') DIED deploy_20261111: $*  -> $L/heartbeat.txt" >> $REG; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }
newdir() { [ -e "$1" ] && die "backup/target $1 already exists"; true; }
envelope() { grep -oE "pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+" $1 | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' '; }
rm -f $L/chain.DONE $L/chain.FAILED

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_SRC && eval $(scramv1 runtime -sh)
set -u

ftest_one() {  # <m> <out dir>; production command line of shellScripts/bkg/fit_bkg.sh (mhLow 95 for m > 1, 100 GOF toys)
  local m=$1 P=$2
  mkdir -p $P/HZAmassInde_fTest
  ( cd $BG && ./bin/fTest_ALP_turnOn -i $MV/mA_M$m/ws/run3.root --saveMultiPdf $P/CMS-HGG_mva_13p6TeV_multipdf.root \
      -D $P/HZAmassInde_fTest --mass_ALP $m -c 1 --isFlashgg 0 --isData 0 -f data, \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $P/ftest.stdout.log 2>&1
  echo "rc=$?" > $P/.ftest_rc
}
export -f ftest_one; export BG MV

if run 1; then
hb "[1] build"
cd $BG
for f in src/PdfModelBuilder.cc lib/libBackgroundProfileFitting.so bin/fTest_ALP_turnOn; do
  [ -e $f.$SFX ] || cp -p $f $f.$SFX || die "backup $f"
done
python3 - <<'EOF' || die "edit PdfModelBuilder.cc"
p = 'src/PdfModelBuilder.cc'; s = open(p).read()
old = 'inline bool openStepGausRanges(int massALP) { return massALP == 0 || (massALP >= 24 && massALP <= 29); }'
new = 'inline bool openStepGausRanges(int massALP) { return massALP == 0 || massALP == 10 || massALP == 15 || (massALP >= 21 && massALP <= 29); }'
if new in s: print('already edited')
else:
    assert s.count(old) == 1, 'openStepGausRanges line not found'
    s = s.replace(old, new)
    s = s.replace('// [PZ 2026-09-27] mA24-29: after the FSR-fix re-production every envelope member had its',
                  '// mA10, 15, 21-23 opened as well (Combine review: wide-truth bias fails there with the per-mass windows,\n// fitted turn-ons on the window edges).\n// [PZ 2026-09-27] mA24-29: after the FSR-fix re-production every envelope member had its', 1)
    open(p, 'w').write(s); print('edited')
EOF
diff $BG/src/PdfModelBuilder.cc.$SFX $BG/src/PdfModelBuilder.cc > $L/PdfModelBuilder.diff
T=$(date +%s)
( cd $BG && make -f makefile lib ) > $L/make_lib.log 2>&1 || die "make lib (see $L/make_lib.log)"
fresh $BG/lib/libBackgroundProfileFitting.so $T || die "library not rebuilt"
python3 $STG/patch_databinned.py $REV $D/src/fTest.databinned.cpp || die "patch_databinned"
python3 $STG/patch_multistart.py $D/src/fTest.databinned.cpp $D/src/fTest.ms.cpp || die "patch_multistart"
python3 $STG/patch_extraenv.py $D/src/fTest.ms.cpp $D/src/fTest_ALP_turnOn.deploy.cpp || die "patch_extraenv"
sed -i "s#\"../../tdrStyle/#\"$FF/tdrStyle/#" $D/src/fTest_ALP_turnOn.deploy.cpp
( cd $BG && make -f makefile vars ) > $L/make_vars.txt
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $L/make_vars.txt); LIBS=$(sed -n 's/^LIBS: *//p' $L/make_vars.txt)
T=$(date +%s)
g++ -o $BG/bin/fTest_ALP_turnOn $CXXF $D/src/fTest_ALP_turnOn.deploy.cpp -L$BG/lib -lBackgroundProfileFitting $LIBS > $L/build_ftest.log 2>&1 || die "g++ fTest (see $L/build_ftest.log)"
fresh $BG/bin/fTest_ALP_turnOn $T || die "fTest binary not rebuilt"
readelf -d $BG/bin/fTest_ALP_turnOn | grep -q "rpath: \[$BG/lib" || die "fTest rpath is not the production library"
hb "    library and fTest rebuilt; regression mA19 (production) and mA10 (staging)"
rm -rf $D/regress/19 $D/regress/10
printf '19 %s\n10 %s\n' $D/regress/19 $D/regress/10 | xargs -P2 -L1 bash -c 'ftest_one $0 $1'
grep -q 'rc=0' $D/regress/19/.ftest_rc && grep -q 'rc=0' $D/regress/10/.ftest_rc || die "regression fTest failed"
M=CMS-HGG_mva_13p6TeV_multipdf.root
python3 $STG/compare_ws.py $FR/19/$M $D/regress/19/$M > $L/regress_19.txt 2>&1; tail -1 $L/regress_19.txt | grep -qx IDENTICAL || die "mA19 differs from production ($L/regress_19.txt)"
python3 $STG/compare_ws.py $STG/fit_results_ms/10/$M $D/regress/10/$M > $L/regress_10.txt 2>&1; tail -1 $L/regress_10.txt | grep -qx IDENTICAL || die "mA10 differs from staging ($L/regress_10.txt)"
hb "    regression: mA19 IDENTICAL to production, mA10 IDENTICAL to staging"
date '+%F %T' > $L/stage1.DONE
fi

if run 2; then
hb "[2] background mA $MASSES"
newdir $FRB; mkdir -p $FRB/AllFitResults_blind $FRB/AllFitResults_unblind
for m in $MASSES; do
  [ -d $FR/$m ] || die "no production directory $FR/$m"
  mv $FR/$m $FRB/$m || die "move $m to backup"
  for b in blind unblind; do
    find $FR/AllFitResults_$b -maxdepth 1 -type f -name "*_${m}.*" -exec mv -t $FRB/AllFitResults_$b/ {} + || die "move AllFitResults_$b $m"
  done
done
T=$(date +%s)
for m in $MASSES; do echo "$m $FR/$m"; done | xargs -P5 -L1 bash -c 'ftest_one $0 $1'
M=CMS-HGG_mva_13p6TeV_multipdf.root
for m in $MASSES; do
  grep -q 'rc=0' $FR/$m/.ftest_rc && fresh $FR/$m/$M $T || die "fTest mA$m failed"
  python3 $STG/compare_ws.py $STG/fit_results_ms/$m/$M $FR/$m/$M > $L/compare_staging_$m.txt 2>&1
  tail -1 $L/compare_staging_$m.txt | grep -qx IDENTICAL || die "mA$m differs from staging ($L/compare_staging_$m.txt)"
  hb "    mA$m $(envelope $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt)| identical to staging"
done
plots_one() {
  local m=$1 P=$FR/$1
  ( cd $BG && ./bin/makeBkgPlots_ALP -b $P/CMS-HGG_mva_13p6TeV_multipdf.root -d $P/BkgPlots_blind --total_OutDir $FR/AllFitResults_blind \
      -o $P/BkgPlots_blind.root --sqrts 13p6TeV --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --maVal $m \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi 172.13 -c 0 --isFlashgg 0 --doBands ) > $P/makeBkgPlots.stdout.log 2>&1
  echo "rc=$?" > $P/.plots_rc
}
export -f plots_one; export FR
T=$(date +%s)
printf '%s\n' $MASSES | xargs -P5 -I{} bash -c 'plots_one {}'
for m in $MASSES; do
  grep -q 'rc=0' $FR/$m/.plots_rc && fresh $FR/$m/BkgPlots_blind.root $T || die "makeBkgPlots mA$m failed"
  fresh $FR/AllFitResults_blind/bkgplot_$m.pdf $T || die "no fresh AllFitResults_blind/bkgplot_$m.pdf"
done
hb "    makeBkgPlots blind done (no unblind plots for the redeployed masses; old ones in $FRB/AllFitResults_unblind)"
date '+%F %T' > $L/stage2.DONE
fi

if run 3; then
hb "[3] signal workspaces renamed"
for c in ele mu; do
  O=$FF/Signal/outdir_$c/signalFit/output_named; newdir $O
  python3 $FF/Signal/scripts/rename_ws.py $O $(ls $FF/Signal/outdir_$c/signalFit/output/[1-9]*_CMS-HGG_sigfit_*_${c}_Hm125.root) > $L/rename_ws_$c.log 2>&1 || die "rename_ws $c ($L/rename_ws_$c.log)"
  n=$(grep -c '^\[OK\]' $L/rename_ws_$c.log); e=$(grep -c '^\[ERROR\]' $L/rename_ws_$c.log)
  [ "$n" = 70 ] && [ "$e" = 0 ] || die "rename_ws $c: $n OK, $e ERROR (expected 70 OK)"
  hb "    $c: 70 files renamed, all checks OK"
done
date '+%F %T' > $L/stage3.DONE
fi

if run 4; then
hb "[4] datacards"
P=${DC}_prenaming
if [ "${STAGE4_CHECKS_ONLY:-0}" = 1 ]; then
# Resume after a stop in the checks: the backup, the regenerated cards and apply_naming are kept from the previous run
# (rerunning them would copy the new cards over the backup); only the checks below are redone.
[ -d ${DC}_$SFX ] && [ -d $P ] && [ -s $L/apply_naming.log ] || die "STAGE4_CHECKS_ONLY: backup, prenaming cards or apply_naming.log missing"
hb "    checks only, on the cards of the previous run (cards $(stat -c %y $DC/1_pruned_datacard_leptons.txt | cut -c1-19))"
T=$(date +%s)
else
if [ "${REDO_CARDS:-0}" = 1 ]; then
# Regenerate after a makeDatacard change (interp_yield rename): the pre-deployment backup is kept as it is; the
# previously deployed cards go to work/cards_v1 for the comparison; prenaming cards and apply_naming.log are redone.
[ -d ${DC}_$SFX ] && [ -d $P ] || die "REDO_CARDS: backup or prenaming dir missing"
V1=$D/work/cards_v1; newdir $V1; mkdir -p $V1
for m in $ALL; do cp -p $DC/${m}_pruned_datacard_leptons.txt $DC/${m}_pruned_datacard_leptons_xsec.txt $V1/ || die "snapshot cards $m"; done
[ -e $L/apply_naming.log.v1 ] && die "apply_naming.log.v1 exists"; mv $L/apply_naming.log $L/apply_naming.log.v1
hb "    REDO_CARDS: pre-deployment backup kept, previously deployed cards in $V1"
else
newdir ${DC}_$SFX; newdir $P
cp -a $DC ${DC}_$SFX || die "backup datacards"
fi
T=$(date +%s)
( cd $FF/Datacard && bash 2_runjob_gen_datacard_makeDatacard.sh ) > $L/makeDatacard.log 2>&1 || die "makeDatacard ($L/makeDatacard.log)"
for m in $ALL; do for s in "" _xsec; do fresh $DC/${m}_pruned_datacard_leptons$s.txt $T || die "card $m$s not regenerated"; done; done
mkdir -p $P
for m in $ALL; do cp -p $DC/${m}_pruned_datacard_leptons.txt $DC/${m}_pruned_datacard_leptons_xsec.txt $P/ || die "copy prenaming $m"; done
for m in $ALL; do for s in "" _xsec; do
  python3 $FF/Datacard/apply_naming.py $P/${m}_pruned_datacard_leptons$s.txt $DC/${m}_pruned_datacard_leptons$s.txt \
    --signal-ws-sibling output_named --drop-unwired-params >> $L/apply_naming.log 2>&1 || die "apply_naming $m$s"
done; done
grep 'dropped' $L/apply_naming.log | sort | uniq -c | tee $L/dropped_params.txt
fi
cmp -s $FF/Datacard/naming/rename_dict.yml ~/doc/HZa/review_comments/Combine_review/rename_dict_draft.yml && \
cmp -s $FF/Datacard/naming/split_dict.yml ~/doc/HZa/review_comments/Combine_review/split_dict_draft.yml || die "production naming dicts differ from the docdir drafts"
python3 ~/doc/HZa/review_comments/Combine_review/check_systematics_draft.py $(ls $P/*_pruned_datacard_leptons*.txt | grep -E '/[0-9]+_') > $L/check_systematics.log 2>&1
cp -p ~/doc/HZa/review_comments/Combine_review/systematics_draft_check.txt $L/
fresh $L/systematics_draft_check.txt $T && grep -qx 'problems: 0' $L/systematics_draft_check.txt || die "systematics.yml check ($L/systematics_draft_check.txt)"
if [ "${REDO_CARDS:-0}" = 1 ]; then
python3 - $DC $D/work/cards_v1 $L/cards_vs_v1.txt <<'EOF' || die "cards differ from the previous deployment beyond the interp_yield rename ($L/cards_vs_v1.txt)"
import os, re, sys
new, v1, out = sys.argv[1:4]; bad = 0; rows = []
E = {'2022preEE': '2022', '2022postEE': '2022EE', '2023preBPix': '2023', '2023postBPix': '2023BPix', '2024': '2024'}
for m in range(1, 31):
    for s in ('', '_xsec'):
        n = '%d_pruned_datacard_leptons%s.txt' % (m, s)
        a = open(os.path.join(v1, n)).read(); b = open(os.path.join(new, n)).read()
        k = len(re.findall(r'\binterpolate_yield_', a))
        a2 = re.sub(r'\binterpolate_yield_(2022preEE|2022postEE|2023preBPix|2023postBPix|2024)_(ele|mu)\b',
                    lambda x: 'CMS_NPS25014_interp_yield_%s_%s' % (E[x.group(1)], x.group(2)), a)
        same = [l.split() for l in a2.splitlines()] == [l.split() for l in b.splitlines()]
        bad += 0 if same and 'interpolate_yield' not in b else 1
        rows.append('%s %-34s old-name occurrences renamed: %d' % ('[OK]   ' if same else '[DIFF] ', n, k))
open(out, 'w').write('\n'.join(rows) + '\n%s\n' % ('ALL OK' if not bad else 'PROBLEMS %d' % bad)); sys.exit(1 if bad else 0)
EOF
hb "    REDO_CARDS: cards identical to the previous deployment up to the interp_yield rename ($L/cards_vs_v1.txt)"
fi
python3 - $DC ${DC}_$SFX $L/card_diff.txt <<'EOF' || die "card diff"
import os, sys
new, old, out = sys.argv[1:4]; res = []; bad = 0
TYPES = ('lnN', 'shape', 'param', 'shapeN', 'lnU', 'gmN')
def parse(p):
    nu, sh = {}, {}
    for l in open(p):
        f = l.split()
        if len(f) > 1 and f[1] in TYPES: nu[f[0]] = f[1]
        if f and f[0] == 'shapes': sh[f[1]] = f[3:]
    return nu, sh
for m in range(1, 31):
    for s in ('', '_xsec'):
        a, b = parse(os.path.join(old, '%d_pruned_datacard_leptons%s.txt' % (m, s))), parse(os.path.join(new, '%d_pruned_datacard_leptons%s.txt' % (m, s)))
        gone = sorted(set(a[0]) - set(b[0])); added = sorted(set(b[0]) - set(a[0]))
        bk = all(a[1][k] == b[1][k] for k in ('bkg_mass', 'data_obs'))
        sig = all(b[1][k][0].startswith('/') and '/signalFit/output_named/' in b[1][k][0] and os.path.isfile(b[1][k][0]) for k in b[1] if k.startswith('ggH_'))
        lumi = sorted(k for k in b[0] if k.startswith('lumi'))
        mc = sum(1 for k in b[0] if '_mcstat_' in k)
        ok = bk and sig and lumi == ['lumi_13p6TeV_2024', 'lumi_13p6TeV_222324', 'lumi_13p6TeV_2324'] and mc == 10
        bad += 0 if ok else 1
        res.append('%s mA%-2d%-5s nuisances %d -> %d, bkg/data shapes unchanged %s, signal -> output_named %s, lumi %s, mcstat %d\n    removed: %s\n    added: %s'
                   % ('[OK]   ' if ok else '[ERROR]', m, s, len(a[0]), len(b[0]), bk, sig, lumi, mc, ' '.join(gone), ' '.join(added)))
open(out, 'w').write('\n'.join(res) + '\n%s\n' % ('ALL OK' if not bad else 'PROBLEMS %d' % bad))
print('card diff:', 'ALL OK' if not bad else 'PROBLEMS %d' % bad); sys.exit(1 if bad else 0)
EOF
hb "    60 cards regenerated, renamed; card diff ALL OK; systematics.yml: $(head -1 $L/systematics_draft_check.txt), problems 0"
date '+%F %T' > $L/stage4.DONE
fi

if run 5; then
hb "[5] text2workspace + blind limits"
if [ "${REDO_LIMITS:-0}" = 1 ]; then
# Rerun after REDO_CARDS: the pre-deployment backups are kept; the previously deployed limits go to work/limits_v1_<k>
# and must be reproduced exactly (the rename does not change the model).
for k in xsec wilson; do
  [ -d $C/root_t2w_${k}_$SFX ] && [ -d $C/output_combine_results_${k}_$SFX ] || die "REDO_LIMITS: backups missing"
  newdir $D/work/limits_v1_$k; cp -a $C/output_combine_results_$k $D/work/limits_v1_$k || die "snapshot limits $k"
done
[ -e $L/limits_compare.txt.v1 ] && die "limits_compare.txt.v1 exists"; mv $L/limits_compare.txt $L/limits_compare.txt.v1
hb "    REDO_LIMITS: pre-deployment backups kept, previously deployed limits in $D/work/limits_v1_{xsec,wilson}"
else
for k in xsec wilson; do
  newdir $C/root_t2w_${k}_$SFX; newdir $C/output_combine_results_${k}_$SFX
  mv $C/root_t2w_$k $C/root_t2w_${k}_$SFX && mkdir $C/root_t2w_$k || die "backup root_t2w_$k"
  mv $C/output_combine_results_$k $C/output_combine_results_${k}_$SFX && mkdir $C/output_combine_results_$k || die "backup output_combine_results_$k"
done
fi
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
for m in $ALL; do printf '%s xsec\n%s wilson\n' $m $m; done | xargs -P5 -L1 bash -c 'lim_one $0 $1'
for m in $ALL; do for k in xsec wilson; do
  fresh $C/root_t2w_$k/${m}_Datacard_leptons.root $T || die "t2w $k mA$m"
  fresh $C/output_combine_results_$k/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "limit $k mA$m"
  tail -1 $L/params_${k}_$m.txt | grep -qx 'ALL OK' || die "check_params $k mA$m ($L/params_${k}_$m.txt)"
done; done
python3 - $C $SFX $L/limits_compare.txt <<'EOF'
import sys, ROOT
C, sfx, out = sys.argv[1:4]; rows = ['card    mA | median before -> after (rel)  | 2.5% 16% 84% 97.5% rel. change']
def q(p):
    f = ROOT.TFile.Open(p); t = f.Get('limit'); d = {}
    for e in t: d[round(e.quantileExpected, 3)] = e.limit
    f.Close(); return d
for k in ('xsec', 'wilson'):
    for m in range(1, 31):
        n = 'higgsCombine%d.AsymptoticLimits.mH125.38.root' % m
        a = q('%s/output_combine_results_%s_%s/%s' % (C, k, sfx, n)); b = q('%s/output_combine_results_%s/%s' % (C, k, n))
        rel = lambda x: (b[x] / a[x] - 1) * 100
        rows.append('%-6s %3d | %.5f -> %.5f (%+.1f%%) | %s' % (k, m, a[0.5], b[0.5], rel(0.5), ' '.join('%+.1f' % rel(x) for x in (0.025, 0.16, 0.84, 0.975))))
open(out, 'w').write('\n'.join(rows) + '\n'); print('\n'.join(rows))
EOF
if [ "${REDO_LIMITS:-0}" = 1 ]; then
python3 - $C $D/work $L/limits_vs_v1.txt <<'EOF' || die "limits differ from the previous deployment ($L/limits_vs_v1.txt)"
import sys, ROOT
C, W, out = sys.argv[1:4]; bad = 0; rows = []
def q(p):
    f = ROOT.TFile.Open(p); t = f.Get('limit'); d = sorted((round(e.quantileExpected, 3), e.limit) for e in t); f.Close(); return d
for k in ('xsec', 'wilson'):
    for m in range(1, 31):
        n = 'higgsCombine%d.AsymptoticLimits.mH125.38.root' % m
        a = q('%s/limits_v1_%s/%s' % (W, k, n)); b = q('%s/output_combine_results_%s/%s' % (C, k, n))
        same = a == b; bad += 0 if same else 1
        rows.append('%s %-6s mA%-2d %s' % ('[OK]   ' if same else '[DIFF] ', k, m, '' if same else '%s -> %s' % (a, b)))
open(out, 'w').write('\n'.join(rows) + '\n%s\n' % ('ALL IDENTICAL' if not bad else 'DIFFERENT %d' % bad)); sys.exit(1 if bad else 0)
EOF
hb "    REDO_LIMITS: all 60 limits identical to the previous deployment ($L/limits_vs_v1.txt)"
fi
hb "    60 workspaces, 60 blind limits, check_params ALL OK on 60 cards; comparison in $L/limits_compare.txt"
date '+%F %T' > $L/stage5.DONE
fi
hb "chain done"
date '+%F %T' > $L/chain.DONE
echo "$(date '+%F %H:%M') DONE deploy_20261111 local stages 1-5 -> $L/heartbeat.txt, limits_compare.txt" >> $REG
