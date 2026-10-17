#!/usr/bin/env bash
# Step F (closure): background-only pseudo-data closure (AN Sec 6.4.2, Fig. sculpt_closure) redone for m_a = 1-3 with
# the low_relE2_w123 BDT, using the production pseudodata_closure scripts exactly as chain_closure30_dyveto.sh did for
# all 30 points (CLOSURE_FSRFIX=1, working points from MVAcut_points_run3.json, real datacard = Wilson card as for the
# other 27 points); the other 27 points are untouched. At most 3 processes.
#   1 backup of the m_a 1-3 products + closure_results.json + plots -> backup_prod/closure_before.tar
#   2 build_pseudodata_ws.py --mA 1,2,3   3 fTest envelopes   4 limits / significance / fitdiag
#   5 extract_results.py + plot_closure.py (all 30 points)    6 unblind background-fit panels of m_a 1-3
# Output: pseudodata_closure/{pseudodata/mA_M{1,2,3},fit_results_run3_pseudodata/{1,2,3},output_*,closure_results.json,plots}
# Logs: logs/closure_*.log; marker logs/closure.DONE | logs/closure.FAILED
set -u
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
PC=$FFIT/pseudodata_closure; S=$PC/scripts
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
L=$D/logs
export CLOSURE_FSRFIX=1
die() { echo "$(date '+%F %T') $*" > $L/closure.FAILED; echo "STOPPING: $*"; exit 1; }
rm -f $L/closure.DONE $L/closure.FAILED
echo "=== [1] backup $(date '+%T')"
[ -s $D/backup_prod/closure_before.tar ] || ( cd $PC && tar -cf $D/backup_prod/closure_before.tar pseudodata/mA_M1 pseudodata/mA_M2 pseudodata/mA_M3 \
  fit_results_run3_pseudodata/1 fit_results_run3_pseudodata/2 fit_results_run3_pseudodata/3 fit_results_run3_pseudodata/AllFitResults_unblind \
  datacards output_limits output_significance output_fitdiag closure_results.json plots ) || die "backup"
cp -p $PC/closure_results.json $L/closure_results_before.json
echo "=== [2] pseudo-data workspaces m_a 1-3 $(date '+%T')"
( cd $S && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 $CPY build_pseudodata_ws.py --mA 1,2,3 ) > $L/closure_build.log 2>&1 || die "build (closure_build.log)"
grep -E "^mA" $L/closure_build.log | sed 's/^/    /'
echo "=== [3] envelope fits $(date '+%T')"
T=$(date +%s)
ALPMASS="1 2 3" NPROC=3 GOF_TOYS=100 bash $S/run_ftest_pseudodata.sh > $L/closure_ftest.log 2>&1
for m in 1 2 3; do f=$PC/fit_results_run3_pseudodata/$m/CMS-HGG_mva_13p6TeV_multipdf.root; [ -s $f ] && [ $(stat -c %Y $f) -ge $T ] || die "multipdf mA$m stale"; done
echo "=== [4] limits / significance / fitdiag $(date '+%T')"
T=$(date +%s)
ALPMASS="1 2 3" NPROC=3 bash $S/run_limits_pseudodata.sh > $L/closure_limits.log 2>&1
n=$(find $PC/output_limits -name '*.root' -newermt "@$T" | wc -l); [ $n -ge 3 ] || die "only $n fresh limit files"
echo "=== [5] results table + figure $(date '+%T')"
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd $FFIT/.. && eval `scramv1 runtime -sh`; set -u
( cd $S && python3 extract_results.py ) > $L/closure_extract.log 2>&1 || die "extract_results"
( cd $S && env -i HOME=$HOME PATH=/usr/bin:/bin $CPY plot_closure.py ) > $L/closure_plot.log 2>&1 || die "plot_closure"
python3 - $L/closure_results_before.json $PC/closure_results.json <<'PY' | tee $L/closure_compare.txt
import json, sys
o, n = (json.load(open(p)) for p in sys.argv[1:3])
def rows(d):
    d = d.get("results", d)
    return {int(k): v for k, v in (d.items() if isinstance(d, dict) else [(r["mA"], r) for r in d])}
o, n = rows(o), rows(n)
for m in sorted(n):
    if o.get(m) != n[m]:
        print("  mA%-2d before %s\n        after  %s" % (m, o.get(m), n[m]))
print("  changed points:", [m for m in sorted(n) if o.get(m) != n[m]])
PY
echo "=== [6] background-fit panels m_a 1-3 $(date '+%T')"
set +u; source $FFIT/setup.sh >/dev/null 2>&1; set -u
U=$PC/fit_results_run3_pseudodata/AllFitResults_unblind; mkdir -p $U
panel() {
  local m=$1 D=$PC/fit_results_run3_pseudodata/$1
  local c=( -b $D/CMS-HGG_mva_13p6TeV_multipdf.root --sqrts 13p6TeV --isMultiPdf --useBinnedData
            --massStep 2.5 --mhVal 125.0 --maVal $m --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135
            --intLumi 172.13 -c 0 --isFlashgg 0 --dataLegend "MC pseudo-data" --ratioYTitle "MC - Fit" --unblind
            -d $D/BkgPlots_unblind --total_OutDir $U -o $D/BkgPlots_unblind.root )
  ( cd $FFIT/Background && timeout 45m ./bin/makeBkgPlots_ALP "${c[@]}" --doBands > $D/makeBkgPlots.stdout.log 2>&1 ) \
    || { echo "[mA=$m] bands timed out/failed -> no bands"; ( cd $FFIT/Background && ./bin/makeBkgPlots_ALP "${c[@]}" >> $D/makeBkgPlots.stdout.log 2>&1 ); }
  echo "[mA=$m] panels done"
}
export -f panel; export PC FFIT U
printf '1\n2\n3\n' | xargs -P 3 -I{} bash -c 'panel {}' > $L/closure_bkgplots.log 2>&1
cat $L/closure_bkgplots.log
date '+%F %T' > $L/closure.DONE
echo "DONE $(date '+%F %T')"
