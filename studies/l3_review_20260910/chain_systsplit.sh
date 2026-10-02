#!/usr/bin/env bash
# Signal systematic datasets rebuilt with the fixed apply_bdt_sig.py (2026-09-27).
#
# Bug: the p2root train/validation/test split is row-index % 10 per file, so each systematic
# sample's 'test' tree was a different random 30% than the nominal one, and apply_bdt_sig.py
# joined it with the nominal pass set -> the Up/Down signal datasets held only the overlap
# (median 48% of the nominal yield, down to 3%), calcPhotonSyst turned that into 5%-capped rate
# terms, and mean/sigma variations were computed on random subsets. Fix (apply_bdt_sig.py, backup
# .bak_systsplit_20260927): systematics read the 'inclusive' tree restricted to the nominal test
# keys and apply their own BDT cut. Tested on mA5/mA20 x 2022preEE/2024: all Up/Down within
# +-0.6% of nominal, cross-channel variations exactly 1, nominal unchanged.
#
# Steps: wait for hza-csev2829 -> apply_bdt_sig (14 mA x 5 eras) -> Tree2WS signal ->
# calcPhotonSyst -> signalFit -> plotter -> effSigma -> datacards (with CSEV) -> t2w -> limits
# -> expected impacts (all 30). The nominal signal model is unchanged, so fTest and bias are not
# rerun.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
EOSMVA=/eos/home-p/pelai/HZa/root_MVAcut
S1=$D/systsplit_single; mkdir -p $S1
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
MASS14="1 2 3 4 5 6 7 8 9 10 15 20 25 30"
TAG=before_systsplit_$(date +%Y%m%d_%H%M)
die() { echo "STOPPING: $*"; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

echo "=============== [0] wait for hza-csev2829 ($(date '+%T')) ==============="
while systemctl --user is-active -q hza-csev2829; do sleep 60; done
grep -q STOPPING $L/csev2829.log && die "csev2829 chain stopped -- check it first"
old_imp=$(grep -oE "impact: [0-9]+ job\(s\) submitted to cluster [0-9]+" $L/csev2829.log | grep -oE "[0-9]+\.?$" | tr -d .)
[ -n "$old_imp" ] && { condor_rm $old_imp >/dev/null 2>&1; echo "  removed superseded impact cluster $old_imp"; }
echo "  go ($(date '+%T'))"

echo "=============== [1] backups ($TAG) ==============="
V=$C/output_combine_results_variants; mkdir -p $V
find $C/output_combine_results -maxdepth 1 -type f ! -regex '.*/higgsCombine[0-9]+\.AsymptoticLimits\.mH125\.38\.root' -exec mv {} $V/ \;
cp -rp $C/output_combine_results $V/output_combine_results_$TAG
cp -rp $F/Datacard/output_Datacard_leptons $F/Datacard/output_Datacard_leptons_$TAG
for ch in ele mu; do
  mkdir -p $F/Signal/sigsyst_variants_$TAG/$ch
  cp -rp $F/Signal/outdir_$ch/calcPhotonSyst $F/Signal/sigsyst_variants_$TAG/$ch/ 2>/dev/null
  cp -p $F/Signal/outdir_$ch/signalFit/output/*_CMS-HGG_sigfit_*.root $F/Signal/sigsyst_variants_$TAG/$ch/ 2>/dev/null
done
echo "  done (the EOS signal trees are rewritten in place; their nominal trees do not change)"

set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv
export PYTHONPATH="${PYTHONPATH:-}:$F/tools:$F/Signal/tools"

echo "=============== [2] apply_bdt_sig (14 mA x 5 eras, 3 in parallel) ==============="
T=$(date +%s)
sed 's/^max_parallel=6$/max_parallel=3/' $F/shellScripts/mva/run_apply_bdt_sig_6jobs.sh > $S1/run_apply_bdt_sig_3jobs.sh
grep -q '^max_parallel=3$' $S1/run_apply_bdt_sig_3jobs.sh || die "parallel edit failed"
bash $S1/run_apply_bdt_sig_3jobs.sh > $L/systsplit_apply_sig.log 2>&1
bad=0; for m in $MASS14; do for y in $YEARS; do fresh $EOSMVA/sig/mA_M$m/output_$y.root $T || { echo "  sig mA_M$m $y stale"; bad=1; }; done; done
[ $bad -eq 0 ] || die "apply_bdt_sig outputs incomplete"
# acceptance: every Up/Down within 3% of nominal, and at least one ratio > 1
env -i HOME=$HOME PATH=/usr/bin:/bin /eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3 - $EOSMVA/sig <<'PY' > $L/systsplit_ratios.txt 2>&1 || { tail -5 $L/systsplit_ratios.txt; die "syst/nominal ratios outside 3%"; }
import glob, sys, uproot
rs = []
for f in sorted(glob.glob(sys.argv[1] + "/mA_M*/output_*.root")):
    t = uproot.open(f)["DiphotonTree"]
    for lep in ("ele", "mu"):
        base = "ggh_125_Za_%s_13p6TeV_cat0" % lep
        nom = t[base].arrays(["weight"], library="np")["weight"].sum()
        for k in t.keys():
            k = k.split(";")[0]
            if k.startswith(base + "_") and nom > 0:
                rs.append((t[k].arrays(["weight"], library="np")["weight"].sum() / nom, f, lep, k))
lo, hi = min(rs), max(rs)
print("n=%d  min %.4f (%s %s %s)  max %.4f (%s %s %s)" % (len(rs), lo[0], *lo[1:], hi[0], *hi[1:]))
sys.exit(0 if (lo[0] > 0.97 and hi[0] < 1.03 and hi[0] > 1.0) else 1)
PY
cat $L/systsplit_ratios.txt | sed 's/^/  /'

echo "=============== [3] Tree2WS (signal only) ==============="
T=$(date +%s)
sed -e 's/^mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)$/mAs_data=()/' \
    $F/Trees2WS/run_tree2ws.sh > $S1/run_tree2ws_sig.sh
grep -q '^mAs_data=()$' $S1/run_tree2ws_sig.sh || die "tree2ws signal-only copy failed"
cd $F/Trees2WS && bash $S1/run_tree2ws_sig.sh > $L/systsplit_tree2ws.log 2>&1
for m in $MASS14; do
  n=$(find $EOSMVA/sig/mA_M$m/ws_Tree2WS -name '*.root' -newermt "@$T" 2>/dev/null | wc -l)
  [ "$n" -ge 10 ] || die "sig ws mA_M$m: only $n fresh files"
done
echo "  signal workspaces 14 x 10 fresh"

echo "=============== [4] signal model (calcPhotonSyst -> signalFit -> plotter -> effSigma) ==============="
T=$(date +%s)
cd $F/shellScripts
for s in 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter 5_runjob_sig_plotEffSigma; do
  echo "  [$(date '+%H:%M:%S')] $s"
  bash $F/shellScripts/sig_sys/$s.sh > $L/systsplit_sig_$s.log 2>&1
done
bad=0
for m in $MASS14; do for y in $YEARS; do for ch in ele mu; do
  fresh $F/Signal/outdir_$ch/signalFit/output/${m}_CMS-HGG_sigfit_${y}_${ch}_Hm125.root $T || { echo "  sigfit mA$m $y $ch stale"; bad=1; }
done; done; done
[ $bad -eq 0 ] || die "signal fits incomplete"
echo "  signal fits: 140/140 fresh"

echo "=============== [5] datacards (with CSEV), t2w, limits ==============="
T=$(date +%s)
cd $F/Datacard
sh 1_runjob_gen_datacard_makeYields.sh  > $L/systsplit_dc_1.log 2>&1
sh 2_runjob_gen_datacard_makeDatacard.sh > $L/systsplit_dc_2.log 2>&1
sh 3_rysn_datacard.sh                    > $L/systsplit_dc_3.log 2>&1
for m in $(seq 1 30); do
  dc=$F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt
  fresh $dc $T || die "datacard mA$m stale"
  grep -q "^CMS_hza_photon_csev" $dc || die "mA$m lost CMS_hza_photon_csev"
  r=$(grep -c rateParam $dc); case " 1 2 3 4 5 6 7 8 9 10 15 20 25 30 " in *" $m "*) want=0;; *) want=10;; esac
  [ "$r" -eq "$want" ] || die "mA$m has $r rateParam lines (want $want)"
done
T=$(date +%s)
cd $C && sh 2_text2ws.sh > $L/systsplit_t2w.log 2>&1 && sh 1_makeLimits.sh > $L/systsplit_limits.log 2>&1
bad=0; for m in $(seq 1 30); do fresh $C/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || bad=1; done
[ $bad -eq 0 ] || die "limits stale"
python3 - $C/output_combine_results $V/output_combine_results_$TAG > $L/systsplit_changed.txt <<'PY'
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def med(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m)); t = f.Get("limit"); v = None
    for e in t:
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
for m in range(1, 31):
    n, o = med(sys.argv[1], m), med(sys.argv[2], m)
    print(m, "%.5f %.5f %.4f" % (o, n, n / o))
PY
echo "  mA old(after CSEV) new ratio:"; sed 's/^/    /' $L/systsplit_changed.txt

echo "=============== [6] condor: expected impacts (all 30) ==============="
bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $(seq 1 30) 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  impact: /'
echo "SUBMITTED $(date '+%F %T')"
