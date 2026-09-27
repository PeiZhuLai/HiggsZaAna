#!/usr/bin/env bash
# flashgg stage of the FSR-fix rerun: working points -> MVA cut -> Tree2WS -> background
# envelope -> signal model -> datacards -> blind expected limits -> limit plots.
#
# Starts only after the dataVmc chain (hza-datavmc2) has PASSED, because the working-point
# optimization (ALP_Optimization.py) reads the dataVmc merged histograms.
#
# Deliberately NOT run here:
#  * RUN_UPDATE_AN: 1_grand.sh does `git add . && git commit && git push` in the shared AN
#    repo. Pushing the AN needs the user's explicit go-ahead, and that repo has unrelated
#    uncommitted files that a blanket `git add .` would sweep in.
#  * observed / unblinded limits and the unblinding procedure.
#
# Rules carried over from earlier incidents (see memory):
#  * background fits run LOCALLY (xargs -P6), never on condor: a condor eviction while
#    saveMultiPdf is writing truncates the multipdf to a ~63-byte corrupt file.
#  * a multipdf is only accepted if it opens and getNumPdfs() > 0 -- "the log has no error"
#    has passed a crashed fit before.
#  * `cmsenv` is an interactive shell function; the flashgg scripts call it (and the
#    datacard ones run under set -e), so it is defined and exported here.
#  * every stage is judged by output mtimes >= the stage start, never by exit codes; the
#    previous outputs are snapshotted first so the before/after comparison survives the
#    overwrite (1_makeLimits.sh starts with `rm -fr output_combine_results`).
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
EOSMVA=/eos/home-p/pelai/HZa/root_MVAcut
SNAP=$D/snapshot_preFSRfix_20260924
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
MASS14="1 2 3 4 5 6 7 8 9 10 15 20 25 30"
MASS30=$(seq 1 30 | tr '\n' ' ')
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"

fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }
# START_STAGE=N skips stages < N (the CMSSW env block always runs). Stages >= 5 are safe to
# re-enter: the background stage only refits mass points whose envelope is missing, broken,
# or older than its input workspace.
START=${START_STAGE:-0}
run() { [ "$START" -le "$1" ]; }
die()   { echo "STOPPING: $*"; exit 1; }

if run 0; then
echo "=============== [0/9] wait for dataVmc2 ==============="
while [ "$(systemctl --user is-active hza-datavmc2 2>/dev/null)" = "active" ]; do sleep 120; done
grep -q "^VERDICT: PASSED" $L/datavmc2.log || die "dataVmc2 did not pass -- see $L/datavmc2.log"
echo "  dataVmc2 passed"
fi

echo
if run 1; then
echo "=============== [1/9] snapshot the pre-FSR-fix results ==============="
mkdir -p "$SNAP"
cp -n $PLOT/output/MVAcut_points_run3.json               "$SNAP/" 2>/dev/null
[ -d "$SNAP/output_combine_results" ] || cp -r $F/Combine/output_combine_results "$SNAP/"
[ -d "$SNAP/fit_results_run3" ]       || cp -r $F/Background/ALP_BkgModel_ReReco/fit_results_run3 "$SNAP/"
[ -d "$SNAP/plot_limits" ]            || cp -r $F/Plots/plot_limits "$SNAP/"
mkdir -p "$SNAP/effSigma"; cp -n $F/Signal/outdir_*/signalFit/*effSigma*.json "$SNAP/effSigma/" 2>/dev/null
cp -n $F/Signal/outdir_ele/*effSigma*.json $F/Signal/outdir_mu/*effSigma*.json "$SNAP/effSigma/" 2>/dev/null
echo "  snapshot: $(du -sh "$SNAP" | cut -f1) in $SNAP"
fi

echo
if run 2; then
echo "=============== [2/9] working points (ALP_Optimization -> collect -> R scan) ==============="
T=$(date +%s)
OPT=$PLOT/plots/optimize_run3UL
# stale copies go OUTSIDE Plot/plots: sync_figures.sh copies that whole tree into the AN
STALE=$PLOT/plots_variants/optimize_run3UL_stale_preFSRfix_20260924; mkdir -p $STALE
mv $OPT/nCat*_all_M*.json $STALE/ 2>/dev/null
cd $PLOT
export PYTHONPATH="${PYTHONPATH:-}:$PLOT/lib"
$ALP_PY scripts/ALP_Optimization.py -y run3 -o $OPT --region 2 -p --sigVSscore -s --doOpt -c 2 --inputTag sideband_rwgt > $L/wp_opt_c2.log 2>&1
$ALP_PY scripts/ALP_Optimization.py -y run3 -o $OPT --region 2 -p --sigVSscore -s --doOpt -c 1 --inputTag sideband_rwgt > $L/wp_opt_c1.log 2>&1
bad=0; for m in $MASS14; do fresh $OPT/nCat1_all_M$m.json $T || { echo "  missing/stale nCat1_all_M$m.json"; bad=1; }; done
[ $bad -eq 0 ] || die "ALP_Optimization did not produce every mass point (logs: wp_opt_c*.log)"
$ALP_PY scripts/collect_MVAcut_points_run3.py > $L/wp_collect.log 2>&1
$ALP_PY scripts/scan_score_R_significance.py --write-json > $L/wp_scan.log 2>&1
fresh $PLOT/output/MVAcut_points_run3.json $T || die "MVAcut_points_run3.json was not rewritten"
$ALP_PY - <<PY
import json
old = {r["mA"]: r for r in json.load(open("$SNAP/MVAcut_points_run3.json"))["results"]}
new = {r["mA"]: r for r in json.load(open("$PLOT/output/MVAcut_points_run3.json"))["results"]}
print("  %4s %10s %10s   %10s %10s" % ("mA", "cut_old", "cut_new", "Z_old", "Z_new"))
for m in sorted(set(old) | set(new)):
    o, n = old.get(m, {}), new.get(m, {})
    print("  %4s %10s %10s   %10s %10s" % (m, o.get("MVAcut"), n.get("MVAcut"), o.get("Significance"), n.get("Significance")))
assert len(new) >= 14, "fewer than 14 working points"
PY
[ $? -eq 0 ] || die "working-point table check failed"
fi

# The datacards scale the 16 interpolated mass points by eff(target)/eff(anchor) read from
# sigEfficiencyVmA_*_quadratic_interp_ma_points.json, which signal_eff_sumw.py writes in the
# post-scan half of Plot/1_runPlot.sh. Run 2 (2026-09-24) skipped that half, so its datacards
# used the June curves (input_base run3_bdt_scored_nominal). Regenerate whenever the JSON is
# older than the working points, and refuse to write datacards from a stale one.
EFFJ_ELE=$PLOT/output/sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json
EFFJ_MU=$PLOT/output/sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json
check_effjson() {
  env -u PYTHONPATH -u PYTHONHOME $ALP_PY - "$EFFJ_ELE" "$EFFJ_MU" "$PLOT/output/MVAcut_points_run3.json" <<'PY'
import json, os, sys
wp = sys.argv[3]; ok = True
for p in sys.argv[1:3]:
    d = json.load(open(p)); meta = d.get("meta", {})
    why = []
    if not meta.get("input_base", "").endswith("run3_bdt_scored_fsrfix"): why.append("input_base=%s" % meta.get("input_base"))
    if os.path.getmtime(p) < os.path.getmtime(wp): why.append("older than MVAcut_points_run3.json")
    vals = d.get("values", {})
    nz = sum(1 for y in vals for m in range(1, 31) if float(vals[y].get(str(m), 0)) > 0)
    if len(vals) != 5 or nz != 150: why.append("%d years, %d/150 non-zero (year,mA) entries" % (len(vals), nz))
    print("  %s: %s" % (os.path.basename(p), "OK" if not why else "STALE (" + "; ".join(why) + ")"))
    ok = ok and not why
sys.exit(0 if ok else 1)
PY
}
if run 7; then
echo "=============== [2b/9] interpolation efficiency JSON (signal_eff_sumw.py) ==============="
if ! check_effjson; then
  cd $PLOT
  export PYTHONPATH="${PYTHONPATH:-}:$PLOT/lib"
  $ALP_PY scripts/signal_eff_sumw.py > $L/wp_signal_eff_sumw.log 2>&1
  check_effjson || die "signal efficiency JSON still stale after signal_eff_sumw.py (wp_signal_eff_sumw.log)"
fi
fi

# ---- from here on: CMSSW environment -------------------------------------------------
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }
export -f cmsenv
cmsenv
export PYTHONPATH="${PYTHONPATH:-}:$F/tools:$F/Signal/tools"
set -u
echo "  [env] CMSSW_BASE=$CMSSW_BASE  python=$(which python3)"

echo
if run 3; then
echo "=============== [3/9] MVA cut: data (30 mA) + signal (14 mA x 5 eras) ==============="
T=$(date +%s)
for d in data sig; do
  [ -d $EOSMVA/${d}_preFSRfix_20260924 ] || mv $EOSMVA/$d $EOSMVA/${d}_preFSRfix_20260924
done
( ok=0; for a in 1 2 3 4 5; do
    if python3 $F/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py > $L/mva_apply_data_try$a.log 2>&1; then ok=1; break; fi
    sleep 30
  done; [ $ok = 1 ] ) &
pd=$!
bash $F/shellScripts/mva/run_apply_bdt_sig_6jobs.sh > $L/mva_apply_sig.log 2>&1
wait $pd || die "apply_bdt_data.py failed 5 times"
bad=0
for m in $MASS30; do fresh $EOSMVA/data/mA_M$m/run3.root $T || { echo "  data mA_M$m stale/missing"; bad=1; }; done
for m in $MASS14; do for y in $YEARS; do fresh $EOSMVA/sig/mA_M$m/output_$y.root $T || { echo "  sig mA_M$m $y stale/missing"; bad=1; }; done; done
[ $bad -eq 0 ] || die "MVA cut outputs incomplete"
echo "  MVA cut outputs: 30 data + 70 signal, all fresh"
# fresh + non-empty is not enough: mA17 data run3.root (2026-09-24) opened fine, reported 90
# entries, and threw an I/O error at entry 41. Read every tree end to end.
python3 $D/audit_mvacut_trees.py > $L/mva_audit.log 2>&1 || { tail -15 $L/mva_audit.log; die "MVA cut trees fail a full read"; }
tail -1 $L/mva_audit.log
fi

echo
if run 4; then
echo "=============== [4/9] Tree2WS ==============="
T=$(date +%s)
cd $F/Trees2WS && bash run_tree2ws.sh > $L/tree2ws.log 2>&1
bad=0
for m in $MASS30; do fresh $EOSMVA/data/mA_M$m/ws/run3.root $T || { echo "  data ws mA_M$m stale/missing"; bad=1; }; done
for m in $MASS14; do
  n=$(find $EOSMVA/sig/mA_M$m/ws_Tree2WS -name '*.root' -newermt "@$T" 2>/dev/null | wc -l)
  [ "$n" -ge 10 ] || { echo "  sig ws mA_M$m: only $n fresh files (expect >= 10 = 5 eras x 2 leptons)"; bad=1; }
done
[ $bad -eq 0 ] || die "Tree2WS outputs incomplete (log: tree2ws.log)"
# a Tree2WS crash leaves a fresh ~600-byte file with no workspace in it (mA17, 2026-09-24),
# which passed the check above; open every one and require a RooWorkspace
python3 - "$EOSMVA" <<'PYWS' || die "some workspaces are empty shells (see above)"
import glob, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B = sys.argv[1]
paths = ["%s/data/mA_M%d/ws/run3.root" % (B, m) for m in range(1, 31)]
paths += sorted(glob.glob("%s/sig/mA_M*/ws_Tree2WS/ws_*.root" % B))
bad = []
for p in paths:
    try:
        f = ROOT.TFile.Open(p)
        # signal ws keep the RooWorkspace inside a TDirectory (CMS_hza_workspace/CMS_hza_workspace)
        def has_ws(d):
            for k in d.GetListOfKeys():
                o = k.ReadObj()
                if o.InheritsFrom("RooWorkspace") or (o.InheritsFrom("TDirectory") and has_ws(o)):
                    return True
            return False
        ok = has_ws(f)
        f.Close()
    except OSError:
        ok = False
    if not ok: bad.append(p)
print("  workspaces opened: %d   without a RooWorkspace: %d" % (len(paths), len(bad)))
for p in bad: print("    EMPTY", p)
sys.exit(1 if bad else 0)
PYWS
echo "  workspaces: 30 data + 14 signal mass points, all fresh and non-empty"
fi

echo
if run 5; then
echo "=============== [5/9] background envelopes (local, 6 in parallel) ==============="
T=$(date +%s)
BK=$F/Background
cd $BK && make -f makefile clean > $L/bkg_build.log 2>&1 && make -f makefile -B >> $L/bkg_build.log 2>&1
[ -x $BK/bin/fTest_ALP_turnOn ] && [ -x $BK/bin/makeBkgPlots_ALP ] || die "background build failed (bkg_build.log)"
TMPD=$D/fit_bkg_single; mkdir -p $TMPD
python3 - "$F/shellScripts/bkg/fit_bkg.sh" "$TMPD" <<'PY'
import re, sys
src, outd = sys.argv[1], sys.argv[2]
s = open(src).read()
s = s.replace("make -f makefile clean && make -f makefile -B", ": # built once by chain_flashgg.sh")
for m in range(1, 31):
    t = re.sub(r'^massList=\( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 \)$',
               'massList=( %d )' % m, s, count=1, flags=re.M)
    assert t != s, "massList line not found"
    t = t.replace('failed_log="failed_mass_points.log"', 'failed_log="failed_mass_points_m%d.log"' % m)
    t = t.replace('summary_log="$path_out_bkg/background_fit_summary.csv"',
                  'summary_log="$path_out_bkg/background_fit_summary_m%d.csv"' % m)
    open("%s/fit_bkg_m%d.sh" % (outd, m), "w").write(t)
print("wrote 30 single-mass copies")
PY
FR=$BK/ALP_BkgModel_ReReco/fit_results_run3
# only refit what is missing, broken, or older than its input workspace -- so re-entering this
# stage after fixing one mass point does not refit the other 29 (mA23 alone takes ~98 min)
todo=$(python3 - "$FR" "$EOSMVA" <<'PYTODO'
import os, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
fr, B = sys.argv[1], sys.argv[2]
order = [23, 4, 1, 2, 3] + [m for m in range(5, 31) if m != 23]
need = []
for m in order:
    p = "%s/%d/CMS-HGG_mva_13p6TeV_multipdf.root" % (fr, m)
    ws = "%s/data/mA_M%d/ws/run3.root" % (B, m)
    good = os.path.exists(p) and os.path.getsize(p) >= 1024 and os.path.getmtime(p) >= os.path.getmtime(ws)
    if good:
        try:
            f = ROOT.TFile.Open(p); n = -1
            for k in f.GetListOfKeys():
                o = k.ReadObj()
                if o.InheritsFrom("RooWorkspace"):
                    it = o.allPdfs().createIterator(); pdf = it.Next()
                    while pdf:
                        if pdf.InheritsFrom("RooMultiPdf"): n = max(n, pdf.getNumPdfs())
                        pdf = it.Next()
            f.Close(); good = n > 0
        except OSError:
            good = False
    if not good: need.append(str(m))
print(" ".join(need))
PYTODO
)
echo "  mass points to (re)fit: ${todo:-none}"
[ -n "$todo" ] && echo $todo | tr ' ' '\n' | xargs -P6 -I{} bash -c "bash $TMPD/fit_bkg_m{}.sh > $L/bkg_m{}.log 2>&1"
FR=$BK/ALP_BkgModel_ReReco/fit_results_run3
head -1 "$(ls $FR/background_fit_summary_m*.csv | head -1)" > $FR/background_fit_summary.csv 2>/dev/null
for m in $MASS30; do tail -n +2 $FR/background_fit_summary_m$m.csv >> $FR/background_fit_summary.csv 2>/dev/null; done
python3 - "$FR" "$T" <<'PY'
import os, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
fr, t0 = sys.argv[1], float(sys.argv[2])
bad = []
for m in range(1, 31):
    p = "%s/%d/CMS-HGG_mva_13p6TeV_multipdf.root" % (fr, m)
    ws = "/eos/home-p/pelai/HZa/root_MVAcut/data/mA_M%d/ws/run3.root" % m
    if not os.path.exists(p) or os.path.getmtime(p) < os.path.getmtime(ws):
        bad.append((m, "missing or older than its input workspace")); continue
    if os.path.getsize(p) < 1024:
        bad.append((m, "size %d B (truncated)" % os.path.getsize(p))); continue
    f = ROOT.TFile.Open(p); n = -1
    for k in f.GetListOfKeys():
        o = k.ReadObj()
        if o.InheritsFrom("RooWorkspace"):
            it = o.allPdfs().createIterator(); pdf = it.Next()
            while pdf:
                if pdf.InheritsFrom("RooMultiPdf"): n = max(n, pdf.getNumPdfs())
                pdf = it.Next()
    f.Close()
    if n <= 0: bad.append((m, "getNumPdfs=%d" % n))
    else: print("  mA%-2d multipdf OK, %d pdfs" % (m, n))
for b in bad: print("  mA%-2d BAD: %s" % b)
sys.exit(1 if bad else 0)
PY
[ $? -eq 0 ] || die "background envelopes incomplete (logs: bkg_m*.log)"
fi

echo
if run 6; then
echo "=============== [6/9] signal model (fTest -> photon syst -> signalFit -> plotter -> effSigma) ==============="
T=$(date +%s)
cd $F/shellScripts
for s in 1_runjob_sig_fTest 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter 5_runjob_sig_plotEffSigma; do
  echo "  [$(date '+%H:%M:%S')] $s"
  bash $F/shellScripts/sig_sys/$s.sh > $L/sig_$s.log 2>&1
done
bad=0
for m in $MASS14; do for y in $YEARS; do for ch in ele mu; do
  fresh $F/Signal/outdir_$ch/signalFit/output/${m}_CMS-HGG_sigfit_${y}_${ch}_Hm125.root $T || { echo "  sigfit mA$m $y $ch stale/missing"; bad=1; }
done; done; done
[ $bad -eq 0 ] || die "signal fits incomplete (logs: sig_*.log)"
echo "  signal fits: 140/140 fresh"
fi

echo
if run 7; then
echo "=============== [7/9] datacards ==============="
T=$(date +%s)
check_effjson || die "interpolation efficiency JSON is stale -- datacards would scale 16 mass points with old curves"
cd $F/Datacard
sh 1_runjob_gen_datacard_makeYields.sh  > $L/dc_1.log 2>&1
sh 2_runjob_gen_datacard_makeDatacard.sh > $L/dc_2.log 2>&1
sh 3_rysn_datacard.sh                    > $L/dc_3.log 2>&1
bad=0; for m in $MASS30; do fresh $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt $T || { echo "  datacard mA$m stale/missing"; bad=1; }; done
[ $bad -eq 0 ] || die "datacards incomplete (logs: dc_*.log)"
echo "  datacards: 30/30 fresh"
# makeDatacard.py compared the string --mass_ALP to an int list until 2026-09-24, so no
# interpolated point ever got its efficiency rateParam. Interpolated: 5 eras x 2 channels.
for m in $MASS30; do
  n=$(grep -c 'rateParam' $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt)
  case " $MASS14 " in *" $m "*) want=0 ;; *) want=10 ;; esac
  [ "$n" -eq "$want" ] || { echo "  datacard mA$m has $n rateParam lines, expected $want"; bad=1; }
done
[ $bad -eq 0 ] || die "interpolation rateParams missing or misplaced"
echo "  rateParams: 16 interpolated x 10, 14 anchors x 0"
fi

echo
if run 8; then
echo "=============== [8/9] blind expected limits ==============="
T=$(date +%s)
cd $F/Combine
sh 2_text2ws.sh    > $L/combine_t2w.log 2>&1
sh 1_makeLimits.sh > $L/combine_limits.log 2>&1
bad=0; for m in $MASS30; do fresh $F/Combine/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || { echo "  limit mA$m stale/missing"; bad=1; }; done
[ $bad -eq 0 ] || die "limits incomplete (logs: combine_*.log)"
python3 - "$F/Combine/output_combine_results" "$SNAP/output_combine_results" <<'PY'
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def med(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m))
    if not f or f.IsZombie(): return None
    t = f.Get("limit"); v = None
    for e in t:
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
new, old = sys.argv[1], sys.argv[2]
print("  %4s %12s %12s %9s" % ("mA", "median_old", "median_new", "new/old"))
for m in range(1, 31):
    o, n = med(old, m), med(new, m)
    r = "%9.3f" % (n / o) if (o and n) else "%9s" % "-"
    print("  %4d %12s %12s %s" % (m, "%.4g" % o if o else "-", "%.4g" % n if n else "-", r))
PY
fi

echo
if run 9; then
echo "=============== [9/9] limit plots ==============="
T=$(date +%s)
cd $F/Plots && sh 2_runLimitsPlot.sh > $L/limit_plots.log 2>&1
n=$(find $F/Plots -name '*.pdf' -newermt "@$T" 2>/dev/null | wc -l)
echo "  fresh limit-plot pdfs: $n"
[ "$n" -gt 0 ] || die "no limit plots produced (limit_plots.log)"
fi

echo
echo "VERDICT: PASSED"
echo "Not run: impacts/bias (condor, fire-and-forget, launched separately) and Update AN (needs the user's go-ahead to push)."
echo "=============== FLASHGG CHAIN DONE $(date '+%F %T') ==============="
