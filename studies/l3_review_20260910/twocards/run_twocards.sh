#!/usr/bin/env bash
# Production switch to two datacards per m_a (user decision 2026-10-15), run once on the current production cards:
#   1 Datacard/make_xsec_datacards.sh            -> <m>_pruned_datacard_leptons_xsec.txt (Wilson card minus 3 theory lines)
#   2 Combine/2_text2ws.sh                         -> root_t2w_{xsec,wilson}, root_t2w -> root_t2w_xsec (old dir kept as *_before_twocards_*)
#   3 Combine/1_makeLimits.sh                      -> output_combine_results_{xsec,wilson}, output_combine_results -> _xsec
#   4 Combine/1b_makeLimits_observed.sh            -> output_combine_results_observed (xsec) and _observed_wilson (old dir kept)
#   5 Plots/2_runLimitsPlot.sh
#   6 comparison table -> twocards/limits_xsec_vs_wilson.txt (Wilson card must reproduce the previous limits exactly)
# Status twocards/logs/status; DONE/STOP line in the varsearch registry.
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/twocards
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/logs/PENDING_RESULTS.txt
ST=$W/logs/status; TAG=$(date +%Y%m%d_%H%M)
stat() { echo "$(date '+%F %T') $*" > $ST; }
stop() { stat "STOP $*"; echo "$(date '+%F %H:%M') STOP  production two datacards  $*" >> $REG; exit 1; }
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && eval $(scramv1 runtime -sh)
cmsenv() { :; }; export -f cmsenv
T0=${T0:-$(date +%s)}; START=${START:-1}; run() { [ "$START" -le "$1" ]; }
fresh() { [ -s "$1" ] && [ "$(command stat -c %Y "$1")" -ge "$T0" ]; }   # command: stat() above is the status writer
if run 1; then
stat "1 xsec cards"
bash $F/Datacard/make_xsec_datacards.sh > $W/logs/1_xsec_cards.log 2>&1 || stop "make_xsec_datacards.sh ($W/logs/1_xsec_cards.log)"
fi
if run 2; then
stat "2 text2workspace (60)"
cd $F/Combine && bash 2_text2ws.sh > $W/logs/2_text2ws.log 2>&1
fi
for c in xsec wilson; do for m in $(seq 1 30); do fresh $F/Combine/root_t2w_$c/${m}_Datacard_leptons.root || stop "t2w $c mA$m missing"; done; done
[ "$(readlink $F/Combine/root_t2w)" = root_t2w_xsec ] || stop "root_t2w is not a link to root_t2w_xsec"
stat "3 expected limits (60)"
OLD=$(ls -d $F/Combine/output_combine_results_before_twocards_* 2>/dev/null | tail -1)
cd $F/Combine && bash 1_makeLimits.sh > $W/logs/3_limits.log 2>&1
[ -z "$OLD" ] && OLD=$(ls -d $F/Combine/output_combine_results_before_twocards_* 2>/dev/null | tail -1)
for c in xsec wilson; do for m in $(seq 1 30); do fresh $F/Combine/output_combine_results_$c/higgsCombine$m.AsymptoticLimits.mH125.38.root || stop "limit $c mA$m missing"; done; done
stat "4 observed limits (60)"
[ -d $F/Combine/output_combine_results_observed ] && cp -rp $F/Combine/output_combine_results_observed $F/Combine/output_combine_results_observed_before_twocards_$TAG
cd $F/Combine && bash 1b_makeLimits_observed.sh > $W/logs/4_observed.log 2>&1
for d in output_combine_results_observed output_combine_results_observed_wilson; do for m in $(seq 1 30); do fresh $F/Combine/$d/higgsCombine$m.AsymptoticLimits.mH125.38.root || stop "observed $d mA$m missing"; done; done
stat "5 plots"
cd $F/Plots && bash 2_runLimitsPlot.sh > $W/logs/5_plots.log 2>&1 || stop "2_runLimitsPlot.sh ($W/logs/5_plots.log)"
stat "6 comparison"
python3 - $F/Combine $OLD > $W/limits_xsec_vs_wilson.txt 2>&1 <<'PY' || stop "comparison failed"
import sys, ROOT
ROOT.gROOT.SetBatch(True)
C, old = sys.argv[1], sys.argv[2]
def q(p):
    f = ROOT.TFile.Open(p); out = {round(e.quantileExpected, 3): e.limit for e in f.Get("limit")}; f.Close(); return out
print("blind expected limit on r (signal at 0.1 pb), median [-1s, +1s]; previous = production before the change (Wilson-type card)")
print("  mA |      previous         |      wilson           |       xsec            | wilson/prev  xsec/wilson")
bad = 0
for m in range(1, 31):
    n = "higgsCombine%d.AsymptoticLimits.mH125.38.root" % m
    p, w, x = q("%s/%s" % (old, n)), q("%s/output_combine_results_wilson/%s" % (C, n)), q("%s/output_combine_results_xsec/%s" % (C, n))
    f = lambda d: "%.4f [%.4f,%.4f]" % (d[0.5], d[0.16], d[0.84])
    bad += abs(w[0.5] / p[0.5] - 1) > 1e-6
    print("  %2d | %s | %s | %s | %10.6f %10.3f" % (m, f(p), f(w), f(x), w[0.5] / p[0.5], x[0.5] / w[0.5]))
print("Wilson card reproduces the previous limits at all 30 points" if not bad else "WARNING: %d points where the Wilson card differs from the previous limit" % bad)
PY
stat "DONE"
echo "$(date '+%F %H:%M') DONE  production two datacards (t2w, expected + observed limits, plots)  -> $W/limits_xsec_vs_wilson.txt" >> $REG
