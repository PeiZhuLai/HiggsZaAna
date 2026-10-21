#!/usr/bin/env bash
# CMS_hza_mva_reweight made symmetric: kappa = exp(max(ln kappa_up, -ln kappa_down)) on both sides, all 30 mass points.
# At m_a = 1-3 every alternative reweighting lowers the efficiency, so the envelope gave kappa_up = 1.000 (one-sided).
# Sandbox (studies/l3_review_20260910/mva_rw_symm_test): expected limit +0.3% / +0.2% / 0% at m_a = 1 / 2 / 3.
# The lnN enters through per-event constants weight_mva_reweight_Up/Down written by apply_bdt_sig.py, so the signal
# MVA-cut trees and Tree2WS workspaces of the 14 generated masses are redone; only those two columns may change.
#   0  wait for deploy_mA20_sigfit_powell (logs/chain.DONE)
#   1  backups: EOS root_MVAcut/sig -> root_MVAcut_before_rwsym (file by file, checksum), AFS products -> backup_prod/,
#      manifest of the production directories
#   2  symmetric lnN JSON (make_s1_lnn_json.py --symmetric, same CSV sources) installed as the production default
#   3  apply_bdt_sig.py, 14 masses x 5 eras; every tree: entries and sum of weights equal to the backup,
#      weight_mva_reweight_Up/Down equal to the new JSON
#   4  Tree2WS signal only; every dataset: entries and sum of weights equal to the backup
#   5  datacards (30 points, both cards); the only changed line is CMS_hza_mva_reweight, = 1/kappa_down on the up side
#   6  text2workspace + expected/observed limits + significance, 30 points x 2 cards
#   7  pseudo-data closure limits/significance/fitdiag, 30 points (datacards are built from the real ones)
#   8  impacts, 30 points x 2 cards (condor; stuck jobs are re-run locally)
#   9  limit plots, regression list
# Signal models and background envelopes are not touched; bias is not redone (the lnN only rescales the signal).
# Markers: logs/chain.DONE or logs/chain.FAILED; heartbeat logs/chain_heartbeat.txt. At most 3 local processes.
# Restart from a stage: START_STAGE=<n> bash chain_prod_rwsym.sh
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_mva_rw_symmetric
PREV=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_mA20_sigfit_powell/logs
S1D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/s1_altsideband_20260927
LNN=$S1D/results/s1_lnn_scheme2_sigrw_dyveto.json
NEWJ=$S1D/results/s1_lnn_scheme2_sigrw_dyveto_symm.json
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
DC=$F/Datacard/output_Datacard_leptons
PC=$F/pseudodata_closure
MV=/eos/home-p/pelai/HZa/root_MVAcut
MVB=/eos/home-p/pelai/HZa/root_MVAcut_before_rwsym
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
CPY=$ALP_PY
SIGM="1 2 3 4 5 6 7 8 9 10 15 20 25 30"
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
TAG=before_rwsym
L=$D/logs; DR=$D/drivers; BK=$D/backup_prod
START=${START_STAGE:-0}
mkdir -p $L $DR $BK
run() { [ "$START" -le "$1" ]; }
hb() { echo "$(date '+%F %T') $*" | tee -a $L/chain_heartbeat.txt; }
die() { hb "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/chain.FAILED; exit 1; }
source /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/eos_fresh.sh
rm -f $L/chain.DONE $L/chain.FAILED
AFSDIRS="Signal/outdir_ele Signal/outdir_mu Datacard/yields_leptons Datacard/output_Datacard_leptons Combine/root_t2w_xsec
 Combine/root_t2w_wilson Combine/output_combine_results_xsec Combine/output_combine_results_wilson
 Combine/output_combine_results_observed Combine/output_combine_results_observed_wilson Combine/output_significance
 Combine/output_impacts Combine/output_impacts_xsec Background/ALP_BkgModel_ReReco/fit_results_run3 Plots/plot_limits"
manifest() { ( cd $F && find $AFSDIRS -type f ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 md5sum ); }

if run 0; then
hb "[0] wait for the m_a = 20 deployment"
while [ ! -e $PREV/chain.DONE ]; do
  [ -e $PREV/chain.FAILED ] && die "m_a = 20 chain failed: $(cat $PREV/chain.FAILED)"
  sleep 300
done
hb "    m_a = 20 chain done ($(cat $PREV/chain.DONE))"
fi

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }   # the Datacard drivers call cmsenv
export -f cmsenv
cmsenv
export PYTHONPATH="$F/tools:$F/Signal/tools"   # same environment as deploy_w123_20261017/chain_prod_m123.sh
set -u

if run 1; then
hb "[1] backups + manifest"
( cd $MV && find sig -maxdepth 3 -type f \( -name 'output_*.root' -o -path '*/ws_Tree2WS/*.root' \) | grep -E "^sig/mA_M($(echo $SIGM | tr ' ' '|'))/" ) > $L/mvacut_backup_list.txt
n=$(wc -l < $L/mvacut_backup_list.txt); [ $n -eq 210 ] || die "expected 210 signal files (14 x 5 trees + 14 x 10 ws), found $n"
while read -r f; do
  x=/eos/user/p/pelai/HZa/root_MVAcut/$f; y=/eos/user/p/pelai/HZa/root_MVAcut_before_rwsym/$f
  ck=$(xrdfs eosuser.cern.ch query checksum $x | awk '{print $2}')
  [ -n "$ck" ] || die "no checksum for $f"
  [ "$(xrdfs eosuser.cern.ch query checksum $y 2>/dev/null | awk '{print $2}')" = "$ck" ] && continue
  xrdfs eosuser.cern.ch mkdir -p $(dirname $y)
  xrdcp -f -s --cksum adler32:$ck root://eosuser.cern.ch/$x root://eosuser.cern.ch/$y || die "backup of $f"
done < $L/mvacut_backup_list.txt
hb "    EOS backup: $n files in $MVB"
[ -s $BK/datacards.tar ] || tar -C $F/Datacard -cf $BK/datacards.tar yields_leptons output_Datacard_leptons || die "backup datacards"
[ -s $BK/combine.tar ] || tar -C $C -cf $BK/combine.tar root_t2w_xsec root_t2w_wilson output_combine_results_xsec output_combine_results_wilson \
  output_combine_results_observed output_combine_results_observed_wilson output_significance output_impacts output_impacts_xsec || die "backup combine"
[ -s $BK/closure.tar ] || tar -C $PC -cf $BK/closure.tar closure_results.json plots datacards root_t2w output_limits output_significance output_fitdiag || die "backup closure"
[ -s $BK/plot_limits.tar ] || tar -C $F/Plots -cf $BK/plot_limits.tar plot_limits || die "backup plot_limits"
mkdir -p $BK/dc_before && cp -p $DC/*_pruned_datacard_leptons*.txt $BK/dc_before/
cp -p $PC/closure_results.json $L/closure_results_before.json
[ -s $L/manifest.before ] || manifest > $L/manifest.before
hb "    manifest: $(wc -l < $L/manifest.before) files"
fi

if run 2; then
hb "[2] symmetric lnN JSON"
SRC=$($ALP_PY -c "import json; print(' '.join(json.load(open('$LNN'))['source']))")
( cd $S1D/results && env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY ../make_s1_lnn_json.py --reference nominal --symmetric $NEWJ $SRC ) > $L/s1_json.log 2>&1 || die "make_s1_lnn_json failed"
$ALP_PY - $NEWJ $LNN <<'PY' | tee -a $L/chain_heartbeat.txt || die "symmetric JSON inconsistent with the current one"
import json, math, sys
new, old = (json.load(open(p)) for p in sys.argv[1:3])
assert new["symmetric"] is True and new["reference"] == old["reference"] == "nominal"
n, o = new["values"], old["values"]; bad = sorted(n) != sorted(o)
for m in sorted(o, key=int):
    lk = max(math.log(o[m]["up"]), -math.log(o[m]["down"]))
    ok = abs(n[m]["up"] - math.exp(lk)) < 2e-5 and abs(n[m]["down"] - math.exp(-lk)) < 2e-5 and n[m]["mvacut"] == o[m]["mvacut"]
    bad |= not ok
    print("    mA%-2s %.3f/%.3f -> %.3f/%.3f%s" % (m, o[m]["down"], o[m]["up"], n[m]["down"], n[m]["up"], "" if ok else "  WRONG"))
sys.exit(1 if bad else 0)
PY
[ -e $LNN.$TAG ] || cp -p $LNN $LNN.$TAG
cp $NEWJ $LNN
fi
export HZA_MVA_REWEIGHT_LNN_JSON=$LNN

if run 3; then
hb "[3] apply_bdt_sig, 14 masses x 5 eras"
for m in $SIGM; do eos_mount_sync /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M$m || die "scored mA_M$m not visible"; done
T=$(date +%s)
sig_one() { local s=$1 y=$2 a; for a in 1 2 3; do
    python3 $F/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples $s --years $y > $L/mva_sig_${s}_$y.log 2>&1 && return 0; sleep 60; done; return 1; }
export -f sig_one; export F L
for m in $SIGM; do for y in $YEARS; do echo mA_M$m $y; done; done | xargs -P3 -L1 bash -c 'sig_one $0 $1 || echo "FAILED $0 $1"' > $L/mva_sig_status.txt
grep -q FAILED $L/mva_sig_status.txt && die "apply_bdt_sig failed: $(grep FAILED $L/mva_sig_status.txt | tr '\n' ' ')"
for m in $SIGM; do for y in $YEARS; do fresh $MV/sig/mA_M$m/output_$y.root $T || die "sig mA$m $y stale"; done; done
python3 - $MV $MVB $LNN "$SIGM" <<'PY' | tee $L/mva_compare.txt | tail -3 | tee -a $L/chain_heartbeat.txt; [ ${PIPESTATUS[0]} -eq 0 ] || die "MVA-cut trees changed beyond the reweight columns (mva_compare.txt)"
import json, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
P, B, J, masses = sys.argv[1], sys.argv[2], json.load(open(sys.argv[3]))["values"], sys.argv[4].split(); bad = 0
def trees(path):
    out = {}; f = ROOT.TFile.Open(path)
    def walk(d, pre):
        for k in d.GetListOfKeys():
            o = k.ReadObj()
            if o.InheritsFrom("TDirectory"): walk(o, pre + o.GetName() + "/")
            elif o.InheritsFrom("TTree"):
                sw = 0.0; ud = set()
                has = bool(o.GetBranch("weight_mva_reweight_Up"))
                for e in o:
                    sw += e.weight
                    if has: ud.add((round(e.weight_mva_reweight_Up, 6), round(e.weight_mva_reweight_Down, 6)))
                out[pre + o.GetName()] = (o.GetEntries(), round(sw, 6), ud)
    walk(f, ""); f.Close(); return out
for m in masses:
    for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
        rel = "sig/mA_M%s/output_%s.root" % (m, y)
        a, b = trees("%s/%s" % (P, rel)), trees("%s/%s" % (B, rel))
        want = {(round(J[m]["up"], 6), round(J[m]["down"], 6))}
        ok = sorted(a) == sorted(b) and all(a[t][:2] == b[t][:2] for t in a) and all(a[t][2] <= want for t in a)
        print("  %-30s %d trees %s" % (rel, len(a), "ok" if ok else "DIFFERS"))
        bad += not ok
print("    %d files compared, %d differ" % (len(masses) * 5, bad))
sys.exit(1 if bad else 0)
PY
fi

if run 4; then
hb "[4] Tree2WS signal"
sed -e 's/^max_parallel=6$/max_parallel=3/' $F/Trees2WS/run_tree2ws.sh > $DR/run_tree2ws_sig.sh
grep -q '^max_parallel=3$' $DR/run_tree2ws_sig.sh || die "tree2ws parallelism not reduced"
T=$(date +%s)
( cd $F/Trees2WS && TREE2WS_PART=sig bash $DR/run_tree2ws_sig.sh ) > $L/tree2ws.log 2>&1
n=$(find $MV/sig -path '*/ws_Tree2WS/*.root' -newermt "@$T" | wc -l); [ $n -eq 140 ] || die "only $n of 140 signal workspaces fresh (tree2ws.log)"
python3 - $MV $MVB "$SIGM" <<'PY' | tee $L/ws_compare.txt | tail -2 | tee -a $L/chain_heartbeat.txt; [ ${PIPESTATUS[0]} -eq 0 ] || die "signal workspaces changed (ws_compare.txt)"
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
P, B, masses = sys.argv[1], sys.argv[2], sys.argv[3].split(); bad = 0
def summ(p):
    f = ROOT.TFile.Open(p); w = f.Get("CMS_hza_workspace/CMS_hza_workspace")
    out = {d.GetName(): (d.numEntries(), round(d.sumEntries(), 6)) for d in w.allData()}; f.Close(); return out
for m in masses:
    for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
        for l in ("ele", "mu"):
            rel = "sig/mA_M%s/ws_Tree2WS/ws_%s_%s.root" % (m, l, y)
            ok = summ("%s/%s" % (P, rel)) == summ("%s/%s" % (B, rel)); bad += not ok
            if not ok: print("  DIFFERS", rel)
print("    140 signal workspaces compared, %d differ" % bad)
sys.exit(1 if bad else 0)
PY
fi

if run 5; then
hb "[5] datacards, 30 points"
T=$(date +%s)
( cd $F/Datacard && sh 1_runjob_gen_datacard_makeYields.sh ) > $L/dc_1.log 2>&1
( cd $F/Datacard && sh 2_runjob_gen_datacard_makeDatacard.sh ) > $L/dc_2.log 2>&1
for m in $(seq 1 30); do for s in "" _xsec; do fresh $DC/${m}_pruned_datacard_leptons$s.txt $T || die "datacard mA$m$s stale (dc_*.log)"; done; done
python3 - $DC $BK/dc_before $LNN <<'PY' | tee $L/dc_compare.txt | tail -1 | tee -a $L/chain_heartbeat.txt; [ ${PIPESTATUS[0]} -eq 0 ] || die "datacards changed beyond CMS_hza_mva_reweight (dc_compare.txt)"
import json, sys
N, O, J = sys.argv[1], sys.argv[2], json.load(open(sys.argv[3]))["values"]
anchors = sorted(int(m) for m in J); bad = 0
for m in range(1, 31):
    a = min(anchors, key=lambda x: (abs(x - m), x)); k = J[str(a)]
    for s in ("", "_xsec"):
        n = open("%s/%d_pruned_datacard_leptons%s.txt" % (N, m, s)).read().split("\n")
        o = open("%s/%d_pruned_datacard_leptons%s.txt" % (O, m, s)).read().split("\n")
        diff = [(x, y) for x, y in zip(n, o) if x != y] + ([("length", "")] if len(n) != len(o) else [])
        ok = len(diff) == 1 and diff[0][0].startswith("CMS_hza_mva_reweight")
        if ok:
            vals = set(v for v in diff[0][0].split()[2:] if v != "-")
            # writeToDatacard.py writes a single kappa_up when kappa_down * kappa_up ~ 1
            ok = vals in ({"%.3f/%.3f" % (k["down"], k["up"])}, {"%.3f" % k["up"]})
        print("  mA%-2d%-5s %s %s" % (m, s, "ok" if ok else "DIFFERS", diff[0][0][:80] if diff else "no change"))
        bad += not ok
print("    60 datacards: only CMS_hza_mva_reweight changed in %d, problems in %d" % (60 - bad, bad))
sys.exit(1 if bad else 0)
PY
fi

lim_one() {   # same as deploy_mA20_sigfit_powell/chain_prod_mA20.sh
  local m=$1 k=$2 sfx="" W
  [ $k = xsec ] && sfx=_xsec
  W=$D/work_combine/${k}_$m; rm -rf $W; mkdir -p $W
  ( cd $C && python3 RunText2Workspace.py --batch local --queue workday --mA $m --card $k ) > $L/t2w_${k}_m$m.log 2>&1
  ( cd $W && combine -M AsymptoticLimits $C/root_t2w_$k/${m}_Datacard_leptons.root --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/limit_${k}_m$m.log 2>&1
    mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_$k/
    combine -M AsymptoticLimits $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons$sfx.txt --cminDefaultMinimizerStrategy 0 -m 125.38 \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/observed_${k}_m$m.log 2>&1
    if [ $k = xsec ]; then mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_observed/
    else mv higgsCombine$m.AsymptoticLimits.mH125.38.root $C/output_combine_results_observed_wilson/; fi
    if [ $k = xsec ]; then
      combine -M Significance $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons_xsec.txt --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -5 --rMax 5 -n $m > $L/significance_m$m.log 2>&1
      mv higgsCombine$m.Significance.mH125.38.root $C/output_significance/
    fi )
}
export -f lim_one; export C D L F

if run 6; then
hb "[6] text2workspace + limits, 30 points x 2 cards"
T=$(date +%s)
for m in $(seq 1 30); do printf "$m xsec\n$m wilson\n"; done | xargs -P3 -L1 bash -c 'lim_one $0 $1'
for m in $(seq 1 30); do
  for k in xsec wilson; do
    fresh $C/root_t2w_$k/${m}_Datacard_leptons.root $T && fresh $C/output_combine_results_$k/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "expected limit $k m$m"
  done
  fresh $C/output_combine_results_observed/higgsCombine$m.AsymptoticLimits.mH125.38.root $T && fresh $C/output_combine_results_observed_wilson/higgsCombine$m.AsymptoticLimits.mH125.38.root $T \
    && fresh $C/output_significance/higgsCombine$m.Significance.mH125.38.root $T || die "observed/significance m$m"
done
mkdir -p $D/old_extract && tar -C $D/old_extract -xf $BK/combine.tar output_combine_results_xsec output_combine_results_wilson
python3 - $D/old_extract $C <<'PY' | tee $L/limit_compare.txt | tail -1 | tee -a $L/chain_heartbeat.txt
import ROOT, sys
O, N = sys.argv[1], sys.argv[2]
def med(p):
    f = ROOT.TFile.Open(p); v = None
    for e in f.Get("limit"):
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
r = []
for k in ("xsec", "wilson"):
    for m in range(1, 31):
        n = "output_combine_results_%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (k, m)
        a, b = med("%s/%s" % (O, n)), med("%s/%s" % (N, n)); r.append(b / a - 1)
        print("  %-6s mA%-2d %.4f -> %.4f (%+.1f%%)" % (k, m, a, b, 100 * (b / a - 1)))
print("    expected median limit change over 60: min %+.1f%%, max %+.1f%%" % (100 * min(r), 100 * max(r)))
PY
fi

if run 7; then
hb "[7] pseudo-data closure, 30 points"
T=$(date +%s)
( export CLOSURE_FSRFIX=1; ALPMASS="$(seq -s ' ' 1 30)" NPROC=3 bash $PC/scripts/run_limits_pseudodata.sh ) > $L/closure_limits.log 2>&1
n=$(find $PC/output_limits -name '*.root' -newermt "@$T" | wc -l); [ $n -ge 30 ] || die "only $n fresh closure limit files"
( cd $PC/scripts && python3 extract_results.py ) > $L/closure_extract.log 2>&1 || die "extract_results"
( cd $PC/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin $CPY plot_closure.py ) > $L/closure_plot.log 2>&1 || die "plot_closure"
python3 - $L/closure_results_before.json $PC/closure_results.json <<'PY' | tee $L/closure_compare.txt | tail -1 | tee -a $L/chain_heartbeat.txt
import json, sys
o, n = (json.load(open(p)) for p in sys.argv[1:3])
def rows(d):
    d = d.get("results", d)
    return {int(k): v for k, v in (d.items() if isinstance(d, dict) else [(r["mA"], r) for r in d])}
o, n = rows(o), rows(n)
for m in sorted(n):
    a, b = o[m]["limit"], n[m]["limit"]
    print("  mA%-2d obs %.4f -> %.4f  exp %.4f -> %.4f  rhat %+.4f -> %+.4f" % (m, a["obs"], b["obs"], a["exp_med"], b["exp_med"], o[m]["rhat"]["r"], n[m]["rhat"]["r"]))
print("    closure: %d points recomputed" % len(n))
PY
fi

if run 8; then
hb "[8] impacts, 30 points x 2 cards (condor)"
T0=$(date +%s)
out=$( (IMPACT_CARD=wilson bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh; IMPACT_CARD=xsec bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh) 2>&1 )
echo "$out" > $L/impact_submit.log
ICL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$out" | sort -u | tr '\n' ' ')
[ -n "${ICL// /}" ] || die "impact submission failed (impact_submit.log)"
hb "    impact clusters $ICL"
start=$(date +%s)
while true; do
  q=$(condor_q $ICL -af ProcId JobStatus 2>/dev/null) || { sleep 300; continue; }
  idle=$(awk '$2==1' <<<"$q" | wc -l); runn=$(awk '$2==2' <<<"$q" | wc -l); held=$(awk '$2==5' <<<"$q" | wc -l)
  hb "    impacts idle=$idle running=$runn held=$held"
  [ $((idle + runn + held)) -eq 0 ] && break
  if [ $runn -eq 0 ] && [ $(( $(date +%s) - start )) -gt 7200 ]; then
    hb "    impact jobs idle/held for 2 h: remove them and run the missing ones locally (one at a time)"
    condor_rm $ICL > /dev/null 2>&1; sleep 30
    for m in $(seq 1 30); do for k in wilson xsec; do
      j=$C/$( [ $k = xsec ] && echo output_impacts_xsec || echo output_impacts )/${m}_impacts.json
      fresh $j $T0 && continue
      IMPACT_CARD=$k bash $F/shellScripts/impact/Condor/runjob_expectedImpact.sh $m $F ${CMSSW_SRC%/src} > $L/impact_local_${k}_m$m.log 2>&1 || hb "    local impact $k m$m rc!=0"
    done; done
    break
  fi
  sleep 1800
done
python3 - $C $T0 <<'PY' | tee $L/impact_check.txt | tail -4 | tee -a $L/chain_heartbeat.txt || die "impacts missing or incomplete (impact_check.txt)"
import json, os, sys
C, T0 = sys.argv[1], float(sys.argv[2]); bad = 0
for k, d in (("wilson", "output_impacts"), ("xsec", "output_impacts_xsec")):
    for m in range(1, 31):
        p = "%s/%s/%d_impacts.json" % (C, d, m)
        if not (os.path.exists(p) and os.path.getmtime(p) >= T0): print("  %s mA%d MISSING/stale" % (k, m)); bad += 1; continue
        j = json.load(open(p)); q = [x for x in j["params"] if x["name"] == "CMS_hza_mva_reweight"]
        if not q: print("  %s mA%d no CMS_hza_mva_reweight" % (k, m)); bad += 1; continue
        lo, hi = q[0]["r"][0] - q[0]["r"][1], q[0]["r"][2] - q[0]["r"][1]
        print("  %-6s mA%-2d mva_reweight impact on r: %+.4f / %+.4f" % (k, m, lo, hi))
print("    impacts: %d of 60 missing" % bad)
sys.exit(bad)
PY
fi

if run 9; then
hb "[9] limit plots + regression"
T=$(date +%s)
bash $F/Plots/2_runLimitsPlot.sh > $L/limitsplot.log 2>&1 || hb "    2_runLimitsPlot.sh rc!=0 (limitsplot.log)"
n=$(find $F/Plots/plot_limits -newermt "@$T" -type f | wc -l); [ $n -ge 1 ] || die "no fresh limit plots"
manifest > $L/manifest.after
python3 - $L/manifest.before $L/manifest.after <<'PY' | tee $L/regression.txt
import re, sys
def rd(p): return {l[34:].strip(): l[:32] for l in open(p)}
a, b = rd(sys.argv[1]), rd(sys.argv[2])
changed = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
pat = re.compile(r"^Datacard/|^Combine/(root_t2w|output_combine_results|output_significance|output_impacts)|^Plots/plot_limits/")
other = [k for k in changed if not pat.search(k)]
print("    changed files: %d, outside datacards/combine/limit plots: %d" % (len(changed), len(other)))
for k in other[:80]: print("      ", k)
sys.exit(1 if other else 0)
PY
[ ${PIPESTATUS[0]} -eq 0 ] || hb "    WARNING: files outside the expected set changed (regression.txt), check by hand"
fi

hb "CHAIN DONE"
date '+%F %T' > $L/chain.DONE
