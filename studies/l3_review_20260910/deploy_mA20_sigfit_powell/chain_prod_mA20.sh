#!/usr/bin/env bash
# m_a = 20, 2023postBPix, electron signal model refitted with --minimizerMethod Powell. The production fit (scipy TNC) never
# left its starting point (dm 0.1, sigma 1/2/3/4, frac 0.5/0.45/0.4; chi2 153 -> 153 on 48 bins, sigma_eff 1.65 GeV vs
# 2.39-2.60 in the other eras). The Powell refit (sandbox studies/l3_review_20260910/sigfit_mA20_postBPix_ele/nRV4_Powell)
# gives chi2 4.6 and sigma_eff 2.54 GeV; only the shape parameters change, the normalization is identical. Of the 140
# signal fits this is the only one stuck at its starting values. The m_a = 20 model is used directly by m_a = 18-22.
#   1  backups of every product that is redone, manifest of the production directories
#   2  signalFit m_a = 20 2023postBPix ele with Powell; the workspace must equal the sandbox refit value by value
#   3  RunPlotter (that era + Run 3 combination) and plot_effisigma.py
#   4  text2workspace + expected/observed limits + significance m_a = 18-22, both cards
#   5  pseudo-data closure m_a = 18-22 (limits/significance/fitdiag only; envelope and pseudo-data are signal-independent)
#   6  bias (10 x 1000 toys, condor) and impacts (both cards, condor) m_a = 18-22; bias plots of all 30 points redrawn
#      (axis titles fixed in plot_bias.py: x = (mu~ - mu)/sigma as computed, y = generated mu)
#   7  limit plots (Plots/2_runLimitsPlot.sh)
#   8  regression: list of every changed production file
# Markers: logs/chain.DONE or logs/chain.FAILED; heartbeat logs/chain_heartbeat.txt. At most 3 local processes.
# Restart from a stage: START_STAGE=<n> bash chain_prod_mA20.sh
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_mA20_sigfit_powell
SB=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/sigfit_mA20_postBPix_ele
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
HS=$BN/bias_outputs_highstat
PC=$F/pseudodata_closure
SIGO=$F/Signal/outdir_ele/signalFit
WS=20_CMS-HGG_sigfit_2023postBPix_ele_Hm125.root
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
MASSES="18 19 20 21 22"
TAG=before_mA20powell
L=$D/logs; BK=$D/backup_prod
START=${START_STAGE:-1}
mkdir -p $L $BK
run() { [ "$START" -le "$1" ]; }
hb() { echo "$(date '+%F %T') $*" | tee -a $L/chain_heartbeat.txt; }
die() { hb "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/chain.FAILED; exit 1; }
source /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/eos_fresh.sh
rm -f $L/chain.DONE $L/chain.FAILED
AFSDIRS="Signal/outdir_ele Signal/outdir_mu Datacard/yields_leptons Datacard/output_Datacard_leptons Combine/root_t2w_xsec
 Combine/root_t2w_wilson Combine/output_combine_results_xsec Combine/output_combine_results_wilson
 Combine/output_combine_results_observed Combine/output_combine_results_observed_wilson Combine/output_significance
 Combine/output_impacts Combine/output_impacts_xsec Background/ALP_BkgModel_ReReco/fit_results_run3
 Combine/Checks/Bias_nominal/bias_outputs_highstat Combine/Checks/Bias_nominal/plots_bias Plots/plot_limits"
manifest() { ( cd $F && find $AFSDIRS -type f ! -path '*/__pycache__/*' ! -path '*/bias_outputs_highstat_variants/*' -print0 | sort -z | xargs -0 md5sum ); }

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_SRC && eval $(scramv1 runtime -sh)
source $F/setup.sh > /dev/null 2>&1
set -u
export PYTHONPATH=$PYTHONPATH:$F/tools:$F/Signal/tools

if run 1; then
hb "[1] backups + manifest"
[ -s $BK/signal_mA20_ele.tar ] || ( cd $F/Signal/outdir_ele && tar -cf $BK/signal_mA20_ele.tar signalFit/output/$WS \
  $(cd $F/Signal/outdir_ele && ls signalFit/Plots | grep -E '(^|_)20_|effSigmaVmA' | sed 's|^|signalFit/Plots/|') ) || die "backup signal"
[ -s $BK/combine_m18-22.tar ] || ( cd $C && tar -cf $BK/combine_m18-22.tar $(for m in $MASSES; do
  ls -d root_t2w_xsec/${m}_* root_t2w_wilson/${m}_* output_combine_results_xsec/higgsCombine$m.* output_combine_results_wilson/higgsCombine$m.* \
     output_combine_results_observed/higgsCombine$m.* output_combine_results_observed_wilson/higgsCombine$m.* \
     output_significance/higgsCombine$m.* output_impacts/${m}_* output_impacts_xsec/${m}_* 2>/dev/null; done) ) || die "backup combine outputs"
[ -s $BK/closure_m18-22.tar ] || ( cd $PC && tar -cf $BK/closure_m18-22.tar closure_results.json plots $(for m in $MASSES; do
  ls -d datacards/${m}_* root_t2w/${m}_* output_limits/*PS$m.* output_significance/*PS$m.* output_fitdiag/*PS$m.* 2>/dev/null; done) ) || die "backup closure"
cp -p $PC/closure_results.json $L/closure_results_before.json
[ -s $BK/plots_bias.tar ] || tar -C $BN -cf $BK/plots_bias.tar plots_bias || die "backup plots_bias"
[ -s $BK/plot_limits.tar ] || tar -C $F/Plots -cf $BK/plot_limits.tar plot_limits || die "backup plot_limits"
[ -s $L/manifest.before ] || manifest > $L/manifest.before
hb "    manifest: $(wc -l < $L/manifest.before) files"
fi

if run 2; then
hb "[2] signalFit m_a = 20 2023postBPix ele, Powell"
T=$(date +%s)
( cd $CMSSW_SRC && python3 $F/Signal/scripts/signalFit.py --mass_ALP 20 --year 2023postBPix --channel ele \
    --inputWSDir /eos/home-p/pelai/HZa/root_MVAcut/sig/mA_M20/ws_Tree2WS --doSystematics --minimizerMethod Powell ) > $L/signalFit.log 2>&1 || die "signalFit rc!=0"
fresh $SIGO/output/$WS $T || die "no fresh workspace"
python3 $SB/cmp_norm.py $SIGO/output/$WS $SB/sandbox/nRV4_Powell/src/flashggFinalFit/Signal/outdir_ele/signalFit/output/$WS > $L/cmp_sandbox.txt 2>&1
hb "    vs sandbox: $(grep 'differing' $L/cmp_sandbox.txt)"
grep -q 'differing 0$' $L/cmp_sandbox.txt || die "production refit differs from the sandbox (cmp_sandbox.txt)"
fi

if run 3; then
hb "[3] RunPlotter + effSigma summary"
T=$(date +%s)
( cd $CMSSW_SRC && python3 $F/Signal/RunPlotter.py --mass_ALP 20 --years 2023postBPix --channel ele \
  && python3 $F/Signal/RunPlotter.py --mass_ALP 20 --channel ele --years '2022preEE,2022postEE,2023preBPix,2023postBPix,2024' ) > $L/runplotter.log 2>&1 || die "RunPlotter rc!=0"
( cd $CMSSW_SRC && python3 $F/Signal/scripts/plot_effisigma.py ) > $L/plot_effisigma.log 2>&1 || die "plot_effisigma rc!=0"
fresh $SIGO/Plots/effSigma_20_2023postBPix_ele.json $T && fresh $SIGO/Plots/effSigmaVmA_ele.pdf $T || die "effSigma products stale"
hb "    effSigma 20 2023postBPix ele: $(tr -d '\n ' < $SIGO/Plots/effSigma_20_2023postBPix_ele.json)"
fi

lim_one() {   # same as deploy_mA1_low100bern4/chain_prod_mA1.sh
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

if run 4; then
hb "[4] text2workspace + limits m_a = $MASSES"
T=$(date +%s)
for m in $MASSES; do printf "$m xsec\n$m wilson\n"; done | xargs -P3 -L1 bash -c 'lim_one $0 $1'
for m in $MASSES; do for k in xsec wilson; do
  fresh $C/root_t2w_$k/${m}_Datacard_leptons.root $T && fresh $C/output_combine_results_$k/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "expected limit $k m$m"
done
  fresh $C/output_combine_results_observed/higgsCombine$m.AsymptoticLimits.mH125.38.root $T && fresh $C/output_combine_results_observed_wilson/higgsCombine$m.AsymptoticLimits.mH125.38.root $T \
    && fresh $C/output_significance/higgsCombine$m.Significance.mH125.38.root $T || die "observed/significance m$m"
  hb "    m$m expected xsec $(grep 'Expected 50.0%' $L/limit_xsec_m$m.log | awk '{print $NF}'), wilson $(grep 'Expected 50.0%' $L/limit_wilson_m$m.log | awk '{print $NF}')"
done
fi

if run 5; then
hb "[5] pseudo-data closure m_a = $MASSES"
T=$(date +%s)
( export CLOSURE_FSRFIX=1; ALPMASS="$MASSES" NPROC=3 bash $PC/scripts/run_limits_pseudodata.sh ) > $L/closure_limits.log 2>&1
n=$(find $PC/output_limits -name '*.root' -newermt "@$T" | wc -l); [ $n -ge 5 ] || die "only $n fresh closure limit files"
( cd $PC/scripts && python3 extract_results.py ) > $L/closure_extract.log 2>&1 || die "extract_results"
( cd $PC/scripts && env -i HOME=$HOME PATH=/usr/bin:/bin $CPY plot_closure.py ) > $L/closure_plot.log 2>&1 || die "plot_closure"
python3 - $L/closure_results_before.json $PC/closure_results.json <<'PY' | tee $L/closure_compare.txt | tee -a $L/chain_heartbeat.txt
import json, sys
o, n = (json.load(open(p)) for p in sys.argv[1:3])
def rows(d):
    d = d.get("results", d)
    return {int(k): v for k, v in (d.items() if isinstance(d, dict) else [(r["mA"], r) for r in d])}
o, n = rows(o), rows(n)
for m in sorted(n):
    if o.get(m) != n[m]:
        print("    mA%-2d before %s\n          after  %s" % (m, o.get(m), n[m]))
print("    closure changed points:", [m for m in sorted(n) if o.get(m) != n[m]])
PY
fi

if run 6; then
hb "[6] bias + impacts m_a = $MASSES (condor)"
mkdir -p $BN/bias_outputs_highstat_variants
for m in $MASSES; do
  [ -d $BN/bias_outputs_highstat_variants/mA_${m}_$TAG ] || mv $HS/mA_$m $BN/bias_outputs_highstat_variants/mA_${m}_$TAG || die "move old bias m$m"
  rm -rf $HS/mA_$m
done
T0=$(date +%s)
out=$( (cd $F && bash shellScripts/bias/Condor/subjob_bias_highstat.sh $MASSES) 2>&1 ); echo "$out" > $L/bias_submit.log
out2=$( (IMPACT_CARD=wilson bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $MASSES; IMPACT_CARD=xsec bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $MASSES) 2>&1 )
echo "$out2" > $L/impact_submit.log
BCL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$out" | sort -u | tr '\n' ' ')
ICL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$out2" | sort -u | tr '\n' ' ')
[ -n "${BCL// /}" ] || die "bias submission failed (bias_submit.log)"
[ -n "${ICL// /}" ] || die "impact submission failed (impact_submit.log)"
hb "    bias clusters $BCL, impact clusters $ICL"
start=$(date +%s); local_impacts=0
while true; do
  qb=$(condor_q $BCL -af ProcId JobStatus 2>/dev/null) || { sleep 300; continue; }
  qi=$(condor_q $ICL -af ProcId JobStatus 2>/dev/null) || { sleep 300; continue; }
  cnt() { awk -v s=$2 '$2==s' <<<"$1" | wc -l; }
  hb "    bias idle=$(cnt "$qb" 1) running=$(cnt "$qb" 2) held=$(cnt "$qb" 5) | impacts idle=$(cnt "$qi" 1) running=$(cnt "$qi" 2) held=$(cnt "$qi" 5)"
  if [ $local_impacts = 0 ] && [ -n "$qi" ] && [ $(cnt "$qi" 2) -eq 0 ] && [ $(( $(date +%s) - start )) -gt 7200 ]; then
    hb "    impact jobs idle/held for 2 h: remove them and run the missing ones locally (one at a time)"
    condor_rm $ICL > /dev/null 2>&1; sleep 30; local_impacts=1
    for m in $MASSES; do for k in wilson xsec; do
      j=$C/$( [ $k = xsec ] && echo output_impacts_xsec || echo output_impacts )/${m}_impacts.json
      fresh $j $T0 && continue
      IMPACT_CARD=$k bash $F/shellScripts/impact/Condor/runjob_expectedImpact.sh $m $F ${CMSSW_SRC%/src} > $L/impact_local_${k}_m$m.log 2>&1 || hb "    local impact $k m$m rc!=0"
    done; done
    qi=""
  fi
  nh=$(cnt "$qb" 5)
  if [ $nh -gt 0 ] && [ $(cnt "$qb" 1) -eq 0 ] && [ $(cnt "$qb" 2) -eq 0 ]; then
    condor_q $BCL -hold -af ProcId HoldReason > $L/bias_held.txt 2>&1
    die "bias: $nh jobs held after the automatic releases (bias_held.txt)"
  fi
  [ -z "$qb" ] && [ -z "$qi" ] && break
  sleep 1800
done
hb "    merge bias"
( cd $F && bash shellScripts/bias/Condor/merge_bias_highstat.sh $MASSES ) > $L/bias_merge.log 2>&1 || die "merge_bias_highstat (bias_merge.log)"
hb "    redraw r_vs_bias of the other 25 points (axis titles)"
for m in $(seq 1 30); do
  case " $MASSES " in *" $m "*) continue;; esac
  ( cd $HS/mA_$m/merged && BIAS_BASE_OVERRIDE=$PWD python3 $BN/plot_bias.py --mA $m ) >> $L/bias_replot.log 2>&1 || hb "    plot_bias m$m rc!=0"
done
( BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh ) > $L/bias_collect.log 2>&1 || hb "    collect_bias_results.sh rc!=0 (bias_collect.log)"
python3 - $HS "$MASSES" $C $T0 <<'PY' | tee -a $L/chain_heartbeat.txt || die "bias or impacts missing"
import json, os, sys
HS, masses, C, T0 = sys.argv[1], sys.argv[2].split(), sys.argv[3], float(sys.argv[4]); bad = 0
for m in masses:
    p = "%s/mA_%s/merged/BiasJson/%s_gaussfit.json" % (HS, m, m)
    if not os.path.exists(p): print("    bias mA%s MISSING" % m); bad += 1; continue
    fr = json.load(open(p))["fit_results"]
    print("    bias mA%s " % m + " ".join("%s:%+.3f" % (k, v["mean"]) for k, v in sorted(fr.items(), key=lambda x: -abs(x[1]["mean"]))))
    for k, d in (("wilson", "output_impacts"), ("xsec", "output_impacts_xsec")):
        q = "%s/%s/%s_impacts.json" % (C, d, m)
        if not (os.path.exists(q) and os.path.getmtime(q) >= T0): print("    impacts %s mA%s MISSING/stale" % (k, m)); bad += 1; continue
        j = json.load(open(q)); r = [x for x in j["POIs"] if x["name"] == "r"][0]["fit"]
        print("    impacts %-6s mA%s r = %.3f -%.3f/+%.3f" % (k, m, r[1], r[1] - r[0], r[2] - r[1]))
sys.exit(bad)
PY
fi

if run 7; then
hb "[7] limit plots"
T=$(date +%s)
bash $F/Plots/2_runLimitsPlot.sh > $L/limitsplot.log 2>&1 || hb "    2_runLimitsPlot.sh rc!=0 (limitsplot.log)"
n=$(find $F/Plots/plot_limits -newermt "@$T" -type f | wc -l); [ $n -ge 1 ] || die "no fresh limit plots"
hb "    $n fresh limit-plot files"
fi

if run 8; then
hb "[8] regression check"
manifest > $L/manifest.after
python3 - $L/manifest.before $L/manifest.after <<'PY' | tee $L/regression.txt
import re, sys
def rd(p): return {l[34:].strip(): l[:32] for l in open(p)}
a, b = rd(sys.argv[1]), rd(sys.argv[2])
changed = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
m = r"(1[89]|2[012])"
pat = re.compile(r"Signal/outdir_ele/signalFit/(output|Plots)/(.*20_.*(2023postBPix|run3)|effSigmaVmA)|(^|/)%s_|higgsCombine%s\.|bias_outputs_highstat/mA_%s/|plots_bias/[0-9]+_r_vs_bias|plots_bias/%s_|bias_summary|Plots/plot_limits/" % (m, m, m, m))
other = [k for k in changed if not pat.search(k)]
print("    changed files: %d, outside the expected set: %d" % (len(changed), len(other)))
for k in other[:80]: print("      ", k)
sys.exit(1 if other else 0)
PY
[ ${PIPESTATUS[0]} -eq 0 ] || hb "    WARNING: files outside the expected set changed (regression.txt), check by hand"
fi

hb "CHAIN DONE"
date '+%F %T' > $L/chain.DONE
