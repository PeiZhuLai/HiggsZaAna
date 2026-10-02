#!/usr/bin/env bash
# Expected impacts (30 mA) + high-stat bias (30 mA x 10 chunks) for the 2024 DY+jets overlap-veto
# retrain (datacards with S1 scheme-2 a_w lnN, CSEV a_w, 2026-09-30 05:53).
# Same flow as chain_impbias.sh / chain_csev_init2829.sh:
#   * bias jobs read a FROZEN copy of the workspaces (each impact job re-runs text2workspace into
#     Combine/root_t2w/<mA> while bias would be reading it)
#   * previous products are kept: bias_outputs_highstat -> *_preDYveto_20260930, plots_bias ->
#     *_preDYveto_20260930, output_impacts copied to output_impacts_preDYveto_20260930
#   * done-markers are products newer than t0, never exit codes or an empty queue
#   * held for memory: RequestMemory 8000 -> 16000 -> 20000, then leave held and report
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FZ=$C/root_t2w_dyveto_20260930
TAG=preDYveto_20260930
HB=$L/dyveto_impbias_heartbeat.txt
hb(){ echo "$(date '+%F %T') $*" | tee -a "$HB"; }

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv() { eval "$(cd $F/.. && scramv1 runtime -sh)"; }
export -f cmsenv
cmsenv
set -u

T0=$(date +%s); echo $T0 > $L/dyveto_impbias.t0
hb "start: freeze workspaces, park previous products"
# workspaces must come from the a_w datacards (t2w of 2026-09-30 05:3x)
for m in $(seq 1 30); do
  w=$C/root_t2w/${m}_Datacard_leptons.root
  [ -s $w ] && [ "$(stat -c %Y $w)" -ge "$(date -d '2026-09-30 05:31' +%s)" ] || { hb "STOPPING: $w missing or older than the a_w rerun"; exit 1; }
done
mkdir -p $FZ && cp -p $C/root_t2w/*_Datacard_leptons.root $FZ/
[ -d $BN/bias_outputs_highstat ] && [ ! -d $BN/bias_outputs_highstat_$TAG ] && mv $BN/bias_outputs_highstat $BN/bias_outputs_highstat_$TAG
[ -d $BN/plots_bias ] && [ ! -d $BN/plots_bias_$TAG ] && mv $BN/plots_bias $BN/plots_bias_$TAG
[ -d $C/output_impacts ] && [ ! -d $C/output_impacts_$TAG ] && cp -rp $C/output_impacts $C/output_impacts_$TAG

hb "submit bias (30 x 10 chunks) + impacts (30)"
bout=$(ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $(seq 1 30) 2>&1)
iout=$(bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $(seq 1 30) 2>&1)
echo "$bout" > $L/dyveto_bias_submit.log; echo "$iout" > $L/dyveto_impact_submit.log
CLUSTERS=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$bout
$iout" | sort -u | tr '\n' ' ')
[ -n "${CLUSTERS// /}" ] || { hb "STOPPING: nothing submitted (dyveto_*_submit.log)"; exit 1; }
hb "clusters: $CLUSTERS"

count_products() {
  imp=0; for m in $(seq 1 30); do f=$C/output_impacts/${m}_impacts.json; [ -s $f ] && [ $(stat -c %Y $f) -ge $T0 ] && imp=$((imp+1)); done
  chunks=0; for m in $(seq 1 30); do for c in $(seq 0 9); do
    compgen -G "$BN/bias_outputs_highstat/mA_${m}/chunk${c}/BiasFits/*_split*_fits.root" >/dev/null && chunks=$((chunks+1)); done; done
}
while true; do
  q=$(condor_q $CLUSTERS -af JobStatus 2>/dev/null)
  if [ $? -ne 0 ]; then sleep 300; continue; fi
  idle=$(grep -c '^1$' <<<"$q"); run=$(grep -c '^2$' <<<"$q"); held=$(grep -c '^5$' <<<"$q")
  # memory ladder for held jobs
  while read -r id mem reason; do
    [ -z "${id:-}" ] && continue
    case "$reason" in *emory*)
      if [ "$mem" -lt 8000 ]; then new=8000; elif [ "$mem" -lt 16000 ]; then new=16000; elif [ "$mem" -lt 20000 ]; then new=20000; else continue; fi
      hb "held $id mem $mem -> $new: condor_qedit $id RequestMemory $new; condor_release $id"
      condor_qedit $id RequestMemory $new >/dev/null; condor_release $id >/dev/null ;;
    esac
  done < <(condor_q $CLUSTERS -constraint 'JobStatus==5' -af:j RequestMemory HoldReason 2>/dev/null | awk '{id=$1; mem=$2; $1="";$2=""; print id, mem, $0}')
  count_products
  hb "idle=$idle running=$run held=$held | impacts $imp/30 | bias chunks $chunks/300"
  [ $((idle + run + held)) -eq 0 ] && break
  [ $((idle + run)) -eq 0 ] && [ "$held" -gt 0 ] && { hb "only held jobs left"; break; }
  sleep 1800
done

count_products
missing_imp=""; for m in $(seq 1 30); do f=$C/output_impacts/${m}_impacts.json; { [ -s $f ] && [ $(stat -c %Y $f) -ge $T0 ]; } || missing_imp="$missing_imp $m"; done
missing_bias=""; for m in $(seq 1 30); do n=0; for c in $(seq 0 9); do
  compgen -G "$BN/bias_outputs_highstat/mA_${m}/chunk${c}/BiasFits/*_split*_fits.root" >/dev/null && n=$((n+1)); done
  [ $n -eq 10 ] || missing_bias="$missing_bias mA$m:$n/10"; done
hb "impacts missing:${missing_imp:- none} | bias incomplete:${missing_bias:- none}"

ROOT_DATACARD_PATH=$FZ bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh $(seq 1 30) > $L/dyveto_bias_merge.log 2>&1
nj=0; bad=""; for m in $(seq 1 30); do
  j=$BN/bias_outputs_highstat/mA_${m}/merged/BiasJson/${m}_gaussfit.json
  p=$BN/bias_outputs_highstat/mA_${m}/merged/plots_bias/${m}_r_vs_bias.pdf
  if [ -s $j ] && [ -s $p ]; then nj=$((nj+1)); else bad="$bad $m"; fi; done
hb "merged: $nj/30${bad:+ missing:$bad}"
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh > $L/dyveto_bias_collect.log 2>&1

# worst |mean pull| per mA (AN threshold 0.2)
python3 - $BN/bias_outputs_highstat > $L/dyveto_bias_summary.txt <<'PY'
import glob, json, os, sys
B = sys.argv[1]
print("mA  worst_function  mean_pull   (all truth functions)")
over = []
for m in range(1, 31):
    p = "%s/mA_%d/merged/BiasJson/%d_gaussfit.json" % (B, m, m)
    if not os.path.exists(p): print("%2d  MISSING" % m); continue
    d = json.load(open(p))
    items = []
    for k, v in d.get("fit_results", d).items():
        mu = v.get("mean", v.get("mu")) if isinstance(v, dict) else None
        if mu is not None: items.append((k, float(mu)))
    if not items: print("%2d  (unparsed) %s" % (m, str(d)[:120])); continue
    k, mu = max(items, key=lambda x: abs(x[1]))
    print("%2d  %-14s %+.3f   %s" % (m, k, mu, " ".join("%s:%+.3f" % x for x in items)))
    over += [(m, a, b) for a, b in items if abs(b) > 0.2]
print("above 0.2:", over if over else "none")
PY
cat $L/dyveto_bias_summary.txt >> $HB
if [ -z "$missing_imp" ] && [ -z "$missing_bias" ] && [ $nj -eq 30 ]; then hb "VERDICT: PASSED"; else hb "VERDICT: FAILED (see above)"; fi
echo "$(date '+%F %H:%M') DONE dyveto impacts/bias (heartbeat dyveto_impbias_heartbeat.txt, summary dyveto_bias_summary.txt)" >> $L/PENDING_RESULTS.txt
