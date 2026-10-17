#!/usr/bin/env bash
# Detached: waits for the chainlow_relE2w123 bias cluster, merges the chunks per mass (merge_bias_highstat.sh of the sandbox copy,
# xsec card), writes chain4/bias_summary.txt (pull mean per envelope member, chain4 vs production) and a DONE/STOP line
# in ../logs/PENDING_RESULTS.txt. Held jobs are left to the cluster's periodic_release; if every remaining job is
# held, it stops and reports. Numbers only.
CL=$1
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chainlow_relE2w123
S=$CH/flashggFinalFit; REG=$(dirname $CH)/logs/PENDING_RESULTS.txt; ST=$CH/logs/watch_bias.status
E=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2w123/bias_outputs_highstat
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal/bias_outputs_highstat
while [ "$(condor_q $CL -af ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
  n=$(condor_q $CL -af JobStatus 2>/dev/null); h=$(echo "$n" | grep -c '^5$'); t=$(echo "$n" | grep -c .)
  echo "$(date '+%F %T') bias cluster $CL in queue $t held $h" > $ST
  [ "$t" -gt 0 ] && [ "$h" -eq "$t" ] && { echo "$(date '+%F %H:%M') STOP  chainlow_relE2w123 bias: all $t remaining jobs held" >> $REG; exit 1; }
  sleep 600
done
echo "$(date '+%F %T') merging" > $ST
HS_OUT_BASE=$E ROOT_DATACARD_PATH=$S/Combine/root_t2w_xsec BASE_DIR=$S bash $S/shellScripts/bias/Condor/merge_bias_highstat.sh 1 2 3 > $CH/logs/bias_merge.log 2>&1
python3 - $E $F > $CH/bias_summary.txt <<'PY'
import json, os, sys
new, prod = sys.argv[1], sys.argv[2]
def get(base, m):
    p = "%s/mA_%d/merged/BiasJson/%d_gaussfit.json" % (base, m, m)
    return json.load(open(p))["fit_results"] if os.path.exists(p) else None
print("bias pull mean per envelope member (xsec card for chain4); |mean| > 0.2 flagged")
for m in (1, 2, 3):
    n, p = get(new, m), get(prod, m)
    fmt = lambda d: " ".join("%s:%+.3f%s" % (k, v["mean"], "!" if abs(v["mean"]) > 0.2 else "") for k, v in sorted(d.items())) if d else "MISSING"
    print("mA%-2d chain4 %-60s | prod %s" % (m, fmt(n), fmt(p)))
PY
miss=$(grep -c "chain4 MISSING" $CH/bias_summary.txt)
echo "$(date '+%F %H:%M') DONE  chainlow_relE2w123 bias merged ($miss/3 missing)  -> $CH/bias_summary.txt" >> $REG
echo "$(date '+%F %T') done" > $ST
