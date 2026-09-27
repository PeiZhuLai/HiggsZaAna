#!/usr/bin/env bash
# Wait for the mA8 / mA11 envelope-adjustment bias + impact jobs, then merge and judge.
#   mA11 (Bern<=2): bias 12748267, impact 12748268
#   mA8  (no Pow1): bias 12748411, impact 12748412
# Verdict per mass: every truth function |mean pull| <= 0.2 (AN Sec 7.5 threshold).
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FZ=$C/root_t2w_fsrfix_20260925
CLUSTERS="${CLUSTERS:-12748267 12748268 12748411 12748412}"
MASSES="${MASSES:-8 11}"
HB=$L/${HBNAME:-envbias}_heartbeat.txt

while true; do
  q=$(condor_q $CLUSTERS -af JobStatus 2>/dev/null) || { sleep 300; continue; }
  idle=$(grep -c '^1$' <<<"$q"); run=$(grep -c '^2$' <<<"$q"); held=$(grep -c '^5$' <<<"$q")
  n=0; NT=$(( $(wc -w <<<"$MASSES") * 10 )); for m in $MASSES; do for c in $(seq 0 9); do
    compgen -G "$BN/bias_outputs_highstat/mA_${m}/chunk${c}/BiasFits/*_split*_fits.root" >/dev/null && n=$((n+1)); done; done
  { echo "$(date '+%F %T') idle=$idle running=$run held=$held  bias chunks $n/$NT"
    [ "$held" -gt 0 ] && condor_q $CLUSTERS -hold -af ClusterId ProcId HoldReason | sed 's/^/  HELD /'; } > $HB
  [ $((idle + run + held)) -eq 0 ] && break
  [ $((idle + run)) -eq 0 ] && [ "$held" -gt 0 ] && break
  sleep 900
done
cat $HB

echo "=============== merge ==============="
ROOT_DATACARD_PATH=$FZ bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh $MASSES > $L/${HBNAME:-envbias}_merge.log 2>&1
echo "  merge rc=$?"
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh

echo "=============== verdict ==============="
python3 - "$BN" "$C" $MASSES <<'PY'
import json, os, sys
BN, C, masses = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3:]]
ok = True
for m in masses:
    p = "%s/bias_outputs_highstat/mA_%d/merged/BiasJson/%d_gaussfit.json" % (BN, m, m)
    if not os.path.exists(p):
        print("  mA%d: NO gaussfit json" % m); ok = False; continue
    fr = json.load(open(p))["fit_results"]
    worst = max(abs(v["mean"]) for v in fr.values())
    print("  mA%-2d %s   worst |mean| %.3f %s" % (m, "  ".join("%s:%+.3f±%.3f" % (k, v["mean"], v["mean_err"]) for k, v in fr.items()),
                                              worst, "PASS" if worst <= 0.2 else "ABOVE 0.2"))
    ok = ok and worst <= 0.2
    j = "%s/output_impacts/%d_impacts.json" % (C, m)
    if os.path.exists(j):
        r = json.load(open(j))["POIs"][0]["fit"]; print("        impact r = %.3f [%.3f, %.3f]" % (r[1], r[0], r[2]))
print("VERDICT: %s" % ("PASSED" if ok else "NEEDS DECISION"))
PY
