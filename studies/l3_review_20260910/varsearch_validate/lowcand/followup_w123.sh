#!/usr/bin/env bash
# Detached follow-up of the c3 validation (user decision 2026-10-15: low-mass = low_swapPt_relE2 features with window
# training at m_a 1-3; m_a = 1 envelope without Laurent, as m_a 24/25 in production). Steps:
#   1 wait for lowcand/logs/drive_c3.status DONE (STOP -> stop)
#   2 gate at the x0.5 cuts (numbers -> lowcand/gate_w123.txt): closure Z < 1 and data GOF min >= 0.05 at m_a 1-3,
#     bootstrap |mean rhat| / sd < 1 at m_a 1 (x0.5), 2, 3; failure -> STOP for the user to decide
#   3 m_a 1-3 signal re-scoring with low_relE2_w123 (51 condor jobs, rescore.sub)
#   4 make_chainlow.py relE2w123 + sandbox fTest without Laurent at m_a = 1 (fTest_noLau_mA1/bin) + build_code_copy.sh
#   5 chainlow.sh under the drive_chainlow watchdog (copy), then watch_bias (copy) -> chainlow_relE2w123/bias_summary.txt
# Registry lines in logs/PENDING_RESULTS.txt.
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
L=$V/lowcand; REG=$V/logs/PENDING_RESULTS.txt; ST=$L/logs/followup_w123.status
N=low_relE2_w123; TAG=relE2w123; CH=$V/chainlow_$TAG; SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
stat() { echo "$(date '+%F %T') $*" > $ST; }
stop() { stat "STOP $*"; echo "$(date '+%F %H:%M') STOP  w123 follow-up  $*" >> $REG; exit 1; }
# 1
while :; do s=$(cat $L/logs/drive_c3.status 2>/dev/null)
  case "$s" in *DONE*) break;; *STOP*) stop "c3 driver stopped: $s";; esac
  stat "waiting for c3: $s"; sleep 300; done
# 2
env -i HOME=$HOME PATH=/usr/bin:/bin $PY - $V $L $N > $L/gate_w123.txt 2>&1 <<'PYG'
import re, sys
import numpy as np
V, L, N = sys.argv[1:4]
cuts = {int(l.split()[0]): l.split()[1:] for l in open("%s/eval/closure_cuts_%s.txt" % (V, N))}
ok = True; wp = []
for m in (1, 2, 3):
    c = cuts[m][4]; wp.append(c)   # x0.5
    line = [l for l in open("%s/scan_%s_mA%d.txt" % (L, N, m)) if abs(float(l.split()[0]) - float(c)) < 1e-6][0]
    z = float(re.search(r"closureZ=(\S+)", line).group(1))
    gof = [float(x) for x in re.findall(r":([0-9.]+)", line.split("| data GOF:")[1].split("| ps GOF")[0])]
    rows = [x.split() for x in open("%s/boot_%s_mA%d_%s.txt" % (L, N, m, c.replace(".", ""))) if "rhat=" in x]
    r = np.array([float(x[3][5:]) for x in rows if x[3][5:] not in ("NA", "")])
    sig = abs(r.mean()) / r.std(ddof=1)
    p = z < 1 and min(gof) >= 0.05 and sig < 1 and r.size >= 25
    ok &= p
    print("mA%d cut %s: closure Z %.2f, data GOF min %.2f, bootstrap %d rep rhat %+.3f +- %.3f (%.2f sd) -> %s"
          % (m, c, z, min(gof), r.size, r.mean(), r.std(ddof=1), sig, "pass" if p else "FAIL"))
print("WP", " ".join(wp)); print("GATE", "PASS" if ok else "FAIL")
PYG
grep -q "^GATE PASS" $L/gate_w123.txt || stop "gate at x0.5 not passed: $(grep -v '^WP' $L/gate_w123.txt | tr '\n' ';')"
read -r _ c1 c2 c3 < <(grep "^WP" $L/gate_w123.txt)
# 3
grep -E "^mA_M[123](_|$)" $V/condor/joblist_sig.txt | awk -v o=$SBE/run3_bdt_scored_$N -v n=$N '{print $1, o, n, "high"}' > $L/list_rescore_sig_w123.txt
[ "$(wc -l < $L/list_rescore_sig_w123.txt)" -eq 51 ] || stop "signal rescore list is not 51 lines"
CL=$(cd $L && condor_submit rescore.sub LIST=list_rescore_sig_w123.txt 2>&1 | grep -oP 'submitted to cluster \K[0-9]+')
[ -n "$CL" ] || stop "signal rescore submission failed"
stat "signal rescore cluster $CL"
# 4
python3 $V/make_chainlow.py $TAG $N $c1 $c2 $c3 > $L/logs/make_chainlow_w123.log 2>&1 || stop "make_chainlow.py failed"
python3 - $CH/chainlow.sh <<'PYF' || stop "fTest patch failed"
import sys
p = sys.argv[1]; s = open(p).read()
a = "( cd $F/Background && ./bin/fTest_ALP_turnOn "
b = "( cd $F/Background && /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/fTest_noLau_mA1/bin/fTest_ALP_turnOn "
assert s.count(a) == 1, s.count(a)
s = s.replace(a, b).replace("#!/usr/bin/env bash\n", "#!/usr/bin/env bash\n# 10-15: background envelopes with the sandbox fTest (no Laurent at m_a = 1; varsearch_validate/fTest_noLau_mA1)\n", 1)
open(p, "w").write(s)
PYF
bash $CH/build_code_copy.sh > $CH/logs/build_code_copy.log 2>&1 || stop "build_code_copy.sh failed ($CH/logs/build_code_copy.log)"
sed -e "s#chainlow_relE2#chainlow_$TAG#g" -e "s#low_swapPt_relE2 x0.5#$N x0.5#" $V/chainlow_relE2/drive_chainlow.sh > $CH/drive_chainlow.sh
sed -e "s#chainlow_relE2#chainlow_$TAG#g" $V/chainlow_relE2/watch_bias.sh > $CH/watch_bias.sh
grep -q "chainlow_relE2\b" $CH/drive_chainlow.sh $CH/watch_bias.sh && stop "path substitution incomplete"
# 5
stat "chainlow running (rescore cluster $CL first)"
bash $CH/drive_chainlow.sh $CL
[ -e $CH/logs/chain4.DONE ] || stop "chainlow did not finish ($(cat $CH/logs/chain4.FAILED 2>/dev/null))"
BC=$(grep -oP 'submitted to cluster \K[0-9]+' $CH/logs/bias_submit.log | head -1)
[ -n "$BC" ] || stop "no bias cluster in bias_submit.log"
stat "bias cluster $BC"
bash $CH/watch_bias.sh $BC
stat "DONE"
