#!/bin/bash
# Writes the condor submit file of the 125.38 bias check (side study, differs from production only in -m 125.38).
# Usage: make_submit.sh <chunks, e.g. "0" or "1 2 3 4 5 6 7 8 9"> <output .submit>
# Truth functions and seeds as in the production high-statistics bias (seed 100000 + 100*mA + chunk), one function
# per job; output split fit files in Bias_nominal/bias_outputs_mh12538/mA_<m>/chunk<c>/BiasFits.
M=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/bias_mh12538
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
CT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4
OUT=$FF/Combine/Checks/Bias_nominal/bias_outputs_mh12538
chunks="$1"; sub="$2"
declare -A FUNCS=( [1]="Bern3 Exp1 Lau1 Pow1" [10]="Bern5 Exp3 Lau1 Pow1" )
cat > "$sub" <<HDR
# Side study (differs from production): production high-statistics bias at mA = 1 and 10 GeV with -m 125.38
# (production limits) instead of the RunBiasStudy.py default -m 125; seeds, toys and fit options as in production.
universe = vanilla
executable = $M/runjob_bias_mh.sh
getenv = False
request_memory = 3000
request_cpus = 1
+JobFlavour = "tomorrow"
batch_name = HZa_bias_mh12538
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
periodic_release = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)
HDR
for m in 1 10; do for f in ${FUNCS[$m]}; do for c in $chunks; do
  tag=mA${m}_${f}_c${c}
  cat >> "$sub" <<JOB

environment = "NTOYS_PER_CHUNK=1000 SELECT_FUNCTION=$f BIAS_MH=125.38 BIAS_HS_ROOT=$OUT"
log = $M/logs/$tag.\$(ClusterId).log
output = $M/logs/$tag.\$(ClusterId).\$(ProcId).out
error = $M/logs/$tag.\$(ClusterId).\$(ProcId).err
arguments = $m $c $FF $CT
queue
JOB
done; done; done
echo "wrote $sub: $(grep -c '^queue' "$sub") jobs"
