#!/usr/bin/env bash
#
# Submit 5 one-off condor jobs that measure copy-then-read vs streaming on real
# worker nodes, one per problem sample.
#
# WHY
#   At the 8 h walltime the production's dominant failure is no longer the wall:
#     exit0 785 | exit1 633 | WALL 386   (recent jobs, MaxRuntime=28800)
#   and every sampled exit-1 carries the same message,
#     OSError: File did not vector_read properly: [ERROR] Operation expired
#   The shortfall tracks job size exactly -- DYGto2LG 770/776 complete, while
#   Data_2024 has 1556/2550 -- which is what a vector-read timeout against
#   high-latency sites looks like. Before changing the shared read path in
#   analysis.py, measure whether xrdcp-then-read is actually faster and whether
#   node-local scratch can hold fpo=4 files.
#
# SCALE: 5 condor jobs. Logs land in this directory. Nothing else is touched.
#
# USAGE
#   bash submit_read_probe.sh            # dry run: print the plan and the .sub
#   DRY_RUN=0 bash submit_read_probe.sh  # submit
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/readtest
S=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana

# One config per sample that is actually short, biggest deficits first.
CFGS=(
  "$S/Data/Data_2024/job_5/Data_2024_config_job5.json"
  "$S/Data/Data_2023preBPix/job_5/Data_2023preBPix_config_job5.json"
  "$S/Bkg_MC/DYJetsTo2Tau_2024/job_5/DYJetsTo2Tau_2024_config_job5.json"
  "$S/Bkg_MC/DYJetsToLL_2022postEE/job_5/DYJetsToLL_2022postEE_config_job5.json"
  "$S/Bkg_MC/DYGto2LG_10to100_2024/job_5/DYGto2LG_10to100_2024_config_job5.json"
)

cat > "$D/probe_wrapper.sh" <<'WEOF'
#!/bin/bash
export X509_USER_PROXY=$PWD/x509up_u175325
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
export PATH="$ENVDIR/bin:$PATH"
export PYTHONNOUSERSITE=1
echo "HOST=$(hostname)  CFG=$1"
"$ENVDIR/bin/python" -s /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/readtest/probe_read_strategy.py "$1" 2
echo "WRAPPER_EXIT=$?"
WEOF
chmod +x "$D/probe_wrapper.sh"

# 8 h so a slow stream cannot be mistaken for a crash; these are 5 jobs.
cat > "$D/probe.sub" <<SEOF
executable   = $D/probe_wrapper.sh
should_transfer_files = Yes
transfer_input_files  = /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/x509up_u175325
output       = $D/probe_\$(ProcId).out
error        = $D/probe_\$(ProcId).err
log          = $D/probe.log
RequestMemory = 4000
RequestDisk   = 30000000
RequestCpus   = 1
JobBatchName  = HZaReadProbe
+JobFlavour   = "workday"
getenv        = True
SEOF
for c in "${CFGS[@]}"; do
    [ -f "$c" ] || { echo "MISSING CONFIG: $c" ; continue ; }
    echo "arguments = $c" >> "$D/probe.sub"
    echo "queue 1"        >> "$D/probe.sub"
done

echo "=== probe.sub ==="; cat "$D/probe.sub"
if [ "$DRY_RUN" = "1" ]; then echo "(dry run: not submitted)"; exit 0; fi
condor_submit "$D/probe.sub"
