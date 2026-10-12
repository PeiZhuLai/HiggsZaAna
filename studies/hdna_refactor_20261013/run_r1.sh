#!/bin/bash
# Regression run R1: one Run 3 signal job (mA_M1_2024 job_1, file 2394cc2e) with the current
# HZa HiggsDNA code, run locally. Usage: run_r1.sh <tag> [template_config] [input_file]
# Output: /eos/home-p/pelai/HZa/test_hdna_regress/<tag>/
set -eo pipefail
TAG=$1
TEMPLATE=${2:-/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC/mA_M1_2024/job_1/mA_M1_2024_config_job1.json}
INPUT=${3:-/eos/home-p/pelai/HZa/test_hdna_regress/inputs/HZa-Zto2L-ato2G_Par-M-1_TuneCP5_13p6TeV_madgraph-pythia8_RunIII2024Summer24NanoAODv15/2394cc2e-fa67-4340-a2ce-cb9bdc55ab98.root}
OUT=/eos/home-p/pelai/HZa/test_hdna_regress/$TAG
HDNA=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
[ -e "$OUT" ] && { echo "output $OUT exists, refusing to overwrite"; exit 1; }
mkdir -p "$OUT"
$PY - "$TEMPLATE" "$INPUT" "$OUT" <<'PYEOF'
import json, os, sys
tmpl, inp, out = sys.argv[1:4]
c = json.load(open(tmpl))
c["files"] = [inp]
c["dir"] = c["output_dir"] = out
c["skimmed_files"] = [os.path.join(out, os.path.basename(inp).replace(".root", "skimmed.root"))]
c["summary_file"] = os.path.join(out, "summary.json")
json.dump(c, open(os.path.join(out, "config.json"), "w"), indent=1)
PYEOF
cd "$HDNA"
export HZA_STAGE_INPUTS=0 X509_USER_PROXY=/tmp/x509up_u$(id -u)
start=$(date +%s)
$PY -c "
import json
from higgs_dna.utils.logger_utils import setup_logger
setup_logger('INFO', None)
from higgs_dna.analysis import run_analysis
run_analysis(json.load(open('$OUT/config.json')))
" > "$OUT/run.log" 2>&1
echo "rc=$? seconds=$(( $(date +%s) - start ))" | tee -a "$OUT/run.log"
