#!/bin/bash
# Local single-job test (ZGG job_1 config copy, outputs redirected here) while the
# condor test jobs wait in the queue. Validates tagger/weights/output on the new sample.
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
export CONDA_PREFIX=$ENVDIR PATH=$ENVDIR/bin:$PATH PYTHONPATH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
export X509_USER_PROXY=/tmp/x509up_u175325 HZA_STAGE_INPUTS=1 PYTHONNOUSERSITE=1
export TMPDIR=/tmp/pelai/claude-175325/-afs-cern-ch-user-p-pelai/4d5b6567-8134-4d1e-8743-829cc5526321/scratchpad
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
python -s -c "
import json
from higgs_dna.utils.logger_utils import setup_logger
setup_logger('INFO', None)
from higgs_dna.analysis import run_analysis
run_analysis(json.load(open('/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v3v8_extrabkg_20260927/localtest/ZGG_job1/config.json')))
"
