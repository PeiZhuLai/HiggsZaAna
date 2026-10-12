#!/bin/bash
# Local validation of the opt-in trigdenom tagger on job 1 of mA_M1_2024 (condor queue is
# slow: user fair-share debt). Same config the condor job would run; stdout = the .out log.
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v2_trigdenom_20260927/localtest
export X509_USER_PROXY=$PWD/x509up_u175325 PYTHONNOUSERSITE=1 HZA_STAGE_INPUTS=1
export TMPDIR=$PWD/scratch; mkdir -p $TMPDIR
export PYTHONPATH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
export PATH=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin:/cvmfs/cms.cern.ch/common:/usr/bin:/bin
echo "START $(date)"
/usr/bin/time -v nice -n 10 python -s /eos/home-p/pelai/HZa/trigdenom_V2_20260927/Sig_MC_trigdenomV2TEST/mA_M1_2024/job_1/mA_M1_2024_executable_job1.py
echo "END rc=$? $(date)"
