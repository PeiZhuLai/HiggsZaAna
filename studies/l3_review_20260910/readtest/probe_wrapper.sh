#!/bin/bash
export X509_USER_PROXY=$PWD/x509up_u175325
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
export PATH="$ENVDIR/bin:$PATH"
export PYTHONNOUSERSITE=1
echo "HOST=$(hostname)  CFG=$1"
"$ENVDIR/bin/python" -s /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/readtest/probe_read_strategy.py "$1" 2
echo "WRAPPER_EXIT=$?"
