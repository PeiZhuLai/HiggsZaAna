#!/bin/bash
# Phase 5: definitive envelope fix. fTest_ALP_turnOn.cpp now excludes Bern4 (all mA) and
# Pow3/Exp3 (mA=1). Re-run background(envelope, recompiles)->datacard->limits->bias.
# Signal/MVAcut/Tree2WS reused from phase3; impact unchanged. AN push held.
set -uo pipefail
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna

RUN_HIGGSZA_CLEAN_OUTPUTS=0 \
RUN_HIGGSZA_MERGE_PARQUET=0 \
RUN_HIGGSZA_P2ROOT=0 \
RUN_HIGGSZA_TRAIN_MVA=0 \
RUN_HIGGSZA_P2ROOT_MVA_SCORE=0 \
RUN_HIGGSZA_DETERMINE_MVA_CUT=0 \
RUN_HIGGSZA_PLOT=0 \
RUN_FLASHGG_ENV=1 \
RUN_FLASHGG_MVA_CUT=0 \
RUN_FLASHGG_TREE2WS=0 \
RUN_FLASHGG_BACKGROUND=1 \
RUN_FLASHGG_SIGNAL=0 \
RUN_FLASHGG_DATACARD=1 \
RUN_FLASHGG_COMBINE_LIMITS=1 \
RUN_FLASHGG_PLOT_LIMITS=1 \
RUN_FLASHGG_IMPACT=0 \
RUN_FLASHGG_BIAS=1 \
RUN_FLASHGG_COLLECT_BKG=1 \
RUN_EXIT_CMSSW_ENV=1 \
RUN_UPDATE_AN=0 \
bash 1_grand.sh
