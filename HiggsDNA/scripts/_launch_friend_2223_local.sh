#!/bin/bash
# friend 生產:submit 用本地 conda-pack env(init 快,local disk)+include fix(ROOT cling 找 assert.h);
# worker 走 B-mode。FPO=50。趕在背景進程被清前送出 tag jobs(jobs 上 condor 即持久)。
LOCALENV=/tmp/pelai/claude-175325/-afs-cern-ch-user-p-pelai/25618aaa-0524-4f8a-86a8-54b12772b793/scratchpad/hza_ana_local
SYS=$LOCALENV/x86_64-conda-linux-gnu/sysroot/usr/include
export C_INCLUDE_PATH=$SYS:/usr/include
export CPLUS_INCLUDE_PATH=$SYS:/usr/include
export CPATH=$SYS:/usr/include
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
ENVDIR=$LOCALENV FPO=50 CATALOG=metadata/samples/za_merged_data_2223_perds.json \
bash scripts/orchestrate_data_ml_production.sh
