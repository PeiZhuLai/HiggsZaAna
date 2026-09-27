#!/usr/bin/env bash
#
# p2root for the FSR-fix production.
#
# Inputs are split across two directories and 1_run_P2Root.sh was updated for it:
#   Bkg_MC, Data -> parquet_DNA_tmp_fsrfix_fpo1   (fpo=1 re-run)
#   Sig_MC       -> parquet_DNA_tmp_fsrfix        (complete, reconciled, not re-run)
# Target is a NEW directory so the AN's current run3_bdt_inputs_nominal is not overwritten.
#
# Known and accepted input gaps (2026-09-19): Data_2022preEE is short 12 DoubleMuon
# Run2022C files vs the previous production; Data_2024 has 93 extra files. Both were
# diagnosed per channel and accepted; see PENDING_RESULTS.txt.
#
# 2026-09-20: the first attempt used the WRONG conda env (hza_ana) and every one of
# the 1209 conversions died, yet 1_run_P2Root.sh printed "completed successfully" and
# exited 0 -- its retry loop swallows failures. RC is worthless here, which is why this
# wrapper counts the produced ROOT files itself.
#
# Two separate faults were behind it, both worth remembering:
#   1. hza_ana has no xgboost at all (Parque2Root_BDT.py line 16 imports XGBClassifier).
#   2. In hza_ana, `import pandas` pulls in the SYSTEM /lib64/libstdc++.so.6 (el9 tops
#      out at GLIBCXX_3.4.30); the later `import ROOT` then cannot dlopen
#      libcppyy_backend.so, which needs GLIBCXX_3.4.31. Importing ROOT alone works fine,
#      so testing `python -c "import ROOT"` gives a FALSE PASS -- the import ORDER in the
#      real script is what breaks it. Symptom is "could not load cppyy_backend library",
#      which names neither pandas nor libstdc++.
# The env that actually has xgboost + ROOT + pandas coexisting is higgs-alp-ana.
#
set -uo pipefail
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna
OUT=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix

# `set -u` must be off across conda activate: activate.d/activate-binutils_linux-64.sh
# dereferences ADDR2LINE et al. unguarded and aborts the whole script under -u.
# An interactive test never has -u set, so it passes there and fails only here.
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana || { echo "FATAL: conda activate higgs-alp-ana failed"; exit 1; }
# Preflight in the SAME import order as the converter, or the check is meaningless.
python -c "
import pandas, uproot
from ROOT import Math, TVector2, TVector3, TLorentzVector
from xgboost import XGBClassifier
" || { echo "FATAL: preflight imports failed, refusing to start"; exit 1; }
echo "[env] python=$(which python)"
set -u

export PYTHONPATH="${PYTHONPATH:-}:$R/HiggsDNA"
mkdir -p "$OUT"
cd "$R/Parquet2Rootfile" && bash 1_run_P2Root.sh
rc=$?

echo "=============== PRODUCT CHECK ==============="
echo "raw exit code from 1_run_P2Root.sh: $rc  (NOT trusted)"
n=$(find "$OUT" -name '*.root' | wc -l)
echo "root files produced: $n"
find "$OUT" -name '*.root' -printf '%TY-%Tm-%Td %TH:%TM  %10s  %p\n' | sort
if [ "$n" -eq 0 ]; then echo "VERDICT: FAILED (no output)"; exit 1; fi
echo "VERDICT: see three-stage reconciliation, file count alone is not sufficient"
