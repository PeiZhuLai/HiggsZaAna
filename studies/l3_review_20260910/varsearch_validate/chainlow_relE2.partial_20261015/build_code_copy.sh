#!/usr/bin/env bash
# Build the code-only sandbox copy of flashggFinalFit for the step-4 full chain of the variable-search BDTs
# (30 m_a points). Production is only read. Copied: source/script/config files of the directories the chain uses
# (no .root/.pdf/.png/logs/outputs). Every production path is then pointed at the sandbox:
#   flashggFinalFit code dir          -> chain4/flashggFinalFit
#   /eos/.../HZa/root_MVAcut          -> /eos/.../sandbox_varsearch_v1/root_MVAcut
#   scored inputs (fsrfix)            -> /eos/.../sandbox_varsearch_v1/run3_bdt_scored
#   Plot/output/MVAcut_points_run3.json, sigEfficiencyVmA_*.json -> chain4/plot_output/
#   S1 lnN default JSON               -> chain4/s1/s1_lnn_scheme2_chain4.json
# 3_rysn_datacard.sh (copies to ~/public/hza_datacard) is replaced by a no-op.
# Output: chain4/flashggFinalFit, chain4/logs/build_residual_paths.txt (must be empty of write targets)
set -euo pipefail
C4=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chain4/flashggFinalFit
V4=$C4/shellScripts/bias/Condor
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chainlow_relE2
S=$CH/flashggFinalFit
SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2   # MVA-cut outputs of this chain
[ -e $S ] && { echo "$S exists; remove it first"; exit 1; }
mkdir -p $S $CH/logs $CH/plot_output $CH/s1
cd $F
DIRS="tools Signal Trees2WS MVAcut/run3_ReReco_Sys/scripts Datacard Combine shellScripts utils tdrStyle"
find $DIRS setup.sh \
  \( -type d \( -name 'output*' -o -name 'logs*' -o -name '*_logs*' -o -name 'logfiles' -o -name 'bias_jobs*' -o -name 'bias_outputs*' \
       -o -name 'condor_logs*' -o -name 'outdir_*' -o -name 'root_t2w*' -o -name 'snapshot*' -o -name 'plots*' -o -name 'yields_*' \
       -o -name '__pycache__' -o -name '_mA*_run' -o -name 'run2_Datacard' -o -name 'Impacts*' -o -name 'impacts*' \) -prune \) -o \
  \( -type f \( -name '*.py' -o -name '*.sh' -o -name '*.C' -o -name '*.cc' -o -name '*.h' -o -name '*.json' -o -name '*.dat' \
       -o -name '*.sub' -o -name '*.submit' -o -name '*.so' -o -name '*.pcm' -o -name 'makefile' -o -name 'Makefile' \) \
       ! -name '*.bak*' -print \) > $CH/logs/build_filelist.txt
echo "files: $(wc -l < $CH/logs/build_filelist.txt)"
rsync -a --files-from=$CH/logs/build_filelist.txt $F/ $S/
cd $S
grep -rlI --exclude=*.so --exclude=*.pcm -e "flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit" -e "/eos/home-p/pelai/HZa/root_" -e "Plot/output/" \
  -e "s1_lnn_scheme2_" . > $CH/logs/build_patched_files.txt || true
while read -r f; do
  sed -i -e "s#/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit#$S#g" \
         -e "s#/eos/home-p/pelai/HZa/root_MVAcut#$SBE/root_MVAcut#g" \
         -e "s#/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix#/eos/home-p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored_low_swapPt_relE2#g" \
         -e "s#/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json#$CH/plot_output/MVAcut_points_run3_chain4.json#g" \
         -e "s#/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sigEfficiencyVmA_#$CH/plot_output/sigEfficiencyVmA_#g" "$f"
done < $CH/logs/build_patched_files.txt
# S1 lnN default (two-line string in apply_bdt_sig.py)
python3 - $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py $CH/s1/s1_lnn_scheme2_chain4.json <<'PY'
import re, sys
p, new = sys.argv[1], sys.argv[2]; s = open(p).read()
t, n = re.subn(r'DEFAULT_MVA_REWEIGHT_LNN_JSON = \([^)]*\)', 'DEFAULT_MVA_REWEIGHT_LNN_JSON = "%s"  # sandbox chain4' % new, s)
assert n == 1, "DEFAULT_MVA_REWEIGHT_LNN_JSON not found"; open(p, "w").write(t)
PY
# commonObjects cwd__
sed -i "s|^cwd__ = .*|cwd__ = \"$S\"  # sandbox chain4 copy, not \$CMSSW_BASE|" $S/tools/commonObjects.py
grep -q "^cwd__ = \"$S\"" $S/tools/commonObjects.py
# setup.sh must put the sandbox tools on PYTHONPATH (scripts source it before adding their own path, so the
# production commonObjects would otherwise be imported first and cwd__ would point at production)
sed -i 's|${CMSSW_BASE}/src/flashggFinalFit/tools|'$S'/tools|' $S/setup.sh
grep -q "$S/tools" $S/setup.sh
# no copy to ~/public
printf '#!/usr/bin/env bash\necho "3_rysn_datacard.sh disabled in the sandbox (would copy to ~/public/hza_datacard)"\n' > $S/Datacard/3_rysn_datacard.sh
# residual production references (any hit here must be read-only, e.g. the background binary)
grep -rnI --exclude=*.so --exclude=*.pcm -e "flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit" -e "/eos/home-p/pelai/HZa/root_" -e "Plot/output" \
  -e "public/hza_datacard" -e "CMSSW_BASE" -e "/eos/project" -e "phys_susy" . > $CH/logs/build_residual_paths.txt || true
echo "patched files: $(wc -l < $CH/logs/build_patched_files.txt); residual lines: $(wc -l < $CH/logs/build_residual_paths.txt)"
# low-mass only: mass lists of the drivers reduced to 1 2 3 (the sig_sys scripts keep theirs: the chain runs them per mass)
sed -i 's/^mAs = \[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30\]/mAs = [1,2,3]/' $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py
grep -q '^mAs = \[1,2,3\]$' $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py
python3 - $S/shellScripts/mva/run_apply_bdt_sig_6jobs.sh <<'PYS'
import re, sys
p = sys.argv[1]; s = open(p).read()
t, n = re.subn(r'samples=\(\n[^)]*\)', 'samples=(\n  "mA_M1" "mA_M2" "mA_M3"\n)', s)
assert n == 1; open(p, "w").write(t)
PYS
sed -i -e 's/^mAs_sig=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)$/mAs_sig=(1 2 3)/' -e 's/^mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)$/mAs_data=(1 2 3)/' $S/Trees2WS/run_tree2ws.sh
grep -q '^mAs_sig=(1 2 3)$' $S/Trees2WS/run_tree2ws.sh && grep -q '^mAs_data=(1 2 3)$' $S/Trees2WS/run_tree2ws.sh
for f in 1_runjob_gen_datacard_makeYields.sh 2_runjob_gen_datacard_makeDatacard.sh; do
  sed -i 's/^mAs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )$/mAs=( 1 2 3 )/' $S/Datacard/$f
  grep -q '^mAs=( 1 2 3 )$' $S/Datacard/$f
done
# bias scripts: the chain4 copies (job work dir in the condor scratch, outputs to HS_OUT_BASE on EOS; 10-15 AFS lesson)
for f in runjob_bias_highstat.sh merge_bias_highstat.sh subjob_bias_highstat.sh; do
  cp $V4/$f $S/shellScripts/bias/Condor/$f
  sed -i "s#$C4#$S#g" $S/shellScripts/bias/Condor/$f
  grep -q 'HS_OUT_BASE' $S/shellScripts/bias/Condor/$f
done
du -sh $S
