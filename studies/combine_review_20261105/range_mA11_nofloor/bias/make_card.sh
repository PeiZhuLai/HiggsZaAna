#!/usr/bin/env bash
# Cross-section card on a staging background model: the production card text with the bkg_mass / data_obs shapes
# pointed to the staging multipdf, then text2workspace with the production options (RunText2Workspace.py:
# -m 125.38 higgsMassRange=115,135, model mu_inclusive = no -P).
# Usage : make_card.sh <m> <staging result dir, e.g. fit_results_ms> <tag>
# Input : Datacard/output_Datacard_leptons/<m>_pruned_datacard_leptons_xsec.txt (read only), ../<resdir>/<m>/CMS-HGG_mva_13p6TeV_multipdf.root
# Output: cards/<tag>/<m>_Datacard_leptons.{txt,root}, cards/<tag>/<m>_t2w.log
set -eo pipefail
m=$1; res=$2; tag=$3
H=$(cd "$(dirname "$0")" && pwd); S=$(dirname "$H")
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
src=$FF/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons_xsec.txt
prod=$FF/Background/ALP_BkgModel_ReReco/fit_results_run3/${m}/CMS-HGG_mva_13p6TeV_multipdf.root
new=$S/$res/$m/CMS-HGG_mva_13p6TeV_multipdf.root
[ -s "$new" ] || { echo "missing $new"; exit 1; }
O=$H/cards/$tag; mkdir -p $O
n=$(grep -c "$prod" $src); [ "$n" = 2 ] || { echo "expected 2 bkg/data shapes lines, found $n"; exit 1; }
sed "s#$prod#$new#" $src > $O/${m}_Datacard_leptons.txt
cd $FF/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
cd $O && text2workspace.py ${m}_Datacard_leptons.txt -o ${m}_Datacard_leptons.root -m 125.38 higgsMassRange=115,135 > ${m}_t2w.log 2>&1
ls -l $O/${m}_Datacard_leptons.root
