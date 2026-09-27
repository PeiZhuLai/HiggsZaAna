#!/usr/bin/env bash
# Revert mA3 to the 0.988 working point (user decision 2026-09-26).
# At 0.940 (R 1.10) the mA3 data grew to 1885 events and no background family described the
# spectrum (GOF: Bern3 1.2e-5, others ~0, envelope entirely fallback), whereas at 0.988 all four
# members have GOF 0.47-0.78 and pass the bias threshold. chain_wp_mA3(_from4).sh had stopped
# before datacards / limits / t2w, so those still hold the 0.988 results and are not touched.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
B=$F/Background
FR=$B/ALP_BkgModel_ReReco/fit_results_run3
BN=$F/Combine/Checks/Bias_nominal
EOSMVA=/eos/home-p/pelai/HZa/root_MVAcut
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
TAG=before_wp0940_20260926_2133
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
die() { echo "STOPPING: $*"; exit 1; }

echo "=============== [1] working point 0.988 restored ==============="
cp -p $PLOT/scripts/scan_score_R_significance.py $PLOT/scripts/scan_score_R_significance.py.wp0940
cp -p $PLOT/scripts/scan_score_R_significance.py.$TAG $PLOT/scripts/scan_score_R_significance.py
grep -n '^FIXED_WP' $PLOT/scripts/scan_score_R_significance.py
cp -p $PLOT/output/MVAcut_points_run3.json $PLOT/output/MVAcut_points_run3.json.wp0940
cd $PLOT && PYTHONPATH=$PLOT/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts $ALP_PY scripts/scan_score_R_significance.py --write-json > $L/wp3_revert_scan.log 2>&1
grep -E '^ +[0-9]+ +0\.' $L/wp3_revert_scan.log
cmp <(python3 -c "import json;print(json.load(open('$PLOT/output/MVAcut_points_run3.json'))['results'])") \
    <(python3 -c "import json;print(json.load(open('$PLOT/output/MVAcut_points_run3.json.$TAG'))['results'])") \
  && echo "  MVAcut_points_run3.json identical to the pre-change backup" || die "JSON differs from backup"

echo "=============== [2] EOS mA_M3 (MVA cut trees + workspaces) restored ==============="
for d in data sig; do
  [ -d $EOSMVA/${d}_variants_mA_M3_$TAG ] || die "backup $EOSMVA/${d}_variants_mA_M3_$TAG missing"
  mv $EOSMVA/$d/mA_M3 $EOSMVA/${d}_variants_mA_M3_wp0940
  cp -rp $EOSMVA/${d}_variants_mA_M3_$TAG $EOSMVA/$d/mA_M3
done
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv
export PYTHONPATH="${PYTHONPATH:-}:$F/tools:$F/Signal/tools"
python3 - $EOSMVA <<'PY' || die "restored trees/workspaces fail a full read"
import glob, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B = sys.argv[1]
f = ROOT.TFile.Open(B + "/data/mA_M3/run3.root"); t = f.Get("DiphotonTree/Data_13p6TeV")
for i in range(t.GetEntries()):
    if t.GetEntry(i) <= 0: print("FAIL data", i); sys.exit(1)
print("  data mA_M3 run3.root: %d entries read" % t.GetEntries())
ws = [B + "/data/mA_M3/ws/run3.root"] + sorted(glob.glob(B + "/sig/mA_M3/ws_Tree2WS/*.root"))
for p in ws:
    g = ROOT.TFile.Open(p)
    if not any(k.ReadObj().InheritsFrom("RooWorkspace") for k in g.GetListOfKeys()): print("EMPTY", p); sys.exit(1)
print("  %d workspaces contain a RooWorkspace" % len(ws))
PY

echo "=============== [3] background fit mA3 restored ==============="
mv $FR/3 $B/fit_results_variants/3_wp0940_20260926
cp -rp $B/fit_results_variants/3_$TAG $FR/3
cat $FR/3/HZAmassInde_fTest/EnvelopeResults.txt | sed 's/^/  /'

echo "=============== [4] signal model mA3 rerun on the restored inputs ==============="
S1=$D/wp3_single
for s in 1_runjob_sig_fTest 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter; do
  echo "  [$(date '+%H:%M:%S')] $s"; bash $S1/$s.sh > $L/wp3_revert_sig_$s.log 2>&1
done
bash $F/shellScripts/sig_sys/5_runjob_sig_plotEffSigma.sh > $L/wp3_revert_sig_5.log 2>&1
python3 - $F/Signal $D/../../../flashgg_run3 $TAG <<'PY' || die "restored signal fit does not match the 0.988 backup"
import glob, os, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
S, tag = sys.argv[1], sys.argv[3]
bk = S + "/sigfit_variants_" + tag
worst = 0.0; n = 0
for new in sorted(glob.glob(S + "/outdir_*/signalFit/output/3_CMS-HGG_sigfit_*_Hm125.root")):
    old = os.path.join(bk, os.path.basename(new))
    if not os.path.exists(old): continue
    def norm(p):
        f = ROOT.TFile.Open(p); w = f.Get("wsig_13p6TeV")
        v = [x for x in w.allFunctions() if x.GetName().endswith("_normThisLumi")]
        val = v[0].getVal() if v else None
        f.Close(); return val
    a, b = norm(new), norm(old)
    if a and b: worst = max(worst, abs(a / b - 1)); n += 1
print("  signal norm vs 0.988 backup: %d files, max rel. diff %.2e" % (n, worst))
sys.exit(0 if n >= 10 and worst < 1e-3 else 1)
PY

echo "=============== [5] bias mA3 restored ==============="
[ -d $BN/bias_outputs_highstat_variants/mA_3_$TAG ] && mv $BN/bias_outputs_highstat_variants/mA_3_$TAG $BN/bias_outputs_highstat/mA_3
ls $BN/bias_outputs_highstat/mA_3/merged/BiasJson/
echo "VERDICT: REVERTED"
