#!/usr/bin/env bash
# chainlow_relE2w123 (make_chainlow.py, 2026-10-15): chain4.sh reduced to m_a = 1, 2, 3 for the low-mass candidate low_relE2_w123
# (cuts [0.98115, 0.98837, 0.99226]); high-mass model = production, m_a >= 4 untouched. Original header:
# Step 4 of the validation of the variable-search BDTs (user decision 2026-10-14: low = model_low with min dR,
# high = model_high, mA1 at the x0.5 cut): full chain for the 30 m_a points in a sandbox.
#
# Inputs:  scored copy /eos/home-p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored (1190 signal + 18 data/DY files)
#          working points plot_output/MVAcut_points_run3_chain4.json (choose_wp.py, from the closure scans)
# Code:    code-only copy flashggFinalFit/ (build_code_copy.sh), every path pointed at the sandbox; background fits
#          use the production fTest_ALP_turnOn / makeBkgPlots_ALP binaries with explicit sandbox paths.
# Stages (START_STAGE=N skips < N; each stage is judged by its products, not exit codes):
#   1 S1 lnN (eval_alt_reweight + make_s1_lnn_json --reference nominal)
#   2 interpolation efficiency JSON (signal_eff_sumw.py)
#   3 MVA cut: data (30 mA) + signal (14 mA x 5 eras), full-read audit
#   4 Tree2WS
#   5 background envelopes (30 mA, 6 in parallel, local)
#   6 signal model (14 mA, sig_sys 1-4, 6 masses in parallel)
#   7 datacards: makeYields + makeDatacard -> "wilson" card (with QCDscale/pdf/alphaS) and "xsec" card (the same
#     without those three lines; model-independent limit on sigma x B at 0.1 pb, user rule 2026-10-14)
#   8 text2workspace + blind expected limits for both cards, compared with production
#   9 closure with the sandbox xsec card (anchors: step-3 pseudo-data fits; interpolated masses: new scans)
#  10 bias (condor, 30 x 10 chunks x 1000 toys, xsec card)
# Production products are fingerprinted before and after; any change, or any file written into production areas
# since the start, is a failure. Markers: logs/chain4.DONE or logs/chain4.FAILED (reason).
set -uo pipefail
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chainlow_relE2w123
V=$(dirname $CH)
S=$CH/flashggFinalFit
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1
MV=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2w123
SC=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored_low_relE2_w123
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
L=$CH/logs
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
WPJ=$CH/plot_output/MVAcut_points_run3_chain4.json
LNN=$CH/s1/s1_lnn_scheme2_chain4.json
EFFJ_ELE=$CH/plot_output/sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json
EFFJ_MU=$CH/plot_output/sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json
MASS14="1 2 3"
MASS30="1 2 3"
INTERP=""
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
FR=$S/Background/ALP_BkgModel_ReReco/fit_results_run3
DC=$S/Datacard/output_Datacard_leptons
START=${START_STAGE:-1}
mkdir -p $L
run() { [ "$START" -le "$1" ]; }
hb() { echo "$(date '+%F %T') $*" | tee -a $L/chain4_heartbeat.txt; }
die() { hb "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/chain4.FAILED; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }
rm -f $L/chain4.DONE $L/chain4.FAILED

fingerprint() {   # production products that must not change
  { for m in $MASS30; do
      ls -l --time-style=+%s $F/Background/ALP_BkgModel_ReReco/fit_results_run3/$m/CMS-HGG_mva_13p6TeV_multipdf.root
      ls -l --time-style=+%s $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt
      ls -l --time-style=+%s $F/Combine/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root
      ls -l --time-style=+%s /eos/home-p/pelai/HZa/root_MVAcut/data/mA_M$m/run3.root
    done
    for m in $MASS14; do for y in $YEARS; do
      ls -l --time-style=+%s /eos/home-p/pelai/HZa/root_MVAcut/sig/mA_M$m/output_$y.root
      for ch in ele mu; do ls -l --time-style=+%s $F/Signal/outdir_$ch/signalFit/output/${m}_CMS-HGG_sigfit_${y}_${ch}_Hm125.root; done
    done; done
    ls -l --time-style=+%s $PLOT/output/MVAcut_points_run3.json $PLOT/output/sigEfficiencyVmA_*_quadratic_interp_ma_points.json
    ls -ld --time-style=+%s /afs/cern.ch/user/p/pelai/public/hza_datacard/mA_M1
  } 2>&1 | awk '{print $6, $7, $5}'
}
[ -s $L/prod_fingerprint.before ] || fingerprint > $L/prod_fingerprint.before
[ -e $L/chain_start.marker ] || touch $L/chain_start.marker

hb "chain4 start (START_STAGE=$START)"
[ -s $WPJ ] || die "no working-point JSON $WPJ (run choose_wp.py)"
python3 - $WPJ <<'PY' || die "working-point JSON does not hold the 14 anchors"
import json, sys
r = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(sys.argv[1]))["results"]}
assert sorted(r) == [1, 2, 3], sorted(r)
print("  WP:", " ".join("%d:%.5f" % kv for kv in sorted(r.items())))
PY
# scored inputs complete: every signal dir x 5 eras plus data/DY, none older than the models
lst=$(xrdfs eosuser.cern.ch ls -R $(echo $SC | sed 's#^/eos/home-p/pelai/#/eos/user/p/pelai/#') 2>/dev/null)
n=$(echo "$lst" | grep -c '\.root$'); p=$(echo "$lst" | grep -c '\.partial$')
[ "$n" -eq 273 ] && [ "$p" -eq 0 ] || die "scored copy incomplete: $n/273 files (255 signal m_a 1-3 + 18 data/DY), $p partial"

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }
export -f cmsenv
cmsenv
export PYTHONPATH="$S/tools:$S/Signal/tools"
set -u
python3 -c "import commonObjects as c; assert c.cwd__ == '$S', c.cwd__; print('  [env] cwd__ =', c.cwd__)" || die "commonObjects does not point at the sandbox"
unset HZA_MVA_REWEIGHT_LNN_JSON

if run 1; then
hb "[1] S1 lnN at the chain4 working points"
T=$(date +%s)
( cd $CH/s1 && env -i HOME=$HOME PATH=/usr/bin:/bin S1_RW_TAG=fsrfix_dyveto $ALP_PY eval_alt_reweight.py --masses 1,2,3 --out yields_chain4.csv ) > $L/s1_eval.log 2>&1 \
  || die "eval_alt_reweight failed (s1_eval.log)"
grep -q "variants: \['nominal', 'zmhsb', 'zmnarrowhsb'\]" $L/s1_eval.log || die "S1 eval did not load the three dyveto reweights (s1_eval.log)"
( cd $CH/s1 && env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY make_s1_lnn_json.py --reference nominal $LNN yields_chain4.csv ) > $L/s1_json.log 2>&1 \
  || die "make_s1_lnn_json failed (s1_json.log)"
fresh $LNN $T || die "S1 JSON not written"
python3 - $LNN $WPJ <<'PY' || die "S1 JSON working points differ from the WP JSON"
import json, sys
j = json.load(open(sys.argv[1])); w = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(sys.argv[2]))["results"]}
v = j.get("values", j)
bad = [m for m in w if abs(float(v[str(m)]["mvacut"]) - w[m]) > 1e-9]
print("  S1 lnN:", " ".join("%s:%s/%s" % (m, v[str(m)].get("down", "?"), v[str(m)].get("up", "?")) for m in sorted(w)))
sys.exit(1 if bad else 0)
PY
fi

if false; then   # no interpolated point among m_a 1-3: interpolation JSON not needed
hb "[2] interpolation efficiency JSON"
T=$(date +%s)
( cd $CH/plot && env -i HOME=$HOME PATH=/usr/bin:/bin PYTHONPATH=$PLOT/lib $ALP_PY signal_eff_sumw.py ) > $L/signal_eff_sumw.log 2>&1
fresh $EFFJ_ELE $T && fresh $EFFJ_MU $T || die "efficiency JSON not written (signal_eff_sumw.log)"
python3 - $EFFJ_ELE $EFFJ_MU $PLOT/output/sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json $WPJ <<'PY' | tee -a $L/chain4_heartbeat.txt || die "efficiency JSON check failed"
import json, sys
e, mu, prod = (json.load(open(p)) for p in sys.argv[1:4])
for d, n in ((e, "ele"), (mu, "mu")):
    v = d["values"]; nz = sum(1 for y in v for m in range(1, 31) if float(v[y].get(str(m), 0)) > 0)
    if len(v) != 5 or nz != 150: print("  %s: %d years %d/150 non-zero" % (n, len(v), nz)); sys.exit(1)
if e["meta"].get("mva_cut_json") != sys.argv[4]: print("  meta.mva_cut_json", e["meta"].get("mva_cut_json")); sys.exit(1)
r = [float(e["values"][y][str(m)]) / float(prod["values"][y][str(m)]) for y in e["values"] for m in range(1, 31)]
print("    ele efficiency chain4/production: %.3f - %.3f" % (min(r), max(r)))
PY
fi

if run 3; then
hb "[3] MVA cut: data (30 mA) + signal (14 mA x 5 eras)"
T=$(date +%s)
for d in data sig; do [ -d $MV/root_MVAcut/$d ] && die "$MV/root_MVAcut/$d exists; move it away first"; done
# 10-14: a transient FUSE read error ("object tag too large") made hadd print the same ROOT error forever and the log
# filled the AFS work quota (16.8 GB). Each try is time-limited and its log capped at 50 MB; rc is taken from the job.
( ok=0; for a in 1 2 3; do
    ( timeout 5400 python3 $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py 2>&1; echo "CHAIN4_RC=$?" ) | head -c 50000000 > $L/mva_data_try$a.log
    grep -q '^CHAIN4_RC=0$' $L/mva_data_try$a.log && { ok=1; break; }; sleep 60; done; [ $ok = 1 ] ) &
pd=$!
bash $S/shellScripts/mva/run_apply_bdt_sig_6jobs.sh > $L/mva_sig.log 2>&1
wait $pd || die "apply_bdt_data failed 3 times (mva_data_try*.log)"
failed=$(grep -oP '^\[ERROR\] Failed \K(mA_M[0-9]+)_(\S+?)(?=;)' $L/mva_sig.log | sed -E 's/^(mA_M[0-9]+)_/\1:/' | tr '\n' ' ')
for job in $failed; do
  ok=0
  for att in 1 2 3; do
    hb "  rerun ${job%%:*} ${job##*:} (attempt $att)"
    python3 $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples ${job%%:*} --years ${job##*:} > $L/mva_sig_rerun_${job%%:*}_${job##*:}.log 2>&1 && { ok=1; break; }
    sleep 60
  done
  [ $ok -eq 1 ] || die "apply_bdt_sig ${job} failed 3 times"
done
bad=0
for m in $MASS30; do fresh $MV/root_MVAcut/data/mA_M$m/run3.root $T || { echo "  data mA$m stale/missing"; bad=1; }; done
python3 - $MV/root_MVAcut/data <<'PYC' || die "data run3.root entries differ from the sum of the era files"
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B = sys.argv[1]; bad = 0
def n(p):
    f = ROOT.TFile.Open(p); t = f.Get("DiphotonTree/Data_13p6TeV"); x = t.GetEntries() if t else -1; f.Close(); return x
for m in (1, 2, 3):
    s = sum(n("%s/mA_M%d/output_%s.root" % (B, m, e)) for e in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"))
    r = n("%s/mA_M%d/run3.root" % (B, m))
    if r != s: print("  mA%d run3 %d vs eras %d" % (m, r, s)); bad = 1
sys.exit(bad)
PYC
for m in $MASS14; do for y in $YEARS; do fresh $MV/root_MVAcut/sig/mA_M$m/output_$y.root $T || { echo "  sig mA$m $y stale/missing"; bad=1; }; done; done
[ $bad -eq 0 ] || die "MVA cut outputs incomplete"
python3 - $MV/root_MVAcut $LNN <<'PY' | tee $L/mva_audit.txt || die "MVA cut trees: wrong tree count, unreadable entries or missing lnN"
import json, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B = sys.argv[1]; J = json.load(open(sys.argv[2])); J = J.get("values", J); bad = 0
def walk(d, acc):
    for k in d.GetListOfKeys():
        o = k.ReadObj()
        if o.InheritsFrom("TDirectory"): walk(o, acc)
        elif o.InheritsFrom("TTree"):
            acc[0] += 1
            for i in range(o.GetEntries()):
                if o.GetEntry(i) <= 0: raise RuntimeError("%s entry %d" % (o.GetName(), i))
            acc[1] += o.GetEntries()
for m in (1, 2, 3):
    f = ROOT.TFile.Open("%s/data/mA_M%d/run3.root" % (B, m)); acc = [0, 0]; walk(f, acc); f.Close()
    print("  data mA%-2d trees %d entries %d" % (m, acc[0], acc[1]))
for m in (1, 2, 3):
    for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
        f = ROOT.TFile.Open("%s/sig/mA_M%d/output_%s.root" % (B, m, y)); acc = [0, 0]; walk(f, acc)
        t = f.Get("DiphotonTree/ggh_125_Za_ele_13p6TeV_cat0")
        lnn_ok = True
        if t and t.GetEntries():
            t.GetEntry(0); up, dn = t.weight_mva_reweight_Up, t.weight_mva_reweight_Down
            lnn_ok = abs(up - float(J[str(m)]["up"])) < 1e-6 and abs(dn - float(J[str(m)]["down"])) < 1e-6
        f.Close()
        if acc[0] != 34 or not lnn_ok:
            bad += 1; print("  BAD sig mA%d %s trees %d lnN %s" % (m, y, acc[0], lnn_ok))
print("  signal files: %d bad" % bad)
sys.exit(1 if bad else 0)
PY
fi

if run 4; then
hb "[4] Tree2WS"
T=$(date +%s)
( cd $S/Trees2WS && bash run_tree2ws.sh ) > $L/tree2ws.log 2>&1
bad=0
for m in $MASS30; do fresh $MV/root_MVAcut/data/mA_M$m/ws/run3.root $T || { echo "  data ws mA$m stale/missing"; bad=1; }; done
for m in $MASS14; do
  n=$(find $MV/root_MVAcut/sig/mA_M$m/ws_Tree2WS -name '*.root' -newermt "@$T" 2>/dev/null | wc -l)
  [ "$n" -ge 10 ] || { echo "  signal ws mA$m: only $n fresh files"; bad=1; }
done
[ $bad -eq 0 ] || die "Tree2WS outputs incomplete (tree2ws.log)"
# the anchors' data workspaces must hold exactly the events the closure scan fitted at the same cut
python3 - $SBE $WPJ $MV <<'PY' | tee -a $L/chain4_heartbeat.txt || die "data workspace counts differ from the closure-scan data"
import json, sys, re, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B = sys.argv[1]; bad = 0
for e in json.load(open(sys.argv[2]))["results"]:
    m, cut = int(e["mA"]), float(e["MVAcut"])
    f = ROOT.TFile.Open("%s/root_MVAcut/data/mA_M%d/ws/run3.root" % (sys.argv[3], m))
    n = f.Get("CMS_hza_workspace").data("Data_13p6TeV").sumEntries()
    log = "%s/closure_scan/cand_low_relE2_w123_drop_mA%d/cut_%s/data.build.log" % (B, m, ("%.5f" % cut).replace(".", ""))
    want = int(re.search(r"data N=(\d+)", open(log).read()).group(1))
    print("  mA%-2d data ws %d, closure-scan data %d %s" % (m, n, want, "OK" if n == want else "MISMATCH")); bad |= (n != want)
sys.exit(bad)
PY
fi

if run 5; then
hb "[5] background envelopes (30 mA, 6 in parallel, production binaries)"
T=$(date +%s)
bkg_one() {
  local m=$1 P=$FR/$1
  mkdir -p $P/HZAmassInde_fTest $FR/AllFitResults_blind $FR/AllFitResults_unblind
  ( cd $F/Background && ./bin/fTest_ALP_turnOn -i $MV/root_MVAcut/data/mA_M$m/ws/run3.root \
      --saveMultiPdf $P/CMS-HGG_mva_13p6TeV_multipdf.root -D $P/HZAmassInde_fTest \
      --mass_ALP $m -c 1 --isFlashgg 0 --isData 0 -f data, --datfile $P/fTest.dat \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $P/ftest.stdout.log 2>&1
  for b in blind unblind; do
    extra=""; [ $b = unblind ] && extra="--unblind"
    ( cd $F/Background && ./bin/makeBkgPlots_ALP -b $P/CMS-HGG_mva_13p6TeV_multipdf.root -d $P/BkgPlots_$b \
        --total_OutDir $FR/AllFitResults_$b -o $P/BkgPlots_$b.root --sqrts 13p6TeV --isMultiPdf --useBinnedData \
        --massStep 2.5 --mhVal 125.0 --maVal $m --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 \
        --intLumi 172.13 -c 0 --isFlashgg 0 --doBands $extra ) >> $P/makeBkgPlots.stdout.log 2>&1
  done
  echo "  [$(date '+%H:%M:%S')] bkg mA$m done"
}
export -f bkg_one; export F FR SBE MV
# slow points first (mA23 took ~98 min in production)
echo 23 4 1 2 3 $(seq 5 30 | grep -vx 23) | tr ' ' '\n' | xargs -P6 -I{} bash -c 'bkg_one {}' | tee -a $L/chain4_heartbeat.txt
python3 - $FR $T <<'PY' | tee $L/bkg_summary.txt || die "background envelopes missing, stale or empty"
import os, re, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
fr, T = sys.argv[1], int(sys.argv[2]); bad = 0
for m in (1, 2, 3):
    p = "%s/%d/CMS-HGG_mva_13p6TeV_multipdf.root" % (fr, m); n = -1
    if os.path.exists(p) and os.path.getmtime(p) >= T and os.path.getsize(p) >= 1024:
        f = ROOT.TFile.Open(p)
        for k in f.GetListOfKeys():
            o = k.ReadObj()
            if o.InheritsFrom("RooWorkspace"):
                it = o.allPdfs().createIterator(); pdf = it.Next()
                while pdf:
                    if pdf.InheritsFrom("RooMultiPdf"): n = max(n, pdf.getNumPdfs())
                    pdf = it.Next()
        f.Close()
    try:
        env = re.findall(r"pdf : ([A-Za-z0-9]+) , gof : ([0-9.e-]+)", open("%s/%d/HZAmassInde_fTest/EnvelopeResults.txt" % (fr, m)).read())
    except OSError:
        env = []
    gmin = min([float(g) for _, g in env]) if env else -1
    flag = "" if gmin > 0.01 else "  <-- GOF <= 0.01"
    print("  mA%-2d npdf %2d  %s%s" % (m, n, " ".join("%s:%s" % e for e in env), flag))
    bad += (n <= 0)
sys.exit(1 if bad else 0)
PY
fi

if run 6; then
hb "[6] signal model (14 mA, sig_sys 1-4, 6 masses in parallel)"
T=$(date +%s)
mkdir -p $CH/sig_single
sig_one() {
  local m=$1
  for s in 1_runjob_sig_fTest 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter; do
    sed -e "s/^mAs=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )$/mAs=( $m )/" $S/shellScripts/sig_sys/$s.sh > $CH/sig_single/${s}_m$m.sh
    grep -q "^mAs=( $m )$" $CH/sig_single/${s}_m$m.sh || { echo "  copy of $s for mA$m failed"; return 1; }
    bash $CH/sig_single/${s}_m$m.sh > $L/sig_${s}_m$m.log 2>&1
  done
  echo "  [$(date '+%H:%M:%S')] signal mA$m done"
}
export -f sig_one; export S CH L
echo $MASS14 | tr ' ' '\n' | xargs -P6 -I{} bash -c 'sig_one {}' | tee -a $L/chain4_heartbeat.txt
bad=0
for m in $MASS14; do for y in $YEARS; do for ch in ele mu; do
  fresh $S/Signal/outdir_$ch/signalFit/output/${m}_CMS-HGG_sigfit_${y}_${ch}_Hm125.root $T || { echo "  sigfit mA$m $y $ch stale/missing"; bad=1; }
done; done; done
[ $bad -eq 0 ] || die "signal fits incomplete (sig_*_m<m>.log)"
hb "    signal fits: 140/140 fresh"
fi

if run 7; then
hb "[7] datacards: wilson (with theory) and xsec (without QCDscale/pdf/alphaS)"
T=$(date +%s)
( cd $S/Datacard && sh 1_runjob_gen_datacard_makeYields.sh ) > $L/dc_1.log 2>&1
( cd $S/Datacard && sh 2_runjob_gen_datacard_makeDatacard.sh ) > $L/dc_2.log 2>&1
for m in $MASS30; do
  W=$DC/${m}_pruned_datacard_leptons.txt; X=$DC/${m}_pruned_datacard_leptons_xsec.txt
  fresh $W $T || die "datacard mA$m stale/missing (dc_*.log)"
  if awk '$1=="shapes"{for(i=1;i<=NF;i++) if($i ~ /\.root$/) print $i}' $W | grep -v "^$S/" | grep -q .; then
    die "datacard mA$m references shapes outside the sandbox"
  fi
  for th in QCDscale_ggH pdf_Higgs_ggH alphaS_ggH; do [ $(awk -v t=$th '$1==t' $W | wc -l) -eq 1 ] || die "mA$m: $th not exactly once in the wilson card"; done
  awk '$1!="QCDscale_ggH" && $1!="pdf_Higgs_ggH" && $1!="alphaS_ggH"' $W > $X
  [ $(( $(wc -l < $W) - $(wc -l < $X) )) -eq 3 ] || die "mA$m: xsec card does not differ by exactly 3 lines"
  case " $INTERP " in *" $m "*) [ $(grep -c "^interpolate_yield_.*rateParam" $W) -eq 10 ] || die "mA$m: interpolation rateParam missing";; esac
done
hb "    30 x 2 datacards OK (shapes in the sandbox, xsec = wilson minus 3 theory lines, rateParam at the 16 interpolated points)"
fi

if run 8; then
hb "[8] text2workspace + blind expected limits (wilson and xsec cards)"
T=$(date +%s)
mkdir -p $S/Combine/root_t2w_wilson $S/Combine/root_t2w_xsec $S/Combine/output_combine_results_wilson $S/Combine/output_combine_results_xsec
lim_one() {
  local m=$1 k=$2 card
  card=$DC/${m}_pruned_datacard_leptons.txt; [ $k = xsec ] && card=$DC/${m}_pruned_datacard_leptons_xsec.txt
  mkdir -p $S/Combine/work_$k/$m
  ( cd $S/Combine/work_$k/$m && text2workspace.py $card -o $S/Combine/root_t2w_$k/${m}_Datacard_leptons.root -m 125.38 higgsMassRange=115,135 > $L/t2w_${k}_m$m.log 2>&1
    combine -M AsymptoticLimits $S/Combine/root_t2w_$k/${m}_Datacard_leptons.root --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/limit_${k}_m$m.log 2>&1
    mv higgsCombine$m.AsymptoticLimits.mH125.38.root $S/Combine/output_combine_results_$k/ )
}
export -f lim_one; export S DC L
for m in $MASS30; do for k in wilson xsec; do echo "$m $k"; done; done | xargs -P6 -L1 bash -c 'lim_one $0 $1'
for m in $MASS30; do for k in wilson xsec; do fresh $S/Combine/output_combine_results_$k/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "limit $k mA$m missing (limit_${k}_m$m.log)"; done; done
python3 - $S/Combine $F/Combine/output_combine_results <<'PY' | tee $CH/limits_vs_production.txt
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def q(path):
    f = ROOT.TFile.Open(path); out = {}
    for e in f.Get("limit"): out[round(e.quantileExpected, 3)] = e.limit
    f.Close(); return out
print("blind expected limit on r (median [-1s, +1s]); production card = with theory uncertainties (wilson type)")
print("  mA |      production       |    chain4 wilson      |     chain4 xsec       | wilson/prod  xsec/prod")
for m in (1, 2, 3):
    p = q("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (sys.argv[2], m))
    w = q("%s/output_combine_results_wilson/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (sys.argv[1], m))
    x = q("%s/output_combine_results_xsec/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (sys.argv[1], m))
    f = lambda d: "%.4f [%.4f,%.4f]" % (d[0.5], d[0.16], d[0.84])
    print("  %2d | %s | %s | %s | %8.3f %9.3f" % (m, f(p), f(w), f(x), w[0.5] / p[0.5], x[0.5] / p[0.5]))
PY
fi

if run 9; then
ln -sfn output_combine_results_xsec $S/Combine/output_combine_results   # the bias script reads the limits there
hb "[9] closure with the sandbox xsec card"
# interpolated masses have no step-3 scan: run the same closure scan (DROP=1) at the cut they inherit
for m in $INTERP; do
  cut=$(python3 -c "
import json; r = {int(e['mA']): float(e['MVAcut']) for e in json.load(open('$WPJ'))['results']}
a = min(r, key=lambda k: (abs(k - $m), k)); print('%.5f' % r[a])")
  echo "$m $cut"
done > $L/closure_interp_cuts.txt
python3 - $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py $L/closure_interp_cuts.txt <<'PY' || die "nearest-anchor rule differs from apply_bdt_data.find_nearest_cut"
import importlib.util, sys
spec = importlib.util.spec_from_file_location("abd", sys.argv[1]); abd = importlib.util.module_from_spec(spec); spec.loader.exec_module(abd)
import json
for l in open(sys.argv[2]):
    m, cut = int(l.split()[0]), float(l.split()[1])
    cuts = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(abd.optimized_BDT_Cut))["results"]}
    a, c = abd.find_nearest_cut(m, cuts)
    assert abs(c - cut) < 1e-9, (m, a, c, cut)
print("  interpolated-mass cuts = apply_bdt_data nearest anchor")
PY
run_interp() { m=$1; cut=$2; [ -s $SBE/closure_scan/cand_low_relE2_w123_drop_mA$m/cut_$(echo $cut | tr -d .)/result.txt ] && return 0
  DROP=1 NPAR=1 bash $V/scan_closure.sh new $m "$cut" > $L/closure_scan_interp_mA$m.log 2>&1; }
export -f run_interp; export SBE V L
xargs -r -P6 -L1 bash -c 'run_interp $0 $1' < $L/closure_interp_cuts.txt
clo_one() {
  local m=$1 cut=$2 tag PSD O
  tag=$(echo $cut | tr -d .); PSD=$SBE/closure_scan/cand_low_relE2_w123_drop_mA$m/cut_$tag/ps/fit; O=$CH/closure_final/mA$m
  [ -s $PSD/CMS-HGG_mva_13p6TeV_multipdf.root ] || { echo "  mA$m: no pseudo-data fit at $PSD"; return 1; }
  mkdir -p $O
  sed "s#$FR/$m/#$PSD/#g" $DC/${m}_pruned_datacard_leptons_xsec.txt > $O/datacard.txt
  grep -q "$PSD/CMS-HGG_mva_13p6TeV_multipdf.root" $O/datacard.txt || { echo "  mA$m: background path not swapped"; return 1; }
  ( cd $O && text2workspace.py datacard.txt -o ws_t2w.root > t2w.log 2>&1
    combine -M Significance ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n PS$m > sig.log 2>&1
    combine -M MultiDimFit ws_t2w.root --algo singles --cminDefaultMinimizerStrategy 0 -m 125.38 --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n PS$m > rhat.log 2>&1
    rm -f ws_t2w.root )
  echo "  mA$m cut $cut closure Z=$(grep -E 'Significance:' $O/sig.log | sed -E 's/.*Significance: //')  r_hat=$(grep -E '^\s+r :' $O/rhat.log | sed -E 's/^\s+r :\s+//; s/ \(68%\)//')"
}
export -f clo_one; export FR DC CH
python3 -c "
import json
for e in json.load(open('$WPJ'))['results']: print(int(e['mA']), '%.5f' % float(e['MVAcut']))" > $L/closure_anchor_cuts.txt
cat $L/closure_anchor_cuts.txt $L/closure_interp_cuts.txt | sort -n | xargs -P6 -L1 bash -c 'clo_one $0 $1' | sort -t A -k2 -n > $CH/closure_final_summary.txt
cat $CH/closure_final_summary.txt | tee -a $L/chain4_heartbeat.txt
[ $(grep -c "closure Z=" $CH/closure_final_summary.txt) -eq 3 ] || die "closure incomplete ($CH/closure_final_summary.txt)"
fi

if run 10; then
hb "[10] bias (condor, 30 x 10 chunks x 1000 toys, xsec card)"
HS_OUT_BASE=$MV/bias_outputs_highstat ROOT_DATACARD_PATH=$S/Combine/root_t2w_xsec NCHUNKS=10 bash $S/shellScripts/bias/Condor/subjob_bias_highstat.sh $MASS30 > $L/bias_submit.log 2>&1
grep -E "submitted to cluster" $L/bias_submit.log | tee -a $L/chain4_heartbeat.txt || die "bias submission failed (bias_submit.log)"
fi

fingerprint > $L/prod_fingerprint.after
diff $L/prod_fingerprint.before $L/prod_fingerprint.after > $L/prod_fingerprint.diff || die "production products changed: see prod_fingerprint.diff"
timeout 1800 find $F /eos/home-p/pelai/HZa/root_MVAcut $PLOT/output -type f -newer $L/chain_start.marker > $L/prod_new_files.txt 2>/dev/null
[ -s $L/prod_new_files.txt ] && die "files written into production areas (prod_new_files.txt)"
hb "production products unchanged; no file in production areas newer than the chain start"
date '+%F %T' > $L/chain4.DONE
hb "CHAIN4 DONE"
