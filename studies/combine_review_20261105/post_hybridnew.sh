#!/bin/bash
# Expected CLs limits from the HybridNew grid (toys) and their ratio to the asymptotic expected limits.
# For each mass: hadd the grid points, read the expected quantiles with --readHybridResults, and
# compare with output_combine_results_xsec. Writes hybridnew_vs_asymptotic.txt next to this script.
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
SUB=${1:-hybridnew}
HN=$FF/Combine/Checks/review_20261105/$SUB
HNA=$FF/Combine/Checks/review_20261105/hybridnew
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` 2>/dev/null
cd $HN || exit 2
for m in 5 20 30; do
  echo "mA$m grid files: $(ls mA$m/higgsCombine_mA${m}_r*.HybridNew.*.root | wc -l)"
  hadd -f grid_mA$m.root mA$m/higgsCombine_mA${m}_r*.HybridNew.*.root > /dev/null
  for q in 0.025 0.16 0.5 0.84 0.975; do
    combine -M HybridNew $FF/Combine/root_t2w_xsec/${m}_Datacard_leptons.root -m 125.38 \
      -D $HNA/higgsCombine_asimov_mA${m}.GenerateOnly.mH125.38.123456.root:toys/toy_asimov \
      --setParameters MH=125.38 --freezeParameters MH --LHCmode LHC-limits \
      --readHybridResults --grid grid_mA$m.root --expectedFromGrid $q -n _read_mA${m}_q$q > read_mA${m}_q$q.log 2>&1
  done
done
python3 - $FF $HN <<'PY' | tee $W/${SUB}_vs_asymptotic.txt
import glob, sys, ROOT
ff, hn = sys.argv[1:3]
qs = (0.025, 0.16, 0.5, 0.84, 0.975)
print("mA  quantile  toys(HybridNew)       asymptotic  ratio")
for m in (5, 20, 30):
    fa = ROOT.TFile.Open("%s/Combine/output_combine_results_xsec/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (ff, m))
    asy = [e.limit for e in fa.Get("limit")]
    for i, q in enumerate(qs):
        fs = glob.glob("%s/higgsCombine_read_mA%d_q%s.HybridNew.mH125.38*.root" % (hn, m, q))
        fr = ROOT.TFile.Open(fs[0]) if fs else None
        t = fr.Get("limit") if fr else None
        lim = [e.limit for e in t][0] if t and t.GetEntries() else float("nan")
        err = [e.limitErr for e in t][0] if t and t.GetEntries() else float("nan")
        print("%-3d %-8s  %.5f +- %.5f  %.5f     %.3f" % (m, q, lim, err, asy[i], lim / asy[i]))
PY
