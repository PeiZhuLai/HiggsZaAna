import sys, numpy as np, uproot, ROOT, types
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from sideband_reweight import load_sideband_reweighter
RW = load_sideband_reweighter("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json", required=True)
# pull the three helpers out of 1_prepare_dataVmc.py without running its argparse
src = open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts/1_prepare_dataVmc.py").read()
seg = src[src.index("class _TrueMassParam:"):src.index("def get_sideband_reweight_uncertainty_weights")]
ns = {"SIDEBAND_REWEIGHTER": RW}; exec(seg, ns)
cfg = types.SimpleNamespace(sig_names=["M5"])
P = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M5/2024.root"
N = 2000
f = ROOT.TFile.Open(P); t = f.Get("inclusive")
loop_w, hash_w = [], []
for i in range(N):
    t.GetEntry(i)
    loop_w.append(RW.weight_for_object(ns["_with_true_mass_param"](t, "M5", cfg), row_index=i))
    hash_w.append(RW.weight_for_object(t, row_index=i))
a = uproot.open(P)["inclusive"].arrays(entry_stop=N, library="pd")
a = a.rename(columns={"pho1PIso_noCorr": "pho1ECALIso", "pho2PIso_noCorr": "pho2ECALIso"})
a["param"] = (a["ALP_m"] - 5.0) / a["H_m"]
df_w = RW.weights_for_dataframe(a)
loop_w, hash_w = np.array(loop_w), np.array(hash_w)
print("max |loop(true mass) - dataframe(true mass)| =", np.max(np.abs(loop_w - df_w)))
print("mean weight: true-mass %.4f  event-hash %.4f  (N=%d)" % (loop_w.mean(), hash_w.mean(), N))
