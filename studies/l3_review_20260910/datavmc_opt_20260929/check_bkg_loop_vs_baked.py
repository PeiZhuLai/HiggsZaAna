import sys, numpy as np, ROOT
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from sideband_reweight import load_sideband_reweighter
RW = load_sideband_reweighter("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json", required=True)
f = ROOT.TFile.Open("/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/DYJetsToLL/2023preBPix.root"); t = f.Get("inclusive")
d = []
for i in range(min(3000, t.GetEntries())):
    t.GetEntry(i)
    d.append(t.weight * RW.weight_for_object(t, row_index=i) - t.weight_sideband_rwgt)
d = np.abs(np.array(d)); print("N", len(d), "max |loop - baked| =", d.max(), " rel max", (d / 1).max())
