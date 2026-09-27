"""Cross-check the FSR-recovery rate directly on NanoAOD.

Applies the AN FSR-photon criteria (pt>2, |eta|<2.4, relIso03<1.8,
dROverEt2<0.012) and requires dR<0.5 to a muon inside the analysis
acceptance, to see how often the recovery should fire.

Input : one signal NanoAODv15 file (xrootd)
Output: printed counts (stdout)
"""
import awkward as ak
import numpy as np
import uproot

PATH = ("root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/"
        "HZa-Zto2L-ato2G_Par-M-1_TuneCP5_13p6TeV_madgraph-pythia8/NANOAODSIM/"
        "150X_mcRun3_2024_realistic_v2-v2/110000/"
        "2394cc2e-fa67-4340-a2ce-cb9bdc55ab98.root")

BR = ["FsrPhoton_pt", "FsrPhoton_eta", "FsrPhoton_phi", "FsrPhoton_relIso03",
      "FsrPhoton_dROverEt2", "FsrPhoton_muonIdx",
      "Muon_pt", "Muon_eta", "Muon_phi", "Muon_looseId", "Muon_pfRelIso03_all"]


def dr(e1, p1, e2, p2):
    de = e1[:, :, None] - e2[:, None, :]
    dp = (p1[:, :, None] - p2[:, None, :] + np.pi) % (2 * np.pi) - np.pi
    return np.sqrt(de ** 2 + dp ** 2)


with uproot.open(PATH + ":Events") as t:
    a = t.arrays(BR, entry_stop=200000, library="ak")

n_ev = len(a["FsrPhoton_pt"])
mu_sel = ((a["Muon_pt"] > 5) & (abs(a["Muon_eta"]) < 2.4)
          & (a["Muon_looseId"] == 1) & (a["Muon_pfRelIso03_all"] < 0.35))
n_mu = ak.num(a["Muon_pt"][mu_sel], axis=1)

fsr_sel = ((a["FsrPhoton_pt"] > 2) & (abs(a["FsrPhoton_eta"]) < 2.4)
           & (a["FsrPhoton_relIso03"] < 1.8) & (a["FsrPhoton_dROverEt2"] < 0.012))

print("events read                      : %d" % n_ev)
print("events with >=2 selected muons   : %d" % int(ak.sum(n_mu >= 2)))
print("events with >=1 raw FsrPhoton    : %d (%.2f%%)"
      % (int(ak.sum(ak.num(a["FsrPhoton_pt"], axis=1) > 0)),
         100.0 * ak.sum(ak.num(a["FsrPhoton_pt"], axis=1) > 0) / n_ev))
print("events with >=1 FsrPhoton passing kinematic+iso+dROverEt2: %d (%.2f%%)"
      % (int(ak.sum(ak.num(a["FsrPhoton_pt"][fsr_sel], axis=1) > 0)),
         100.0 * ak.sum(ak.num(a["FsrPhoton_pt"][fsr_sel], axis=1) > 0) / n_ev))

# now require dR < 0.5 to a selected muon
base = (n_mu >= 2) & (ak.num(a["FsrPhoton_pt"][fsr_sel], axis=1) > 0)
sub = {k: a[k][base] for k in BR}
fs = fsr_sel[base]
ms = mu_sel[base]
d = dr(ak.to_numpy(ak.fill_none(ak.pad_none(sub["FsrPhoton_eta"][fs], 4, clip=True), 999.0)),
       ak.to_numpy(ak.fill_none(ak.pad_none(sub["FsrPhoton_phi"][fs], 4, clip=True), 999.0)),
       ak.to_numpy(ak.fill_none(ak.pad_none(sub["Muon_eta"][ms], 4, clip=True), -999.0)),
       ak.to_numpy(ak.fill_none(ak.pad_none(sub["Muon_phi"][ms], 4, clip=True), -999.0)))
match = (d < 0.5).any(axis=(1, 2))
print("  of those, >=1 pairing with dR(mu,gamma)<0.5: %d (%.2f%% of all events)"
      % (int(match.sum()), 100.0 * match.sum() / n_ev))
print("  muonIdx of passing FSR photons  :",
      np.unique(ak.to_numpy(ak.flatten(sub["FsrPhoton_muonIdx"][fs])), return_counts=True))
