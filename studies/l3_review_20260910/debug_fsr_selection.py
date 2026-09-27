"""Run the analysis FSR-photon selection and assignment standalone on NanoAOD.

Purpose: find out why the production dresses almost no muon even though
n_fsr > 0 in 12-28% of muon-channel events.

Input : one signal NanoAODv15 file (xrootd)
Output: printed diagnostics (stdout)
"""
import sys

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA")
from higgs_dna.selections import photon_selections  # noqa: E402

PATH = ("root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/"
        "HZa-Zto2L-ato2G_Par-M-1_TuneCP5_13p6TeV_madgraph-pythia8/NANOAODSIM/"
        "150X_mcRun3_2024_realistic_v2-v2/110000/"
        "2394cc2e-fa67-4340-a2ce-cb9bdc55ab98.root")

OPTS = {"iso": 1.8, "eta": 2.4, "pt": 2.0, "dROverEt2": 0.012}

BR = ["FsrPhoton_pt", "FsrPhoton_eta", "FsrPhoton_phi", "FsrPhoton_relIso03",
      "FsrPhoton_dROverEt2",
      "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_looseId",
      "Muon_pfRelIso03_all",
      "Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass",
      "Photon_pt", "Photon_eta", "Photon_phi"]

with uproot.open(PATH + ":Events") as t:
    a = t.arrays(BR, entry_stop=60000, library="ak")


def coll(prefix, extra=()):
    d = {"pt": a[prefix + "_pt"], "eta": a[prefix + "_eta"],
         "phi": a[prefix + "_phi"]}
    d["mass"] = a[prefix + "_mass"] if prefix + "_mass" in a.fields \
        else ak.zeros_like(a[prefix + "_pt"])
    for e in extra:
        d[e] = a[prefix + "_" + e]
    return ak.Array(ak.zip(d), with_name="Momentum4D")


fsr = coll("FsrPhoton", ("relIso03", "dROverEt2"))
mu_all = coll("Muon", ("looseId", "pfRelIso03_all"))
ele = coll("Electron")
pho = coll("Photon")

mu = mu_all[(mu_all.pt > 5) & (abs(mu_all.eta) < 2.4)
            & (mu_all.looseId == 1) & (mu_all.pfRelIso03_all < 0.35)]
ele_sel = ele[(ele.pt > 7) & (abs(ele.eta) < 2.5)]
pho_sel = pho[(pho.pt > 10) & (abs(pho.eta) < 2.5)]

sel = photon_selections.select_resolved_fsr_photons(
    FSRphotons=fsr, electrons=ele_sel, muons=mu, photons=pho_sel,
    options=OPTS)

print("selection type          :", str(ak.type(sel)))
print("selection has None      :", bool(ak.any(ak.is_none(sel, axis=1))))
sel_fsr = fsr[sel]
print("selected-FSR type       :", str(ak.type(sel_fsr)))
n = ak.num(sel_fsr, axis=1)
print("events with n_fsr>0     : %d / %d (%.2f%%)"
      % (int(ak.sum(n > 0)), len(n), 100.0 * ak.sum(n > 0) / len(n)))
n_mu = ak.num(mu, axis=1)
both = (n > 0) & (n_mu >= 2)
print("  and >=2 selected muons: %d" % int(ak.sum(both)))

# how many of those have a selected FSR photon within dR<0.5 of a selected muon
sub_fsr = sel_fsr[both]
sub_mu = mu[both]
pairs = ak.cartesian({"lep": sub_mu, "ph": sub_fsr}, axis=1, nested=True)
dr = pairs.lep.deltaR(pairs.ph)
has = ak.any(ak.any(dr < 0.5, axis=2), axis=1)
print("  with dR(mu,gamma)<0.5 : %d (%.2f%% of the %d two-muon+FSR events)"
      % (int(ak.sum(has)), 100.0 * ak.sum(has) / max(int(ak.sum(both)), 1),
         int(ak.sum(both))))
