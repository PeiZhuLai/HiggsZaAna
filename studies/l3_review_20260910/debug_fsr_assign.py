"""Replay the tagger's FSR sequence step by step to localize the lost dressing.

Mirrors za_tagger_resolved.py lines ~692-860 verbatim (selection -> with_field
mass -> Momentum4D -> costheta -> add_object_fields -> ak.num -> assign) and
prints the array type and the surviving counts after each step.

Input : one signal NanoAODv15 file (xrootd)
Output: printed diagnostics (stdout)
"""
import sys

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA")
from higgs_dna.selections import photon_selections  # noqa: E402
from higgs_dna.utils import awkward_utils  # noqa: E402

PATH = ("root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/"
        "HZa-Zto2L-ato2G_Par-M-1_TuneCP5_13p6TeV_madgraph-pythia8/NANOAODSIM/"
        "150X_mcRun3_2024_realistic_v2-v2/110000/"
        "2394cc2e-fa67-4340-a2ce-cb9bdc55ab98.root")
OPTS = {"iso": 1.8, "eta": 2.4, "pt": 2.0, "dROverEt2": 0.012}
DUMMY = -999.0

BR = ["FsrPhoton_pt", "FsrPhoton_eta", "FsrPhoton_phi", "FsrPhoton_relIso03",
      "FsrPhoton_dROverEt2", "FsrPhoton_muonIdx",
      "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_looseId",
      "Muon_pfRelIso03_all",
      "Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass",
      "Photon_pt", "Photon_eta", "Photon_phi"]

with uproot.open(PATH + ":Events") as t:
    a = t.arrays(BR, entry_stop=60000, library="ak")


def coll(prefix, extra=()):
    d = {"pt": a[prefix + "_pt"], "eta": a[prefix + "_eta"], "phi": a[prefix + "_phi"]}
    d["mass"] = a[prefix + "_mass"] if prefix + "_mass" in a.fields \
        else ak.zeros_like(a[prefix + "_pt"])
    for e in extra:
        d[e] = a[prefix + "_" + e]
    return ak.Array(ak.zip(d), with_name="Momentum4D")


fsr_raw = coll("FsrPhoton", ("relIso03", "dROverEt2", "muonIdx"))
mu_all = coll("Muon", ("looseId", "pfRelIso03_all"))
ele = coll("Electron")
pho = coll("Photon")

muons = mu_all[(mu_all.pt > 5) & (abs(mu_all.eta) < 2.4)
               & (mu_all.looseId == 1) & (mu_all.pfRelIso03_all < 0.35)]
electrons = ele[(ele.pt > 7) & (abs(ele.eta) < 2.5)]
photons = pho[(pho.pt > 10) & (abs(pho.eta) < 2.5)]

sel = photon_selections.select_resolved_fsr_photons(
    FSRphotons=fsr_raw, electrons=electrons, muons=muons, photons=photons,
    options=OPTS)
print("step1 selection type :", ak.type(sel))

FSRphotons = fsr_raw[sel]
print("step2 after mask     :", ak.type(FSRphotons),
      "| n>0:", int(ak.sum(ak.num(FSRphotons, axis=1) > 0)))

FSRphotons = ak.with_field(FSRphotons, ak.ones_like(FSRphotons.pt) * 0.0, "mass")
print("step3 after mass     :", ak.type(FSRphotons))
FSRphotons = ak.Array(FSRphotons, with_name="Momentum4D")
FSRphotons = ak.with_field(FSRphotons, FSRphotons.costheta, "costheta")
print("step5 after costheta :", ak.type(FSRphotons))

events = ak.Array({"dummy": np.zeros(len(FSRphotons))})
awkward_utils.add_object_fields(events=events, name="gamma_fsr",
                                objects=FSRphotons, n_objects=1,
                                dummy_value=DUMMY)
n_fsr = ak.num(FSRphotons, axis=1)
gpt = ak.to_numpy(events["gamma_fsr_pt"])
n_np = ak.to_numpy(n_fsr)
print("step6 n_fsr>0 = %d, gamma_fsr_pt real = %d, BOTH = %d"
      % (int((n_np > 0).sum()), int((gpt > -100).sum()),
         int(((n_np > 0) & (gpt > -100)).sum())))

# ---- assignment, verbatim from za_tagger_resolved.assign_fsr_photon ----
zero_ph = ak.zip({"pt": 0.0, "eta": 0.0, "phi": 0.0, "mass": 0.0,
                  "dROverEt2": np.inf}, with_name="Momentum4D")
fsr_padded = ak.fill_none(ak.pad_none(FSRphotons, 1, clip=False), zero_ph)
print("step7 fsr_padded     :", ak.type(fsr_padded))

pairs = ak.cartesian({"lep": muons, "ph": fsr_padded}, axis=1, nested=True)
dR = pairs.lep.deltaR(pairs.ph)
valid = dR < 0.5
has_valid = ak.any(valid, axis=2)
n_mu = ak.num(muons, axis=1)
base = (n_np > 0) & (ak.to_numpy(n_mu) >= 2)
print("step8 events with >=2 mu and n_fsr>0 : %d" % int(base.sum()))
print("      of those, >=1 lepton has_valid : %d"
      % int(ak.sum(ak.any(has_valid, axis=1)[base])))
