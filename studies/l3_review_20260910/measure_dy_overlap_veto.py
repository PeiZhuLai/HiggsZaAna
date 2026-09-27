"""Measure the DY+jets / DY+gamma overlap-removal veto fraction.

Reproduces higgs_dna/taggers/mc_overlap_tagger.py:get_n_iso_photon exactly and
reports the fraction of events removed from each sample, to answer the L3
question on how DY+jets and DY+gamma are separated.

Input : one NanoAODv15 file per sample (xrootd)
Output: printed table (stdout)
"""
import sys
import awkward as ak
import numpy as np
import uproot

BRANCHES = ["GenPart_pdgId", "GenPart_pt", "GenPart_eta", "GenPart_phi",
            "GenPart_statusFlags"]


def delta_r(a, b):
    da = a.eta[:, :, None] - b.eta[:, None, :]
    dp = a.phi[:, :, None] - b.phi[:, None, :]
    dp = (dp + np.pi) % (2 * np.pi) - np.pi
    return np.sqrt(da ** 2 + dp ** 2)


def n_iso_photons(gp, pt_thresh=10.0, iso_dr=0.05):
    iso_cut = ((gp.pdgId == 22) & (gp.pt > pt_thresh)
               & (((gp.statusFlags & 0x1) != 0) | ((gp.statusFlags & 0x100) != 0)))
    photons = gp[iso_cut]
    truth_cut = (gp.pdgId != 22) & (gp.pt > 5) & ((gp.statusFlags & 0x100) != 0)
    truth = gp[truth_cut]
    dr = delta_r(photons, truth)
    keep = ak.all(dr > iso_dr, axis=-1)
    keep = ak.fill_none(keep, True)
    return ak.num(photons[keep], axis=1)


def run(label, path, entries):
    with uproot.open(path + ":Events") as t:
        arr = t.arrays(BRANCHES, entry_stop=entries, library="ak", how="zip")
    n = n_iso_photons(arr.GenPart)
    tot = len(n)
    with_iso = int(ak.sum(n > 0))
    print(f"{label:<22} events={tot:>8}  with >=1 iso gen photon={with_iso:>8} "
          f"({100.0 * with_iso / tot:5.2f}%)  without={tot - with_iso:>8} "
          f"({100.0 * (tot - with_iso) / tot:5.2f}%)")


if __name__ == "__main__":
    n_ev = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    for label, path in [
        ("DY+jets (DYto2Mu 2024)",
         "root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/"
         "DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/"
         "150X_mcRun3_2024_realistic_v2-v6/100000/"
         "005c5bf4-d3d3-4fc1-9df3-99c83183c5f5.root"),
        ("DY+gamma (DYGto2LG 2024)",
         "root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/"
         "DYGto2LG-1Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/"
         "150X_mcRun3_2024_realistic_v2-v2/100000/"
         "01dbd91c-e594-4ab6-9f75-baab9718be3d.root"),
    ]:
        run(label, path, n_ev)
