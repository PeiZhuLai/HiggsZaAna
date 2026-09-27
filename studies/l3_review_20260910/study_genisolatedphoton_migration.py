"""Compare the analysis' GenPart-based overlap-removal tag with GenIsolatedPhoton.

The L3 reviewer asked whether the DY+jets / DY+gamma split should be made with
the particle-level GenIsolatedPhoton collection instead of our GenPart-based
prompt-photon count. This builds the 2x2 migration matrix between the two
definitions so the size of the change can be judged before adopting it.

Current definition (higgs_dna/taggers/mc_overlap_tagger.py):
  GenPart, pdgId == 22, pt > 10 GeV, statusFlags isPrompt OR fromHardProcess,
  and dR > 0.05 from every prompt (fromHardProcess) GenPart with pt > 5 GeV
  and pdgId != 22.

Alternative:
  at least one GenIsolatedPhoton with pt > 10 GeV.

Input : NanoAODv15 files of DYto2Mu-2Jets and DYGto2LG (2024), via xrootd
Output: printed migration matrices (stdout) + a small text summary
"""
import sys
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

OUT = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910"
           "/genisolatedphoton_migration.txt")

SAMPLES = {
    "DY+jets (DYto2Mu-2Jets 2024)": [
        "root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/"
        "DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/"
        "150X_mcRun3_2024_realistic_v2-v6/100000/"
        "005c5bf4-d3d3-4fc1-9df3-99c83183c5f5.root",
    ],
    "DY+gamma (DYGto2LG 2024)": [
        "root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/"
        "DYGto2LG-1Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/"
        "150X_mcRun3_2024_realistic_v2-v2/100000/"
        "01dbd91c-e594-4ab6-9f75-baab9718be3d.root",
    ],
}

BRANCHES = ["GenPart_pdgId", "GenPart_pt", "GenPart_eta", "GenPart_phi",
            "GenPart_statusFlags", "GenIsolatedPhoton_pt"]

PT_MIN = 10.0
ISO_DR = 0.05


def delta_r(a, b):
    de = a.eta[:, :, None] - b.eta[:, None, :]
    dp = a.phi[:, :, None] - b.phi[:, None, :]
    dp = (dp + np.pi) % (2 * np.pi) - np.pi
    return np.sqrt(de ** 2 + dp ** 2)


def n_genpart_iso_photons(gp):
    """Reproduce mc_overlap_tagger.get_n_iso_photon exactly."""
    photon_cut = ((gp.pdgId == 22) & (gp.pt > PT_MIN)
                  & (((gp.statusFlags & 0x1) != 0)
                     | ((gp.statusFlags & 0x100) != 0)))
    photons = gp[photon_cut]
    truth_cut = (gp.pdgId != 22) & (gp.pt > 5) & ((gp.statusFlags & 0x100) != 0)
    truth = gp[truth_cut]
    keep = ak.fill_none(ak.all(delta_r(photons, truth) > ISO_DR, axis=-1), True)
    return ak.num(photons[keep], axis=1)


def run(label, paths, entries, lines):
    n_cur = []
    n_alt = []
    for path in paths:
        with uproot.open(path + ":Events") as tree:
            arr = tree.arrays(BRANCHES, entry_stop=entries, library="ak",
                              how="zip")
        n_cur.append(ak.to_numpy(n_genpart_iso_photons(arr.GenPart)))
        pt = arr["GenIsolatedPhoton"]["pt"]
        n_alt.append(ak.to_numpy(ak.num(pt[pt > PT_MIN], axis=1)))
    cur = np.concatenate(n_cur) > 0
    alt = np.concatenate(n_alt) > 0
    total = len(cur)

    both = int((cur & alt).sum())
    only_cur = int((cur & ~alt).sum())
    only_alt = int((~cur & alt).sum())
    neither = int((~cur & ~alt).sum())

    lines.append("")
    lines.append("%s  (%d events)" % (label, total))
    lines.append("                       GenIsolatedPhoton>=1   none")
    lines.append("  GenPart tag >=1      %8d (%5.2f%%)  %8d (%5.2f%%)"
                 % (both, 100.0 * both / total, only_cur,
                    100.0 * only_cur / total))
    lines.append("  GenPart tag none     %8d (%5.2f%%)  %8d (%5.2f%%)"
                 % (only_alt, 100.0 * only_alt / total, neither,
                    100.0 * neither / total))
    lines.append("  disagreement         %8d (%5.2f%%)"
                 % (only_cur + only_alt, 100.0 * (only_cur + only_alt) / total))
    for line in lines[-6:]:
        print(line, flush=True)


if __name__ == "__main__":
    n_ev = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    lines = ["GenPart-based vs GenIsolatedPhoton overlap-removal tag",
             "photon pT > %.0f GeV in both definitions" % PT_MIN]
    for label, paths in SAMPLES.items():
        run(label, paths, n_ev, lines)
    OUT.write_text("\n".join(lines) + "\n")
    print("\nwrote", OUT)
