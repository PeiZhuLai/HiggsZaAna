"""Locate which FSR-selection term introduces None into the mask.

The tagger's FSR-photon mask comes back as ?bool, so events.FsrPhoton[mask]
is option-typed: ak.num counts the None placeholders (n_fsr > 0 in 12-19 % of
muon events) while almost none of them carry a real photon, and the dressing
never happens. This script re-implements select_resolved_fsr_photons term by
term inside a real tagger run and prints the None count of each term.

Input : fsr_rerun/config.json (one signal NanoAODv15 file)
Output: [TERM] lines on stdout
"""
import json
import sys

import awkward as ak

from higgs_dna.utils.logger_utils import setup_logger
from higgs_dna.selections import object_selections, photon_selections
from higgs_dna.taggers import za_tagger_resolved

setup_logger("INFO")
_done = {"n": 0}


def probe_select(FSRphotons, electrons, muons, photons, options,
                 name="none", tagger=None):
    terms = {}
    terms["pt"] = FSRphotons.pt > options["pt"]
    terms["eta"] = abs(FSRphotons.eta) < options["eta"]
    terms["iso"] = FSRphotons.relIso03 < options["iso"]
    terms["dROverEt2"] = FSRphotons.dROverEt2 < options["dROverEt2"]
    terms["clean_ele"] = object_selections.delta_R(FSRphotons, electrons, 0.001)
    far_e = object_selections.delta_R(FSRphotons, electrons, 0.5)
    far_m = object_selections.delta_R(FSRphotons, muons, 0.5)
    terms["far_from_electrons"] = far_e
    terms["far_from_muons"] = far_m
    terms["lep_indR0p5"] = ~(far_e & far_m)
    photons_sorted = photons[ak.argsort(photons.pt, ascending=False)]
    lead_photons = photons_sorted[:, :1]
    terms["clean_lead_photon"] = object_selections.delta_R(
        FSRphotons, lead_photons, 0.2)

    if _done["n"] == 0:
        for label, coll in (("electrons", electrons), ("muons", muons),
                            ("photons", photons)):
            print("[TERM] n %-10s: %d events with 0, type=%s"
                  % (label, int(ak.sum(ak.num(coll, axis=1) == 0)),
                     ak.type(coll)), flush=True)
        for key, val in terms.items():
            try:
                n_none = int(ak.sum(ak.is_none(val, axis=1)))
            except Exception as exc:  # noqa: BLE001
                n_none = "<%s>" % exc
            print("[TERM] %-20s type=%s  None entries=%s"
                  % (key, ak.type(val), n_none), flush=True)
        _done["n"] = 1

    return (terms["pt"] & terms["eta"] & terms["iso"] & terms["dROverEt2"]
            & terms["clean_ele"] & terms["lep_indR0p5"]
            & terms["clean_lead_photon"])


photon_selections.select_resolved_fsr_photons = probe_select
za_tagger_resolved.photon_selections = photon_selections

from higgs_dna.analysis import run_analysis  # noqa: E402

run_analysis(json.load(open(sys.argv[1])))
