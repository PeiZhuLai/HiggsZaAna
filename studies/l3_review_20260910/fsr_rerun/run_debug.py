"""Run the tagger on one file with probes inside the FSR selection/assignment."""
import json
import sys

import awkward as ak
import numpy as np

from higgs_dna.utils.logger_utils import setup_logger
from higgs_dna.taggers.za_tagger_resolved import ZaTaggerRun3

setup_logger("INFO")


def _describe(tag, arr):
    try:
        print(f"[PROBE] {tag}: type={ak.type(arr)}", flush=True)
    except Exception as exc:
        print(f"[PROBE] {tag}: <no type> {exc}", flush=True)


_orig_sel = ZaTaggerRun3.select_FSRphotons
_orig_assign = ZaTaggerRun3.assign_fsr_photon


def sel_probe(self, FSRphotons, electrons, muons, photons, options):
    out = _orig_sel(self, FSRphotons, electrons, muons, photons, options)
    _describe("select_FSRphotons INPUT FsrPhoton", FSRphotons)
    _describe("select_FSRphotons MASK", out)
    try:
        print("[PROBE] mask has None:", bool(ak.any(ak.is_none(out, axis=1))),
              "| n_true events:", int(ak.sum(ak.sum(out, axis=1) > 0)), flush=True)
    except Exception as exc:
        print("[PROBE] mask None check failed:", exc, flush=True)
    return out


def assign_probe(self, leptons, fsr_photons):
    _describe("assign_fsr_photon leptons", leptons)
    _describe("assign_fsr_photon fsr_photons", fsr_photons)
    n = ak.num(fsr_photons, axis=1)
    print("[PROBE] assign: events with n_fsr>0 =", int(ak.sum(n > 0)), flush=True)
    try:
        first = ak.firsts(ak.pad_none(fsr_photons, 1, clip=True))
        print("[PROBE] assign: first photon is None for",
              int(ak.sum(ak.is_none(first))), "of", len(first), flush=True)
    except Exception as exc:
        print("[PROBE] assign firsts failed:", exc, flush=True)
    out = _orig_assign(self, leptons, fsr_photons)
    try:
        moved = ak.sum(ak.sum(np.abs(out.pt - leptons.pt) > 1e-4, axis=1) > 0)
        print("[PROBE] assign: events with >=1 dressed lepton =", int(moved), flush=True)
    except Exception as exc:
        print("[PROBE] assign moved check failed:", exc, flush=True)
    return out


ZaTaggerRun3.select_FSRphotons = sel_probe
ZaTaggerRun3.assign_fsr_photon = assign_probe

from higgs_dna.analysis import run_analysis  # noqa: E402
run_analysis(json.load(open(sys.argv[1])))
