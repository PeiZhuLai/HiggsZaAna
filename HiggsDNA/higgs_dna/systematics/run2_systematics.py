"""Run 2 UL (NanoAODv9) corrections and uncertainties that have no Run 3 counterpart in HZa.

- l1_prefiring_sf: L1 ECAL/muon pre-firing weight stored in NanoAODv9
  (L1PreFiringWeight_{Nom,Up,Dn}); event weight, as pu_reweight_sf.
- egm_scale_unc_run2_{electron,photon}: energy-scale uncertainty from the EGM UL payload
  (data/<era>_UL/EGM_ScaleUnc.json, correction "UL-EGM_ScaleUnc", multiplicative up/down
  factors as a function of supercluster eta and seed gain). In NanoAODv9 the scale is
  applied to data only and Electron/Photon_dEscaleUp/Down are zero in MC, so they cannot
  be used for the MC variation.
- egm_smear_unc_run2_{electron,photon}: smearing uncertainty from the NanoAODv9 branches
  dEsigmaUp/Down. In CMSSW 10_6 (PhysicsTools/NanoAOD/python/{electrons,photons}_cff.py)
  they are energy differences in GeV, dEsigmaUp = E(post-correction) - E(sigma up), so the
  varied transverse momentum is pt * (1 - dEsigma / E) with E = pt * cosh(eta).

The Run 2 scale and smearing corrections themselves are already applied in NanoAODv9, so
these functions only provide up/down variations and never return a nominal one.
"""

import os

import awkward
import numpy
import correctionlib

import logging
logger = logging.getLogger(__name__)

from higgs_dna.utils import awkward_utils

RUN2_YEARS = ("2016preVFP", "2016postVFP", "2017", "2018")
EGM_SCALE_UNC_FILE = os.path.join(os.path.dirname(__file__), "data", "%s_UL", "EGM_ScaleUnc.json")


def _check_run2(year, name):
    if year not in RUN2_YEARS:
        raise ValueError("[run2_systematics : %s] only for Run 2 UL years %s, got %s" % (name, RUN2_YEARS, year))


def l1_prefiring_sf(events, year, central_only):
    _check_run2(year, "l1_prefiring_sf")
    required = ["L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Dn"]
    missing = awkward_utils.missing_fields(events, required)
    if missing:
        raise ValueError("[run2_systematics : l1_prefiring_sf] missing fields %s" % missing)
    variations = {"central": awkward.values_astype(events.L1PreFiringWeight_Nom, numpy.float64)}
    if not central_only:
        variations["up"] = awkward.values_astype(events.L1PreFiringWeight_Up, numpy.float64)
        variations["down"] = awkward.values_astype(events.L1PreFiringWeight_Dn, numpy.float64)
    return variations


def _egm_scale_unc(objects, year, sc_eta, nominal_only, modify_nominal, name):
    _check_run2(year, name)
    if nominal_only:
        raise ValueError("[run2_systematics : %s] provides variations only; the Run 2 scale is "
                         "already applied in NanoAODv9" % name)
    # No "nominal" key is returned, so the nominal collection is never modified.
    counts = awkward.num(objects.pt)
    pt = awkward.to_numpy(awkward.flatten(objects.pt))
    eta = awkward.to_numpy(awkward.flatten(sc_eta))
    gain = awkward.to_numpy(awkward.flatten(objects.seedGain)).astype(numpy.int32)
    unknown = ~numpy.isin(gain, [1, 6, 12])
    if unknown.any():
        # The payload only has gains 12/6/1; treat anything else as gain 12 and say so.
        logger.warning("[run2_systematics : %s] %d objects with seedGain not in {1,6,12}; using 12",
                       name, int(unknown.sum()))
        gain = numpy.where(unknown, 12, gain)
    corr = correctionlib.CorrectionSet.from_file(EGM_SCALE_UNC_FILE % year)["UL-EGM_ScaleUnc"]
    up = corr.evaluate(year, "scaleup", eta, gain)
    down = corr.evaluate(year, "scaledown", eta, gain)
    return {"up": awkward.unflatten(pt * up, counts), "down": awkward.unflatten(pt * down, counts)}


def egm_scale_unc_run2_electron(events, year, nominal_only, modify_nominal):
    ele = events.Electron
    return _egm_scale_unc(ele, year, ele.eta + ele.deltaEtaSC, nominal_only, modify_nominal,
                          "egm_scale_unc_run2_electron")


def egm_scale_unc_run2_photon(events, year, nominal_only, modify_nominal):
    # NanoAODv9 photons have no supercluster eta; their eta is computed at the supercluster
    # position, so it is used directly.
    pho = events.Photon
    return _egm_scale_unc(pho, year, pho.eta, nominal_only, modify_nominal, "egm_scale_unc_run2_photon")


def _egm_smear_unc(objects, year, nominal_only, modify_nominal, name):
    _check_run2(year, name)
    if nominal_only:
        raise ValueError("[run2_systematics : %s] provides variations only; the Run 2 smearing is "
                         "already applied in NanoAODv9" % name)
    energy = objects.pt * numpy.cosh(objects.eta)
    safe = energy > 0
    return {
        "up": awkward.where(safe, objects.pt * (1.0 - objects.dEsigmaUp / energy), objects.pt),
        "down": awkward.where(safe, objects.pt * (1.0 - objects.dEsigmaDown / energy), objects.pt),
    }


def egm_smear_unc_run2_electron(events, year, nominal_only, modify_nominal):
    return _egm_smear_unc(events.Electron, year, nominal_only, modify_nominal, "egm_smear_unc_run2_electron")


def egm_smear_unc_run2_photon(events, year, nominal_only, modify_nominal):
    return _egm_smear_unc(events.Photon, year, nominal_only, modify_nominal, "egm_smear_unc_run2_photon")
