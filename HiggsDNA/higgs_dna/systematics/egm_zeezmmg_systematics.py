"""
Two-stage (Zee + Zmmg) photon energy scale and smearing for Run 3.

Background
----------
The default photon scale and smearing in this framework
(``photon_systematics.photon_scale_smear_run3``) reads the central EGM
``photonSS_EtDependent.json`` payloads.  Those are derived from Z->ee only:
the electron-based calibration is transferred to photons, and the residual
electron/photon difference is covered by a hard-coded flat scale uncertainty.

The two-stage calibration adds a second step measured directly on photons with
Z->mumu+gamma (FSR photons), on top of the Z->ee stage.  Both stages live in a
single correctionlib payload,
``JSONs/scaleAndSmearing/Hgg/EGMScalesSmearing_ZeeZmmg_<RunID>.v1.json``:

  compound ``Scale_Zee``  : stage-1 scale, Z->ee only        (data)
  compound ``Scale_Zmmg`` : stage-2 residual scale, Z->mumu+gamma
  compound ``Scale``      : the product of both stages       (data)
  ``SmearAndSyst_Zee``    : smearing + its systematics       (MC)
  ``SmearAndSyst_Zmmg``   : residual photon scale systematic (MC)
  ``Smear1G``             : single-Gaussian smearing, used for sigma_E/E
  ``EGMRandomGenerator``  : deterministic random numbers keyed on
                            (event, seediEtaOriX, seediPhiOriY)

Run 3 uses a double-Gaussian smearing; Run 2 would use ``Smear1G``.  Only the
Run 3 eras are shipped here, because the Run 2 payloads are derived for
NanoAODv15 and this analysis reads NanoAODv9 for Run 2 (feeding v9 into them
would double count the scale and smearing already applied at the NanoAOD step).

Branch contract
---------------
The function is a drop-in replacement for ``photon_scale_smear_run3`` and
writes exactly the same branches, so ``metadata/*.json``
``independent_collections`` needs no change:

  data : Photon.pt, Photon.energyErr, Photon.corrected_pt,
         Photon.corrected_energyErr
  MC   : the same, plus
         Photon.dEsigmaUp / dEsigmaDown                (additive, smearing)
         Photon.energyErr_dEsigmaUp / _dEsigmaDown
         Photon.pt_ScaleUp / pt_ScaleDown              (absolute, scale)
         Photon.energyErr_ScaleUp / _ScaleDown

``pt_ScaleUp/Down`` combines the Zee and the Zmmg scale systematics in
quadrature so that the existing single ``Photon_scale`` nuisance keeps working.
The two stages are additionally written out separately as
``pt_ScaleZeeUp/Down`` and ``pt_ScaleZmmgUp/Down`` (plus the matching
``energyErr_*``), so that they can be declared as two decorrelated nuisances
later without rerunning: the Zmmg stage is correlated between 2022 and 2023,
which is why the payloads for those eras are a combined 2022+2023 fit.

Required NanoAOD branches
-------------------------
On top of what the default path needs, the deterministic smearing requires
``Photon_seediEtaOriX`` and ``Photon_seediPhiOriY``, and ``event``.  If they are
absent the function refuses to run and returns the events untouched, so that a
missing branch can never silently degrade into an uncorrected sample.
"""

import os

import awkward
import numpy

import correctionlib

from higgs_dna.utils import awkward_utils
from higgs_dna.utils.logger_utils import simple_logger

logger = simple_logger(__name__)


# --------------------------------------------------------------------------
# Payloads
# --------------------------------------------------------------------------
# RunID naming follows the upstream Hgg ingredients: the 2022 and 2023 eras are
# a combined "2223" fit (the Zmmg stage is correlated across them), 2024 and
# 2025 are standalone.
ZEEZMMG_FILE = {
    "2022preEE": "EGMScalesSmearing_ZeeZmmg_22232022preEE.v1.json",
    "2022postEE": "EGMScalesSmearing_ZeeZmmg_22232022postEE.v1.json",
    "2023preBPix": "EGMScalesSmearing_ZeeZmmg_22232023preBPIX.v1.json",
    "2023postBPix": "EGMScalesSmearing_ZeeZmmg_22232023postBPIX.v1.json",
    "2024": "EGMScalesSmearing_ZeeZmmg_2024.v1.json",
    "2025": "EGMScalesSmearing_ZeeZmmg_2025.v1.json",
}

# Double Gaussian for Run 3; a Run 2 era would use "1G" (Smear1G).
GAUSSIANS = {
    "2022preEE": "2G",
    "2022postEE": "2G",
    "2023preBPix": "2G",
    "2023postBPix": "2G",
    "2024": "2G",
    "2025": "2G",
}

# Extra photon scale systematics folded into the Zmmg stage, following the
# upstream implementation: the muon momentum scale limits how well the Zmmg
# scale can be measured, and the ECAL non-linearity is not constrained by the
# FSR photons, which are soft.
MUON_SCALE_SYST = 0.0005
NONLINEARITY_SYST_EB = 0.0015
NONLINEARITY_SYST_EE = 0.0025
NONLINEARITY_PT_THRESHOLD = 80.0

# The 2024 and 2025 calibrations have a known non-closure for photons at high
# eta; the smearing uncertainty is inflated there until it is understood.
HIGH_EE_ETA = 2.1
HIGH_EE_SMEAR_INFLATION = 2.0
HIGH_EE_INFLATED_YEARS = {"2024", "2025"}
# The reference implementation masks on the *signed* supercluster eta, so only
# the +eta endcap gets the inflated uncertainty.  That asymmetry has no physical
# motivation and looks like an upstream oversight, but it is reproduced here so
# that this analysis stays numerically identical to the reference.  Set this to
# True to inflate both endcaps instead.
HIGH_EE_SYMMETRIC = False

PHOTON_PT_MIN = 10.0
MAX_ABS_ETA = 3.0

_JSON_SUBDIR = os.path.join("JSONs", "scaleAndSmearing", "Hgg")


def zeezmmg_supported(year):
    """True if a two-stage payload is shipped for this era."""
    return year in ZEEZMMG_FILE


def zeezmmg_json_path(year):
    return os.path.join(os.path.dirname(__file__), _JSON_SUBDIR, ZEEZMMG_FILE[year])


def double_smearing(std_normal, std_flat, mu, sigma, sigma_scale, frac):
    """
    Double Gaussian smearing factor.

    ``sigma_scale`` is the width of the tail Gaussian relative to the core one,
    ``frac`` selects between the two, and ``mu`` shifts the tail. This is the
    convention used by the Zee+Zmmg payloads (upstream ``old_convention=False``).
    """
    scale_core = 1.0 + sigma * std_normal
    scale_tail = mu * (1.0 + sigma_scale * sigma * std_normal)
    return numpy.where(std_flat > frac, scale_tail, scale_core)


def _out_of_acceptance(abs_eta, pt):
    """Photons the payloads do not cover: corrections are left at unity there."""
    return (abs_eta > MAX_ABS_ETA) | (pt < PHOTON_PT_MIN)


def _energy_err(energy_err, pt, sc_eta, rho, factor):
    """sigma_E propagated through a relative smearing rho and a pt factor."""
    return numpy.sqrt(energy_err ** 2 + (pt * numpy.cosh(sc_eta) * rho) ** 2) * factor


def photon_scale_smear_zeezmmg_run3(events, year, is_data):
    """
    Apply the two-stage (Zee + Zmmg) photon scale and smearing.

    :param events: awkward array of events, modified in place
    :param year: era tag, e.g. "2023preBPix"
    :param is_data: True for data (scale only), False for MC (smearing + systematics)
    :return: the events, with the branches listed in the module docstring
    """
    if not zeezmmg_supported(year):
        logger.warning(
            "[ZeeZmmg SaS] No two-stage payload for year %s, photons are left uncorrected. "
            "Supported years: %s", year, sorted(ZEEZMMG_FILE)
        )
        return events

    required_fields = [
        ("Photon", "eta"), ("Photon", "pt"), ("Photon", "r9"),
        ("Photon", "energyErr"), ("Photon", "isScEtaEB"), ("Photon", "isScEtaEE"),
        ("Photon", "seedGain"), ("run"),
    ]
    if not is_data:
        # The smearing random numbers are keyed on the seed crystal, so that a
        # given photon always gets the same random number regardless of how the
        # sample is split into chunks.
        required_fields += [
            ("Photon", "seediEtaOriX"), ("Photon", "seediPhiOriY"), ("event"),
        ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        # Deliberately fatal rather than a warning: the two-stage method is
        # opt-in, so silently handing back uncorrected photons would leave a
        # miscalibrated sample that looks like a successful run.  Add the
        # branches to the "branches" list of the analysis config.
        raise ValueError(
            "[ZeeZmmg SaS] Cannot apply the two-stage photon scale and smearing for %s (%s): "
            "the following fields are missing from the input: %s"
            % (year, "data" if is_data else "MC", str(missing_fields))
        )

    path_json = zeezmmg_json_path(year)
    try:
        cset = correctionlib.CorrectionSet.from_file(path_json)
    except OSError:
        logger.error("[ZeeZmmg SaS] Could not open the payload %s. Exiting.", path_json)
        raise

    logger.info(
        "[ZeeZmmg SaS] Applying two-stage (Zee+Zmmg) photon scale and smearing for %s "
        "(%s, %s) from %s",
        year, "data" if is_data else "MC", GAUSSIANS[year], os.path.basename(path_json),
    )

    photons = events["Photon"]
    n_photons = awkward.num(photons.pt)
    flat = awkward.flatten(photons)

    # For photons the NanoAOD eta is already the supercluster eta.
    sc_eta = awkward.to_numpy(flat.eta)
    abs_sc_eta = numpy.abs(sc_eta)
    pt_raw = awkward.to_numpy(flat.pt)
    r9 = awkward.to_numpy(flat.r9)
    energy_err = awkward.to_numpy(flat.energyErr)
    seed_gain = awkward.to_numpy(flat.seedGain)
    is_eb = awkward.to_numpy(flat.isScEtaEB)
    run = awkward.to_numpy(
        awkward.flatten(awkward.broadcast_arrays(events["run"], photons.pt)[0])
    )

    skip = _out_of_acceptance(abs_sc_eta, pt_raw)
    ones = numpy.ones_like(pt_raw, dtype=float)
    zeros = numpy.zeros_like(pt_raw, dtype=float)

    # ---------------------------------------------------------------- data --
    if is_data:
        # Both stages at once: the compound "Scale" is Scale_Zee * Scale_Zmmg.
        scale = numpy.where(
            skip, ones,
            cset.compound["Scale"].evaluate("scale", run, sc_eta, r9, pt_raw, seed_gain),
        )
        corrected_pt = pt_raw * scale

        # sigma_E/E needs the single-Gaussian width even in the 2G eras.
        rho = numpy.where(
            skip, zeros,
            cset["Smear1G"].evaluate("smear", corrected_pt, r9, sc_eta),
        )
        corrected_energy_err = _energy_err(energy_err, pt_raw, sc_eta, rho, scale)

        events["Photon", "corrected_pt"] = awkward.unflatten(corrected_pt, n_photons)
        events["Photon", "corrected_energyErr"] = awkward.unflatten(corrected_energy_err, n_photons)
        events["Photon", "pt"] = awkward.unflatten(corrected_pt, n_photons)
        events["Photon", "energyErr"] = awkward.unflatten(corrected_energy_err, n_photons)
        return events

    # ------------------------------------------------------------------ MC --
    event_number = awkward.to_numpy(
        awkward.flatten(awkward.broadcast_arrays(events["event"], photons.pt)[0])
    )
    seed_ieta = awkward.to_numpy(awkward.flatten(photons.seediEtaOriX))
    seed_iphi = awkward.to_numpy(awkward.flatten(photons.seediPhiOriY))

    rng = cset["EGMRandomGenerator"]
    std_normal = rng.evaluate("stdnormal", event_number, seed_ieta, seed_iphi)
    std_flat = rng.evaluate("stdflat", event_number, seed_ieta, seed_iphi)

    smear_zee = cset["SmearAndSyst_Zee"]
    two_gaussian = GAUSSIANS[year] == "2G"

    # Shape parameters of the tail Gaussian: they do not change between the
    # nominal and the varied smearing, so evaluate them once.
    if two_gaussian:
        tail_mu = smear_zee.evaluate("mu", pt_raw, r9, sc_eta)
        tail_reso_scale = smear_zee.evaluate("reso_scale", pt_raw, r9, sc_eta)
        tail_frac = smear_zee.evaluate("frac", pt_raw, r9, sc_eta)

    def _smear_factor(sigma):
        """Turn a smearing width into the multiplicative pt factor."""
        if not two_gaussian:
            return 1.0 + sigma * std_normal
        return double_smearing(
            std_normal, std_flat, tail_mu, sigma, tail_reso_scale, tail_frac
        )

    sigma_nominal = numpy.where(
        skip, zeros, smear_zee.evaluate("smear", pt_raw, r9, sc_eta)
    )
    correction = numpy.where(skip, ones, _smear_factor(sigma_nominal))
    corrected_pt = pt_raw * correction

    # sigma_E/E always from the single-Gaussian width, as recommended upstream.
    rho = numpy.where(skip, zeros, cset["Smear1G"].evaluate("smear", pt_raw, r9, sc_eta))
    corrected_energy_err = _energy_err(energy_err, pt_raw, sc_eta, rho, correction)

    events["Photon", "corrected_pt"] = awkward.unflatten(corrected_pt, n_photons)
    events["Photon", "corrected_energyErr"] = awkward.unflatten(corrected_energy_err, n_photons)

    # --- smearing systematics -------------------------------------------
    if two_gaussian and year in HIGH_EE_INFLATED_YEARS:
        esmear = smear_zee.evaluate("esmear", pt_raw, r9, sc_eta)
        high_ee = numpy.abs(sc_eta) > HIGH_EE_ETA if HIGH_EE_SYMMETRIC \
            else sc_eta > HIGH_EE_ETA
        inflation = numpy.where(high_ee, HIGH_EE_SMEAR_INFLATION, 1.0)
        sigma_up = sigma_nominal + inflation * esmear
        sigma_down = sigma_nominal - inflation * esmear
    else:
        sigma_up = smear_zee.evaluate("smear_up", pt_raw, r9, sc_eta)
        sigma_down = smear_zee.evaluate("smear_down", pt_raw, r9, sc_eta)
    sigma_down = numpy.maximum(0.0, sigma_down)

    for label, sigma in (("Up", sigma_up), ("Down", sigma_down)):
        sigma = numpy.where(skip, zeros, sigma)
        factor = numpy.where(skip, ones, _smear_factor(sigma))
        pt_varied = pt_raw * factor
        # The framework consumes the smearing variation as an additive shift.
        events["Photon", "dEsigma" + label] = awkward.unflatten(
            pt_varied - corrected_pt, n_photons
        )
        events["Photon", "energyErr_dEsigma" + label] = awkward.unflatten(
            _energy_err(energy_err, pt_raw, sc_eta, sigma, factor), n_photons
        )

    # --- scale systematics ----------------------------------------------
    # Stage 1: Z->ee.
    escale_zee = numpy.where(
        skip, zeros, smear_zee.evaluate("escale", pt_raw, r9, sc_eta)
    )

    # Stage 2: Z->mumu+gamma. Note the input order is (ScEta, r9, pt) in this
    # payload, the transpose of SmearAndSyst_Zee.
    escale_zmmg = numpy.where(
        skip, zeros,
        cset["SmearAndSyst_Zmmg"].evaluate("escale", sc_eta, r9, pt_raw),
    )
    escale_zmmg = numpy.sqrt(escale_zmmg ** 2 + MUON_SCALE_SYST ** 2)
    nonlinearity = numpy.where(is_eb, NONLINEARITY_SYST_EB, NONLINEARITY_SYST_EE)
    # The threshold is evaluated on the smeared pt, matching the reference
    # implementation (the uncertainty itself is computed from the raw pt).
    nonlinearity = numpy.where(corrected_pt > NONLINEARITY_PT_THRESHOLD, nonlinearity, 0.0)
    escale_zmmg = numpy.where(
        skip, zeros, numpy.sqrt(escale_zmmg ** 2 + nonlinearity ** 2)
    )

    # Combined, for the single Photon_scale nuisance the metadata declares.
    escale_total = numpy.sqrt(escale_zee ** 2 + escale_zmmg ** 2)

    # Scale variations act on the smeared pt.  "Scale" is the combined nuisance
    # the metadata already declares; the per-stage ones are extra outputs.
    for tag, escale in (
        ("Scale", escale_total),
        ("ScaleZee", escale_zee),
        ("ScaleZmmg", escale_zmmg),
    ):
        for label, sign in (("Up", 1.0), ("Down", -1.0)):
            factor = 1.0 + sign * escale
            events["Photon", "pt_%s%s" % (tag, label)] = awkward.unflatten(
                corrected_pt * factor, n_photons
            )
            events["Photon", "energyErr_%s%s" % (tag, label)] = awkward.unflatten(
                corrected_energy_err * factor, n_photons
            )

    events["Photon", "pt"] = awkward.unflatten(corrected_pt, n_photons)
    events["Photon", "energyErr"] = awkward.unflatten(corrected_energy_err, n_photons)

    return events
