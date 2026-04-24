import awkward
import numpy

from correctionlib import _core
import correctionlib

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.systematics.utils import systematic_from_bins, ic_systematic_from_bins

def _clib_inputs(obj):
    """Return list of input names for correctionlib highlevel objects (Correction/CompoundCorrection)."""
    return [v.name for v in getattr(obj, "inputs", [])]

def _clib_eval(obj, args_by_name, override_order=None):
    """
    Evaluate a correctionlib object using its declared input order.
    args_by_name: dict(name -> numpy array/scalar)
    override_order: optional list of names to force a specific order.
    """
    order = override_order if override_order is not None else _clib_inputs(obj)
    missing = [k for k in order if k not in args_by_name]
    if missing:
        raise ValueError(f"Missing inputs for correctionlib evaluate: {missing} (needed order={order})")
    return obj.evaluate(*[args_by_name[k] for k in order])

def _pick_2024_mc_scale_correction(evaluator):
    """
    For 2024 MC, avoid compound['Scale'] because it requires run/seedGain.
    Try to find a non-compound scale correction that only needs photon-level inputs.
    """
    # Highlevel CorrectionSet supports dict-like access by key
    for name in getattr(evaluator, "keys", lambda: [])():
        try:
            corr = evaluator[name]
            ins = _clib_inputs(corr)
        except Exception:
            continue

        # Heuristic: must accept syst and pt, and must NOT require run/seedGain (data-only)
        if ("syst" in ins) and ("pt" in ins) and ("run" not in ins) and ("seedGain" not in ins):
            # Prefer names that look like scale
            if "Scale" in name or "scale" in name:
                return name
    # fallback: any correction matching heuristic
    for name in getattr(evaluator, "keys", lambda: [])():
        try:
            corr = evaluator[name]
            ins = _clib_inputs(corr)
        except Exception:
            continue
        if ("syst" in ins) and ("pt" in ins) and ("run" not in ins) and ("seedGain" not in ins):
            return name
    return None

def _pick_2024_mc_scale_unc_correction(evaluator, scale_corr_name, required_inputs):
    """
    Find a matching uncertainty correction for the chosen 2024 MC scale correction.
    Heuristic: name contains 'unc'/'Unc'/'err' and inputs match required_inputs exactly.
    """
    for name in getattr(evaluator, "keys", lambda: [])():
        if name == scale_corr_name:
            continue
        if not (("unc" in name.lower()) or ("err" in name.lower())):
            continue
        try:
            corr = evaluator[name]
            ins = _clib_inputs(corr)
        except Exception:
            continue
        if ins == required_inputs:
            return name
    return None

########################
##### Photon ID SF #####
########################

PHOTON_ID_SF_FILE = {
    "2016" : ["higgs_dna/systematics/data/2016postVFP_UL/photon.json", "higgs_dna/systematics/data/2016preVFP_UL/hzg_phidvalidate_2016_scalefactors.json"],
    "2016preVFP" : ["higgs_dna/systematics/data/2016preVFP_UL/photon.json", "higgs_dna/systematics/data/2016preVFP_UL/hzg_phidvalidate_2016APV_scalefactors.json"],
    "2016postVFP" : ["higgs_dna/systematics/data/2016postVFP_UL/photon.json", "higgs_dna/systematics/data/2016postVFP_UL/hzg_phidvalidate_2016_scalefactors.json"],
    "2017" : ["higgs_dna/systematics/data/2017_UL/photon.json", "higgs_dna/systematics/data/2017_UL/hzg_phidvalidate_2017_scalefactors.json"],
    "2018" : ["higgs_dna/systematics/data/2018_UL/photon.json", "higgs_dna/systematics/data/2018_UL/hzg_phidvalidate_2018_scalefactors.json"],
    "2022preEE" : ["higgs_dna/systematics/data/2022preEE_UL/photon.json", "higgs_dna/systematics/data/2022preEE_UL/hzg_phidvalidate_2022_scalefactors.json"],
    "2022postEE" : ["higgs_dna/systematics/data/2022postEE_UL/photon.json", "higgs_dna/systematics/data/2022postEE_UL/hzg_phidvalidate_2022EE_scalefactors.json"],
    "2023preBPix" : ["higgs_dna/systematics/data/2023preBPix_UL/photon.json"],
    "2023postBPix" : ["higgs_dna/systematics/data/2023postBPix_UL/photon.json"],
    "2024" : ["higgs_dna/systematics/data/2023postBPix_UL/photon.json"], # FIXME
}

PHOTON_ID_SF = {
    "2016" : "2016postVFP",
    "2016preVFP" : "2016preVFP",
    "2016postVFP" : "2016postVFP",
    "2017" : "2017",
    "2018" : "2018",
    "2022preEE" : "2022Re-recoBCD",
    "2022postEE" : "2022Re-recoE+PromptFG",
    "2023preBPix" : "2023PromptC",
    "2023postBPix" : "2023PromptD",
    "2024" : "2023PromptD" # FIXME
}

PHOTON_ID_EVAL = {
    "2016" : "UL-Photon-ID-SF",
    "2016preVFP" : "UL-Photon-ID-SF",
    "2016postVFP" : "UL-Photon-ID-SF",
    "2017" : "UL-Photon-ID-SF",
    "2018" : "UL-Photon-ID-SF",
    "2022preEE" : "Photon-ID-SF",
    "2022postEE" : "Photon-ID-SF",
    "2023preBPix" : "Photon-ID-SF",
    "2023postBPix" : "Photon-ID-SF",
    "2024" : "Photon-ID-SF" # FIXME
}

def photon_id_sf(events, year, central_only, input_collection, working_point = "wp80"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "phi")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluators = []
    for file in PHOTON_ID_SF_FILE[year]:
        evaluators.append(_core.CorrectionSet.from_file(misc_utils.expand_path(file)))

    photons = events[input_collection]

    # Flatten photons then convert to numpy for compatibility with correctionlib
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    variations_list = []

    evaluator = evaluators[0]

    pho_eta = numpy.clip(
        awkward.to_numpy(photons_flattened.eta),
        -2.49999,
        2.49999 # SFs only valid up to eta 2.5
    )

    pho_pt = numpy.clip(
        awkward.to_numpy(photons_flattened.pt),
        20.0, # SFs only valid for pT >= 20.0
        499.999 # and pT < 500.
    )

    pho_phi = numpy.clip(
        awkward.to_numpy(photons_flattened.phi),
        -999.0,
        999.0 # SFs only valid for phi in [-pi, pi]
    )

    variations = {}
    if int(year[:4]) < 2023:
        sf = evaluator[PHOTON_ID_EVAL[year]].evalv(
                PHOTON_ID_SF[year],
                "sf",
                working_point,
                pho_eta,
                pho_pt
        )
    else:
        # For 2023, need to handle pt < 20 photons separately
        pho_pt_orig = awkward.to_numpy(photons_flattened.pt)
        pt_below_20_mask = pho_pt_orig < 20.0
        
        # Create separate pt arrays for below/above 20 GeV
        pho_pt_below_20 = numpy.clip(pho_pt_orig, 15.0, 19.999)
        pho_pt_above_20 = numpy.clip(pho_pt_orig, 20.0, 499.999)
        
        # Calculate SF for pt >= 20
        sf_above_20 = evaluator[PHOTON_ID_EVAL[year]].evalv(
                PHOTON_ID_SF[year],
                "sf",
                working_point,
                pho_eta,
                pho_pt_above_20,
                pho_phi
        )
        
        # Calculate SF for pt < 20 with modified working point
        sf_below_20 = evaluator[PHOTON_ID_EVAL[year]].evalv(
                PHOTON_ID_SF[year],
                "sf",
                f"{working_point}Below20",
                pho_eta,
                pho_pt_below_20,
                pho_phi
        )
        
        # Combine results based on pt threshold
        sf = numpy.where(pt_below_20_mask, sf_below_20, sf_above_20)
    variations["central"] = awkward.unflatten(sf, n_photons)

    if not central_only:
        syst_vars = ["sfup", "sfdown"] 
        for syst_var in syst_vars:
            if int(year[:4]) < 2023:
                syst = evaluator[PHOTON_ID_EVAL[year]].evalv(
                        PHOTON_ID_SF[year],
                        syst_var,
                        working_point,
                        pho_eta,
                        pho_pt
                )
            else:
                # For 2023, need to handle pt < 20 photons separately
                pho_pt_orig = awkward.to_numpy(photons_flattened.pt)
                pt_below_20_mask = pho_pt_orig < 20.0
                
                # Create separate pt arrays for below/above 20 GeV
                pho_pt_below_20 = numpy.clip(pho_pt_orig, 15.0, 19.999)
                pho_pt_above_20 = numpy.clip(pho_pt_orig, 20.0, 499.999)
                
                # Calculate syst for pt >= 20
                syst_above_20 = evaluator[PHOTON_ID_EVAL[year]].evalv(
                        PHOTON_ID_SF[year],
                        syst_var,
                        working_point,
                        pho_eta,
                        pho_pt_above_20,
                        pho_phi
                )
                
                # Calculate syst for pt < 20 with modified working point
                syst_below_20 = evaluator[PHOTON_ID_EVAL[year]].evalv(
                        PHOTON_ID_SF[year],
                        syst_var,
                        f"{working_point}Below20",
                        pho_eta,
                        pho_pt_below_20,
                        pho_phi
                )
                
                # Combine results based on pt threshold
                syst = numpy.where(pt_below_20_mask, syst_below_20, syst_above_20)
                
            if "up" in syst_var:
                syst_var_name = "up"
            elif "down" in syst_var:
                syst_var_name = "down"
            variations[syst_var_name] = awkward.unflatten(syst, n_photons)
    variations_list.append(variations)

    if len(evaluators) == 2:
        evaluator = evaluators[1]

        pho_eta = numpy.clip(
            awkward.to_numpy(photons_flattened.eta),
            -2.49999,
            2.49999 # SFs only valid up to eta 2.5
        )

        pho_pt = numpy.clip(
            awkward.to_numpy(photons_flattened.pt),
            15.0, # SFs only valid for pT >= 15.0
            20.0,
        )

        variations = {}
        sf = evaluator["sf_pass"].evalv(
                pho_pt,
                pho_eta
        )
        variations["central"] = awkward.unflatten(sf, n_photons)

        if not central_only:
            syst = evaluator["unc_pass"].evalv(
                    pho_pt,
                    pho_eta
            )
            variations["up"] = awkward.unflatten(sf + syst, n_photons)
            variations["down"] = awkward.unflatten(sf - syst, n_photons)
        variations_list.append(variations)

    if len(variations_list) > 1:
        for var in variations_list[0].keys():
            variations[var] = awkward.where(
                (photons.pt < 15.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
                awkward.ones_like(variations_list[0][var], dtype=float),
                awkward.where(photons.pt < 20.0, variations_list[1][var], variations_list[0][var])
            )
    else:
        for var in variations_list[0].keys():
            variations[var] = awkward.where(
                (photons.pt < 15.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
                awkward.ones_like(variations_list[0][var], dtype=float),
                variations_list[0][var]
            )

    return variations


############################
##### HZa Photon ID SF #####
############################

PHOTON_ZAID_SF_FILE = {
    "2022preEE" : ["higgs_dna/systematics/data/2022preEE_UL/hza_phid_2022_scalefactors.json"],
    "2022postEE" : ["higgs_dna/systematics/data/2022postEE_UL/hza_phid_2022EE_scalefactors.json"],
    "2023preBPix" : ["higgs_dna/systematics/data/2023preBPix_UL/hza_phid_2023_scalefactors.json"],
    "2023postBPix" : ["higgs_dna/systematics/data/2023postBPix_UL/hza_phid_2023BPix_scalefactors.json"],
    "2024" : ["higgs_dna/systematics/data/2024_UL/hza_phid_2024_scalefactors.json"],
}

PHOTON_ZAID_SF_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hza_phid_2023BPixHole_scalefactors.json"
}


def photon_zaid_sf(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """

    required_fields = [
        (input_collection, "eta"), (input_collection, "pt")
    ]
    if year == "2023postBPix":
        required_fields.append((input_collection, "phi"))

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluators = []
    for file in PHOTON_ZAID_SF_FILE[year]:
        evaluators.append(_core.CorrectionSet.from_file(misc_utils.expand_path(file)))
    evaluator_hole = None
    if year == "2023postBPix":
        evaluator_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(PHOTON_ZAID_SF_FILE_HOLE[year]))

    photons = events[input_collection]

    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    pho_eta = numpy.clip(
        awkward.to_numpy(photons_flattened.eta),
        -2.49999,
        2.49999,
    )
    pho_pt_orig = awkward.to_numpy(photons_flattened.pt)
    pho_pt = numpy.clip(pho_pt_orig, 10.0, 499.999)

    is_hole_region = None
    if year == "2023postBPix":
        pho_phi = awkward.to_numpy(photons_flattened.phi)
        is_hole_region = (pho_eta > -1.566) & (pho_eta < 0.0) & (pho_phi > -1.2) & (pho_phi < -0.8)

    def _eval_pass(evaluator, var):
        if var == "central":
            return evaluator["sf_pass"].evalv(pho_pt, pho_eta)
        sf0 = evaluator["sf_pass"].evalv(pho_pt, pho_eta)
        unc = evaluator["unc_pass"].evalv(pho_pt, pho_eta)
        if var == "up":
            return sf0 + unc
        if var == "down":
            return sf0 - unc
        raise ValueError(f"Unknown var={var}")

    variations = {}
    for var in (["central"] if central_only else ["central", "up", "down"]):
        sf = _eval_pass(evaluators[0], var)

        if evaluator_hole is not None:
            sf_hole = _eval_pass(evaluator_hole, var)
            sf = numpy.where(is_hole_region, sf_hole, sf)

        variations[var] = awkward.unflatten(sf, n_photons)

    for var in variations.keys():
        variations[var] = awkward.where(
            (photons.pt < 10.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
            awkward.ones_like(variations[var], dtype=float),
            variations[var],
        )

    return variations

#################################
### Custom Photon e-veto SF #####
#################################

PHOTON_CSEV_SF_FILE = {
    "2022preEE" : ["higgs_dna/systematics/data/2022preEE_UL/hza_phcsev_2022_scalefactors.json"],
    "2022postEE" : ["higgs_dna/systematics/data/2022postEE_UL/hza_phcsev_2022EE_scalefactors.json"],
    "2023preBPix" : ["higgs_dna/systematics/data/2023preBPix_UL/hza_phcsev_2023_scalefactors.json"],
    "2023postBPix" : ["higgs_dna/systematics/data/2023postBPix_UL/hza_phcsev_2023BPix_scalefactors.json"],
    "2024" : ["higgs_dna/systematics/data/2024_UL/hza_phcsev_2024_scalefactors.json"],
}

PHOTON_CSEV_SF_FILE_HOLE = {
    "2023postBPix" : "higgs_dna/systematics/data/2023postBPix_UL/hza_phcsev_2023BPixHole_scalefactors.json"
}

def photon_CSEV_sf(events, year, central_only, input_collection, working_point = "MVA"):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """
    
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "r9")
    ]
    if year == "2023postBPix":
        required_fields.append((input_collection, "phi"))

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(PHOTON_CSEV_SF_FILE[year][0]))
    evaluator_hole = None
    if year == "2023postBPix":
        evaluator_hole = _core.CorrectionSet.from_file(misc_utils.expand_path(PHOTON_CSEV_SF_FILE_HOLE[year]))

    photons = events[input_collection]

    # Flatten photons then convert to numpy for compatibility with correctionlib
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    pho_eta = numpy.clip(
        awkward.to_numpy(photons_flattened.eta),
        -2.49999,
        2.49999 # SFs only valid up to eta 2.5
    )
    pho_abs_eta = numpy.abs(pho_eta)
    pho_pt = numpy.clip(
        awkward.to_numpy(photons_flattened.pt),
        10.0,
        499.999
    )

    pho_r9 = numpy.clip(
        awkward.to_numpy(photons_flattened.r9),
        0.0, # SFs only valid for R9 >= 0.0
        99999.0
    )

    is_hole_region = None
    if year == "2023postBPix":
        pho_phi = awkward.to_numpy(photons_flattened.phi)
        is_hole_region = (pho_eta > -1.566) & (pho_eta < 0.0) & (pho_phi > -1.2) & (pho_phi < -0.8)

    # Calculate SF and syst
    variations = {}

    sf = evaluator["sf_pass"].evalv(pho_abs_eta, pho_pt, pho_r9)
    if evaluator_hole is not None:
        sf_hole = evaluator_hole["sf_pass"].evalv(pho_abs_eta, pho_pt, pho_r9)
        sf = numpy.where(is_hole_region, sf_hole, sf)
    variations["central"] = awkward.unflatten(sf, n_photons)

    if not central_only:
        syst = evaluator["unc_pass"].evalv(pho_abs_eta, pho_pt, pho_r9)
        if evaluator_hole is not None:
            syst_hole = evaluator_hole["unc_pass"].evalv(pho_abs_eta, pho_pt, pho_r9)
            syst = numpy.where(is_hole_region, syst_hole, syst)
        variations["up"] = awkward.unflatten(sf + syst, n_photons)
        variations["down"] = awkward.unflatten(sf - syst, n_photons)

    for var in variations.keys():
        # Set SFs = 1 for photons which are not applicable
        variations[var] = awkward.where(
                (photons.pt < 10.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
                awkward.ones_like(variations[var], dtype=float),
                variations[var]
        )

    return variations

########################
##### Fake Photon SF ###
########################

PHOTON_Fake_Photon_SF_FILE = {
    "2016preVFP" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2016postVFP" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2017" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2018" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2022preEE" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2022postEE" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2023preBPix" : "higgs_dna/systematics/data/fakephoton/fakephoton.json",
    "2023postBPix" : "higgs_dna/systematics/data/fakephoton/fakephoton.json"
}

def photon_fake_photon_sf(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/electronExample.py
    """
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "genPartFlav")
    ]
    
    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        logger.warning(f"Missing fields for photon_fake_photon_sf: {missing_fields}")

    is_run2 = not (year.startswith("2022") or year.startswith("2023"))
    run_str = "run2" if is_run2 else "run3"
    
    photons = events[input_collection]
    n_photons = awkward.num(photons)
    
    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(PHOTON_Fake_Photon_SF_FILE[year]))
    photons_flattened = awkward.flatten(photons)

    photons_abseta = numpy.clip(
            numpy.abs(awkward.to_numpy(photons_flattened.eta)),
            0.0,
            2.4999
    )
    photons_pt = numpy.clip(
            awkward.to_numpy(photons_flattened.pt),
            15.0,
            499.999
    )

    isjet = awkward.where(photons_flattened.genPartFlav == 1, 0, 1)    
    photons_isjet = awkward.to_numpy(isjet).astype(int)
    
    sf = evaluator["fakephoton_corrections"].evalv(
        run_str,
        photons_isjet,
        photons_pt,
        photons_abseta
    )

    variations = {}
    variations["central"] = awkward.unflatten(sf, n_photons)

    if not central_only:
        # No systematics for fake photon SF, so up/down are the same as central
        variations["up"] = variations["central"]
        variations["down"] = variations["central"]

    for var in variations.keys():
        variations[var] = awkward.where(
            (photons.pt < 15.0) | (photons.pt >= 500.0) | (abs(photons.eta) >= 2.5),
            awkward.ones_like(variations[var], dtype=float),
            variations[var]
        )

    return variations

##################
### Trigger SF ###
##################
# Note: since the trigger sf applies separate SF for the lead/sublead photons,
# it's easiest to just cast this as an EventWeightSystematic (rather than ObjectWeightSystematic as we would typically do)
# and just multiply the lead/sublead variations manually by hand

from higgs_dna.systematics.data.trigger_sf import LEAD_TRIGGER_SF_2016, SUBLEAD_TRIGGER_SF_2016, LEAD_TRIGGER_SF_2017, SUBLEAD_TRIGGER_SF_2017, LEAD_TRIGGER_SF_2018, SUBLEAD_TRIGGER_SF_2018 
lead_trigger_sf_bins = {
    "2016" : LEAD_TRIGGER_SF_2016,
    "2016preVFP" : LEAD_TRIGGER_SF_2016,
    "2016postVFP" : LEAD_TRIGGER_SF_2016,
    "2017" : LEAD_TRIGGER_SF_2017,
    "2018" : LEAD_TRIGGER_SF_2018,
    "2022preEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : LEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : LEAD_TRIGGER_SF_2018 #FIXME
}
sublead_trigger_sf_bins = {
    "2016" : SUBLEAD_TRIGGER_SF_2016,
    "2016preVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2016postVFP" : SUBLEAD_TRIGGER_SF_2016,
    "2017" : SUBLEAD_TRIGGER_SF_2017,
    "2018" : SUBLEAD_TRIGGER_SF_2018, 
    "2022preEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2022postEE" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023preBPix" : SUBLEAD_TRIGGER_SF_2018, #FIXME
    "2023postBPix" : SUBLEAD_TRIGGER_SF_2018 #FIXME
}

def trigger_sf(events, central_only, year):
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9"), ("Photon", "pt")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message)

    variations_lead = systematic_from_bins(
        bins = lead_trigger_sf_bins[year], 
        variables = {
            "photon_r9" : events.LeadPhoton.r9,
            "photon_eta" : abs(events.LeadPhoton.eta),
            "photon_pt" : events.LeadPhoton.pt
        },
        central_only = central_only
    )

    variations_sublead = systematic_from_bins(
        bins = sublead_trigger_sf_bins[year],
        variables = {
            "photon_r9" : events.SubleadPhoton.r9,
            "photon_eta" : abs(events.SubleadPhoton.eta),
            "photon_pt" : events.SubleadPhoton.pt
        },
        central_only = central_only
    )
 
    # Multiply up/down/central variations together, following this treatment in flashgg:
    # https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Systematics/interface/DiPhotonFromSinglePhotonViewBase.h#L85-L87
    variations = {}
    for key in variations_lead.keys():
        variations[key] = variations_lead[key] * variations_sublead[key]

    return variations

############
### FNUF ###
############

from higgs_dna.systematics.data.fnuf import FNUF_2016, FNUF_2017, FNUF_2018
fnuf_bins = {
    "2016" : FNUF_2016,
    "2016preVFP" : FNUF_2016,
    "2016postVFP" : FNUF_2016,
    "2017" : FNUF_2017,
    "2018" : FNUF_2018,
    "2022preEE" : FNUF_2018,
    "2022postEE" : FNUF_2018,
    "2023preBPix" : FNUF_2018,
    "2023postBPix" : FNUF_2018,
    "2024" : FNUF_2018,
}

def fnuf_unc(events, year, nominal_only, modify_nominal):
    """

    """
    required_fields = [
        ("Photon", "eta"), ("Photon", "r9")
    ]

    missing_fields = awkward_utils.missing_fields(events, required_fields)

    if missing_fields:
        message = "[photon_systematics : photon_preselection_sf] The events array is missing the following fields: %s which are needed as inputs." % (str(missing_fields))
        logger.exception(message)
        raise ValueError(message) 

    photons = events.Photon

    loc = "all"  # fixed default; previously a function argument

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "eb":
        mask = photons.isScEtaEB == True
    elif loc == "ee":
        mask = photons.isScEtaEE == True

    variations = ic_systematic_from_bins(
        bins = fnuf_bins[year], 
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations

################
### Material ###
################

from higgs_dna.systematics.data.material import MATERIAL_2016, MATERIAL_2017, MATERIAL_2018
material_bins = {
    "2016" : MATERIAL_2016,
    "2016preVFP" : MATERIAL_2016,
    "2016postVFP" : MATERIAL_2016,
    "2017" : MATERIAL_2017,
    "2018" : MATERIAL_2018,
    "2022preEE" : MATERIAL_2018,
    "2022postEE" : MATERIAL_2018,
    "2023preBPix" : MATERIAL_2018,
    "2023postBPix" : MATERIAL_2018,
    "2024" : MATERIAL_2018
}

def material_unc(events, year, nominal_only, modify_nominal):
    photons = events.Photon

    loc = "all"  # fixed default; previously a function argument

    if loc == "all":
        mask = photons.pt > 0
    elif loc == "central_barrel":
        mask = abs(photons.eta) <= 1.0
    elif loc == "outer_barrel":
        mask = (abs(photons.eta) > 1.0) & (abs(photons.eta) <= 1.5)
    elif loc == "forward":
        mask = abs(photons.eta) > 1.5
    
    variations = ic_systematic_from_bins(
        bins = material_bins[year],
        variables = {
            "photon_eta" : abs(events.Photon.eta),
            "photon_r9" : events.Photon.r9
        },
        branch = photons.pt,
        nominal_only = nominal_only,
        modify_nominal = modify_nominal,
        mask = mask
    )

    return variations  


def dummy_photon_pt_syst(events):
    photons = events.Photon

    variations = {}
    variations["up"] = photons.pt + awkward.ones_like(photons.pt)
    variations["down"] = photons.pt - awkward.ones_like(photons.pt)
    return variations


#######################
### Photon MC Smear ###
#######################

def photon_mc_smear(events, r9, loc):
    photons = events.Photon

    # Split into high/low r9 and EE/EB
    if r9 == "high":
        r9_cut = photons.r9 > 0.94
    elif r9 == "low":
        r9_cut = photons.r9 <= 0.94

    if loc == "eb":
        loc_cut = photons.isScEtaEB == True
    elif loc == "ee":
        loc_cut = photons.isScEtaEE == True

    mask = r9_cut & loc_cut

    variations = {}
    variations["up"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaUp,
            photons.pt
    )
    variations["down"] = awkward.where(
            mask,
            photons.pt + photons.dEsigmaDown,
            photons.pt
    )

    return variations

photon_scale_FILE = {
    "2022preEE" : "jsonpog-integration/POG/EGM/2022_Summer22/photonSS_EtDependent.json",
    "2022postEE" : "jsonpog-integration/POG/EGM/2022_Summer22EE/photonSS_EtDependent.json",
    "2023preBPix" : "jsonpog-integration/POG/EGM/2023_Summer23/photonSS_EtDependent.json",
    "2023postBPix" : "jsonpog-integration/POG/EGM/2023_Summer23BPix/photonSS_EtDependent.json",
    "2024" : "jsonpog-integration/POG/EGM/2024_Summer24/photonSS_EtDependent.json"
}

smear_names = {
    "2022preEE" : "EGMSmearAndSyst_PhoPTsplit_2022preEE",
    "2022postEE" : "EGMSmearAndSyst_PhoPTsplit_2022postEE",
    "2023preBPix" : "EGMSmearAndSyst_PhoPTsplit_2023preBPIX",
    "2023postBPix" : "EGMSmearAndSyst_PhoPTsplit_2023postBPIX",
    "2024" : "EGMSmearAndSyst_PhoPT_2024"
}

scale_names = {
    "2022preEE" : "EGMScale_Compound_Pho_2022preEE",
    "2022postEE" : "EGMScale_Compound_Pho_2022postEE",
    "2023preBPix" : "EGMScale_Compound_Pho_2023preBPIX",
    "2023postBPix" : "EGMScale_Compound_Pho_2023postBPIX",
    "2024" : "Scale",  # fix: compound_corrections.name in 2024 JSON is "Scale"
}

def photon_scale_smear_run3(events, year, is_data):
    """
    Align implementation with hzg_photon_systematics.py:
      - data: apply scale -> corrected_pt, corrected_energyErr then return
      - MC: apply smear central -> corrected_pt, corrected_energyErr
            + smear_up/down -> dEsigmaUp/Down and energyErr_dEsigmaUp/Down
            + scale_up/down -> pt_ScaleUp/Down and energyErr_ScaleUp/Down
    """
    logger.info("[Photon Systematics] Applying photon scale and smear corrections for year %s", year)

    required_fields = [
        ("Photon", "eta"), ("Photon", "pt"), ("Photon", "r9"), ("Photon", "energyErr")
    ]
    if is_data:
        required_fields.append(("run"))
        required_fields.append(("Photon", "seedGain"))

    missing_fields = awkward_utils.missing_fields(events, required_fields)
    if missing_fields:
        logger.warning(
            "[Photon Systematics] Photon scale and smear corrections will not be applied, because the following fields are missing: %s",
            str(missing_fields),
        )
        return events

    # Use the same evaluator choice pattern as hzg_photon_systematics.py
    if is_data:
        evaluator = correctionlib.CorrectionSet.from_file(misc_utils.expand_path(photon_scale_FILE[year]))
    else:
        evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(photon_scale_FILE[year]))

    photons = events["Photon"]
    n_photons = awkward.num(photons)
    photons_flattened = awkward.flatten(photons)

    photons_scEta = awkward.to_numpy(photons_flattened.eta)
    photons_AbsScEta = numpy.abs(photons_scEta)
    photons_pt = awkward.to_numpy(photons_flattened.pt)
    photons_r9 = awkward.to_numpy(photons_flattened.r9)
    photons_energyerr = awkward.to_numpy(photons_flattened.energyErr)
    photons_smear_eta = photons_scEta if year == "2024" else photons_AbsScEta
    photon_pt_min = 10.0

    # NEW: log first 10 photon pt before applying corrections
    _nprint = int(min(10, len(photons_pt)))
    logger.debug(
        "[Photon Systematics] Photon pt BEFORE scale+smear (first %d): %s",
        _nprint,
        photons_pt[:_nprint],
    )

    if is_data:
        photons_seedGain = awkward.to_numpy(photons_flattened.seedGain)
        run_arr_flattened = numpy.repeat(awkward.to_numpy(events.run), n_photons)

        if year == "2024":
            scale = awkward.where(
                (photons_AbsScEta > 3.0) | (photons_pt < photon_pt_min),
                awkward.ones_like(photons_pt, dtype=float),
                evaluator.compound[scale_names[year]].evaluate(
                    "scale",
                    run_arr_flattened,
                    photons_scEta,
                    photons_r9,
                    photons_pt,
                    photons_seedGain,
                ),
            )
        else:
            scale = awkward.where(
                (photons_AbsScEta > 3.0) | (photons_pt < photon_pt_min),
                awkward.ones_like(photons_pt, dtype=float),
                evaluator.compound[scale_names[year]].evaluate(
                    "scale",
                    run_arr_flattened,
                    photons_scEta,
                    photons_r9,
                    photons_AbsScEta,
                    photons_pt,
                    photons_seedGain,
                ),
            )
                        
        corrected_pt = awkward.to_numpy(photons_pt * scale)
        events["Photon", "corrected_pt"] = awkward.unflatten(corrected_pt, n_photons)

        # NEW: log first 10 photon pt after applying corrections (data: scale only)
        logger.debug(
            "[Photon Systematics] Photon pt AFTER scale (data) (first %d): %s",
            _nprint,
            corrected_pt[:_nprint],
        )

        # NOTE: match hzg: smear is evaluated on corrected_pt and energyErr updated accordingly, scaled by scale
        evaluator_smear = _core.CorrectionSet.from_file(misc_utils.expand_path(photon_scale_FILE[year]))
        smear = awkward.where(
            (photons_AbsScEta > 3.0) | (photons_pt < photon_pt_min),
            awkward.zeros_like(photons_pt, dtype=float),
            evaluator_smear[smear_names[year]].evalv("smear", corrected_pt, photons_r9, photons_smear_eta),
        )
        corrected_energyErr = numpy.sqrt((photons_energyerr) ** 2 + (photons_pt * numpy.cosh(photons_scEta) * smear) ** 2) * scale
        events["Photon", "corrected_energyErr"] = awkward.unflatten(corrected_energyErr, n_photons)
        events["Photon", "pt"] = awkward.unflatten(corrected_pt, n_photons)
        events["Photon", "energyErr"] = awkward.unflatten(corrected_energyErr, n_photons)
        return events

    # --- MC branch: apply central smear ---
    smear = awkward.where(
        (photons_AbsScEta > 3.0) | (photons_pt < photon_pt_min),
        awkward.zeros_like(photons_pt, dtype=float),
        evaluator[smear_names[year]].evalv("smear", photons_pt, photons_r9, photons_smear_eta),
    )
    rng = numpy.random.default_rng(seed=123)
    smear_random = rng.normal(loc=0.0, scale=1.0, size=len(photons_pt))
    smear_val = awkward.where(
        (photons_AbsScEta > 3.0) | (photons_pt < photon_pt_min),
        awkward.ones_like(photons_pt, dtype=float),
        1.0 + numpy.abs(smear) * smear_random,
    )

    corrected_pt = photons_pt * smear_val
    corrected_energyErr = numpy.sqrt((photons_energyerr) ** 2 + (photons_pt * numpy.cosh(photons_scEta) * smear) ** 2) * smear_val

    # NEW: log first 10 photon pt after applying corrections (MC: central smear)
    logger.debug(
        "[Photon Systematics] Photon pt AFTER smear (MC) (first %d): %s",
        _nprint,
        awkward.to_numpy(corrected_pt)[:_nprint],
    )

    # --- MC branch: smear systematics ---
    for syst in ["smear_up", "smear_down"]:
        smear_syst = awkward.where(
            (photons_AbsScEta > 3.0) | (photons_pt < photon_pt_min),
            awkward.zeros_like(photons_pt, dtype=float),
            evaluator[smear_names[year]].evalv(syst, photons_pt, photons_r9, photons_smear_eta),
        )

        smear_syst_val = awkward.where(
            (photons_AbsScEta > 3.0) | (photons_pt < photon_pt_min),
            awkward.ones_like(photons_pt, dtype=float),
            1.0 + numpy.abs(smear_syst) * smear_random,
        )
        pt_syst = photons_pt * smear_syst_val
        events["Photon", "dEsigma" + syst.replace("smear_", "").capitalize()] = awkward.unflatten(
            pt_syst - corrected_pt,
            n_photons,
        )
        events["Photon", "energyErr_dEsigma" + syst.replace("smear_", "").capitalize()] = awkward.unflatten(
            numpy.sqrt((photons_energyerr) ** 2 + (photons_pt * numpy.cosh(photons_scEta) * smear_syst) ** 2) * smear_syst_val,
            n_photons,
        )

    # write central
    events["Photon", "corrected_pt"] = awkward.unflatten(corrected_pt, n_photons)
    events["Photon", "corrected_energyErr"] = awkward.unflatten(corrected_energyErr, n_photons)
    photons = events["Photon"]

    # --- MC branch: scale systematics (match hzg: evaluated from the same correction name as smear) ---
    for syst in ["scale_up", "scale_down"]:
        scale = awkward.where(
            (abs(photons.eta) > 3.0) | (photons.pt < photon_pt_min),
            awkward.ones_like(photons.pt, dtype=float),
            awkward.unflatten(
                evaluator[smear_names[year]].evalv(syst, photons_pt, photons_r9, photons_smear_eta),
                n_photons,
            ),
        )
        events["Photon", "pt_Scale" + syst.replace("scale_", "").capitalize()] = photons.corrected_pt * scale
        events["Photon", "energyErr_Scale" + syst.replace("scale_", "").capitalize()] = photons.corrected_energyErr * scale

    events["Photon", "pt"] = awkward.unflatten(corrected_pt, n_photons)
    events["Photon", "energyErr"] = awkward.unflatten(corrected_energyErr, n_photons)

    # NEW: log systematics branches (match photon_systematics_copy.py style; first 10)
    if ("Photon" in events.fields) and ("dEsigmaUp" in events["Photon"].fields):
        logger.debug(
            "[Photon Systematics] Photon pt Smear Up (dEsigmaUp) (first %d): %s",
            _nprint,
            awkward.to_numpy(awkward.flatten(events["Photon", "dEsigmaUp"]))[:_nprint],
        )
    if ("Photon" in events.fields) and ("dEsigmaDown" in events["Photon"].fields):
        logger.debug(
            "[Photon Systematics] Photon pt Smear Down (dEsigmaDown) (first %d): %s",
            _nprint,
            awkward.to_numpy(awkward.flatten(events["Photon", "dEsigmaDown"]))[:_nprint],
        )
    if ("Photon" in events.fields) and ("pt_ScaleUp" in events["Photon"].fields):
        logger.debug(
            "[Photon Systematics] Photon pt Scale Up (pt_ScaleUp) (first %d): %s",
            _nprint,
            awkward.to_numpy(awkward.flatten(events["Photon", "pt_ScaleUp"]))[:_nprint],
        )
    if ("Photon" in events.fields) and ("pt_ScaleDown" in events["Photon"].fields):
        logger.debug(
            "[Photon Systematics] Photon pt Scale Down (pt_ScaleDown) (first %d): %s",
            _nprint,
            awkward.to_numpy(awkward.flatten(events["Photon", "pt_ScaleDown"]))[:_nprint],
        )

    return events
