import awkward

import numpy

from correctionlib import _core

from higgs_dna.utils import misc_utils, awkward_utils

###################################
### b-tag continuous reshape SF ###
###################################

BTAG_RESHAPE_SF_FILE = {
    "2016" : "higgs_dna/systematics/data/2016postVFP_UL/btag_mceff.json", 
    "2016preVFP" : "higgs_dna/systematics/data/2016preVFP_UL/btagging.json", 
    "2016postVFP" : "higgs_dna/systematics/data/2016postVFP_UL/btagging.json", 
    "2017" : "higgs_dna/systematics/data/2017_UL/btagging.json",
    "2018" : "higgs_dna/systematics/data/2018_UL/btagging.json"
}

DEEPJET_RESHAPE_SF = {
    "2016" : "deepJet_shape", 
    "2016preVFP" : "deepJet_shape",
    "2016postVFP" : "deepJet_shape",
    "2017" : "deepJet_shape",
    "2018" : "deepJet_shape"
}

DEEPJET_VARIATIONS = { # b, c, light
    "up_correlated" : [5, 4, 0], 
    "down_correlated" : [5, 4, 0],
    "up_uncorrelated" : [5, 4, 0],
    "down_uncorrelated" : [5, 4, 0],
    # "up_jes" : [5, 0], # applicable to b (5) and light (0) jets, but not charm (4)
    # "up_lf" : [5],
    # "up_hfstats1" : [5],
    # "up_hfstats2" : [5],
    # "up_cferr1" : [4],
    # "up_cferr2" : [4],
    # "up_hf" : [0],
    # "up_lfstats1" : [0],
    # "up_lfstats2" : [0],
    # "down_jes" : [5, 0], # applicable to b (5) and light (0) jets, but not charm(4)
    # "down_lf" : [5],
    # "down_hfstats1" : [5],
    # "down_hfstats2" : [5],
    # "down_cferr1" : [4],
    # "down_cferr2" : [4],
    # "down_hf" : [0],
    # "down_lfstats1" : [0],
    # "down_lfstats2" : [0],
}


def btag_deepjet_wp_sf(events, year, central_only, input_collection):
    """
    See:
        - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/BTV_bjets_Run2_UL/
        - https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/btvExample.py

    Note: application of SFs should not change the overall normalization of a sample (before any b-tagging selection) and each sample should be adjusted by an overall weight derived in a phase space with no requirements on b-jets such that the normalization is unchanged. TODO: link BTV TWiki that describes this.
    """
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "hadronFlavour"), (input_collection, "btagDeepFlavB") 
    ]
    missing_fields = awkward_utils.missing_fields(events, required_fields)

    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_RESHAPE_SF_FILE[year]))
   
    jets = events[input_collection]
    jets["flavor"] = jets.hadronFlavour
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end 
    jets_flattened = awkward.flatten(jets)

    jet_flavor = awkward.to_numpy(jets_flattened.flavor)
    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.0,
        2.49999 # SFs only valid up to eta 2.5
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        99999999.
    )

    variations_list = ["central"]
    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()

    variations = {}

    central_sf = numpy.ones_like(jet_flavor)
    for f in [0, 4, 5]:
        central_sf = numpy.where(
            jet_flavor == f,
            evaluator["deepJet_comb" if f > 0 else "deepJet_incl"].evalv(
                "central",
                "M",
                numpy.ones_like(jet_flavor) * f,
                jet_abs_eta,
                jet_pt
            ),
            central_sf
        )

    variations["central"] = awkward.unflatten(central_sf, n_jets)

    for var in variations_list:
        if var == "central":
            continue
        applicable_flavors = DEEPJET_VARIATIONS[var] # the up/down variations are only applicable to specific flavors of jet
        var_sf = central_sf 
        for f in applicable_flavors:
            var_sf = numpy.where(
                jet_flavor == f,
                evaluator["deepJet_comb" if f > 0 else "deepJet_incl"].evalv(
                    var,
                    "M",
                    numpy.ones_like(jet_flavor) * f,
                    jet_abs_eta,
                    jet_pt
                ),
                var_sf
            )

        variations[var] = awkward.unflatten(var_sf, n_jets) # make jagged again

    for var in variations.keys():
        # Set SFs = 1 for jets which are not applicable (pt <= 20 or |eta| >= 2.5)
        variations[var] = awkward.where(
                (jets.pt <= 20.0) | (abs(jets.eta) >= 2.5),
                awkward.ones_like(variations[var]),
                variations[var]
        )
    # Swap 'up' or 'down' with the rest of the variation name
    swapped_variations = {
        (var.replace("up_", "") + "_up") if "up" in var else (var.replace("down_", "") + "_down") if "down" in var else var: val
        for var, val in variations.items()
    }

    return swapped_variations
def btag_deepjet_mujet_sf(events, year, central_only, input_collection, working_point ="M"):
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "hadronFlavour")
    ]    
    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_RESHAPE_SF_FILE[year]))
    
    jets = events[input_collection]   
    jets["flavor"] = jets.hadronFlavour
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end 
    jets_flattened = awkward.flatten(jets)

    jet_flavor = awkward.to_numpy(jets_flattened.flavor)
    is_light = jet_flavor == 0
    
    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.0,
        2.49999 # SFs only valid up to eta 2.5
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        999.99
    )
    variations_list = ["central"]

    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()

    variations = {}
    print(type(jet_flavor))
    print(type(numpy.where(is_light,4,jet_flavor)))
    central_sf = evaluator["deepJet_mujets"].evalv(
            "central",
            working_point,
            numpy.where(is_light,4,jet_flavor),
            jet_abs_eta,
            jet_pt
    )    
    central_sf = numpy.where(is_light,1,central_sf)
    variations["central"] = awkward.unflatten(central_sf, n_jets)
    if not central_only:
        syst_vars = ["up_correlated", "down_correlated", "up_uncorrelated", "down_uncorrelated"]
        for syst_var in syst_vars:
            syst = evaluator["deepJet_mujets"].evalv(
                    syst_var,
                    working_point,
                    numpy.where(is_light,4,jet_flavor),
                    jet_abs_eta,
                    jet_pt,
            )
            if "up_correlated" in syst_var:
                syst_var_name = "up_correlated"
            elif "down_correlated" in syst_var:
                syst_var_name = "down_correlated"
            elif "up_uncorrelated" in syst_var:
                syst_var_name = "up_uncorrelated"
            elif "down_uncorrelated" in syst_var:
                syst_var_name = "down_uncorrelated"
            variations[syst_var_name] = awkward.unflatten(syst, n_jets)
    return variations

def btag_deepjet_incl_sf(events, year, central_only, input_collection, working_point ="M"):
    required_fields = [
        (input_collection, "eta"), (input_collection, "pt"), (input_collection, "hadronFlavour")
    ]    
    evaluator = _core.CorrectionSet.from_file(misc_utils.expand_path(BTAG_RESHAPE_SF_FILE[year]))
    
    jets = events[input_collection]   
    jets["flavor"] = jets.hadronFlavour
    n_jets = awkward.num(jets) # save n_jets to convert back to jagged format at the end 
    jets_flattened = awkward.flatten(jets)

    jet_flavor = awkward.to_numpy(jets_flattened.flavor)
    is_light = jet_flavor == 0

    jet_abs_eta = numpy.clip(
        awkward.to_numpy(abs(jets_flattened.eta)),
        0.0,
        2.49999 # SFs only valid up to eta 2.5
    )
    jet_pt = numpy.clip(
        awkward.to_numpy(jets_flattened.pt),
        20.0, # SFs only valid for pT > 20.
        999.99
    )
    variations_list = ["central"]

    if not central_only:
        variations_list += DEEPJET_VARIATIONS.keys()
    
    variations = {}
    central_sf = evaluator["deepJet_incl"].evalv(
            "central",
            working_point,
            numpy.where(~is_light,0,jet_flavor),
            jet_abs_eta,
            jet_pt
    )   
    central_sf = numpy.where(~is_light,1,central_sf)
    variations["central"] = awkward.unflatten(central_sf, n_jets)
    if not central_only:
        syst_vars = ["up_correlated", "down_correlated", "up_uncorrelated", "down_uncorrelated"]
        for syst_var in syst_vars:
            syst = evaluator["deepJet_incl"].evalv(
                    syst_var,
                    working_point,
                    numpy.where(is_light,4,jet_flavor),
                    jet_abs_eta,
                    jet_pt,
            )
            if "up_correlated" in syst_var:
                syst_var_name = "up_correlated"
            elif "down_correlated" in syst_var:
                syst_var_name = "down_correlated"
            elif "up_uncorrelated" in syst_var:
                syst_var_name = "up_uncorrelated"
            elif "down_uncorrelated" in syst_var:
                syst_var_name = "down_uncorrelated"
            variations[syst_var_name] = awkward.unflatten(syst, n_jets)
    return variations
