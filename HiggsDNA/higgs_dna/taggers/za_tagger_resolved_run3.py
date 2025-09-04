import awkward as ak
import time
import numpy
import numba
import vector

from pdb import set_trace

vector.register_awkward()

import logging
# logger = logging.getLogger(__name__)
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.taggers.tagger import Tagger, NOMINAL_TAG 
from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.selections import object_selections, lepton_selections, jet_selections, tau_selections, physics_utils
from higgs_dna.selections import gen_selections

# Run3 Photon cut-based ID
# H/E Loose
        # "hoe_barrel": 0.129991,
        # "hoe_endcap": 0.153428,
# H/E Medium
        # "hoe_barrel": 0.0583054,
        # "hoe_endcap": 0.00518075,
# H/E Tight
        # "hoe_barrel": 0.0417588,
        # "hoe_endcap": 0.00254267,
#-------------------------------------
# Charged Hadron Iso Loose
        # "PFChIso_barrel": 1.88518,
        # "PFChIso_endcap": 1.65396,
# Charged Hadron Iso Medium
        # "PFChIso_barrel": 0.939289,
        # "PFChIso_endcap": 0.970286,
# Charged Hadron Iso Tight
        # "PFChIso_barrel": 0.316306,
        # "PFChIso_endcap": 0.292664,
#-------------------------------------
# HCal Iso Loose
        # "PFHCalIso_barrel": [6.34397,   0.0100547, 5.78332e-05],
        # "PFHCalIso_endcap": [1.85881,   0.0116989, 7.47603-05],
# HCal Iso Medium
        # "PFHCalIso_barrel": [2.18903,   0.0100547, 5.78332e-05],
        # "PFHCalIso_endcap": [0.0336699, 0.0116989, 7.47603e-05],
# HCal Iso Tight
        # "PFHCalIso_barrel": [0.39057,   0.0100547, 5.78332e-05],
        # "PFHCalIso_endcap": [0.0292617, 0.0116989, 7.47603e-05],
#-------------------------------------
# Discard 
# sieie Loose
        # "sieie_barrel": 0.0114521,
        # "sieie_endcap": 0.0276744,
# sieie Medium
        # "sieie_barrel": 0.0100086,
        # "sieie_endcap": 0.0268736,
# sieie Tight
        # "sieie_barrel": 0.00999299,
        # "sieie_endcap": 0.0268702,
#-------------------------------------
# ECal Iso Loose
        # "PFECalIso_barrel": [0.703789, 0.000652035],
        # "PFECalIso_endcap": [6.61585, 0.000195486],
# ECal Iso Medium
        # "PFECalIso_barrel": [0.227697, 0.000652035],
        # "PFECalIso_endcap": [1.124, 0.000195486],
# ECal Iso Tight
        # "PFECalIso_barrel": [0.14189, 0.000652035],
        # "PFECalIso_endcap": [1.04269, 0.000195486],

DUMMY_VALUE = -999.
DEFAULT_OPTIONS = {
    "photons" : {
        "use_central_nano" : True,
        "pt" : 10.0,
        "eta" : [
            [0.0, 1.4442],
            [1.566, 2.5]
        ],
        "mvaID_barrel" : -0.4,
        "mvaID_endcap" : -0.58,
        # Loose
        # "hoe_barrel": 0.129991,
        # "hoe_endcap": 0.153428,
        # "PFChIso_barrel": 1.88518,
        # "PFChIso_endcap": 1.65396,
        # "PFHCalIso_barrel": [6.34397,   0.0100547, 5.78332e-05],
        # "PFHCalIso_endcap": [1.85881,   0.0116989, 7.47603-05],
        # "sieie_barrel": 0.0114521,
        # "sieie_endcap": 0.0276744,
        # "PFECalIso_barrel": [0.703789, 0.000652035],
        # "PFECalIso_endcap": [6.61585, 0.000195486],
        # # Medium
        # "hoe_barrel": 0.0583054,
        # "hoe_endcap": 0.00518075,
        # "PFChIso_barrel": 0.939289,
        # "PFChIso_endcap": 0.970286,
        # "PFHCalIso_barrel": [2.18903,   0.0100547, 5.78332e-05],
        # "PFHCalIso_endcap": [0.0336699, 0.0116989, 7.47603e-05],
        # "sieie_barrel": 0.0100086,
        # "sieie_endcap": 0.0268736,
        # "PFECalIso_barrel": [0.227697, 0.000652035],
        # "PFECalIso_endcap": [1.124, 0.000195486],
        # # Tight
        "hoe_barrel": 0.0417588,
        "hoe_endcap": 0.00254267,
        "PFChIso_barrel": 0.316306,
        "PFChIso_endcap": 0.292664,
        "PFHCalIso_barrel": [0.39057,   0.0100547, 5.78332e-05],
        "PFHCalIso_endcap": [0.0292617, 0.0116989, 7.47603e-05],
        "sieie_barrel": 0.00999299,
        "sieie_endcap": 0.0268702,
        "PFECalIso_barrel": [0.14189, 0.000652035],
        "PFECalIso_endcap": [1.04269, 0.000195486],
        #------------------------------------------
        "PFHCalIso_EA_EB_1": [0.17005,  -0.000835],
        "PFHCalIso_EA_EB_2": [0.208571, -0.000905],
        "PFHCalIso_EA_EE_1": [0.246494, -0.000722],
        "PFHCalIso_EA_EE_2": [0.306529, -0.000608],
        "PFHCalIso_EA_EE_3": [0.322673, -0.000750],
        "PFHCalIso_EA_EE_4": [0.315793, -0.000795],
        "PFHCalIso_EA_EE_5": [0.36531,  -0.000439],
        "e_veto" : 0.5,
        "PFECalIso_EA_EB_1": [0.0866519, -0.0002296],
        "PFECalIso_EA_EB_2": [0.0730397, -0.0002134],
        "PFECalIso_EA_EE_1": [0.0542479, -0.00010934],
        "PFECalIso_EA_EE_2": [0.0486181, -6.2098e-05],
        "PFECalIso_EA_EE_3": [0.0412923, -9.63732e-06],
        "PFECalIso_EA_EE_4": [0.03555,   -5.79549e-05],
        "PFECalIso_EA_EE_5": [0.0360895, -7.28546e-06]
    },
    "zgammas" : {
        "relative_pt_gamma" : 15.0/110.,
        "mass_h" : [95., 180.],
        "mass_sum" : 185,
        "select_highest_pt_sum" : True
    },
    "single_muon_trigger":{
        "2016":["HLT_IsoMu24", "HLT_IsoTkMu24"],
        "2017":["HLT_IsoMu27"],
        "2018":["HLT_IsoMu24"],
        "2022":["HLT_IsoMu24"],
        "2023":["HLT_IsoMu24"]
    },
    "single_ele_trigger":{
        "2016":["HLT_Ele27_WPTight_Gsf"],
        "2017":["HLT_Ele32_WPTight_Gsf_L1DoubleEG"],
        "2018":["HLT_Ele32_WPTight_Gsf"],
        "2022":["HLT_Ele30_WPTight_Gsf"],
        "2023":["HLT_Ele30_WPTight_Gsf"]
    },
    "double_muon_trigger":{
        "2016":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"],
        "2017":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
        "2018":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
        "2022":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
        "2023":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"]
    },
    "double_ele_trigger":{
        "2016":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
        "2017":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2018":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2022":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2023":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]
    },
    # "trigger" : {
    #     "2016" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele27_WPTight_Gsf"],
    #     "2016UL_preVFP" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele27_WPTight_Gsf"],
    #     "2016UL_postVFP" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele27_WPTight_Gsf"],
    #     "2017" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Ele27_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG", "HLT_Ele35_WPTight_Gsf", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_IsoMu24_eta2p1"],
    #     "2018" : ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Mu37_TkMu27", "HLT_IsoMu20", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55", "HLT_IsoMu24_eta2p1", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_DoubleEle25_CaloIdL_MW", "HLT_Ele27_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG", "HLT_Ele35_WPTight_Gsf", "HLT_Ele20_WPLoose_Gsf"],
    #     "2022" : ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Ele30_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf", "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"]
    # }, 
    "electrons" : {
        "pt" : 7.0
    },
    "muons" : {
        "pt" : 5.0
    },
    "lead_ele_pt":{
        "2016": 30,
        "2017": 35,
        "2018": 35,
        "2022": 35,
        "2023": 35
    },
    "lead_mu_pt":{
        "2016": 25,
        "2017": 28,
        "2018": 25,
        "2022": 25,
        "2023": 25
    },
    "jets" : {
        "pt" : 30.0,
        "eta" : 4.7,
        "eta" : 4.7,
        "dr_photons" : 0.4,
        "dr_electrons" : 0.4,
        "dr_muons" : 0.4,
        "jets_horn" : {
            "pt" : 50.0,
            "eta" : [2.5, 3.0]
        }
    },
    "btag_med": {
        "2016preVFP": 0.2598,
        "2016postVFP": 0.2489,
        "2017": 0.3040,
        "2018": 0.2783,
        "2022preEE": 0.3086,
        "2022postEE": 0.3196,
        "2023preBPix": 0.2431,
        "2023postBPix": 0.2435
    },
    "FSRphotons" : {
        "iso" : 1.8,
        "eta" : 2.4,
        "pt" : 2.0,
        "dROverEt2" : 0.012
    },
    "gen_info" : {
        "calculate" : False,
        "max_dr" : 0.2,
        "max_pt_diff" : 15.
    }  
}


# Diphoton preselection below synced with flashgg, see details in:
#   - https://indico.cern.ch/event/1071721/contributions/4551056/attachments/2320292/3950844/HiggsDNA_DiphotonPreselectionAndSystematics_30Sep2021.pdf
# EA from 
#   - https://indico.cern.ch/event/1204277/contributions/5064356/attachments/2538496/4369369/CutBasedPhotonID_20221031.pdf


def to_momentum4d(obj):
    out = ak.zip({
        "pt": obj.pt,
        "eta": obj.eta,
        "phi": obj.phi,
        "mass": obj.mass if "mass" in obj.fields else ak.zeros_like(obj.pt),
    }, with_name="Momentum4D")
    # 强制 materialize layout 以触发行为绑定
    out = ak.Array(out.layout)  # 👈 核心所在
    return out

class ZaTaggerRun3(Tagger):
    def __init__(self, name = "default_zgamma_tagger", options = {}, is_data = None, year = None):
        super(ZaTaggerRun3, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )


    def calculate_selection(self, events):
        """
        Select photons and create diphoton pairs.
        Add a record "Diphoton" to events array with relevant information about each diphoton pair.
        In principle, there can be more than one Diphoton pair per event.
        """

        # Determine what type of rho variable is available in nanoAOD
        # To be deleted once a standard rho variable is added to central nanoAOD
        if self.options["photons"]["use_central_nano"]:
            if "fixedGridRhoAll" in events.fields:
                rho = events.fixedGridRhoAll
            elif "fixedGridRhoFastjetAll" in events.fields:
                rho = events.fixedGridRhoFastjetAll
            elif "Rho_fixedGridRhoAll" in events.fields:
                rho = events.Rho_fixedGridRhoAll
            else:
                raise RuntimeError("Rho not found in central nanoAOD! Cannot apply PU correction.")
        else:
            rho = ak.ones_like(events.Photon)

        if not self.is_data:
            self.overlap_removal(events=events)

        zgamma_selection, zgammas = self.produce_and_select_zgammas(
                events = events,
                rho = rho,
                options = self.options["zgammas"]
        )

        if not self.is_data and self.options["gen_info"]["calculate"]:
            zgammas = self.calculate_gen_info(zgammas, self.options["gen_info"])

        return zgamma_selection, zgammas 

    def overlap_removal(self, events):
        """
        Select isolation photons in events
        Add number of isolation photons (n_iso_photons) in output .parquet file to indetify thé overlap events
        """
        
        """
        statusFlags usage: (events.GenPart.statusFlags // numpy.power(2, i)) % 2 == 1
        "statusFlags" is a number with 14 bits. 
        Filling "1" on corresponding digit when the particle meets one of the 14 conditions, else remaining "0".
        Echo paticles can meet more than one kind of condition, thus, more than one digit in "statusFlags" is "1".
        """
        iso_photons_cut = (events.GenPart.pdgId == 22) & (events.GenPart.pt > 15) & (abs(events.GenPart.eta) < 2.6) & (( (events.GenPart.statusFlags // numpy.power(2, 0)) % 2 == 1 ) | ( (events.GenPart.statusFlags // numpy.power(2, 8)) % 2 == 1 ))
        iso_photons = events.GenPart[iso_photons_cut]

        truth_objects_cut =  (events.GenPart.pdgId != 22) & (events.GenPart.pt > 5) & ( (events.GenPart.statusFlags // numpy.power(2, 8)) % 2 == 1 ) 
        truth_objects = events.GenPart[truth_objects_cut]

        iso_cut = object_selections.delta_R(iso_photons, truth_objects, 0.05)
        iso_photons = iso_photons[iso_cut]

        n_iso_photons = ak.num(iso_photons)
        awkward_utils.add_field(events, "n_iso_photons", n_iso_photons, overwrite=True)

    def produce_and_select_zgammas(self, events, rho, options):
        """
        Perform diphoton preselection.
        For events with more than 2 photons, more than 1 diphoton candidate
        per event is possible.

        :param events: events array to calculate diphoton candidates from
        :type events: ak.highlevel.Array
        :param photons: array of selected photons from events
        :type photons: ak.highlevel.Array
        :param options: dictionary containing configurable options for the diphoton preselection
        :type options: dict
        :return: boolean array indicating which events pass the diphoton preselection
        :rtype: ak.highlevel.Array
        """

        start = time.time()

        # events = events[(events.run == 356077) & (events.luminosityBlock == 158) & (events.event == 196518696)]
        # Electrons
        electron_cut = lepton_selections.select_electrons(
            electrons = events.Electron,
            options = self.options["electrons"],
            clean = {
            },
            name = "SelectedElectron",
            tagger = self,
            year = self.year[:4]
        )
        
        electrons = awkward_utils.add_field(
            events = events,
            name = "SelectedElectron",
            data = events.Electron[electron_cut]
        )
        
        # generate the index in the original array and add to electrons
        arr = ak.local_index(events.Electron["pt"], axis=1)[electron_cut]
        electron_idx = ak.mask(arr, ak.num(arr) > 0)
        awkward_utils.add_field(events = electrons, name = "Idx", data = electron_idx)

        # Muons
        muon_cut = lepton_selections.select_muons(
            muons = events.Muon,
            options = self.options["muons"],
            clean = {
            },
            name = "SelectedMuon",
            tagger = self
        )

        muons = awkward_utils.add_field(
            events = events,
            name = "SelectedMuon",
            data = events.Muon[muon_cut]
        )

        if not self.is_data:
            # gen mu not reco
            gen_muons = events.GenPart[(abs(events.GenPart.pdgId) == 13)]
            reco_muons = events.Muon

            unmatched_gen_mask = ak.fill_none(object_selections.delta_R(gen_muons, reco_muons, 0.4), True)
            
            unmatched_gen_muons = gen_muons[unmatched_gen_mask]
            unmatched_gen_muons = unmatched_gen_muons[unmatched_gen_muons.pt > self.options["muons"]["pt"]]

            unmatched_gen_muons = awkward_utils.add_field(
                events = events,
                name = "SelectedGenNoRecoMuon",
                data = unmatched_gen_muons
            )

            # gen ele not reco
            gen_eles = events.GenPart[(abs(events.GenPart.pdgId) == 11)]
            reco_eles = events.Electron

            unmatched_gen_mask = ak.fill_none(object_selections.delta_R(gen_eles, reco_eles, 0.4), True)

            unmatched_gen_eles = gen_eles[unmatched_gen_mask]
            unmatched_gen_eles = unmatched_gen_eles[unmatched_gen_eles.pt > self.options["electrons"]["pt"]]

            unmatched_gen_eles = awkward_utils.add_field(
                events = events,
                name = "SelectedGenNoRecoElectron",
                data = unmatched_gen_eles
            )

        # Photons
        photon_selection = self.select_photons(
                photons = events.Photon,
                options = self.options["photons"],
                electrons = electrons,
                rho = rho,
                year = self.year[:4]
        )
        photons = events.Photon[photon_selection]
        
        # lepton-photon overlap removal 
        clean_photon_mask = ak.fill_none(object_selections.delta_R(photons, muons, 0.3), True) & ak.fill_none(object_selections.delta_R(photons, electrons, 0.3), True)
        # object_selections.delta_R(photons, muons, 0.3) & object_selections.delta_R(photons, electrons, 0.3)
        photons = photons[clean_photon_mask]
        photons = ak.with_field(photons, ak.ones_like(photons.pt) * 0.0, "mass")

        # photons = awkward_utils.add_field(
        #     events = events,
        #     name = "SelectedPhoton",
        #     data = photons[:, :2]
        # )

        # # Jets
        # jet_cut = jet_selections.select_jets(
        #     jets = events.Jet,
        #     options = self.options["jets"],
        #     clean = {
        #         "photons" : {
        #             "objects" : photons,
        #             "min_dr" : self.options["jets"]["dr_photons"]
        #         },
        #         "electrons" : {
        #             "objects" : electrons,
        #             "min_dr" : self.options["jets"]["dr_electrons"]
        #         },
        #         "muons" : {
        #             "objects" : muons,
        #             "min_dr" : self.options["jets"]["dr_muons"]
        #         }
        #     },
        #     year = self.year,
        #     name = "SelectedJet",
        #     tagger = self
        # )
        # jets = awkward_utils.add_field(
        #     events = events,
        #     name = "SelectedJet",
        #     data = events.Jet[jet_cut]
        # )

        FSRphoton_selection = self.select_FSRphotons(
                FSRphotons = events.FsrPhoton,
                electrons = electrons,
                muons = muons,
                photons = photons,
                options = self.options["FSRphotons"]
        )
        FSRphotons = awkward_utils.add_field(
            events = events,
            name = "SelectedFSRPhotons",
            data = events.FsrPhoton[FSRphoton_selection]
        )
        FSRphotons = ak.with_field(FSRphotons, ak.ones_like(FSRphotons.pt) * 0.0, "mass")

        # if "2017" in self.year or "2018" in self.year:
        #     year = self.year[:4]
        #     b_jet_cut = jets.btagDeepFlavB > self.options["btag_med"][year]
        # else:
        #     b_jet_cut = jets.btagDeepFlavB > self.options["btag_med"][self.year]
        # jets = ak.with_field(jets, b_jet_cut, "is_med_bjet") 

        # Add object fields to events array
        for objects, name in zip([electrons, muons], ["electron", "muon"]):
            awkward_utils.add_object_fields(
                events = events,
                name = name,
                objects = objects,
                n_objects = 4,
                dummy_value = DUMMY_VALUE
            )
        
        if not self.is_data:
            dZ = events.GenVtx_z - events.PV_z
            awkward_utils.add_field(events, "dZ", dZ, overwrite=True)

        n_electrons = ak.fill_none(ak.num(electrons), 0)
        # N_e_cut = n_electrons>=2
        awkward_utils.add_field(events, "n_electrons", n_electrons, overwrite=True)

        n_muons = ak.num(muons)
        # N_mu_cut = n_muons>=2
        awkward_utils.add_field(events, "n_muons", n_muons, overwrite=True)

        n_leptons = n_electrons + n_muons
        # N_e_mu_cut = N_e_cut | N_mu_cut
        awkward_utils.add_field(events, "n_leptons", n_leptons, overwrite=True)

        # n_jets = ak.num(jets)
        # # logger.debug(f"Number of jets(tagger): {n_jets[:10]}")
        # awkward_utils.add_field(events, "n_jets", n_jets, overwrite=True)

        # n_b_jets = ak.sum(b_jet_cut, axis=1)
        # awkward_utils.add_field(events, "n_b_jets", n_b_jets, overwrite=True)

        n_photons = ak.num(photons)
        logger.debug(f"Number of photons(tagger): {n_photons[:10]}")
        awkward_utils.add_field(events, "n_photons", n_photons, overwrite=True)

        # PDG ID
        electrons = ak.with_field(electrons, ak.ones_like(electrons.pt) * 11, "id")
        # electrons = ak.with_field(electrons, ak.ones_like(electrons.pt) * 0.00051099895, "mass")
        muons = ak.with_field(muons, ak.ones_like(muons.pt) * 13, "id")

        # leptons ptE_error
        electrons = ak.with_field(electrons, electrons.energyErr, "ptE_error")
        muons = ak.with_field(muons, muons.ptErr, "ptE_error")

        # Sort objects by pt
        photons = photons[ak.argsort(photons.pt, ascending=False, axis=1)]
        electrons = electrons[ak.argsort(electrons.pt, ascending=False, axis=1)]
        muons = muons[ak.argsort(muons.pt, ascending=False, axis=1)]

        # self.select_fake_and_medium_photons(events=events, photons=photons)

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        electrons = ak.Array(electrons, with_name = "Momentum4D")
        muons = ak.Array(muons, with_name = "Momentum4D")
        photons = ak.Array(photons, with_name = "Momentum4D")
        FSRphotons = ak.Array(FSRphotons, with_name = "Momentum4D")

        awkward_utils.add_field(
                events = events,
                name = "SelectedPhoton",
                data = photons,
        )

        awkward_utils.add_object_fields(
            events = events,
            name = "gamma_fsr",
            objects = FSRphotons,
            n_objects = 1,
            dummy_value = DUMMY_VALUE
        )

        awkward_utils.add_field(events, "n_fsr", ak.num(FSRphotons), overwrite=True)
        logger.debug(f"Number of FSR photons: {ak.num(FSRphotons)}")
        logger.debug(f"Total number of FSR photons: {sum(ak.num(FSRphotons)>0)}")


        # 未修正的 Z boson 重建（不包含 FSR）
        ee_pairs_noFSR = ak.combinations(electrons, 2, fields=["LeadLepton", "SubleadLepton"])
        os_cut_noFSR = (ee_pairs_noFSR.LeadLepton.charge * ee_pairs_noFSR.SubleadLepton.charge) == -1
        ee_pairs_noFSR = ee_pairs_noFSR[os_cut_noFSR]
        mm_pairs_noFSR = ak.combinations(muons, 2, fields=["LeadLepton", "SubleadLepton"])
        os_cut_noFSR = (mm_pairs_noFSR.LeadLepton.charge * mm_pairs_noFSR.SubleadLepton.charge) == -1
        mm_pairs_noFSR = mm_pairs_noFSR[os_cut_noFSR]
        z_cands_noFSR = ak.concatenate([ee_pairs_noFSR, mm_pairs_noFSR], axis=1)
        z_cands_noFSR["ZCand"] = z_cands_noFSR.LeadLepton + z_cands_noFSR.SubleadLepton
        z_cands_noFSR = z_cands_noFSR[ak.argsort(abs(z_cands_noFSR.ZCand.mass - 91.1876), axis=1)]
        z_ee_cut_noFSR = ak.fill_none(ak.firsts(z_cands_noFSR).LeadLepton.id == 11, False)
        z_mumu_cut_noFSR = ak.fill_none(ak.firsts(z_cands_noFSR).LeadLepton.id == 13, False)

        # 修正 electrons 和 muons（包含 FSR）
        # electrons_withFSR = self.assign_fsr_photon(electrons, FSRphotons)
        electrons_withFSR = electrons
        muons_withFSR = self.assign_fsr_photon(muons, FSRphotons)

        # 修正的 Z boson 重建（包含 FSR）
        ee_pairs = ak.combinations(electrons_withFSR, 2, fields = ["LeadLepton", "SubleadLepton"])
        os_cut = (ee_pairs.LeadLepton.charge * ee_pairs.SubleadLepton.charge) == -1
        ee_pairs = ee_pairs[os_cut]
        mm_pairs = ak.combinations(muons_withFSR, 2, fields = ["LeadLepton", "SubleadLepton"])
        os_cut = (mm_pairs.LeadLepton.charge * mm_pairs.SubleadLepton.charge) == -1
        mm_pairs = mm_pairs[os_cut]
        z_cands = ak.concatenate([ee_pairs, mm_pairs], axis = 1)
        z_cands["ZCand"] = z_cands.LeadLepton + z_cands.SubleadLepton
        z_cands = z_cands[ak.argsort(abs(z_cands.ZCand.mass - 91.1876), axis = 1)]
        z_ee_cut = ak.fill_none(ak.firsts(z_cands).LeadLepton.id == 11, False)
        z_mumu_cut = ak.fill_none(ak.firsts(z_cands).LeadLepton.id == 13, False)

        # Make trigger cuts 
        if self.year is not None:
            year = self.year[:4]
            single_ele_trigger_cut = ak.num(events.Photon) < 0 # dummy cut, all False
            double_ele_trigger_cut = ak.num(events.Photon) < 0
            single_mu_trigger_cut = ak.num(events.Photon) < 0
            double_mu_trigger_cut = ak.num(events.Photon) < 0 
            for hlt in self.options["single_ele_trigger"][year]:
                if hasattr(events, hlt):
                    single_ele_trigger_cut = (single_ele_trigger_cut) | (events[hlt] == True)
                else:
                    logger.debug("[ZGammaTagger] %s is not in these event" % (hlt))
            for hlt in self.options["double_ele_trigger"][year]: # logical OR of all triggers
                if hasattr(events, hlt):
                    double_ele_trigger_cut = (double_ele_trigger_cut) | (events[hlt] == True)
                else:
                    logger.debug("[ZGammaTagger] %s is not in these event" % (hlt))
            for hlt in self.options["single_muon_trigger"][year]: # logical OR of all triggers
                if hasattr(events, hlt):
                    single_mu_trigger_cut = (single_mu_trigger_cut) | (events[hlt] == True)
                else:
                    logger.debug("[ZGammaTagger] %s is not in these event" % (hlt))
            for hlt in self.options["double_muon_trigger"][year]: # logical OR of all triggers
                if hasattr(events, hlt):
                    double_mu_trigger_cut = (double_mu_trigger_cut) | (events[hlt] == True)
                else:
                    logger.debug("[ZGammaTagger] %s is not in these event" % (hlt))
        else:
            single_ele_trigger_cut = ak.num(events.Photon) >= 0 # dummy cut, all True
            double_ele_trigger_cut = ak.num(events.Photon) >= 0
            single_mu_trigger_cut = ak.num(events.Photon) >= 0
            double_mu_trigger_cut = ak.num(events.Photon) >= 0
            
        if "2017" in self.year:
            single_ele_trigger_cut = ak.any((events.TrigObj.id == 11) & ((events.TrigObj.filterBits & 0x400) != 0), axis=1) & single_ele_trigger_cut

        trigger_cut = single_ele_trigger_cut | double_ele_trigger_cut | single_mu_trigger_cut  | double_mu_trigger_cut
        ele_trigger_cut = single_ele_trigger_cut | double_ele_trigger_cut
        mu_trigger_cut = single_mu_trigger_cut  | double_mu_trigger_cut
        # HLT lepton status: [FIXME]
        # 0: failed all triggers ->  double_ele_trigger_cut = False, single_ele_trigger_cut = False
        # 1: passed lower dilepton trigger -> double_ele_trigger_cut = False, single_ele_trigger_cut = True
        # 2: passed upper dilepton trigger -> double_ele_trigger_cut = True, single_ele_trigger_cut = False
        # 3: passed single lepton trigger ->  single_ele_trigger_cut = True
        # HLT_ele_cat0 = (~double_ele_trigger_cut) & (~single_ele_trigger_cut) 
        # HLT_mu_cat0 = (~double_mu_trigger_cut) & (~single_mu_trigger_cut)
        # HLT_ele_cat1 = (~double_ele_trigger_cut) & single_ele_trigger_cut
        # HLT_mu_cat1 = (~double_mu_trigger_cut) & single_mu_trigger_cut
        # HLT_ele_cat2 = double_ele_trigger_cut & (~single_ele_trigger_cut)
        # HLT_mu_cat2 = double_mu_trigger_cut & (~single_mu_trigger_cut)
        # HLT_ele_cat3 = single_ele_trigger_cut
        # HLT_mu_cat3 = single_mu_trigger_cut
        # def GetLeptonProbability(lepton_pt, lepton_eta, is_data, is_electron, trigger_leg) :
        #     if (is_electron):
        #         if (trigger_leg == pass_lowerdilep):
        #             if is_data:
        #                 prob_map = diele12_correction_data
        #                 unc_map = diele12_uncertainty_data
        #             else :
        #                 prob_map = diele12_correction_mc
        #                 unc_map = diele12_uncertainty_mc
        #         elif (trigger_leg == pass_upperdilep):
        #             if is_data:
        #                 prob_map = diele23_correction_data
        #                 unc_map = diele23_uncertainty_data
        #             else :
        #                 prob_map = diele23_correction_mc
        #                 unc_map = diele23_uncertainty_mc
        #         else :  
        #             if is_data:
        #                 prob_map = singleele_correction_data
        #                 unc_map = singleele_uncertainty_data
        #             else :
        #                 prob_map = singleele_correction_mc
        #                 unc_map = singleele_uncertainty_mc
        #     else:
        #         if (trigger_leg == pass_lowerdilep):
        #             if is_data:
        #                 prob_map = dimu12_correction_data
        #                 unc_map = dimu12_uncertainty_data
        #             else :
        #                 prob_map = dimu12_correction_mc
        #                 unc_map = dimu12_uncertainty_mc
        #         elif (trigger_leg == pass_upperdilep):
        #             if is_data:
        #                 prob_map = dimu23_correction_data
        #                 unc_map = dimu23_uncertainty_data
        #             else :
        #                 prob_map = dimu23_correction_mc
        #                 unc_map = dimu23_uncertainty_mc
        #         else :  
        #             if is_data:
        #                 prob_map = singlemu_correction_data
        #                 unc_map = singlemu_uncertainty_data
        #             else :
        #                 prob_map = singlemu_correction_mc
        #                 unc_map = singlemu_uncertainty_mc
        #     prob = prob_map
        #     uncr = unc_map
        #     return prob, uncr

        # def GetFlavorProbability(lepton_pt, lepton_eta, pass_singlelep, pass_dilep, is_data, is_electron):
        #     n_pass_lower = ak.num(1*HLT_ele_cat1[1*HLT_ele_cat1>0])+ak.num(1*HLT_mu_cat1[1*HLT_mu_cat1>0])
        #     n_pass_upper = ak.num(1*HLT_ele_cat2[1*HLT_ele_cat2>0])+ak.num(1*HLT_mu_cat2[1*HLT_mu_cat2>0])
        #     n_pass_single = ak.num(1*HLT_ele_cat3[1*HLT_ele_cat3>0])+ak.num(1*HLT_mu_cat3[1*HLT_mu_cat3>0])
        #     relevant_cat = ak.ones_like(lepton_pt) * True
        #     relevant_cat = ak.where(((n_pass_single>0)!=pass_singlelep), ak.ones_like(relevant_cat) * False, relevant_cat)
        #     relevant_cat = ak.where((((n_pass_upper>0)&(n_pass_lower>1))!=pass_dilep), ak.ones_like(relevant_cat) * False, relevant_cat)
        #     lep_prob = ak.zero_like(lepton_pt) 
        #     lep_unc = ak.zero_like(lepton_pt)

        # def GetTotalProbability(electron_pt,muon_pt, electron_eta,muon_eta, pass_singleel,pass_singlemu,pass_diel,pass_dimu, is_data):
        #     electron_prob = GetFlavorProbability(electron_pt, electron_eta,pass_singleel, pass_diel, is_data, True)
        #     muon_prob = GetFlavorProbability(muon_pt, muon_eta,pass_singlemu, pass_dimu, is_data, False)
        
        # mc_prob = GetTotalProbability(electron_pt, muon_pt, electron_eta, muon_eta, pass_singleel, pass_singlemu, pass_diel, pass_dimu, False)


        # 0: lep_prob = 1 - probability(1), lep_unc = uncertainty(1)
        # 1: lep_prob = prob(1) - prob(2), lep_unc = sqrt(uncertainty(1)**2 + uncertainty(2)**2)
        # 2: lep_prob = prob(2) - prob(3), lep_unc = sqrt(uncertainty(2)**2 + uncertainty(3)**2)
        # 3: lep_prob = prob(3), lep_unc = uncertainty(3)
        
        #if mc_prob < 0.001: mc_prob = 0.
        # sf = 1.0;
        # unc = 0.0;
        # propagate_uncertainty_ratio(data_prob, data_unc, mc_prob, mc_unc, sf, unc);
        # sf,unc=1
        # if mc_prob !=0: sf = data_prob/mc_prob; unc = sqrt((data_unc/mc_prob)**2 + ((mc_unc*data_prob)/(mc_prob*mc*prob))**2)
        # elif mc_prob == 0 sf = 1; unc = 0
        #[FIXME]
        # a trick that give 0 to the empty array    
        if self.year is not None:
            year = self.year[:4]
            e_cut = ak.fill_none(ak.pad_none(electrons.pt, 1, axis=1)[:, 0], 0) > self.options["lead_ele_pt"][year]
            m_cut = ak.fill_none(ak.pad_none(muons.pt, 1, axis=1)[:, 0], 0) > self.options["lead_mu_pt"][year]
        else:
            e_cut = ak.fill_none(ak.pad_none(electrons.pt, 1, axis=1)[:, 0], 0) > 25
            m_cut = ak.fill_none(ak.pad_none(muons.pt, 1, axis=1)[:, 0], 0) > 20
        ee_cut = (ak.fill_none(ak.pad_none(electrons.pt, 1, axis=1)[:, 0], 0) > 25) & (ak.fill_none(ak.pad_none(electrons.pt, 2, axis=1)[:, 1], 0) > 15)
        mm_cut = (ak.fill_none(ak.pad_none(muons.pt, 1, axis=1)[:, 0], 0) > 20) & (ak.fill_none(ak.pad_none(muons.pt, 2, axis=1)[:, 1], 0) > 10)
        
        ele_trigger_pt_cut = (single_ele_trigger_cut & e_cut) | (double_ele_trigger_cut & ee_cut)
        mu_trigger_pt_cut = (single_mu_trigger_cut & m_cut) | (double_mu_trigger_cut & mm_cut)
        trigger_pt_cut = (single_ele_trigger_cut & e_cut) | (double_ele_trigger_cut & ee_cut) | (single_mu_trigger_cut & m_cut) | (double_mu_trigger_cut & mm_cut)
        
        z_mumu = z_mumu_cut 
        z_ee = z_ee_cut
        awkward_utils.add_field(events, "z_mumu", z_mumu, overwrite=True)
        awkward_utils.add_field(events, "z_ee", z_ee, overwrite=True)

        # mass_cut = (z_cands.ZCand.mass > 80.) & (z_cands.ZCand.mass < 100.)
        mass_cut = z_cands.ZCand.mass > 50.
        z_cands = z_cands[mass_cut] # OSSF lepton pairs with m_ll > 50.
        z_cands_noFSR = z_cands_noFSR[mass_cut]

        # # HEM cut
        # if self.year=="2018" and self.is_data:
        #     hem_run=events.run > 319077        
        #     # checked 65.15623538907509% events in data could pass this run cut
        #     hem_jet=ak.num(events.Jet[(events.Jet.phi>-1.57) & (events.Jet.phi<-0.87) & (events.Jet.eta>-3) & (events.Jet.eta<-1.3)])>0
        #     hem_fatjet=ak.num(events.FatJet[(events.FatJet.phi>-1.57) & (events.FatJet.phi<-0.87) & (events.FatJet.eta>-3) & (events.FatJet.eta<-1.3)])>0
        #     hem_cut=~((hem_run & hem_jet) | (hem_run & hem_fatjet))        
        # elif self.year=="2018" and not self.is_data:
        #     #random number generator from 0 to 1
        #     fraction=0.6515623538907509
        #     events['random'] = numpy.random.rand(len(events))
        #     hem_run=events.random < fraction
        #     hem_jet=ak.num(events.Jet[(events.Jet.phi>-1.57) & (events.Jet.phi<-0.87) & (events.Jet.eta>-3) & (events.Jet.eta<-1.3)])>0
        #     hem_fatjet=ak.num(events.FatJet[(events.FatJet.phi>-1.57) & (events.FatJet.phi<-0.87) & (events.FatJet.eta>-3) & (events.FatJet.eta<-1.3)])>0
        #     hem_cut=~((hem_run & hem_jet) | (hem_run & hem_fatjet))
        # else:
        #     hem_cut=ak.num(events.Photon) >= 0 
        # events = events[hem_cut]
        
        # # Construct di-electron/di-muon pairs
        # ee_pairs = ak.combinations(electrons, 2, fields = ["LeadLepton", "SubleadLepton"])
        # if self.year is not None:
        #     year = self.year[:4]
        #     e_cut = ee_pairs.LeadLepton.pt > self.options["lead_ele_pt"][year]
        # else:
        #     e_cut = ee_pairs.LeadLepton.pt > 25
        # ee_cut = (ee_pairs.LeadLepton.pt > 25) & (ee_pairs.SubleadLepton.pt > 15)
        # ee_pairs = ak.concatenate([ee_pairs[double_ele_trigger_cut & ee_cut], ee_pairs[single_ele_trigger_cut & e_cut]], axis = 1)
        # ele_pair_cut = ak.num(ee_pairs) >= 1

        # mm_pairs = ak.combinations(muons, 2, fields = ["LeadLepton", "SubleadLepton"])
        # if self.year is not None:
        #     year = self.year[:4]
        #     m_cut = mm_pairs.LeadLepton.pt > self.options["lead_mu_pt"][year]
        # else:
        #     m_cut = mm_pairs.LeadLepton.pt > 20
        # mm_cut = (mm_pairs.LeadLepton.pt > 20) & (mm_pairs.SubleadLepton.pt > 10)
        # mm_pairs = ak.concatenate([mm_pairs[double_mu_trigger_cut & mm_cut], mm_pairs[single_mu_trigger_cut & m_cut]], axis = 1)
        # muon_pair_cut = ak.num(mm_pairs) >= 1

        # # Concatenate these together
        # z_cands = ak.concatenate([ee_pairs, mm_pairs], axis = 1)

        # # Make Z candidate-level cuts
        # # print('!!!!charge of lead lepton{}, charge of sublead lepton {}, flavour of lepton {}'.format(ak.flatten(z_cands[os_cut]).LeadLepton.charge, ak.flatten(z_cands[os_cut]).SubleadLepton.charge, ak.flatten(z_cands[os_cut]).LeadLepton.id))
        # os_cut = (z_cands.LeadLepton.charge * z_cands.SubleadLepton.charge) == -1
        # z_cands = z_cands[os_cut]
        # os_cut = ak.num(z_cands) >= 1

        # z_cands["ZCand"] = z_cands.LeadLepton + z_cands.SubleadLepton # these add as 4-vectors since we registered them as "Momentum4D" objects
        # mass_cut = (z_cands.ZCand.mass > 80.) & (z_cands.ZCand.mass < 100.)
        # # mass_cut = (z_cands.ZCand.mass > 50.)
        # z_cands = z_cands[mass_cut] # OSSF lepton pairs with m_ll > 50.

        # z_cands = z_cands[ak.argsort(abs(z_cands.ZCand.mass - 91.1876), axis = 1)] # take the one with mass closest to mZ

        has_z_cand = ak.num(z_cands) >= 1
        z_cand = ak.firsts(z_cands)
        z_cand_noFSR = ak.firsts(z_cands_noFSR)

        # Add Z-related fields to array
        for field in ["pt", "eta", "phi", "mass", "charge", "id", "ptE_error"]:
            if not field in ["charge", "id", "ptE_error"]:
                # Include FSR
                awkward_utils.add_field(
                        events,
                        "Z_%s" % field,
                        ak.fill_none(getattr(z_cand.ZCand, field), DUMMY_VALUE)
                )
                # Without FSR
                awkward_utils.add_field(
                    events,
                    f"Z_noFSR_{field}",
                    ak.fill_none(getattr(z_cand_noFSR.ZCand, field), DUMMY_VALUE)
                )
            awkward_utils.add_field(
                    events,
                    "Z_lead_lepton_%s" % field,
                    ak.fill_none(z_cand.LeadLepton[field], DUMMY_VALUE)
            )
            awkward_utils.add_field(
                    events,
                    "Z_sublead_lepton_%s" % field,
                    ak.fill_none(z_cand.SubleadLepton[field], DUMMY_VALUE)
            )
            awkward_utils.add_field(
                events,
                f"Z_noFSR_lead_lepton_{field}",
                ak.fill_none(z_cand_noFSR.LeadLepton[field], DUMMY_VALUE)
            )
            awkward_utils.add_field(
                events,
                f"Z_noFSR_sublead_lepton_{field}",
                ak.fill_none(z_cand_noFSR.SubleadLepton[field], DUMMY_VALUE)
            )

        # Make gamma candidate-level cuts
        has_2gamma_cand = (ak.num(photons) >= 2) #& (events.n_iso_photons == 0) # only for dy samples

        gamma_pairs = ak.combinations(photons, 2, fields=["LeadPhoton", "SubleadPhoton"])
        gamma_pairs = gamma_pairs[ak.argsort(gamma_pairs.LeadPhoton.pt, ascending=False, axis=1)]
        
        gamma_pairs["LeadPhoton"] = ak.with_name(gamma_pairs.LeadPhoton, "Momentum4D")
        gamma_pairs["SubleadPhoton"] = ak.with_name(gamma_pairs.SubleadPhoton, "Momentum4D")

        alp_cand = ak.firsts(gamma_pairs)

        # alp_cand["ALPCand"] = alp_cand.LeadPhoton + alp_cand.SubleadPhoton
        alp_cand["ALPCand"] = ak.with_name(alp_cand.LeadPhoton + alp_cand.SubleadPhoton, "Momentum4D")

        events = self.calculate_alp_photon_isolation(events, alp_cand, photons)

        # Add ALP-related fields
        for field in ["pt", "eta", "phi", "mass", "energyErr", "r9", "sieie", "hoe_PUcorr", "hcalPFClusterIso", "ecalPFClusterIso"]:
            if not field in ["energyErr", "r9", "sieie", "hoe_PUcorr", "hcalPFClusterIso", "ecalPFClusterIso"]:
                awkward_utils.add_field(
                    events,
                    "ALP_%s" % field,
                    ak.fill_none(getattr(alp_cand.ALPCand, field), DUMMY_VALUE)
                )
            awkward_utils.add_field(
                events,
                "ALP_lead_photon_%s" % field,
                ak.fill_none(alp_cand.LeadPhoton[field], DUMMY_VALUE)
            )
            awkward_utils.add_field(
                events,
                "ALP_sublead_photon_%s" % field,
                ak.fill_none(alp_cand.SubleadPhoton[field], DUMMY_VALUE)
            )
        if int(self.year[:4]) < 2020:
            awkward_utils.add_field(events, "ALP_lead_photon_chiso",  alp_cand.LeadPhoton.pfRelIso03_chg) #run2
            awkward_utils.add_field(events, "ALP_lead_photon_alliso", alp_cand.LeadPhoton.pfRelIso03_all) #run2
            awkward_utils.add_field(events, "ALP_sublead_photon_chiso",  alp_cand.SubleadPhoton.pfRelIso03_chg) #run2
            awkward_utils.add_field(events, "ALP_sublead_photon_alliso", alp_cand.SubleadPhoton.pfRelIso03_all) #run2
        elif int(self.year[:4]) > 2020:
            awkward_utils.add_field(events, "ALP_lead_photon_chiso",  alp_cand.LeadPhoton.pfRelIso03_chg_quadratic) #run3
            awkward_utils.add_field(events, "ALP_lead_photon_alliso",  alp_cand.LeadPhoton.pfRelIso03_all_quadratic) #run3
            awkward_utils.add_field(events, "ALP_sublead_photon_chiso",  alp_cand.SubleadPhoton.pfRelIso03_chg_quadratic) #run3
            awkward_utils.add_field(events, "ALP_sublead_photon_alliso",  alp_cand.SubleadPhoton.pfRelIso03_all_quadratic) #run3
        
        # Gamma candidate 
        gamma_cand = ak.firsts(photons)
        gamma_mvaID_WPL = ((gamma_cand.isScEtaEB & (gamma_cand.mvaID > self.options["photons"]["mvaID_barrel"])) | (gamma_cand.isScEtaEE & (gamma_cand.mvaID > self.options["photons"]["mvaID_endcap"])))
        gamma_e_veto = gamma_cand.electronVeto > self.options["photons"]["e_veto"]

        awkward_utils.add_field(gamma_cand, "mass", ak.ones_like(gamma_cand.pt) * 0) #TODO: run3 BUG

        # Add gamma-related fields to array
        for field in ["pt", "eta", "phi", "mass", "mvaID", "energyErr", "sieie", "hoe", "r9", "mvaID_WP80", "mvaID_WP90"]:
            awkward_utils.add_field(
                events,
                "gamma_%s" % field,
                ak.fill_none(getattr(gamma_cand, field), DUMMY_VALUE)
            )
        awkward_utils.add_field(events, "gamma_mvaID_WPL",  gamma_mvaID_WPL)
        awkward_utils.add_field(events, "gamma_e_veto",  gamma_e_veto)
        if int(self.year[:4]) < 2020:
            awkward_utils.add_field(events, "gamma_chiso",  gamma_cand.pfRelIso03_chg) #run2
            awkward_utils.add_field(events, "gamma_alliso",  gamma_cand.pfRelIso03_all) #run2
        elif int(self.year[:4]) > 2020:
            awkward_utils.add_field(events, "gamma_chiso",  gamma_cand.pfRelIso03_chg_quadratic) #run3
            awkward_utils.add_field(events, "gamma_alliso",  gamma_cand.pfRelIso03_all_quadratic) #run3
        #awkward_utils.add_field(events, "gamma_mvaID_17",  gamma_cand.mvaID_Fall17V2) #run3


        # Make Higgs candidate-level cuts
        # Use FSR Z boson 
        h_cand = (z_cand.ZCand + alp_cand.ALPCand)
        sel_h_1 = (z_cand.ZCand.mass + h_cand.mass) > options["mass_sum"]
        sel_h_2 = (h_cand.mass > options["mass_h"][0]) & (h_cand.mass < options["mass_h"][1])
        sel_h_1 = ak.fill_none(sel_h_1, value = False)
        sel_h_2 = ak.fill_none(sel_h_2, value = False)
        
        # Use No FSR Z boson 
        h_cand_noFSR = z_cand_noFSR.ZCand + alp_cand.ALPCand
        sel_h_1_noFSR = (z_cand_noFSR.ZCand.mass + h_cand_noFSR.mass) > options["mass_sum"]
        sel_h_2_noFSR = (h_cand_noFSR.mass > options["mass_h"][0]) & (h_cand_noFSR.mass < options["mass_h"][1])
        sel_h_1_noFSR = ak.fill_none(sel_h_1_noFSR, value=False)
        sel_h_2_noFSR = ak.fill_none(sel_h_2_noFSR, value=False)

        print(f'!!!!has H: {sum(has_z_cand & has_2gamma_cand & z_mumu_cut)} | {sum(has_z_cand & has_2gamma_cand & z_ee_cut)}')

        # Add Higgs-related fields to array
        for field in ["pt", "eta", "phi", "mass"]:
            awkward_utils.add_field(
                events,
                "H_%s" % field,
                ak.fill_none(getattr(h_cand, field), DUMMY_VALUE)
            )
            # 不包含 FSR
            awkward_utils.add_field(
                events,
                f"H_noFSR_{field}",
                ak.fill_none(getattr(h_cand_noFSR, field), DUMMY_VALUE)
            )

        # additional leptons
        leptons = ak.concatenate([electrons, muons], axis = 1)
        max_I_mini = ak.fill_none(ak.max(leptons.miniPFRelIso_all, axis = 1), 9999)
        awkward_utils.add_field(events, "max_I_mini", max_I_mini)
        
        veto_Z_leptons = (leptons.pt != events.Z_lead_lepton_pt) & (leptons.pt != events.Z_sublead_lepton_pt)
        additional_leptons = leptons[veto_Z_leptons]       
        additional_leptons = additional_leptons[ak.argsort(additional_leptons.pt, ascending=False, axis=1)]

        for objects, name in zip([additional_leptons], ["additional_lepton"]):
            awkward_utils.add_object_fields(
                events = events,
                name = name,
                objects = objects,
                n_objects = 2,
                dummy_value = DUMMY_VALUE
            )

        event_filter = (events.Flag_goodVertices & 
                        events.Flag_globalSuperTightHalo2016Filter & 
                        ((ak.num(events.Photon) >= 0) if "202" in self.year else events.Flag_HBHENoiseFilter) & 
                        ((ak.num(events.Photon) >= 0) if "202" in self.year else events.Flag_HBHENoiseIsoFilter) & 
                        events.Flag_EcalDeadCellTriggerPrimitiveFilter & 
                        events.Flag_BadPFMuonFilter & 
                        events.Flag_BadPFMuonDzFilter & 
                        events.Flag_hfNoisyHitsFilter & 
                        events.Flag_eeBadScFilter & 
                        ((ak.num(events.Photon) >= 0) if "2016" in self.year else events.Flag_ecalBadCalibFilter) # 2016 dummy cut, all True
                        )
        
        all_cuts = trigger_pt_cut & has_z_cand & has_2gamma_cand & sel_h_1 & sel_h_2 & event_filter #& ak.fill_none((h_cand.mass>80) & (h_cand.mass < options["mass_h"][1]), False)

        for cut_type in ["zgammas", "zgammas_ele", "zgammas_mu", "zgammas_w", "zgammas_ele_w", "zgammas_mu_w"]:
            if "_w" in cut_type:
                if hasattr(events, 'Generator_weight'):
                    weighted = True
                else:
                    continue
            else:
                weighted = False

            cut0 = ak.num(events.Photon) >= 0

            cut1 = z_ee_cut | z_mumu_cut
            if "ele" in cut_type:
                cut1 = z_ee_cut
            elif "mu" in cut_type:
                cut1 = z_mumu_cut
            cut2 = cut1 & trigger_cut
            if "ele" in cut_type:
                cut2 = cut1 & ele_trigger_cut
            elif "mu" in cut_type:
                cut2 = cut1 & mu_trigger_cut
            cut3 = cut2 & trigger_pt_cut
            if "ele" in cut_type:
                cut3 = cut2 & ele_trigger_pt_cut
            elif "mu" in cut_type:
                cut3 = cut2 & mu_trigger_pt_cut
            cut4 = cut3 & has_z_cand
            cut5 = cut4 & has_2gamma_cand
            cut6 = cut5 & sel_h_1
            cut7 = cut6 & sel_h_2
            cut8 = cut7 & event_filter
            
            if cut_type == "zgammas_ele":
                ee_all_cut = cut8
            if cut_type == "zgammas_mu":
                mm_all_cut = cut8
            
            # if cut_type == "zgammas_ele":
            #     print(f"!!!start check events tag({cut_type})!!!")
            #     for i in events[cut1]:
            #         print(f"{i.run} {i.luminosityBlock} {i.event}")
            #     print(f"!!!end check events tag({cut_type})!!!")

            self.register_event_cuts(
                names = ["all", "N_lep_sel", "trig_cut", "lep_pt_cut", "has_z_cand", "has_2g_cand", "sel_h_1", "sel_h_2", "event", "all cuts"],
                results = [cut0, cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8, all_cuts],
                events = events,
                cut_type = cut_type,
                weighted = weighted
            )
            # cut_names = ["N_lep_sel", "trig_cut", "lep_pt_cut", "has_g_cand", "os_cut", "has_z_cand", "sel_h_1", "sel_h_2", "sel_h_3"]
            # for cut, cut_name in zip([cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8, cut9], cut_names):
            #     awkward_utils.add_field(events, f"{cut_type}_{cut_name}", cut)

        all_cuts = ee_all_cut | mm_all_cut
        
        # print(f"Sum of all_cuts: {sum(all_cuts)}")
        # all_cuts = ee_all_cut | mm_all_cut
        # print(f"Sum of all_cuts: {sum(all_cuts)}")

        # checked_cut = (z_ee_cut | z_mumu_cut) & pair_cut
        # checked_events = events[checked_cut]
        # print("!!!start check events tag(inclusive)!!!")
        # for event in checked_events:
        #     print(event.run, event.luminosityBlock, event.event, sep=" ")
        # print("!!!end check events tag(inclusive)!!!")

        # checked_cut = z_ee_cut & ee_trigger_pt_cut
        # checked_events = events[checked_cut]
        # print("!!!start check events tag(electron)!!!")
        # for event in checked_events:
        #     print(event.run, event.luminosityBlock, event.event, sep=" ")
        # print("!!!end check events tag(electron)!!!")

        # checked_cut = z_mumu_cut & mm_trigger_pt_cut
        # checked_events = events[checked_cut]
        # print("!!!start check events tag(muon)!!!")
        # for event in checked_events:
        #     print(event.run, event.luminosityBlock, event.event, sep=" ")
        # print("!!!end check events tag(muon)!!!")

        # self.register_cuts(
        #     names = ["has_z_cand", "has_2gamma_cand", "sel_h_1", "sel_h_2", "sel_h_3", "all cuts"],
        #     results = [has_z_cand, has_2gamma_cand, sel_h_1, sel_h_2, sel_h_3, all_cuts],
        #     cut_type = "zgammas_unweighted"
        # )


        elapsed_time = time.time() - start
        logger.debug("[ZGammaTagger] %s, syst variation : %s, total time to execute select_zgammas: %.6f s" % (self.name, self.current_syst, elapsed_time))

        #dummy_cut =  ak.num(events.Photon) >= 0
        return all_cuts, events 


    # def select_fake_and_medium_photons(self, events, photons):
    #     # | pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho | 
    #     # listed from the right side
    #     mask1 = 0b10101010101010  # full medium ID
    #     mask2 = 0b00101010101010  # remove Isopho
    #     mask3 = 0b10001010101010  # remove IsoNeu
    #     mask4 = 0b10100010101010  # remove Isoch 
    #     mask5 = 0b10101000101010  # remove sigma ieie
    #     mask6 = 0b10100000101010  # remove the Isoch and sigma ieie

    #     # photons = photons[photons.pixelSeed]
    #     bitmap = photons.vidNestedWPBitmap

    #     # select medium and control photons
    #     # after adding the photons that pass the full ID, add the photons that pass the inverted ID
    #     # select control photons that don't pass the full ID but pass ID that one of cut inverted, which means this cut is inverted
    #     # also the fake photon enriched region
    #     # use photon_selection to identified the type of events
    #     medium_and_control_cut = ((bitmap & mask1) == mask1) | ((bitmap & (mask2 + (3<<12))) == mask2) | ((bitmap & (mask3 + (3<<10))) == mask3) | ((bitmap & (mask4 + (3<<8))) == mask4) | ((bitmap & mask6) == mask6)
    #     selected_medium_or_control_photons = photons[medium_and_control_cut] # append the medium and control photons

    #     pass_selection1 = ak.num(selected_medium_or_control_photons) >= 1  # select medium and control photons without fake photon
    #     awkward_utils.add_field(events, "pass_selection1", pass_selection1)  # has selected medium and control photons

    #     med_cand = ak.firsts(selected_medium_or_control_photons) #similar to the gamma_cand

    #     bitmap = med_cand.vidNestedWPBitmap 
    #     photon_selection = (
    #         (((bitmap & mask1) == mask1) << 0) + 
    #         (((bitmap & (mask2 + (3<<12))) == mask2) << 1) + 
    #         ((((bitmap & (mask3 + (3<<10)))) == mask3) << 2) + 
    #         (((bitmap & (mask4 + (3<<8))) == mask4) << 3) + 
    #         (((bitmap & (mask5 + (3<<6))) == mask5) << 4) + 
    #         (((bitmap & mask5) == mask5) << 5) + 
    #         ((((bitmap & mask6) == mask6) & ((bitmap & mask5) != mask5)) << 6)
    #     )
    #     awkward_utils.add_field(events, "photon_selection", ak.fill_none(photon_selection, -1))
    #     #pass_selection1 && photon_selection==1 or 5 or 7 -> build ture template from MC and data template from data
    #     #pass_selection1 && photon_selection==4 or 6 or 8  && chiso side band -> build fake tempalte from data
    #     #pass_selection1 && ((photon_selection!=1 && photon_selection==2) || (!1 && ==3) || (!1 && ==4) || (!=1 && ==5)) -> build non-prompt photon sample from data

    #     awkward_utils.add_field(events, "photon_is_barrel", ak.fill_none(med_cand.isScEtaEB, -1))
    #     awkward_utils.add_field(events, "photon_is_endcap", ak.fill_none(med_cand.isScEtaEE, -1))
    #     for field in ["pt", "eta", "phi", "sieie"]:
    #         awkward_utils.add_field(
    #             events,
    #             "photon_%s" % field,
    #             ak.fill_none(getattr(med_cand, field), DUMMY_VALUE)
    #         )
    #     awkward_utils.add_field(events, "photon_chiso",  ak.fill_none(med_cand.pfRelIso03_chg, DUMMY_VALUE))
    
    def calculate_gen_info(self, zgammas, options):
        """
        Calculate gen info, adding the following fields to the events array:
            GenHggHiggs : [pt, eta, phi, mass, dR]
            GenHggLeadPhoton : [pt, eta, phi, mass, dR, pt_diff]
            GenHggSubleadPhoton : [pt, eta, phi, mass, dR, pt_diff]
            LeadPhoton : [gen_dR, gen_pt_diff]
            SubleadPhoton : [gen_dR, gen_pt_diff]

        Perform both matching of
            - closest gen photons from Higgs to reco lead/sublead photons from diphoton candidate
            - closest reco photons to gen photons from Higgs

        If no match is found for a given reco/gen photon, it will be given values of -999. 
        """
        gen_hzg = gen_selections.select_x_to_yz(zgammas.GenPart, 25, 23, 22)

        print("DEBUG: bing", gen_hzg.fields)
        
        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgHiggs",
                objects = gen_hzg.GenParent,
                n_objects = 1
        )

        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgLeadGenChild",
                objects = gen_hzg.LeadGenChild,
                n_objects = 1
        )

        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgSubleadGenChild",
                objects = gen_hzg.SubleadGenChild,
                n_objects = 1
        )

        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgLeadGenChildChild1",
                objects = gen_hzg.LeadGenChildChild1,
                n_objects = 1
        )

        awkward_utils.add_object_fields(
                events = zgammas,
                name = "GenHzgLeadGenChildChild2",
                objects = gen_hzg.LeadGenChildChild2,
                n_objects = 1
        )

        return zgammas 
        

    def select_photons(self, photons, options, electrons, rho, year):
        """
        Enforces all photon cuts that are commmon to both
        leading and subleading photons in the diphoton preselection.
        Cuts specific to a diphoton pair are not enforced here.

        :param photons: input collection of photons to use for calculating selected photons
        :type photons: ak.highlevel.Array
        :param rho: energy density in each event, used for corrections to photon isolation cuts
        :type rho: ak.highlevel.Array
        :param options: dictionary containing configurable options for the photon selection
        :type options: dict
        :return: boolean array indicating which photons pass the photon selection
        :rtype: ak.highlevel.Array
        """
        tagger_name = "none" if self is None else self.name

        logger.debug("[select_objects] : Tagger '%s', selecting objects '%s', with the following requirements:" % (tagger_name, "SelectedPhoton"))
        for cut, value in options.items():
            logger.debug("\t '%s' : %s" % (cut, str(value)))

        # nominal
        no_cut = photons.pt > 0

        # pt
        pt_cut = photons.pt > options["pt"]

        # eta
        #eta_cut = Tagger.get_range_cut(abs(photons.eta), options["eta"]) | (photons.isScEtaEB | photons.isScEtaEE)
        eta_cut = (photons.isScEtaEB | photons.isScEtaEE) 
        # eta_cut = ((photons.isScEtaEB & (photons.mvaID > options["mvaID_barrel"])) | (photons.isScEtaEE & (photons.mvaID > options["mvaID_endcap"])))

        customized_id_cut = []
        official_id_cut = []
        sieie_cut = []
        PFECalIso_cut = []
        rho_broadcasted, _ = ak.broadcast_arrays(rho, photons.pt)
        rho = rho_broadcasted
        photon_abs_eta = numpy.abs(photons.eta)
        if int(year) < 2020:
            customized_id_cut = ak.ones_like(photons.pt) # all true, dummy TODO
            official_id_cut = ak.ones_like(photons.pt) # all true, dummy TODO
        elif int(year) > 2020:
            # ID
            # id_cut = photons.mvaID_WP80
            ''' 1. Rho corrected H/E '''
            hoe_barrel_cut = photons.hoe_PUcorr < options["hoe_barrel"]
            hoe_endcap_cut = photons.hoe_PUcorr < options["hoe_endcap"]

            ''' 2. Rho corrected PF charged hadron isolation '''
            PFChIso_barrel_cut = photons.pfRelIso03_chg_quadratic < options["PFChIso_barrel"]
            PFChIso_endcap_cut = photons.pfRelIso03_chg_quadratic < options["PFChIso_endcap"]

            ''' 3. Rho corrected PF HCal isolation '''
            # quadratic EA corrections in Run3 : https://indico.cern.ch/event/1204277/contributions/5064356/attachments/2538496/4369369/CutBasedPhotonID_20221031.pdf
            coef_1, coef_2, coef_3 = options["PFHCalIso_barrel"]
            parabola_cut = coef_1 + coef_2 * photons.pt + coef_3 * photons.pt ** 2
            # PFHCalIso_barrel_cut = photons.hcalPFClusterIso < parabola_cut
            PFHCalIso_barrel_cut = (
                ((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0))
                & (
                    (photons.hcalPFClusterIso
                    - rho * options["PFHCalIso_EA_EB_1"][0]
                    - rho**2 * options["PFHCalIso_EA_EB_1"][1])
                    < parabola_cut
                )
            ) | (
                ((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442))
                & (
                    (photons.hcalPFClusterIso
                    - rho * options["PFHCalIso_EA_EB_2"][0]
                    - rho**2 * options["PFHCalIso_EA_EB_2"][1])
                    < parabola_cut
                )
            )

            coef_1, coef_2, coef_3 = options["PFHCalIso_endcap"]
            parabola_cut = coef_1 + coef_2 * photons.pt + coef_3 * photons.pt ** 2
            # PFHCalIso_endcap_cut = photons.hcalPFClusterIso < parabola_cut
            PFHCalIso_endcap_cut = (
                (
                    ((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0))
                    & (
                        photons.hcalPFClusterIso
                        - (rho * options["PFHCalIso_EA_EE_1"][0])
                        - (rho**2 * options["PFHCalIso_EA_EE_1"][1])
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2))
                    & (
                        photons.hcalPFClusterIso
                        - (rho * options["PFHCalIso_EA_EE_2"][0])
                        - (rho**2 * options["PFHCalIso_EA_EE_2"][1])
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3))
                    & (
                        photons.hcalPFClusterIso
                        - (rho * options["PFHCalIso_EA_EE_3"][0])
                        - (rho**2 * options["PFHCalIso_EA_EE_3"][1])
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4))
                    & (
                        photons.hcalPFClusterIso
                        - (rho * options["PFHCalIso_EA_EE_4"][0])
                        - (rho**2 * options["PFHCalIso_EA_EE_4"][1])
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5))
                    & (
                        photons.hcalPFClusterIso
                        - (rho * options["PFHCalIso_EA_EE_5"][0])
                        - (rho**2 * options["PFHCalIso_EA_EE_5"][1])
                        < parabola_cut
                    )
                )
            )

            ''' 4. Sigma IEtaIEta (Discarded) '''
            sieie_barrel_cut = photons.sieie < options["sieie_barrel"]
            sieie_endcap_cut = photons.sieie < options["sieie_endcap"]

            ''' 5. Rho corrected PF ECal isolation (Discarded) '''
            # quadratic EA corrections in Run3 : https://indico.cern.ch/event/1204277/contributions/5064356/attachments/2538496/4369369/CutBasedPhotonID_20221031.pdf
            coef_1, coef_2 = options["PFECalIso_barrel"]
            parabola_cut = coef_1 + coef_2 * photons.pt
            # PFECalIso_barrel_cut = photons.ecalPFClusterIso < parabola_cut
            PFECalIso_barrel_cut = (
                ((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0))
                & (
                    (photons.ecalPFClusterIso
                    - rho * options["PFECalIso_EA_EB_1"][0]
                    - rho**2 * options["PFECalIso_EA_EB_1"][1])
                    < parabola_cut
                )
            ) | (
                ((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442))
                & (
                    (photons.ecalPFClusterIso
                    - rho * options["PFECalIso_EA_EB_2"][0]
                    - rho**2 * options["PFECalIso_EA_EB_2"][1])
                    < parabola_cut
                )
            )

            coef_1, coef_2 = options["PFECalIso_endcap"]
            parabola_cut = coef_1 + coef_2 * photons.pt
            # PFECalIso_endcap_cut = photons.ecalPFClusterIso < parabola_cut
            PFECalIso_endcap_cut = (
                (
                    ((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0))
                    & (
                        photons.ecalPFClusterIso
                        - (rho * options["PFECalIso_EA_EE_1"][0])
                        - (rho**2 * options["PFECalIso_EA_EE_1"][1])
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2))
                    & (
                        photons.ecalPFClusterIso
                        - (rho * options["PFECalIso_EA_EE_2"][0])
                        - (rho**2 * options["PFECalIso_EA_EE_2"][1])
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3))
                    & (
                        photons.ecalPFClusterIso
                        - (rho * options["PFECalIso_EA_EE_3"][0])
                        - (rho**2 * options["PFECalIso_EA_EE_3"][1])
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4))
                    & (
                        photons.ecalPFClusterIso
                        - (rho * options["PFECalIso_EA_EE_4"][0])
                        - (rho**2 * options["PFECalIso_EA_EE_4"][1])
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5))
                    & (
                        photons.ecalPFClusterIso
                        - (rho * options["PFECalIso_EA_EE_5"][0])
                        - (rho**2 * options["PFECalIso_EA_EE_5"][1])
                        < parabola_cut
                    )
                )
            )

            # 1, 2, 3
            hoe_cut = (photons.isScEtaEB & hoe_barrel_cut) | (photons.isScEtaEE & hoe_endcap_cut)
            PFChIso_cut = (photons.isScEtaEB & PFChIso_barrel_cut) | (photons.isScEtaEE & PFChIso_endcap_cut)
            PFHCalIso_cut = (photons.isScEtaEB & PFHCalIso_barrel_cut) | (photons.isScEtaEE & PFHCalIso_endcap_cut)

            customized_id_cut = hoe_cut & PFChIso_cut & PFHCalIso_cut

            # 4, 5
            sieie_cut = (photons.isScEtaEB & sieie_barrel_cut) | (photons.isScEtaEE & sieie_endcap_cut)
            PFECalIso_cut = (photons.isScEtaEB & PFECalIso_barrel_cut) | (photons.isScEtaEE & PFECalIso_endcap_cut)

            official_id_cut = hoe_cut & PFChIso_cut & PFHCalIso_cut & sieie_cut & PFECalIso_cut

        # Custom Photon ID
        id_cut = customized_id_cut
        
        # sieie
        # id_cut = sieie_cut

        # PFECalIso 
        # id_cut = PFECalIso_cut

        # Official Photon ID 
        # id_cut = official_id_cut

        # electron veto
        e_veto_cut = (photons.electronVeto > options["e_veto"])
        
        photon_ele_idx = ak.where(ak.num(photons.electronIdx, axis=1) == 0, ak.ones_like(photons.pt)*-1, photons.electronIdx)
        new_pho = ak.unflatten(ak.unflatten(ak.flatten(photon_ele_idx), [1]*ak.sum(ak.num(photon_ele_idx))), ak.num(photon_ele_idx, axis=1))
        new_ele = ak.broadcast_arrays(electrons.Idx[:,None], new_pho, depth_limit=2)[0]
        eg_overlap_cut = ~ak.where(
            ak.is_none(electrons.Idx),
            ak.broadcast_arrays(photons.electronIdx, False)[1],
            ak.flatten(ak.any(new_pho[:, :, None] == new_ele, axis=-2), axis=-1)
        ) # some events may have no electrons, so we need to replace None with False

        # use_central_nano = options["use_central_nano"] # indicates whether we are using central nanoAOD (with some branches that are necessary for full diphoton preselection missing) or custom nanoAOD (with these branches added)

        # all_cuts = pt_cut & eta_cut & id_cut & e_veto_cut & eg_overlap_cut
        # # all_cuts = pt_cut & eta_cut & e_veto_cut & eg_overlap_cut # bing for CR selection

        # self.register_cuts(
        #         names = ["no_cut", "pt", "eta", "id", "e_veto", "ele_pho_overlap", "all"], #"pt", "eta", "id", "e_veto", "ele_pho_overlap", "all"
        #         results = [no_cut, pt_cut, eta_cut, id_cut, e_veto_cut, eg_overlap_cut, all_cuts], #pt_cut, eta_cut, id_cut, e_veto_cut, eg_overlap_cut, all_cuts
        #         cut_type = "SelectedPhoton"
        # )

        cut_names = ["no_cut", "pt", "eta", "id", "e_veto", "ele_pho_overlap"]
        cut_results = [no_cut, pt_cut, eta_cut, id_cut, e_veto_cut, eg_overlap_cut]
        # Make all_cuts and perform N-1 cut
        all_cuts = photons.pt > 0
        for i, cut in enumerate(cut_results):
            all_cuts = (all_cuts) & cut
            if i == 0:
                # In an event, at least 2 photons pass selections
                cut_results[i] = (ak.sum(cut, axis=1) > 1)
            else:
                cut_results[i] = (ak.sum(cut, axis=1) > 1) & cut_results[i-1]

        # Print out cut flow results
        self.register_cuts(
                names = cut_names,
                results = cut_results,
                cut_type = "SelectedPhoton"
        )

        return all_cuts

    # def calculate_alp_photon_isolation(self, events, alp_cand, photons, delta_r_cone=0.3, delta_r_self_match=0.08):
    #     """
    #     Calculate PF photon isolation for the ALP candidate (ALPCand).
    #     Sums the pt of reconstructed photons within a cone of delta_r_cone around ALPCand,
    #     excluding the lead and sublead photons that form the ALP candidate.

    #     :param events: Input events array to add isolation field to
    #     :type events: ak.highlevel.Array
    #     :param alp_cand: ALP candidate with LeadPhoton, SubleadPhoton, and ALPCand fields
    #     :type alp_cand: ak.highlevel.Array
    #     :param photons: Reconstructed photons (events.Photon after selection)
    #     :type photons: ak.highlevel.Array
    #     :param delta_r_cone: Delta R cone size for isolation (default: 0.3)
    #     :type delta_r_cone: float
    #     :param delta_r_self_match: Delta R to identify lead/sublead photons (default: 0.08)
    #     :type delta_r_self_match: float
    #     :return: Events array with added isolation field
    #     :rtype: ak.highlevel.Array
    #     """

    #     # Initialize isolation array with zeros
    #     iso_pt = ak.full_like(alp_cand["ALPCand"].pt, 0.0, dtype=float)

    #     # Ensure arrays have Momentum4D behavior
    #     alp_cand["ALPCand"] = ak.with_name(alp_cand["ALPCand"], "Momentum4D")
    #     photons = ak.with_name(photons, "Momentum4D")
    #     lead_photon = ak.with_name(alp_cand.LeadPhoton, "Momentum4D")
    #     sublead_photon = ak.with_name(alp_cand.SubleadPhoton, "Momentum4D")

    #     print(len(photons), len(alp_cand))

    #     print(photons.type.show())
    #     print(alp_cand["ALPCand"].type.show())
    #     assert hasattr(photons, "delta_r")
    #     assert hasattr(alp_cand["ALPCand"], "delta_r")

    #     # Compute DeltaR between ALPCand and all photons
    #     alp_broadcasted = ak.broadcast_arrays(alp_cand["ALPCand"], photons)[0]
    #     dr = alp_broadcasted.delta_r(photons)

    #     # Exclude lead and sublead photons by checking DeltaR with them
    #     lead_photon_broadcasted = ak.broadcast_arrays(lead_photon, photons)[0]
    #     lead_photon_vec_broadcasted = ak.with_name(lead_photon_broadcasted, "Momentum4D")
    #     dr_lead = lead_photon_vec_broadcasted.delta_r(photons)

    #     sublead_photon_broadcasted = ak.broadcast_arrays(sublead_photon, photons)[0]
    #     sublead_photon_vec_broadcasted = ak.with_name(sublead_photon_broadcasted, "Momentum4D")
    #     dr_sublead = sublead_photon_vec_broadcasted.delta_r(photons)
    #     is_self_matched = (dr_lead < delta_r_self_match) | (dr_sublead < delta_r_self_match)

    #     # Isolation cone condition (ΔR ≤ 0.3) and exclude self-matched photons
    #     iso_mask = (dr <= delta_r_cone) & (~is_self_matched)

    #     # Sum photon pt within the isolation cone
    #     iso_pt = ak.sum(photons.pt[iso_mask], axis=-1)

    #     # Fill None values with 0.0 for events with no photons in the cone
    #     iso_pt = ak.fill_none(iso_pt, 0.0)

    #     # Add isolation field to events
    #     events = ak.with_field(events, iso_pt, "ALP_PhotonIso")

    #     return events

    def calculate_alp_photon_isolation(self, events, alp_cand, photons, delta_r_cone=0.3, delta_r_self_match=0.08):

        # logger.debug(f"len(photons) = {len(photons)}, len(alp_cand) = {len(alp_cand)}")
        # print(len(photons), len(alp_cand))
        # photons_P4 = to_momentum4d(photons)
        # print(f"photons_P4: {type(photons_P4)}, {ak.type(photons_P4)}")
        # ph0 = photons_P4[0][0]
        # print(ph0)  # 这应该是一个 Momentum4D object
        # print(hasattr(ph0, "deltaR"))  # ✅ True！
        # print(ph0.deltaR(ph0))         # ✅ 应该是 0.0
        # print(photons.fields)
        # print(f"photons: {type(photons)}, {ak.type(photons)}")
        # print(hasattr(photons, "deltaR"))
        # print(photons[0].deltaR(photons[0]))
        # print(alp_b.fields)
        # print(f"alp_b: {type(alp_b)}, {ak.type(alp_b)}")
        # print(f"pho_b: {type(pho_b)}, {ak.type(pho_b)}")

        # Ensure momentum behavior
        photons = ak.with_name(photons, "Momentum4D")
        lead = ak.with_name(alp_cand.LeadPhoton, "Momentum4D")
        sublead = ak.with_name(alp_cand.SubleadPhoton, "Momentum4D")
        alp_vec = ak.with_name(alp_cand.ALPCand, "Momentum4D")

        # Broadcast to photons
        alp_b, pho_b = ak.broadcast_arrays(alp_vec, photons)
        alp_b = ak.with_name(alp_b, "Momentum4D")
        pho_b = ak.with_name(pho_b, "Momentum4D")
        
        dr = alp_b.deltaR(pho_b)

        # Exclude self photons (lead/sublead)
        lead_b = ak.broadcast_arrays(lead, photons)[0]
        sublead_b = ak.broadcast_arrays(sublead, photons)[0]
        dr_lead = lead_b.deltaR(photons)
        dr_sublead = sublead_b.deltaR(photons)
        is_self = (dr_lead < delta_r_self_match) | (dr_sublead < delta_r_self_match)

        # Mask and sum
        iso_mask = (dr <= delta_r_cone) & (~is_self)
        iso_pt = ak.sum(photons.pt[iso_mask], axis=-1)
        iso_pt = ak.fill_none(iso_pt, 0.0)

        # Add to events
        return ak.with_field(events, iso_pt, "ALP_PhotonIso")

    def select_FSRphotons(self, FSRphotons, electrons, muons, photons, options):
        
        FSR_pt_cut = FSRphotons.pt > options["pt"]
        FSR_eta_cut = abs(FSRphotons.eta) < options["eta"]

        FSR_iso_cut = FSRphotons.relIso03 < options["iso"]
        FSR_dROverEt2_cut = FSRphotons.dROverEt2 < options["dROverEt2"]
        
        # Clean FSR photons by checking dR with electrons
        FSRphoton_clean_electrons = object_selections.delta_R(FSRphotons, electrons, 0.001)

        # Discard FSRphoton if dR(FSRphoton, lep) > 0.5
        dR0p5_FSRphoton_lep = object_selections.delta_R(FSRphotons, electrons, 0.5) & object_selections.delta_R(FSRphotons, muons, 0.5)
        FSRphoton_lep_indR0p5 = ~dR0p5_FSRphoton_lep

        photons_sorted = photons[ak.argsort(photons.pt, ascending=False)]

        lead_photons = photons_sorted[:, :1]

        FSRphoton_clean_photons = object_selections.delta_R(FSRphotons, lead_photons, 0.2)

        FSR_all_cuts = FSR_pt_cut & FSR_eta_cut & FSR_iso_cut & FSR_dROverEt2_cut & FSRphoton_clean_electrons & FSRphoton_lep_indR0p5 & FSRphoton_clean_photons

        return FSR_all_cuts

    def assign_fsr_photon(self, leptons, fsr_photons):
        
        # ---- 0. 若所有事件都没有 photon，直接返回 ----
        if ak.max(ak.num(fsr_photons)) == 0:
            return leptons                    # nothing to do

        # ---- 1. 确保每个 event 至少有 1 个 “占位 photon” ----
        #    （避免在完全空列表上索引越界）
        zero_ph = ak.zip({"pt":   0.0,
                        "eta":  0.0,
                        "phi":  0.0,
                        "mass": 0.0,
                        "dROverEt2": numpy.inf,   # 填一个 ∞，保证不会被选中
                        },
                        with_name="Momentum4D")

        fsr_padded = ak.fill_none(ak.pad_none(fsr_photons, 1, clip=False),
                                zero_ph)     # 长度 <1 时补 None→zero_ph

        # ---- 2. 列出 (event, lepton, photon) 组合 ----
        pairs = ak.cartesian({"lep": leptons, "ph": fsr_padded},
                            axis=1,
                            nested=True)      # shape: (ev, lep, pho)

        dR = pairs.lep.deltaR(pairs.ph)

        # ---- 3. 只允许 dR < 0.5 的配对 ----
        valid_mask       = dR < 0.5
        dROverEt2        = pairs.ph.dROverEt2
        dROverEt2_masked = ak.where(valid_mask, dROverEt2, numpy.inf)

        # ---- 4. 每个 (ev, lep) 取最小值所在的 photon 下标 ----
        best_idx  = ak.argmin(dROverEt2_masked, axis=2)    # (ev, lep)
        has_valid = ak.any(valid_mask,          axis=2)    # (ev, lep)

        # ---- 5. 把无合法 photon 的 lepton 标成 None（用 ak.mask！）----
        best_idx_opt = ak.mask(best_idx, has_valid)        # OptionArray[int64]

        # ---- 6. 取出最佳 photon；None → None ----
        best_ph = fsr_padded[best_idx_opt]                 # OptionArray[Momentum4D]

        # ---- 7. None → 零动量向量 ----
        best_ph = ak.fill_none(best_ph, zero_ph)

        # ==== 8. 四動量相加，但只更新 pt/eta/phi/mass，其餘欄位保留 ====
        vec_sum  = leptons + best_ph
        corrected = leptons
        for fld in ("pt", "eta", "phi", "mass"):
            corrected = ak.with_field(corrected, getattr(vec_sum, fld), fld)

        # 最後修正 mass（這時才是真的寫到 corrected 上）
        corrected = ak.with_field(
            corrected,
            ak.where((corrected.mass < 0) | ak.is_none(corrected.mass), 0.0, corrected.mass),
            "mass"
        )

        # print("\n✅ 以下是有成功匹配 FSR photon 的 lepton：")

        # for i in range(len(leptons)):  # 所有 event
        #     if leptons[i] is None:
        #         continue

        #     for j in range(len(leptons[i])):  # 每個 lepton
        #         if not has_valid[i][j]:
        #             continue  # 忽略沒匹配成功的

        #         l0 = leptons[i][j]
        #         l1 = corrected[i][j]
        #         ph = best_ph[i][j]
        #         print(f"\nEvent {i}, Lepton {j} matched:")
        #         print(f"  Before: pt={l0.pt:.2f}, eta={l0.eta:.2f}, phi={l0.phi:.2f}, mass={l0.mass:.2f}")
        #         print(f"  After : pt={l1.pt:.2f}, eta={l1.eta:.2f}, phi={l1.phi:.2f}, mass={l1.mass:.2f}")
        #         print(f"  FSR   : pt={ph.pt:.2f},  eta={ph.eta:.2f},  phi={ph.phi:.2f},  mass={ph.mass:.2f}")

        return corrected


# Below is an example of how the diphoton preselection could be performed with an explicit loop (C++ style) 
# that is compiled with numba for increased performance.

# NOTE: pre-compiled numba functions should be defined outside the class,
# as numba does not like when a class instance is passed to a function
@numba.njit
def produce_diphotons(photons, n_photons, lead_pt_cut, lead_pt_mgg_cut, sublead_pt_mgg_cut):
    n_events = len(photons)

    diphotons_offsets = numpy.zeros(n_events + 1, numpy.int64)
    diphotons_contents = []
    lead_photons_idx = []
    sublead_photons_idx = []

    # Loop through events and select diphoton pairs
    for i in range(n_events):
        n_diphotons_event = 0
        # Enumerate photon pairs
        if n_photons[i] >= 2: # only try if there are at least 2 photons
            sum_pt = 0
            for j in range(n_photons[i]):
                for k in range(j+1, n_photons[i]):
                    # Choose who is the leading photon
                    lead_idx = j if photons[i][j].pt > photons[i][k].pt else k
                    sublead_idx = k if photons[i][j].pt > photons[i][k].pt else j
                    lead_photon = photons[i][lead_idx]
                    sublead_photon = photons[i][sublead_idx]

                    # Lead photon must satisfy lead pt cut
                    if lead_photon.pt < lead_pt_cut:
                        continue

                    # Construct four vectors
                    lead_photon_vec = vector.obj(
                            pt = lead_photon.pt,
                            eta = lead_photon.eta,
                            phi = lead_photon.phi,
                            mass = lead_photon.mass
                    )
                    sublead_photon_vec = vector.obj(
                            pt = sublead_photon.pt,
                            eta = sublead_photon.eta,
                            phi = sublead_photon.phi,
                            mass = sublead_photon.mass
                    )
                    diphoton = vector.obj(px = 0., py = 0., pz = 0., E = 0.) # IMPORTANT NOTE: you need to initialize this to an empty vector first. Otherwise, you will get ZeroDivisionError exceptions for like 1 out of a million events (seemingly only with numba). 
                    diphoton = diphoton + lead_photon_vec + sublead_photon_vec

                    if (diphoton.mass < 100) | (diphoton.mass > 180):
                        continue

                    if lead_photon.pt / diphoton.mass < lead_pt_mgg_cut:
                        continue

                    if sublead_photon.pt / diphoton.mass < sublead_pt_mgg_cut:
                        continue

                    # This diphoton candidate passes
                    n_diphotons_event += 1

                    diphotons_contents.append([
                        diphoton.pt,
                        diphoton.eta,
                        diphoton.phi,
                        diphoton.mass
                    ])

                    lead_photons_idx.append(lead_idx)
                    sublead_photons_idx.append(sublead_idx)

        diphotons_offsets[i+1] = diphotons_offsets[i] + n_diphotons_event

    return diphotons_offsets, numpy.array(diphotons_contents), numpy.array(lead_photons_idx), numpy.array(sublead_photons_idx)
