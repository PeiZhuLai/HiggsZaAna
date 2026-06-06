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
from higgs_dna.tools.zmass_constraint_refit import refit_lepton_pts

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
        # # Loose
        "loose_hoe_barrel": 0.129991,
        "loose_hoe_endcap": 0.153428,
        "loose_PFChIso_barrel": 1.88518,
        "loose_PFChIso_endcap": 1.65396,
        "loose_PFHCalIso_barrel": [6.34397,   0.0100547, 5.78332e-05],
        "loose_PFHCalIso_endcap": [1.85881,   0.0116989, 7.47603e-05],
        "loose_sieie_barrel": 0.0114521,
        "loose_sieie_endcap": 0.0276744,
        "loose_PFECalIso_barrel": [0.703789, 0.000652035],
        "loose_PFECalIso_endcap": [6.61585, 0.000195486],
        # # Medium
        "medium_hoe_barrel": 0.0583054,
        "medium_hoe_endcap": 0.00518075,
        "medium_PFChIso_barrel": 0.939289,
        "medium_PFChIso_endcap": 0.970286,
        "medium_PFHCalIso_barrel": [2.18903,   0.0100547, 5.78332e-05],
        "medium_PFHCalIso_endcap": [0.0336699, 0.0116989, 7.47603e-05],
        "medium_sieie_barrel": 0.0100086,
        "medium_sieie_endcap": 0.0268736,
        "medium_PFECalIso_barrel": [0.227697, 0.000652035],
        "medium_PFECalIso_endcap": [1.124, 0.000195486],
        # # Tight
        "tight_hoe_barrel": 0.0417588,
        "tight_hoe_endcap": 0.00254267,
        "tight_PFChIso_barrel": 0.316306,
        "tight_PFChIso_endcap": 0.292664,
        "tight_PFHCalIso_barrel": [0.39057,   0.0100547, 5.78332e-05],
        "tight_PFHCalIso_endcap": [0.0292617, 0.0116989, 7.47603e-05],
        "tight_sieie_barrel": 0.00999299,
        "tight_sieie_endcap": 0.0268702,
        "tight_PFECalIso_barrel": [0.14189, 0.000652035],
        "tight_PFECalIso_endcap": [1.04269, 0.000195486],
        #------------------------------------------
        # Effective Area Correction for hoe, PFChargedHadronIso, PFHCalIso, PFECalIso (EA * rho)
        #------------------------------------------
        "hoe_EA_EB_1": [0.00198598,  -1.15014e-05],
        "hoe_EA_EB_2": [0.00249958 , -1.15014e-05],

        "hoe_EA_EE_1": [0.00302416, -1.51973e-05],
        "hoe_EA_EE_2": [0.00349404, -1.49651e-05],
        "hoe_EA_EE_3": [0.00334674, -1.47232e-05],
        "hoe_EA_EE_4": [0.00416787, -2.13958e-05],
        "hoe_EA_EE_5": [0.00458106, -2.80795e-05],
        #------------------------------------------
        "PFCHIso_EA_EB_1": [0.0342898, -0.000103508],
        "PFCHIso_EA_EB_2": [0.0281424, -3.1494e-05],
        "PFCHIso_EA_EE_1": [0.0288533, -6.66148e-05],
        "PFCHIso_EA_EE_2": [0.028789, -6.84993e-05],
        "PFCHIso_EA_EE_3": [0.0264064, -8.89189e-05],
        "PFCHIso_EA_EE_4": [0.025587, -5.90178e-05],
        "PFCHIso_EA_EE_5": [0.0224817, -4.22712e-05],
        #------------------------------------------
        "PFECalIso_EA_EB_1": [0.0866519, -0.0002296],
        "PFECalIso_EA_EB_2": [0.0730397, -0.0002134],
        "PFECalIso_EA_EE_1": [0.0542479, -0.00010934],
        "PFECalIso_EA_EE_2": [0.0486181, -6.2098e-05],
        "PFECalIso_EA_EE_3": [0.0412923, -9.63732e-06],
        "PFECalIso_EA_EE_4": [0.03555,   -5.79549e-05],
        "PFECalIso_EA_EE_5": [0.0360895, -7.28546e-06],
        #------------------------------------------
        "PFHCalIso_EA_EB_1": [0.17005,  -0.000835],
        "PFHCalIso_EA_EB_2": [0.208571, -0.000905],
        "PFHCalIso_EA_EE_1": [0.246494, -0.000722],
        "PFHCalIso_EA_EE_2": [0.306529, -0.000608],
        "PFHCalIso_EA_EE_3": [0.322673, -0.000750],
        "PFHCalIso_EA_EE_4": [0.315793, -0.000795],
        "PFHCalIso_EA_EE_5": [0.36531,  -0.000439],
        "e_veto" : 0.5
    },
    "zgammas" : {
        "relative_pt_gamma" : 15.0/110.,
        "mass_h" : [95., 180.],
        "mass_sum" : 185,
        "select_highest_pt_sum" : True,
        # 新增：ALP isolation 除錯開關與最多印出的事件數
        "debug_alp_iso": False,
        "debug_alp_iso_max_events": 10,
        # NEW: 啟用 15 種 photonID 情境 cutflow
        "study_photon_id_scenarios": True,
        # （可選）若只想跑子集就填 list；None=全跑
        # 例如: ["phid_custom_tight", "phid_official_medium", ...]
        "study_photon_id_scenarios_subset": None,

        # NEW: electron sip3d scenario study
        "study_ele_ip3d_scenarios": True,

        # Study replacing nominal muon isolation with pfRelIso04_all < 0.2
        "study_mu_iso04_replacement": True,
        "study_mu_iso04_cut": 0.2,

        # Study removing nominal muon sip3d cut
        "study_mu_no_sip3d_replacement": True,

        # NEW: lepton pT-binned trigger efficiency study (phid-like)
        "study_lep_trigger_eff_ptbins": True,
        # bins: [low0, low1, ..., last_low, +inf)
        "ele_trigger_eff_ptbins": list(range(8, 102, 2)),
        "mu_trigger_eff_ptbins":  list(range(8, 102, 2)),
    },
    "single_muon_trigger":{
        "2016":["HLT_IsoMu24", "HLT_IsoTkMu24"],
        "2017":["HLT_IsoMu27"],
        "2018":["HLT_IsoMu24"],
        "2022":["HLT_IsoMu24"],
        "2023":["HLT_IsoMu24"],
        "2024":["HLT_IsoMu24"]
    },
    "single_ele_trigger":{
        "2016":["HLT_Ele27_WPTight_Gsf"],
        "2017":["HLT_Ele32_WPTight_Gsf_L1DoubleEG"],
        "2018":["HLT_Ele32_WPTight_Gsf"],
        "2022":["HLT_Ele30_WPTight_Gsf"],
        "2023":["HLT_Ele30_WPTight_Gsf"],
        "2024":["HLT_Ele30_WPTight_Gsf"]
    },
    "double_muon_trigger":{
        "2016":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"],
        "2017":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
        "2018":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
        "2022":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
        "2023":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
        "2024":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"]
    },
    "double_ele_trigger":{
        "2016":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
        "2017":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2018":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2022":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2023":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
        "2024":["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]
    },
    "electrons" : {
        "pt" : 7.0,
        # NOTE: keep nominal behavior unchanged by default
        # "ip3d_max": 4.0,  # enable only if you want nominal electron selection to include it
    },
    "muons" : {
        "pt" : 5.0
    },
    "lead_ele_pt":{
        "2016": 30,
        "2017": 35,
        "2018": 35,
        "2022": 35,
        "2023": 35,
        "2024": 35
    },
    "lead_mu_pt":{
        "2016": 25,
        "2017": 28,
        "2018": 25,
        "2022": 25,
        "2023": 25,
        "2024": 25
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
        "2023postBPix": 0.2435,
        "2024": 0.2435
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
    out = ak.Array(out.layout) 
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

        # NEW: 存 rho 到 events（event-level），方便後續 write_events 輸出
        awkward_utils.add_field(events, "fixedGridRhoAll", ak.fill_none(rho, 0.0), overwrite=True)

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

    def _build_best_z_candidates(self, electrons, muons):
        ee_pairs = ak.combinations(electrons, 2, fields=["LeadLepton", "SubleadLepton"])
        ee_pairs = ee_pairs[(ee_pairs.LeadLepton.charge * ee_pairs.SubleadLepton.charge) == -1]

        mm_pairs = ak.combinations(muons, 2, fields=["LeadLepton", "SubleadLepton"])
        mm_pairs = mm_pairs[(mm_pairs.LeadLepton.charge * mm_pairs.SubleadLepton.charge) == -1]

        z_cands = ak.concatenate([ee_pairs, mm_pairs], axis=1)
        z_cands["ZCand"] = z_cands.LeadLepton + z_cands.SubleadLepton
        z_cands = z_cands[ak.argsort(abs(z_cands.ZCand.mass - 91.1876), axis=1)]

        best_z = ak.firsts(z_cands)
        z_ee_cut = ak.fill_none(best_z.LeadLepton.id == 11, False)
        z_mumu_cut = ak.fill_none(best_z.LeadLepton.id == 13, False)
        return z_cands, best_z, z_ee_cut, z_mumu_cut

    def _refit_z_candidate(self, z_cand):
        """對選中的 Z 候選做 Z mass-constraint 逐輕子 pT refit，回傳 dressed refit Z 四動量。

        z_cand 為 event-level（option type）的最佳 Z 候選；其 LeadLepton / SubleadLepton
        皆帶有 bare_* 與 fsr_* 欄位（見 assign_fsr_photon / _attach_refit_fields）：
        bare_* 為 dressed 前的裸輕子四動量，fsr_* 為這顆輕子被指派到的 FSR 光子
        （沒有時為 0）。eta/phi/mass 固定，僅浮動兩顆裸輕子的 pT；FSR 光子固定併入
        m_ll 約束，refit 後再加回對應輕子組成 dressed 的 refit Z。沒有合法 Z 的事件
        回傳 None（讓下游 fill_none 寫入 DUMMY_VALUE）。
        """
        lead = z_cand.LeadLepton
        sub = z_cand.SubleadLepton

        def _flat(arr, fill=0.0):
            return ak.to_numpy(ak.fill_none(arr, fill))

        pt1_refit, pt2_refit = refit_lepton_pts(
            _flat(lead.bare_pt), _flat(lead.bare_eta), _flat(lead.bare_phi),
            _flat(lead.bare_mass), _flat(lead.ptE_error),
            _flat(sub.bare_pt), _flat(sub.bare_eta), _flat(sub.bare_phi),
            _flat(sub.bare_mass), _flat(sub.ptE_error),
            _flat(lead.fsr_pt), _flat(lead.fsr_eta), _flat(lead.fsr_phi),
            _flat(sub.fsr_pt), _flat(sub.fsr_eta), _flat(sub.fsr_phi),
        )

        def _dressed(obj, pt_new):
            lep = ak.zip(
                {
                    "pt": ak.Array(pt_new),
                    "eta": ak.fill_none(obj.bare_eta, 0.0),
                    "phi": ak.fill_none(obj.bare_phi, 0.0),
                    "mass": ak.fill_none(obj.bare_mass, 0.0),
                },
                with_name="Momentum4D",
            )
            # FSR 光子（無質量）；pt==0 表示沒有 FSR → 貢獻零向量
            fsr_pt = ak.fill_none(obj.fsr_pt, 0.0)
            fsr = ak.zip(
                {
                    "pt": fsr_pt,
                    "eta": ak.fill_none(obj.fsr_eta, 0.0),
                    "phi": ak.fill_none(obj.fsr_phi, 0.0),
                    "mass": ak.zeros_like(fsr_pt),
                },
                with_name="Momentum4D",
            )
            return lep + fsr

        z_refit_p4 = _dressed(lead, pt1_refit) + _dressed(sub, pt2_refit)

        # 沒有合法 Z 候選的事件標成 None，讓下游 fill_none 寫入 DUMMY_VALUE
        valid = ~ak.is_none(z_cand.LeadLepton.bare_pt)
        z_refit_p4 = ak.mask(z_refit_p4, valid)
        # per-lepton 的 refit pT（裸輕子，未含 FSR），同樣對無效事件 mask
        lead_pt_refit = ak.mask(ak.Array(pt1_refit), valid)
        sub_pt_refit = ak.mask(ak.Array(pt2_refit), valid)
        return z_refit_p4, lead_pt_refit, sub_pt_refit

    def _build_trigger_pt_cuts(
        self,
        electrons,
        muons,
        single_ele_trigger_cut,
        double_ele_trigger_cut,
        single_mu_trigger_cut,
        double_mu_trigger_cut,
    ):
        lead_ele_pt = ak.fill_none(ak.pad_none(electrons.pt, 1, axis=1)[:, 0], 0.0)
        sublead_ele_pt = ak.fill_none(ak.pad_none(electrons.pt, 2, axis=1)[:, 1], 0.0)
        lead_mu_pt = ak.fill_none(ak.pad_none(muons.pt, 1, axis=1)[:, 0], 0.0)
        sublead_mu_pt = ak.fill_none(ak.pad_none(muons.pt, 2, axis=1)[:, 1], 0.0)

        if self.year is not None:
            year = self.year[:4]
            e_cut = lead_ele_pt > self.options["lead_ele_pt"][year]
            m_cut = lead_mu_pt > self.options["lead_mu_pt"][year]
        else:
            e_cut = lead_ele_pt > 25
            m_cut = lead_mu_pt > 20

        ee_cut = (lead_ele_pt > 25) & (sublead_ele_pt > 15)
        mm_cut = (lead_mu_pt > 20) & (sublead_mu_pt > 10)

        ele_trigger_pt_cut = (single_ele_trigger_cut & e_cut) | (double_ele_trigger_cut & ee_cut)
        mu_trigger_pt_cut = (single_mu_trigger_cut & m_cut) | (double_mu_trigger_cut & mm_cut)
        trigger_pt_cut = ele_trigger_pt_cut | mu_trigger_pt_cut

        return {
            "lead_ele_pt": lead_ele_pt,
            "sublead_ele_pt": sublead_ele_pt,
            "lead_mu_pt": lead_mu_pt,
            "sublead_mu_pt": sublead_mu_pt,
            "ele_trigger_pt_cut": ele_trigger_pt_cut,
            "mu_trigger_pt_cut": mu_trigger_pt_cut,
            "trigger_pt_cut": trigger_pt_cut,
        }

    def _register_sequential_event_cutflow(self, events, cut_type, steps, weighted=False):
        names = ["all"]
        results = [ak.num(events.Photon) >= 0]
        running = results[0]

        for name, mask in steps:
            running = running & ak.fill_none(mask, False)
            names.append(name)
            results.append(running)

        names.append("all cuts")
        results.append(running)

        self.register_event_cuts(
            names=names,
            results=results,
            events=events,
            cut_type=cut_type,
            weighted=weighted,
        )
        return running

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
            clean = {},
            name = "SelectedElectron",
            tagger = self,
            year = self.year[:4]
        )
        
        electrons = awkward_utils.add_field(
            events = events,
            name = "SelectedElectron",
            data = events.Electron[electron_cut]
        )

        # NEW: event-level flag: 是否「選到的 electrons」全都通過 sip3d<4（或至少無 electron 時視為 True）
        pass_sel_ele_sip3d = ak.fill_none(ak.all(electrons.sip3d < 4, axis=1), True)
        awkward_utils.add_field(events, "pass_sel_ele_sip3d", ak.fill_none(pass_sel_ele_sip3d, False), overwrite=True)

        # generate the index in the original array and add to electrons
        arr = ak.local_index(events.Electron["pt"], axis=1)[electron_cut]
        electron_idx = ak.mask(arr, ak.num(arr) > 0)
        awkward_utils.add_field(events = electrons, name = "Idx", data = electron_idx)

        # Muons
        muon_cut = lepton_selections.select_muons(
            muons = events.Muon,
            options = self.options["muons"],
            clean = {},
            name = "SelectedMuon",
            tagger = self
        )

        muons = awkward_utils.add_field(
            events = events,
            name = "SelectedMuon",
            data = events.Muon[muon_cut]
        )

        muons_iso04 = None
        if options.get("study_mu_iso04_replacement", False):
            # Study path only: keep nominal muon cut unchanged, rebuild an alternative
            # selection where pfRelIso03_all is replaced by pfRelIso04_all.
            muon_iso04_options = dict(self.options["muons"])
            muon_iso04_options["pfRelIso03_all"] = None
            muon_iso04_options["pfRelIso04_all"] = options.get("study_mu_iso04_cut", 0.2)

            muon_iso04_cut = lepton_selections.select_muons(
                muons = events.Muon,
                options = muon_iso04_options,
                clean = {},
                name = "SelectedMuonIso04",
                tagger = self
            )

            muons_iso04 = awkward_utils.add_field(
                events = events,
                name = "SelectedMuonIso04",
                data = events.Muon[muon_iso04_cut]
            )

        muons_nosip3d = None
        if options.get("study_mu_no_sip3d_replacement", False):
            # Study path only: keep nominal muon cut unchanged, rebuild an alternative
            # selection where the muon sip3d requirement is removed.
            muon_nosip3d_options = dict(self.options["muons"])
            muon_nosip3d_options["sip3d"] = None

            muon_nosip3d_cut = lepton_selections.select_muons(
                muons = events.Muon,
                options = muon_nosip3d_options,
                clean = {},
                name = "SelectedMuonNoSip3d",
                tagger = self
            )

            muons_nosip3d = awkward_utils.add_field(
                events = events,
                name = "SelectedMuonNoSip3d",
                data = events.Muon[muon_nosip3d_cut]
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
        photon_selection, photons_with_flags = self.select_photons(
                photons = events.Photon,
                options = self.options["photons"],
                electrons = electrons,
                rho = rho,
                year = self.year[:4]
        )

        # 重要：把 pass_* flags 存回 events（保留「所有 photon」的 flags，方便後續研究）
        awkward_utils.add_field(
            events=events,
            name="Photon",
            data=photons_with_flags,
            overwrite=True
        )

        # 後續真正用來組 ALP / 做分析的 photons，才套用 selection mask
        photons = events.Photon[photon_selection]

        # lepton-photon overlap removal 
        clean_photon_mask = ak.fill_none(object_selections.delta_R(photons, muons, 0.3), True) & ak.fill_none(object_selections.delta_R(photons, electrons, 0.3), True)
        photons = photons[clean_photon_mask]
        photons = ak.with_field(photons, ak.ones_like(photons.pt) * 0.0, "mass")

        # 加入在原始 events.Photon 中的索引，方便稍後做全體 photon 的排除
        pho_orig_idx = ak.local_index(events.Photon.pt, axis=1)
        photons = ak.with_field(photons, pho_orig_idx[photon_selection][clean_photon_mask], "origIndex")


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

        FSRphotons_iso04 = None
        if muons_iso04 is not None:
            # Keep the photon side nominal and only rebuild the muon-dependent study path.
            FSRphoton_selection_iso04 = self.select_FSRphotons(
                    FSRphotons = events.FsrPhoton,
                    electrons = electrons,
                    muons = muons_iso04,
                    photons = photons,
                    options = self.options["FSRphotons"]
            )
            FSRphotons_iso04 = events.FsrPhoton[FSRphoton_selection_iso04]
            FSRphotons_iso04 = ak.with_field(FSRphotons_iso04, ak.ones_like(FSRphotons_iso04.pt) * 0.0, "mass")

        FSRphotons_nosip3d = None
        if muons_nosip3d is not None:
            FSRphoton_selection_nosip3d = self.select_FSRphotons(
                    FSRphotons = events.FsrPhoton,
                    electrons = electrons,
                    muons = muons_nosip3d,
                    photons = photons,
                    options = self.options["FSRphotons"]
            )
            FSRphotons_nosip3d = events.FsrPhoton[FSRphoton_selection_nosip3d]
            FSRphotons_nosip3d = ak.with_field(FSRphotons_nosip3d, ak.ones_like(FSRphotons_nosip3d.pt) * 0.0, "mass")

        # Add object fields to events array
        object_fields = [
            (electrons, "electron"),
            (muons, "muon"),
        ]
        if muons_iso04 is not None:
            object_fields.append((muons_iso04, "muon_iso04"))
        if muons_nosip3d is not None:
            object_fields.append((muons_nosip3d, "muon_nosip3d"))

        for objects, name in object_fields:
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

        if muons_iso04 is not None:
            awkward_utils.add_field(events, "n_muons_iso04", ak.fill_none(ak.num(muons_iso04), 0), overwrite=True)
        if muons_nosip3d is not None:
            awkward_utils.add_field(events, "n_muons_nosip3d", ak.fill_none(ak.num(muons_nosip3d), 0), overwrite=True)

        n_leptons = n_electrons + n_muons
        # N_e_mu_cut = N_e_cut | N_mu_cut
        awkward_utils.add_field(events, "n_leptons", n_leptons, overwrite=True)

        n_photons = ak.num(photons)
        logger.debug(f"Number of photons(tagger): {n_photons[:10]}")
        awkward_utils.add_field(events, "n_photons", n_photons, overwrite=True)

        # PDG ID
        electrons = ak.with_field(electrons, ak.ones_like(electrons.pt) * 11, "id")
        # electrons = ak.with_field(electrons, ak.ones_like(electrons.pt) * 0.00051099895, "mass")
        muons = ak.with_field(muons, ak.ones_like(muons.pt) * 13, "id")
        if muons_iso04 is not None:
            muons_iso04 = ak.with_field(muons_iso04, ak.ones_like(muons_iso04.pt) * 13, "id")
        if muons_nosip3d is not None:
            muons_nosip3d = ak.with_field(muons_nosip3d, ak.ones_like(muons_nosip3d.pt) * 13, "id")

        # leptons ptE_error
        electrons = ak.with_field(electrons, electrons.energyErr, "ptE_error")
        muons = ak.with_field(muons, muons.ptErr, "ptE_error")
        if muons_iso04 is not None:
            muons_iso04 = ak.with_field(muons_iso04, muons_iso04.ptErr, "ptE_error")
        if muons_nosip3d is not None:
            muons_nosip3d = ak.with_field(muons_nosip3d, muons_nosip3d.ptErr, "ptE_error")

        # Sort objects by pt
        photons = photons[ak.argsort(photons.pt, ascending=False, axis=1)]
        electrons = electrons[ak.argsort(electrons.pt, ascending=False, axis=1)]
        muons = muons[ak.argsort(muons.pt, ascending=False, axis=1)]
        if muons_iso04 is not None:
            muons_iso04 = muons_iso04[ak.argsort(muons_iso04.pt, ascending=False, axis=1)]
        if muons_nosip3d is not None:
            muons_nosip3d = muons_nosip3d[ak.argsort(muons_nosip3d.pt, ascending=False, axis=1)]

        # self.select_fake_and_medium_photons(events=events, photons=photons)

        # Register as `vector.Momentum4D` objects so we can do four-vector operations with them
        electrons = ak.Array(electrons, with_name = "Momentum4D")
        muons = ak.Array(muons, with_name = "Momentum4D")
        if muons_iso04 is not None:
            muons_iso04 = ak.Array(muons_iso04, with_name = "Momentum4D")
        if muons_nosip3d is not None:
            muons_nosip3d = ak.Array(muons_nosip3d, with_name = "Momentum4D")
        photons = ak.Array(photons, with_name = "Momentum4D")
        FSRphotons = ak.Array(FSRphotons, with_name = "Momentum4D")
        if FSRphotons_iso04 is not None:
            FSRphotons_iso04 = ak.Array(FSRphotons_iso04, with_name = "Momentum4D")
        if FSRphotons_nosip3d is not None:
            FSRphotons_nosip3d = ak.Array(FSRphotons_nosip3d, with_name = "Momentum4D")
        FSRphotons = ak.with_field(FSRphotons, FSRphotons.costheta, "costheta")
        if FSRphotons_iso04 is not None:
            FSRphotons_iso04 = ak.with_field(FSRphotons_iso04, FSRphotons_iso04.costheta, "costheta")
        if FSRphotons_nosip3d is not None:
            FSRphotons_nosip3d = ak.with_field(FSRphotons_nosip3d, FSRphotons_nosip3d.costheta, "costheta")

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
        z_cands_noFSR, z_cand_noFSR, z_ee_cut_noFSR, z_mumu_cut_noFSR = self._build_best_z_candidates(
            electrons,
            muons,
        )

        # 修正 electrons 和 muons（包含 FSR）
        # electrons_withFSR = self.assign_fsr_photon(electrons, FSRphotons)
        # 電子不做 FSR recovery（FSR 已大致併入 supercluster），但仍附上 Z-refit
        # 需要的 bare 四動量與 FSR=0 欄位，使下游 z_cand 的兩種 channel 欄位一致。
        _ele_zero = ak.zeros_like(electrons.pt)
        electrons_withFSR = self._attach_refit_fields(
            electrons, electrons, _ele_zero, _ele_zero, _ele_zero
        )
        muons_withFSR = self.assign_fsr_photon(muons, FSRphotons)
        z_cands, _, z_ee_cut, z_mumu_cut = self._build_best_z_candidates(
            electrons_withFSR,
            muons_withFSR,
        )

        z_cands_iso04 = None
        z_mumu_cut_iso04 = ak.num(events.Photon) < 0
        if muons_iso04 is not None and FSRphotons_iso04 is not None:
            muons_iso04_withFSR = self.assign_fsr_photon(muons_iso04, FSRphotons_iso04)
            z_cands_iso04, _, _, z_mumu_cut_iso04 = self._build_best_z_candidates(
                electrons_withFSR,
                muons_iso04_withFSR,
            )

        z_cands_nosip3d = None
        z_mumu_cut_nosip3d = ak.num(events.Photon) < 0
        if muons_nosip3d is not None and FSRphotons_nosip3d is not None:
            muons_nosip3d_withFSR = self.assign_fsr_photon(muons_nosip3d, FSRphotons_nosip3d)
            z_cands_nosip3d, _, _, z_mumu_cut_nosip3d = self._build_best_z_candidates(
                electrons_withFSR,
                muons_nosip3d_withFSR,
            )

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

        pass_only_dimu = ak.fill_none(double_mu_trigger_cut, False) & ~ak.fill_none(single_mu_trigger_cut, False)
        pass_only_diel = ak.fill_none(double_ele_trigger_cut, False) & ~ak.fill_none(single_ele_trigger_cut, False)
        awkward_utils.add_field(events, "pass_only_dimu", pass_only_dimu, overwrite=True)
        awkward_utils.add_field(events, "pass_simu_or_dimu", ak.fill_none(mu_trigger_cut, False), overwrite=True)
        awkward_utils.add_field(events, "pass_only_diel", pass_only_diel, overwrite=True)
        awkward_utils.add_field(events, "pass_siel_or_diel", ak.fill_none(ele_trigger_cut, False), overwrite=True)

        trigger_pt_info = self._build_trigger_pt_cuts(
            electrons,
            muons,
            single_ele_trigger_cut,
            double_ele_trigger_cut,
            single_mu_trigger_cut,
            double_mu_trigger_cut,
        )
        ele_trigger_pt_cut = trigger_pt_info["ele_trigger_pt_cut"]
        mu_trigger_pt_cut = trigger_pt_info["mu_trigger_pt_cut"]
        trigger_pt_cut = trigger_pt_info["trigger_pt_cut"]

        mu_trigger_pt_cut_iso04 = ak.num(events.Photon) < 0
        if muons_iso04 is not None:
            trigger_pt_info_iso04 = self._build_trigger_pt_cuts(
                electrons,
                muons_iso04,
                single_ele_trigger_cut,
                double_ele_trigger_cut,
                single_mu_trigger_cut,
                double_mu_trigger_cut,
            )
            mu_trigger_pt_cut_iso04 = trigger_pt_info_iso04["mu_trigger_pt_cut"]

        mu_trigger_pt_cut_nosip3d = ak.num(events.Photon) < 0
        if muons_nosip3d is not None:
            trigger_pt_info_nosip3d = self._build_trigger_pt_cuts(
                electrons,
                muons_nosip3d,
                single_ele_trigger_cut,
                double_ele_trigger_cut,
                single_mu_trigger_cut,
                double_mu_trigger_cut,
            )
            mu_trigger_pt_cut_nosip3d = trigger_pt_info_nosip3d["mu_trigger_pt_cut"]

        # NEW: lepton pT-binned trigger efficiency study (phid-like cutflows)
        if self.options.get("zgammas", {}).get("study_lep_trigger_eff_ptbins", False):
            # 4 種 trigger 定義（你指定要研究的）
            trig_defs = {
                "double_ele": double_ele_trigger_cut,
                "double_mu":  double_mu_trigger_cut,
                "OR_ele":     (single_ele_trigger_cut | double_ele_trigger_cut),
                "OR_mu":      (single_mu_trigger_cut  | double_mu_trigger_cut),
            }

            def _register_ptbins(lep: str, ptvals, bins, trig_key: str, trig_mask, ord_label: str):
                # bins: list of lows; auto add +inf as last high
                lows = list(bins)
                highs = lows[1:] + [numpy.inf]

                for lo, hi in zip(lows, highs):
                    in_bin = (ptvals >= lo) & (ptvals < hi)
                    passed = in_bin & ak.fill_none(trig_mask, False)

                    bin_label = f"pt{int(lo)}to{('Inf' if hi == numpy.inf else int(hi))}"
                    cut_type = f"trigeff_{lep}_{ord_label}_{trig_key}_{bin_label}"

                    self.register_event_cuts(
                        names=["in_bin", "pass_trigger"],
                        results=[ak.fill_none(in_bin, False), ak.fill_none(passed, False)],
                        events=events,
                        cut_type=cut_type,
                        weighted=False,
                    )

            ele_bins = self.options.get("zgammas", {}).get("ele_trigger_eff_ptbins", list(range(8, 102, 2)))
            mu_bins  = self.options.get("zgammas", {}).get("mu_trigger_eff_ptbins",  list(range(8, 102, 2)))

            # electron: lead/sublead 各研究 double_ele 與 OR_ele
            _register_ptbins("ele", trigger_pt_info["lead_ele_pt"],    ele_bins, "double_ele", trig_defs["double_ele"], "lead")
            _register_ptbins("ele", trigger_pt_info["lead_ele_pt"],    ele_bins, "OR_ele",     trig_defs["OR_ele"],     "lead")
            _register_ptbins("ele", trigger_pt_info["sublead_ele_pt"], ele_bins, "double_ele", trig_defs["double_ele"], "sublead")
            _register_ptbins("ele", trigger_pt_info["sublead_ele_pt"], ele_bins, "OR_ele",     trig_defs["OR_ele"],     "sublead")

            # muon: lead/sublead 各研究 double_mu 與 OR_mu
            _register_ptbins("mu", trigger_pt_info["lead_mu_pt"],    mu_bins, "double_mu", trig_defs["double_mu"], "lead")
            _register_ptbins("mu", trigger_pt_info["lead_mu_pt"],    mu_bins, "OR_mu",     trig_defs["OR_mu"],     "lead")
            _register_ptbins("mu", trigger_pt_info["sublead_mu_pt"], mu_bins, "double_mu", trig_defs["double_mu"], "sublead")
            _register_ptbins("mu", trigger_pt_info["sublead_mu_pt"], mu_bins, "OR_mu",     trig_defs["OR_mu"],     "sublead")

        z_mumu = z_mumu_cut 
        z_ee = z_ee_cut
        awkward_utils.add_field(events, "z_mumu", z_mumu, overwrite=True)
        awkward_utils.add_field(events, "z_ee", z_ee, overwrite=True)
        if muons_iso04 is not None:
            awkward_utils.add_field(events, "z_mumu_iso04", z_mumu_cut_iso04, overwrite=True)
        if muons_nosip3d is not None:
            awkward_utils.add_field(events, "z_mumu_nosip3d", z_mumu_cut_nosip3d, overwrite=True)

        # mass_cut = (z_cands.ZCand.mass > 80.) & (z_cands.ZCand.mass < 100.)
        mass_cut = z_cands.ZCand.mass > 50.
        z_cands = z_cands[mass_cut] # OSSF lepton pairs with m_ll > 50.
        z_cands_noFSR = z_cands_noFSR[mass_cut]

        has_z_cand = ak.num(z_cands) >= 1
        # awkward_utils.add_field(events, "pass_1_Zcand", ak.fill_none(has_z_cand, False), overwrite=True)

        z_cand = ak.firsts(z_cands)
        z_cand_noFSR = ak.firsts(z_cands_noFSR)

        has_z_cand_iso04 = ak.num(events.Photon) < 0
        z_cand_iso04 = None
        if z_cands_iso04 is not None:
            mass_cut_iso04 = z_cands_iso04.ZCand.mass > 50.
            z_cands_iso04 = z_cands_iso04[mass_cut_iso04]
            has_z_cand_iso04 = ak.num(z_cands_iso04) >= 1
            z_cand_iso04 = ak.firsts(z_cands_iso04)

        has_z_cand_nosip3d = ak.num(events.Photon) < 0
        z_cand_nosip3d = None
        if z_cands_nosip3d is not None:
            mass_cut_nosip3d = z_cands_nosip3d.ZCand.mass > 50.
            z_cands_nosip3d = z_cands_nosip3d[mass_cut_nosip3d]
            has_z_cand_nosip3d = ak.num(z_cands_nosip3d) >= 1
            z_cand_nosip3d = ak.firsts(z_cands_nosip3d)

        # --- Z mass-constraint kinematic refit (per-lepton pT likelihood fit) ---
        # 對選中的 Z 候選，固定 eta/phi/mass，僅浮動兩顆「裸輕子」的 pT 以最大化
        # Gauss(pt1)*Gauss(pt2)*ZLineshape(m_ll) 的似然，藉此提升 m_ll → m_Za 的
        # 解析度（見 higgs_dna/tools/zmass_constraint_refit.py）。muon channel 的
        # FSR 光子會被併入 m_ll 約束（固定在 reco 值），refit 後再加回對應輕子上以
        # 組成 dressed 的 refit Z；電子不做 FSR（已大致併入 supercluster）。
        z_refit_p4, lead_pt_refit, sub_pt_refit = self._refit_z_candidate(z_cand)

        # Add Z-related fields to array
        # NOTE: sip3d only (no ip3d fallback)
        for field in ["pt", "eta", "costheta", "phi", "mass", "charge", "id", "ptE_error", "sip3d"]:
            if field == "sip3d":
                awkward_utils.add_field(events, "Z_lead_lepton_sip3d", ak.fill_none(z_cand.LeadLepton.sip3d, DUMMY_VALUE))
                awkward_utils.add_field(events, "Z_sublead_lepton_sip3d", ak.fill_none(z_cand.SubleadLepton.sip3d, DUMMY_VALUE))
                awkward_utils.add_field(events, "Z_noFSR_lead_lepton_sip3d", ak.fill_none(z_cand_noFSR.LeadLepton.sip3d, DUMMY_VALUE))
                awkward_utils.add_field(events, "Z_noFSR_sublead_lepton_sip3d", ak.fill_none(z_cand_noFSR.SubleadLepton.sip3d, DUMMY_VALUE))
                continue

            if not field in ["charge", "id", "ptE_error", "sip3d"]:
                awkward_utils.add_field(
                        events,
                        "Z_%s" % field,
                        ak.fill_none(getattr(z_cand.ZCand, field), DUMMY_VALUE)
                )
                awkward_utils.add_field(
                    events,
                    f"Z_noFSR_{field}",
                    ak.fill_none(getattr(z_cand_noFSR.ZCand, field), DUMMY_VALUE)
                )

            awkward_utils.add_field(
                    events,
                    "Z_lead_lepton_%s" % field,
                    ak.fill_none(getattr(z_cand.LeadLepton, field), DUMMY_VALUE)
            )
            awkward_utils.add_field(
                    events,
                    "Z_sublead_lepton_%s" % field,
                    ak.fill_none(getattr(z_cand.SubleadLepton, field), DUMMY_VALUE)
            )
            awkward_utils.add_field(
                events,
                f"Z_noFSR_lead_lepton_{field}",
                ak.fill_none(getattr(z_cand_noFSR.LeadLepton, field), DUMMY_VALUE)
            )
            awkward_utils.add_field(
                events,
                f"Z_noFSR_sublead_lepton_{field}",
                ak.fill_none(getattr(z_cand_noFSR.SubleadLepton, field), DUMMY_VALUE)
            )

        # Z mass-constraint refit 輸出（leptons + muon FSR；ALP/photon 不變）
        for field in ["pt", "eta", "phi", "mass"]:
            awkward_utils.add_field(
                events,
                "Z_%s_refit" % field,
                ak.fill_none(getattr(z_refit_p4, field), DUMMY_VALUE),
                overwrite=True,
            )
        # per-lepton 的 refit pT（裸輕子，未含 FSR）
        awkward_utils.add_field(
            events, "Z_lead_lepton_pt_refit",
            ak.fill_none(lead_pt_refit, DUMMY_VALUE), overwrite=True,
        )
        awkward_utils.add_field(
            events, "Z_sublead_lepton_pt_refit",
            ak.fill_none(sub_pt_refit, DUMMY_VALUE), overwrite=True,
        )
        # m_Z refit 別名
        awkward_utils.add_field(
            events, "mZ_refit",
            ak.fill_none(z_refit_p4.mass, DUMMY_VALUE), overwrite=True,
        )

        # Make gamma candidate-level cuts
        has_2gamma_cand = (ak.num(photons) >= 2) #& (events.n_iso_photons == 0) # only for dy samples
        # awkward_utils.add_field(events, "pass_1_ALP", ak.fill_none(has_2gamma_cand, False), overwrite=True)

        # ------------------------------------------------------------
        # 重新以 index 方式構建 gamma pairs（保留原本邏輯但加入索引）
        # ------------------------------------------------------------
        idx_photons = ak.local_index(photons.pt, axis=1)
        pairs_idx = ak.combinations(idx_photons, 2, fields=["i1", "i2"])
        lead_idx = ak.where(photons.pt[pairs_idx.i1] >= photons.pt[pairs_idx.i2], pairs_idx.i1, pairs_idx.i2)
        sub_idx  = ak.where(photons.pt[pairs_idx.i1] >= photons.pt[pairs_idx.i2], pairs_idx.i2, pairs_idx.i1)

        lead_ph = photons[lead_idx]
        sub_ph  = photons[sub_idx]

        # 同時保存原始 Photon 索引，後續隔離時用來排除 ALP 的兩顆 photon
        lead_orig = photons.origIndex[lead_idx]
        sub_orig  = photons.origIndex[sub_idx]

        gamma_pairs = ak.zip({
            "LeadPhoton": lead_ph,
            "SubleadPhoton": sub_ph,
            "LeadIndex": lead_idx,
            "SubleadIndex": sub_idx,
            "LeadOrigIndex": lead_orig,
            "SubleadOrigIndex": sub_orig
        })

        # 依照 leadPhoton pt 排序並取第一個（最高 lead pt）
        gamma_pairs = gamma_pairs[ak.argsort(gamma_pairs.LeadPhoton.pt, ascending=False, axis=1)]
        alp_cand = ak.firsts(gamma_pairs)

        # 建立 ALP 四動量
        alp_cand = ak.with_field(alp_cand,
                                 ak.with_name(alp_cand.LeadPhoton + alp_cand.SubleadPhoton, "Momentum4D"),
                                 "ALPCand")
        # ------------------------------------------------------------
        # 以 ALP 四動量對「所有未經過篩選的 photon」做 ΔR<0.3 的 pt 加總，
        # 並排除組成 ALP 的那兩顆 photon（用原始索引判斷）
        # ------------------------------------------------------------
        alp_iso = self._compute_alp_iso(events.Photon, alp_cand, dr_cone=0.3)
        events = ak.with_field(events, alp_iso, "ALP_PhotonIso")

        # 只打印非 0 的 ALP_PhotonIso（除錯）
        if self.options.get("zgammas", {}).get("debug_alp_iso", False):
            try:
                iso_np = ak.to_numpy(alp_iso)
                nz_idx = numpy.flatnonzero(iso_np > 0)
                n_all = iso_np.shape[0]
                n_nz = nz_idx.shape[0]
                logger.info(f"[ALP ISO] nonzero_count={n_nz}/{n_all}")
                K = int(min(self.options["zgammas"].get("debug_alp_iso_max_events", 20), n_nz))
                for j in range(K):
                    i = int(nz_idx[j])
                    run = getattr(events, "run", None)
                    ls  = getattr(events, "luminosityBlock", None)
                    evt = getattr(events, "event", None)
                    run_i = int(run[i]) if run is not None else -1
                    ls_i  = int(ls[i])  if ls  is not None else -1
                    evt_i = int(evt[i]) if evt is not None else -1
                    logger.info(f"[ALP ISO][nonzero #{j}] idx={i} run/lumi/evt={run_i}/{ls_i}/{evt_i} iso={float(iso_np[i]):.3f}")
            except Exception as e:
                logger.debug(f"[ALP ISO] nonzero print failed: {e}")
        # ------------------------------------------------------------

        # Add ALP-related fields
        for field in ["pt", "eta", "costheta", "phi", "mass", "electronVeto", "energyRaw", "energyErr", "r9", "sieie", "hoe", "hoe_PUcorr", "hcalPFClusterIso", "hcalPFClusterIso_PUcorr", "ecalPFClusterIso", "ecalPFClusterIso_PUcorr", "chiso_PUcorr", "sieip", "etaWidth", "phiWidth", "s4", "trkSumPtHollowConeDR03", "trkSumPtSolidConeDR04", "pfChargedIso", "pfChargedIsoWorstVtx", "esEffSigmaRR", "esEnergyOverRawE", "mvaID"]:
            if field in ["pt","eta","costheta","phi","mass"]:
                awkward_utils.add_field(
                    events,
                    "ALP_%s" % field,
                    ak.fill_none(getattr(alp_cand.ALPCand, field), DUMMY_VALUE)
                )
            name = f"ALP_lead_photon_{field}"
            if name not in events.fields:
                awkward_utils.add_field(events, name, ak.fill_none(getattr(alp_cand.LeadPhoton, field), DUMMY_VALUE))

            name = f"ALP_sublead_photon_{field}"
            if name not in events.fields:
                awkward_utils.add_field(events, name, ak.fill_none(getattr(alp_cand.SubleadPhoton, field), DUMMY_VALUE))


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
        
        # NEW: ALP lead/sublead photon 與 Z lead/sublead lepton 的最近 dR（取兩個 lepton 中較小者）
        # 用 event-level branch 組 Momentum4D，避免處理複雜的 jagged/索引情況
        alp_lead_pho_p4 = ak.zip(
            {
                "pt":   ak.fill_none(events.ALP_lead_photon_pt, 0.0),
                "eta":  ak.fill_none(events.ALP_lead_photon_eta, 0.0),
                "phi":  ak.fill_none(events.ALP_lead_photon_phi, 0.0),
                "mass": ak.zeros_like(ak.fill_none(events.ALP_lead_photon_pt, 0.0)),
            },
            with_name="Momentum4D",
        )
        alp_sub_pho_p4 = ak.zip(
            {
                "pt":   ak.fill_none(events.ALP_sublead_photon_pt, 0.0),
                "eta":  ak.fill_none(events.ALP_sublead_photon_eta, 0.0),
                "phi":  ak.fill_none(events.ALP_sublead_photon_phi, 0.0),
                "mass": ak.zeros_like(ak.fill_none(events.ALP_sublead_photon_pt, 0.0)),
            },
            with_name="Momentum4D",
        )
        z_lead_lep_p4 = ak.zip(
            {
                "pt":   ak.fill_none(events.Z_lead_lepton_pt, 0.0),
                "eta":  ak.fill_none(events.Z_lead_lepton_eta, 0.0),
                "phi":  ak.fill_none(events.Z_lead_lepton_phi, 0.0),
                "mass": ak.zeros_like(ak.fill_none(events.Z_lead_lepton_pt, 0.0)),
            },
            with_name="Momentum4D",
        )
        z_sub_lep_p4 = ak.zip(
            {
                "pt":   ak.fill_none(events.Z_sublead_lepton_pt, 0.0),
                "eta":  ak.fill_none(events.Z_sublead_lepton_eta, 0.0),
                "phi":  ak.fill_none(events.Z_sublead_lepton_phi, 0.0),
                "mass": ak.zeros_like(ak.fill_none(events.Z_sublead_lepton_pt, 0.0)),
            },
            with_name="Momentum4D",
        )

        # 有效性：需要 Z cand 與 ALP cand (用 pt != DUMMY_VALUE 當 guard)
        has_valid_z   = (ak.fill_none(events.Z_pt, DUMMY_VALUE) != DUMMY_VALUE)
        has_valid_alp = (ak.fill_none(events.ALP_pt, DUMMY_VALUE) != DUMMY_VALUE)
        valid_dr = ak.fill_none(has_valid_z & has_valid_alp, False)

        lead_dr_min = ak.where(
            valid_dr,
            ak.min(
                ak.concatenate(
                    [
                        alp_lead_pho_p4.deltaR(z_lead_lep_p4)[:, None],
                        alp_lead_pho_p4.deltaR(z_sub_lep_p4)[:, None],
                    ],
                    axis=1,
                ),
                axis=1,
            ),
            DUMMY_VALUE,
        )
        sub_dr_min = ak.where(
            valid_dr,
            ak.min(
                ak.concatenate(
                    [
                        alp_sub_pho_p4.deltaR(z_lead_lep_p4)[:, None],
                        alp_sub_pho_p4.deltaR(z_sub_lep_p4)[:, None],
                    ],
                    axis=1,
                ),
                axis=1,
            ),
            DUMMY_VALUE,
        )

        awkward_utils.add_field(events, "ALP_lead_photon_lep_near_dR", ak.fill_none(lead_dr_min, DUMMY_VALUE), overwrite=True)
        awkward_utils.add_field(events, "ALP_sublead_photon_lep_near_dR", ak.fill_none(sub_dr_min, DUMMY_VALUE), overwrite=True)

        # Gamma candidate 
        gamma_cand = ak.firsts(photons)
        gamma_mvaID_WPL = ((gamma_cand.isScEtaEB & (gamma_cand.mvaID > self.options["photons"]["mvaID_barrel"])) | (gamma_cand.isScEtaEE & (gamma_cand.mvaID > self.options["photons"]["mvaID_endcap"])))
        gamma_e_veto = gamma_cand.electronVeto > self.options["photons"]["e_veto"]

        # awkward_utils.add_field(gamma_cand, "mass", ak.ones_like(gamma_cand.pt) * 0) #TODO: run3 BUG
        if "mass" not in gamma_cand.fields:
            awkward_utils.add_field(gamma_cand, "mass", ak.zeros_like(gamma_cand.pt))

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
        sel_h = (h_cand.mass > options["mass_h"][0]) & (h_cand.mass < options["mass_h"][1])
        sel_h = ak.fill_none(sel_h, value = False)
        
        # Use No FSR Z boson 
        h_cand_noFSR = z_cand_noFSR.ZCand + alp_cand.ALPCand
        sel_h_noFSR = (h_cand_noFSR.mass > options["mass_h"][0]) & (h_cand_noFSR.mass < options["mass_h"][1])
        sel_h_noFSR = ak.fill_none(sel_h_noFSR, value=False)

        sel_h_iso04 = ak.num(events.Photon) < 0
        if z_cand_iso04 is not None:
            h_cand_iso04 = z_cand_iso04.ZCand + alp_cand.ALPCand
            sel_h_iso04 = (h_cand_iso04.mass > options["mass_h"][0]) & (h_cand_iso04.mass < options["mass_h"][1])
            sel_h_iso04 = ak.fill_none(sel_h_iso04, value=False)

        sel_h_nosip3d = ak.num(events.Photon) < 0
        if z_cand_nosip3d is not None:
            h_cand_nosip3d = z_cand_nosip3d.ZCand + alp_cand.ALPCand
            sel_h_nosip3d = (h_cand_nosip3d.mass > options["mass_h"][0]) & (h_cand_nosip3d.mass < options["mass_h"][1])
            sel_h_nosip3d = ak.fill_none(sel_h_nosip3d, value=False)

        logger.debug(
            f'!!!!has H: {sum(has_z_cand & has_2gamma_cand & z_mumu_cut)} | '
            f'{sum(has_z_cand & has_2gamma_cand & z_ee_cut)}'
        )

        # Add Higgs-related fields to array
        for field in ["pt", "eta", "costheta", "phi", "mass"]:
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

        # Z mass-constraint refit 後的 Higgs(=Za) 候選：refit 後的 Z 加上不變的 ALP
        h_cand_refit = z_refit_p4 + alp_cand.ALPCand
        for field in ["pt", "eta", "phi", "mass"]:
            awkward_utils.add_field(
                events,
                "H_%s_refit" % field,
                ak.fill_none(getattr(h_cand_refit, field), DUMMY_VALUE),
                overwrite=True,
            )
        # m_H refit 別名
        awkward_utils.add_field(
            events, "mH_refit",
            ak.fill_none(h_cand_refit.mass, DUMMY_VALUE), overwrite=True,
        )

        all_cuts = trigger_pt_cut & has_z_cand & has_2gamma_cand & sel_h #& ak.fill_none((h_cand.mass>80) & (h_cand.mass < options["mass_h"][1]), False)

        ee_all_cut = ak.num(events.Photon) < 0
        mm_all_cut = ak.num(events.Photon) < 0
        nominal_cutflows = {
            "zgammas": (
                [("N_lep_sel", z_ee_cut | z_mumu_cut), ("trig_cut", trigger_cut), ("lep_pt_cut", trigger_pt_cut)],
                False,
            ),
            "zgammas_ele": (
                [("N_lep_sel", z_ee_cut), ("trig_cut", ele_trigger_cut), ("lep_pt_cut", ele_trigger_pt_cut)],
                False,
            ),
            "zgammas_mu": (
                [("N_lep_sel", z_mumu_cut), ("trig_cut", mu_trigger_cut), ("lep_pt_cut", mu_trigger_pt_cut)],
                False,
            ),
            "zgammas_w": (
                [("N_lep_sel", z_ee_cut | z_mumu_cut), ("trig_cut", trigger_cut), ("lep_pt_cut", trigger_pt_cut)],
                True,
            ),
            "zgammas_ele_w": (
                [("N_lep_sel", z_ee_cut), ("trig_cut", ele_trigger_cut), ("lep_pt_cut", ele_trigger_pt_cut)],
                True,
            ),
            "zgammas_mu_w": (
                [("N_lep_sel", z_mumu_cut), ("trig_cut", mu_trigger_cut), ("lep_pt_cut", mu_trigger_pt_cut)],
                True,
            ),
        }
        trailing_nominal_steps = [
            ("has_z_cand", has_z_cand),
            ("has_2g_cand", has_2gamma_cand),
            ("sel_h", sel_h),
        ]

        for cut_type, (leading_steps, weighted) in nominal_cutflows.items():
            if weighted and not hasattr(events, "Generator_weight"):
                continue

            final_cut = self._register_sequential_event_cutflow(
                events=events,
                cut_type=cut_type,
                steps=leading_steps + trailing_nominal_steps,
                weighted=weighted,
            )

            if cut_type == "zgammas_ele":
                ee_all_cut = final_cut
            elif cut_type == "zgammas_mu":
                mm_all_cut = final_cut

        # electron sip3d scenario cutflow
        if self.options.get("zgammas", {}).get("study_ele_ip3d_scenarios", False):
            final_ele_sip3d_cut = ak.num(events.Photon) < 0
            ele_sip3d_steps = [
                ("N_lep_sel", z_ee_cut),
                ("ele_sip3d_cut", events.pass_sel_ele_sip3d),
                ("trig_cut", ele_trigger_cut),
                ("lep_pt_cut", ele_trigger_pt_cut),
                ("has_z_cand", has_z_cand),
                ("has_2g_cand", has_2gamma_cand),
                ("sel_h", sel_h),
            ]

            for weighted in (False, True):
                if weighted and not hasattr(events, "Generator_weight"):
                    continue

                final_cut = self._register_sequential_event_cutflow(
                    events=events,
                    cut_type="zgammas_ele_elesip3d_w" if weighted else "zgammas_ele_elesip3d",
                    steps=ele_sip3d_steps,
                    weighted=weighted,
                )
                if not weighted:
                    final_ele_sip3d_cut = final_cut

            awkward_utils.add_field(events, "pass_allcuts_elesip3d", ak.fill_none(final_ele_sip3d_cut, False), overwrite=True)

        if muons_iso04 is not None:
            final_mu_iso04_cut = ak.num(events.Photon) < 0
            mu_iso04_steps = [
                ("N_lep_sel", z_mumu_cut_iso04),
                ("trig_cut", mu_trigger_cut),
                ("lep_pt_cut", mu_trigger_pt_cut_iso04),
                ("has_z_cand", has_z_cand_iso04),
                ("has_2g_cand", has_2gamma_cand),
                ("sel_h", sel_h_iso04),
            ]

            for weighted in (False, True):
                if weighted and not hasattr(events, "Generator_weight"):
                    continue

                final_cut = self._register_sequential_event_cutflow(
                    events=events,
                    cut_type="zgammas_mu_muiso04_w" if weighted else "zgammas_mu_muiso04",
                    steps=mu_iso04_steps,
                    weighted=weighted,
                )
                if not weighted:
                    final_mu_iso04_cut = final_cut

            awkward_utils.add_field(events, "pass_allcuts_muiso04", ak.fill_none(final_mu_iso04_cut, False), overwrite=True)

        if muons_nosip3d is not None:
            final_mu_nosip3d_cut = ak.num(events.Photon) < 0
            mu_nosip3d_steps = [
                ("N_lep_sel", z_mumu_cut_nosip3d),
                ("trig_cut", mu_trigger_cut),
                ("lep_pt_cut", mu_trigger_pt_cut_nosip3d),
                ("has_z_cand", has_z_cand_nosip3d),
                ("has_2g_cand", has_2gamma_cand),
                ("sel_h", sel_h_nosip3d),
            ]

            for weighted in (False, True):
                if weighted and not hasattr(events, "Generator_weight"):
                    continue

                final_cut = self._register_sequential_event_cutflow(
                    events=events,
                    cut_type="zgammas_mu_munosip3d_w" if weighted else "zgammas_mu_munosip3d",
                    steps=mu_nosip3d_steps,
                    weighted=weighted,
                )
                if not weighted:
                    final_mu_nosip3d_cut = final_cut

            awkward_utils.add_field(events, "pass_allcuts_munosip3d", ak.fill_none(final_mu_nosip3d_cut, False), overwrite=True)

        # --- NEW: 一次註冊 15 種 photon ID 情境對 cutflow("all cuts") 的影響 ---
        def _scenario_has_2gamma(events_, id_key: str):
            """
            以「所有 photon」上的 pass_phid_* flag 定義該情境的 has_2gamma_cand。
            注意：這裡避免重跑 ALP build；目標是研究 photonID 對 cutflow 的影響（尤其 all cuts）
            """
            flag_name = f"pass_{id_key}"
            if flag_name not in events_.Photon.fields:
                return ak.fill_none(ak.num(events_.Photon) < 0, False)  # all False
            return ak.fill_none(ak.sum(events_.Photon[flag_name], axis=1) >= 2, False)

        g_kin_cut = ak.fill_none(ak.sum(events.Photon["pass_ph_kinematic"], axis=1) >= 2, False)

        if self.options.get("zgammas", {}).get("study_photon_id_scenarios", False):
            scenario_keys = [
                # tight
                "phid_custom_tight",
                "phid_custom_extend_tight",
                "phid_sieie_tight",
                "phid_PFECalIso_tight",
                "phid_official_tight",
                # medium
                "phid_custom_medium",
                "phid_custom_extend_medium",
                "phid_sieie_medium",
                "phid_PFECalIso_medium",
                "phid_official_medium",
                # loose
                "phid_custom_loose",
                "phid_custom_extend_loose",
                "phid_sieie_loose",
                "phid_PFECalIso_loose",
                "phid_official_loose",
            ]
            subset = self.options.get("zgammas", {}).get("study_photon_id_scenarios_subset", None)
            if subset:
                scenario_keys = [k for k in scenario_keys if k in set(subset)]

            # 共用既有 cut 定義，只替換 has_2gamma_cand
            for id_key in scenario_keys:
                
                for weighted in (False, True):
                    if weighted and not hasattr(events, "Generator_weight"):
                        continue

                    has_2g_s = _scenario_has_2gamma(events, id_key)
                    final_cut = self._register_sequential_event_cutflow(
                        events=events,
                        cut_type=f"zgammas_{id_key}",
                        steps=[
                            ("N_lep_sel", z_ee_cut | z_mumu_cut),
                            ("trig_cut", trigger_cut),
                            ("lep_pt_cut", trigger_pt_cut),
                            ("has_z_cand", has_z_cand),
                            ("g_kin_cut", g_kin_cut),
                            ("has_2g_cand", has_2g_s),
                            ("sel_h", sel_h),
                        ],
                        weighted=weighted,
                    )

                    if not weighted:
                        awkward_utils.add_field(events, f"pass_allcuts_{id_key}", ak.fill_none(final_cut, False), overwrite=True)

        all_cuts = ee_all_cut | mm_all_cut

        elapsed_time = time.time() - start
        logger.debug("[ZGammaTagger] %s, syst variation : %s, total time to execute select_zgammas: %.6f s" % (self.name, self.current_syst, elapsed_time))

        # Nominal
        return all_cuts, events

        # Motive Study
        # keep_all = ak.num(events.Photon) >= 0
        # return keep_all, events


    def calculate_gen_info(self, zgammas, options):
        """
        Store gen-level info for H -> Z a, Z -> ll, a -> gamma gamma .
        Expected fields from gen_selections.select_x_to_yz:
          ['LeadGenChild','SubleadGenChild','LeadGenChildChild1','LeadGenChildChild2','GenParent']
        """
        # H(25) -> Z(23) + a(9000005), and Z -> l l (l = e/mu/t)
        gen_hza = gen_selections.select_x_to_yz(zgammas.GenPart, 25, 23, 9000005)

        # Robust mapping:
        # - For (23,9000005) with abs(pdgId) sort desc: LeadGenChild=ALP (9000005) OR Z? actually 9000005>23 so LEAD would be ALP.
        #   Your comment says "leading for Z" but current sorting makes ALP lead. So enforce by pdgId here.
        z_cand  = ak.where(abs(gen_hza.LeadGenChild.pdgId) == 23, gen_hza.LeadGenChild, gen_hza.SubleadGenChild)
        a_cand  = ak.where(abs(gen_hza.LeadGenChild.pdgId) == 9000005, gen_hza.LeadGenChild, gen_hza.SubleadGenChild)

        # Z -> ll children: sort by pt to define lead/sublead leptons
        l1 = gen_hza.LeadGenChildChild1
        l2 = gen_hza.LeadGenChildChild2
        lep_lead = ak.where(l1.pt >= l2.pt, l1, l2)
        lep_sub  = ak.where(l1.pt >= l2.pt, l2, l1)
        # PDG convention: e-/mu-/tau- have positive pdgId; anti-leptons have negative pdgId.
        lep_minus = ak.where(l1.pdgId > 0, l1, l2)
        lep_plus = ak.where(l1.pdgId > 0, l2, l1)

        awkward_utils.add_object_fields(
            events=zgammas,
            name="GenHzaHiggs",
            objects=gen_hza.GenParent,
            n_objects=1
        )
        awkward_utils.add_object_fields(
            events=zgammas,
            name="GenHzaZ",
            objects=z_cand,
            n_objects=1
        )
        awkward_utils.add_object_fields(
            events=zgammas,
            name="GenHzaALP",
            objects=a_cand,
            n_objects=1
        )
        awkward_utils.add_object_fields(
            events=zgammas,
            name="GenHzaZLeadLep",
            objects=lep_lead,
            n_objects=1
        )
        awkward_utils.add_object_fields(
            events=zgammas,
            name="GenHzaZSubleadLep",
            objects=lep_sub,
            n_objects=1
        )
        awkward_utils.add_object_fields(
            events=zgammas,
            name="GenHzaZLepMinus",
            objects=lep_minus,
            n_objects=1
        )
        awkward_utils.add_object_fields(
            events=zgammas,
            name="GenHzaZLepPlus",
            objects=lep_plus,
            n_objects=1
        )

        # Optional: store useful dR's (event-level scalars)
        awkward_utils.add_field(
            zgammas,
            "GenHza_dR_ZALP",
            ak.fill_none(z_cand.deltaR(a_cand), DUMMY_VALUE),
            overwrite=True
        )
        awkward_utils.add_field(
            zgammas,
            "GenHza_dR_ll",
            ak.fill_none(lep_lead.deltaR(lep_sub), DUMMY_VALUE),
            overwrite=True
        )

        # (Optional) a(9000005) -> γγ : try to select ALP->γγ separately; if absent, branches will be dummy
        try:
            gen_agam = gen_selections.select_x_to_yz(zgammas.GenPart, 9000005, 22, 22)

            # NEW: robust 取 child（支援 direct 2-body 回傳；legacy 4-body 也仍可用）
            g1 = gen_agam.LeadGenChild
            g2 = gen_agam.SubleadGenChild

            # Store an unordered-by-pt convention in addition to pT-ordered aliases.
            # For a spin-0 a -> gamma gamma decay either photon is equivalent up to
            # cos(theta_a) -> -cos(theta_a); abs(cos(theta_a)) is convention-safe.
            awkward_utils.add_object_fields(
                events=zgammas,
                name="GenALPPhoton1",
                objects=g1,
                n_objects=1
            )
            awkward_utils.add_object_fields(
                events=zgammas,
                name="GenALPPhoton2",
                objects=g2,
                n_objects=1
            )

            # 若 event 沒有 pair，g1/g2 會是 None；排序前先保護
            g1_pt = ak.fill_none(g1.pt, -1.0)
            g2_pt = ak.fill_none(g2.pt, -1.0)

            pho_lead = ak.where(g1_pt >= g2_pt, g1, g2)
            pho_sub  = ak.where(g1_pt >= g2_pt, g2, g1)

            awkward_utils.add_object_fields(
                events=zgammas,
                name="GenALPLeadPho",
                objects=pho_lead,
                n_objects=1
            )
            awkward_utils.add_object_fields(
                events=zgammas,
                name="GenALPSubleadPho",
                objects=pho_sub,
                n_objects=1
            )
            awkward_utils.add_field(
                zgammas,
                "GenALP_dR_gg",
                ak.fill_none(pho_lead.deltaR(pho_sub), DUMMY_VALUE),
                overwrite=True
            )
        except Exception:
            # keep backward compatible; no crash if decay not present / selection logic doesn't find it
            pass

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

        # logger.debug("[select_objects] : Tagger '%s', selecting objects '%s', with the following requirements:" % (tagger_name, "SelectedPhoton"))
        # for cut, value in options.items():
            # logger.debug("\t '%s' : %s" % (cut, str(value)))

        # nominal
        no_cut = photons.pt > 0

        # pt
        pt_cut = photons.pt > options["pt"]

        # eta
        #eta_cut = Tagger.get_range_cut(abs(photons.eta), options["eta"]) | (photons.isScEtaEB | photons.isScEtaEE)
        eta_cut = (photons.isScEtaEB | photons.isScEtaEE) 
        # eta_cut = ((photons.isScEtaEB & (photons.mvaID > options["mvaID_barrel"])) | (photons.isScEtaEE & (photons.mvaID > options["mvaID_endcap"])))

        ph_kinemtaic = pt_cut & eta_cut
        photons = ak.with_field(photons, ph_kinemtaic, "pass_ph_kinematic")

        # photon ID
        phid_custom_tight = []
        phid_custom_extend_tight = []
        phid_sieie_tight = []
        phid_PFECalIso_tight = []
        phid_official_tight = []

        phid_custom_mediate = []
        phid_custom_extend_mediate = []
        phid_sieie_mediate = []
        phid_PFECalIso_mediate = []
        phid_official_mediate = []

        phid_custom_loose = []
        phid_custom_extend_loose = []
        phid_sieie_loose = []
        phid_PFECalIso_loose = []
        phid_official_loose = []

        rho_broadcasted, _ = ak.broadcast_arrays(rho, photons.pt)
        rho = rho_broadcasted
        photon_abs_eta = numpy.abs(photons.eta)

        # NEW: 若兩個 photon dR < 0.3，計算 PFECalIso 時從 ecalPFClusterIso 扣掉彼此 pt
        # 只修正 ecalPFClusterIso（因為需求點名它），且只在 Run3 ID 計算中使用
        ecalPFClusterIso_forID = photons.ecalPFClusterIso
        if int(year) > 2020:
            try:
                pho_v4 = to_momentum4d(photons)
                pairs = ak.combinations(pho_v4, 2, fields=["p1", "p2"])
                dr12 = pairs.p1.deltaR(pairs.p2)
                close = dr12 < 0.3

                # 對每一對 (i,j)，若 close，則 i 扣 j.pt、j 扣 i.pt
                pt_other_to_p1 = ak.where(close, pairs.p2.pt, 0.0)
                pt_other_to_p2 = ak.where(close, pairs.p1.pt, 0.0)

                # 聚合回每顆 photon：sum(pt of other photons within cone)
                # 對 ak.combinations 的輸出，用 local_index 重新定位回原本 photon index
                idx = ak.local_index(photons.pt, axis=1)
                idx_pairs = ak.combinations(idx, 2, fields=["i1", "i2"])

                npho = ak.num(photons.pt, axis=1)
                sum_other = ak.zeros_like(photons.pt, dtype=float)

                # scatter-add: 將 pt_other_to_p1 加到 i1；pt_other_to_p2 加到 i2
                sum_other = ak.where(sum_other >= 0, sum_other, 0.0)  # keep shape, no-op guard
                sum_other = ak.flatten(sum_other, axis=None)  # flatten for add-at
                i1_flat = ak.to_numpy(ak.flatten(idx_pairs.i1))
                i2_flat = ak.to_numpy(ak.flatten(idx_pairs.i2))
                add1 = ak.to_numpy(ak.flatten(pt_other_to_p1))
                add2 = ak.to_numpy(ak.flatten(pt_other_to_p2))

                # event-wise offset for flattened indexing
                counts = ak.to_numpy(npho)
                offsets = numpy.zeros_like(counts)
                offsets[1:] = numpy.cumsum(counts[:-1])
                off_rep = numpy.repeat(offsets, ak.to_numpy(ak.num(idx_pairs.i1)))

                numpy.add.at(sum_other, off_rep + i1_flat, add1)
                numpy.add.at(sum_other, off_rep + i2_flat, add2)

                sum_other = ak.unflatten(sum_other, counts)
                ecalPFClusterIso_forID = ak.where(
                    ecalPFClusterIso_forID - sum_other > 0,
                    ecalPFClusterIso_forID - sum_other,
                    0.0,
                )
            except Exception as _e:
                # 回退：不修正，避免因少數奇怪 event 造成整體崩潰
                ecalPFClusterIso_forID = photons.ecalPFClusterIso

        def _apply_quadratic_ea_corr(iso_values, ea_prefix):
            corrected = iso_values
            eta_regions = [
                (((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0)),   f"{ea_prefix}_EA_EB_1"),
                (((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442)), f"{ea_prefix}_EA_EB_2"),
                (((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0)), f"{ea_prefix}_EA_EE_1"),
                (((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2)),   f"{ea_prefix}_EA_EE_2"),
                (((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3)),   f"{ea_prefix}_EA_EE_3"),
                (((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4)),   f"{ea_prefix}_EA_EE_4"),
                (((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5)),   f"{ea_prefix}_EA_EE_5"),
            ]
            for eta_mask, option_key in eta_regions:
                corrected = ak.where(
                    eta_mask,
                    iso_values
                    - rho * options[option_key][0]
                    - rho**2 * options[option_key][1],
                    corrected,
                )
            return corrected

        hcalPFClusterIso_PUcorr = photons.hcalPFClusterIso
        ecalPFClusterIso_PUcorr = photons.ecalPFClusterIso
        if "pfRelIso03_chg" in photons.fields:
            chiso_PUcorr = photons.pfRelIso03_chg
        elif "pfRelIso03_chg_quadratic" in photons.fields:
            chiso_PUcorr = photons.pfRelIso03_chg_quadratic
        else:
            raise AttributeError("Photon is missing both 'pfRelIso03_chg' and 'pfRelIso03_chg_quadratic'.")

        if int(year) > 2020:
            photons = ak.with_field(photons, photons.pfRelIso03_chg_quadratic, "pfRelIso03_chg")
            if "pfRelIso03_all_quadratic" in photons.fields:
                photons = ak.with_field(photons, photons.pfRelIso03_all_quadratic, "pfRelIso03_all")
            hcalPFClusterIso_PUcorr = _apply_quadratic_ea_corr(
                photons.hcalPFClusterIso, "PFHCalIso"
            )
            ecalPFClusterIso_PUcorr = _apply_quadratic_ea_corr(
                ecalPFClusterIso_forID, "PFECalIso"
            )
            chiso_PUcorr = photons.pfRelIso03_chg_quadratic

        photons = ak.with_field(photons, hcalPFClusterIso_PUcorr, "hcalPFClusterIso_PUcorr")
        photons = ak.with_field(photons, ecalPFClusterIso_PUcorr, "ecalPFClusterIso_PUcorr")
        photons = ak.with_field(photons, chiso_PUcorr, "chiso_PUcorr")

        if int(year) < 2020:
            phid_custom_tight = ak.ones_like(photons.pt) # all true, dummy TODO
            phid_official_tight = ak.ones_like(photons.pt) # all true, dummy TODO

            # NEW: 補齊 medium/loose（Run2 先 dummy all-true，避免欄位缺失）
            phid_custom_mediate = ak.ones_like(photons.pt)
            phid_custom_extend_mediate = ak.ones_like(photons.pt)
            phid_sieie_mediate = ak.ones_like(photons.pt)
            phid_PFECalIso_mediate = ak.ones_like(photons.pt)
            phid_official_mediate = ak.ones_like(photons.pt)

            phid_custom_loose = ak.ones_like(photons.pt)
            phid_custom_extend_loose = ak.ones_like(photons.pt)
            phid_sieie_loose = ak.ones_like(photons.pt)
            phid_PFECalIso_loose = ak.ones_like(photons.pt)
            phid_official_loose = ak.ones_like(photons.pt)

        elif int(year) > 2020:
            # ID
            # id_cut = photons.mvaID_WP80
            ''' 1. Rho corrected H/E '''
            hoe_barrel_cut = photons.hoe_PUcorr < options["tight_hoe_barrel"]
            hoe_endcap_cut = photons.hoe_PUcorr < options["tight_hoe_endcap"]

            ''' 2. Rho corrected PF charged hadron isolation '''
            PFChIso_barrel_cut = photons.chiso_PUcorr < options["tight_PFChIso_barrel"]
            PFChIso_endcap_cut = photons.chiso_PUcorr < options["tight_PFChIso_endcap"]

            ''' 3. Rho corrected PF HCal isolation '''
            # quadratic EA corrections in Run3 : https://indico.cern.ch/event/1204277/contributions/5064356/attachments/2538496/4369369/CutBasedPhotonID_20221031.pdf
            coef_1, coef_2, coef_3 = options["tight_PFHCalIso_barrel"]
            parabola_cut = coef_1 + coef_2 * photons.pt + coef_3 * photons.pt ** 2
            # PFHCalIso_barrel_cut = photons.hcalPFClusterIso < parabola_cut
            PFHCalIso_barrel_cut = (
                ((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0))
                & (
                    photons.hcalPFClusterIso_PUcorr
                    < parabola_cut
                )
            ) | (
                ((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442))
                & (
                    photons.hcalPFClusterIso_PUcorr
                    < parabola_cut
                )
            )

            coef_1, coef_2, coef_3 = options["tight_PFHCalIso_endcap"]
            parabola_cut = coef_1 + coef_2 * photons.pt + coef_3 * photons.pt ** 2
            # PFHCalIso_endcap_cut = photons.hcalPFClusterIso < parabola_cut
            PFHCalIso_endcap_cut = (
                (
                    ((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0))
                    & (
                        photons.hcalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2))
                    & (
                        photons.hcalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3))
                    & (
                        photons.hcalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4))
                    & (
                        photons.hcalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5))
                    & (
                        photons.hcalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
            )

            ''' 4. Sigma IEtaIEta (Discarded) '''
            sieie_barrel_cut = photons.sieie < options["tight_sieie_barrel"]
            sieie_endcap_cut = photons.sieie < options["tight_sieie_endcap"]

            ''' 5. Rho corrected PF ECal isolation (Discarded) '''
            # quadratic EA corrections in Run3 : https://indico.cern.ch/event/1204277/contributions/5064356/attachments/2538496/4369369/CutBasedPhotonID_20221031.pdf
            coef_1, coef_2 = options["tight_PFECalIso_barrel"]
            parabola_cut = coef_1 + coef_2 * photons.pt
            # PFECalIso_barrel_cut = photons.ecalPFClusterIso < parabola_cut
            PFECalIso_barrel_cut = (
                ((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0))
                & (
                    photons.ecalPFClusterIso_PUcorr
                    < parabola_cut
                )
            ) | (
                ((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442))
                & (
                    photons.ecalPFClusterIso_PUcorr
                    < parabola_cut
                )
            )

            coef_1, coef_2 = options["tight_PFECalIso_endcap"]
            parabola_cut = coef_1 + coef_2 * photons.pt
            # PFECalIso_endcap_cut = photons.ecalPFClusterIso < parabola_cut
            PFECalIso_endcap_cut = (
                (
                    ((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0))
                    & (
                        photons.ecalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2))
                    & (
                        photons.ecalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3))
                    & (
                        photons.ecalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4))
                    & (
                        photons.ecalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
                | (
                    ((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5))
                    & (
                        photons.ecalPFClusterIso_PUcorr
                        < parabola_cut
                    )
                )
            )

            # 1, 2, 3
            hoe_cut = (photons.isScEtaEB & hoe_barrel_cut) | (photons.isScEtaEE & hoe_endcap_cut)
            PFChIso_cut = (photons.isScEtaEB & PFChIso_barrel_cut) | (photons.isScEtaEE & PFChIso_endcap_cut)
            PFHCalIso_cut = (photons.isScEtaEB & PFHCalIso_barrel_cut) | (photons.isScEtaEE & PFHCalIso_endcap_cut)
            # 4, 5
            sieie_cut = (photons.isScEtaEB & sieie_barrel_cut) | (photons.isScEtaEE & sieie_endcap_cut)
            PFECalIso_cut = (photons.isScEtaEB & PFECalIso_barrel_cut) | (photons.isScEtaEE & PFECalIso_endcap_cut)

            phid_custom_tight = hoe_cut & PFChIso_cut & PFHCalIso_cut
            phid_custom_extend_tight = hoe_cut & PFChIso_cut & PFHCalIso_cut & sieie_cut
            phid_sieie_tight = sieie_cut
            phid_PFECalIso_tight = PFECalIso_cut
            phid_official_tight = hoe_cut & PFChIso_cut & PFHCalIso_cut & sieie_cut & PFECalIso_cut

            # NEW: 依樣畫葫蘆算 medium / loose（用 options["medium_*"] / ["loose_*"]）
            def _calc_wp(wp: str):
                # H/E
                hoe_barrel = photons.hoe_PUcorr < options[f"{wp}_hoe_barrel"]
                hoe_endcap = photons.hoe_PUcorr < options[f"{wp}_hoe_endcap"]

                # PFChIso
                ch_barrel = photons.chiso_PUcorr < options[f"{wp}_PFChIso_barrel"]
                ch_endcap = photons.chiso_PUcorr < options[f"{wp}_PFChIso_endcap"]

                # PFHCalIso (quadratic threshold + EA corrections)
                c1, c2, c3 = options[f"{wp}_PFHCalIso_barrel"]
                thr_b = c1 + c2 * photons.pt + c3 * photons.pt**2
                pfh_b = (
                    ((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0))
                    & (photons.hcalPFClusterIso_PUcorr < thr_b)
                ) | (
                    ((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442))
                    & (photons.hcalPFClusterIso_PUcorr < thr_b)
                )

                c1, c2, c3 = options[f"{wp}_PFHCalIso_endcap"]
                thr_e = c1 + c2 * photons.pt + c3 * photons.pt**2
                pfh_e = (
                    ((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0))
                    & (photons.hcalPFClusterIso_PUcorr < thr_e)
                ) | (
                    ((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2))
                    & (photons.hcalPFClusterIso_PUcorr < thr_e)
                ) | (
                    ((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3))
                    & (photons.hcalPFClusterIso_PUcorr < thr_e)
                ) | (
                    ((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4))
                    & (photons.hcalPFClusterIso_PUcorr < thr_e)
                ) | (
                    ((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5))
                    & (photons.hcalPFClusterIso_PUcorr < thr_e)
                )

                # sieie
                sie_b = photons.sieie < options[f"{wp}_sieie_barrel"]
                sie_e = photons.sieie < options[f"{wp}_sieie_endcap"]

                # PFECalIso (linear threshold + EA corrections)
                c1, c2 = options[f"{wp}_PFECalIso_barrel"]
                thr_ec_b = c1 + c2 * photons.pt
                pfe_b = (
                    ((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0))
                    & (photons.ecalPFClusterIso_PUcorr < thr_ec_b)
                ) | (
                    ((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442))
                    & (photons.ecalPFClusterIso_PUcorr < thr_ec_b)
                )

                c1, c2 = options[f"{wp}_PFECalIso_endcap"]
                thr_ec_e = c1 + c2 * photons.pt
                pfe_e = (
                    ((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0))
                    & (photons.ecalPFClusterIso_PUcorr < thr_ec_e)
                ) | (
                    ((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2))
                    & (photons.ecalPFClusterIso_PUcorr < thr_ec_e)
                ) | (
                    ((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3))
                    & (photons.ecalPFClusterIso_PUcorr < thr_ec_e)
                ) | (
                    ((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4))
                    & (photons.ecalPFClusterIso_PUcorr < thr_ec_e)
                ) | (
                    ((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5))
                    & (photons.ecalPFClusterIso_PUcorr < thr_ec_e)
                )

                hoe = (photons.isScEtaEB & hoe_barrel) | (photons.isScEtaEE & hoe_endcap)
                ch  = (photons.isScEtaEB & ch_barrel) | (photons.isScEtaEE & ch_endcap)
                pfh = (photons.isScEtaEB & pfh_b) | (photons.isScEtaEE & pfh_e)
                sie = (photons.isScEtaEB & sie_b) | (photons.isScEtaEE & sie_e)
                pfe = (photons.isScEtaEB & pfe_b) | (photons.isScEtaEE & pfe_e)

                return {
                    f"phid_custom_{wp}": hoe & ch & pfh,
                    f"phid_custom_extend_{wp}": hoe & ch & pfh & sie,
                    f"phid_sieie_{wp}": sie,
                    f"phid_PFECalIso_{wp}": pfe,
                    f"phid_official_{wp}": hoe & ch & pfh & sie & pfe,
                }

            # 產生 medium/loose
            m = _calc_wp("medium")
            l = _calc_wp("loose")

            phid_custom_mediate = m["phid_custom_medium"]
            phid_custom_extend_mediate = m["phid_custom_extend_medium"]
            phid_sieie_mediate = m["phid_sieie_medium"]
            phid_PFECalIso_mediate = m["phid_PFECalIso_medium"]
            phid_official_mediate = m["phid_official_medium"]

            phid_custom_loose = l["phid_custom_loose"]
            phid_custom_extend_loose = l["phid_custom_extend_loose"]
            phid_sieie_loose = l["phid_sieie_loose"]
            phid_PFECalIso_loose = l["phid_PFECalIso_loose"]
            phid_official_loose = l["phid_official_loose"]

        # NEW: 把 15 種 photon-ID 結果存成 photon-level flags（供上游 scenario study 使用）
        photons = ak.with_field(photons, phid_custom_tight, "pass_phid_custom_tight")
        photons = ak.with_field(photons, phid_custom_extend_tight, "pass_phid_custom_extend_tight")
        photons = ak.with_field(photons, phid_sieie_tight, "pass_phid_sieie_tight")
        photons = ak.with_field(photons, phid_PFECalIso_tight, "pass_phid_PFECalIso_tight")
        photons = ak.with_field(photons, phid_official_tight, "pass_phid_official_tight")

        photons = ak.with_field(photons, phid_custom_mediate, "pass_phid_custom_medium")
        photons = ak.with_field(photons, phid_custom_extend_mediate, "pass_phid_custom_extend_medium")
        photons = ak.with_field(photons, phid_sieie_mediate, "pass_phid_sieie_medium")
        photons = ak.with_field(photons, phid_PFECalIso_mediate, "pass_phid_PFECalIso_medium")
        photons = ak.with_field(photons, phid_official_mediate, "pass_phid_official_medium")

        photons = ak.with_field(photons, phid_custom_loose, "pass_phid_custom_loose")
        photons = ak.with_field(photons, phid_custom_extend_loose, "pass_phid_custom_extend_loose")
        photons = ak.with_field(photons, phid_sieie_loose, "pass_phid_sieie_loose")
        photons = ak.with_field(photons, phid_PFECalIso_loose, "pass_phid_PFECalIso_loose")
        photons = ak.with_field(photons, phid_official_loose, "pass_phid_official_loose")

        # Custom Photon ID（nominal 用哪個，仍照你現在的策略；study 用 flags 另外算）
        id_cut = phid_custom_tight


        # electron veto
        e_veto_cut = (photons.electronVeto > options["e_veto"])
        # photons = ak.with_field(photons, e_veto_cut, "pass_photon_e_veto")

        photon_ele_idx = ak.where(ak.num(photons.electronIdx, axis=1) == 0, ak.ones_like(photons.pt)*-1, photons.electronIdx)
        new_pho = ak.unflatten(ak.unflatten(ak.flatten(photon_ele_idx), [1]*ak.sum(ak.num(photon_ele_idx))), ak.num(photon_ele_idx, axis=1))
        new_ele = ak.broadcast_arrays(electrons.Idx[:,None], new_pho, depth_limit=2)[0]
        eg_overlap_cut = ~ak.where(
            ak.is_none(electrons.Idx),
            ak.broadcast_arrays(photons.electronIdx, False)[1],
            ak.flatten(ak.any(new_pho[:, :, None] == new_ele, axis=-2), axis=-1)
        ) # some events may have no electrons, so we need to replace None with False

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

        # 回傳：最終 selection mask + 帶 flags 的 photons
        return all_cuts, photons

    def _compute_alp_iso(self, all_photons, alp_cand, dr_cone=0.3):
        """
        使用 ALP 四動量與「所有未經過篩選的 photon（events.Photon）」計算隔離變量：
          sum_{photons} pt，條件：
            - ΔR(ALP, photon) < dr_cone
            - 排除組成 ALP 的兩顆 photon（以原始索引 LeadOrigIndex/SubleadOrigIndex 排除）
        對於沒有有效 ALP 的事件回傳 0。
        """
        n_all = ak.num(all_photons)
        has_any_pho = n_all > 0
        # 事件是否有 ALP 候選（gamma pair）
        has_alp = ~ak.is_none(alp_cand) & ~ak.is_none(alp_cand.ALPCand)
        valid = has_any_pho & has_alp

        dbg = bool(self.options.get("zgammas", {}).get("debug_alp_iso", False))
        if dbg:
            try:
                n_valid = int(ak.sum(valid))
                logger.info(f"[ALP ISO] total_events={len(all_photons)}, valid_events={n_valid}, dr_cone={0.3}")
            except Exception as e:
                logger.info(f"[ALP ISO] debug summary failed: {e}")

        if not ak.any(valid):
            return ak.zeros_like(n_all, dtype=float)

        # 取出有效事件的物件
        all_ph_v = all_photons[valid]
        acand_v = alp_cand[valid]

        # 轉為 Momentum4D
        all_ph_v4 = to_momentum4d(all_ph_v)
        alp_vec = ak.with_name(acand_v.ALPCand, "Momentum4D")

        # ΔR 計算（ALP 與所有 photon）
        alp_b, pho_b = ak.broadcast_arrays(alp_vec, all_ph_v4)
        alp_b = ak.with_name(alp_b, "Momentum4D")
        pho_b = ak.with_name(pho_b, "Momentum4D")
        dr = alp_b.deltaR(pho_b)  # (ev, n_ph)

        # 事件內 photon 的原始索引，以及 ALP 兩顆 photon 在原始 Photon 陣列的索引
        idx_all = ak.local_index(all_ph_v.pt, axis=1)              # (ev, n_ph)
        lead_orig = ak.fill_none(acand_v.LeadOrigIndex, -1)        # (ev,)
        sub_orig  = ak.fill_none(acand_v.SubleadOrigIndex, -1)     # (ev,)

        # 顯式廣播至 (ev, n_ph)
        lead_orig_b = lead_orig[:, None]
        sub_orig_b  = sub_orig[:, None]

        # 排除 ALP 本身的兩顆 photon
        not_parts = (idx_all != lead_orig_b) & (idx_all != sub_orig_b)

        # cone 內且非 ALP 自身兩顆 photon
        iso_mask = (dr < dr_cone) & not_parts

        # 加總 pt
        iso_valid = ak.sum(all_ph_v.pt[iso_mask], axis=1)
        iso_valid = ak.fill_none(iso_valid, 0.0)

        if dbg:
            try:
                # 統計 cone 內（不排除 ALP 兩顆）與排除後的數量
                near_any = ak.sum(dr < dr_cone, axis=1)
                near_after = ak.sum(iso_mask, axis=1)

                n_valid = int(len(iso_valid))
                frac_any = float(ak.sum(near_any > 0)) / max(n_valid, 1)
                frac_after = float(ak.sum(near_after > 0)) / max(n_valid, 1)
                logger.info(f"[ALP ISO] frac(any_in_cone)={frac_any:.3f}, frac(after_excl)={frac_after:.3f}")

                K = int(min(self.options["zgammas"].get("debug_alp_iso_max_events", 5), n_valid))
                for i in range(K):
                    try:
                        npho_i = int(ak.num(all_ph_v[i]))
                        li = int(lead_orig[i])
                        si = int(sub_orig[i])
                        na = int(near_any[i])
                        nb = int(near_after[i])
                        iso_i = float(iso_valid[i])

                        # 逐 event 檢查匹配，避免 axis 錯誤
                        has_lead_match_i = False
                        has_sub_match_i = False
                        if npho_i > 0 and li >= 0 and li < npho_i:
                            has_lead_match_i = True
                        if npho_i > 0 and si >= 0 and si < npho_i:
                            has_sub_match_i = True

                        # 計算 dr(ALP, lead/sub)（若 index 合法）
                        dr_lead = None
                        dr_sub = None
                        if has_lead_match_i:
                            dr_lead = float(dr[i][li])
                        if has_sub_match_i:
                            dr_sub = float(dr[i][si])

                        logger.info(f"[ALP ISO][ev#{i}] nPhoAll={npho_i} leadOrig={li} subOrig={si} "
                                    f"nInCone={na} nInConeExcl={nb} iso={iso_i:.3f} "
                                    f"hasLeadMatch={has_lead_match_i} hasSubMatch={has_sub_match_i} "
                                    f"dr(ALP,lead)={dr_lead} dr(ALP,sub)={dr_sub}")

                        # 印出最小的幾個 ΔR 幫助理解幾何關係（若 npho_i>0）
                        if npho_i > 0:
                            dr_sorted = ak.to_list(ak.sort(dr[i]))
                            logger.debug(f"[ALP ISO][ev#{i}] minDRs={dr_sorted[:min(5, len(dr_sorted))]}")

                    except Exception as ie:
                        logger.debug(f"[ALP ISO][ev#{i}] debug print failed: {ie}")
            except Exception as e:
                logger.info(f"[ALP ISO] debug block failed: {e}")

        # 回填到完整事件長度
        valid_np = ak.to_numpy(valid)
        iso_valid_np = ak.to_numpy(iso_valid)
        iso_full_np = numpy.zeros(valid_np.shape[0], dtype=iso_valid_np.dtype)
        iso_full_np[valid_np] = iso_valid_np
        return ak.Array(iso_full_np)

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
            # 仍附上 refit 需要的「裸輕子」與「FSR=0」欄位，確保下游 z_cand 一致
            zero = ak.zeros_like(leptons.pt)
            return self._attach_refit_fields(leptons, leptons, zero, zero, zero)

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

        # 附上 Z mass-constraint refit 需要的輸入：dressed 之前的「裸輕子」四動量
        # 與這顆輕子實際被指派到的 FSR 光子（沒有時為 0）。如此下游選到的 z_cand
        # 不論哪一對都能取得 refit 所需的 bare/FSR 資訊。
        corrected = self._attach_refit_fields(
            corrected, leptons, best_ph.pt, best_ph.eta, best_ph.phi
        )

        return corrected

    def _attach_refit_fields(self, leptons_dressed, bare_leptons, fsr_pt, fsr_eta, fsr_phi):
        """把 Z-refit 所需的 bare 輕子 4-momentum 與指派的 FSR 光子記在輕子上。"""
        out = ak.with_field(leptons_dressed, bare_leptons.pt,   "bare_pt")
        out = ak.with_field(out,             bare_leptons.eta,  "bare_eta")
        out = ak.with_field(out,             bare_leptons.phi,  "bare_phi")
        out = ak.with_field(out,             bare_leptons.mass, "bare_mass")
        out = ak.with_field(out, fsr_pt,  "fsr_pt")
        out = ak.with_field(out, fsr_eta, "fsr_eta")
        out = ak.with_field(out, fsr_phi, "fsr_phi")
        return out
