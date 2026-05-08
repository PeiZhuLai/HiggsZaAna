import copy
import logging
import numpy
import awkward as ak
import vector

vector.register_awkward()

logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger
from higgs_dna.utils import awkward_utils, misc_utils

Z_MASS = 91.1876
DUMMY_VALUE = -999.0

MUON_TNP_FIELDS = [
    "pt",
    "eta",
    "phi",
    "mass",
    "charge",
    "ptErr",
    "ptE_error",
    "pfRelIso03_all",
    "miniPFRelIso_all",
    "dxy",
    "dz",
    "sip3d",
    "tightId",
    "mediumPromptId",
    "mediumId",
    "looseId",
    "highPtId",
    "isGlobal",
    "isTracker",
    "nTrackerLayers",
    "id",
    "Idx",
]

PHOTON_TNP_FIELDS = [
    "pt",
    "eta",
    "phi",
    "mass",
    "mvaID",
    "mvaID_WP80",
    "mvaID_WP90",
    "energyRaw",
    "energyErr",
    "r9",
    "sieie",
    "hoe",
    "hoe_PUcorr",
    "hcalPFClusterIso",
    "ecalPFClusterIso",
    "pfRelIso03_chg",
    "pfRelIso03_all",
    "pfRelIso03_chg_quadratic",
    "pfRelIso03_all_quadratic",
    "sieip",
    "etaWidth",
    "phiWidth",
    "s4",
    "trkSumPtHollowConeDR03",
    "trkSumPtSolidConeDR04",
    "pfChargedIso",
    "pfChargedIsoWorstVtx",
    "esEffSigmaRR",
    "esEnergyOverRawE",
    "electronVeto",
    "pixelSeed",
    "isScEtaEB",
    "isScEtaEE",
    "electronIdx",
    "jetIdx",
    "genPartFlav",
    "seedGain",
    "origIndex",
    "pass_ph_kinematic",
    "pass_phid_custom_tight",
]

P4_FIELDS = ["pt", "eta", "phi", "mass"]

DEFAULT_OPTIONS = {
    "dimuon": {
        "single_muon_lead_pt": 25.0,
        "single_muon_sublead_pt": 4.0,
        "double_muon_lead_pt": 20.0,
        "double_muon_sublead_pt": 10.0,
    },
    "photons": {
        "use_central_nano": True,
        "pt": 10.0,
        "eta": [
            [0.0, 1.4442],
            [1.566, 2.5],
        ],
        "tight_hoe_barrel": 0.0417588,
        "tight_hoe_endcap": 0.00254267,
        "tight_PFChIso_barrel": 0.316306,
        "tight_PFChIso_endcap": 0.292664,
        "tight_PFHCalIso_barrel": [0.39057, 0.0100547, 5.78332e-05],
        "tight_PFHCalIso_endcap": [0.0292617, 0.0116989, 7.47603e-05],
        "tight_sieie_barrel": 0.00999299,
        "tight_sieie_endcap": 0.0268702,
        "tight_PFECalIso_barrel": [0.14189, 0.000652035],
        "tight_PFECalIso_endcap": [1.04269, 0.000195486],
        "hoe_EA_EB_1": [0.00198598, -1.15014e-05],
        "hoe_EA_EB_2": [0.00249958, -1.15014e-05],
        "hoe_EA_EE_1": [0.00302416, -1.51973e-05],
        "hoe_EA_EE_2": [0.00349404, -1.49651e-05],
        "hoe_EA_EE_3": [0.00334674, -1.47232e-05],
        "hoe_EA_EE_4": [0.00416787, -2.13958e-05],
        "hoe_EA_EE_5": [0.00458106, -2.80795e-05],
        "PFCHIso_EA_EB_1": [0.0342898, -0.000103508],
        "PFCHIso_EA_EB_2": [0.0281424, -3.1494e-05],
        "PFCHIso_EA_EE_1": [0.0288533, -6.66148e-05],
        "PFCHIso_EA_EE_2": [0.028789, -6.84993e-05],
        "PFCHIso_EA_EE_3": [0.0264064, -8.89189e-05],
        "PFCHIso_EA_EE_4": [0.025587, -5.90178e-05],
        "PFCHIso_EA_EE_5": [0.0224817, -4.22712e-05],
        "PFECalIso_EA_EB_1": [0.0866519, -0.0002296],
        "PFECalIso_EA_EB_2": [0.0730397, -0.0002134],
        "PFECalIso_EA_EE_1": [0.0542479, -0.00010934],
        "PFECalIso_EA_EE_2": [0.0486181, -6.2098e-05],
        "PFECalIso_EA_EE_3": [0.0412923, -9.63732e-06],
        "PFECalIso_EA_EE_4": [0.03555, -5.79549e-05],
        "PFECalIso_EA_EE_5": [0.0360895, -7.28546e-06],
        "PFHCalIso_EA_EB_1": [0.17005, -0.000835],
        "PFHCalIso_EA_EB_2": [0.208571, -0.000905],
        "PFHCalIso_EA_EE_1": [0.246494, -0.000722],
        "PFHCalIso_EA_EE_2": [0.306529, -0.000608],
        "PFHCalIso_EA_EE_3": [0.322673, -0.00075],
        "PFHCalIso_EA_EE_4": [0.315793, -0.000795],
        "PFHCalIso_EA_EE_5": [0.36531, -0.000439],
        "e_veto": 0.5,
    },
        "fsr": {
            "min_muon_gamma_dr": 0.8,
            "near_muon_gamma_dr_min": 0.4,
        "photon_mva_min": -0.7,
        "zmmg_mass": [80.0, 100.0],
        "mass_sum_max": 180.0,
    },
        "muons": {
            "pt": 4.0,
            "eta": 2.4,
        },
    "trigger": {
        "single_muon": [
            "HLT_IsoMu24",
        ],
        "double_muon": [
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        ],
    },
}


class TnPZmmgTagger(Tagger):
    def __init__(self, name="TnPZmmgTagger", options=None, is_data=None, year=None):
        super(TnPZmmgTagger, self).__init__(name, options, is_data, year)

        if options is None:
            options = {}
        if not options:
            self.options = copy.deepcopy(DEFAULT_OPTIONS)
        else:
            self.options = misc_utils.update_dict(
                original=DEFAULT_OPTIONS,
                new=options,
            )

    def calculate_selection(self, events):
        year = self.get_analysis_year(events)
        rho = self.get_rho(events)
        event_npv = self.get_event_npv(events)
        awkward_utils.add_field(
            events,
            "fixedGridRhoAll",
            ak.fill_none(rho, 0.0),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "rho",
            ak.fill_none(rho, 0.0),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "event_nPV",
            ak.fill_none(event_npv, DUMMY_VALUE),
            overwrite=True,
        )
        if "PV_z" in events.fields:
            awkward_utils.add_field(
                events,
                "PV",
                ak.fill_none(events.PV_z, DUMMY_VALUE),
                overwrite=True,
            )
            awkward_utils.add_field(
                events,
                "PV_z",
                ak.fill_none(events.PV_z, DUMMY_VALUE),
                overwrite=True,
            )

        muon_selection = self.select_muons(events.Muon, self.options["muons"])
        muons = events.Muon[muon_selection]
        muon_idx = ak.local_index(events.Muon.pt, axis=1)[muon_selection]
        muons = ak.with_field(muons, muon_idx, "Idx")
        muons = ak.with_field(muons, ak.ones_like(muons.pt) * 13, "id")
        muons = ak.with_field(muons, muons.ptErr, "ptE_error")
        muons = muons[ak.argsort(muons.pt, ascending=False, axis=1)]
        muons = ak.Array(muons, with_name="Momentum4D")
        awkward_utils.add_field(events, "SelectedMuon", muons, overwrite=True)

        photon_selection, photons_with_flags = self.select_photons(
            photons=events.Photon,
            rho=rho,
            year=year,
        )
        awkward_utils.add_field(events, "Photon", photons_with_flags, overwrite=True)

        photons = events.Photon[photon_selection]
        photon_idx = ak.local_index(events.Photon.pt, axis=1)[photon_selection]
        photons = ak.with_field(photons, photon_idx, "origIndex")
        photons = photons[ak.argsort(photons.pt, ascending=False, axis=1)]
        if "mass" not in photons.fields:
            photons = ak.with_field(photons, ak.zeros_like(photons.pt), "mass")
        photons = ak.Array(photons, with_name="Momentum4D")
        awkward_utils.add_field(events, "SelectedPhoton", photons, overwrite=True)

        awkward_utils.add_field(events, "n_muons", ak.num(muons), overwrite=True)
        awkward_utils.add_field(events, "n_photons", ak.num(photons), overwrite=True)

        if (not self.is_data) and "GenVtx_z" in events.fields and "PV_z" in events.fields:
            awkward_utils.add_field(
                events,
                "dZ",
                events.GenVtx_z - events.PV_z,
                overwrite=True,
            )

        dimuon_candidates = ak.combinations(
            muons,
            2,
            fields=["LeadMuon", "SubleadMuon"],
        )
        os_cut = (dimuon_candidates.LeadMuon.charge * dimuon_candidates.SubleadMuon.charge) == -1
        dimuon_candidates = dimuon_candidates[os_cut]
        dimuon_lead_pt = ak.where(
            dimuon_candidates.LeadMuon.pt >= dimuon_candidates.SubleadMuon.pt,
            dimuon_candidates.LeadMuon.pt,
            dimuon_candidates.SubleadMuon.pt,
        )
        dimuon_sublead_pt = ak.where(
            dimuon_candidates.LeadMuon.pt >= dimuon_candidates.SubleadMuon.pt,
            dimuon_candidates.SubleadMuon.pt,
            dimuon_candidates.LeadMuon.pt,
        )
        single_trigger_cut, double_trigger_cut, trigger_cut = self.get_dimuon_trigger_decisions(events)
        single_trigger_broadcast = ak.broadcast_arrays(single_trigger_cut, dimuon_lead_pt)[0]
        double_trigger_broadcast = ak.broadcast_arrays(double_trigger_cut, dimuon_lead_pt)[0]
        single_muon_pt_cut = (
            (dimuon_lead_pt > self.options["dimuon"]["single_muon_lead_pt"])
            & (dimuon_sublead_pt > self.options["dimuon"]["single_muon_sublead_pt"])
        )
        double_muon_pt_cut = (
            (dimuon_lead_pt > self.options["dimuon"]["double_muon_lead_pt"])
            & (dimuon_sublead_pt > self.options["dimuon"]["double_muon_sublead_pt"])
        )
        trigger_dependent_muon_pt_cut = (
            (single_trigger_broadcast & single_muon_pt_cut)
            | (double_trigger_broadcast & double_muon_pt_cut)
        )
        dimuon_candidates = dimuon_candidates[trigger_dependent_muon_pt_cut]
        dimuon_candidates["Dimuon"] = (
            dimuon_candidates.LeadMuon + dimuon_candidates.SubleadMuon
        )

        zmmg_candidates_all = ak.cartesian(
            {
                "DimuonCand": dimuon_candidates,
                "ZPhoton": photons,
            },
            axis=1,
        )
        zmmg_candidates_all["LeadMuon"] = zmmg_candidates_all.DimuonCand.LeadMuon
        zmmg_candidates_all["SubleadMuon"] = zmmg_candidates_all.DimuonCand.SubleadMuon
        zmmg_candidates_all["Dimuon"] = zmmg_candidates_all.DimuonCand.Dimuon
        zmmg_candidates_all["Zmmg"] = zmmg_candidates_all.Dimuon + zmmg_candidates_all.ZPhoton

        lead_muon_dr = zmmg_candidates_all.LeadMuon.deltaR(zmmg_candidates_all.ZPhoton)
        sublead_muon_dr = zmmg_candidates_all.SubleadMuon.deltaR(zmmg_candidates_all.ZPhoton)
        lead_is_near = lead_muon_dr <= sublead_muon_dr

        zmmg_candidates_all["NearMuon"] = self.choose_object(
            lead_is_near,
            zmmg_candidates_all.LeadMuon,
            zmmg_candidates_all.SubleadMuon,
        )
        zmmg_candidates_all["FarMuon"] = self.choose_object(
            lead_is_near,
            zmmg_candidates_all.SubleadMuon,
            zmmg_candidates_all.LeadMuon,
        )
        zmmg_candidates_all["TagMuon1"] = zmmg_candidates_all.LeadMuon
        zmmg_candidates_all["TagMuon2"] = zmmg_candidates_all.SubleadMuon
        zmmg_candidates_all["ProbePhoton"] = zmmg_candidates_all.ZPhoton
        zmmg_candidates_all["LeadMuonGammaDR"] = lead_muon_dr
        zmmg_candidates_all["SubleadMuonGammaDR"] = sublead_muon_dr
        zmmg_candidates_all["TagMuon1IsNear"] = lead_is_near
        zmmg_candidates_all["TagMuon2IsNear"] = ~lead_is_near
        zmmg_candidates_all["minMuonGammaDR"] = ak.where(
            lead_is_near,
            lead_muon_dr,
            sublead_muon_dr,
        )
        zmmg_candidates_all["farMuonGammaDR"] = ak.where(
            lead_is_near,
            sublead_muon_dr,
            lead_muon_dr,
        )

        near_dr_max_cut = (
            zmmg_candidates_all.minMuonGammaDR < self.options["fsr"]["min_muon_gamma_dr"]
        )
        near_dr_min_cut = (
            zmmg_candidates_all.minMuonGammaDR > self.options["fsr"]["near_muon_gamma_dr_min"]
        )
        zmmg_mass_cut = (
            zmmg_candidates_all.Zmmg.mass > self.options["fsr"]["zmmg_mass"][0]
        ) & (
            zmmg_candidates_all.Zmmg.mass < self.options["fsr"]["zmmg_mass"][1]
        )
        mass_sum_cut = (
            zmmg_candidates_all.Dimuon.mass + zmmg_candidates_all.Zmmg.mass
        ) < self.options["fsr"]["mass_sum_max"]

        fsr_candidate_cut = (
            near_dr_max_cut
            & near_dr_min_cut
            & zmmg_mass_cut
            & mass_sum_cut
        )
        zmmg_candidates = zmmg_candidates_all[fsr_candidate_cut]
        # If an event has multiple mu-mu-gamma candidates, keep the one closest to the PDG Z mass.
        zmmg_candidates = zmmg_candidates[
            ak.argsort(
                numpy.abs(zmmg_candidates.Zmmg.mass - Z_MASS),
                ascending=True,
                axis=1,
            )
        ]

        has_dimuon = ak.num(dimuon_candidates) >= 1
        has_photon = ak.num(photons) >= 1
        candidate_cut = ak.num(zmmg_candidates) >= 1
        dimuon_cut = trigger_cut & has_dimuon
        photon_cut = dimuon_cut & has_photon
        presel_cut = trigger_cut & candidate_cut
        near_dr_max_event_cut = photon_cut & (ak.num(zmmg_candidates_all[near_dr_max_cut]) >= 1)
        near_dr_min_event_cut = photon_cut & (
            ak.num(zmmg_candidates_all[near_dr_max_cut & near_dr_min_cut]) >= 1
        )
        zmmg_mass_event_cut = photon_cut & (
            ak.num(
                zmmg_candidates_all[
                    near_dr_max_cut
                    & near_dr_min_cut
                    & zmmg_mass_cut
                ]
            ) >= 1
        )
        mass_sum_event_cut = photon_cut & (
            ak.num(zmmg_candidates_all[fsr_candidate_cut]) >= 1
        )

        candidate_trigger_cut = ak.broadcast_arrays(
            trigger_cut,
            zmmg_candidates.Zmmg.mass,
        )[0]
        best_candidate = ak.firsts(zmmg_candidates[candidate_trigger_cut])

        for field in [
            "LeadMuon",
            "SubleadMuon",
            "NearMuon",
            "FarMuon",
            "TagMuon1",
            "TagMuon2",
            "Dimuon",
            "ZPhoton",
            "ProbePhoton",
            "Zmmg",
        ]:
            events[field] = best_candidate[field]
        awkward_utils.add_field(
            events,
            "minMuonGammaDR",
            ak.fill_none(best_candidate.minMuonGammaDR, DUMMY_VALUE),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "farMuonGammaDR",
            ak.fill_none(best_candidate.farMuonGammaDR, DUMMY_VALUE),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "passing_dimuon_trigger",
            ak.fill_none(trigger_cut, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "has_dimuon_candidate",
            ak.fill_none(has_dimuon, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "has_photon_candidate",
            ak.fill_none(has_photon, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "has_zmmg_candidate",
            ak.fill_none(candidate_cut, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "pass_tnp_presel",
            ak.fill_none(presel_cut, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "tag_muon_1_is_near",
            ak.fill_none(best_candidate.TagMuon1IsNear, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "tag_muon_2_is_near",
            ak.fill_none(best_candidate.TagMuon2IsNear, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "tag_muon_1_probe_photon_dr",
            ak.fill_none(best_candidate.LeadMuonGammaDR, DUMMY_VALUE),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "tag_muon_2_probe_photon_dr",
            ak.fill_none(best_candidate.SubleadMuonGammaDR, DUMMY_VALUE),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "zmmg_minus_dimuon_mass",
            ak.fill_none(best_candidate.Zmmg.mass - best_candidate.Dimuon.mass, DUMMY_VALUE),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "abs_zmmg_mass_minus_z",
            ak.fill_none(numpy.abs(best_candidate.Zmmg.mass - Z_MASS), DUMMY_VALUE),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "z_mumu",
            ak.fill_none(presel_cut, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "z_ee",
            ak.zeros_like(events.run, dtype=bool),
            overwrite=True,
        )

        probe_photon_e_veto_value = best_candidate.ProbePhoton.electronVeto
        probe_photon_pass_csev = (
            best_candidate.ProbePhoton.electronVeto > self.options["photons"]["e_veto"]
        )
        probe_photon_pass_pixel_veto = best_candidate.ProbePhoton.pixelSeed < 0.5
        awkward_utils.add_field(
            events,
            "probe_photon_e_veto",
            ak.fill_none(probe_photon_e_veto_value, DUMMY_VALUE),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "probe_photon_pass_csev",
            ak.fill_none(probe_photon_pass_csev, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "probe_photon_pass_pixel_veto",
            ak.fill_none(probe_photon_pass_pixel_veto, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "probe_photon_pass_mva_min",
            ak.fill_none(
                best_candidate.ProbePhoton.mvaID > self.options["fsr"]["photon_mva_min"],
                False,
            ),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "probe_photon_lep_near_dR",
            ak.fill_none(best_candidate.minMuonGammaDR, DUMMY_VALUE),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "probe_photon_lep_far_dR",
            ak.fill_none(best_candidate.farMuonGammaDR, DUMMY_VALUE),
            overwrite=True,
        )
        if int(year) < 2020:
            awkward_utils.add_field(
                events,
                "probe_photon_chiso",
                ak.fill_none(best_candidate.ProbePhoton.pfRelIso03_chg, DUMMY_VALUE),
                overwrite=True,
            )
            awkward_utils.add_field(
                events,
                "probe_photon_alliso",
                ak.fill_none(best_candidate.ProbePhoton.pfRelIso03_all, DUMMY_VALUE),
                overwrite=True,
            )
        else:
            awkward_utils.add_field(
                events,
                "probe_photon_chiso",
                ak.fill_none(best_candidate.ProbePhoton.pfRelIso03_chg_quadratic, DUMMY_VALUE),
                overwrite=True,
            )
            awkward_utils.add_field(
                events,
                "probe_photon_alliso",
                ak.fill_none(best_candidate.ProbePhoton.pfRelIso03_all_quadratic, DUMMY_VALUE),
                overwrite=True,
            )

        self.add_flat_object_fields(events, "tag_muon_1", best_candidate.TagMuon1, MUON_TNP_FIELDS)
        self.add_flat_object_fields(events, "tag_muon_2", best_candidate.TagMuon2, MUON_TNP_FIELDS)
        self.add_flat_object_fields(events, "probe_photon", best_candidate.ProbePhoton, PHOTON_TNP_FIELDS)
        self.add_flat_object_fields(events, "dimuon", best_candidate.Dimuon, P4_FIELDS)
        self.add_flat_object_fields(events, "zmmg", best_candidate.Zmmg, P4_FIELDS)

        self.register_cuts(
            names=[
                "passing dimuon HLT",
                "OSSF dimuon with trigger-dependent mu pT",
                "photon",
                "min dR(mu,gamma) < 0.8",
                "dR(gamma, near muon) > 0.4",
                "80 < m_mumugamma < 100",
                "m_mumu + m_mumugamma < 180",
                "all",
            ],
            results=[
                trigger_cut,
                dimuon_cut,
                photon_cut,
                near_dr_max_event_cut,
                near_dr_min_event_cut,
                zmmg_mass_event_cut,
                mass_sum_event_cut,
                presel_cut,
            ],
            cut_type="event",
        )

        return presel_cut, events

    def get_analysis_year(self, events):
        if self.year is not None:
            return str(self.year)[:4]

        photon_fields = events.Photon.fields
        if "hoe_PUcorr" in photon_fields or "pfRelIso03_chg_quadratic" in photon_fields:
            return "2022"

        return "2018"

    def get_rho(self, events):
        if self.options["photons"]["use_central_nano"]:
            if "fixedGridRhoAll" in events.fields:
                return events.fixedGridRhoAll
            if "fixedGridRhoFastjetAll" in events.fields:
                return events.fixedGridRhoFastjetAll
            if "Rho_fixedGridRhoAll" in events.fields:
                return events.Rho_fixedGridRhoAll

            message = (
                "[TnPZmmgTagger : calculate_selection] Rho not found in central nanoAOD."
            )
            logger.error(message)
            raise RuntimeError(message)

        return ak.ones_like(events.run, dtype=numpy.float64)

    def get_event_npv(self, events):
        if "PV_npvsGood" in events.fields:
            return events.PV_npvsGood
        if "PV_npvs" in events.fields:
            return events.PV_npvs
        if "PV" in events.fields and hasattr(events.PV, "fields"):
            if "npvsGood" in events.PV.fields:
                return events.PV.npvsGood
            if "npvs" in events.PV.fields:
                return events.PV.npvs
        return ak.ones_like(events.run, dtype=numpy.float64) * DUMMY_VALUE

    def select_muons(self, muons, options):
        pt_cut = muons.pt > options["pt"]
        eta_cut = numpy.abs(muons.eta) < options["eta"]
        if "mediumPromptId" in muons.fields:
            id_cut = muons.mediumPromptId == True
            id_name = "mediumPrompt ID"
        elif "mediumId" in muons.fields:
            logger.warning(
                "[TnPZmmgTagger] Muon.mediumPromptId is unavailable; falling back to Muon.mediumId."
            )
            id_cut = muons.mediumId == True
            id_name = "medium ID (fallback)"
        else:
            logger.error("[TnPZmmgTagger] Neither Muon.mediumPromptId nor Muon.mediumId is available.")
            raise RuntimeError("Muon.mediumPromptId is required for the Zmmg TnP muon selection.")

        all_cuts = pt_cut & eta_cut & id_cut

        self.register_cuts(
            names=["pt", "eta", id_name, "all"],
            results=[pt_cut, eta_cut, id_cut, all_cuts],
            cut_type="muon",
        )

        return all_cuts

    def select_photons(self, photons, rho, year):
        if int(year) > 2020:
            if (
                "pfRelIso03_chg" not in photons.fields
                and "pfRelIso03_chg_quadratic" in photons.fields
            ):
                photons = ak.with_field(
                    photons,
                    photons.pfRelIso03_chg_quadratic,
                    "pfRelIso03_chg",
                )
            if (
                "pfRelIso03_all" not in photons.fields
                and "pfRelIso03_all_quadratic" in photons.fields
            ):
                photons = ak.with_field(
                    photons,
                    photons.pfRelIso03_all_quadratic,
                    "pfRelIso03_all",
                )

        _, photons_with_flags = self._select_photons_with_flags(
            photons=photons,
            options=self.options["photons"],
            rho=rho,
            year=year,
        )

        photon_selection = (
            photons_with_flags.pass_ph_kinematic
            & photons_with_flags.pass_phid_custom_tight
        )
        return photon_selection, photons_with_flags

    def _select_photons_with_flags(self, photons, options, rho, year):
        no_cut = photons.pt > 0
        pt_cut = photons.pt > options["pt"]
        eta_cut = photons.isScEtaEB | photons.isScEtaEE

        ph_kinemtaic = pt_cut & eta_cut
        photons = ak.with_field(photons, ph_kinemtaic, "pass_ph_kinematic")

        phid_custom_tight = []

        rho_broadcasted, _ = ak.broadcast_arrays(rho, photons.pt)
        rho = rho_broadcasted
        photon_abs_eta = numpy.abs(photons.eta)

        def _apply_quadratic_ea_corr(iso_values, ea_prefix):
            corrected = iso_values
            eta_regions = [
                (((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0)), f"{ea_prefix}_EA_EB_1"),
                (((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442)), f"{ea_prefix}_EA_EB_2"),
                (((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0)), f"{ea_prefix}_EA_EE_1"),
                (((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2)), f"{ea_prefix}_EA_EE_2"),
                (((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3)), f"{ea_prefix}_EA_EE_3"),
                (((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4)), f"{ea_prefix}_EA_EE_4"),
                (((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5)), f"{ea_prefix}_EA_EE_5"),
            ]
            for eta_mask, option_key in eta_regions:
                corrected = ak.where(
                    eta_mask,
                    iso_values - rho * options[option_key][0] - rho**2 * options[option_key][1],
                    corrected,
                )
            return corrected

        hcalPFClusterIso_PUcorr = photons.hcalPFClusterIso
        if "pfRelIso03_chg" in photons.fields:
            chiso_PUcorr = photons.pfRelIso03_chg
        elif "pfRelIso03_chg_quadratic" in photons.fields:
            chiso_PUcorr = photons.pfRelIso03_chg_quadratic
        else:
            raise AttributeError(
                "Photon is missing both 'pfRelIso03_chg' and 'pfRelIso03_chg_quadratic'."
            )

        if int(year) > 2020:
            photons = ak.with_field(photons, photons.pfRelIso03_chg_quadratic, "pfRelIso03_chg")
            if "pfRelIso03_all_quadratic" in photons.fields:
                photons = ak.with_field(photons, photons.pfRelIso03_all_quadratic, "pfRelIso03_all")
            hcalPFClusterIso_PUcorr = _apply_quadratic_ea_corr(
                photons.hcalPFClusterIso, "PFHCalIso"
            )
            chiso_PUcorr = photons.pfRelIso03_chg_quadratic

        photons = ak.with_field(photons, hcalPFClusterIso_PUcorr, "hcalPFClusterIso_PUcorr")
        photons = ak.with_field(photons, chiso_PUcorr, "chiso_PUcorr")

        if int(year) < 2020:
            phid_custom_tight = ak.ones_like(photons.pt)
        elif int(year) > 2020:
            hoe_barrel_cut = photons.hoe_PUcorr < options["tight_hoe_barrel"]
            hoe_endcap_cut = photons.hoe_PUcorr < options["tight_hoe_endcap"]

            PFChIso_barrel_cut = photons.chiso_PUcorr < options["tight_PFChIso_barrel"]
            PFChIso_endcap_cut = photons.chiso_PUcorr < options["tight_PFChIso_endcap"]

            coef_1, coef_2, coef_3 = options["tight_PFHCalIso_barrel"]
            parabola_cut = coef_1 + coef_2 * photons.pt + coef_3 * photons.pt**2
            PFHCalIso_barrel_cut = (
                ((photon_abs_eta > 0.0) & (photon_abs_eta < 1.0))
                & (photons.hcalPFClusterIso_PUcorr < parabola_cut)
            ) | (
                ((photon_abs_eta > 1.0) & (photon_abs_eta < 1.4442))
                & (photons.hcalPFClusterIso_PUcorr < parabola_cut)
            )

            coef_1, coef_2, coef_3 = options["tight_PFHCalIso_endcap"]
            parabola_cut = coef_1 + coef_2 * photons.pt + coef_3 * photons.pt**2
            PFHCalIso_endcap_cut = (
                ((photon_abs_eta > 1.566) & (photon_abs_eta < 2.0))
                & (photons.hcalPFClusterIso_PUcorr < parabola_cut)
            ) | (
                ((photon_abs_eta > 2.0) & (photon_abs_eta < 2.2))
                & (photons.hcalPFClusterIso_PUcorr < parabola_cut)
            ) | (
                ((photon_abs_eta > 2.2) & (photon_abs_eta < 2.3))
                & (photons.hcalPFClusterIso_PUcorr < parabola_cut)
            ) | (
                ((photon_abs_eta > 2.3) & (photon_abs_eta < 2.4))
                & (photons.hcalPFClusterIso_PUcorr < parabola_cut)
            ) | (
                ((photon_abs_eta > 2.4) & (photon_abs_eta < 2.5))
                & (photons.hcalPFClusterIso_PUcorr < parabola_cut)
            )

            hoe_cut = (photons.isScEtaEB & hoe_barrel_cut) | (photons.isScEtaEE & hoe_endcap_cut)
            PFChIso_cut = (photons.isScEtaEB & PFChIso_barrel_cut) | (
                photons.isScEtaEE & PFChIso_endcap_cut
            )
            PFHCalIso_cut = (photons.isScEtaEB & PFHCalIso_barrel_cut) | (
                photons.isScEtaEE & PFHCalIso_endcap_cut
            )

            phid_custom_tight = hoe_cut & PFChIso_cut & PFHCalIso_cut

        photons = ak.with_field(photons, phid_custom_tight, "pass_phid_custom_tight")

        id_cut = phid_custom_tight

        cut_results = [no_cut, pt_cut, eta_cut, id_cut]
        all_cuts = photons.pt > 0
        for i, cut in enumerate(cut_results):
            all_cuts = all_cuts & cut
            if i == 0:
                cut_results[i] = ak.sum(cut, axis=1) > 1
            else:
                cut_results[i] = (ak.sum(cut, axis=1) > 1) & cut_results[i - 1]

        return all_cuts, photons

    def choose_object(self, condition, first, second):
        fields = first.fields
        data = {
            field: ak.where(condition, first[field], second[field])
            for field in fields
        }
        return ak.zip(data, with_name="Momentum4D")

    def add_flat_object_fields(self, events, name, objects, fields):
        available_fields = [field for field in fields if field in objects.fields]
        if available_fields:
            awkward_utils.add_object_fields(
                events=events,
                name=name,
                objects=ak.singletons(objects),
                n_objects=1,
                dummy_value=DUMMY_VALUE,
                fields=available_fields,
                overwrite=True,
            )

        derived_fields = [
            field for field in fields if field not in available_fields and hasattr(objects, field)
        ]
        for field in derived_fields:
            awkward_utils.add_field(
                events,
                f"{name}_{field}",
                ak.fill_none(getattr(objects, field), DUMMY_VALUE),
                overwrite=True,
            )

    def passing_dimuon_trigger(self, events):
        return self.get_dimuon_trigger_decisions(events)[2]

    def get_dimuon_trigger_decisions(self, events):
        trigger_paths = self.options["trigger"]
        if isinstance(trigger_paths, dict) and not (
            "single_muon" in trigger_paths or "double_muon" in trigger_paths
        ):
            year = self.get_analysis_year(events)
            trigger_paths = trigger_paths.get(year, trigger_paths.get("default", trigger_paths))

        if isinstance(trigger_paths, dict):
            single_muon_paths = trigger_paths.get("single_muon", [])
            double_muon_paths = trigger_paths.get("double_muon", [])
        else:
            single_muon_paths = [hlt for hlt in trigger_paths if "IsoMu" in hlt]
            double_muon_paths = [hlt for hlt in trigger_paths if hlt not in single_muon_paths]

        single_trigger_cut = self._evaluate_trigger_paths(events, single_muon_paths)
        double_trigger_cut = self._evaluate_trigger_paths(events, double_muon_paths)
        trigger_cut = single_trigger_cut | double_trigger_cut
        return single_trigger_cut, double_trigger_cut, trigger_cut

    def _evaluate_trigger_paths(self, events, trigger_paths):
        trigger_cut = ak.zeros_like(events.run, dtype=bool)

        for hlt in trigger_paths:
            if hlt in events.fields:
                trigger_cut = trigger_cut | (events[hlt] == True)
            else:
                logger.debug("[TnPZmmgTagger] %s is not available in these events.", hlt)

        return trigger_cut
