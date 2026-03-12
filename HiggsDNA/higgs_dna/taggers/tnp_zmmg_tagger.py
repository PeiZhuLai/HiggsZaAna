import copy
import logging
import numpy
import awkward as ak
import vector

vector.register_awkward()

logger = logging.getLogger(__name__)

from higgs_dna.taggers.tagger import Tagger
from higgs_dna.taggers.za_tagger_resolved import (
    DEFAULT_OPTIONS as ZA_DEFAULT_OPTIONS,
    ZaTaggerRun3,
)
from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.selections import lepton_selections, object_selections

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
    "pfRelIso03_chg",
    "pfRelIso03_all",
    "miniPFRelIso_all",
    "dxy",
    "dz",
    "sip3d",
    "tightId",
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
    "pass_phid_custom_extend_tight",
    "pass_phid_sieie_tight",
    "pass_phid_PFECalIso_tight",
    "pass_phid_official_tight",
    "pass_phid_custom_medium",
    "pass_phid_custom_extend_medium",
    "pass_phid_sieie_medium",
    "pass_phid_PFECalIso_medium",
    "pass_phid_official_medium",
    "pass_phid_custom_loose",
    "pass_phid_custom_extend_loose",
    "pass_phid_sieie_loose",
    "pass_phid_PFECalIso_loose",
    "pass_phid_official_loose",
]

P4_FIELDS = ["pt", "eta", "phi", "mass"]

DEFAULT_OPTIONS = {
    "dimuon": {
        "lead_muon_pt": 20.0,
    },
    "photons": copy.deepcopy(ZA_DEFAULT_OPTIONS["photons"]),
    "electrons": {
        "pt": 7.0,
    },
    "fsr": {
        "min_muon_gamma_dr": 0.8,
        "near_muon_gamma_dr_min": 0.4,
        "far_muon_pt": 20.0,
        "dimuon_mass_min": 35.0,
        "photon_mva_min": -0.7,
        "zmmg_mass": [60.0, 120.0],
        "mass_sum_max": 180.0,
    },
    "muons": {
        "pt": 10.0,
        "eta": 2.4,
        "pfRelIso03_chg": 0.2,
    },
    "trigger": [
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8",
    ],
}


class TnPZmmgTagger(Tagger):
    def __init__(self, name="TnPZmmgTagger", options={}, is_data=None, year=None):
        super(TnPZmmgTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                original=DEFAULT_OPTIONS,
                new=options,
            )

    def calculate_selection(self, events):
        year = self.get_analysis_year(events)
        rho = self.get_rho(events)
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

        electrons = self.select_electrons_for_photon_id(events, year)
        awkward_utils.add_field(events, "SelectedElectron", electrons, overwrite=True)
        awkward_utils.add_field(events, "n_electrons", ak.num(electrons), overwrite=True)

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
            electrons=electrons,
            rho=rho,
            year=year,
        )
        awkward_utils.add_field(events, "Photon", photons_with_flags, overwrite=True)

        photons = events.Photon[photon_selection]
        # Keep photons close to muons for FSR; only retain the electron cleaning.
        clean_photon_mask = ak.fill_none(object_selections.delta_R(photons, electrons, 0.3), True)
        photons = photons[clean_photon_mask]
        photon_idx = ak.local_index(events.Photon.pt, axis=1)[photon_selection][clean_photon_mask]
        photons = ak.with_field(photons, photon_idx, "origIndex")
        photons = photons[ak.argsort(photons.pt, ascending=False, axis=1)]
        if "mass" not in photons.fields:
            photons = ak.with_field(photons, ak.zeros_like(photons.pt), "mass")
        photons = ak.Array(photons, with_name="Momentum4D")
        awkward_utils.add_field(events, "SelectedPhoton", photons, overwrite=True)

        awkward_utils.add_field(events, "n_muons", ak.num(muons), overwrite=True)
        awkward_utils.add_field(events, "n_photons", ak.num(photons), overwrite=True)
        awkward_utils.add_field(
            events,
            "n_leptons",
            ak.num(electrons) + ak.num(muons),
            overwrite=True,
        )

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
        lead_muon_pt_cut = dimuon_lead_pt > self.options["dimuon"]["lead_muon_pt"]
        dimuon_candidates = dimuon_candidates[lead_muon_pt_cut]
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
        far_pt_cut = zmmg_candidates_all.FarMuon.pt > self.options["fsr"]["far_muon_pt"]
        dimuon_mass_cut = zmmg_candidates_all.Dimuon.mass > self.options["fsr"]["dimuon_mass_min"]
        photon_mva_cut = zmmg_candidates_all.ZPhoton.mvaID > self.options["fsr"]["photon_mva_min"]
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
            & far_pt_cut
            & dimuon_mass_cut
            & photon_mva_cut
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

        trigger_cut = self.passing_dimuon_trigger(events)
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
        far_pt_event_cut = photon_cut & (
            ak.num(zmmg_candidates_all[near_dr_max_cut & near_dr_min_cut & far_pt_cut]) >= 1
        )
        dimuon_mass_event_cut = photon_cut & (
            ak.num(
                zmmg_candidates_all[
                    near_dr_max_cut & near_dr_min_cut & far_pt_cut & dimuon_mass_cut
                ]
            ) >= 1
        )
        photon_mva_event_cut = photon_cut & (
            ak.num(
                zmmg_candidates_all[
                    near_dr_max_cut
                    & near_dr_min_cut
                    & far_pt_cut
                    & dimuon_mass_cut
                    & photon_mva_cut
                ]
            ) >= 1
        )
        zmmg_mass_event_cut = photon_cut & (
            ak.num(
                zmmg_candidates_all[
                    near_dr_max_cut
                    & near_dr_min_cut
                    & far_pt_cut
                    & dimuon_mass_cut
                    & photon_mva_cut
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

        probe_photon_e_veto = best_candidate.ProbePhoton.electronVeto > self.options["photons"]["e_veto"]
        probe_photon_pass_pixel_veto = best_candidate.ProbePhoton.pixelSeed < 0.5
        awkward_utils.add_field(
            events,
            "probe_photon_e_veto",
            ak.fill_none(probe_photon_e_veto, False),
            overwrite=True,
        )
        awkward_utils.add_field(
            events,
            "probe_photon_pass_csev",
            ak.fill_none(probe_photon_e_veto, False),
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
                "OSSF dimuon with lead mu pT > 20",
                "photon",
                "min dR(mu,gamma) < 0.8",
                "dR(gamma, near muon) > 0.4",
                "pT_far > 20",
                "m_mumu > 35",
                "photon ID MVA > -0.7",
                "60 < m_mumugamma < 120",
                "m_mumu + m_mumugamma < 180",
                "all",
            ],
            results=[
                trigger_cut,
                dimuon_cut,
                photon_cut,
                near_dr_max_event_cut,
                near_dr_min_event_cut,
                far_pt_event_cut,
                dimuon_mass_event_cut,
                photon_mva_event_cut,
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

    def select_electrons_for_photon_id(self, events, year):
        electron_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electrons"],
            clean={},
            name="SelectedElectron",
            tagger=self,
            year=year,
        )

        electrons = events.Electron[electron_cut]
        electron_idx = ak.local_index(events.Electron.pt, axis=1)[electron_cut]
        electron_idx = ak.mask(electron_idx, ak.num(electron_idx) > 0)
        electrons = ak.with_field(electrons, electron_idx, "Idx")
        return electrons

    def select_muons(self, muons, options):
        pt_cut = muons.pt > options["pt"]
        eta_cut = numpy.abs(muons.eta) < options["eta"]
        id_cut = muons.tightId == True
        iso_cut = muons.pfRelIso03_chg < options["pfRelIso03_chg"]

        all_cuts = pt_cut & eta_cut & id_cut & iso_cut

        self.register_cuts(
            names=["pt", "eta", "tight ID", "pfRelIso03_chg", "all"],
            results=[pt_cut, eta_cut, id_cut, iso_cut, all_cuts],
            cut_type="muon",
        )

        return all_cuts

    def select_photons(self, photons, electrons, rho, year):
        return ZaTaggerRun3.select_photons(
            self,
            photons=photons,
            options=self.options["photons"],
            electrons=electrons,
            rho=rho,
            year=year,
        )

    def choose_object(self, condition, first, second):
        fields = first.fields
        data = {
            field: ak.where(condition, first[field], second[field])
            for field in fields
        }
        return ak.zip(data, with_name="Momentum4D")

    def add_flat_object_fields(self, events, name, objects, fields):
        available_fields = [field for field in fields if field in objects.fields]
        if not available_fields:
            return

        awkward_utils.add_object_fields(
            events=events,
            name=name,
            objects=ak.singletons(objects),
            n_objects=1,
            dummy_value=DUMMY_VALUE,
            fields=available_fields,
            overwrite=True,
        )

    def passing_dimuon_trigger(self, events):
        trigger_cut = ak.zeros_like(events.run, dtype=bool)

        trigger_paths = self.options["trigger"]
        if isinstance(trigger_paths, dict):
            year = self.get_analysis_year(events)
            trigger_paths = trigger_paths.get(year, trigger_paths.get("default", []))

        for hlt in trigger_paths:
            if hlt in events.fields:
                trigger_cut = trigger_cut | (events[hlt] == True)
            else:
                logger.debug("[TnPZmmgTagger] %s is not available in these events.", hlt)

        return trigger_cut
