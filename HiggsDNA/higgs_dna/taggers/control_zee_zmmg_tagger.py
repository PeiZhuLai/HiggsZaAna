import awkward as ak
import numpy
import vector

vector.register_awkward()

import copy
import logging
from higgs_dna.utils.logger_utils import simple_logger
logger = simple_logger(__name__)

from higgs_dna.utils import awkward_utils, misc_utils
from higgs_dna.selections import object_selections, lepton_selections
from higgs_dna.taggers.za_tagger_resolved import ZaTaggerRun3, DEFAULT_OPTIONS as ZA_DEFAULT_OPTIONS

Z_MASS = 91.1876
DUMMY_VALUE = -999.0

# 在 HZa 既有的 ZaTaggerRun3 預設選項上，再加上兩個 control region 專屬的視窗設定。
# - zee  ：Z->ee  control，用來驗證 electron reco/ID scale factor。
# - zmmg ：Z->mumu+gamma (FSR) control，用來驗證 muon SF 與 photon ID/CSEV SF。
DEFAULT_OPTIONS = copy.deepcopy(ZA_DEFAULT_OPTIONS)
DEFAULT_OPTIONS["zee"] = {
    "mass": [80.0, 100.0],   # m_ee 視窗
    "lead_pt": 25.0,
    "sublead_pt": 15.0,
}
DEFAULT_OPTIONS["zmmg"] = {
    "lead_pt": 25.0,         # 觸發 muon (single-muon HLT) 的 pT 門檻
    "sublead_pt": 4.0,       # FSR 端 muon 可以很軟
    "min_muon_gamma_dr": 0.8,    # 最近 muon 與 gamma 的 dR 必須 < 0.8 (FSR)
    "near_muon_gamma_dr_min": 0.4,  # 但不能太近 (避免重疊)
    "mmg_mass": [80.0, 100.0],   # m_mumugamma 視窗 (拉回 Z peak)
    "mass_sum_max": 180.0,       # m_mumu + m_mumugamma < 180 -> 壓非 FSR
}


class ControlZeeZmmgTagger(ZaTaggerRun3):
    """Z->ee 與 Z->mumu+gamma 合併 control-sample tagger。

    重用 ZaTaggerRun3 的 lepton / photon 選取（與正式分析一致），但事件層級
    改成挑 Z->ee 或 Z->mumugamma(FSR) 兩個 control region，任一通過即保留，
    並以 z_ee / z_mumu flag 區分。

    為了讓生產時的 scale-factor 權重「各區只作用在該區的物件」，輸出的
    SelectedElectron / SelectedMuon / SelectedPhoton 都是 region-scoped：
      - Z->ee   事件 : 只填 SelectedElectron (-> electron reco/ID SF)
      - Z->mmg  事件 : 只填 SelectedMuon + SelectedPhoton(probe) (-> muon SF, photon SF)
    其他區的集合為空 list，對應的 object SF 自動回到 1。
    """

    def __init__(self, name="control_zee_zmmg_tagger", options={}, is_data=None, year=None):
        # 繞過 ZaTaggerRun3.__init__ 用的是它自己的 DEFAULT_OPTIONS；這裡改用本檔擴充版。
        super(ZaTaggerRun3, self).__init__(name, options, is_data, year)
        if not options:
            self.options = copy.deepcopy(DEFAULT_OPTIONS)
        else:
            self.options = misc_utils.update_dict(
                original=copy.deepcopy(DEFAULT_OPTIONS),
                new=options,
            )

    # ------------------------------------------------------------------
    # 主入口
    # ------------------------------------------------------------------
    def calculate_selection(self, events):
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

        awkward_utils.add_field(events, "fixedGridRhoAll", ak.fill_none(rho, 0.0), overwrite=True)

        if not self.is_data:
            self.overlap_removal(events=events)

        selection, events = self.produce_control_regions(events, rho)
        return selection, events

    # ------------------------------------------------------------------
    # trigger helper
    # ------------------------------------------------------------------
    def _or_triggers(self, events, hlt_list):
        cut = ak.zeros_like(events.run, dtype=bool)
        for hlt in hlt_list:
            if hasattr(events, hlt):
                cut = cut | (events[hlt] == True)
            else:
                logger.debug("[ControlZeeZmmgTagger] %s not in events." % hlt)
        return cut

    # ------------------------------------------------------------------
    # control-region 選取
    # ------------------------------------------------------------------
    def produce_control_regions(self, events, rho):
        year = self.year[:4] if self.year is not None else "2022"

        # --- Electrons ---
        electron_cut = lepton_selections.select_electrons(
            electrons=events.Electron,
            options=self.options["electrons"],
            clean={},
            name="SelectedElectron",
            tagger=self,
            year=year,
        )
        electrons = awkward_utils.add_field(
            events=events, name="SelectedElectron", data=events.Electron[electron_cut]
        )
        # select_photons 內的 e/gamma overlap 需要 electrons.Idx
        arr = ak.local_index(events.Electron["pt"], axis=1)[electron_cut]
        electron_idx = ak.mask(arr, ak.num(arr) > 0)
        awkward_utils.add_field(events=electrons, name="Idx", data=electron_idx)

        # --- Muons ---
        muon_cut = lepton_selections.select_muons(
            muons=events.Muon,
            options=self.options["muons"],
            clean={},
            name="SelectedMuon",
            tagger=self,
        )
        muons = awkward_utils.add_field(
            events=events, name="SelectedMuon", data=events.Muon[muon_cut]
        )

        # --- Photons (與正式分析相同的 cut-based tight ID) ---
        photon_selection, photons_with_flags = self.select_photons(
            photons=events.Photon,
            options=self.options["photons"],
            electrons=electrons,
            rho=rho,
            year=year,
        )
        awkward_utils.add_field(events, "Photon", photons_with_flags, overwrite=True)
        photons = events.Photon[photon_selection]
        # 注意：Zmmg 的 probe 是貼著 muon 的 FSR photon，故「不」做 photon-muon 重疊移除，
        #       select_photons 內部已用 electronIdx 移除與 electron 重疊的 photon。
        photons = ak.with_field(photons, ak.zeros_like(photons.pt), "mass")

        # PDG id / ptE_error
        electrons = ak.with_field(electrons, ak.ones_like(electrons.pt) * 11, "id")
        muons = ak.with_field(muons, ak.ones_like(muons.pt) * 13, "id")
        electrons = ak.with_field(electrons, electrons.energyErr, "ptE_error")
        muons = ak.with_field(muons, muons.ptErr, "ptE_error")

        # 依 pt 排序
        electrons = electrons[ak.argsort(electrons.pt, ascending=False, axis=1)]
        muons = muons[ak.argsort(muons.pt, ascending=False, axis=1)]
        photons = photons[ak.argsort(photons.pt, ascending=False, axis=1)]

        electrons = ak.Array(electrons, with_name="Momentum4D")
        muons = ak.Array(muons, with_name="Momentum4D")
        photons = ak.Array(photons, with_name="Momentum4D")

        n_electrons = ak.fill_none(ak.num(electrons), 0)
        n_muons = ak.fill_none(ak.num(muons), 0)
        n_photons = ak.fill_none(ak.num(photons), 0)
        awkward_utils.add_field(events, "n_electrons", n_electrons, overwrite=True)
        awkward_utils.add_field(events, "n_muons", n_muons, overwrite=True)
        awkward_utils.add_field(events, "n_photons", n_photons, overwrite=True)
        awkward_utils.add_field(events, "n_leptons", n_electrons + n_muons, overwrite=True)

        if (not self.is_data) and "GenVtx_z" in events.fields and "PV_z" in events.fields:
            awkward_utils.add_field(events, "dZ", events.GenVtx_z - events.PV_z, overwrite=True)

        # triggers
        ele_trigger_cut = (
            self._or_triggers(events, self.options["single_ele_trigger"].get(year, []))
            | self._or_triggers(events, self.options["double_ele_trigger"].get(year, []))
        )
        mu_trigger_cut = (
            self._or_triggers(events, self.options["single_muon_trigger"].get(year, []))
            | self._or_triggers(events, self.options["double_muon_trigger"].get(year, []))
        )

        # ==============================================================
        # Z -> ee  control
        # ==============================================================
        ee_pairs = ak.combinations(electrons, 2, fields=["LeadLepton", "SubleadLepton"])
        os_cut = (ee_pairs.LeadLepton.charge * ee_pairs.SubleadLepton.charge) == -1
        ee_pairs = ee_pairs[os_cut]
        ee_pairs["ZCand"] = ee_pairs.LeadLepton + ee_pairs.SubleadLepton

        zee_mass_lo, zee_mass_hi = self.options["zee"]["mass"]
        ee_mass_cut = (ee_pairs.ZCand.mass > zee_mass_lo) & (ee_pairs.ZCand.mass < zee_mass_hi)
        ee_pairs = ee_pairs[ee_mass_cut]
        ee_pairs = ee_pairs[ak.argsort(abs(ee_pairs.ZCand.mass - Z_MASS), axis=1)]
        zee_cand = ak.firsts(ee_pairs)

        e_lead_pt = ak.fill_none(
            ak.where(zee_cand.LeadLepton.pt >= zee_cand.SubleadLepton.pt,
                     zee_cand.LeadLepton.pt, zee_cand.SubleadLepton.pt), 0.0)
        e_sub_pt = ak.fill_none(
            ak.where(zee_cand.LeadLepton.pt >= zee_cand.SubleadLepton.pt,
                     zee_cand.SubleadLepton.pt, zee_cand.LeadLepton.pt), 0.0)
        ee_kin_cut = (e_lead_pt > self.options["zee"]["lead_pt"]) & (e_sub_pt > self.options["zee"]["sublead_pt"])

        has_ee_cand = ak.num(ee_pairs) >= 1
        pass_zee = ak.fill_none(has_ee_cand & ee_kin_cut & ele_trigger_cut, False)

        # ==============================================================
        # Z -> mu mu gamma (FSR) control
        # ==============================================================
        mm_pairs = ak.combinations(muons, 2, fields=["LeadMuon", "SubleadMuon"])
        os_cut = (mm_pairs.LeadMuon.charge * mm_pairs.SubleadMuon.charge) == -1
        mm_pairs = mm_pairs[os_cut]
        mm_pairs["Dimuon"] = mm_pairs.LeadMuon + mm_pairs.SubleadMuon

        zmmg = ak.cartesian({"mm": mm_pairs, "ph": photons}, axis=1)
        zmmg["LeadMuon"] = zmmg.mm.LeadMuon
        zmmg["SubleadMuon"] = zmmg.mm.SubleadMuon
        zmmg["Dimuon"] = zmmg.mm.Dimuon
        zmmg["Photon"] = zmmg.ph
        zmmg["Zmmg"] = zmmg.Dimuon + zmmg.ph

        lead_dr = zmmg.LeadMuon.deltaR(zmmg.ph)
        sub_dr = zmmg.SubleadMuon.deltaR(zmmg.ph)
        min_dr = ak.where(lead_dr <= sub_dr, lead_dr, sub_dr)
        far_dr = ak.where(lead_dr <= sub_dr, sub_dr, lead_dr)
        zmmg["minMuonGammaDR"] = min_dr
        zmmg["farMuonGammaDR"] = far_dr

        opt = self.options["zmmg"]
        near_dr_max_cut = min_dr < opt["min_muon_gamma_dr"]
        near_dr_min_cut = min_dr > opt["near_muon_gamma_dr_min"]
        mmg_mass_cut = (zmmg.Zmmg.mass > opt["mmg_mass"][0]) & (zmmg.Zmmg.mass < opt["mmg_mass"][1])
        mass_sum_cut = (zmmg.Dimuon.mass + zmmg.Zmmg.mass) < opt["mass_sum_max"]
        fsr_cand_cut = near_dr_max_cut & near_dr_min_cut & mmg_mass_cut & mass_sum_cut

        zmmg = zmmg[fsr_cand_cut]
        # 同一事件多個候選 -> 取最接近 Z peak 的 m_mumugamma
        zmmg = zmmg[ak.argsort(abs(zmmg.Zmmg.mass - Z_MASS), axis=1)]
        zmmg_cand = ak.firsts(zmmg)

        m_lead_pt = ak.fill_none(
            ak.where(zmmg_cand.LeadMuon.pt >= zmmg_cand.SubleadMuon.pt,
                     zmmg_cand.LeadMuon.pt, zmmg_cand.SubleadMuon.pt), 0.0)
        m_sub_pt = ak.fill_none(
            ak.where(zmmg_cand.LeadMuon.pt >= zmmg_cand.SubleadMuon.pt,
                     zmmg_cand.SubleadMuon.pt, zmmg_cand.LeadMuon.pt), 0.0)
        mm_kin_cut = (m_lead_pt > opt["lead_pt"]) & (m_sub_pt > opt["sublead_pt"])

        has_zmmg_cand = ak.num(zmmg) >= 1
        pass_zmmg = ak.fill_none(has_zmmg_cand & mm_kin_cut & mu_trigger_cut, False)

        # ==============================================================
        # region-scoped SF 集合
        # ==============================================================
        sel_electrons = electrons[ak.broadcast_arrays(pass_zee, electrons.pt)[0]]
        awkward_utils.add_field(events, "SelectedElectron", sel_electrons, overwrite=True)

        sel_muons = muons[ak.broadcast_arrays(pass_zmmg, muons.pt)[0]]
        awkward_utils.add_field(events, "SelectedMuon", sel_muons, overwrite=True)

        probe_photon = ak.mask(zmmg_cand.Photon, pass_zmmg)
        sel_photon = ak.singletons(probe_photon)
        awkward_utils.add_field(events, "SelectedPhoton", sel_photon, overwrite=True)

        # ==============================================================
        # 輸出 branch
        # ==============================================================
        awkward_utils.add_field(events, "z_ee", pass_zee, overwrite=True)
        awkward_utils.add_field(events, "z_mumu", pass_zmmg, overwrite=True)

        # Z->ee
        for field in ["pt", "eta", "phi", "mass"]:
            awkward_utils.add_field(events, "Zee_%s" % field,
                                    ak.fill_none(getattr(zee_cand.ZCand, field), DUMMY_VALUE), overwrite=True)
        for field in ["pt", "eta", "phi", "mass", "charge"]:
            awkward_utils.add_field(events, "Zee_lead_lepton_%s" % field,
                                    ak.fill_none(zee_cand.LeadLepton[field], DUMMY_VALUE), overwrite=True)
            awkward_utils.add_field(events, "Zee_sublead_lepton_%s" % field,
                                    ak.fill_none(zee_cand.SubleadLepton[field], DUMMY_VALUE), overwrite=True)

        # Z->mumugamma
        for field in ["pt", "eta", "phi", "mass"]:
            awkward_utils.add_field(events, "Zmmg_%s" % field,
                                    ak.fill_none(getattr(zmmg_cand.Zmmg, field), DUMMY_VALUE), overwrite=True)
            awkward_utils.add_field(events, "dimuon_%s" % field,
                                    ak.fill_none(getattr(zmmg_cand.Dimuon, field), DUMMY_VALUE), overwrite=True)
        for field in ["pt", "eta", "phi", "mass", "charge"]:
            awkward_utils.add_field(events, "zmmg_lead_muon_%s" % field,
                                    ak.fill_none(zmmg_cand.LeadMuon[field], DUMMY_VALUE), overwrite=True)
            awkward_utils.add_field(events, "zmmg_sublead_muon_%s" % field,
                                    ak.fill_none(zmmg_cand.SubleadMuon[field], DUMMY_VALUE), overwrite=True)
        for field in ["pt", "eta", "phi", "mass", "mvaID", "r9"]:
            awkward_utils.add_field(events, "probe_photon_%s" % field,
                                    ak.fill_none(zmmg_cand.Photon[field], DUMMY_VALUE), overwrite=True)
        awkward_utils.add_field(events, "probe_photon_lep_near_dR",
                                ak.fill_none(zmmg_cand.minMuonGammaDR, DUMMY_VALUE), overwrite=True)
        awkward_utils.add_field(events, "probe_photon_lep_far_dR",
                                ak.fill_none(zmmg_cand.farMuonGammaDR, DUMMY_VALUE), overwrite=True)
        awkward_utils.add_field(events, "zmmg_minus_dimuon_mass",
                                ak.fill_none(zmmg_cand.Zmmg.mass - zmmg_cand.Dimuon.mass, DUMMY_VALUE), overwrite=True)

        # 方便畫圖的合併 Z 變數（ee 取 m_ee，mmg 取 m_mumugamma）
        zee_mass = ak.fill_none(zee_cand.ZCand.mass, DUMMY_VALUE)
        zmmg_mass = ak.fill_none(zmmg_cand.Zmmg.mass, DUMMY_VALUE)
        zee_pt = ak.fill_none(zee_cand.ZCand.pt, DUMMY_VALUE)
        zmmg_pt = ak.fill_none(zmmg_cand.Zmmg.pt, DUMMY_VALUE)
        awkward_utils.add_field(events, "Z_mass",
                                ak.where(pass_zee, zee_mass, ak.where(pass_zmmg, zmmg_mass, DUMMY_VALUE)), overwrite=True)
        awkward_utils.add_field(events, "Z_pt",
                                ak.where(pass_zee, zee_pt, ak.where(pass_zmmg, zmmg_pt, DUMMY_VALUE)), overwrite=True)

        # ==============================================================
        # cutflow
        # ==============================================================
        presel_cut = pass_zee | pass_zmmg
        self.register_cuts(
            names=["Z->ee region", "Z->mumugamma region", "all"],
            results=[pass_zee, pass_zmmg, presel_cut],
            cut_type="event",
        )

        return presel_cut, events
