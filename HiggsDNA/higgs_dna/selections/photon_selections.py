import awkward as ak
import numpy
import vector

vector.register_awkward()

from higgs_dna.selections import object_selections

"""
Resolved H->Za photon selections.

從 higgs_dna/taggers/za_tagger_resolved.py 抽出的 photon selection 邏輯，
內容與原 ZaTaggerRun3.select_photons / select_FSRphotons 逐字一致（cut values、
變數命名、array 操作、回傳型態皆未更動）。tagger 端保留同名 method 作為
thin wrapper 委派到這裡，維持既有 interface 與子類相容性。
"""


def to_momentum4d(obj):
    out = ak.zip({
        "pt": obj.pt,
        "eta": obj.eta,
        "phi": obj.phi,
        "mass": obj.mass if "mass" in obj.fields else ak.zeros_like(obj.pt),
    }, with_name="Momentum4D")
    out = ak.Array(out.layout) 
    return out


def select_resolved_photons(photons, options, electrons, rho, year, name="none", tagger=None):
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
    tagger_name = "none" if tagger is None else tagger.name

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

    # NEW: per-photon ID sub-cut flags used only by the fine a-candidate cutflow
    # decomposition (see za_tagger_resolved.produce_and_select_zgammas). Default to
    # all-True (matches the Run2 dummy ID, which is all-true); the Run3 branch below
    # overrides them with the real H/E, PFchIso, PFHCalIso sub-cuts. This does NOT
    # change any selection: the nominal id_cut is still phid_custom_tight.
    ph_id_hoe_flag = no_cut
    ph_id_chiso_flag = no_cut
    ph_id_hcaliso_flag = no_cut

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
        # NEW: capture the 3 nominal ID sub-cuts (their AND == phid_custom_tight == id_cut)
        ph_id_hoe_flag = hoe_cut
        ph_id_chiso_flag = PFChIso_cut
        ph_id_hcaliso_flag = PFHCalIso_cut
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

    # NEW: per-photon boolean flags for the fine a-candidate cutflow decomposition.
    # These only ADD fields to the returned photons collection; they do not alter any
    # cut value, mask, ordering, or the returned all_cuts selection. Note that
    #   pass_ph_id_hoe & pass_ph_id_chiso & pass_ph_id_hcaliso == id_cut (phid_custom_tight)
    # so an upstream cumulative "N_gamma>=2" count using these reproduces the nominal
    # photon selection exactly (see za_tagger_resolved.produce_and_select_zgammas).
    photons = ak.with_field(photons, pt_cut, "pass_ph_pt")
    photons = ak.with_field(photons, eta_cut, "pass_ph_eta")
    photons = ak.with_field(photons, ph_id_hoe_flag, "pass_ph_id_hoe")
    photons = ak.with_field(photons, ph_id_chiso_flag, "pass_ph_id_chiso")
    photons = ak.with_field(photons, ph_id_hcaliso_flag, "pass_ph_id_hcaliso")
    photons = ak.with_field(photons, e_veto_cut, "pass_ph_eveto")
    photons = ak.with_field(photons, eg_overlap_cut, "pass_ph_eg_overlap")

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

def select_resolved_fsr_photons(FSRphotons, electrons, muons, photons, options, name="none", tagger=None):

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
