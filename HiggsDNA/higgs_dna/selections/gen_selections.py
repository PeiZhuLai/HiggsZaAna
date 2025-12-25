import awkward
import vector

vector.register_awkward()

from higgs_dna.utils import awkward_utils

def select_x_to_yz(gen_part, x_pdgId, y_pdgId, z_pdgId):
    """
    Return all x->yz decays, sorted by x_pt.
    Values of `None` for any of `x_pdgId`, `y_pdgId`, or `z_pdgId` will result in no requirement on pdgId.
    """
    if not isinstance(gen_part, vector.Vector4D):
        gen_part = awkward.Array(gen_part, with_name = "Momentum4D")

    if x_pdgId is None:
        x_pdgId = abs(gen_part.pdgId)
    if y_pdgId is None:
        y_pdgId = abs(gen_part.pdgId)
    if z_pdgId is None:
        z_pdgId = abs(gen_part.pdgId)

    # ------------------------------------------------------------
    # NEW: 先做「標準兩體衰變」選擇：child(y/z) 的 mother 是 x
    # 這條路徑會涵蓋 ALP->γγ 這種沒有 "child child" 的情況
    # ------------------------------------------------------------
    child_is_yz = (abs(gen_part.pdgId) == y_pdgId) | (abs(gen_part.pdgId) == z_pdgId)
    mother_is_x = abs(gen_part.pdgId[gen_part.genPartIdxMother]) == x_pdgId
    gen_yz_from_x_direct = gen_part[child_is_yz & mother_is_x]
    gen_yz_from_x_direct = gen_yz_from_x_direct[
        awkward.argsort(abs(gen_yz_from_x_direct.pdgId), ascending=False, axis=1)
    ]

    # ------------------------------------------------------------
    # 兼容舊邏輯：處理 x->(y,z) 且 y/z 再衰變（例如 Z->ll）的特殊輸出格式
    # 只有在確實找到「孫代 lepton」時才走這條路徑
    # ------------------------------------------------------------
    gen_yz_from_x_legacy = gen_part[
        (
            ((abs(gen_part.pdgId) == y_pdgId) | (abs(gen_part.pdgId) == z_pdgId))
            & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == x_pdgId)
        )
        | (
            (
                ((abs(gen_part.pdgId) == 11) | (abs(gen_part.pdgId) == 13) | (abs(gen_part.pdgId) == 15))
                & (abs(gen_part.pdgId[gen_part.genPartIdxMother]) == y_pdgId)
                & (abs(gen_part.pdgId[gen_part.genPartIdxMother[gen_part.genPartIdxMother]]) == x_pdgId)
            )
        )
    ]
    gen_yz_from_x_legacy = gen_yz_from_x_legacy[
        awkward.argsort(abs(gen_yz_from_x_legacy.pdgId), ascending=False, axis=1)
    ]

    # FIX: use event-level mask (not scalar bool) for switching
    use_legacy_evt = awkward.num(gen_yz_from_x_legacy, axis=1) >= 4

    # -------------------------
    # legacy: 4-body 輸出（維持原本欄位）
    # NOTE: even when <4 objects, awkward.combinations will yield empty per-event lists (good)
    # -------------------------
    gen_child_pairs_legacy = awkward.combinations(
        gen_yz_from_x_legacy, 4,
        fields=["LeadGenChild", "SubleadGenChild", "LeadGenChildChild1", "LeadGenChildChild2"]
    )
    gen_child_pairs_legacy = gen_child_pairs_legacy[
        gen_child_pairs_legacy.LeadGenChild.genPartIdxMother == gen_child_pairs_legacy.SubleadGenChild.genPartIdxMother
    ]
    gen_child_pairs_legacy["GenParent"] = gen_part[gen_child_pairs_legacy.LeadGenChild.genPartIdxMother]
    gen_child_pairs_legacy[("GenParent", "dR")] = gen_child_pairs_legacy.LeadGenChild.deltaR(gen_child_pairs_legacy.SubleadGenChild)
    gen_child_pairs_legacy[("GenParent", "Child_1_Id")] = abs(gen_child_pairs_legacy.LeadGenChild.pdgId)
    gen_child_pairs_legacy[("GenParent", "Child_2_Id")] = abs(gen_child_pairs_legacy.SubleadGenChild.pdgId)
    if awkward.any(awkward.num(gen_child_pairs_legacy) >= 2):
        gen_child_pairs_legacy = gen_child_pairs_legacy[
            awkward.argsort(gen_child_pairs_legacy.GenParent.pt, ascending=False, axis=1)
        ]

    # -------------------------
    # direct: 2-body 輸出（補齊 legacy 期望的欄位：ChildChild1/2 給空）
    # -------------------------
    gen_child_pairs_direct = awkward.combinations(
        gen_yz_from_x_direct, 2, fields=["LeadGenChild", "SubleadGenChild"]
    )
    gen_child_pairs_direct = gen_child_pairs_direct[
        gen_child_pairs_direct.LeadGenChild.genPartIdxMother == gen_child_pairs_direct.SubleadGenChild.genPartIdxMother
    ]
    gen_child_pairs_direct["GenParent"] = gen_part[gen_child_pairs_direct.LeadGenChild.genPartIdxMother]
    gen_child_pairs_direct[("GenParent", "dR")] = gen_child_pairs_direct.LeadGenChild.deltaR(gen_child_pairs_direct.SubleadGenChild)
    gen_child_pairs_direct[("GenParent", "Child_1_Id")] = abs(gen_child_pairs_direct.LeadGenChild.pdgId)
    gen_child_pairs_direct[("GenParent", "Child_2_Id")] = abs(gen_child_pairs_direct.SubleadGenChild.pdgId)

    # 兼容：補上 legacy 會用到的欄位（對 2-body 沒有意義 -> 空/None）
    empty_like = awkward.pad_none(gen_child_pairs_direct.LeadGenChild, 0, axis=1)  # empty list per event
    gen_child_pairs_direct["LeadGenChildChild1"] = empty_like
    gen_child_pairs_direct["LeadGenChildChild2"] = empty_like

    if awkward.any(awkward.num(gen_child_pairs_direct) >= 2):
        gen_child_pairs_direct = gen_child_pairs_direct[
            awkward.argsort(gen_child_pairs_direct.GenParent.pt, ascending=False, axis=1)
        ]

    # FIX: select per-event, not scalar
    return awkward.where(use_legacy_evt, gen_child_pairs_legacy, gen_child_pairs_direct)
