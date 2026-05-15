#!/usr/bin/env python3
"""
Angular distributions for H -> Z a -> l+l- gamma gamma.

Physics summary:
  * cos(theta_Z) probes the helicity/polarization of the spin-1 Z boson.
    A dominantly longitudinal Z gives an approximate 1 - cos^2(theta_Z)
    shape before detector acceptance and analysis selections.
  * cos(theta_a) checks the spin-0 ALP decay.  At generator level,
    a -> gamma gamma should be isotropic in the ALP rest frame, so
    cos(theta_a) should be approximately flat.
  * The angle phi between the Z and ALP decay planes is a more CP-sensitive
    observable than cos(theta_a) alone.
  * Non-flat reconstructed-level shapes should not be interpreted as
    intrinsic spin/helicity effects without comparison to generator-level
    shapes, because photon/lepton acceptance, ID, pT thresholds, and photon
    ordering can sculpt them.

The reusable functions below use numpy four-vectors with fields
px, py, pz, E.  They intentionally do not depend on ROOT so the same code can
run on HiggsDNA parquet, converted ROOT, or nano-like ROOT files.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import awkward as ak

_CACHE_DIR = "/private/tmp/hza_angular_matplotlib_cache"
os.makedirs(_CACHE_DIR, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", _CACHE_DIR)
os.environ.setdefault("XDG_CACHE_HOME", _CACHE_DIR)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DUMMY_VALUE = -999.0
DEFAULT_MASSES = (1, 5, 10, 20, 30)
DEFAULT_INPUT_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal"


Vec = Dict[str, np.ndarray]


def _as_np(x) -> np.ndarray:
    if isinstance(x, np.ndarray):
        return x
    try:
        return ak.to_numpy(ak.flatten(x, axis=None))
    except Exception:
        return np.asarray(x)


def make_p4_from_pt_eta_phi_m(pt, eta, phi, mass) -> Vec:
    pt = _as_np(pt).astype(float)
    eta = _as_np(eta).astype(float)
    phi = _as_np(phi).astype(float)
    mass = _as_np(mass).astype(float)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(np.maximum(px * px + py * py + pz * pz + mass * mass, 0.0))
    return {"px": px, "py": py, "pz": pz, "E": energy}


def sanitize_pt_eta_phi_m(pt, eta, phi, mass_value) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pt = _as_np(pt).astype(float)
    eta = _as_np(eta).astype(float)
    phi = _as_np(phi).astype(float)
    mass_value = _as_np(mass_value).astype(float)
    invalid = (
        ~np.isfinite(pt)
        | ~np.isfinite(eta)
        | ~np.isfinite(phi)
        | ~np.isfinite(mass_value)
        | (pt <= 0.0)
        | np.isclose(pt, DUMMY_VALUE)
        | np.isclose(eta, DUMMY_VALUE)
        | np.isclose(phi, DUMMY_VALUE)
        | (np.abs(eta) > 20.0)
    )
    pt = np.where(invalid, np.nan, pt)
    eta = np.where(invalid, np.nan, eta)
    phi = np.where(invalid, np.nan, phi)
    mass_value = np.where(invalid, np.nan, mass_value)
    return pt, eta, phi, mass_value


def add_p4(a: Vec, b: Vec) -> Vec:
    return {k: a[k] + b[k] for k in ("px", "py", "pz", "E")}


def spatial(vec: Vec) -> np.ndarray:
    return np.stack([vec["px"], vec["py"], vec["pz"]], axis=-1)


def mass(vec: Vec) -> np.ndarray:
    m2 = vec["E"] * vec["E"] - vec["px"] * vec["px"] - vec["py"] * vec["py"] - vec["pz"] * vec["pz"]
    return np.sqrt(np.maximum(m2, 0.0))


def finite_p4_mask(*vecs: Vec) -> np.ndarray:
    if not vecs:
        return np.asarray([], dtype=bool)
    mask = np.ones_like(vecs[0]["E"], dtype=bool)
    for vec in vecs:
        for key in ("px", "py", "pz", "E"):
            mask &= np.isfinite(vec[key])
        mask &= vec["E"] > 0.0
    return mask


def boost_to_rest_frame(vec: Vec, parent: Vec) -> Vec:
    """Boost vec into the rest frame of parent."""
    beta = -spatial(parent) / parent["E"][:, None]
    beta2 = np.sum(beta * beta, axis=1)
    beta2_safe = np.where(beta2 > 0.0, beta2, 1.0)
    gamma = np.ones_like(beta2)
    good = (beta2 > 0.0) & (beta2 < 1.0)
    gamma[good] = 1.0 / np.sqrt(1.0 - beta2[good])

    p = spatial(vec)
    bp = np.sum(beta * p, axis=1)
    factor = np.zeros_like(bp)
    factor[good] = ((gamma[good] - 1.0) * bp[good] / beta2_safe[good]) + gamma[good] * vec["E"][good]

    p_prime = p + factor[:, None] * beta
    e_prime = gamma * (vec["E"] + bp)
    invalid = beta2 >= 1.0
    if np.any(invalid):
        p_prime[invalid, :] = np.nan
        e_prime[invalid] = np.nan
    return {"px": p_prime[:, 0], "py": p_prime[:, 1], "pz": p_prime[:, 2], "E": e_prime}


def unit_vector(vec3: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vec3, axis=-1)
    out = np.full_like(vec3, np.nan, dtype=float)
    good = norm > 0.0
    out[good] = vec3[good] / norm[good, None]
    return out


def cos_angle(vec3_a: np.ndarray, vec3_b: np.ndarray) -> np.ndarray:
    ua = unit_vector(vec3_a)
    ub = unit_vector(vec3_b)
    cos_value = np.sum(ua * ub, axis=-1)
    return np.clip(cos_value, -1.0, 1.0)


def compute_cos_theta_Z(lep_minus: Vec, lep_plus: Vec, gamma1: Vec, gamma2: Vec) -> np.ndarray:
    z_vec = add_p4(lep_minus, lep_plus)
    a_vec = add_p4(gamma1, gamma2)
    h_vec = add_p4(z_vec, a_vec)
    z_h = boost_to_rest_frame(z_vec, h_vec)
    lep_minus_h = boost_to_rest_frame(lep_minus, h_vec)
    lep_minus_z = boost_to_rest_frame(lep_minus_h, z_h)
    return cos_angle(spatial(lep_minus_z), spatial(z_h))


def compute_cos_theta_a(lep_minus: Vec, lep_plus: Vec, gamma1: Vec, gamma2: Vec) -> np.ndarray:
    z_vec = add_p4(lep_minus, lep_plus)
    a_vec = add_p4(gamma1, gamma2)
    h_vec = add_p4(z_vec, a_vec)
    a_h = boost_to_rest_frame(a_vec, h_vec)
    gamma1_h = boost_to_rest_frame(gamma1, h_vec)
    gamma1_a = boost_to_rest_frame(gamma1_h, a_h)
    return cos_angle(spatial(gamma1_a), spatial(a_h))


def compute_decay_plane_phi(
    lep_minus: Vec,
    lep_plus: Vec,
    gamma1: Vec,
    gamma2: Vec,
    signed: bool = False,
) -> np.ndarray:
    z_vec = add_p4(lep_minus, lep_plus)
    a_vec = add_p4(gamma1, gamma2)
    h_vec = add_p4(z_vec, a_vec)

    lep_minus_h = boost_to_rest_frame(lep_minus, h_vec)
    lep_plus_h = boost_to_rest_frame(lep_plus, h_vec)
    gamma1_h = boost_to_rest_frame(gamma1, h_vec)
    gamma2_h = boost_to_rest_frame(gamma2, h_vec)
    z_h = boost_to_rest_frame(z_vec, h_vec)

    n_z = unit_vector(np.cross(spatial(lep_minus_h), spatial(lep_plus_h)))
    n_a = unit_vector(np.cross(spatial(gamma1_h), spatial(gamma2_h)))
    cos_phi = np.clip(np.sum(n_z * n_a, axis=-1), -1.0, 1.0)
    phi = np.arccos(cos_phi)

    if signed:
        sign_value = np.sign(np.sum(spatial(z_h) * np.cross(n_z, n_a), axis=-1))
        sign_value = np.where(sign_value == 0.0, 1.0, sign_value)
        phi = sign_value * phi
    return phi


def compute_cosTheta_Z(lep_minus: Vec, lep_plus: Vec, gamma1: Vec, gamma2: Vec) -> np.ndarray:
    z_vec = add_p4(lep_minus, lep_plus)
    a_vec = add_p4(gamma1, gamma2)
    h_vec = add_p4(z_vec, a_vec)
    z_h = boost_to_rest_frame(z_vec, h_vec)
    beam = {
        "px": np.zeros_like(h_vec["E"]),
        "py": np.zeros_like(h_vec["E"]),
        "pz": np.ones_like(h_vec["E"]),
        "E": np.ones_like(h_vec["E"]),
    }
    beam_h = boost_to_rest_frame(beam, h_vec)
    return cos_angle(spatial(z_h), spatial(beam_h))


def compute_cosTheta_a(lep_minus: Vec, lep_plus: Vec, gamma1: Vec, gamma2: Vec) -> np.ndarray:
    z_vec = add_p4(lep_minus, lep_plus)
    a_vec = add_p4(gamma1, gamma2)
    h_vec = add_p4(z_vec, a_vec)
    a_h = boost_to_rest_frame(a_vec, h_vec)
    beam = {
        "px": np.zeros_like(h_vec["E"]),
        "py": np.zeros_like(h_vec["E"]),
        "pz": np.ones_like(h_vec["E"]),
        "E": np.ones_like(h_vec["E"]),
    }
    beam_h = boost_to_rest_frame(beam, h_vec)
    return cos_angle(spatial(a_h), spatial(beam_h))


def _valid_values(values: np.ndarray, xmin: float, xmax: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values) & (values >= xmin) & (values <= xmax)]


@dataclass
class SampleArrays:
    label: str
    mass: Optional[int]
    arrays: ak.Array
    source_files: List[Path]


def infer_mass_from_path(path: Path) -> Optional[int]:
    text = str(path)
    patterns = (
        r"mA[_-]?M(\d+)",
        r"mA[_-]?(\d+)",
        r"ALP[_-]?M(\d+)",
        r"ALP[_-]?(\d+)",
        r"\bM(\d+)\b",
    )
    for pattern in patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def discover_files(inputs: Sequence[str], masses: Sequence[int], years: Sequence[str]) -> Dict[str, List[Path]]:
    files: List[Tuple[Path, bool]] = []
    for item in inputs:
        path = Path(item).expanduser()
        if path.is_file() and path.suffix.lower() in {".root", ".parquet"}:
            files.append((path, True))
        elif path.is_dir():
            files.extend((p, False) for p in sorted(path.rglob("*")) if p.suffix.lower() in {".root", ".parquet"})

    selected: Dict[str, List[Path]] = {}
    for path, is_explicit_file in files:
        mass_value = infer_mass_from_path(path)
        if masses and mass_value not in masses and not (is_explicit_file and mass_value is None):
            continue
        if years and not any(year in str(path) for year in years):
            continue
        label = f"m_{{a}} = {mass_value} GeV" if mass_value is not None else path.stem
        selected.setdefault(label, []).append(path)
    return dict(sorted(selected.items(), key=lambda kv: (infer_mass_from_path(kv[1][0]) is None, infer_mass_from_path(kv[1][0]) or 9999, kv[0])))


def _root_tree(file_obj, preferred: str):
    if preferred in file_obj:
        return file_obj[preferred]
    for key in ("inclusive", "Events", "tree"):
        if key in file_obj:
            return file_obj[key]
    for key in file_obj.keys():
        obj = file_obj[key]
        if hasattr(obj, "arrays"):
            return obj
    raise RuntimeError("No TTree found")


def read_files(paths: Sequence[Path], tree_name: str) -> ak.Array:
    chunks = []
    root_paths = [p for p in paths if p.suffix.lower() == ".root"]
    parquet_paths = [p for p in paths if p.suffix.lower() == ".parquet"]

    if root_paths:
        import uproot

        for path in root_paths:
            with uproot.open(path) as root_file:
                tree = _root_tree(root_file, tree_name)
                chunks.append(tree.arrays(library="ak"))

    for path in parquet_paths:
        try:
            chunks.append(ak.from_parquet(path))
        except AttributeError as exc:
            if "PyExtensionType" not in str(exc):
                raise
            import pyarrow.parquet as pq

            table = pq.read_table(path)
            chunks.append(ak.Array(table.to_pydict()))

    if not chunks:
        return ak.Array([])
    if len(chunks) == 1:
        return chunks[0]
    return ak.concatenate(chunks, axis=0)


def has_fields(arrays: ak.Array, fields: Sequence[str]) -> bool:
    available = set(arrays.fields)
    return all(field in available for field in fields)


def flat_branch_p4(arrays: ak.Array, prefix: str, default_mass: float = 0.0) -> Optional[Vec]:
    fields = arrays.fields
    need = [f"{prefix}_pt", f"{prefix}_eta", f"{prefix}_phi"]
    if not all(field in fields for field in need):
        return None
    mass_field = f"{prefix}_mass"
    pt = _as_np(arrays[f"{prefix}_pt"])
    eta = _as_np(arrays[f"{prefix}_eta"])
    phi = _as_np(arrays[f"{prefix}_phi"])
    if mass_field in fields:
        m = _as_np(arrays[mass_field])
    else:
        m = np.zeros_like(pt) + default_mass
    pt, eta, phi, m = sanitize_pt_eta_phi_m(pt, eta, phi, m)
    return make_p4_from_pt_eta_phi_m(pt, eta, phi, m)


def vectors_from_flat_gen(arrays: ak.Array) -> Optional[Tuple[Vec, Vec, Vec, Vec]]:
    lep_minus = flat_branch_p4(arrays, "GenHzaZLepMinus", 0.0)
    lep_plus = flat_branch_p4(arrays, "GenHzaZLepPlus", 0.0)
    if lep_minus is None or lep_plus is None:
        lep1 = flat_branch_p4(arrays, "GenHzaZLeadLep", 0.0)
        lep2 = flat_branch_p4(arrays, "GenHzaZSubleadLep", 0.0)
        if lep1 is None or lep2 is None:
            return None
        if "GenHzaZLeadLep_pdgId" in arrays.fields and "GenHzaZSubleadLep_pdgId" in arrays.fields:
            pdg1 = _as_np(arrays["GenHzaZLeadLep_pdgId"])
            pdg2 = _as_np(arrays["GenHzaZSubleadLep_pdgId"])
            lep_minus = _select_vec(np.asarray(pdg1) > 0, lep1, lep2)
            lep_plus = _select_vec(np.asarray(pdg1) > 0, lep2, lep1)
            missing_charge = ~np.isfinite(pdg1) | ~np.isfinite(pdg2)
            if np.any(missing_charge):
                lep_minus = _select_vec(missing_charge, lep1, lep_minus)
                lep_plus = _select_vec(missing_charge, lep2, lep_plus)
        else:
            lep_minus, lep_plus = lep1, lep2

    gamma1 = flat_branch_p4(arrays, "GenALPPhoton1", 0.0)
    gamma2 = flat_branch_p4(arrays, "GenALPPhoton2", 0.0)
    if gamma1 is None or gamma2 is None:
        gamma1 = flat_branch_p4(arrays, "GenALPLeadPho", 0.0)
        gamma2 = flat_branch_p4(arrays, "GenALPSubleadPho", 0.0)
    if gamma1 is None or gamma2 is None:
        return None
    return lep_minus, lep_plus, gamma1, gamma2


def _select_vec(mask: np.ndarray, a: Vec, b: Vec) -> Vec:
    return {key: np.where(mask, a[key], b[key]) for key in ("px", "py", "pz", "E")}


def vectors_from_reco(arrays: ak.Array) -> Optional[Tuple[Vec, Vec, Vec, Vec]]:
    lep1 = flat_branch_p4(arrays, "Z_lead_lepton", 0.0)
    lep2 = flat_branch_p4(arrays, "Z_sublead_lepton", 0.0)
    gamma1 = flat_branch_p4(arrays, "ALP_lead_photon", 0.0)
    gamma2 = flat_branch_p4(arrays, "ALP_sublead_photon", 0.0)
    if any(vec is None for vec in (lep1, lep2, gamma1, gamma2)):
        return None

    q1_name = "Z_lead_lepton_charge"
    q2_name = "Z_sublead_lepton_charge"
    if q1_name in arrays.fields and q2_name in arrays.fields:
        q1 = _as_np(arrays[q1_name])
        q2 = _as_np(arrays[q2_name])
        lep_minus = _select_vec(q1 < 0.0, lep1, lep2)  # type: ignore[arg-type]
        lep_plus = _select_vec(q1 < 0.0, lep2, lep1)  # type: ignore[arg-type]
        bad_charge = ~(((q1 < 0.0) & (q2 > 0.0)) | ((q2 < 0.0) & (q1 > 0.0)))
        if np.any(bad_charge):
            lep_minus = _select_vec(bad_charge, lep1, lep_minus)
            lep_plus = _select_vec(bad_charge, lep2, lep_plus)
    else:
        lep_minus, lep_plus = lep1, lep2  # type: ignore[assignment]
    return lep_minus, lep_plus, gamma1, gamma2  # type: ignore[return-value]


def _children_of(gen: ak.Array, mother_idx: ak.Array, target_pdg: int) -> ak.Array:
    local_idx = ak.local_index(gen.pt, axis=1)
    return gen[(abs(gen.pdgId) == abs(target_pdg)) & (gen.genPartIdxMother == mother_idx)]


def vectors_from_nano_genpart(arrays: ak.Array) -> Optional[Tuple[Vec, Vec, Vec, Vec]]:
    required = ("GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass", "GenPart_pdgId", "GenPart_genPartIdxMother")
    if not has_fields(arrays, required):
        return None

    gen = ak.zip(
        {
            "pt": arrays.GenPart_pt,
            "eta": arrays.GenPart_eta,
            "phi": arrays.GenPart_phi,
            "mass": arrays.GenPart_mass,
            "pdgId": arrays.GenPart_pdgId,
            "genPartIdxMother": arrays.GenPart_genPartIdxMother,
        }
    )
    idx = ak.local_index(gen.pt, axis=1)
    h_idx = idx[abs(gen.pdgId) == 25]
    h_idx = ak.firsts(h_idx)
    has_h = ~ak.is_none(h_idx)

    z = ak.firsts(gen[(abs(gen.pdgId) == 23) & (gen.genPartIdxMother == h_idx)])
    alp = ak.firsts(gen[(abs(gen.pdgId) == 9000005) & (gen.genPartIdxMother == h_idx)])
    z_idx = idx[(abs(gen.pdgId) == 23) & (gen.genPartIdxMother == h_idx)]
    alp_idx = idx[(abs(gen.pdgId) == 9000005) & (gen.genPartIdxMother == h_idx)]
    z_idx = ak.firsts(z_idx)
    alp_idx = ak.firsts(alp_idx)

    leptons = gen[((abs(gen.pdgId) == 11) | (abs(gen.pdgId) == 13) | (abs(gen.pdgId) == 15)) & (gen.genPartIdxMother == z_idx)]
    photons = gen[(abs(gen.pdgId) == 22) & (gen.genPartIdxMother == alp_idx)]
    leptons = leptons[ak.argsort(leptons.pdgId, axis=1)]
    photons = photons[ak.argsort(photons.pt, ascending=False, axis=1)]

    lep_plus = ak.firsts(leptons[leptons.pdgId < 0])
    lep_minus = ak.firsts(leptons[leptons.pdgId > 0])
    gamma1 = ak.firsts(photons)
    gamma2 = ak.firsts(photons[:, 1:])

    valid = has_h & ~ak.is_none(z) & ~ak.is_none(alp) & ~ak.is_none(lep_minus) & ~ak.is_none(lep_plus) & ~ak.is_none(gamma1) & ~ak.is_none(gamma2)
    if not bool(ak.any(valid)):
        return None

    def obj_to_p4(obj) -> Vec:
        return make_p4_from_pt_eta_phi_m(
            ak.fill_none(obj.pt, np.nan),
            ak.fill_none(obj.eta, np.nan),
            ak.fill_none(obj.phi, np.nan),
            ak.fill_none(obj.mass, 0.0),
        )

    return obj_to_p4(lep_minus), obj_to_p4(lep_plus), obj_to_p4(gamma1), obj_to_p4(gamma2)


def apply_mask_to_vecs(vecs: Tuple[Vec, Vec, Vec, Vec], mask: np.ndarray) -> Tuple[Vec, Vec, Vec, Vec]:
    return tuple({key: vec[key][mask] for key in ("px", "py", "pz", "E")} for vec in vecs)  # type: ignore[return-value]


def base_valid_mask(vecs: Tuple[Vec, Vec, Vec, Vec]) -> np.ndarray:
    lep_minus, lep_plus, gamma1, gamma2 = vecs
    z_vec = add_p4(lep_minus, lep_plus)
    a_vec = add_p4(gamma1, gamma2)
    h_vec = add_p4(z_vec, a_vec)
    mask = finite_p4_mask(lep_minus, lep_plus, gamma1, gamma2, z_vec, a_vec, h_vec)
    mask &= mass(z_vec) > 0.0
    mask &= mass(a_vec) >= 0.0
    mask &= mass(h_vec) > 0.0
    for vec in (lep_minus, lep_plus, gamma1, gamma2):
        mask &= ~np.isclose(vec["E"], DUMMY_VALUE)
    return mask


def final_selection_mask(arrays: ak.Array, mass_value: Optional[int], args: argparse.Namespace) -> Optional[np.ndarray]:
    n = len(arrays)
    if n == 0:
        return None

    if args.final_cut_expression:
        env = {field: _as_np(arrays[field]) for field in arrays.fields}
        env.update({"np": np, "abs": np.abs})
        try:
            return np.asarray(eval(args.final_cut_expression, {"__builtins__": {}}, env), dtype=bool)
        except Exception as exc:
            print(f"[final] Could not evaluate --final-cut-expression: {exc}")
            return None

    if args.final_score_branch:
        candidates = [args.final_score_branch]
    else:
        candidates = []
        if mass_value is not None:
            candidates.extend(
                [
                    f"BDT_Score_mA_M{mass_value}",
                    f"BDT_Score_mA_M{mass_value:02d}",
                    f"mvaVal_larger_M{mass_value}",
                    f"mvaVal_M{mass_value}",
                    f"mvaVal_larger_ALP_M{mass_value}",
                    f"mvaVal_ALP_M{mass_value}",
                ]
            )
        candidates.extend(("pass_final_selection", "pass_final", "final_selection"))

    for branch in candidates:
        if branch not in arrays.fields:
            continue
        values = _as_np(arrays[branch])
        if values.dtype == np.bool_:
            return values.astype(bool)
        if args.final_score_cut is not None:
            return np.asarray(values >= args.final_score_cut, dtype=bool)
        if branch.startswith(("pass_", "final_")):
            return np.asarray(values, dtype=bool)
    return None


def compute_observables(vecs: Tuple[Vec, Vec, Vec, Vec]) -> Dict[str, np.ndarray]:
    lep_minus, lep_plus, gamma1, gamma2 = vecs
    phi = compute_decay_plane_phi(lep_minus, lep_plus, gamma1, gamma2, signed=False)
    signed_phi = compute_decay_plane_phi(lep_minus, lep_plus, gamma1, gamma2, signed=True)
    return {
        "costhetaZ": compute_cos_theta_Z(lep_minus, lep_plus, gamma1, gamma2),
        "costhetaA": compute_cos_theta_a(lep_minus, lep_plus, gamma1, gamma2),
        "abscosthetaA": np.abs(compute_cos_theta_a(lep_minus, lep_plus, gamma1, gamma2)),
        "phi_decayplanes": phi,
        "cosphi_decayplanes": np.cos(phi),
        "signedphi_decayplanes": signed_phi,
        "cosThetaZ": compute_cosTheta_Z(lep_minus, lep_plus, gamma1, gamma2),
        "cosThetaA": compute_cosTheta_a(lep_minus, lep_plus, gamma1, gamma2),
    }


PLOT_SPECS = {
    "costhetaZ": (r"$\cos\theta_Z$", -1.0, 1.0, 24, "hza_costhetaZ_{level}"),
    "costhetaA": (r"$\cos\theta_a$", -1.0, 1.0, 24, "hza_costhetaA_{level}"),
    "abscosthetaA": (r"$|\cos\theta_a|$", 0.0, 1.0, 20, "hza_abscosthetaA_{level}"),
    "phi_decayplanes": (r"$\phi$ between decay planes", 0.0, math.pi, 24, "hza_phi_decayplanes_{level}"),
    "cosphi_decayplanes": (r"$\cos\phi$ between decay planes", -1.0, 1.0, 24, "hza_cosphi_decayplanes_{level}"),
    "signedphi_decayplanes": (r"signed $\phi$ between decay planes", -math.pi, math.pi, 24, "hza_signedphi_decayplanes_{level}"),
    "cosThetaZ": (r"$\cos\Theta_Z$", -1.0, 1.0, 24, "hza_cosThetaZ_{level}"),
    "cosThetaA": (r"$\cos\Theta_a$", -1.0, 1.0, 24, "hza_cosThetaA_{level}"),
}


def plot_observable(
    values_by_label: Mapping[str, np.ndarray],
    observable: str,
    level: str,
    output_dir: Path,
    title_suffix: str,
) -> None:
    xlabel, xmin, xmax, nbins, name_template = PLOT_SPECS[observable]
    bins = np.linspace(xmin, xmax, nbins + 1)

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    any_drawn = False
    for label, values in values_by_label.items():
        clean = _valid_values(values, xmin, xmax)
        if clean.size == 0:
            continue
        ax.hist(clean, bins=bins, histtype="step", density=True, linewidth=1.8, label=f"{label} ({clean.size})")
        any_drawn = True

    if not any_drawn:
        plt.close(fig)
        return

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Unit-area density")
    ax.set_xlim(xmin, xmax)
    ax.grid(True, alpha=0.25)
    ax.set_title(f"H #rightarrow Za angular distribution ({title_suffix})")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()

    stem = name_template.format(level=level)
    for ext in ("pdf", "png"):
        fig.savefig(output_dir / f"{stem}.{ext}")
    plt.close(fig)


def process_sample(sample: SampleArrays, args: argparse.Namespace) -> Dict[str, Dict[str, np.ndarray]]:
    out: Dict[str, Dict[str, np.ndarray]] = {}

    gen_vecs = vectors_from_nano_genpart(sample.arrays)
    gen_source = "GenPart"
    if gen_vecs is None:
        gen_vecs = vectors_from_flat_gen(sample.arrays)
        gen_source = "GenHza flat branches"
    if gen_vecs is not None:
        mask = base_valid_mask(gen_vecs)
        if np.any(mask):
            out["gen"] = compute_observables(apply_mask_to_vecs(gen_vecs, mask))
            print(f"[gen] {sample.label}: {int(np.sum(mask))} events from {gen_source}")

    reco_vecs = vectors_from_reco(sample.arrays)
    if reco_vecs is not None:
        mask = base_valid_mask(reco_vecs)
        if np.any(mask):
            out["reco"] = compute_observables(apply_mask_to_vecs(reco_vecs, mask))
            print(f"[reco] {sample.label}: {int(np.sum(mask))} baseline events")

            final_mask = final_selection_mask(sample.arrays, sample.mass, args)
            if final_mask is not None:
                final_mask = np.asarray(final_mask, dtype=bool)
                if final_mask.shape[0] == mask.shape[0]:
                    combined = mask & final_mask
                    if np.any(combined):
                        out["reco_final"] = compute_observables(apply_mask_to_vecs(reco_vecs, combined))
                        print(f"[final] {sample.label}: {int(np.sum(combined))} final/BDT-selected events")
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot normalized H -> Za angular observables from ROOT/parquet/nano-like signal files."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="inputs",
        action="append",
        default=[],
        help="Input ROOT/parquet file or directory. Can be repeated. Defaults to the Run 3 P2Root signal input base.",
    )
    parser.add_argument("--tree", default="inclusive", help="ROOT tree name, with automatic fallback to inclusive/Events.")
    parser.add_argument(
        "-o",
        "--output-dir",
        default="Plot/plots/hza_angular_distributions",
        help="Directory for PDF/PNG outputs.",
    )
    parser.add_argument(
        "--masses",
        default=",".join(str(m) for m in DEFAULT_MASSES),
        help="Comma-separated ALP masses to overlay. Use 'all' to keep all discovered masses.",
    )
    parser.add_argument(
        "--years",
        default="",
        help="Optional comma-separated year/era strings used to filter input paths, e.g. 2022preEE,2022postEE.",
    )
    parser.add_argument("--final-score-branch", default=None, help="Score branch for reco_final selection.")
    parser.add_argument("--final-score-cut", type=float, default=None, help="Cut value for --final-score-branch or auto score branch.")
    parser.add_argument(
        "--final-cut-expression",
        default=None,
        help="Numpy expression using branch names, e.g. '(BDT_Score_mA_M5 > 0.95) & (H_mass > 115) & (H_mass < 135)'.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    inputs = args.inputs or [DEFAULT_INPUT_BASE]
    if args.masses.strip().lower() == "all":
        masses: List[int] = []
    else:
        masses = [int(x) for x in re.split(r"[, ]+", args.masses.strip()) if x]
    years = [x for x in re.split(r"[, ]+", args.years.strip()) if x]

    files_by_label = discover_files(inputs, masses, years)
    if not files_by_label:
        raise SystemExit(f"No ROOT/parquet files found under: {', '.join(inputs)}")

    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    by_level: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {"gen": {}, "reco": {}, "reco_final": {}}
    for label, paths in files_by_label.items():
        mass_value = infer_mass_from_path(paths[0])
        print(f"[input] {label}: {len(paths)} file(s)")
        arrays = read_files(paths, args.tree)
        sample = SampleArrays(label=label, mass=mass_value, arrays=arrays, source_files=paths)
        sample_obs = process_sample(sample, args)
        for level, obs_map in sample_obs.items():
            for observable, values in obs_map.items():
                by_level[level].setdefault(observable, {})[label] = values

    level_titles = {
        "gen": "gen level, before reco selection when full GenPart is available",
        "reco": "reco level, baseline selected candidates",
        "reco_final": "reco level, final/BDT selection",
    }
    for level, observable_map in by_level.items():
        if not observable_map:
            continue
        for observable in PLOT_SPECS:
            values_by_label = observable_map.get(observable, {})
            if values_by_label:
                plot_observable(values_by_label, observable, level, output_dir, level_titles[level])

    print(f"[output] Wrote plots to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
