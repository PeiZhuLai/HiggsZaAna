"""Double-Sided Crystal Ball (DSCB) lepton pT response model.

The H->Z(ll)gamma kinematic refit (``tools/zmass_constraint_refit.py``) floats
the two lepton pT against the Z line shape.  Each floated pT is held near its
reconstructed value by a *constraint* term, which historically is a Gaussian

    log P(pt) = -0.5 * ( (pt - pt_reco) / sigma_pt )^2

with ``sigma_pt`` = the muon pT error (Run 3: ``Muon_bsConstrainedPtErr``, the
error of the beam-spot-constrained pT; Run 2: ``Muon_ptErr``) for muons and
``Electron_energyErr * sin(theta)`` for electrons.  That is a statement about
the per-lepton response

    r = (pt_reco - pt_true) / sigma_pt

being a unit Gaussian.  It is not: measured against the dressed generator
lepton the pull is asymmetric (a long low-side tail from bremsstrahlung /
un-recovered FSR) and roughly 40-60% too wide, so the Gaussian both
under-weights the tails and mis-states the core width.

This module provides the alternative constraint: a Double-Sided Crystal Ball in
the SAME pull variable ``r``, i.e. a Gaussian core with independent power-law
tails on each side,

    t = (r - mean) / sigma

    f(t) = A_L * (B_L - t)^-nL              t < -alphaL
           exp(-t^2 / 2)                    -alphaL <= t <= alphaR
           A_R * (B_R + t)^-nR              t > alphaR

    A_X = (nX / alphaX)^nX * exp(-alphaX^2 / 2)
    B_X = nX / alphaX - alphaX

with ``A_X``/``B_X`` fixed by continuity and continuity of the first derivative
at the two transition points.  ``mean = 0, sigma = 1, alphaL = alphaR = inf``
reproduces the Gaussian exactly, so the DSCB is a strict generalization and the
parameters are dimensionless (they multiply the existing ``sigma_pt``, they do
not replace it).

Two things live here:

* the PDF itself (:func:`dscb_pdf`, :func:`dscb_log_pdf`), analytically
  normalized and fully protected against the pathological corners
  (``sigma -> 0``, ``alpha -> 0``, ``n -> 1``, ``|t| -> inf``);
* the calibration lookup (:class:`LeptonPtDSCBCalibration`) that maps
  (era, flavor, reco pT, reco eta) onto a fitted parameter set.

The parameters themselves are measured by ``studies/lepton_pt_dscb`` from MC
truth; nothing in this module fits anything.
"""

from __future__ import annotations

import json
import os
from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np

# ── numerical guards ─────────────────────────────────────────────────────────
# Every one of these is a HARD clamp applied inside the PDF, not just a fit
# boundary: the refit evaluates the constraint on a (N, G, G) grid where a
# single inf/NaN would silently poison an event's argmin.
MIN_SIGMA = 1e-3      # pull units; sigma -> 0 would make the PDF a delta
MIN_ALPHA = 1e-2      # alpha -> 0 puts the transition at the peak and A_X -> inf
MIN_N = 1.0 + 1e-3    # n <= 1 makes the power-law tail non-integrable
MAX_N = 200.0         # n -> inf is a Gaussian tail; clamping keeps A_X finite
_TINY = 1e-300

PARAM_NAMES: Tuple[str, ...] = ("mean", "sigma", "alphaL", "nL", "alphaR", "nR")

# The parameter set that makes the DSCB numerically identical to the current
# Gaussian constraint.  Used by the unit tests and as the ultimate fallback for
# a bin with no usable calibration.
GAUSSIAN_LIMIT_PARAMS: Dict[str, float] = {
    "mean": 0.0, "sigma": 1.0,
    "alphaL": 50.0, "nL": 10.0,
    "alphaR": 50.0, "nR": 10.0,
}

FLAVORS: Tuple[str, ...] = ("electron", "muon")

# NanoAOD |pdgId| -> flavor label used throughout the calibration file.
PDGID_TO_FLAVOR: Dict[int, str] = {11: "electron", 13: "muon"}


# ── error function (scipy when available, else a pure-numpy fallback) ────────
try:  # pragma: no cover - exercised implicitly by whichever branch is taken
    from scipy.special import erf as _erf
except ImportError:  # pragma: no cover
    def _erf(x):
        """Abramowitz & Stegun 7.1.26; |error| < 1.5e-7, ample for a norm."""
        x = np.asarray(x, dtype=np.float64)
        sign = np.sign(x)
        z = np.abs(x)
        t = 1.0 / (1.0 + 0.3275911 * z)
        poly = t * (
            0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))
        )
        return sign * (1.0 - poly * np.exp(-z * z))


_SQRT_HALF_PI = np.sqrt(0.5 * np.pi)
_SQRT_HALF = np.sqrt(0.5)


def sanitize_params(mean, sigma, alphaL, nL, alphaR, nR):
    """Clamp a parameter set into the region where the DSCB is a proper PDF.

    Returns float64 arrays broadcastable together.  Non-finite inputs fall back
    to the Gaussian-limit value for that parameter rather than propagating NaN,
    because a NaN constraint silently removes an event from the refit.
    """
    def _clean(value, default, lo=None, hi=None):
        arr = np.asarray(value, dtype=np.float64)
        arr = np.where(np.isfinite(arr), arr, default)
        if lo is not None:
            arr = np.maximum(arr, lo)
        if hi is not None:
            arr = np.minimum(arr, hi)
        return arr

    mean = _clean(mean, GAUSSIAN_LIMIT_PARAMS["mean"])
    sigma = _clean(sigma, GAUSSIAN_LIMIT_PARAMS["sigma"], lo=MIN_SIGMA)
    alphaL = _clean(alphaL, GAUSSIAN_LIMIT_PARAMS["alphaL"], lo=MIN_ALPHA)
    alphaR = _clean(alphaR, GAUSSIAN_LIMIT_PARAMS["alphaR"], lo=MIN_ALPHA)
    nL = _clean(nL, GAUSSIAN_LIMIT_PARAMS["nL"], lo=MIN_N, hi=MAX_N)
    nR = _clean(nR, GAUSSIAN_LIMIT_PARAMS["nR"], lo=MIN_N, hi=MAX_N)
    return mean, sigma, alphaL, nL, alphaR, nR


def dscb_norm(sigma, alphaL, nL, alphaR, nR):
    """Integral of the unnormalized DSCB over the whole real line.

    Closed form -- core is an erf, each tail is a power-law integral:

        N / sigma = sqrt(pi/2) * [erf(alphaR/sqrt2) + erf(alphaL/sqrt2)]
                  + exp(-alphaL^2/2) * nL / (alphaL * (nL - 1))
                  + exp(-alphaR^2/2) * nR / (alphaR * (nR - 1))

    Inputs are assumed already sanitized (``nX > 1``, ``alphaX > 0``,
    ``sigma > 0``), which is what makes both tail terms finite.
    """
    core = _SQRT_HALF_PI * (_erf(alphaR * _SQRT_HALF) + _erf(alphaL * _SQRT_HALF))
    left = np.exp(-0.5 * alphaL * alphaL) * nL / (alphaL * (nL - 1.0))
    right = np.exp(-0.5 * alphaR * alphaR) * nR / (alphaR * (nR - 1.0))
    return sigma * (core + left + right)


def dscb_log_pdf(x, mean, sigma, alphaL, nL, alphaR, nR, normalized=True):
    """log of the (normalized) DSCB density evaluated at ``x``.

    All arguments broadcast against each other, so per-event parameter arrays
    can be evaluated against a per-event (or per-grid-point) ``x`` in one call.

    The tails are computed in log space, which is what keeps them finite for
    arbitrarily large ``|x|``: the power law becomes ``-n * log(...)``, a number
    of order a few hundred at worst, instead of an underflow to 0 followed by
    ``log(0) = -inf``.

    ``normalized=False`` drops the ``log(N)`` term.  The refit does not need it
    (a per-event constant cannot move the per-event minimum) but keeping the
    default normalized means an NLL from this function is directly comparable
    across parameter sets, which the calibration fit does need.
    """
    mean, sigma, alphaL, nL, alphaR, nR = sanitize_params(
        mean, sigma, alphaL, nL, alphaR, nR
    )
    x = np.asarray(x, dtype=np.float64)
    t = (x - mean) / sigma

    # Core.  Evaluated everywhere; np.where picks it only in the core region.
    log_core = -0.5 * t * t

    # Left tail.  B_L - t >= n/alpha > 0 wherever this branch is selected; the
    # maximum() only protects the (discarded) evaluation outside it.
    log_a_left = nL * np.log(nL / alphaL) - 0.5 * alphaL * alphaL
    u_left = np.maximum(nL / alphaL - alphaL - t, _TINY)
    log_left = log_a_left - nL * np.log(u_left)

    log_a_right = nR * np.log(nR / alphaR) - 0.5 * alphaR * alphaR
    u_right = np.maximum(nR / alphaR - alphaR + t, _TINY)
    log_right = log_a_right - nR * np.log(u_right)

    out = np.where(t < -alphaL, log_left, np.where(t > alphaR, log_right, log_core))
    if normalized:
        out = out - np.log(dscb_norm(sigma, alphaL, nL, alphaR, nR))
    # A non-finite result can only come from a non-finite x; map it to a
    # vanishing (not undefined) probability so the caller's argmin stays sane.
    return np.where(np.isfinite(out), out, -np.inf)


def dscb_pdf(x, mean, sigma, alphaL, nL, alphaR, nR, normalized=True):
    """Normalized DSCB density; ``exp`` of :func:`dscb_log_pdf`, never negative."""
    return np.exp(dscb_log_pdf(x, mean, sigma, alphaL, nL, alphaR, nR,
                               normalized=normalized))


# ── calibration file ─────────────────────────────────────────────────────────
# Schema (written by studies/lepton_pt_dscb/derive_parameters.py):
#
# {
#   "metadata": {"version": ..., "pull_definition": ..., ...},
#   "binning": {
#       "pt_edges": [...],                      # GeV, last edge may be "inf";
#                                               # or {"electron": [...], "muon": [...]}
#       "eta_edges": {"electron": [...], "muon": [...]},  # |eta|
#       "aux_var":   {"muon": "abs_dxy", ...},            # OPTIONAL 3rd axis
#       "aux_edges": {"muon": [...], ...},                # OPTIONAL, per flavor
#       "era_overrides": {                                # OPTIONAL, per era
#           "2024": {"muon": {"pt_edges": [...], "eta_edges": [...],
#                             "aux_var": "abs_dxy", "aux_edges": [...]}}}
#   },
#   "era_fallback": {"2026": "2025", ...},      # optional
#   "eras": {
#     "2018": {
#       "electron": {
#         "pt_20_30": {"eta_0p0_0p8": {"mean": ..., "sigma": ..., "alphaL": ...,
#                                      "nL": ..., "alphaR": ..., "nR": ...,
#                                      "errors": {...}, "fit_status": ...,
#                                      "n_events": ..., "pt_range": [...],
#                                      "eta_range": [...]}, ...}, ...}, ...}, ...}
# }
#
# The nested dict is the human-readable form.  On load it is flattened into
# (n_pt, n_eta, n_aux) numpy arrays per (era, flavor) so the per-event lookup is
# a few searchsorted calls and a fancy index, with no Python loop over events.
#
# The third axis is OPTIONAL and per flavor.  A flavor without "aux_edges" gets
# a single bin spanning (-inf, inf) and its entries stay two levels deep, so a
# calibration written before the axis existed loads unchanged.  Where it IS
# declared the entries gain a third level keyed like the others, e.g.
# "pt_20_30" -> "eta_0p0_0p9" -> "abs_dxy_0_0p0005".
#
# "era_overrides" lets one (era, flavor) use its own binning, e.g. Run 3 muons
# re-binned for Muon_bsConstrainedPtErr while Run 2 muons keep theirs. An
# override is COMPLETE for that flavor: pt_edges and eta_edges are required and
# the third axis exists only if it names aux_var and aux_edges itself; nothing
# is inherited from the top-level binning. Overrides are looked up by the era
# the table is actually read from (after era_fallback).

CALIBRATION_SCHEMA_VERSION = 1

# Default location of the calibration produced by the study.  Overridable via
# the ``HZG_LEPTON_PT_DSCB_JSON`` environment variable so a production run can
# point at a versioned copy on EOS without a code change.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_CALIBRATION_PATH = os.path.join(
    _REPO_ROOT, "studies", "lepton_pt_dscb", "results",
    "lepton_pt_dscb_parameters.json",
)
CALIBRATION_PATH_ENV = "HZG_LEPTON_PT_DSCB_JSON"


def _as_float_edges(edges: Sequence) -> np.ndarray:
    """Edge list -> float64 array; the string "inf" is allowed as the last edge."""
    return np.array(
        [np.inf if isinstance(e, str) and e.lower() in ("inf", "+inf") else float(e)
         for e in edges],
        dtype=np.float64,
    )


def format_bin_key(prefix: str, lo: float, hi: float) -> str:
    """Canonical bin key, e.g. ``pt_20_30`` / ``eta_0p0_0p8``.

    ``.`` is not legal in the flat key style used elsewhere in the repo's JSON,
    so it is spelled ``p``; ``inf`` stays literal.
    """
    def _fmt(value: float) -> str:
        if not np.isfinite(value):
            return "inf"
        text = ("%g" % value)
        return text.replace(".", "p").replace("-", "m")

    return f"{prefix}_{_fmt(lo)}_{_fmt(hi)}"


class LeptonPtDSCBCalibration:
    """Fitted DSCB parameters, looked up by (era, flavor, reco pT, reco eta).

    Out-of-range handling is explicit and never silent:

    * ``pt`` or ``|eta|`` below the first edge or above the last edge is
      CLAMPED into the edge bin, and the returned ``in_range`` mask is False
      there.  Clamping (rather than falling back to the Gaussian) is the right
      default because the edge bins are exactly the extrapolation the response
      is smoothest in, but the caller can always use ``in_range`` to do
      something else.
    * an unknown era raises ``KeyError`` unless the file declares a fallback
      for it in ``era_fallback`` -- production must never silently use another
      year's resolution.
    * a bin whose fit did not produce usable parameters is stored with
      ``fit_status`` recorded and the parameters it fell back to; see the study
      README.  Nothing is invented at lookup time.
    """

    def __init__(self, payload: Dict, source: str = "<dict>"):
        self.source = source
        self.metadata = dict(payload.get("metadata", {}))
        self.era_fallback = dict(payload.get("era_fallback", {}))

        binning = payload["binning"]
        self.eta_edges = {
            flavor: _as_float_edges(edges)
            for flavor, edges in binning["eta_edges"].items()
        }
        # pT edges are either shared by all flavors (a list) or per flavor (a
        # dict): electrons are split further in |eta| and sc_et_over_pt, so they
        # need coarser pT bins than muons to keep each cell populated.
        raw_pt = binning["pt_edges"]
        if isinstance(raw_pt, dict):
            missing = set(self.eta_edges) - set(raw_pt)
            if missing:
                raise KeyError(f"{source}: binning.pt_edges has no entry for {sorted(missing)}")
            self.pt_edges_by_flavor = {f: _as_float_edges(raw_pt[f]) for f in self.eta_edges}
            self.pt_edges = dict(self.pt_edges_by_flavor)
        else:
            shared = _as_float_edges(raw_pt)
            self.pt_edges_by_flavor = {f: shared for f in self.eta_edges}
            self.pt_edges = shared
        # Optional third axis. Absent -> one bin spanning everything, which
        # makes the 3-D lookup below reduce exactly to the old 2-D one.
        aux_edges = binning.get("aux_edges") or {}
        aux_var = binning.get("aux_var") or {}
        self.aux_var = {f: str(aux_var[f]) for f in aux_var}
        self.aux_edges = {}
        for flavor in self.eta_edges:
            if flavor in aux_edges:
                self.aux_edges[flavor] = _as_float_edges(aux_edges[flavor])
            else:
                self.aux_edges[flavor] = np.array([-np.inf, np.inf])

        # Per-(era, flavor) binning: the top-level one unless era_overrides
        # replaces it (see the schema note above).
        self.era_overrides = {}
        for era, per_flavor in (binning.get("era_overrides") or {}).items():
            for flavor, spec in per_flavor.items():
                if flavor not in self.eta_edges:
                    raise KeyError(f"{source}: era_overrides[{era!r}] names unknown "
                                   f"flavor {flavor!r}")
                missing = {"pt_edges", "eta_edges"} - set(spec)
                if missing:
                    raise KeyError(f"{source}: era_overrides[{era!r}][{flavor!r}] "
                                   f"lacks {sorted(missing)} (overrides are complete)")
                has_var, has_edges = spec.get("aux_var"), spec.get("aux_edges")
                if bool(has_var) != bool(has_edges):
                    raise KeyError(f"{source}: era_overrides[{era!r}][{flavor!r}] "
                                   "needs both aux_var and aux_edges, or neither")
                self.era_overrides[(era, flavor)] = (
                    _as_float_edges(spec["pt_edges"]),
                    _as_float_edges(spec["eta_edges"]),
                    str(has_var) if has_var else None,
                    _as_float_edges(has_edges) if has_edges else np.array([-np.inf, np.inf]),
                )

        self._tables: Dict[Tuple[str, str], Dict[str, np.ndarray]] = {}
        self._status: Dict[Tuple[str, str], np.ndarray] = {}
        self._n_events: Dict[Tuple[str, str], np.ndarray] = {}
        for era, per_flavor in payload["eras"].items():
            for flavor, per_pt in per_flavor.items():
                self._tables[(era, flavor)] = self._build_table(flavor, per_pt, era)
        unused = set(self.era_overrides) - set(self._tables)
        if unused:
            raise KeyError(f"{source}: era_overrides for {sorted(unused)} have no table")
        self.eras = sorted({era for era, _ in self._tables})

    def binning_for(self, era: Optional[str], flavor: str):
        """(pt_edges, eta_edges, aux_var or None, aux_edges) used for ``era``.

        ``era`` is the era the table is read from (already fallback-resolved);
        None gives the top-level binning.
        """
        if era is not None and (era, flavor) in self.era_overrides:
            return self.era_overrides[(era, flavor)]
        return (self.pt_edges_by_flavor[flavor], self.eta_edges[flavor],
                self.aux_var.get(flavor), self.aux_edges[flavor])

    # -- construction ---------------------------------------------------------
    def _build_table(self, flavor: str, per_pt: Dict,
                     era: Optional[str] = None) -> Dict[str, np.ndarray]:
        pt_edges, eta_edges, aux_name, aux_edges = self.binning_for(era, flavor)
        n_pt = len(pt_edges) - 1
        n_eta = len(eta_edges) - 1
        n_aux = len(aux_edges) - 1
        shape = (n_pt, n_eta, n_aux)
        table = {
            name: np.full(shape, GAUSSIAN_LIMIT_PARAMS[name], dtype=np.float64)
            for name in PARAM_NAMES
        }
        table["_filled"] = np.zeros(shape, dtype=bool)
        table["_n_events"] = np.zeros(shape, dtype=np.float64)
        table["_status"] = np.empty(shape, dtype=object)
        table["_status"][:] = "missing"

        def _store(entry, i, j, k):
            for name in PARAM_NAMES:
                table[name][i, j, k] = float(entry[name])
            table["_filled"][i, j, k] = True
            table["_n_events"][i, j, k] = float(entry.get("n_events", 0.0))
            table["_status"][i, j, k] = str(entry.get("fit_status", "unknown"))

        for i in range(n_pt):
            pt_key = format_bin_key("pt", pt_edges[i], pt_edges[i + 1])
            per_eta = per_pt.get(pt_key)
            if per_eta is None:
                continue
            for j in range(n_eta):
                eta_key = format_bin_key("eta", eta_edges[j], eta_edges[j + 1])
                entry = per_eta.get(eta_key)
                if entry is None:
                    continue
                if aux_name is None:
                    # No third axis: the entry IS the bin (old two-level form).
                    _store(entry, i, j, 0)
                    continue
                for k in range(n_aux):
                    aux_key = format_bin_key(aux_name, aux_edges[k], aux_edges[k + 1])
                    sub = entry.get(aux_key)
                    if sub is not None:
                        _store(sub, i, j, k)
        return table

    @classmethod
    def from_json(cls, path: Optional[str] = None) -> "LeptonPtDSCBCalibration":
        """Load from ``path``, ``$HZG_LEPTON_PT_DSCB_JSON`` or the study default."""
        if path is None:
            path = os.environ.get(CALIBRATION_PATH_ENV, DEFAULT_CALIBRATION_PATH)
        with open(path) as handle:
            payload = json.load(handle)
        version = payload.get("metadata", {}).get("schema_version")
        if version is not None and int(version) != CALIBRATION_SCHEMA_VERSION:
            raise ValueError(
                f"{path}: calibration schema version {version} != "
                f"{CALIBRATION_SCHEMA_VERSION} expected by this code"
            )
        return cls(payload, source=path)

    @property
    def electron_pterr_definition(self) -> str:
        """Electron sigma_pt the table was derived with: "raw" or "n2p".

        Tables written before this was recorded used the raw NanoAOD
        ``Electron_energyErr * sin(theta)``.
        """
        return str(self.metadata.get("electron_pterr_definition", "raw"))

    @property
    def muon_pterr_definition(self) -> str:
        """Muon sigma_pt the table was derived with: "raw" or "n2p".

        "raw" is ``Muon_ptErr`` in every era; "n2p" is ``Muon_bsConstrainedPtErr``
        in Run 3 and ``Muon_ptErr`` in Run 2 (nano2pico mu_producer.cpp). Tables
        written before this was recorded used the raw one.
        """
        return str(self.metadata.get("muon_pterr_definition", "raw"))

    # -- lookup ---------------------------------------------------------------
    def resolve_era(self, era: str) -> str:
        """The era whose table is actually used for ``era`` (may be a fallback)."""
        if any(key[0] == era for key in self._tables):
            return era
        target = self.era_fallback.get(era)
        if target is not None and any(key[0] == target for key in self._tables):
            return target
        raise KeyError(
            f"{self.source}: no DSCB calibration for era {era!r} "
            f"(available: {self.eras}; declared fallbacks: {sorted(self.era_fallback)})"
        )

    def bin_indices(self, flavor: str, pt, eta, aux=None, era: Optional[str] = None):
        """(i_pt, i_eta, i_aux, in_range) for reco ``pt`` [GeV] and ``eta``.

        ``eta`` is used as ``|eta|``; the binning is symmetric by construction.
        ``aux`` is the flavor's third-axis variable (see ``aux_var``); it is
        ignored for flavors without one, and a flavor that HAS one but is
        called without it falls into bin 0 and is reported out of range rather
        than silently taking the first bin's parameters.  ``era`` (the era the
        table is read from) selects its era_overrides binning, if any.
        """
        pt_edges, eta_edges, aux_name, aux_edges = self.binning_for(era, flavor)
        pt = np.asarray(pt, dtype=np.float64)
        abs_eta = np.abs(np.asarray(eta, dtype=np.float64))

        n_pt = len(pt_edges) - 1
        n_eta = len(eta_edges) - 1
        n_aux = len(aux_edges) - 1
        i_pt = np.clip(np.searchsorted(pt_edges, pt, side="right") - 1, 0, n_pt - 1)
        i_eta = np.clip(
            np.searchsorted(eta_edges, abs_eta, side="right") - 1, 0, n_eta - 1
        )
        in_range = (
            np.isfinite(pt) & np.isfinite(abs_eta)
            & (pt >= pt_edges[0]) & (pt < pt_edges[-1])
            & (abs_eta >= eta_edges[0]) & (abs_eta < eta_edges[-1])
        )

        if n_aux == 1:
            return i_pt, i_eta, np.zeros(pt.shape, dtype=np.int64), in_range
        if aux is None:
            raise ValueError(
                f"{self.source}: flavor {flavor!r} is binned in "
                f"{aux_name!r}; lookup needs that value"
            )
        aux_arr = np.asarray(aux, dtype=np.float64)
        i_aux = np.clip(
            np.searchsorted(aux_edges, aux_arr, side="right") - 1, 0, n_aux - 1
        )
        in_range &= (
            np.isfinite(aux_arr)
            & (aux_arr >= aux_edges[0]) & (aux_arr < aux_edges[-1])
        )
        return i_pt, i_eta, i_aux, in_range

    def lookup(self, era: str, flavor: str, pt, eta, aux=None) -> Dict[str, np.ndarray]:
        """Per-entry DSCB parameters for one era and one lepton flavor.

        Returns a dict with the six :data:`PARAM_NAMES` arrays plus
        ``in_range`` (bool) and ``n_events`` (the MC statistics behind the bin
        each entry landed in), all shaped like ``pt``.
        """
        if flavor not in self.eta_edges:
            raise KeyError(f"{self.source}: unknown flavor {flavor!r}")
        resolved = self.resolve_era(era)
        table = self._tables[(resolved, flavor)]
        i_pt, i_eta, i_aux, in_range = self.bin_indices(flavor, pt, eta, aux, era=resolved)
        out = {name: table[name][i_pt, i_eta, i_aux] for name in PARAM_NAMES}
        out["in_range"] = in_range
        out["n_events"] = table["_n_events"][i_pt, i_eta, i_aux]
        return out

    def lookup_by_pdgid(self, era: str, pdg_id, pt, eta, aux=None,
                        flavors=None) -> Dict[str, np.ndarray]:
        """Mixed-flavor lookup driven by ``|pdgId|`` (11 = e, 13 = mu).

        Rows with any other pdgId keep the Gaussian-limit parameters and are
        reported as out of range, so a mis-labelled row degrades to exactly the
        current behaviour instead of silently taking the wrong flavor.

        ``flavors`` restricts which flavors actually receive measured
        parameters; the rest get the Gaussian limit, which the DSCB reproduces
        exactly, so they behave as if the option had never been switched on.
        They are reported ``in_range=True``: running Gaussian there is a
        deliberate configuration, not a bin that fell off the table.

        This exists because the measured gain is entirely in the electron
        channel (+5.4% on m_ll out of sample) while muons come out at -0.1 to
        -0.2%, i.e. nothing -- their pull is already close to Gaussian and the
        Z mass constraint, not the per-lepton width, sets the refit resolution.
        ``flavors=("electron",)`` takes the gain without paying for the muon
        side.
        """
        if flavors is not None:
            unknown = set(flavors) - set(FLAVORS)
            if unknown:
                raise ValueError(
                    f"{self.source}: unknown flavor(s) {sorted(unknown)}; "
                    f"choose from {list(FLAVORS)}"
                )
            flavors = set(flavors)
        pdg_id = np.abs(np.asarray(pdg_id, dtype=np.int64))
        pt = np.asarray(pt, dtype=np.float64)
        eta = np.asarray(eta, dtype=np.float64)
        out = {name: np.full(pt.shape, GAUSSIAN_LIMIT_PARAMS[name], dtype=np.float64)
               for name in PARAM_NAMES}
        out["in_range"] = np.zeros(pt.shape, dtype=bool)
        out["n_events"] = np.zeros(pt.shape, dtype=np.float64)
        for code, flavor in PDGID_TO_FLAVOR.items():
            sel = pdg_id == code
            if not np.any(sel):
                continue
            if flavors is not None and flavor not in flavors:
                # Deliberately Gaussian: the arrays already hold the Gaussian
                # limit, so only the reporting needs setting.
                out["in_range"][sel] = True
                continue
            aux_sel = None
            if aux is not None:
                aux_sel = np.asarray(aux, dtype=np.float64)[sel]
            sub = self.lookup(era, flavor, pt[sel], eta[sel], aux_sel)
            for key, value in sub.items():
                out[key][sel] = value
        return out

    # -- introspection --------------------------------------------------------
    def status_table(self, era: str, flavor: str) -> np.ndarray:
        """(n_pt, n_eta, n_aux) array of per-bin ``fit_status`` strings."""
        return self._tables[(self.resolve_era(era), flavor)]["_status"]

    def coverage(self, era: str, flavor: str) -> float:
        """Fraction of (pt, eta) bins that carry a real fit for this era/flavor."""
        filled = self._tables[(self.resolve_era(era), flavor)]["_filled"]
        return float(filled.mean())

    def describe(self) -> str:
        lines = [f"DSCB calibration from {self.source}",
                 f"  version   : {self.metadata.get('version', '?')}",
                 f"  eras      : {self.eras}"]
        for flavor, edges in sorted(self.eta_edges.items()):
            lines.append(f"  pT ({flavor}): {self.pt_edges_by_flavor[flavor].tolist()}")
            lines.append(f"  |eta| ({flavor}): {edges.tolist()}")
        for flavor, name in sorted(self.aux_var.items()):
            lines.append(f"  {name} ({flavor}): {self.aux_edges[flavor].tolist()}")
        for (era, flavor), (pt_e, eta_e, aux_n, aux_e) in sorted(self.era_overrides.items()):
            aux_txt = f", {aux_n} {aux_e.tolist()}" if aux_n else ""
            lines.append(f"  override {era} {flavor}: pT {pt_e.tolist()}, "
                         f"|eta| {eta_e.tolist()}{aux_txt}")
        for era in self.eras:
            for flavor in sorted(self.eta_edges):
                if (era, flavor) in self._tables:
                    lines.append("  coverage %-12s %-8s %.0f%%"
                                 % (era, flavor, 100.0 * self.coverage(era, flavor)))
        return "\n".join(lines)


_CALIBRATION_CACHE: Dict[str, LeptonPtDSCBCalibration] = {}


def load_calibration(path: Optional[str] = None) -> LeptonPtDSCBCalibration:
    """Cached :meth:`LeptonPtDSCBCalibration.from_json`.

    The processor calls this once per variation per chunk; parsing the JSON
    every time would be pure waste, and the object is read-only.
    """
    key = path or os.environ.get(CALIBRATION_PATH_ENV, DEFAULT_CALIBRATION_PATH)
    if key not in _CALIBRATION_CACHE:
        _CALIBRATION_CACHE[key] = LeptonPtDSCBCalibration.from_json(key)
    return _CALIBRATION_CACHE[key]


def clear_calibration_cache() -> None:
    """Drop the cached calibrations (tests that write temporary JSON need this)."""
    _CALIBRATION_CACHE.clear()


def stack_params(params: Dict[str, np.ndarray]) -> Tuple[np.ndarray, ...]:
    """(mean, sigma, alphaL, nL, alphaR, nR) tuple from a lookup dict."""
    return tuple(np.asarray(params[name], dtype=np.float64) for name in PARAM_NAMES)


def params_from_iterable(values: Iterable) -> Dict[str, float]:
    """Dict form of a 6-tuple in :data:`PARAM_NAMES` order."""
    values = list(values)
    if len(values) != len(PARAM_NAMES):
        raise ValueError(f"expected {len(PARAM_NAMES)} parameters, got {len(values)}")
    return {name: float(value) for name, value in zip(PARAM_NAMES, values)}
