# Copied from the HZgamma HiggsDNA (higgs_dna/tools/zmass_constraint_refit.py, repo HEAD 54303273) as the
# opt-in "hzgamma" Z refit of the HZa tagger (option z_refit_method / env HZA_Z_REFIT). Not used by default.

"""Z mass-constraint kinematic refit for H->Z(ll)gamma.

Vectorized NumPy port of the per-lepton pT likelihood fit used in nano2pico
(``src/KinZfitter.cpp``) and HiggsZGammaAna/hzgml (``scripts/z_refit.py``).

For each event the two lepton pT are floated (eta/phi/mass held fixed) to
maximize the likelihood

    L = Constraint(pt1; pt1_reco, sigma1)
      * Constraint(pt2; pt2_reco, sigma2)
      * Shape(m_ll)

The per-lepton ``Constraint`` is selected by the ``lepton_constraint=``
argument of :func:`refit_lepton_pts`:

  ``gaussian``  the NOMINAL model, ``Gauss(pt; pt_reco, sigma)`` -- i.e. the
                per-lepton pull ``r = (pt_reco - pt)/sigma`` is a unit Gaussian
  ``dscb``      a Double-Sided Crystal Ball in the same pull ``r``, with
                (mean, sigma, alphaL, nL, alphaR, nR) measured per era /
                flavor / pT / eta by ``studies/lepton_pt_dscb``. See
                :mod:`higgs_dna.tools.lepton_pt_dscb`.

The DSCB is a strict generalization: mean = 0, sigma = 1 and both alpha large
reproduce the Gaussian, so switching the option with those parameters is a
no-op. ``gaussian`` remains the default so existing production is unchanged.

``Shape(m_ll)`` is the Z lineshape.  Two models are available via the
``lineshape=`` argument of :func:`refit_lepton_pts`:

  ``legacy_cb3gauss``  the NOMINAL empirical shape (Crystal Ball + 3 Gaussians,
                       identical electron/muon parameters in z_refit.py)
  ``pythia_q``         the PYTHIA 8 relativistic running-width Breit-Wigner
                       q(m) used to generate the off-shell Z in H -> Z gamma
  ``pythia_q_smeared`` q(m) (x) Gauss(sigma_res): the truth shape turned into an
                       effective reco-level template. ``sigma_res`` is the only
                       parameter that may differ between Z->ee and Z->mumu; see
                       LINESHAPE_PARAMS_BY_CHANNEL and the ``channel`` argument.

The negative log-likelihood

    NLL = -[ logConstraint1 + logConstraint2 + log(Shape(m_ll)) ]

is minimized on a per-event 2D grid spanning +-``n_sigma`` * sigma around the
reconstructed pT (matching z_refit's +-2 sigma bounds and >=5 GeV floor).

FSR-photon recovery (muon channel) IS supported: a muon's associated FSR photon
is folded into m(ll) for the lineshape constraint and added back onto the parent
muon to build the dressed refit Z. Unlike KinZfitter's NLL_1 / NLL_2 the FSR
photon pT is held fixed at its reconstructed value (not independently floated);
this captures the FSR energy recovery while keeping the fit a fast 2-D grid.
Electrons get no FSR (their FSR is largely already in the supercluster energy),
matching nano2pico. The signal photon is added by the caller afterwards to build
the refit H candidate.
"""

import logging

import numpy as np
from math import erf as _erf_scalar


def erf(x):
    """Vectorized erf; scipy is not a hard dependency of this module."""
    return np.vectorize(_erf_scalar)(x)

from higgs_dna.tools.lepton_pt_dscb import (
    PARAM_NAMES as DSCB_PARAM_NAMES,
    dscb_log_pdf,
)

# ── per-lepton pT constraint models ──────────────────────────────────────────
# The nominal (and default) model is the Gaussian the analysis has always used;
# ``dscb`` swaps in the measured Double-Sided Crystal Ball response. Changing
# the default here changes production, so it is deliberately a named constant.
LEPTON_CONSTRAINT_GAUSSIAN = "gaussian"
LEPTON_CONSTRAINT_DSCB = "dscb"
LEPTON_CONSTRAINTS = (LEPTON_CONSTRAINT_GAUSSIAN, LEPTON_CONSTRAINT_DSCB)
DEFAULT_LEPTON_CONSTRAINT = LEPTON_CONSTRAINT_GAUSSIAN

logger = logging.getLogger(__name__)

# Z lineshape parameters.  In hzgml/z_refit.py the electron and muon parameter
# sets are identical, so a single set is used here.  Source:
# ELECTRON_LINESHAPE_PARAMS / MUON_LINESHAPE_PARAMS in scripts/z_refit.py.
LINESHAPE_PARAMS = {
    "meanCB": 90.8919,
    "sigmaCB": 4.00007,
    "alphaCB": 1.1981,
    "nCB": 3.25604,
    "meanGauss1": 96.4278,
    "sigmaGauss1": 6.17509,
    "meanGauss2": 91.1649,
    "sigmaGauss2": 0.856305,
    "meanGauss3": 91.1513,
    "sigmaGauss3": 1.77148,
    "f1": 0.86449,
    "f2": 0.514371,
    "f3": 0.648225,
}

# PYTHIA-style relativistic Breit-Wigner with a running width -- the Z lineshape
# PYTHIA 8 uses when generating the off-shell Z in H -> Z gamma. Unlike
# LINESHAPE_PARAMS (an empirical, already-smeared reco-level template) this is a
# TRUTH-level shape: it carries no detector resolution at all.
#
#   q(m) = [ 2 m^2 Gamma(m) / ( (m^2 - mZ^2)^2 + (m Gamma(m))^2 ) ]
#          * sqrt(1 - m^2 / mH^2)
#   Gamma(m) = GammaZ * sqrt( (m^2 - mthr^2) / (mZ^2 - mthr^2) )
#
# The sqrt(1 - m^2/mH^2) factor is the two-body H -> Z gamma phase space, so q(m)
# vanishes at the kinematic endpoint m = mH.
PYTHIA_LINESHAPE_PARAMS = {
    "mZ": 91.1876,      # GeV, Z pole mass
    "gammaZ": 2.4978,   # GeV, Z total width at the pole
    "mH": 125.0,        # GeV, Higgs mass (phase-space endpoint)
    "mthr": 1.2536,     # GeV, running-width threshold
}

# Grid-search configuration (faithful to z_refit's bounds).
DEFAULT_N_GRID = 11      # points scanned per lepton in the coarse pass (z_refit uses 11)
DEFAULT_N_SIGMA = 2.0    # half-window in units of sigma (z_refit bound: reco +- 2 sigma)
DEFAULT_PT_FLOOR = 5.0   # GeV, minimum allowed pT (z_refit bound)

# Second, refined pass around the coarse winner. The single 11-point pass alone
# quantizes m(ll)^refit badly enough to make it BIMODAL: with a 0.4*sigma step
# the grid cannot land on the Z pole, so events pile up on the two points
# straddling it and the pole bin sits ~11% below its neighbours. Measured on
# ggH 2018 (deficit of the 91.12 GeV bin relative to its neighbours):
#
#     n_grid  step     deficit
#     11      0.40 s   -11.5 %     <- production
#     21      0.20 s   -25.9 %     <- simply doubling makes it WORSE
#     31      0.13 s   -18.9 %
#     41      0.10 s    -7.5 %
#     51      0.08 s    -0.5 %     <- converged
#
# The artifact is not monotonic in the point count, so it has to be driven to a
# step well below the resolution scale. Brute force at 51 points costs 51^2/11^2
# = 21x; refining once around the coarse winner reaches the same 0.08*sigma
# effective step for 2*11^2 = 242 evaluations, i.e. 10.7x cheaper than that.
DEFAULT_N_REFINE = 11    # points in the refinement pass; <= 1 disables it

# Extra zooms after the refinement and the edge rescue: each re-scans n_refine
# points per lepton over +- the previous step around the current minimum, so the
# step shrinks by 2/(n_refine-1) = 5x per round (0.08 -> 0.016 -> 0.0032 ->
# 0.00064 sigma). This takes the result to the continuous minimum of the NLL,
# i.e. what a Minuit fit started in the same valley finds. 0 gives the
# two-pass (coarse + refinement) result.
DEFAULT_N_EXTRA_REFINE = 3

# Edge rescue. The coarse scan spans reco +- n_sigma; an event whose optimum
# lies outside simply gets the boundary point, which is not a fit result but the
# statement "the optimum is somewhere I was not allowed to look". Measured on
# ggF 2018, the fraction of events pinned to the edge:
#
#                       gaussian   dscb
#     legacy shape        0.00 %   0.09 %
#     pythia_q            0.10 %   2.57 %
#
# Both new options push it up and together they multiply: pythia_q is much
# narrower than the empirical shape so it pulls m(ll) harder, and the DSCB tails
# make the per-lepton term cheaper to violate. Widening the window for everyone
# would change the nominal scan definition (and cost time everywhere), so
# instead the events that actually hit the edge are re-scanned with the window
# doubled, repeatedly, until they come off it or the cap is reached.
EDGE_RESCUE_MAX_ROUNDS = 3      # 2 -> 4 -> 8 -> 16 sigma; 0 disables
EDGE_RESCUE_GROWTH = 2.0
_EDGE_TOL = 0.999               # matches the study's definition of "on the edge"

# Grid of the early pythia-lineshape production, still used by the
# zlineshape_pythia closure study and the determinism tests. It is NOT what the
# processor runs: that is DEFAULT_N_GRID / DEFAULT_N_REFINE plus
# DEFAULT_N_EXTRA_REFINE zooms (HZG_ZREFIT_EXTRA_REFINE), so reproducing stored
# lepton_*_pt_refit values needs those, and n_extra_refine=0 for productions
# made before the extra zooms.
PRODUCTION_N_GRID = 11
PRODUCTION_N_REFINE = 1

_EVENT_BATCH = 20000     # events per slice at the reference (11, 11) grid; the
                         # actual slice is scaled below so peak memory is flat

# Numeric-normalization grid for pdf_over_range(); 20001 points over a 20 GeV
# window is ~1 MeV spacing, far finer than any structure in either lineshape.
_NORM_N_POINTS = 20001


def _theta(eta):
    return 2.0 * np.arctan(np.exp(-eta))


def _crystal_ball(x, mean, sigma, alpha, n):
    """NORMALIZED Crystal Ball, matching nano2pico KinZfitter::evaluateShape.

    The normalization is not cosmetic. This shape is mixed with three
    NORMALIZED Gaussians by the fractions f1/f2/f3, so leaving the CB
    unnormalized does not rescale the total -- it rescales the CB RELATIVE to
    the Gaussians, which changes the shape and therefore moves the refit
    minimum. As coded before, the CB integrated to sigma*(C+D) = 11.2 instead
    of 1, giving it an effective weight of 0.82 where f1/f2/f3 intend 0.29
    (2.8x over-weighted). Since the CB is the broad component (sigma 4.0 GeV,
    against 0.86 and 1.77 for Gauss2/3), that made the whole prior far too wide
    and the Z constraint correspondingly too loose.

    n2p (src/KinZfitter.cpp):
        CBC = (n/alpha)*(1/(n-1))*exp(-alpha^2/2)      tail integral
        CBD = sqrt(pi/2)*(1+erf(alpha/sqrt(2)))        core integral
        CBt = CBt/(sigma*(CBC+CBD))

    Assumes ``alpha > 0`` and ``n > 1`` (both true for LINESHAPE_PARAMS).
    """
    t = (x - mean) / sigma
    a = (n / alpha) ** n * np.exp(-0.5 * alpha * alpha)
    b = n / alpha - alpha
    core = np.exp(-0.5 * t * t)
    tail = a * np.power(np.maximum(b - t, 1e-9), -n)
    shape = np.where(t >= -alpha, core, tail)
    c = (n / alpha) * (1.0 / (n - 1.0)) * np.exp(-0.5 * alpha * alpha)
    d = np.sqrt(np.pi / 2.0) * (1.0 + erf(alpha / np.sqrt(2.0)))
    return shape / (sigma * (c + d))


def _gauss(x, mean, sigma):
    return np.exp(-0.5 * ((x - mean) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))


def _z_lineshape(mz, p=LINESHAPE_PARAMS):
    cb = _crystal_ball(mz, p["meanCB"], p["sigmaCB"], p["alphaCB"], p["nCB"])
    g1 = _gauss(mz, p["meanGauss1"], p["sigmaGauss1"])
    g2 = _gauss(mz, p["meanGauss2"], p["sigmaGauss2"])
    g3 = _gauss(mz, p["meanGauss3"], p["sigmaGauss3"])
    a = p["f1"] * cb + (1.0 - p["f1"]) * g1
    b = p["f2"] * a + (1.0 - p["f2"]) * g2
    return p["f3"] * b + (1.0 - p["f3"]) * g3


def pythia_running_width(m, mZ=91.1876, gammaZ=2.4978, mthr=1.2536):
    """Relativistic running Z width Gamma(m) [GeV]; 0 at or below ``mthr``.

    Gamma(m) = gammaZ * sqrt( (m^2 - mthr^2) / (mZ^2 - mthr^2) )
    """
    m = np.asarray(m, dtype=np.float64)
    num = m * m - mthr * mthr
    den = mZ * mZ - mthr * mthr
    if den <= 0:
        raise ValueError(f"mZ ({mZ}) must exceed mthr ({mthr})")
    return gammaZ * np.sqrt(np.maximum(num, 0.0) / den)


def pythia_z_lineshape(
    m,
    mZ=91.1876,
    gammaZ=2.4978,
    mH=125.0,
    mthr=1.2536,
    norm_range=None,
):
    """PYTHIA-style relativistic running-width Z lineshape q(m).

        q(m) = [ 2 m^2 Gamma(m) / ( (m^2 - mZ^2)^2 + (m Gamma(m))^2 ) ]
               * sqrt(1 - m^2 / mH^2)

    with ``Gamma(m)`` from :func:`pythia_running_width`. This is the truth-level
    shape PYTHIA 8 uses for the off-shell Z in H -> Z gamma; it contains NO
    detector resolution (contrast :func:`_z_lineshape`, whose sigmaCB = 4 GeV
    makes it an already-smeared reco-level template).

    Fully protected: q(m) = 0 for m <= mthr, q(m) = 0 for m >= mH, no negative
    sqrt arguments, and any residual NaN/inf is mapped to 0. The result is
    guaranteed finite and non-negative.

    The raw expression is NOT normalized. Pass ``norm_range=(lo, hi)`` to divide
    by the numeric integral over that window, giving a proper PDF on [lo, hi].
    """
    m = np.asarray(m, dtype=np.float64)
    # Support region: strictly between the running-width threshold and the
    # H -> Z gamma phase-space endpoint.
    good = np.isfinite(m) & (m > mthr) & (m < mH)
    # Evaluate on a sanitized copy so out-of-support entries can never produce
    # a warning or a NaN that would then have to be cleaned up downstream.
    m_safe = np.where(good, m, mZ)

    m2 = m_safe * m_safe
    width = pythia_running_width(m_safe, mZ=mZ, gammaZ=gammaZ, mthr=mthr)
    m_width = m_safe * width

    denom = (m2 - mZ * mZ) ** 2 + m_width * m_width
    phase_space = np.sqrt(np.maximum(1.0 - m2 / (mH * mH), 0.0))

    value = np.where(
        denom > 0.0,
        2.0 * m2 * width / np.where(denom > 0.0, denom, 1.0) * phase_space,
        0.0,
    )
    value = np.where(good & np.isfinite(value), value, 0.0)
    value = np.maximum(value, 0.0)

    if norm_range is not None:
        value = value / _lineshape_integral(
            "pythia_q", tuple(float(x) for x in norm_range),
            _params_key({"mZ": mZ, "gammaZ": gammaZ, "mH": mH, "mthr": mthr}),
        )
    if np.ndim(m) == 0:
        return float(value)
    return value


def pythia_smeared_z_lineshape(
    m,
    mZ=91.1876,
    gammaZ=2.4978,
    mH=125.0,
    mthr=1.2536,
    sigma_res=1.7307,
    norm_range=None,
):
    """PYTHIA q(m) convolved with a Gaussian detector resolution ``sigma_res``.

    q(m) itself is truth level and carries no resolution, whereas the nominal
    CB+3Gauss template is measurably pre-smeared (its fitted sigma_res runs to
    the lower bound). This model is the middle ground: the PYTHIA truth line
    shape turned into an EFFECTIVE reco-level template by one explicit Gaussian
    smearing, whose width is the only channel-dependent parameter.

    The convolution is done once on a fine uniform grid and cached; evaluation
    is a linear interpolation, so this stays usable inside the per-event refit's
    (N, G, G) grid search. ``sigma_res <= 0`` falls back to the bare q(m).

    Unlike q(m) the result is NOT forced to zero above ``mH``: smearing a
    distribution that ends at the phase-space endpoint genuinely leaks past it,
    and truncating there would put a hard edge into the constraint.
    """
    if sigma_res is None or sigma_res <= 0:
        return pythia_z_lineshape(m, mZ=mZ, gammaZ=gammaZ, mH=mH, mthr=mthr,
                                  norm_range=norm_range)

    grid, values = _smeared_grid(
        _params_key({"mZ": mZ, "gammaZ": gammaZ, "mH": mH, "mthr": mthr,
                     "sigma_res": sigma_res})
    )
    m_arr = np.asarray(m, dtype=np.float64)
    finite = np.isfinite(m_arr)
    out = np.interp(np.where(finite, m_arr, grid[0]), grid, values, left=0.0, right=0.0)
    out = np.where(finite & np.isfinite(out), out, 0.0)
    out = np.maximum(out, 0.0)

    if norm_range is not None:
        out = out / _lineshape_integral(
            PYTHIA_SMEARED_LINESHAPE, tuple(float(x) for x in norm_range),
            _params_key({"mZ": mZ, "gammaZ": gammaZ, "mH": mH, "mthr": mthr,
                         "sigma_res": sigma_res}),
        )
    if np.ndim(m) == 0:
        return float(out)
    return out


# Grid used for the cached q(m) (x) Gauss convolution: 2 MeV spacing, extended
# 6 sigma past the phase-space endpoint so the smeared tail is not clipped.
_SMEAR_GRID_STEP = 0.01    # GeV; 50x finer than the 0.5 GeV analysis bins
_SMEAR_PAD_NSIGMA = 6.0
_SMEAR_KERNEL_NSIGMA = 6.0
_SMEARED_GRID_CACHE = {}


def _smeared_grid(params_key):
    """(grid, convolved values) for one parameter set; built once, then cached."""
    if params_key in _SMEARED_GRID_CACHE:
        return _SMEARED_GRID_CACHE[params_key]

    params = dict(params_key)
    sigma = params.pop("sigma_res")
    pad = _SMEAR_PAD_NSIGMA * sigma
    lo = max(0.0, params["mthr"] - pad)
    hi = params["mH"] + pad
    n_points = int(np.ceil((hi - lo) / _SMEAR_GRID_STEP)) + 1
    grid = np.linspace(lo, hi, n_points)
    step = grid[1] - grid[0]

    truth = pythia_z_lineshape(grid, **params)
    half = int(np.ceil(_SMEAR_KERNEL_NSIGMA * sigma / step))
    offsets = np.arange(-half, half + 1) * step
    kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
    kernel /= kernel.sum()
    values = np.convolve(truth, kernel, mode="same")
    values = np.maximum(np.where(np.isfinite(values), values, 0.0), 0.0)

    _SMEARED_GRID_CACHE[params_key] = (grid, values)
    return grid, values


# ── Lineshape model registry ─────────────────────────────────────────────────
# Each entry maps a model name to (default parameter dict, evaluator). The
# evaluator takes (m, params) and returns the UNNORMALIZED shape value.
LEGACY_LINESHAPE = "legacy_cb3gauss"   # the nominal model -- do not change
PYTHIA_LINESHAPE = "pythia_q"
PYTHIA_SMEARED_LINESHAPE = "pythia_q_smeared"
DEFAULT_LINESHAPE = LEGACY_LINESHAPE

LINESHAPE_MODELS = {
    LEGACY_LINESHAPE: (LINESHAPE_PARAMS, lambda m, p: _z_lineshape(m, p)),
    PYTHIA_LINESHAPE: (PYTHIA_LINESHAPE_PARAMS, lambda m, p: pythia_z_lineshape(m, **p)),
    PYTHIA_SMEARED_LINESHAPE: (
        dict(PYTHIA_LINESHAPE_PARAMS, sigma_res=1.7307),
        lambda m, p: pythia_smeared_z_lineshape(m, **p),
    ),
}

# ── Per-channel parameter sets ───────────────────────────────────────────────
# The TRUTH-level Z line shape cannot depend on the decay channel: mZ, GammaZ,
# mH and mthr are physics constants. What differs between Z->ee and Z->mumu is
# the detector resolution, so the only channel-dependent parameter here is
# sigma_res of the smeared model.
#
# The two values below are the sigma_res fitted per channel to the signal m_ll
# spectrum (plot/scripts/zlineshape_pythia_fit.py --resolution gauss, full
# Run2+Run3 signal): the electron channel needs nearly twice the muon width,
# which is the bremsstrahlung low-side tail showing up as extra smearing.
#
# For the legacy model both channels share one set, mirroring upstream
# hzgml/z_refit.py where ELECTRON_LINESHAPE_PARAMS and MUON_LINESHAPE_PARAMS
# are numerically identical.
CHANNELS = ("ee", "mumu")

PYTHIA_SIGMA_RES_BY_CHANNEL = {"ee": 2.3587, "mumu": 1.2955}

LINESHAPE_PARAMS_BY_CHANNEL = {
    LEGACY_LINESHAPE: {channel: dict(LINESHAPE_PARAMS) for channel in CHANNELS},
    PYTHIA_LINESHAPE: {channel: dict(PYTHIA_LINESHAPE_PARAMS) for channel in CHANNELS},
    PYTHIA_SMEARED_LINESHAPE: {
        channel: dict(PYTHIA_LINESHAPE_PARAMS,
                      sigma_res=PYTHIA_SIGMA_RES_BY_CHANNEL[channel])
        for channel in CHANNELS
    },
}


def resolve_lineshape_params(name=DEFAULT_LINESHAPE, params=None, channel=None):
    """Full parameter dict for ``name``, with ``params`` overriding the defaults.

    ``channel`` ("ee" or "mumu") selects the per-channel defaults from
    :data:`LINESHAPE_PARAMS_BY_CHANNEL` before ``params`` is applied on top.
    """
    if name not in LINESHAPE_MODELS:
        raise ValueError(
            f"unknown lineshape model {name!r}; choose from {sorted(LINESHAPE_MODELS)}"
        )
    if channel is None:
        merged = dict(LINESHAPE_MODELS[name][0])
    else:
        if channel not in CHANNELS:
            raise ValueError(f"unknown channel {channel!r}; choose from {CHANNELS}")
        merged = dict(LINESHAPE_PARAMS_BY_CHANNEL[name][channel])
    if params:
        merged.update(params)
    return merged


def build_lineshape(name=DEFAULT_LINESHAPE, params=None, norm_range=None,
                    channel=None):
    """Return a vectorized ``f(m_ll) -> shape value`` for the named model.

    ``norm_range=(lo, hi)`` returns a properly normalized PDF on [lo, hi]
    (numeric integral, so it works for every model). Left as None the shape is
    unnormalized -- which is all the per-event refit needs, since a constant
    factor shifts the NLL by a constant and cannot move the minimum.

    ``channel`` picks the per-channel parameter set; see
    :func:`resolve_lineshape_params`.
    """
    merged = resolve_lineshape_params(name, params, channel)
    _, evaluate = LINESHAPE_MODELS[name]
    if norm_range is None:
        return lambda m: evaluate(m, merged)
    integral = _lineshape_integral(
        name, tuple(float(x) for x in norm_range), _params_key(merged)
    )
    return lambda m: evaluate(m, merged) / integral


def _params_key(params):
    """Hashable, order-independent key for a parameter dict."""
    return tuple(sorted((k, float(v)) for k, v in params.items()))


_INTEGRAL_CACHE = {}


def _lineshape_integral(name, norm_range, params_key):
    """Numeric integral of an unnormalized lineshape over ``norm_range``.

    Composite Simpson on a fixed fine grid, cached per (model, range, params).
    Raises if the integral is not strictly positive -- a silent divide by ~0
    would poison every downstream number.
    """
    cache_key = (name, norm_range, params_key)
    if cache_key in _INTEGRAL_CACHE:
        return _INTEGRAL_CACHE[cache_key]

    lo, hi = norm_range
    if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
        raise ValueError(f"invalid normalization range {norm_range!r}")

    params = dict(params_key)
    _, evaluate = LINESHAPE_MODELS[name]
    # Odd point count so Simpson's rule applies over the whole grid.
    n_points = _NORM_N_POINTS if _NORM_N_POINTS % 2 else _NORM_N_POINTS + 1
    grid = np.linspace(lo, hi, n_points)
    values = np.asarray(evaluate(grid, params), dtype=np.float64)
    values = np.where(np.isfinite(values), values, 0.0)

    step = (hi - lo) / (n_points - 1)
    weights = np.ones(n_points)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    integral = float(step / 3.0 * np.dot(weights, values))

    if not np.isfinite(integral) or integral <= 0.0:
        raise ValueError(
            f"lineshape {name!r} integrates to {integral!r} over {norm_range!r}; "
            "cannot normalize"
        )
    _INTEGRAL_CACHE[cache_key] = integral
    return integral


def _log_lepton_constraint(pt_grid, pt_reco, err, dscb):
    """log P(pt_reco | pt) for one lepton, on the (N, G, G) candidate grid.

    ``pt_reco``/``err`` are (N,) and are broadcast onto the grid here. The pull
    is defined EXACTLY as the calibration measures it,

        r = (pt_reco - pt) / err

    i.e. reco minus truth, with the floated ``pt`` playing the role of truth.
    The sign matters for the DSCB and only for the DSCB: the Gaussian is even,
    but the measured response has a longer LOW-side tail (a reconstructed
    lepton loses energy, it does not gain it), so flipping the sign would put
    the tail on the wrong side.

    ``dscb`` is None for the Gaussian model, else the 6-tuple of per-event
    parameter arrays in :data:`DSCB_PARAM_NAMES` order.

    The Gaussian branch is kept as the bare ``-0.5 r^2`` the analysis has always
    used (its normalization is a per-event constant and cannot move the
    minimum); the DSCB branch is normalized, which likewise cannot move the
    per-event minimum but keeps NLL values comparable between parameter sets.
    """
    residual = (pt_reco[:, None, None] - pt_grid) / err[:, None, None]
    if dscb is None:
        return -0.5 * residual * residual
    return dscb_log_pdf(residual, *(p[:, None, None] for p in dscb))


def _refit_batch(pt1, eta1, phi1, m1, err1, pt2, eta2, phi2, m2, err2,
                 fsr_e, fsr_px, fsr_py, fsr_pz,
                 n_grid, n_sigma, pt_floor, shape_fn,
                 centre1=None, centre2=None, half1=None, half2=None,
                 dscb1=None, dscb2=None):
    """Core grid search for one batch of events. Returns (pt1_refit, pt2_refit).

    ``fsr_{e,px,py,pz}`` are the per-event SUM of the (muon-channel) FSR photon
    4-momenta to fold into the m(ll) used by the Z-lineshape constraint. They are
    held fixed at their reconstructed values (only the two lepton pT are floated),
    so FSR enters the mass constraint but is not itself re-fitted. Pass zeros for
    no FSR.

    ``dscb1``/``dscb2`` select the per-lepton constraint model: None gives the
    nominal Gaussian, a 6-tuple of (N,) parameter arrays gives the DSCB.
    """
    n = pt1.shape[0]
    # Default: the coarse pass, centred on the reconstructed pT with a
    # +-n_sigma*err half-window. The refinement pass passes its own centres and
    # (much smaller) half-widths instead.
    if centre1 is None:
        centre1, centre2 = pt1, pt2
        half1, half2 = n_sigma * err1, n_sigma * err2
    frac = np.linspace(-1.0, 1.0, n_grid)  # (G,), in units of the half-window

    cand1 = np.maximum(centre1[:, None] + frac[None, :] * half1[:, None], pt_floor)
    cand2 = np.maximum(centre2[:, None] + frac[None, :] * half2[:, None], pt_floor)

    # (N, G, G): axis 1 -> lepton-1 grid, axis 2 -> lepton-2 grid.
    p1 = cand1[:, :, None]
    p2 = cand2[:, None, :]

    th1 = _theta(eta1)[:, None, None]
    th2 = _theta(eta2)[:, None, None]
    sin1 = np.sin(th1)
    sin2 = np.sin(th2)
    cot1 = np.cos(th1) / sin1
    cot2 = np.cos(th2) / sin2

    e1 = np.sqrt(p1 ** 2 / sin1 ** 2 + m1[:, None, None] ** 2)
    e2 = np.sqrt(p2 ** 2 / sin2 ** 2 + m2[:, None, None] ** 2)
    # Add the fixed FSR-photon 4-momentum (per event) to the dilepton system.
    tot_e = e1 + e2 + fsr_e[:, None, None]
    px = p1 * np.cos(phi1)[:, None, None] + p2 * np.cos(phi2)[:, None, None] + fsr_px[:, None, None]
    py = p1 * np.sin(phi1)[:, None, None] + p2 * np.sin(phi2)[:, None, None] + fsr_py[:, None, None]
    pz = p1 * cot1 + p2 * cot2 + fsr_pz[:, None, None]
    mz2 = tot_e ** 2 - px ** 2 - py ** 2 - pz ** 2
    mz = np.sqrt(np.maximum(mz2, 0.0))

    log_c1 = _log_lepton_constraint(p1, pt1, err1, dscb1)
    log_c2 = _log_lepton_constraint(p2, pt2, err2, dscb2)
    log_shape = np.log(np.maximum(shape_fn(mz), 1e-10))

    nll = -(log_c1 + log_c2 + log_shape)  # (N, G, G); minimize
    nll = np.where(np.isfinite(nll), nll, np.inf)

    best = np.argmin(nll.reshape(n, -1), axis=1)
    gi, gj = np.unravel_index(best, (n_grid, n_grid))
    rows = np.arange(n)

    # POST-FIT uncertainty on m(ll) from the profile likelihood: the spread of
    # m(ll) over the delta-NLL <= 0.5 region. This is what the mass error of the
    # refit Z has to be built from -- propagating the raw lepton resolutions
    # would ignore the constraint, which is the whole point of the refit.
    #
    # Only m(ll) is returned, not per-lepton errors: profiling one lepton out
    # lets the other compensate along the constraint, so each lepton's profiled
    # width comes back equal to its input resolution and carries no information.
    sigma_mz, _sat = _profile_sigma(mz.reshape(n, -1), nll.reshape(n, -1))

    return cand1[rows, gi], cand2[rows, gj], sigma_mz


def _profile_sigma(values, nll_profile):
    """1-sigma half-width of ``values`` over the delta-NLL <= 0.5 region.

    ``values`` and ``nll_profile`` are (N, K) over whatever grid was scanned;
    the K axis is treated as an unordered set of points, so this works for the
    flattened 2-D grid as well as for a 1-D profile.

    Returns (sigma, saturated). ``saturated`` marks events whose 1-sigma region
    reaches both extremes of the scanned range -- there the scan is too narrow
    to contain the interval and the width is a LOWER BOUND, not a measurement.
    Where the region collapses onto a single grid point the width is below the
    grid resolution and half a step is used instead.
    """
    n, k = values.shape
    if k < 2:
        return np.zeros(n), np.ones(n, dtype=bool)
    delta = nll_profile - np.min(nll_profile, axis=1)[:, None]
    inside = delta <= 0.5
    masked = np.where(inside, values, np.nan)
    lo = np.nanmin(masked, axis=1)
    hi = np.nanmax(masked, axis=1)
    sigma = 0.5 * (hi - lo)

    full_lo = np.nanmin(values, axis=1)
    full_hi = np.nanmax(values, axis=1)
    # Half a grid step when only the minimum survives (width below resolution).
    half_step = 0.5 * (full_hi - full_lo) / (k - 1)
    sigma = np.where(sigma > 0, sigma, half_step)
    saturated = (lo <= full_lo) & (hi >= full_hi)
    return np.where(np.isfinite(sigma), sigma, 0.0), saturated


def _fsr_cartesian(pt, eta, phi):
    """Cartesian (e, px, py, pz) of a massless FSR photon; pt<=0 -> all zeros."""
    pt = np.asarray(pt, dtype=np.float64)
    eta = np.asarray(eta, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    has = np.isfinite(pt) & (pt > 0) & np.isfinite(eta) & np.isfinite(phi)
    pt = np.where(has, pt, 0.0)
    eta = np.where(has, eta, 0.0)
    phi = np.where(has, phi, 0.0)
    theta = _theta(eta)
    sin = np.sin(theta)
    e = np.where(has, pt / sin, 0.0)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = np.where(has, pt * np.cos(theta) / sin, 0.0)
    return e, px, py, pz


def _prepare_dscb(lepton_constraint, params, n, label):
    """Normalize a user-supplied DSCB parameter set to a 6-tuple of (N,) arrays.

    Returns None for the Gaussian model. Raising here rather than falling back
    to the Gaussian is deliberate: a DSCB run that quietly reverted to the
    nominal constraint would look like a null result instead of a broken job.
    """
    if lepton_constraint != LEPTON_CONSTRAINT_DSCB:
        if params is not None:
            raise ValueError(
                f"{label} was given but lepton_constraint={lepton_constraint!r}; "
                f"pass lepton_constraint={LEPTON_CONSTRAINT_DSCB!r} to use it"
            )
        return None
    if params is None:
        raise ValueError(
            f"lepton_constraint={LEPTON_CONSTRAINT_DSCB!r} requires {label}"
        )
    if isinstance(params, dict):
        missing = [k for k in DSCB_PARAM_NAMES if k not in params]
        if missing:
            raise ValueError(f"{label} is missing DSCB parameters {missing}")
        values = [params[k] for k in DSCB_PARAM_NAMES]
    else:
        values = list(params)
        if len(values) != len(DSCB_PARAM_NAMES):
            raise ValueError(
                f"{label} must have {len(DSCB_PARAM_NAMES)} entries "
                f"({', '.join(DSCB_PARAM_NAMES)}), got {len(values)}"
            )
    out = []
    for name, value in zip(DSCB_PARAM_NAMES, values):
        arr = np.asarray(value, dtype=np.float64)
        if arr.ndim == 0:
            arr = np.full(n, float(arr))
        elif arr.shape != (n,):
            raise ValueError(
                f"{label}[{name!r}] has shape {arr.shape}, expected () or ({n},)"
            )
        out.append(arr)
    return tuple(out)


def refit_lepton_pts(
    pt1, eta1, phi1, m1, err1,
    pt2, eta2, phi2, m2, err2,
    fsr1_pt=None, fsr1_eta=None, fsr1_phi=None,
    fsr2_pt=None, fsr2_eta=None, fsr2_phi=None,
    n_grid=DEFAULT_N_GRID,
    n_sigma=DEFAULT_N_SIGMA,
    pt_floor=DEFAULT_PT_FLOOR,
    params=None,
    event_batch=None,
    n_refine=DEFAULT_N_REFINE,
    return_errors=False,
    edge_rescue_rounds=EDGE_RESCUE_MAX_ROUNDS,
    n_extra_refine=DEFAULT_N_EXTRA_REFINE,
    lineshape=DEFAULT_LINESHAPE,
    channel=None,
    lepton_constraint=DEFAULT_LEPTON_CONSTRAINT,
    dscb_params1=None,
    dscb_params2=None,
):
    """Refit the two lepton pT under the Z mass constraint.

    All inputs are 1-D array-likes of equal length (one entry per event).
    ``fsr{1,2}_{pt,eta,phi}`` are optional per-lepton FSR photons (muon channel);
    pass pt<=0 / None where a lepton has no FSR photon. The FSR photons are folded
    into m(ll) for the lineshape constraint (held at their reco values); the
    returned pT are the lepton *core* pT only -- the caller adds the FSR photon
    4-momentum back onto each lepton to build the dressed refit Z.

    ``lineshape`` selects the Z-constraint model from :data:`LINESHAPE_MODELS`
    (default :data:`LEGACY_LINESHAPE`, the nominal empirical CB+3Gauss);
    ``params`` overrides that model's parameters and defaults to the model's own.
    ``channel`` ("ee"/"mumu") picks the per-channel parameter set -- all rows
    passed in one call are assumed to be that channel, so a channel-split refit
    means calling this once per channel on the corresponding event subset.
    The shape is used unnormalized on purpose: a constant factor only shifts the
    NLL by a constant and cannot move the per-event minimum.

    ``lepton_constraint`` selects the per-lepton pT constraint:
    ``"gaussian"`` (default, the nominal analysis) or ``"dscb"``. For
    ``"dscb"`` both ``dscb_params1`` and ``dscb_params2`` are required and are
    per-lepton parameter sets -- either a dict keyed by
    :data:`DSCB_PARAM_NAMES` or a 6-tuple in that order. Each entry may be a
    scalar (same response for every row) or a length-N array (the usual case:
    the per-event lookup from
    :class:`~higgs_dna.tools.lepton_pt_dscb.LeptonPtDSCBCalibration`).

    ``n_refine`` > 1 adds a second pass that re-scans one coarse step either
    side of the coarse winner, giving an effective step of
    ``(2*n_sigma/(n_grid-1)) * 2/(n_refine-1)`` sigma. This is what removes the
    bimodal m(ll)^refit the single coarse pass produces; set ``n_refine=1`` (with
    ``n_grid=PRODUCTION_N_GRID``) to reproduce the stored production values.

    ``n_extra_refine`` further zooms (each 5x finer, see
    :data:`DEFAULT_N_EXTRA_REFINE`) run after the refinement and the edge
    rescue and bring the result to the continuous minimum; 0 disables them.

    Returns ``(pt1_refit, pt2_refit)`` as float64 NumPy arrays. Events with
    non-finite kinematics or non-positive errors fall back to the reconstructed pT
    (errors default to 1% of pT when missing).

    With ``return_errors=True`` the return becomes
    ``(pt1_refit, pt2_refit, m_ll_err)``. ``m_ll_err`` is the POST-FIT
    uncertainty on m(ll): the delta-NLL <= 0.5 half-width of the profile
    likelihood, i.e. what the Z constraint leaves rather than what the leptons
    alone would give. It is taken from the coarse scan, which spans +-n_sigma
    and so contains the interval. Rows that were not refitted get NaN.
    """
    shape_fn = build_lineshape(lineshape, params, channel=channel)
    if lepton_constraint not in LEPTON_CONSTRAINTS:
        raise ValueError(
            f"unknown lepton_constraint {lepton_constraint!r}; "
            f"choose from {LEPTON_CONSTRAINTS}"
        )
    # Keep the (N, G, G) working set flat whatever the grid size, so a finer
    # scan costs time rather than memory.
    if event_batch is None:
        event_batch = max(
            2000,
            int(_EVENT_BATCH * (DEFAULT_N_GRID / float(max(n_grid, 1))) ** 2),
        )
    pt1 = np.asarray(pt1, dtype=np.float64)
    pt2 = np.asarray(pt2, dtype=np.float64)
    eta1 = np.asarray(eta1, dtype=np.float64)
    eta2 = np.asarray(eta2, dtype=np.float64)
    phi1 = np.asarray(phi1, dtype=np.float64)
    phi2 = np.asarray(phi2, dtype=np.float64)
    m1 = np.asarray(m1, dtype=np.float64)
    m2 = np.asarray(m2, dtype=np.float64)
    err1 = np.asarray(err1, dtype=np.float64)
    err2 = np.asarray(err2, dtype=np.float64)

    n = pt1.shape[0]
    dscb1 = _prepare_dscb(lepton_constraint, dscb_params1, n, "dscb_params1")
    dscb2 = _prepare_dscb(lepton_constraint, dscb_params2, n, "dscb_params2")
    if n == 0:
        if return_errors:
            return pt1.copy(), pt2.copy(), np.zeros(0)
        return pt1.copy(), pt2.copy()

    # Total (fixed) FSR 4-momentum per event = sum over the two leptons' FSR photons.
    zeros = np.zeros(n)
    e1f, px1f, py1f, pz1f = _fsr_cartesian(
        zeros if fsr1_pt is None else fsr1_pt,
        zeros if fsr1_eta is None else fsr1_eta,
        zeros if fsr1_phi is None else fsr1_phi,
    )
    e2f, px2f, py2f, pz2f = _fsr_cartesian(
        zeros if fsr2_pt is None else fsr2_pt,
        zeros if fsr2_eta is None else fsr2_eta,
        zeros if fsr2_phi is None else fsr2_phi,
    )
    fsr_e = e1f + e2f
    fsr_px = px1f + px2f
    fsr_py = py1f + py2f
    fsr_pz = pz1f + pz2f

    # Sanitize errors: fall back to 1% of pT (floored) when missing/non-positive.
    err1 = np.where(np.isfinite(err1) & (err1 > 0), err1, np.maximum(0.01 * pt1, 1e-3))
    err2 = np.where(np.isfinite(err2) & (err2 > 0), err2, np.maximum(0.01 * pt2, 1e-3))

    pt1_refit = pt1.copy()
    pt2_refit = pt2.copy()
    # Unrefitted rows keep the input resolution; their m(ll) error is undefined.
    mll_err = np.full(n, np.nan)
    n_rescued = 0

    # Only refit rows with finite kinematics; others keep the reco pT.
    valid = (
        np.isfinite(pt1) & np.isfinite(pt2)
        & np.isfinite(eta1) & np.isfinite(eta2)
        & np.isfinite(phi1) & np.isfinite(phi2)
        & (pt1 > 0) & (pt2 > 0)
    )
    idx = np.nonzero(valid)[0]
    if idx.size == 0:
        if return_errors:
            return pt1_refit, pt2_refit, mll_err
        return pt1_refit, pt2_refit

    for start in range(0, idx.size, event_batch):
        sl = idx[start:start + event_batch]
        common = (fsr_e[sl], fsr_px[sl], fsr_py[sl], fsr_pz[sl],
                  n_grid, n_sigma, pt_floor, shape_fn)
        sub1 = None if dscb1 is None else tuple(p[sl] for p in dscb1)
        sub2 = None if dscb2 is None else tuple(p[sl] for p in dscb2)
        r1, r2, smz = _refit_batch(
            pt1[sl], eta1[sl], phi1[sl], m1[sl], err1[sl],
            pt2[sl], eta2[sl], phi2[sl], m2[sl], err2[sl], *common,
            dscb1=sub1, dscb2=sub2,
        )
        if n_refine and n_refine > 1:
            # One coarse step either side of the winner, sampled n_refine times.
            # Effective step = (2*n_sigma/(n_grid-1)) * (2/(n_refine-1)) sigma.
            step1 = (2.0 * n_sigma / (n_grid - 1)) * err1[sl]
            step2 = (2.0 * n_sigma / (n_grid - 1)) * err2[sl]
            # Keep the COARSE m(ll) error: the refined scan spans only +-1
            # coarse step, which is narrower than the 1-sigma interval and would
            # truncate it. The refinement is for the central value.
            r1, r2, _smz_fine = _refit_batch(
                pt1[sl], eta1[sl], phi1[sl], m1[sl], err1[sl],
                pt2[sl], eta2[sl], phi2[sl], m2[sl], err2[sl], *common[:4],
                n_refine, n_sigma, pt_floor, shape_fn,
                centre1=r1, centre2=r2, half1=step1, half2=step2,
                dscb1=sub1, dscb2=sub2,
            )
        # --- edge rescue -------------------------------------------------
        # Re-scan only the events pinned to the scan boundary, with a window
        # grown geometrically. Everything else keeps the nominal scan exactly,
        # so this cannot move an event that already had an interior optimum.
        if edge_rescue_rounds:
            width = n_sigma
            coarse_step_sigma = 2.0 * n_sigma / (n_grid - 1)
            for _round in range(int(edge_rescue_rounds)):
                on_edge = (
                    (np.abs(r1 - pt1[sl]) >= width * err1[sl] * _EDGE_TOL)
                    | (np.abs(r2 - pt2[sl]) >= width * err2[sl] * _EDGE_TOL)
                )
                # A leg already sitting on the pT floor cannot move outward on
                # the low side; rescuing it would loop forever.
                on_edge &= ~((r1 <= pt_floor) | (r2 <= pt_floor))
                if not np.any(on_edge):
                    break
                width *= EDGE_RESCUE_GROWTH
                # Grow the POINT COUNT with the window so the step stays put.
                # Widening alone would coarsen the scan and the rescue would
                # return a worse optimum than the one it replaced -- measured:
                # agreement with an equal-step wide reference DROPPED from 52%
                # to 22% on the rescued events before this was fixed.
                n_grid_e = int(round(2.0 * width / coarse_step_sigma)) + 1
                idx_e = np.nonzero(on_edge)[0]
                sub1_e = None if sub1 is None else tuple(p[idx_e] for p in sub1)
                sub2_e = None if sub2 is None else tuple(p[idx_e] for p in sub2)
                e1, e2, esmz = _refit_batch(
                    pt1[sl][idx_e], eta1[sl][idx_e], phi1[sl][idx_e],
                    m1[sl][idx_e], err1[sl][idx_e],
                    pt2[sl][idx_e], eta2[sl][idx_e], phi2[sl][idx_e],
                    m2[sl][idx_e], err2[sl][idx_e],
                    fsr_e[sl][idx_e], fsr_px[sl][idx_e],
                    fsr_py[sl][idx_e], fsr_pz[sl][idx_e],
                    n_grid_e, width, pt_floor, shape_fn,
                    dscb1=sub1_e, dscb2=sub2_e,
                )
                if n_refine and n_refine > 1:
                    step1_e = (2.0 * width / (n_grid_e - 1)) * err1[sl][idx_e]
                    step2_e = (2.0 * width / (n_grid_e - 1)) * err2[sl][idx_e]
                    e1, e2, _ = _refit_batch(
                        pt1[sl][idx_e], eta1[sl][idx_e], phi1[sl][idx_e],
                        m1[sl][idx_e], err1[sl][idx_e],
                        pt2[sl][idx_e], eta2[sl][idx_e], phi2[sl][idx_e],
                        m2[sl][idx_e], err2[sl][idx_e],
                        fsr_e[sl][idx_e], fsr_px[sl][idx_e],
                        fsr_py[sl][idx_e], fsr_pz[sl][idx_e],
                        n_refine, width, pt_floor, shape_fn,
                        centre1=e1, centre2=e2, half1=step1_e, half2=step2_e,
                        dscb1=sub1_e, dscb2=sub2_e,
                    )
                ok_e = np.isfinite(e1) & np.isfinite(e2)
                r1[idx_e[ok_e]] = e1[ok_e]
                r2[idx_e[ok_e]] = e2[ok_e]
                smz[idx_e[ok_e]] = esmz[ok_e]
                n_rescued += int(ok_e.sum())

        # --- extra zooms to the continuous minimum ---------------------------
        # Every event (rescued or not) ends the scan above at the same
        # 0.08-sigma step, so one schedule serves all of them.
        if n_refine and n_refine > 1 and n_extra_refine:
            step_sigma = (2.0 * n_sigma / (n_grid - 1)) * (2.0 / (n_refine - 1))
            for _zoom in range(int(n_extra_refine)):
                r1, r2, _ = _refit_batch(
                    pt1[sl], eta1[sl], phi1[sl], m1[sl], err1[sl],
                    pt2[sl], eta2[sl], phi2[sl], m2[sl], err2[sl], *common[:4],
                    n_refine, n_sigma, pt_floor, shape_fn,
                    centre1=r1, centre2=r2,
                    half1=step_sigma * err1[sl], half2=step_sigma * err2[sl],
                    dscb1=sub1, dscb2=sub2,
                )
                step_sigma *= 2.0 / (n_refine - 1)

        good = np.isfinite(r1) & np.isfinite(r2)
        pt1_refit[sl[good]] = r1[good]
        pt2_refit[sl[good]] = r2[good]
        mll_err[sl[good]] = smz[good]

    if n_rescued and edge_rescue_rounds:
        logger.info(
            "zmass refit: %d/%d legs-pairs were pinned to the +-%.1f sigma scan "
            "edge and re-scanned with a widened window (up to %.0f sigma).",
            n_rescued, idx.size, n_sigma,
            n_sigma * EDGE_RESCUE_GROWTH ** int(edge_rescue_rounds),
        )

    if return_errors:
        return pt1_refit, pt2_refit, mll_err
    return pt1_refit, pt2_refit
