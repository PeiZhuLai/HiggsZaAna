"""Z mass-constraint kinematic refit for H->Z(ll)gamma.

Vectorized NumPy port of the per-lepton pT likelihood fit used in nano2pico
(``src/KinZfitter.cpp``) and HiggsZGammaAna/hzgml (``scripts/z_refit.py``).

For each event the two lepton pT are floated (eta/phi/mass held fixed) to
maximize the likelihood

    L = Gauss(pt1; pt1_reco, sigma1)
      * Gauss(pt2; pt2_reco, sigma2)
      * Shape(m_ll)

where ``Shape(m_ll)`` is the empirical Z lineshape (Crystal Ball + 3 Gaussians,
identical electron/muon parameters in z_refit.py).  The negative log-likelihood

    NLL = -[ logGauss1 + logGauss2 + log(Shape(m_ll)) ]

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

import numpy as np

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

# Grid-search configuration (faithful to z_refit's bounds).
DEFAULT_N_GRID = 11      # points scanned per lepton (z_refit uses 11)
DEFAULT_N_SIGMA = 2.0    # half-window in units of sigma (z_refit bound: reco +- 2 sigma)
DEFAULT_PT_FLOOR = 5.0   # GeV, minimum allowed pT (z_refit bound)
_EVENT_BATCH = 20000     # process events in slices to cap the (N, G, G) peak memory


def _theta(eta):
    return 2.0 * np.arctan(np.exp(-eta))


def _crystal_ball(x, mean, sigma, alpha, n):
    """Crystal Ball value (unnormalized core, matching z_refit.crystal_ball_pdf).

    Assumes ``alpha > 0`` (true for LINESHAPE_PARAMS).
    """
    t = (x - mean) / sigma
    a = (n / alpha) ** n * np.exp(-0.5 * alpha * alpha)
    b = n / alpha - alpha
    core = np.exp(-0.5 * t * t)
    tail = a * np.power(np.maximum(b - t, 1e-9), -n)
    return np.where(t >= -alpha, core, tail)


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


def _refit_batch(pt1, eta1, phi1, m1, err1, pt2, eta2, phi2, m2, err2,
                 fsr_e, fsr_px, fsr_py, fsr_pz,
                 n_grid, n_sigma, pt_floor, params):
    """Core grid search for one batch of events. Returns (pt1_refit, pt2_refit).

    ``fsr_{e,px,py,pz}`` are the per-event SUM of the (muon-channel) FSR photon
    4-momenta to fold into the m(ll) used by the Z-lineshape constraint. They are
    held fixed at their reconstructed values (only the two lepton pT are floated),
    so FSR enters the mass constraint but is not itself re-fitted. Pass zeros for
    no FSR.
    """
    n = pt1.shape[0]
    frac = np.linspace(-n_sigma, n_sigma, n_grid)  # (G,)

    cand1 = np.maximum(pt1[:, None] + frac[None, :] * err1[:, None], pt_floor)  # (N, G)
    cand2 = np.maximum(pt2[:, None] + frac[None, :] * err2[:, None], pt_floor)  # (N, G)

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

    log_gauss1 = -0.5 * ((p1 - pt1[:, None, None]) / err1[:, None, None]) ** 2
    log_gauss2 = -0.5 * ((p2 - pt2[:, None, None]) / err2[:, None, None]) ** 2
    log_shape = np.log(np.maximum(_z_lineshape(mz, params), 1e-10))

    nll = -(log_gauss1 + log_gauss2 + log_shape)  # (N, G, G); minimize
    nll = np.where(np.isfinite(nll), nll, np.inf)

    best = np.argmin(nll.reshape(n, -1), axis=1)
    gi, gj = np.unravel_index(best, (n_grid, n_grid))
    rows = np.arange(n)
    return cand1[rows, gi], cand2[rows, gj]


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


def refit_lepton_pts(
    pt1, eta1, phi1, m1, err1,
    pt2, eta2, phi2, m2, err2,
    fsr1_pt=None, fsr1_eta=None, fsr1_phi=None,
    fsr2_pt=None, fsr2_eta=None, fsr2_phi=None,
    n_grid=DEFAULT_N_GRID,
    n_sigma=DEFAULT_N_SIGMA,
    pt_floor=DEFAULT_PT_FLOOR,
    params=LINESHAPE_PARAMS,
    event_batch=_EVENT_BATCH,
):
    """Refit the two lepton pT under the Z mass constraint.

    All inputs are 1-D array-likes of equal length (one entry per event).
    ``fsr{1,2}_{pt,eta,phi}`` are optional per-lepton FSR photons (muon channel);
    pass pt<=0 / None where a lepton has no FSR photon. The FSR photons are folded
    into m(ll) for the lineshape constraint (held at their reco values); the
    returned pT are the lepton *core* pT only -- the caller adds the FSR photon
    4-momentum back onto each lepton to build the dressed refit Z.

    Returns ``(pt1_refit, pt2_refit)`` as float64 NumPy arrays. Events with
    non-finite kinematics or non-positive errors fall back to the reconstructed pT
    (errors default to 1% of pT when missing).
    """
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
    if n == 0:
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

    # Only refit rows with finite kinematics; others keep the reco pT.
    valid = (
        np.isfinite(pt1) & np.isfinite(pt2)
        & np.isfinite(eta1) & np.isfinite(eta2)
        & np.isfinite(phi1) & np.isfinite(phi2)
        & (pt1 > 0) & (pt2 > 0)
    )
    idx = np.nonzero(valid)[0]
    if idx.size == 0:
        return pt1_refit, pt2_refit

    for start in range(0, idx.size, event_batch):
        sl = idx[start:start + event_batch]
        r1, r2 = _refit_batch(
            pt1[sl], eta1[sl], phi1[sl], m1[sl], err1[sl],
            pt2[sl], eta2[sl], phi2[sl], m2[sl], err2[sl],
            fsr_e[sl], fsr_px[sl], fsr_py[sl], fsr_pz[sl],
            n_grid, n_sigma, pt_floor, params,
        )
        good = np.isfinite(r1) & np.isfinite(r2)
        pt1_refit[sl[good]] = r1[good]
        pt2_refit[sl[good]] = r2[good]

    return pt1_refit, pt2_refit
