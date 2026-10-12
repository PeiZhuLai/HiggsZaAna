# Run-2 (2016preVFP/2016postVFP/2017/2018 UL) muon momentum corrections = Rochester
# (RoccoR). The Run-3 muon POG "muonscarekit" (muon_systematics.py:muon_pt_scare)
# only ships Run-3 JSONs, so Run-2 muons would otherwise get NO momentum
# correction. This module is a pure numpy/awkward port of the official RoccoR.cc
# (systematics/roccor/RoccoR.cc + RoccoR20XXUL.txt), central set only (s=0, m=0);
# the RoccoR uncertainty replicas are negligible vs other HZg systematics (see
# systematics/roccor/README) and are deferred.
#
# Correction (matches PhysicsTools/NanoAODTools muonScaleResProducer.py Run-2 path):
#   data              : pt_corr = pt * kScaleDT(Q, pt, eta, phi)
#   MC (gen-matched)  : pt_corr = pt * kSpreadMC(Q, pt, eta, phi, genPt)   [recommended]
#   MC (not matched)  : pt_corr = pt * kSmearMC(Q, pt, eta, phi, nTrkLayers, u)
#   pt > 200 GeV      : pt_corr = pt (data and MC; nano2pico mu_producer.cpp:126-136)

import os
import logging

import numpy as np
import awkward as ak
from scipy.special import erfinv

logger = logging.getLogger(__name__)

_PI = 3.14159  # RoccoR uses this truncated value for DPHI/MPHI, keep identical

# Muons with pt above this (GeV) are not Rochester-corrected (nano2pico
# mu_producer.cpp:126-136); see muon_pt_rochester_run2.
ROCCOR_MAX_PT = 200.0

# year -> RoccoR UL text file (staged in systematics/roccor/)
_ROCCOR_FILES = {
    "2016preVFP": "RoccoR2016aUL.txt",
    "2016postVFP": "RoccoR2016bUL.txt",
    "2017": "RoccoR2017UL.txt",
    "2018": "RoccoR2018UL.txt",
}

_CACHE = {}


# --------------------------------------------------------------------------- #
# Parser: read the central set (sys=0, mem=0) of a RoccoR .txt into arrays.
# --------------------------------------------------------------------------- #
def _parse_roccor(path):
    NPHI = None
    DPHI = None
    NETA = None            # CP eta bins (signed)
    etabin = None          # CP eta boundaries (len NETA+1)
    RMIN = RTRK = RETA = None
    BETA = None            # resolution eta boundaries (len RETA+1, abs eta)

    # CP scale params, [type][etaBin][phiBin]; type 0=MC, 1=DT
    M = A = X = None
    kRes = None            # [type][resEtaBin]
    rsPar = None           # [3][resEtaBin][trkBin]
    cb_s = cb_a = cb_n = None  # [resEtaBin][trkBin]
    nTrk = None            # [type][resEtaBin][trkBin+1]

    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if not s.strip():
                continue
            head = s[:4]
            tok = s.split()
            if s[:7] == "VERSION":
                continue
            if head == "NSET":
                continue
            if head == "NMEM":
                continue
            if head == "TVAR":
                continue
            if head == "RMIN":
                RMIN = int(tok[1]); continue
            if head == "RTRK":
                RTRK = int(tok[1]); continue
            if head == "RETA":
                RETA = int(tok[1]); BETA = [float(x) for x in tok[2:2 + RETA + 1]]; continue
            if head == "CPHI":
                NPHI = int(tok[1]); DPHI = 2 * _PI / NPHI; continue
            if head == "CETA":
                NETA = int(tok[1]); etabin = [float(x) for x in tok[2:2 + NETA + 1]]; continue

            # data line: sys mem tag ...  -- keep only central set/member
            sys_, mem_, tag = int(tok[0]), int(tok[1]), tok[2]
            if sys_ != 0 or mem_ != 0:
                continue

            # lazy-allocate now that header dims are known
            if M is None:
                M = np.ones((2, NETA, NPHI))
                A = np.zeros((2, NETA, NPHI))
                X = np.zeros((2, NETA, NPHI))
                kRes = np.ones((2, RETA))
                rsPar = np.zeros((3, RETA, RTRK))
                cb_s = np.ones((RETA, RTRK))
                cb_a = np.full((RETA, RTRK), 10.0)
                cb_n = np.full((RETA, RTRK), 10.0)
                nTrk = np.zeros((2, RETA, RTRK + 1))

            if tag == "C":
                typ, var, b = int(tok[3]), int(tok[4]), int(tok[5])
                vals = [float(x) for x in tok[6:6 + NPHI]]
                for i, v in enumerate(vals):
                    if var == 0:
                        M[typ, b, i] = 1.0 + v / 100.0
                    elif var == 1:
                        A[typ, b, i] = v / 100.0
            elif tag == "X":
                typ = int(tok[3]); dKdX = float(tok[4])
                X[typ, :, :] = dKdX
            elif tag == "F":
                typ = int(tok[3])
                vals = [float(x) for x in tok[4:4 + RETA]]
                for i, v in enumerate(vals):
                    kRes[typ, i] = v
            elif tag == "R":
                var, b = int(tok[3]), int(tok[4])
                vals = [float(x) for x in tok[5:5 + RTRK]]
                for i, v in enumerate(vals):
                    if var == 0:
                        rsPar[0, b, i] = v
                    elif var == 1:
                        rsPar[1, b, i] = v
                    elif var == 2:
                        rsPar[2, b, i] = v / 100.0
                    elif var == 3:
                        cb_s[b, i] = v
                    elif var == 4:
                        cb_a[b, i] = v
                    elif var == 5:
                        cb_n[b, i] = v
            elif tag == "T":
                typ, b = int(tok[3]), int(tok[4])
                vals = [float(x) for x in tok[5:5 + RTRK + 1]]
                for i, v in enumerate(vals):
                    nTrk[typ, b, i] = v

    # Precompute CrystalBall init constants per (resEtaBin, trkBin) -- see RoccoR.h
    _cberr = np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore")
    _cberr.__enter__()
    fa = np.abs(cb_a)
    ex = np.exp(-fa * fa / 2.0)
    A_ = np.power(cb_n / fa, cb_n) * ex
    C1 = cb_n / fa / (cb_n - 1.0) * ex
    D1 = 2.0 * np.sqrt(_PI / 2.0) * np.vectorize(_erf)(fa / np.sqrt(2.0))
    cb = {
        "m": np.zeros_like(cb_s),  # CrystalBall.m defaults to 0
        "s": cb_s, "a": cb_a, "n": cb_n,
        "B": cb_n / fa - fa,
        "N": 1.0 / cb_s / (D1 + 2 * C1),
        "k": 1.0 / (cb_n - 1.0),
    }
    cb["C"] = (D1 + 2 * C1) / C1
    cb["D"] = (D1 + 2 * C1) / 2.0
    cb["NC"] = cb["N"] * cb_s * C1
    cb["Ns"] = cb["N"] * cb_s
    cb["F"] = 1.0 - fa * fa / cb_n
    cb["G"] = cb_s * cb_n / fa
    cb["cdfMa"] = _cb_cdf_scalarbins(cb, cb["m"] - cb_a * cb_s)
    cb["cdfPa"] = _cb_cdf_scalarbins(cb, cb["m"] + cb_a * cb_s)
    _cberr.__exit__(None, None, None)

    return dict(
        NPHI=NPHI, DPHI=DPHI, NETA=NETA, etabin=np.array(etabin),
        RETA=RETA, RTRK=RTRK, RMIN=RMIN, BETA=np.array(BETA),
        M=M, A=A, X=X, kRes=kRes, rsPar=rsPar, nTrk=nTrk, cb=cb,
    )


def _erf(x):
    # scalar erf via erfinv's sibling; use math.erf to avoid scipy import cost here
    import math
    return math.erf(x)


def _cb_cdf_scalarbins(cb, x):
    # cdf evaluated per-bin (x has bin shape); mirrors CrystalBall::cdf. The
    # unused/default (a=n=10) bins produce harmless inf/nan in the non-selected
    # np.where branches -> silence those; the tail branches are masked out.
    m, s, a = cb["m"], cb["s"], cb["a"]
    d = (x - m) / s
    with np.errstate(divide="ignore", invalid="ignore"):
        out = cb["Ns"] * (cb["D"] - np.sqrt(_PI / 2.0) * np.vectorize(_erf)(-d / np.sqrt(2.0)))
        lo = d < -a
        hi = d > a
        out = np.where(lo, cb["NC"] / np.power(cb["F"] - s * d / cb["G"], cb["n"] - 1), out)
        out = np.where(hi, cb["NC"] * (cb["C"] - np.power(cb["F"] + s * d / cb["G"], 1 - cb["n"])), out)
    return out


def _load_roccor(year):
    if year not in _CACHE:
        fname = _ROCCOR_FILES.get(year)
        if fname is None:
            raise ValueError(f"RoccoR (Run-2) not available for year '{year}'")
        path = os.path.join(os.path.dirname(__file__), "roccor", fname)
        _CACHE[year] = _parse_roccor(path)
    return _CACHE[year]


# --------------------------------------------------------------------------- #
# Vectorized bin lookups (flat numpy arrays in)
# --------------------------------------------------------------------------- #
def _cp_eta_bin(eta, rc):
    # RoccoR::etaBin -- signed eta against CETA boundaries; returns 0..NETA-1
    idx = np.searchsorted(rc["etabin"][1:-1], eta, side="right")
    return np.clip(idx, 0, rc["NETA"] - 1)


def _phi_bin(phi, rc):
    ibin = ((phi - (-_PI)) / rc["DPHI"]).astype(np.int64)
    return np.clip(ibin, 0, rc["NPHI"] - 1)


def _res_eta_bin(eta, rc):
    # RocRes::etaBin -- |eta| against BETA boundaries; returns 0..RETA-1
    idx = np.searchsorted(rc["BETA"][1:-1], np.abs(eta), side="right")
    return np.clip(idx, 0, rc["RETA"] - 1)


def _k_scale(typ, Q, pt, eta, phi, rc):
    H = _cp_eta_bin(eta, rc)
    F = _phi_bin(phi, rc)
    M = rc["M"][typ, H, F]
    A = rc["A"][typ, H, F]
    X = rc["X"][typ, H, F]
    return 1.0 / (M + Q * A * pt + X / pt)


def kScaleDT(Q, pt, eta, phi, rc):
    return _k_scale(1, Q, pt, eta, phi, rc)


def kSpreadMC(Q, pt, eta, phi, gpt, rc):
    k = _k_scale(0, Q, pt, eta, phi, rc)          # kScaleMC
    Hr = _res_eta_bin(eta, rc)
    x = gpt / (k * pt)
    kData = rc["kRes"][1, Hr]
    kMC = rc["kRes"][0, Hr]
    kspread = x / (1.0 + (x - 1.0) * kData / kMC)  # RocRes::kSpread (3-arg)
    return k * kspread


def _cb_invcdf(cb, H, F, u):
    # CrystalBall::invcdf per-muon; cb params selected by (H, F)
    m = cb["m"][H, F]; s = cb["s"][H, F]
    NC = cb["NC"][H, F]; k = cb["k"][H, F]; C = cb["C"][H, F]
    G = cb["G"][H, F]; Fc = cb["F"][H, F]; Ns = cb["Ns"][H, F]; D = cb["D"][H, F]
    cdfMa = cb["cdfMa"][H, F]; cdfPa = cb["cdfPa"][H, F]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        mid = m - np.sqrt(2.0) * s * erfinv((D - u / Ns) / np.sqrt(_PI / 2.0))
        lo = m + G * (Fc - np.power(NC / np.where(u <= 0, 1e-300, u), k))
        hi = m - G * (Fc - np.power(np.clip(C - u / NC, 1e-300, None), -k))
    out = np.where(u < cdfMa, lo, mid)
    out = np.where(u > cdfPa, hi, out)
    return out


def _sigma(pt, Hr, F, rc):
    dpt = pt - 45.0
    return rc["rsPar"][0, Hr, F] + rc["rsPar"][1, Hr, F] * dpt + rc["rsPar"][2, Hr, F] * dpt * dpt


def kSmearMC(Q, pt, eta, phi, nl, u, rc):
    k = _k_scale(0, Q, pt, eta, phi, rc)          # kScaleMC
    Hr = _res_eta_bin(eta, rc)
    F = np.clip(nl.astype(np.int64) - rc["RMIN"], 0, rc["RTRK"] - 1)
    d = rc["kRes"][1, Hr]
    m = rc["kRes"][0, Hr]
    ptk = k * pt
    sig = _sigma(ptk, Hr, F, rc)
    fac = np.where(d > m, np.sqrt(np.clip(d * d - m * m, 0, None)), 0.0)
    x = fac * sig * _cb_invcdf(rc["cb"], Hr, F, u)
    kextra = np.where(x <= -1, 1.0, 1.0 / (1.0 + x))
    return k * kextra


def _det_uniform(eta, phi, pt):
    # deterministic, reproducible u in [0,1) for the non-gen-matched smear
    h = np.sin(phi * 12.9898 + eta * 78.233 + pt * 37.719) * 43758.5453
    return h - np.floor(h)


# --------------------------------------------------------------------------- #
# Main entry -- mirrors muon_systematics.muon_pt_scare(is_correction=True)
# --------------------------------------------------------------------------- #
def muon_pt_rochester_run2(pt, events, year="2018", unc_type=None, is_correction=True):
    """Apply Run-2 Rochester (RoccoR) muon momentum correction; sets events.Muon.pt."""
    if not is_correction:
        # RoccoR uncertainty replicas are negligible vs other HZg systematics -> no-op.
        return events

    rc = _load_roccor(year)
    is_data = not hasattr(events, "genWeight")

    counts = ak.num(events.Muon.pt)
    muons = ak.flatten(events.Muon)

    Q = ak.to_numpy(muons.charge).astype(np.float64)
    ptn = ak.to_numpy(muons.pt).astype(np.float64)
    eta = ak.to_numpy(muons.eta).astype(np.float64)
    phi = ak.to_numpy(muons.phi).astype(np.float64)

    if is_data:
        sf = kScaleDT(Q, ptn, eta, phi, rc)
        pt_corr = ptn * sf
    else:
        nl = ak.to_numpy(muons.nTrackerLayers).astype(np.float64)
        # gen-matched muon pt (Muon.genPartIdx -> GenPart.pt); -1 => not matched
        gidx = ak.to_numpy(muons.genPartIdx).astype(np.int64)
        gpt_all = ak.to_numpy(ak.flatten(events.GenPart.pt)).astype(np.float64)
        # build per-muon gen offset to index the flattened GenPart array
        gp_counts = ak.to_numpy(ak.num(events.GenPart.pt)).astype(np.int64)
        mu_counts = ak.to_numpy(counts).astype(np.int64)
        gp_offset = np.repeat(np.concatenate([[0], np.cumsum(gp_counts)[:-1]]), mu_counts)
        # `>= 0` is deliberate: GenPart index 0 is a valid match. nano2pico uses
        # `Muon_genPartIdx[imu] > 0` (mu_producer.cpp:147), which drops index 0 --
        # an n2p bug that we do not copy.
        matched = (gidx >= 0) & (gidx < np.repeat(gp_counts, mu_counts))
        gpt = np.where(matched, gpt_all[np.where(matched, gp_offset + gidx, 0)], 0.0)

        sf_spread = kSpreadMC(Q, ptn, eta, phi, np.where(matched, gpt, ptn), rc)
        u = _det_uniform(eta, phi, ptn)
        sf_smear = kSmearMC(Q, ptn, eta, phi, nl, u, rc)
        sf = np.where(matched, sf_spread, sf_smear)
        pt_corr = ptn * sf

    # Muons above 200 GeV are left uncorrected, data and MC alike, as in
    # nano2pico mu_producer.cpp:126-136 (`if(pt>200.f){ ... push_back(pt); continue; }`;
    # the muon POG high-pT recommendation is separate and RoccoR is not meant
    # for it). Strictly greater than 200, so pt == 200 is still corrected. The
    # Run-3 ScaRe path does the same (muon_systematics.py, `pt > 200`).
    pt_corr = np.where(ptn > ROCCOR_MAX_PT, ptn, pt_corr)

    muons = ak.with_field(muons, ak.Array(pt_corr), "pt")
    events["Muon"] = ak.unflatten(muons, counts)
    logger.debug("[ Muon RoccoR ] Run-2 Rochester muon pt correction applied (%s)", year)
    return events
