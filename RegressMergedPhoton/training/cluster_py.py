"""Offline (numpy) reimplementation of the MLPhoton clustering + image builder.

This is a deliberately literal port of the online C++ used to produce the
MLPhoton collection:

    RecoEgamma/EgammaMLPhotonProducers/src/Cluster.cc          (Cluster, makeImage, compute_En)
    RecoEgamma/EgammaMLPhotonProducers/src/MLPhotonProducer.cc (DoPairings, FindNearest,
                                                                calculateLorentzVector)

Why a port and not a rewrite: the network inputs are whatever `makeImage()`
happens to produce, quirks included. If the training preprocessing differs from
the online preprocessing by even a shifted pixel, the retrained model is being
fed a different detector. So the quirks below are reproduced ON PURPOSE and must
NOT be "fixed":

  * makeImage() does `ieta += ieta_min + isize/2` -- it ADDS the minimum where
    subtracting it would be the natural choice. Harmless in the end (the very
    next step subtracts `cxval`, the first crystal's shifted coordinate) but it
    does change nothing only because of that cancellation, so leave it alone.
  * the image is centred on the FIRST crystal in the member list, not on the
    highest-energy one. Member order therefore matters and is preserved exactly:
    combine() appends C2's crystals after C1's.
  * x is clamped into [0, isize-1] while y is first shifted (while y_min < 0)
    then shifted back (if y_max > isize) and only afterwards clamped -- the two
    axes are genuinely treated differently.
  * the classifier is fed the RAW energy image, the regressor the SUM-NORMALISED
    one. (In C++ this is obscured by `img` being a 900x900 vector of which only
    row 0 is ever used.)

Known non-reproducible corner: std::sort is not stable, so clusters with exactly
equal |p|^2 may be ordered differently here (numpy stable sort) than online. With
float energies this is rare; closure_test.py reports any event where it bites.

Usage:
    from cluster_py import build_clusters, MLPhotonModels
    clusters = build_clusters(rh_ieta, rh_iphi, rh_eta, rh_phi, rh_energy)
    models = MLPhotonModels(classifier_onnx, regressor_onnx)
    out = models.predict(clusters, pv_z)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

ISIZE = 30
MATCH_DR = 0.15   # MLPhotonProducer::MATCH_DR
ISO_DR = 0.3      # MLPhotonProducer::ISO_DR (isolation needs PF candidates; not ported)
ECAL_R = 129.0    # calculateLorentzVector


# --------------------------------------------------------------------------
# Cluster
# --------------------------------------------------------------------------
@dataclass
class Cluster:
    """One-to-one with the C++ Cluster class.

    `vec` is the C++ math::RhoEtaPhiVectorF held as cartesian (px, py, pz) so
    that combine() is plain vector addition, exactly as `vec += C.vec` is.
    """

    px: float
    py: float
    pz: float
    ietas: list = field(default_factory=list)
    iphis: list = field(default_factory=list)
    Es: list = field(default_factory=list)
    xcoords: list = field(default_factory=list)
    ycoords: list = field(default_factory=list)
    image: np.ndarray | None = None  # (ISIZE, ISIZE) float32, raw energies

    # -- kinematics -------------------------------------------------------
    @property
    def mag2(self) -> float:
        """C++ vec.Mag2(): squared 3-momentum, the sort key in the pairing loop."""
        return self.px * self.px + self.py * self.py + self.pz * self.pz

    def eta(self) -> float:
        pt = math.hypot(self.px, self.py)
        if pt == 0.0:
            return math.inf if self.pz > 0 else -math.inf
        return math.asinh(self.pz / pt)

    def phi(self) -> float:
        return math.atan2(self.py, self.px)

    def total_energy(self) -> float:
        """C++ getTotalE(): plain sum of member crystal energies (NOT |p|)."""
        return float(np.sum(self.Es, dtype=np.float64))

    def combine(self, other: "Cluster") -> None:
        self.px += other.px
        self.py += other.py
        self.pz += other.pz
        self.ietas += other.ietas
        self.iphis += other.iphis
        self.Es += other.Es

    def copy(self) -> "Cluster":
        """C++ copy constructor: copies vec/Es/ietas/iphis only."""
        return Cluster(self.px, self.py, self.pz,
                       list(self.ietas), list(self.iphis), list(self.Es))

    # -- image ------------------------------------------------------------
    def make_image(self) -> np.ndarray:
        """Literal port of Cluster::makeImage(). Mutates ietas/iphis in place
        and fills xcoords/ycoords (compute_En depends on them)."""
        ietas = self.ietas
        iphis = self.iphis

        # 1) unwrap clusters straddling the iphi = 360/1 seam
        iphi_max = max(iphis)
        iphi_min = min(iphis)
        if iphi_max - iphi_min > ISIZE and iphi_max > 340 and iphi_min < 20:
            for ii in range(len(iphis)):
                if iphis[ii] < ISIZE:
                    iphis[ii] = 360 + iphis[ii] + 1

        # 2) the "+= min" quirk (see module docstring) -- cancels at step 3
        ieta_min = min(ietas)
        iphi_min = min(iphis)
        for ii in range(len(ietas)):
            ietas[ii] += ieta_min + (ISIZE // 2)
        for ii in range(len(iphis)):
            iphis[ii] += iphi_min + (ISIZE // 2)

        # 3) centre on the FIRST crystal of the member list
        cxval = ietas[0]
        cyval = iphis[0]
        self.xcoords = [ietas[ii] - cxval + (ISIZE // 2) for ii in range(len(ietas))]
        self.ycoords = [iphis[ii] - cyval + (ISIZE // 2) for ii in range(len(iphis))]

        # 4) y-only shifting (x gets no equivalent treatment)
        y_min = min(self.ycoords)
        while y_min < 0:
            shift = abs(y_min)
            self.ycoords = [d + shift for d in self.ycoords]
            y_min = min(self.ycoords)

        y_max = max(self.ycoords)
        if y_max > ISIZE:
            shift = y_max - ISIZE
            self.ycoords = [d - shift for d in self.ycoords]

        # 5) fill, clamping both axes into range
        img = np.zeros((ISIZE, ISIZE), dtype=np.float32)
        for ii in range(len(self.xcoords)):
            x = self.xcoords[ii]
            y = self.ycoords[ii]
            if x >= ISIZE:
                x = ISIZE - 1
            if y >= ISIZE:
                y = ISIZE - 1
            if x < 0:
                x = 0
            if y < 0:
                y = 0
            img[x][y] += np.float32(self.Es[ii])

        self.image = img
        return img

    def compute_en(self, n: float) -> float:
        """Cluster::compute_En(n). Requires make_image() to have run first
        (xcoords/ycoords are filled there)."""
        total_e = float(np.sum(self.Es, dtype=np.float64))
        num = 0.0
        for xx, yy in zip(self.xcoords, self.ycoords):
            num += pow(xx * xx + yy * yy, n / 2.0)
        return num / total_e


# --------------------------------------------------------------------------
# Clustering
# --------------------------------------------------------------------------
def _delta_r(c1: Cluster, c2: Cluster) -> float:
    deta = c1.eta() - c2.eta()
    dphi = c1.phi() - c2.phi()
    while dphi > math.pi:
        dphi -= 2.0 * math.pi
    while dphi <= -math.pi:
        dphi += 2.0 * math.pi
    return math.sqrt(deta * deta + dphi * dphi)


def _find_nearest(c1: Cluster, in_c: list) -> tuple:
    """MLPhotonProducer::FindNearest -- scans from index 1, defaults to index 1."""
    min_r = 999.0
    index = 1
    for i in range(1, len(in_c)):
        dr = _delta_r(c1, in_c[i])
        if dr < min_r:
            min_r = dr
            index = i
    return index, min_r


def _do_pairings(in_c: list, radius: float) -> list:
    """MLPhotonProducer::DoPairings. One pass: each cluster is merged with at
    most one partner, hence the outer convergence loop in build_clusters()."""
    out_c = []
    in_c = list(in_c)
    while len(in_c) > 0:
        c1 = in_c[0].copy()
        if len(in_c) == 1:
            out_c.append(c1)
            in_c.pop(0)
        else:
            idx, best_dr = _find_nearest(c1, in_c)
            c2 = in_c[idx].copy()
            if best_dr > radius:
                out_c.append(c1)
                out_c.append(c2)
            else:
                c1.combine(c2)
                out_c.append(c1)
            in_c.pop(idx)
            in_c.pop(0)
    return out_c


def build_clusters(rh_ieta, rh_iphi, rh_eta, rh_phi, rh_energy,
                   match_dr: float = MATCH_DR) -> list:
    """Reproduce the cluster list MLPhotonProducer::produce() builds.

    Inputs are the per-event rh_* vectors from EBRecHitDumper, IN THE ORDER THE
    DUMPER WROTE THEM (which is the order the producer iterates the RecHit
    collection). Only hits with E > 0 are used -- the producer's zero
    suppression is `> 0.`, so a dumper run with minHitEnergy = 0.0 must still be
    filtered here.
    """
    rh_ieta = np.asarray(rh_ieta)
    rh_iphi = np.asarray(rh_iphi)
    rh_eta = np.asarray(rh_eta, dtype=np.float64)
    rh_phi = np.asarray(rh_phi, dtype=np.float64)
    rh_energy = np.asarray(rh_energy, dtype=np.float64)

    keep = rh_energy > 0.0
    clusters = []
    for ieta, iphi, eta, phi, e in zip(rh_ieta[keep], rh_iphi[keep],
                                       rh_eta[keep], rh_phi[keep], rh_energy[keep]):
        pt = e / math.cosh(eta)
        clusters.append(Cluster(
            px=pt * math.cos(phi),
            py=pt * math.sin(phi),
            pz=pt * math.sinh(eta),
            ietas=[int(ieta)], iphis=[int(iphi)], Es=[float(e)],
        ))

    size_before, size_after = 1, 0
    while size_before != size_after:
        # descending |p|^2; stable here vs unstable std::sort online (see docstring)
        order = np.argsort([-c.mag2 for c in clusters], kind="stable")
        clusters = [clusters[i] for i in order]
        size_before = len(clusters)
        clusters = _do_pairings(clusters, match_dr)
        size_after = len(clusters)

    return clusters


# --------------------------------------------------------------------------
# Lorentz vector
# --------------------------------------------------------------------------
def calculate_lorentz_vector(moe: float, energy: float, eta: float, phi: float,
                             zpv: float) -> dict:
    """MLPhotonProducer::calculateLorentzVector -- corrects eta for the primary
    vertex z assuming the shower sits on the ECAL barrel surface (R = 129 cm)."""
    theta = 2.0 * math.atan(math.exp(-1.0 * abs(eta)))
    sign = 1 if eta >= 0 else -1
    z = ECAL_R / math.tan(theta) * sign
    zp = abs(zpv) + abs(z)
    thetaprime = math.atan(ECAL_R / zp)
    etaprime = math.log(math.tan(thetaprime / 2.0))

    if (zpv < 0 and z >= 0) or (zpv >= 0 and z >= 0 and z > zpv) or (zpv < 0 and z < 0 and z >= zpv):
        etaprime *= -1.0

    m = moe * energy
    p = math.sqrt(max(energy * energy - m * m, 0.0))
    pt = p / math.cosh(etaprime)
    return {"pt": pt, "eta": etaprime, "phi": phi, "mass": m}


# --------------------------------------------------------------------------
# ONNX inference
# --------------------------------------------------------------------------
class MLPhotonModels:
    """Runs the two ONNX models exactly as MLPhotonProducer does.

    Note the asymmetry in the inputs, which is easy to miss in the C++:
      regressor  <- image normalised to unit sum, plus cluster eta
      classifier <- RAW (unnormalised) energy image
    """

    def __init__(self, classifier_path: str, regressor_path: str,
                 providers=("CPUExecutionProvider",)):
        import onnxruntime as ort

        so = ort.SessionOptions()
        so.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        self.classifier = ort.InferenceSession(classifier_path, so, providers=list(providers))
        self.regressor = ort.InferenceSession(regressor_path, so, providers=list(providers))
        self.cls_input = self.classifier.get_inputs()[0].name
        self.reg_inputs = [i.name for i in self.regressor.get_inputs()]

    def predict_one(self, cluster: Cluster, pv_z: float) -> dict:
        img = cluster.image if cluster.image is not None else cluster.make_image()
        flat = img.reshape(-1).astype(np.float32)

        img_sum = float(np.sum(flat, dtype=np.float32))
        r_img = (flat / img_sum) if img_sum != 0.0 else flat.copy()

        eta = cluster.eta()
        energy = cluster.total_energy()

        cls_out = self.classifier.run(
            None, {self.cls_input: flat.reshape(1, 1, ISIZE, ISIZE)})[0].reshape(-1)
        reg_out = self.regressor.run(None, {
            self.reg_inputs[0]: r_img.reshape(1, 1, ISIZE, ISIZE),
            self.reg_inputs[1]: np.array([[eta]], dtype=np.float32),
        })[0].reshape(-1)

        # The C++ softmaxes the classifier output even though the graph already
        # ends in LogSoftmax; softmax(log_softmax(x)) == softmax(x), so this
        # matches.
        expo = np.exp(cls_out.astype(np.float64))
        denom = float(np.sum(expo))

        moe = float(reg_out[0])
        p4 = calculate_lorentz_vector(moe, energy, eta, cluster.phi(), pv_z)

        return {
            "pt": p4["pt"], "eta": p4["eta"], "phi": p4["phi"], "mass": p4["mass"],
            "moe": moe,
            "energy": energy,
            "cluster_eta": eta,
            "monophotonScore": float(expo[0] / denom),
            "diphotonScore": float(expo[1] / denom),
            "hadronScore": float(expo[2] / denom),
            "r1": cluster.compute_en(1.0) / cluster.compute_en(0.0),
            "r2": cluster.compute_en(2.0) / cluster.compute_en(0.0),
            "r3": cluster.compute_en(3.0) / cluster.compute_en(0.0),
            "n_crystals": len(cluster.Es),
        }

    def predict(self, clusters: list, pv_z: float) -> list:
        return [self.predict_one(c, pv_z) for c in clusters]
