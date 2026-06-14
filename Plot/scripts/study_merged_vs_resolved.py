#!/usr/bin/env python3
"""
Study where to use the MERGED-photon category vs the RESOLVED-photon category
(vs both) as a function of m_a.

Physics driver: the two photons from a -> gamma gamma have an opening angle that
scales ~ 2*m_a / pT(a).  When dR(gamma,gamma) drops below the ECAL supercluster
merging scale (~3-4 crystals, dR ~ 0.05-0.07) the two photons are reconstructed
as ONE object -> the resolved (N_gamma >= 2) selection loses them, and the merged
single-photon reconstruction is required.

For every available m_a point we read the signal MC and report, at generator
acceptance (both gen photons in |eta| < 2.5), the distribution of the *true*
opening angle dR(g1,g2), plus the reconstructed-object multiplicity actually
seen by the two taggers.

Inputs (2024 era is representative for the gen opening angle, which depends only
on the a boost, not the data-taking era):
  RESOLVED tagger output : parquet_DNA/Sig_MC/mA_M<X>_2024/merged_nominal.parquet      (m_a = 1..30)
  MERGED  tagger output  : parquet_merged_DNA_tmp/.../mA_MLNANO_M<X>_2024/merged_nominal.parquet (m_a = 0.1..1)
"""
import glob
import numpy as np
import pyarrow.parquet as pq

# ---- ECAL merging scale (barrel crystal ~0.0174 in eta-phi) -----------------
DR_MERGE   = 0.07   # below this the two photons share one supercluster -> merged regime
DR_RESOLVE = 0.30   # above this they are comfortably two separate photons

RES_BASE = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA/Sig_MC"
MRG_BASES = [
    "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all",
    "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_merged_DNA_tmp/Sig_MC_MLNANO",
    "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_merged_DNA_tmp/Sig_MC_MLNANO_M1",
]

# m_a points: (label, gen m_a in GeV, source)
SUBGEV = [(f"M0p{i}", i/10.0) for i in range(1, 10)]
RESOLV = [("M1", 1.0), ("M2", 2.0), ("M3", 3.0), ("M4", 4.0), ("M5", 5.0),
          ("M6", 6.0), ("M7", 7.0), ("M8", 8.0), ("M9", 9.0), ("M10", 10.0),
          ("M15", 15.0), ("M20", 20.0), ("M25", 25.0), ("M30", 30.0)]


def find_resolved(tag):
    p = f"{RES_BASE}/mA_{tag}_2024/merged_nominal.parquet"
    return p if glob.glob(p) else None


def find_merged(tag):
    for b in MRG_BASES:
        p = f"{b}/mA_MLNANO_{tag}_2024/merged_nominal.parquet"
        if glob.glob(p):
            return p
    return None


def gen_dr(path):
    """True dR(g1,g2) for events with both gen photons in |eta|<2.5."""
    t = pq.read_table(path, columns=[
        "GenALPPhoton1_eta", "GenALPPhoton1_phi",
        "GenALPPhoton2_eta", "GenALPPhoton2_phi"]).to_pandas()
    e1, p1 = t["GenALPPhoton1_eta"].to_numpy(), t["GenALPPhoton1_phi"].to_numpy()
    e2, p2 = t["GenALPPhoton2_eta"].to_numpy(), t["GenALPPhoton2_phi"].to_numpy()
    acc = (np.abs(e1) < 2.5) & (np.abs(e2) < 2.5)
    de = e1 - e2
    dp = np.arctan2(np.sin(p1 - p2), np.cos(p1 - p2))
    dr = np.sqrt(de**2 + dp**2)
    return dr[acc & np.isfinite(dr)]


def reco_npho(path):
    """Reconstructed sub-leading ALP photon pt; >0 means a 2nd photon object exists."""
    f = pq.ParquetFile(path)
    if "ALP_sublead_photon_pt" not in f.schema.names:
        return None
    return pq.read_table(path, columns=["ALP_sublead_photon_pt"]).to_pandas()["ALP_sublead_photon_pt"].to_numpy()


HDR = f"{'m_a':>6} {'src':>4} {'N_acc':>7} {'dR med':>7} {'dR p10':>7} {'dR p90':>7} " \
      f"{'<%.2f' % DR_MERGE:>7} {'mid%':>6} {'>%.2f' % DR_RESOLVE:>7}"
print(HDR)
print("-" * len(HDR))


def row(label, ma, path, src):
    dr = gen_dr(path)
    n = len(dr)
    if n == 0:
        print(f"{label:>6} {src:>4}  (no events in acceptance)")
        return
    frac_m = np.mean(dr < DR_MERGE) * 100
    frac_r = np.mean(dr > DR_RESOLVE) * 100
    frac_mid = 100 - frac_m - frac_r
    print(f"{ma:>6.1f} {src:>4} {n:>7d} "
          f"{np.median(dr):>7.3f} {np.percentile(dr,10):>7.3f} {np.percentile(dr,90):>7.3f} "
          f"{frac_m:>6.1f} {frac_mid:>6.1f} {frac_r:>6.1f}")


for tag, ma in SUBGEV:
    p = find_merged(tag)
    if p:
        row(tag, ma, p, "mrg")
for tag, ma in RESOLV:
    p = find_resolved(tag)
    if p:
        row(tag, ma, p, "res")

print()
print(f"Columns: dR percentiles of TRUE dR(g1,g2) at gen acceptance (|eta|<2.5 both photons).")
print(f"  <{DR_MERGE}  : both photons share one supercluster -> MERGED category needed")
print(f"  mid     : transition ({DR_MERGE}-{DR_RESOLVE}) -> both categories have acceptance")
print(f"  >{DR_RESOLVE} : cleanly separated -> RESOLVED category")
