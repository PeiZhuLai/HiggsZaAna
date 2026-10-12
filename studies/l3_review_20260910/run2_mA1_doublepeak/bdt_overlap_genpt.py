"""Overlap events: is the BDT preference from the true ALP boost or from the overestimated energy?
pho1Pt_oHm of the overlap events is replaced by its value with the generator-level photon pT
(leading gen ALP photon pT / reco m_llgg x reco-to-gen pT ratio removed), keeping everything else."""
import glob, pickle
import numpy as np, uproot
SIG = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M1"
MODEL = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_lowmass_run3.pkl"
FEATS = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm", "param"]
BR = FEATS[:-1] + ["ALP_m", "H_m", "factor", "GenALP_dR_gg", "GenALPLeadPho_pt", "pho1Pt"]
a = {}
for p in sorted(glob.glob(SIG + "/*.root")):
    for k, v in uproot.open(p)["inclusive"].arrays(BR, library="np").items():
        a.setdefault(k, []).append(v)
a = {k: np.concatenate(v) for k, v in a.items()}
a["param"] = (a["ALP_m"] - 1.0) / a["H_m"]
sel = (a["H_m"] > 95) & (a["H_m"] < 180)
a = {k: v[sel] for k, v in a.items()}
clf = pickle.load(open(MODEL, "rb"))
w, ov = a["factor"], a["GenALP_dR_gg"] < 0.05
def rate(X, m):
    s = clf.predict_proba(X[m])[:, 1]
    return w[m][s > 0.963].sum() / w[m].sum()
X = np.column_stack([a[f] for f in FEATS])
def q(x, m): return np.percentile(x[m], 50)
print("median pho1Pt / gen lead-photon pT: overlap %.2f  normal %.2f" % (
    q(a["pho1Pt"] / a["GenALPLeadPho_pt"], ov), q(a["pho1Pt"] / a["GenALPLeadPho_pt"], ~ov)))
print("median gen lead-photon pT / 125: overlap %.3f  normal %.3f" % (
    q(a["GenALPLeadPho_pt"] / 125., ov), q(a["GenALPLeadPho_pt"] / 125., ~ov)))
Xg = X.copy(); Xg[:, 3] = a["GenALPLeadPho_pt"] / a["H_m"]
print("WP pass rate, overlap: reco %.3f  with gen photon pT %.3f" % (rate(X, ov), rate(Xg, ov)))
print("WP pass rate, normal : reco %.3f  with gen photon pT %.3f" % (rate(X, ~ov), rate(Xg, ~ov)))
