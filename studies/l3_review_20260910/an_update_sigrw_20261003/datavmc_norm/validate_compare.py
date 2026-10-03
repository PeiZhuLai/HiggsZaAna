#!/usr/bin/env python3
"""Validation of the signal normalization N (and the --optimizeBranches fix) in 1_prepare_dataVmc.py,
mA5, one era, SR region. Files (this directory):
  ..._val_new_<era>.root        current script (--optimizeBranches): w * R * N
  ..._val_oldnoopt_<era>.root   pre-change script without --optimizeBranches: w * R (correct R)
  ..._val_old.root              pre-change script with --optimizeBranches (= the 09-29 condor setup)
Independent emulation with uproot + SidebandReweighter.weights_for_dataframe (true-mass param),
channel = n_electrons==2 / n_muons==2 (apply_bdt_sig.py), N from signal_rw_norms.txt."""
import sys, os, numpy as np, uproot
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from sideband_reweight import SidebandReweighter
W = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/an_update_sigrw_20261003"
D = W + "/datavmc_norm"
ERA = sys.argv[1] if len(sys.argv) > 1 else "2022preEE"
MASS = 5
norms = {}
for line in open(f"{W}/signal_rw_norms.txt"):
    p = line.split(); norms[(p[0], p[1], p[2])] = float(p[4])
Ne, Nm = norms[(f"mA_M{MASS}", ERA, "ele")], norms[(f"mA_M{MASS}", ERA, "mu")]
op = lambda f: uproot.open(f"{D}/{f}")["raw_plots"]
new, onop = op(f"ALP_plot_run3_UL_SR_sigrwnorm_val_new_{ERA}.root"), op(f"ALP_plot_run3_UL_SR_sigrwnorm_val_oldnoopt_{ERA}.root")
old = op("ALP_plot_run3_UL_SR_sigrwnorm_val_old.root") if ERA == "2022preEE" else None

rw = SidebandReweighter.from_json("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json")
fr = uproot.open(f"/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M{MASS}/{ERA}.root")["test"].arrays(library="pd")
fr["param"] = (fr["ALP_m"] - MASS) / fr["H_m"]
R = np.asarray(rw.weights_for_dataframe(fr), dtype=float)
fz = fr.copy()
for c in ("pho1IetaIeta55", "pho2IetaIeta55", "pho1R9", "pho2R9", "Z_m"): fz[c] = 0.0
Rz = np.asarray(rw.weights_for_dataframe(fz), dtype=float)   # what the --optimizeBranches loop computed before the fix
w = fr["weight"].to_numpy(float); H = fr["H_m"].to_numpy(float)
ele = fr["n_electrons"].to_numpy() == 2; mu = fr["n_muons"].to_numpy() == 2
N = Ne * ele + Nm * mu
sel = (H > 115) & (H < 135)
score = fr[f"MVA_Score_mA_M{MASS}"].to_numpy(float)
edges = np.linspace(-0.1, 1.1, 241)
def h(mask, wt):
    return np.histogram(score[mask], bins=edges, weights=wt[mask])[0]
def maxrel(a, b):
    m = np.abs(b) > 1e-9
    return float(np.max(np.abs(a - b)[m] / np.abs(b)[m])) if m.any() else 0.0
key = f"mvaVal_M{MASS}_M{MASS}"
hn, hon = new[key].values(), onop[key].values()
print(f"mA{MASS} {ERA} SR: N_ele={Ne} N_mu={Nm}")
print(f"[1] new      vs emulated w*R*N : max rel bin diff {maxrel(hn, h(sel, w*R*N)):.1e}  totals {hn.sum():.6f} / {h(sel, w*R*N).sum():.6f}")
print(f"[1] oldnoopt vs emulated w*R   : max rel bin diff {maxrel(hon, h(sel, w*R)):.1e}  totals {hon.sum():.6f} / {h(sel, w*R).sum():.6f}")
if old is not None:
    ho = old[key].values()
    print(f"[1] old(09-29 setup) vs emulated w*R(5 vars=0): max rel bin diff {maxrel(ho, h(sel, w*Rz)):.1e}  totals {ho.sum():.6f} / {h(sel, w*Rz).sum():.6f}")
Se, Sm = (w*R)[sel & ele & ~mu].sum(), (w*R)[sel & mu & ~ele].sum()
Sb, S0 = (w*R)[sel & ele & mu].sum(), (w*R)[sel & ~ele & ~mu].sum()
expected = (Ne*Se + Nm*Sm + (Ne+Nm)*Sb) / (Se + Sm + Sb + S0)
print(f"[2] total yield ratio new/oldnoopt = {hn.sum()/hon.sum():.7f}; expected N mix = {expected:.7f}  (diff {hn.sum()/hon.sum()-expected:+.1e}; "
      f"neither-channel share {S0/(Se+Sm+Sb+S0):.2e}, both {Sb:.3g})")
cum = lambda x: x[::-1].cumsum()[::-1] / x.sum()
thr = edges[:-1]; picks = [np.searchsorted(thr, c - 1e-9) for c in (0.5, 0.9, 0.95, 0.975, 0.99, 0.995)]
en, eon = cum(hn), cum(hon)
print("[3] efficiency above threshold, new vs oldnoopt (both correct R):")
for i in picks: print(f"      score > {thr[i]:.3f}: {eon[i]:.7f} {en[i]:.7f}  {en[i]-eon[i]:+.1e}")
print(f"    inclusive: max |delta eff| over 240 thresholds {np.max(np.abs(en-eon)):.1e}")
for name, m in (("ele", ele & ~mu), ("mu", mu & ~ele)):
    a, b = cum(h(sel & m, w*R)), cum(h(sel & m, w*R*N))
    print(f"    per channel {name}: max |delta eff| {np.nanmax(np.abs(a-b)):.1e}")
if old is not None:
    eo = cum(ho)
    print("[5] effect of the --optimizeBranches bug on the 09-29 signal histograms (old vs oldnoopt, same N=1):")
    print(f"    total yield old/oldnoopt = {ho.sum()/hon.sum():.4f}")
    for i in picks: print(f"      score > {thr[i]:.3f}: eff old {eo[i]:.5f}  correct {eon[i]:.5f}  diff {eo[i]-eon[i]:+.4f}")
