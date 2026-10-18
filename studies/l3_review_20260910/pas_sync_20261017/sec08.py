#!/usr/bin/env python3
"""NPS-25-014 Sec 8 synced to AN-25-172 Sec 8.1 (m_a = 1 second population: overlapping showers, energy
overestimated by 40-60%, correctly paired; 0.3% of the selected signal with the low-mass BDT; window 94.5-97.2%).
sigma_eff ranges from the per-era effSigma JSONs of the current signal fits (flashggFinalFit/Signal/outdir_<ch>/
signalFit/Plots/effSigma_<m>_<era>_<ch>.json, 14 masses x 5 eras): muon 1.35-1.76 GeV, electron 1.65-2.87 GeV;
the effSigmaVmA plots were redrawn (they still held the pre-relE2 m_a = 1-3 values) and show one curve per era."""
PAS = "/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/"
fn = "Sec-08-00-SignalExtraction.tex"
s = open(PAS + fn).read()
pairs = [
    (r"""At $\mass{a}=1~\GeV$ the two photons are nearly collinear and the pairing of the
reconstructed photon candidates becomes ambiguous, which produces a second
population of events at $134<\mllgg<137~\GeV$, separated from the signal peak by
a clear minimum. Its kinematics are incompatible with the decay of a 1\GeV\
pseudoscalar, the median leading-photon \pt\ being 33.7\GeV\ against 19.1\GeV\
in the peak, so these mis-paired events are excluded by setting the upper edge
of the fit window to 133\GeV\ at this mass point. The efficiency times
acceptance used to normalize the model is evaluated in the same restricted
window, so that the model is normalized only to the events it describes.""",
     r"""At $\mass{a}=1~\GeV$ the two photons are nearly collinear. For the most
collimated pairs the two electromagnetic showers overlap and each photon cluster
collects part of the energy of its partner, so the reconstructed photon-pair
\pt\ exceeds the generated one by about 40--60\%. This shifts \mllgg\ upward and
produces a second population of events at $134<\mllgg<137~\GeV$, separated from
the signal peak by a clear minimum; the photons in these events are correctly
paired with the generated ones. Since the energy response of overlapping photon
showers has not been validated with data, these events are excluded by setting
the upper edge of the fit window to 133\GeV\ at this mass point. After the
low-mass BDT selection they amount to 0.3\% of the selected signal, and between
94.5 and 97.2\% of the simulated yield, depending on era and channel, falls
inside the window. The efficiency times acceptance used to normalize the model
is evaluated in the same restricted window, so that the model is normalized only
to the events it describes."""),
    (r"""channels. The muon channel has a better resolution than the electron channel
across the entire \mass{a}\ range, $\sigma_{\mathrm{eff}}$ ranging from about
1.5 to 2.2\GeV\ against 2.4 to 3.3\GeV; both channels show the broadening at
very low \mass{a}\ associated with the merged-photon regime.""",
     r"""channels. The muon channel has a better resolution than the electron channel
across the entire \mass{a}\ range, $\sigma_{\mathrm{eff}}$ ranging from about
1.4 to 1.8\GeV\ against 1.7 to 2.9\GeV\ depending on \mass{a}\ and era, with
little dependence on \mass{a}."""),
    (r"""    \caption{Signal $\sigma_{\mathrm{eff}}$ as a function of \mass{a}\ for the
      electron (left) and muon (right) channels, combined across the Run~3 eras.}""",
     r"""    \caption{Signal $\sigma_{\mathrm{eff}}$ as a function of \mass{a}\ for the
      electron (left) and muon (right) channels, shown separately for each Run~3
      era.}"""),
]
for a, b in pairs:
    assert s.count(a) == 1, a[:80]
    s = s.replace(a, b)
open(PAS + fn, "w").write(s); print("edited", fn)
