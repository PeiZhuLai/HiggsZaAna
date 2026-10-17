#!/usr/bin/env python3
"""AN-25-172 Sec 8.1 (m_a = 1 GeV signal model) for the low-mass BDT with the sub-leading-photon sigma_E/E input.
Numbers: Plot/plots/mA1_doublepeak/summary.txt (BDT > 0.980): 134-137 GeV population 0.16/46.41 = 0.3% of the weight,
44% gen dR < 0.05, median response 1.21; peak 120-130: 2%, 1.00; matched > 95%; second m_gg peak 1.31/46.41 = 2.8%.
signalFit log (an_numbers/mA1_sigfit.txt): in-window 94.5-97.2%, chi2/ndf 0.21-0.93, nGauss ele 2,2,2,1,3 mu 2,2,2,2,2.
effSigma (an_numbers/effsigma.txt): 1.43-2.74 GeV. Tail below 123 GeV (an_numbers/tail123.txt): model vs simulation within
6.2 points, largest ele 2023preBPix 18.0% vs 24.2%. Backup: Sec-08-01-Parameterizedsignalmodel.tex.bak_relE2_20261016"""
AN = "/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172/"
fn = "Sec-08-01-Parameterizedsignalmodel.tex"
s = open(AN + fn).read()
pairs = [
    (r"""channel combinations. After the BDT selection, 94\% of the events at
$134<m_{\ell\ell\gamma\gamma}<137~\GeV$ have $\Delta R(\gamma,\gamma)<0.05$ and their
median photon-pair $p_{T}$ response is 1.43, compared with 26\% and 1.01 for the events in
the peak ($120<m_{\ell\ell\gamma\gamma}<130~\GeV$). In both populations the two""",
     r"""channel combinations before the BDT selection. The current low-mass BDT, which does not use the
photon transverse momenta (Table~\ref{tab:bdt_input_lowmass}), selects few of these events: after
its working point the events at $134<m_{\ell\ell\gamma\gamma}<137~\GeV$ are 0.3\% of the selected
signal, 44\% of them have $\Delta R(\gamma,\gamma)<0.05$ and their median photon-pair $p_{T}$
response is 1.21, compared with 2\% and 1.00 for the events in the peak
($120<m_{\ell\ell\gamma\gamma}<130~\GeV$). With the previous low-mass BDT, which used the
leading-photon $p_{T}/m_{\ell\ell\gamma\gamma}$, these events were about two thirds of the selected
signal. In both populations the two"""),
    (r"""photons with overestimated energy rather than mis-pairings. The same events form the
second peak of the $m_{\gamma\gamma}$ distribution at 1.2--1.9~\GeV
(Figures~\ref{fig:ma_true_eras} and~\ref{fig:mA1_doublepeak}). Because""",
     r"""photons with overestimated energy rather than mis-pairings. The same events form the
second peak of the $m_{\gamma\gamma}$ distribution at 1.2--1.9~\GeV
(Figures~\ref{fig:ma_true_eras} and~\ref{fig:mA1_doublepeak}), which holds 2.8\% of the selected
signal. Because"""),
    (r"""the events it describes. Between 26.3\% and 46.8\% of the simulated yield,
depending on era and channel, falls inside the window; the remainder is
discarded rather than concentrated into the narrow core, which would give an
optimistic limit. With this window all ten fits converge, with $\chi^{2}$ per
degree of freedom between $0.25$ and $0.84$. The effective resolution ranges
from $1.75$ to $3.86~\GeV$ depending on era and channel, the spread reflecting
the residual tail discussed below rather than a difference in detector
performance.""",
     r"""the events it describes. Between 94.5\% and 97.2\% of the simulated yield,
depending on era and channel, falls inside the window (between 26.3\% and 46.8\% with the
previous low-mass BDT); the remainder is
discarded rather than concentrated into the narrow core, which would give an
optimistic limit. With this window all ten fits converge, with $\chi^{2}$ per
degree of freedom between $0.21$ and $0.93$. The effective resolution ranges
from $1.43$ to $2.74~\GeV$ depending on era and channel, the spread reflecting
the residual tail discussed below rather than a difference in detector
performance."""),
]
for a, b in pairs:
    assert s.count(a) == 1, a[:90]
    s = s.replace(a, b)
i = s.index(r"""Even inside this window the $m_{a}=1~\GeV$ shape is not a pure resolution peak,""")
j = s.index(r"""effect on the combined limit at this mass point is negligible.""", i) + len(r"""effect on the combined limit at this mass point is negligible.""")
new = r"""Even inside this window the $m_{a}=1~\GeV$ shape is not a pure resolution peak: in nine of
the ten era and channel combinations at least one additional Gaussian component (two in the
electron channel of 2024) is retained alongside the narrow core to absorb the residual tail, and a
single Gaussian is sufficient only in the electron channel of 2023postBPix. In all ten
combinations the model reproduces the fraction of the in-window yield below $123~\GeV$ to within
6.2 percentage points; the largest difference is in the electron channel of 2023preBPix, where the
model places $18.0\%$ of its probability below $123~\GeV$ and the simulated sample $24.2\%$."""
s = s[:i] + new + s[j:]
open(AN + fn, "w").write(s)
print("Sec 8.1 edited")
