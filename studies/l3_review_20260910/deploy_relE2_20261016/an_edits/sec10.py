#!/usr/bin/env python3
"""AN-25-172 Sec 10 for the low-mass BDT with the sub-leading-photon sigma_E/E input (low_swapPt_relE2).
Expected cross-section-card limits (r x 100 fb): mA1 28.28, mA2 9.53, mA3 5.74 fb (Combine/output_combine_results_xsec);
Run 2 17.8 / 4.7 / 2.8 fb -> ratios 1.59 / 2.03 / 2.05. Signal efficiency at the working point (tab:boundary,
optimize_run3UL/nCat1_all_M*.json): 0.110 / 0.339 / 0.448. m_a = 1 GeV composition after the working point
(Plot/plots/mA1_doublepeak/summary.txt): generator dR < 0.05 4%, second m_gg peak 2.8%; in-window fraction of the
signal model 94.5-97.2% (signalFit log). Backup: Sec-10-00-ResultAndInterpretation.tex.bak_relE2_20261016"""
AN = "/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172/"
fn = "Sec-10-00-ResultAndInterpretation.tex"
s = open(AN + fn).read()
pairs = [
    (r"""rising to about 55~fb at $m_{a}=1~\GeV$ where the two photons are most collimated""",
     r"""rising to about 28~fb at $m_{a}=1~\GeV$ where the two photons are most collimated"""),
    (r"""where they are weaker than the Run 2 ones by factors of about 3.1, 2.2 and 2.3 at $m_{a}=1$, 2 and 3~\GeV.""",
     r"""where they are weaker than the Run 2 ones by factors of about 1.6, 2.0 and 2.0 at $m_{a}=1$, 2 and 3~\GeV."""),
    (r"""about 55~fb at $m_{a}=1~\GeV$, with the best expected sensitivity at""",
     r"""about 28~fb at $m_{a}=1~\GeV$, with the best expected sensitivity at"""),
    (r"""expected limits are instead weaker: 54.5, 10.4 and 6.5~fb at $m_{a}=1$, 2 and 3~\GeV,
compared with 17.8, 4.7 and 2.8~fb in Run 2, i.e.\ larger by factors of about 3.1, 2.2
and 2.3.""",
     r"""expected limits are instead weaker: 28.3, 9.5 and 5.7~fb at $m_{a}=1$, 2 and 3~\GeV,
compared with 17.8, 4.7 and 2.8~fb in Run 2, i.e.\ larger by factors of about 1.6, 2.0
and 2.0."""),
    (r"""spectrum. The low-mass BDT therefore uses only five inputs, without shower-shape variables
(Table~\ref{tab:bdt_input_lowmass}), and its working points at $m_{a}=1$, 2 and 3~\GeV are
chosen with the sculpting ratio and the background-only closure test rather than at the maximum
expected significance (Sections~\ref{sec:sculpt_R} and~\ref{sec:sculpt_closure}). They keep
9.5\%, 40\% and 59\% of the preselected signal, compared with 56--77\% at $m_{a}\geq4~\GeV$
(Table~\ref{tab:boundary}).""",
     r"""spectrum. The low-mass BDT therefore uses only five inputs, without shower-shape variables
and without the photon transverse momenta (Table~\ref{tab:bdt_input_lowmass}), and its working
points at $m_{a}=1$, 2 and 3~\GeV are chosen with the background-only closure test rather than
at the maximum expected significance (Sections~\ref{sec:sculpt_R} and~\ref{sec:sculpt_closure}).
They keep 11\%, 34\% and 45\% of the preselected signal, compared with 56--77\% at
$m_{a}\geq4~\GeV$ (Table~\ref{tab:boundary})."""),
]
for a, b in pairs:
    assert s.count(a) == 1, a[:80]
    s = s.replace(a, b)
# the m_a = 1 GeV overlapping-shower paragraph: from "At $m_{a}=1~\GeV$ this BDT also favors" to "...rather than from
# the detector."
i = s.index(r"""(Table~\ref{tab:boundary}). At $m_{a}=1~\GeV$ this BDT also favors the events with overlapping""")
j = s.index(r"""the detector.""", i) + len(r"""the detector.""")
old = s[i:j]
assert "26--47\\%" in old and "two thirds" in old, old[:200]
new = r"""(Table~\ref{tab:boundary}). At $m_{a}=1~\GeV$ the previous low-mass BDT, which used the
leading-photon $p_{T}/m_{\ell\ell\gamma\gamma}$ instead of the relative energy uncertainty of the
sub-leading photon, favored the events with overlapping photon showers
(Figure~\ref{fig:mA1_doublepeak}), whose overestimated photon energy raises that input. These
events are about a quarter of the signal after the preselection; they formed two thirds of the
signal selected by that BDT and were excluded from the window used to build the signal model,
because their energy response is not validated with data (Section~\ref{sec:signal_ngauss}), so that
only 26--47\% of the selected signal entered the model. With the current low-mass BDT they are 4\%
of the selected signal, and 94.5--97.2\% of it, depending on era and lepton channel, falls inside
the window of the signal model. The overlapping showers themselves are not specific to Run~3: in
the 2018 simulation of the Run~2 analysis they occur at the same rate (8.0\% of the events have a
generator-level $\Delta R(\gamma,\gamma)<0.05$, compared with 8.5\% in the 2024 simulation) and
with the same energy overestimate, so the difference arose from the selection rather than from
the detector."""
s = s[:i] + new + s[j:]
open(AN + fn, "w").write(s)
print("Sec 10 edited")
