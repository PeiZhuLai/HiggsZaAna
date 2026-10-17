#!/usr/bin/env python3
"""AN-25-172 Sec 7.4 / 7.5 for the low-mass BDT with the sub-leading-photon sigma_E/E input (low_swapPt_relE2).
Numbers: production bias (Bias_nominal/bias_outputs_highstat, 10 x 1000 toys, xsec card) and envelope fits
(fit_results_run3/{1,2,3}); alternatives for m_a = 1 from varsearch_validate/chainlow_relE2/redo_mA1_{noLau,hiord}.txt.
Backups: Sec-07-0{4,5}-*.tex.bak_relE2_20261016"""
AN = "/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172/"


def sub(fn, pairs):
    p = AN + fn
    s = open(p).read()
    for a, b in pairs:
        assert s.count(a) == 1, (fn, a[:80])
        s = s.replace(a, b)
    open(p, "w").write(s)


sub("Sec-07-04-Envelopselection.tex", [
    (r"""$m_{a}=3~\GeV$, it is $-0.089$. At $m_{a}=2$, 4 and 27~\GeV the""",
     r"""$m_{a}=3~\GeV$, it was $-0.089$; with the current low-mass BDT, which uses the relative energy
uncertainty of the sub-leading photon (Table~\ref{tab:bdt_input_lowmass}), and an uncertainty of
$0.965/1.000$, it is $-0.081$ (Bern1, goodness-of-fit probability 0.45). At $m_{a}=2$, 4 and 27~\GeV the"""),
    (r"""lowest order, and removing the power-law member was found not to reduce the pulls
of the remaining functions.
""",
     r"""lowest order, and removing the power-law member was found not to reduce the pulls
of the remaining functions.

At $m_{a}=1~\GeV$, with the current low-mass BDT, the Laurent member Lau1 has a mean pull of
$-0.261\pm0.010$ and is the only function of the whole scan above the threshold; the other two
members give $-0.173$ (Exp1) and $+0.004$ (Pow1), with goodness-of-fit probabilities of 0.43
(Lau1), 0.49 (Exp1) and 0.35 (Pow1). All three are already at the lowest order, so the rule
used at the other mass points, falling back to the next lower order, cannot be applied. Two
alternatives were tested with the same 10000 toys per function. Removing the Laurent family
leaves Exp1 and Pow1 with goodness-of-fit probabilities of 0.49 and 0.35 and moves the excess
to Exp1 ($-0.267$, with $-0.035$ for Pow1), while the expected limit becomes 3\% weaker.
Allowing higher orders of the exponential and power-law families selects Lau1, Exp1 and Pow3,
and Pow3 gives $-0.317$ (Lau1 $-0.259$). Neither alternative removes the excess, so the
envelope of Table~\ref{tab:bkg_true} is kept at $m_{a}=1~\GeV$, and the residual bias is
reported in Section~\ref{sec:bias_studies}.
"""),
])

sub("Sec-07-05-Biasstudy.tex", [
    (r"""described in Section~\ref{sec:envelope_selection}, the mean of the pull
distribution stays within 0.2 for every one of the 108 truth functions of the
scan (one truth function per family and mass point, i.e.\ four families over the
30 mass points minus the ten points at which the Bernstein family and the two points
at which the Laurent family is absent). The largest values are $-0.182\pm0.011$ (Exp3) at
$m_{a}=10~\GeV$, $-0.180\pm0.010$ (Lau1) at $m_{a}=30~\GeV$ and $-0.178\pm0.011$ (Exp3) at
$m_{a}=9~\GeV$; the mean pull of every
individual truth function is listed in
Tables~\ref{tab:bkg_diagnostics_mA1to10}--\ref{tab:bkg_diagnostics_mA21to30}.
The largest of these corresponds to an increase of the total uncertainty of only
$\sqrt{1+0.182^{2}}\approx1.016$, i.e.\ a 2\% effect on the signal strength
measurement. These results indicate that the impact of the residual bias on the
analysis is expected to be negligible.""",
     r"""described in Section~\ref{sec:envelope_selection}, the mean of the pull
distribution stays within 0.2 for 107 of the 108 truth functions of the
scan (one truth function per family and mass point, i.e.\ four families over the
30 mass points minus the ten points at which the Bernstein family and the two points
at which the Laurent family is absent). The exception is Lau1 at $m_{a}=1~\GeV$, with
$-0.261\pm0.010$; at this mass point Exp1 gives $-0.173$ and Pow1 $+0.004$. Among the other
mass points the largest values are $-0.182\pm0.011$ (Exp3) at
$m_{a}=10~\GeV$, $-0.180\pm0.010$ (Lau1) at $m_{a}=30~\GeV$ and $-0.178\pm0.011$ (Exp3) at
$m_{a}=9~\GeV$; the mean pull of every
individual truth function is listed in
Tables~\ref{tab:bkg_diagnostics_mA1to10}--\ref{tab:bkg_diagnostics_mA21to30}.
Outside $m_{a}=1~\GeV$ the largest of these corresponds to an increase of the total
uncertainty of only $\sqrt{1+0.182^{2}}\approx1.016$, i.e.\ a 2\% effect on the signal strength
measurement, and the impact of the residual bias on the analysis is expected to be negligible.
At $m_{a}=1~\GeV$ the pull of Lau1 is negative, i.e.\ when the background follows this function
the fit underestimates the injected signal by 0.26 of its uncertainty. For an upper limit this
is not a conservative direction: the expected limit at this mass point can be too strong by up
to a quarter of the uncertainty on the signal strength. The envelope cannot be changed to remove
the excess (Section~\ref{sec:envelope_selection}): removing the Laurent family moves it to Exp1
($-0.267$), and higher orders add Pow3 ($-0.317$)."""),
    (r"""brings the largest bias of that mass point to $-0.106$. With the final envelopes and the
final signal model, in which the simulated signal is also reweighted
(Section~\ref{ssec:mcweights}), the largest bias at $m_{a}=3$, 11 and $20~\GeV$ is
$-0.089$, $-0.154$ and $-0.107$.""",
     r"""brings the largest bias of that mass point to $-0.106$. With the final envelopes and the
final signal model, in which the simulated signal is also reweighted
(Section~\ref{ssec:mcweights}), the largest bias at $m_{a}=3$, 11 and $20~\GeV$ is
$-0.081$, $-0.154$ and $-0.107$, the value at $m_{a}=3~\GeV$ being obtained with the current
low-mass BDT."""),
    (r"""the truth function and the vertical lines mark the $\pm0.2$ acceptance threshold. Every
  function of the scan lies inside the threshold, the largest deviation being $-0.182$ at
  $m_{a}=10~\GeV$, which inflates the total uncertainty on the signal strength by only
  $\sqrt{1+0.182^{2}}\approx1.016$.""",
     r"""the truth function and the vertical lines mark the $\pm0.2$ acceptance threshold. Every
  function of the scan lies inside the threshold except Lau1 at $m_{a}=1~\GeV$ ($-0.261$,
  Section~\ref{sec:bias_studies}); elsewhere the largest deviation is $-0.182$ at
  $m_{a}=10~\GeV$, which inflates the total uncertainty on the signal strength by only
  $\sqrt{1+0.182^{2}}\approx1.016$."""),
])
print("Sec 7.4 / 7.5 edited")
