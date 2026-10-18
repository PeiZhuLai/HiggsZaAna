#!/usr/bin/env python3
"""NPS-25-014 Sec 10 synced to AN-25-172 Sec 10 (cross-section card, current master):
2-10 fb for m_a = 2-30, best 1.8 fb at m_a = 5, 28 fb at m_a = 1; Run 3 vs Run 2: up to ~36% better, weaker at
m_a = 1-3 (28.3/9.5/5.7 vs 17.8/4.7/2.8 fb, x1.6/2.0/2.0) and at 6, 26, 28 (26/28/17%) and 20, 27, 30 (<=6%);
low-mass reason = sculpting protection (11/34/45% of preselected signal vs 56-77%); S/sqrt(B) bracket 37-55%;
lumi +25%; Wilson improvement up to ~22%; theory uncertainties only in the B and C_ZH/Lambda interpretations.
AN is blinded (observed-limit sentence commented out in AN Sec 10), so the PAS observed-limit statements are removed.
Limit figures moved to plot_limits/2_runLimitsPlot/ as in the AN."""
PAS = "/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/"
fn = "Sec-10-00-ResultAndInterpretation.tex"
s = open(PAS + fn).read()
i = s.index(r"""Figure~\ref{fig:limit} (left) shows the expected upper limits on""")
j = s.index(r"""observed limits are consistent with the expected limits within uncertainties.""") + \
    len(r"""observed limits are consistent with the expected limits within uncertainties.""")
new = r"""Figure~\ref{fig:limit} (left) shows the expected upper limits on
$\sigma(\pp\to\PH)\,\mathcal{B}(\PH\to\Za\to\ell\ell\gamma\gamma)$ as a
function of \mass{a}. The expected limits are in the range 2--10\unit{fb} over
$2\leq\mass{a}\leq30~\GeV$, with the best expected sensitivity of 1.8\unit{fb}
at $\mass{a}=5~\GeV$, and about 28\unit{fb} at $\mass{a}=1~\GeV$, where the two
photons are most collimated. The kinematics of the final-state particles change
with \mass{a}, so the accessible background phase space changes as well, leading
to a mass-dependent sensitivity. Figure~\ref{fig:limit} (right) compares the
Run~3 expected limit with the Run~2 result of Ref.~\cite{CMS:2021wlt}, under the
same reference signal cross section of 0.1\unit{pb}. Above $\mass{a}=3~\GeV$ the
Run~3 expected limits are stronger at most mass points, by up to about 36\%;
they are weaker by 26, 28, and 17\% at $\mass{a}=6$, 26, and 28\GeV\ and by at
most 6\% at $\mass{a}=20$, 27, and 30\GeV. At the lowest masses they are
weaker: 28.3, 9.5, and 5.7\unit{fb} at $\mass{a}=1$, 2, and 3\GeV, compared with
17.8, 4.7, and 2.8\unit{fb} in Run~2. These limits are weakened mainly by the
protection against background sculpting, which the Run~2 analysis did not
quantify: the low-mass classifier uses neither shower-shape variables nor photon
transverse momenta, and its thresholds are chosen with the background-only
closure test rather than at the maximum expected significance, so that they keep
11, 34, and 45\% of the preselected signal at $\mass{a}=1$, 2, and 3\GeV,
compared with 56--77\% at $\mass{a}\geq4~\GeV$. The gain at higher masses comes
mainly from lowering the subleading-lepton \pt\ threshold, which raises the
preselection efficiency from about 2\% in Run~2 to about 4\% in Run~3, together
with the larger integrated luminosity (172 versus 138\fbinv, an increase of about
25\%), the updated event classifier, and the use of dedicated $\PZ+\gamma$
simulated samples in the training. Scaling the expected sensitivity of the first
two effects with $S/\sqrt{B}$ would predict an improvement of between about 37
and 55\%; the smaller achieved value is attributed to the roughly 10\% lower
multivariate signal efficiency, which follows from the degraded \mllgg\
resolution in Run~3."""
s = s[:i] + new + s[j:]
pairs = [
    (r"""directly, its QCD-scale, PDF, and $\alpha_{s}$ uncertainties propagate to the
limit on \cZhLambda\ as well as to the limit on $\sigma\,\mathcal{B}$.
Figure~\ref{fig:limit_cZH} shows the expected and observed limits on
\cZhLambda\ and the comparison with the Run~2 result. The relative improvement
on \cZhLambda\ with respect to Ref.~\cite{CMS:2021wlt} is about 30\%, smaller""",
     r"""directly, its QCD-scale, PDF, and $\alpha_{s}$ uncertainties are included in
the limits on $\mathcal{B}(\hToZa)$ and \cZhLambda, whereas the limit on
$\sigma\,\mathcal{B}$ does not include them; they weaken the expected limits on
\cZhLambda\ by 0--4\%. Figure~\ref{fig:limit_cZH} shows the expected limits on
\cZhLambda\ and the comparison with the Run~2 result. The relative improvement
on \cZhLambda\ with respect to Ref.~\cite{CMS:2021wlt} is up to about 22\%, smaller"""),
    (r"""    \caption{Expected and observed 95\% CL upper limits on \cZhLambda\ as a""",
     r"""    \caption{Expected 95\% CL upper limits on \cZhLambda\ as a"""),
    ("figure_ALP/run3/flashggFinalFit/Plots/plot_limits/Limits_XS}", "figure_ALP/run3/flashggFinalFit/Plots/plot_limits/2_runLimitsPlot/Limits_XS}"),
    ("figure_ALP/run3/flashggFinalFit/Plots/plot_limits/compare_xs_run2.pdf}", "figure_ALP/run3/flashggFinalFit/Plots/plot_limits/2_runLimitsPlot/compare_xs_run2.pdf}"),
    ("figure_ALP/run3/flashggFinalFit/Plots/plot_limits/Limits_Wilson}", "figure_ALP/run3/flashggFinalFit/Plots/plot_limits/2_runLimitsPlot/Limits_Wilson}"),
    ("figure_ALP/run3/flashggFinalFit/Plots/plot_limits/compare_wilson_run2.pdf}", "figure_ALP/run3/flashggFinalFit/Plots/plot_limits/2_runLimitsPlot/compare_wilson_run2.pdf}"),
]
for a, b in pairs:
    assert s.count(a) == 1, a[:80]
    s = s.replace(a, b)
open(PAS + fn, "w").write(s); print("edited", fn)
