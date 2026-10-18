#!/usr/bin/env python3
"""NPS-25-014 Summary and abstract synced to AN-25-172 Sec 10/11 (blinded: no observed-limit statements;
2-10 fb for m_a = 2-30, 28 fb at m_a = 1; up to ~36% better than Run 2 above 3 GeV at most masses, weaker at
m_a = 1-3 mainly because of the sculpting protection; C_ZH/Lambda up to ~22% better)."""
PAS = "/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/"

def edit(fn, pairs):
    s = open(PAS + fn).read()
    for a, b in pairs:
        assert s.count(a) == 1, (fn, a[:80])
        s = s.replace(a, b)
    open(PAS + fn, "w").write(s); print("edited", fn)

edit("Sec-11-00-Summary.tex", [
    (r"""with $\ell=\Pe,\Pgm$. The expected limits are in the range 2--10\unit{fb} over
$2\leq\mass{a}\leq30~\GeV$ and improve by more than 40\% on those from the
preceding Run~2 search~\cite{CMS:2021wlt}. No significant excess over the
background expectation is observed at any mass point. The cross-section limits
are also interpreted in terms of the ALP effective coupling \cZhLambda, which
is constrained about 30\% more tightly than in Run~2.""",
     r"""with $\ell=\Pe,\Pgm$. The expected limits are in the range 2--10\unit{fb} over
$2\leq\mass{a}\leq30~\GeV$ and about 28\unit{fb} at $\mass{a}=1~\GeV$. Above
$\mass{a}=3~\GeV$ they improve by up to about 36\% on those from the preceding
Run~2 search~\cite{CMS:2021wlt} at most mass points, while at
$\mass{a}=1$--3\GeV\ they are weaker, mainly because of the protection against
background sculpting by the event classifier. The cross-section limits are also
interpreted in terms of the ALP effective coupling \cZhLambda, which is
constrained up to about 22\% more tightly than in Run~2."""),
])

edit("NPS-25-014.tex", [
    (r"""No significant excess over the standard model background expectation is observed. Upper limits at 95\% confidence level are set on the product of the Higgs boson production cross section and the branching fraction into this final state; the expected limits range from 2 to 10~fb for $2 \leq m_{a} \leq 30$ GeV and improve by more than 40\% on those of the preceding Run 2 search.""",
     r"""Expected upper limits at 95\% confidence level are set on the product of the Higgs boson production cross section and the branching fraction into this final state; they range from 2 to 10~fb for $2 \leq m_{a} \leq 30$ GeV, reach about 28~fb at $m_{a} = 1$ GeV, and improve by up to about 36\% on those of the preceding Run 2 search at most masses above 3 GeV."""),
])
