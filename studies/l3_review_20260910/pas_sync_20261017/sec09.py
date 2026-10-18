#!/usr/bin/env python3
"""NPS-25-014 Sec 9 synced to AN-25-172 Sec 9 (tab:uncertainty_summary and Sec 9.1-9.3):
pileup 2-3%; conversion-safe electron veto 2-3% (yield-weighted mean), uncorrelated across eras (per-era 5-12% in
2022-2023, <1% in 2024); shape reweighting = asymmetric log-normal from three control regions, 1-4% (mean 2%);
QCD scale/PDF/alpha_s only in the B and C_ZH/Lambda datacard (weaken those limits by 0-4%, AN Sec 10);
shape (scale/resolution) effects on yield <1%, on peak/width typically 0-0.5 GeV.
The previous closing claim (QCD scale + background modeling degrade the cross-section limit by ~10%) has no
counterpart in the AN and is not compatible with the cross-section datacard, which has no QCD-scale parameter."""
PAS = "/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/"
fn = "Sec-09-00-SystematicUncertainty.tex"
s = open(PAS + fn).read()
pairs = [
    (r"""categories. Their effect on the signal yield is below 1\%, while the signal peak
position and resolution shift by at most 0.5\GeV\ across the \mass{a}\ range.""",
     r"""categories. Their effect on the signal yield is below 1\%, while the signal peak
position and resolution typically shift by up to 0.5\GeV\ across the \mass{a}\ range."""),
    (r"""  parameter. The effect on the signal yield is 1--4\% across the \mass{a}\ range.""",
     r"""  parameter. The effect on the signal yield is 2--3\% across the \mass{a}\ range."""),
    (r"""  separate uncertainty is assigned to the photon conversion-safe electron veto
  and, because the pixel-detector running conditions differ from year to year,
  is implemented with one nuisance parameter per year.""",
     r"""  separate uncertainty is assigned to the photon conversion-safe electron veto
  and, because the pixel-detector running conditions differ between the
  data-taking periods, is treated as uncorrelated across the eras; it changes the
  yield of the corresponding era by 4--12\% in 2022 and 2023 and by less than 1\%
  in 2024, i.e.\ by 2--3\% on the total signal yield."""),
    (r"""  \item \textbf{Training-variable reweighting.} The iterative data-to-simulation
  reweighting applied to the BDT input variables is propagated by varying each
  correction by 100\% of its value, i.e.\ taking the full difference between the
  reweighted and unreweighted templates; the impact on the signal yield is
  1--2\%.""",
     r"""  \item \textbf{Training-variable reweighting.} The data-to-simulation
  reweighting of the BDT input variables, which is also applied to the signal,
  is rederived in two alternative control regions that in addition require the
  dilepton mass to lie away from the Z boson peak. The largest upward and
  downward changes of the signal efficiency after the BDT selection, among the
  two alternative reweightings and no reweighting, define an asymmetric
  log-normal uncertainty for each \mass{a}; its average effect on the signal
  yield is 1--4\%."""),
    (r"""  section used to normalize the signal~\cite{LHCHiggsXSHandbook4,Karlberg2024};
  the effects on the yield are 6\% (QCD scale), 3\% ($\alpha_{s}$), and 2\%
  (PDF).""",
     r"""  section used to normalize the signal~\cite{LHCHiggsXSHandbook4,Karlberg2024};
  the effects on the yield are 6\% (QCD scale), 3\% ($\alpha_{s}$), and 2\%
  (PDF). These uncertainties are included only in the interpretations in terms
  of $\mathcal{B}(\hToZa)$ and \cZhLambda; the limits on
  $\sigma(\pp\to\PH)\,\mathcal{B}$ do not include them."""),
    (r"""    signal normalization, with their type, correlation across the Run~3 eras,
    and $1\sigma$ effect on the signal yield.}""",
     r"""    signal normalization, with their type, correlation across the Run~3 eras,
    and $1\sigma$ effect on the signal yield, averaged over the eras and lepton
    channels and given as a range over the \mass{a}\ hypotheses. The QCD scale,
    $\alpha_{s}$, and PDF uncertainties enter only the interpretations in terms
    of $\mathcal{B}(\hToZa)$ and \cZhLambda.}"""),
    (r"""    Pileup modeling                            & Norm       & Correlated & 1--4 \\""",
     r"""    Pileup modeling                            & Norm       & Correlated & 2--3 \\"""),
    (r"""    Training-variable reweighting              & Norm       & Correlated & 1--2 \\""",
     r"""    Training-variable reweighting              & Norm (asymmetric) & Correlated & 1--4 \\
    Photon conversion-safe electron veto       & Norm       & Uncorrelated & 2--3 \\"""),
    (r"""The inclusion of the systematic uncertainties has a moderate impact on the
expected limits. The dominant contributions are the perturbative truncation
(QCD scale) uncertainty in the signal cross section and the background modeling,
which together degrade the expected cross-section limit by about 10\% across the
probed \mass{a}\ range.""",
     r"""The expected limits are dominated by the statistical uncertainty. Among the
systematic effects, the background modeling enters through the discrete
profiling of the envelope of Section~\ref{sec:bkg}; the theoretical
uncertainties in the H boson production cross section, of which the QCD scale
is the largest, weaken the expected limits on $\mathcal{B}(\hToZa)$ and
\cZhLambda\ by 0--4\%."""),
]
for a, b in pairs:
    assert s.count(a) == 1, a[:80]
    s = s.replace(a, b)
open(PAS + fn, "w").write(s); print("edited", fn)

# second pass: drop the unsupported "dominated by the statistical uncertainty" sentence (the section opening
# already quotes the AN statement that the signal-yield systematics are small compared with the statistical one)
s = open(PAS + fn).read()
a = r"""The expected limits are dominated by the statistical uncertainty. Among the
systematic effects, the background modeling enters"""
b = r"""Among the systematic effects, the background modeling enters"""
if s.count(a) == 1:
    open(PAS + fn, "w").write(s.replace(a, b)); print("edited", fn, "(second pass)")
