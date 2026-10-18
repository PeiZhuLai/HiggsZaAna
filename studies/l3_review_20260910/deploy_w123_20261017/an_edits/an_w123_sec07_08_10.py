#!/usr/bin/env python3
"""AN-25-172 Sec 7.4, 7.5, 8.1, 10 for the low-mass BDT low_relE2_w123 (production since 2026-10-17).
Numbers (deploy_w123_20261017/an_numbers/):
  bias (Bias_nominal/bias_outputs_highstat, 10 x 1000 toys, xsec card) + envelope GOF (fit_results_run3/{1,2,3}):
    mA1 Lau1 -0.246 +- 0.010 (GOF 0.33), Exp1 -0.159 (0.45), Pow1 +0.001 (0.24); mA3 Bern1 -0.100 (0.32)
    the m_a = 1 alternatives (no Laurent, higher orders) were measured on the previous low-mass model (Lau1 -0.261,
    varsearch_validate/chainlow_relE2/redo_mA1_{noLau,hiord}.txt) and are quoted as such
  Plot/plots/mA1_doublepeak/summary.txt (BDT > 0.981): 134-137 GeV 0.09 / 45.98 = 0.2 %, gen dR < 0.05 in 10 %, response
    1.02; peak 1 % and 1.00; one-to-one matched >= 97 %; second m_gg peak 0.5 %; overlapping showers 1 % of the selection
  signalFit log (w123_numbers.txt): in-window 95.9-98.8 %, chi2/ndf 0.24-0.69, one Gaussian only in ele 2023postBPix,
    two elsewhere; effSigma m_a = 1: 1.33-2.70 GeV
  tail123_w123.txt: model vs simulation below 123 GeV within 5.0 points, largest ele 2023preBPix 19.8 % vs 24.8 %
  limits_w123.txt: xsec card 30.47 / 9.69 / 5.47 fb at m_a 1-3 (Run 2 17.82 / 4.74 / 2.78 fb, ratios 1.71 / 2.05 / 1.97);
    m_a >= 4 unchanged; Plot/output/latexTable_MVAcut_eff_bkg.txt: efficiency 0.117 / 0.362 / 0.463
Backups: <file>.bak_w123_20261017"""
import shutil

AN = "/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172/"
TAG = ".bak_w123_20261017"


def edit(fn, pairs):
    p = AN + fn
    s = open(p).read()
    for a, b in pairs:
        assert s.count(a) == 1, (fn, s.count(a), a[:100])
        s = s.replace(a, b)
    shutil.copy2(p, p + TAG)
    open(p, "w").write(s)
    print("edited", fn, len(pairs))


edit("Sec-07-04-Envelopselection.tex", [
    (r"""$m_{a}=3~\GeV$, it was $-0.089$; with the current low-mass BDT, which uses the relative energy
uncertainty of the sub-leading photon (Table~\ref{tab:bdt_input_lowmass}), and an uncertainty of
$0.965/1.000$, it is $-0.081$ (Bern1, goodness-of-fit probability 0.45).""",
     r"""$m_{a}=3~\GeV$, it was $-0.089$; with the current low-mass BDT, which uses the relative energy
uncertainty of the sub-leading photon (Table~\ref{tab:bdt_input_lowmass}) and trains the signal of
all three low masses inside the window, and an uncertainty of
$0.965/1.000$, it is $-0.100$ (Bern1, goodness-of-fit probability 0.32)."""),
    (r"""At $m_{a}=1~\GeV$, with the current low-mass BDT, the Laurent member Lau1 has a mean pull of
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
reported in Section~\ref{sec:bias_studies}.""",
     r"""At $m_{a}=1~\GeV$, with the current low-mass BDT, the Laurent member Lau1 has a mean pull of
$-0.246\pm0.010$ and is the only function of the whole scan above the threshold; the other two
members give $-0.159$ (Exp1) and $+0.001$ (Pow1), with goodness-of-fit probabilities of 0.33
(Lau1), 0.45 (Exp1) and 0.24 (Pow1). All three are already at the lowest order, so the rule
used at the other mass points, falling back to the next lower order, cannot be applied. Two
alternatives were tested, with the same 10000 toys per function, on the previous version of the
low-mass model, in which only the $m_{a}=1~\GeV$ signal was restricted to the training window and
Lau1 gave $-0.261$. Removing the Laurent family left Exp1 and Pow1 and moved the excess to Exp1
($-0.267$, with $-0.035$ for Pow1), while the expected limit became 3\% weaker. Allowing higher
orders of the exponential and power-law families selected Lau1, Exp1 and Pow3, and Pow3 gave
$-0.317$ (Lau1 $-0.259$). Neither alternative removed the excess, so the
envelope of Table~\ref{tab:bkg_true} is kept at $m_{a}=1~\GeV$, and the residual bias is
reported in Section~\ref{sec:bias_studies}."""),
])

edit("Sec-07-05-Biasstudy.tex", [
    (r"""$-0.261\pm0.010$; at this mass point Exp1 gives $-0.173$ and Pow1 $+0.004$.""",
     r"""$-0.246\pm0.010$; at this mass point Exp1 gives $-0.159$ and Pow1 $+0.001$."""),
    (r"""function of the scan lies inside the threshold except Lau1 at $m_{a}=1~\GeV$ ($-0.261$,""",
     r"""function of the scan lies inside the threshold except Lau1 at $m_{a}=1~\GeV$ ($-0.246$,"""),
])

edit("Sec-08-01-Parameterizedsignalmodel.tex", [
    (r"""its working point the events at $134<m_{\ell\ell\gamma\gamma}<137~\GeV$ are 0.3\% of the selected
signal, 44\% of them have $\Delta R(\gamma,\gamma)<0.05$ and their median photon-pair $p_{T}$
response is 1.21, compared with 2\% and 1.00 for the events in the peak""",
     r"""its working point the events at $134<m_{\ell\ell\gamma\gamma}<137~\GeV$ are 0.2\% of the selected
signal, 10\% of them have $\Delta R(\gamma,\gamma)<0.05$ and their median photon-pair $p_{T}$
response is 1.02, compared with 1\% and 1.00 for the events in the peak"""),
    (r"""(Figures~\ref{fig:ma_true_eras} and~\ref{fig:mA1_doublepeak}), which holds 2.8\% of the selected
signal.""",
     r"""(Figures~\ref{fig:ma_true_eras} and~\ref{fig:mA1_doublepeak}), which holds 0.5\% of the selected
signal."""),
    (r"""high side of $m_{\ell\ell\gamma\gamma}$; with the current low-mass BDT they are 4\% of the
    selected signal.}""",
     r"""high side of $m_{\ell\ell\gamma\gamma}$; with the current low-mass BDT they are 1\% of the
    selected signal.}"""),
    (r"""the events it describes. Between 94.5\% and 97.2\% of the simulated yield,""",
     r"""the events it describes. Between 95.9\% and 98.8\% of the simulated yield,"""),
    (r"""degree of freedom between $0.21$ and $0.93$. The effective resolution ranges
from $1.43$ to $2.74~\GeV$ depending on era and channel, the spread reflecting""",
     r"""degree of freedom between $0.24$ and $0.69$. The effective resolution ranges
from $1.33$ to $2.70~\GeV$ depending on era and channel, the spread reflecting"""),
    (r"""Even inside this window the $m_{a}=1~\GeV$ shape is not a pure resolution peak: in nine of
the ten era and channel combinations at least one additional Gaussian component (two in the
electron channel of 2024) is retained alongside the narrow core to absorb the residual tail, and a
single Gaussian is sufficient only in the electron channel of 2023postBPix. In all ten
combinations the model reproduces the fraction of the in-window yield below $123~\GeV$ to within
6.2 percentage points; the largest difference is in the electron channel of 2023preBPix, where the
model places $18.0\%$ of its probability below $123~\GeV$ and the simulated sample $24.2\%$.""",
     r"""Even inside this window the $m_{a}=1~\GeV$ shape is not a pure resolution peak: in nine of
the ten era and channel combinations a second Gaussian component is retained alongside the
narrow core to absorb the residual tail, and a single Gaussian is sufficient only in the electron
channel of 2023postBPix. In all ten combinations the model reproduces the fraction of the in-window
yield below $123~\GeV$ to within 5.0 percentage points; the largest difference is in the electron
channel of 2023preBPix, where the model places $19.8\%$ of its probability below $123~\GeV$ and
the simulated sample $24.8\%$."""),
])

edit("Sec-10-00-ResultAndInterpretation.tex", [
    (r"""The expected limits are in the range of 2--10~fb over $2 \le m_{a} \le 30~\GeV$, rising to about 28~fb at $m_{a}=1~\GeV$""",
     r"""The expected limits are in the range of 2--10~fb over $2 \le m_{a} \le 30~\GeV$, rising to about 30~fb at $m_{a}=1~\GeV$"""),
    (r"""where they are weaker than the Run 2 ones by factors of about 1.6, 2.0 and 2.0 at $m_{a}=1$, 2 and 3~\GeV.""",
     r"""where they are weaker than the Run 2 ones by factors of about 1.7, 2.0 and 2.0 at $m_{a}=1$, 2 and 3~\GeV."""),
    (r"""about 28~fb at $m_{a}=1~\GeV$, with the best expected sensitivity at""",
     r"""about 30~fb at $m_{a}=1~\GeV$, with the best expected sensitivity at"""),
    (r"""expected limits are instead weaker: 28.3, 9.5 and 5.7~fb at $m_{a}=1$, 2 and 3~\GeV,
compared with 17.8, 4.7 and 2.8~fb in Run 2, i.e.\ larger by factors of about 1.6, 2.0
and 2.0.""",
     r"""expected limits are instead weaker: 30.5, 9.7 and 5.5~fb at $m_{a}=1$, 2 and 3~\GeV,
compared with 17.8, 4.7 and 2.8~fb in Run 2, i.e.\ larger by factors of about 1.7, 2.0
and 2.0."""),
    (r"""spectrum. The low-mass BDT therefore uses only five inputs, without shower-shape variables
and without the photon transverse momenta (Table~\ref{tab:bdt_input_lowmass}), and its working""",
     r"""spectrum. The low-mass BDT therefore uses only five inputs, without shower-shape variables
and without the photon transverse momenta (Table~\ref{tab:bdt_input_lowmass}), its signal is
trained only inside the window $115<m_{\ell\ell\gamma\gamma}<133~\GeV$, and its working"""),
    (r"""They keep 11\%, 34\% and 45\% of the preselected signal, compared with 56--77\% at""",
     r"""They keep 12\%, 36\% and 46\% of the preselected signal, compared with 56--77\% at"""),
    (r"""only 26--47\% of the selected signal entered the model. With the current low-mass BDT they are 4\%
of the selected signal, and 94.5--97.2\% of it, depending on era and lepton channel, falls inside""",
     r"""only 26--47\% of the selected signal entered the model. With the current low-mass BDT they are 1\%
of the selected signal, and 95.9--98.8\% of it, depending on era and lepton channel, falls inside"""),
])
