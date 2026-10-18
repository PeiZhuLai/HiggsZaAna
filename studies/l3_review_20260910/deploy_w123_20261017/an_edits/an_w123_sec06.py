#!/usr/bin/env python3
"""AN-25-172 Sec 6.0, 6.1, 6.3, 6.4, 6.5 for the low-mass BDT low_relE2_w123 (relE2 inputs, signal of m_a = 1, 2, 3 all
trained inside 115-133 GeV; production since 2026-10-17). Numbers (deploy_w123_20261017/an_numbers/):
  model meta: test AUC 0.9375, KS p sig 0.502 / bkg 0.243
  Plot/output/latexTable_sculptR.txt: R 1.30 / 1.27 / 1.26 (relE2 1.19 / 1.23 / 1.30)
  pseudodata_closure/closure_results.json (w123_numbers.txt): 20/30 in +-1 sigma, 30/30 in +-2 sigma, below: 5 7 12 22 23,
    above: 2 25 28 29 30; mA1 rhat -0.150 -0.105/+0.129 Z 0 (obs 0.234, -1 sigma 0.217); mA2 +0.062 Z 1.32; mA3 +0.004 Z 0.16
  closure without training DY, own signal model, xsec card (chainlow_relE2w123/logs/chainlow.log): mA1 Z 0.24 rhat +0.040 -0.157/+0.181
  working-point scan (lowcand/summary_c3.txt, training DY removed, previous signal model) and bootstrap (lowcand/boot_low_relE2_w123_*):
    mA1 x1 +0.67 +- 1.01, x0.5 +0.35 +- 0.99; mA2 x0.5 +0.02 +- 0.09; mA3 x0.5 -0.01 +- 0.05
  Plot/plots/mAmigratedMatrix/mAmigration_matrix_Run3.pdf diagonal: mA1 99.8, mA2 99.3, mA3 99.1 %
  Plot/output/latexTable_significanceVcat.txt: mA1 +2.0 %, at most 0.4 % elsewhere (mA2)
  Plot/plots/mA1_doublepeak/summary.txt: second m_gg peak 0.22 / 45.98 = 0.5 %
Sec 6.0 also corrects the training-range sentence: the high-mass signal is trained on 95-180 GeV, the low-mass signal on
115-133 GeV only. Backups: <file>.bak_w123_20261017"""
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


edit("Sec-06-00-BDTOptimisation.tex", [
    (r"""The signal training uses the signal region, $115~\GeV < m_{\ell^{+}\ell^{-}\gamma\gamma} < 135~\GeV$, while the background training uses the full mass range, $95~\GeV < m_{\ell^{+}\ell^{-}\gamma\gamma} < 180~\GeV$.""",
     r"""The background training uses the full mass range, $95~\GeV < m_{\ell^{+}\ell^{-}\gamma\gamma} < 180~\GeV$. The signal training of the high-mass model (see below) uses the same range, while the low-mass model uses the signal only within $115~\GeV < m_{\ell^{+}\ell^{-}\gamma\gamma} < 133~\GeV$."""),
])

edit("Sec-06-01-Inputvariable.tex", [
    (r"""in its training the $m_{a}=1~\GeV$ signal enters only within $115 < m_{\ell\ell\gamma\gamma} < 133~\GeV$, with its weight rescaled to the full $m_{a}=1~\GeV$ weight;""",
     r"""in its training the signal of each of the three mass points enters only within $115 < m_{\ell\ell\gamma\gamma} < 133~\GeV$, with its weight rescaled to the full weight of that mass point;"""),
])

edit("Sec-06-03-BDTevaluation.tex", [
    (r"""while the low-mass model reaches $0.940$ (both evaluated with the nominal simulation weights, without the sideband reweighting, and for the low-mass model with the $m_{a}=1~\GeV$ signal not restricted to the training window);""",
     r"""while the low-mass model reaches $0.938$ (both evaluated with the nominal simulation weights, without the sideband reweighting, and for the low-mass model with the signal not restricted to the training window);"""),
    (r"""and for the low-mass model the $m_{a}=1~\GeV$ signal enters only within $115<m_{\ell\ell\gamma\gamma}<133~\GeV$. The resulting p-values are $0.46$ (signal) and $0.49$ (background) for the high-mass model, and $0.62$ (signal) and $0.30$ (background) for the low-mass model.""",
     r"""and for the low-mass model the signal at $m_{a}=1$, 2 and $3~\GeV$ enters only within $115<m_{\ell\ell\gamma\gamma}<133~\GeV$. The resulting p-values are $0.46$ (signal) and $0.49$ (background) for the high-mass model, and $0.50$ (signal) and $0.24$ (background) for the low-mass model."""),
])

edit("Sec-06-04-BDTperformancecheck.tex", [
    (r"""$m_{a}=1~\GeV$ signal; it replaces the leading-photon $p_{T}/m_{\ell\ell\gamma\gamma}$ by the
relative energy uncertainty of the sub-leading photon. At its working points it yields
$R(m_{a})=1.19$, $1.23$, and $1.30$ at $m_{a}=1$, $2$, and $3~\GeV$
(Table~\ref{tab:sculpt_R}), compared with 1.00, 1.18 and 1.16 for the previous model. The""",
     r"""$m_{a}=1~\GeV$ signal; it replaces the leading-photon $p_{T}/m_{\ell\ell\gamma\gamma}$ by the
relative energy uncertainty of the sub-leading photon. The signal of all three mass points enters
its training only within $115<m_{\ell\ell\gamma\gamma}<133~\GeV$, so that the three points follow
the same training rule. At its working points it yields
$R(m_{a})=1.30$, $1.27$, and $1.26$ at $m_{a}=1$, $2$, and $3~\GeV$
(Table~\ref{tab:sculpt_R}), compared with 1.19, 1.23 and 1.30 for the same inputs with only the
$m_{a}=1~\GeV$ signal restricted to that window, and 1.00, 1.18 and 1.16 for the previous model. The"""),
    ("      1 & 1.19 & 4  & 0.64 \\\\\n", "      1 & 1.30 & 4  & 0.64 \\\\\n"),
    ("      2 & 1.23 & 5  & 0.38 \\\\\n", "      2 & 1.27 & 5  & 0.38 \\\\\n"),
    ("      3 & 1.30 & 6  & 0.77 \\\\\n", "      3 & 1.26 & 6  & 0.77 \\\\\n"),
    (r"""(Table~\ref{tab:closure_lowmass_scan}); the resulting values are $R=1.19$, $1.23$, and $1.30$""",
     r"""(Table~\ref{tab:closure_lowmass_scan}); the resulting values are $R=1.30$, $1.27$, and $1.26$"""),
    (r"""$\pm1\sigma$ band. Of the ten points outside the $\pm1\sigma$ band, six ($m_{a}=1$, 5, 7, 12, 22
and $23~\GeV$) have a stronger-than-expected limit, i.e.\ a deficit, and four
($m_{a}=25$, 28, 29 and $30~\GeV$) a weaker one.""",
     r"""$\pm1\sigma$ band. Of the ten points outside the $\pm1\sigma$ band, five ($m_{a}=5$, 7, 12, 22
and $23~\GeV$) have a stronger-than-expected limit, i.e.\ a deficit, and five
($m_{a}=2$, 25, 28, 29 and $30~\GeV$) a weaker one."""),
    (r"""low-mass regime most susceptible to sculpting shows no fake signal: the local significances
are $0.00\sigma$, $0.97\sigma$ and $0.46\sigma$ at $m_{a}=1$, 2 and $3~\GeV$, with
$\hat{\mu}=+0.044$ and $+0.013$ at $m_{a}=2$ and $3~\GeV$, the working points being those
selected with this test (see below). At $m_{a}=1~\GeV$ the pseudo-data show a deficit in the
signal window, $\hat{\mu}=-0.32^{+0.11}_{-0.10}$, with the observed limit below the
$-1\sigma$ band. The simulated background used here includes the events used to train the
BDT, which the classifier rejects more efficiently than independent events; repeated without
the training events, the test gives $\hat{\mu}=-0.16^{+0.16}_{-0.14}$ and a local significance
of zero at this mass point.""",
     r"""low-mass regime most susceptible to sculpting shows no significant fake signal: the local
significances are $0.00\sigma$, $1.32\sigma$ and $0.16\sigma$ at $m_{a}=1$, 2 and $3~\GeV$, with
$\hat{\mu}=+0.062^{+0.053}_{-0.047}$ and $+0.004$ at $m_{a}=2$ and $3~\GeV$, the working points being
those selected with this test (see below). At $m_{a}=1~\GeV$ the pseudo-data show a small deficit
in the signal window, $\hat{\mu}=-0.15^{+0.13}_{-0.11}$, with the observed limit inside the
$\pm1\sigma$ band. The simulated background used here includes the events used to train the
BDT, which the classifier rejects more efficiently than independent events; repeated without
the training events, the test gives $\hat{\mu}=+0.04^{+0.18}_{-0.16}$ and a local significance of
$0.24\sigma$ at this mass point."""),
    (r"""At $m_{a}=1~\GeV$ only the cuts at 1 and 0.5 times the data yield pass; the number of simulated
background events in the Higgs-peak window is small at these cuts, so 30 bootstrap replicas of the
simulation, in which every simulated event is weighted by a Poisson(1) random number and the full
closure procedure is repeated, were used to compare the two: the mean $\hat{\mu}$ of the replicas is
$+0.28$ with a spread of $0.47$ at the looser cut and $+0.03$ with a spread of $0.69$ at the tighter
one, in units of the signal yield of the previous model. At $m_{a}=2~\GeV$ four of the five cuts
pass, and at $m_{a}=3~\GeV$ the two loosest cuts fail the goodness-of-fit requirement, the
first-order Bernstein member giving a probability below 0.01. At all three mass points the adopted
working point is the tightest cut of the scan, 0.98036, 0.98661 and 0.99194, which keeps 156, 267
and 176 data events.""",
     r"""At $m_{a}=1~\GeV$ only the tightest cut, at 0.5 times the data yield, passes. The number of
simulated background events in the Higgs-peak window is small at the tight cuts, so 30 bootstrap
replicas of the simulation, in which every simulated event is weighted by a Poisson(1) random
number and the full closure procedure is repeated, were used to check the result: the mean
$\hat{\mu}$ of the replicas is $+0.35$ with a spread of $0.99$ at the adopted cut, and $+0.67$ with
a spread of $1.01$ at the cut with the data yield of the previous working point, in units of the
signal yield of the previous model. At $m_{a}=2~\GeV$ four of the five cuts pass, the cut at the
previous data yield failing the closure requirement, and at $m_{a}=3~\GeV$ the loosest cut fails
the goodness-of-fit requirement, the first-order Bernstein member giving a probability below 0.01;
at the adopted cuts the bootstrap gives $+0.02\pm0.09$ and $-0.01\pm0.05$ at $m_{a}=2$ and
$3~\GeV$. At all three mass points the adopted working point is the tightest cut of the scan,
0.98115, 0.98837 and 0.99226, which keeps 167, 267 and 169 data events."""),
    (r"""      1 & $\times2$   & 0.96410 & 552 & 0.03 & $+0.49$ & 1.11 \\
        & $\times1.5$ & 0.96837 & 428 & 0.18 & $+0.45$ & 1.17 \\
        & $\times1$   & 0.97363 & 296 & 0.15 & $+0.25$ & 0.75 \\
        & $\times0.7$ & 0.97699 & 215 & 0.25 & $+0.81$ & 2.34 \\
        & $\mathbf{\times0.5}$ & \textbf{0.98036} & \textbf{156} & \textbf{0.35} & $\mathbf{-0.16}$ & \textbf{0.00} \\ \hline
      2 & $\times2$   & 0.96381 & 997 & 0.30 & $+0.032$ & 0.50 \\
        & $\times1.5$ & 0.97130 & 776 & 0.35 & $+0.029$ & 0.49 \\
        & $\times1$   & 0.97830 & 537 & 0.37 & $+0.040$ & 0.80 \\
        & $\times0.7$ & 0.98300 & 375 & 0.70 & $+0.056$ & 1.25 \\
        & $\mathbf{\times0.5}$ & \textbf{0.98661} & \textbf{267} & \textbf{0.83} & $\mathbf{+0.021}$ & \textbf{0.55} \\ \hline
      3 & $\times2$   & 0.96921 & 695 & $<0.01$ & $+0.069$ & 1.48 \\
        & $\times1.5$ & 0.97714 & 529 & $<0.01$ & $+0.003$ & 0.08 \\
        & $\times1$   & 0.98421 & 362 & 0.10 & $-0.017$ & 0.00 \\
        & $\times0.7$ & 0.98911 & 250 & 0.32 & $-0.016$ & 0.00 \\
        & $\mathbf{\times0.5}$ & \textbf{0.99194} & \textbf{176} & \textbf{0.45} & $\mathbf{+0.013}$ & \textbf{0.62} \\ \hline""",
     r"""      1 & $\times2$   & 0.96548 & 562 & 0.01 & $+1.37$ & 2.85 \\
        & $\times1.5$ & 0.96914 & 442 & 0.08 & $+0.58$ & 1.38 \\
        & $\times1$   & 0.97518 & 303 & 0.16 & $+0.61$ & 1.59 \\
        & $\times0.7$ & 0.97847 & 223 & 0.13 & $+0.54$ & 1.58 \\
        & $\mathbf{\times0.5}$ & \textbf{0.98115} & \textbf{167} & \textbf{0.24} & $\mathbf{+0.28}$ & \textbf{0.91} \\ \hline
      2 & $\times2$   & 0.96623 & 995 & 0.36 & $+0.016$ & 0.25 \\
        & $\times1.5$ & 0.97414 & 756 & 0.34 & $+0.019$ & 0.34 \\
        & $\times1$   & 0.98041 & 526 & 0.44 & $+0.090$ & 1.84 \\
        & $\times0.7$ & 0.98509 & 371 & 0.70 & $-0.004$ & 0.00 \\
        & $\mathbf{\times0.5}$ & \textbf{0.98837} & \textbf{267} & \textbf{0.70} & $\mathbf{+0.018}$ & \textbf{0.47} \\ \hline
      3 & $\times2$   & 0.97028 & 689 & $<0.01$ & $-0.026$ & 0.00 \\
        & $\times1.5$ & 0.97783 & 530 & 0.02 & $+0.004$ & 0.15 \\
        & $\times1$   & 0.98516 & 352 & 0.22 & $-0.010$ & 0.00 \\
        & $\times0.7$ & 0.98907 & 253 & 0.17 & $-0.028$ & 0.00 \\
        & $\mathbf{\times0.5}$ & \textbf{0.99226} & \textbf{169} & \textbf{0.32} & $\mathbf{+0.001}$ & \textbf{0.04} \\ \hline"""),
])

edit("Sec-06-05-Eventcategorization.tex", [
    (r"""a gain of at most about 2\% (1.9\% at $m_{a}=1~\GeV$, at most 0.5\% elsewhere)""",
     r"""a gain of at most about 2\% (2.0\% at $m_{a}=1~\GeV$, at most 0.4\% elsewhere)"""),
    (r"""The correct-$m_{a}$ probability ranges from 55\% at $m_{a}=8~\GeV$ to 99\% at
$m_{a}=3~\GeV$: it is above 90\% at $m_{a}=2$, 3 and $15$--$30~\GeV$, and between
55\% and 78\% at $m_{a}=1$ and 4--10~\GeV, where the signal selected by a given""",
     r"""The correct-$m_{a}$ probability ranges from 55\% at $m_{a}=8~\GeV$ to 99.8\% at
$m_{a}=1~\GeV$: it is above 90\% at $m_{a}=1$, 2, 3 and $15$--$30~\GeV$, and between
55\% and 78\% at $m_{a}=4$--10~\GeV, where the signal selected by a given"""),
    (r"""showers and an overestimated photon energy form a small second, displaced peak, with 2.8\% of
the selected signal""",
     r"""showers and an overestimated photon energy form a small second, displaced peak, with 0.5\% of
the selected signal"""),
    (r"""there a small second peak at 1.2--1.9~\GeV, with 2.8\% of the selected signal, is formed""",
     r"""there a small second peak at 1.2--1.9~\GeV, with 0.5\% of the selected signal, is formed"""),
])
