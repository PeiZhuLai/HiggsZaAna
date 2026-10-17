#!/usr/bin/env python3
"""AN-25-172 Sec 6.4 for the low-mass BDT with the sub-leading-photon sigma_E/E input (low_swapPt_relE2).
Numbers:
  sculpting ratio at the working points (Plot/output/latexTable_sculptR.txt): 1.19 / 1.23 / 1.30 (previous 1.00 / 1.18 / 1.16)
  production closure, all 30 points (pseudodata_closure/closure_results.json, logs/closure_extract.log): mA1 obs 0.171 below
    -1sigma (inside -2sigma 0.144), rhat -0.320 -0.097/+0.113, Z 0; mA2 Z 0.97 rhat +0.044; mA3 Z 0.46 rhat +0.013;
    20/30 within +-1sigma, 30/30 within +-2sigma; down: 1,5,7,12,22,23; up: 25,28,29,30
  closure without the training events, current signal model, xsec card (chainlow_relE2/closure_final_summary.txt): mA1
    Z 0, rhat -0.157 -0.139/+0.163
  working-point scan (sandbox closure_scan/cand_low_swapPt_relE2_drop_mA*, training DY removed, previous signal model)
  bootstrap at mA1 (30 replicas; docdir): x1 +0.28 +- 0.47, x0.5 +0.03 +- 0.69 (production at 0.963: -0.46 +- 0.39)
Also corrects the statement that the goodness-of-fit and bias study establish the absence of sculpting (the bias study
generates its toys from the background functions and cannot see a peak sculpted into the true background).
Backup: Sec-06-04-BDTperformancecheck.tex.bak_relE2_20261016"""
AN = "/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172/"
fn = "Sec-06-04-BDTperformancecheck.tex"
s = open(AN + fn).read()


def rep(a, b):
    global s
    assert s.count(a) == 1, a[:100]
    s = s.replace(a, b)


def rep_span(start, end, new):
    global s
    i = s.index(start)
    j = s.index(end, i) + len(end)
    s = s[:i] + new + s[j:]


rep(r"""bias on the result. Whether any residual sculpting actually biases the measurement is
established by the background-function goodness-of-fit at each working point and by the
signal-injection bias study over the function envelope (Section~\ref{sec:bias_studies}).
The sculpting ratio is used to steer the working-point choice toward $R\approx1$ as a
conservative precaution.""",
    r"""bias on the result. Whether any residual sculpting fakes a signal is tested directly by
the background-only pseudo-data closure test of Section~\ref{sec:sculpt_closure}; the
signal-injection bias study of Section~\ref{sec:bias_studies} generates its toys from the
background functions themselves and therefore cannot reveal a peak sculpted into the true
background. The sculpting ratio was used to steer the working points of the previous
low-mass model toward $R\approx1$ as a conservative precaution.""")

rep_span(r"""The sculpting ratio is the figure of merit used to select the reduced input set of the
low-mass model.""", r"""$R=1.18$.""",
         r"""The sculpting ratio was the figure of merit used to select the input set of the previous
low-mass model. When the full 16-variable set is used at $m_{a}=1,2,3~\GeV$, the
sub-leading-photon observables and the $\sigma_{i\eta i\eta}^{5\times 5}$ shower-shape
variables, which are poorly modeled and strongly correlated with $m_{\ell\ell\gamma\gamma}$
in the merged-photon regime, drive $R(m_{a})$ well above unity. The previous low-mass model
used five inputs including the leading-photon $p_{T}/m_{\ell\ell\gamma\gamma}$, with working
points at $R\approx1$. The current model (Table~\ref{tab:bdt_input_lowmass}) was selected in a
dedicated search over input sets, requiring no excess of the simulated background in the
Higgs-peak window on a simulated test sample not used in the training and a single-peaked
$m_{a}=1~\GeV$ signal; it replaces the leading-photon $p_{T}/m_{\ell\ell\gamma\gamma}$ by the
relative energy uncertainty of the sub-leading photon. At its working points it yields
$R(m_{a})=1.19$, $1.23$, and $1.30$ at $m_{a}=1$, $2$, and $3~\GeV$
(Table~\ref{tab:sculpt_R}), compared with 1.00, 1.18 and 1.16 for the previous model. The
working points are no longer chosen at $R\approx1$ but with the background-only closure test
(Section~\ref{sec:sculpt_closure}), which tests directly whether the residual $R>1$ fakes a
signal; at these working points it does not. $R$ is computed on the inclusive simulation,
which includes the events used to train the BDT.""")

rep(r"""Figures~\ref{fig:bdt_sculpt} and~\ref{fig:bdt_mass_shapes}. As noted above, the sculpting
ratio serves as a conservative diagnostic to guide the working-point choice; the definitive
demonstration that the residual sculpting does not bias the final signal extraction is
provided by the background-function goodness-of-fit and the signal-injection bias study of
Section~\ref{sec:bias_studies}, which use the data-driven background model rather than the
simulation.""",
    r"""Figures~\ref{fig:bdt_sculpt} and~\ref{fig:bdt_mass_shapes}. As noted above, the sculpting
ratio is a diagnostic; whether the residual sculpting fakes a signal is tested directly by
the background-only closure test of Section~\ref{sec:sculpt_closure}.""")

rep_span(r"""The BDT score threshold of the low-mass model is fixed from a joint scan of the""",
         r"""of Section~\ref{sec:sculpt_closure}.""",
         r"""Figure~\ref{fig:scan_score_R} shows the sculpting ratio $R(m_{a})$ and the expected Asimov
significance $Z$ as a function of the applied BDT score cut for the low-mass model at
$m_{a}=1,2,3~\GeV$. Each point corresponds to one candidate cut, with the color encoding the
associated significance. Tightening the cut increases $Z$ but, in the merged-photon regime,
also drives $R(m_{a})$ above unity. The working points (red stars) are chosen with the
background-only closure test of Section~\ref{sec:sculpt_closure}
(Table~\ref{tab:closure_lowmass_scan}).""")

rep_span(r"""    and the red star marks the selected working point, chosen to keep $R(m_{a})$ as""",
         r"""    expected significance. The""",
         r"""    and the red star marks the adopted working point, chosen with the background-only
    pseudo-data closure test of Section~\ref{sec:sculpt_closure}
    (Table~\ref{tab:closure_lowmass_scan}); the resulting values are $R=1.19$, $1.23$, and $1.30$
    at $m_{a}=1$, $2$, and $3~\GeV$. The""")

rep(r"""$1<m_{a}<30~\GeV$ the observed limit obtained on the background-only pseudo-data tracks
the Asimov expectation: all 30 mass points lie within the $\pm2\sigma$ band and 21 within the
$\pm1\sigma$ band. Of the nine points outside the $\pm1\sigma$ band, five ($m_{a}=5$, 7, 12, 22
and $23~\GeV$) have a stronger-than-expected limit, i.e.\ a downward fluctuation, and four
($m_{a}=25$, 28, 29 and $30~\GeV$) a weaker one.""",
    r"""$1<m_{a}<30~\GeV$ the observed limit obtained on the background-only pseudo-data tracks
the Asimov expectation: all 30 mass points lie within the $\pm2\sigma$ band and 20 within the
$\pm1\sigma$ band. Of the ten points outside the $\pm1\sigma$ band, six ($m_{a}=1$, 5, 7, 12, 22
and $23~\GeV$) have a stronger-than-expected limit, i.e.\ a deficit, and four
($m_{a}=25$, 28, 29 and $30~\GeV$) a weaker one.""")

rep(r"""pseudo-data show a deficit in the signal window; the local significance there is zero. The
low-mass regime most susceptible to sculpting is clean: the local significances are
$0.00\sigma$, $1.12\sigma$ and $0.96\sigma$ at $m_{a}=1$, 2 and $3~\GeV$, the
$m_{a}=3~\GeV$ working point being the one selected with this test (see below). Any residual enhancement appears only in the simulation, whereas the final
background is the data-driven fit to the data over the full
$95<m_{\ell\ell\gamma\gamma}<180~\GeV$ range, including the signal window, and it is in any
case covered by the background-function goodness-of-fit and the signal-injection bias study
of Section~\ref{sec:bias_studies}, so it does not bias the
measurement.""",
    r"""pseudo-data show a deficit in the signal window; the local significance there is zero. The
low-mass regime most susceptible to sculpting shows no fake signal: the local significances
are $0.00\sigma$, $0.97\sigma$ and $0.46\sigma$ at $m_{a}=1$, 2 and $3~\GeV$, with
$\hat{\mu}=+0.044$ and $+0.013$ at $m_{a}=2$ and $3~\GeV$, the working points being those
selected with this test (see below). At $m_{a}=1~\GeV$ the pseudo-data show a deficit in the
signal window, $\hat{\mu}=-0.32^{+0.11}_{-0.10}$, with the observed limit below the
$-1\sigma$ band. The simulated background used here includes the events used to train the
BDT, which the classifier rejects more efficiently than independent events; repeated without
the training events, the test gives $\hat{\mu}=-0.16^{+0.16}_{-0.14}$ and a local significance
of zero at this mass point. Any residual enhancement appears only in the simulation, whereas
the final background is the data-driven fit to the data over the full
$95<m_{\ell\ell\gamma\gamma}<180~\GeV$ range, including the signal window.""")

rep_span(r"""\paragraph{Working point at $m_{a}=3~\GeV$.}""",
         r"""      0.996 & 1.57 & 23.5 &  75 & 0.34--0.52 & $+0.032$ & 2.2 \\ \hline
    \end{tabular}
  \end{center}
\end{table}""",
         r"""\paragraph{Working points at $m_{a}=1$, 2 and $3~\GeV$.}
The working points of the low-mass model are chosen with this closure test. For each of the three
mass points, five BDT score cuts are scanned, defined by the number of selected data events in the
$m_{\ell\ell\gamma\gamma}$ sidebands, which is 2, 1.5, 1, 0.7 and 0.5 times the number selected at
the working point of the previous low-mass model. At each cut the background envelope is fitted to
the real data (100 toys for the goodness-of-fit, as in the nominal fit) and to the background-only
pseudo-data, built as above but without the simulated events used in the training of the BDT. A cut
is accepted if the closure significance is below $1\sigma$ and every member of the data envelope has
a goodness-of-fit probability above 0.01, and among the accepted cuts the one with the largest
expected sensitivity is adopted. The results are listed in Table~\ref{tab:closure_lowmass_scan}.
At $m_{a}=1~\GeV$ only the cuts at 1 and 0.5 times the data yield pass; the number of simulated
background events in the Higgs-peak window is small at these cuts, so 30 bootstrap replicas of the
simulation, in which every simulated event is weighted by a Poisson(1) random number and the full
closure procedure is repeated, were used to compare the two: the mean $\hat{\mu}$ of the replicas is
$+0.28$ with a spread of $0.47$ at the looser cut and $+0.03$ with a spread of $0.69$ at the tighter
one, in units of the signal yield of the previous model. At $m_{a}=2~\GeV$ four of the five cuts
pass, and at $m_{a}=3~\GeV$ the two loosest cuts fail the goodness-of-fit requirement, the
first-order Bernstein member giving a probability below 0.01. At all three mass points the adopted
working point is the tightest cut of the scan, 0.98036, 0.98661 and 0.99194, which keeps 156, 267
and 176 data events.

\begin{table}[htbp]
  \begin{center}
    \topcaption{Working-point scan of the low-mass model at $m_{a}=1$, 2 and $3~\GeV$. For each
    BDT score cut: the number of selected data events relative to the working point of the
    previous low-mass model and in absolute terms, the smallest goodness-of-fit probability among
    the members of the background envelope fitted to the data, and the background-only
    pseudo-data closure result (without the simulated events used in the training; best-fit
    $\hat{\mu}$ in units of the signal yield of the previous model, and local significance).
    The adopted working points are in bold.}
    \label{tab:closure_lowmass_scan}
    \small
    \begin{tabular}{|c|c|c|c|c|c|c|} \hline
      $m_{a}$ [\GeV] & data yield & BDT cut & $N_{\mathrm{data}}$ & min.\ data GOF $p$-value & closure $\hat{\mu}$ & closure $Z$ \\ \hline
      1 & $\times2$   & 0.96410 & 552 & 0.03 & $+0.49$ & 1.11 \\
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
        & $\mathbf{\times0.5}$ & \textbf{0.99194} & \textbf{176} & \textbf{0.45} & $\mathbf{+0.013}$ & \textbf{0.62} \\ \hline
    \end{tabular}
  \end{center}
\end{table}
""")


open(AN + fn, "w").write(s)
print("Sec 6.4 edited")
