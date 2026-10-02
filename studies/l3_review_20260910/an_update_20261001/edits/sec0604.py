PAIRS = [
# --- low-/high-mass R paragraph
(r"""yields $R(m_{a})=0.87$, $0.96$, and $1.24$ at $m_{a}=1$, $2$, and $3~\GeV$, as listed in
Table~\ref{tab:sculpt_R}. The values at $m_{a}=1$ and $2~\GeV$ are close to unity. At
$m_{a}=3~\GeV$ no working point both brings $R$ close to unity and keeps a background model
that describes the data: $R$ reaches about 1.1 only at cuts of 0.94 or looser, where the
selected data sample has about 1900 events or more and no background function passes the
goodness-of-fit test. The $m_{a}=3~\GeV$ working point is therefore set by the
background-only pseudo-data closure scan of Section~\ref{sec:sculpt_closure}, which tests
directly whether the residual $R>1$ fakes a signal.""",
r"""yields $R(m_{a})=1.00$, $1.18$, and $1.16$ at $m_{a}=1$, $2$, and $3~\GeV$, as listed in
Table~\ref{tab:sculpt_R}. The value at $m_{a}=1~\GeV$ is at unity, and those at $m_{a}=2$ and
$3~\GeV$ exceed it by less than 20\%. At $m_{a}=3~\GeV$ the cuts scanned between 0.980 and
0.996 give $R$ between 1.16 and 1.57 (Table~\ref{tab:closure_mA3_scan}), the lowest value being
at the adopted working point. The $m_{a}=3~\GeV$ working point is set by the
background-only pseudo-data closure scan of Section~\ref{sec:sculpt_closure}, which tests
directly whether the residual $R>1$ fakes a signal."""),
(r"""criterion adopted for the low-mass points; on the final samples the same threshold gives
$R=0.96$. For the high-mass model the full input set yields
$R(m_{a})$ in the range $0.34$--$1.30$ across the generated grid. The values well below unity
at $m_{a}=4$--$6~\GeV$ are not significant: after the tight working points (0.995) the effective
number of simulated background events in the peak window is only 1--5, and a bootstrap of the
simulation gives an uncertainty of 0.22--0.32 on $R$ there. Where the simulation constrains $R$,
the largest enhancements are $R=1.12\pm0.08$ and $1.30\pm0.07$ at $m_{a}=25$ and $30~\GeV$;
at these two masses $R$ exceeds unity over the whole range of scanned BDT thresholds, so it is a
property of the high-score region rather than of the chosen working point""",
r"""criterion adopted for the low-mass points; on the final samples the same threshold gives
$R=1.18$. For the high-mass model the full input set yields
$R(m_{a})$ in the range $0.38$--$1.25$ across the generated grid. The values well below unity
at $m_{a}=4$--$7~\GeV$ are not significant: after the tight working points (0.990--0.995) only a
few effective simulated background events remain in the peak window, and a bootstrap of the
simulation at these masses, performed in the previous iteration of the study, gave an
uncertainty of 0.17--0.32 on $R$. The largest enhancements are $R=1.25$ at $m_{a}=30~\GeV$ and
$R=1.04$ at $m_{a}=25~\GeV$; at these two masses $R$ exceeds unity over the whole range of
scanned BDT thresholds (1.06--1.24 and 1.23--1.40 on the scan grid), so it is a
property of the high-score region rather than of the chosen working point"""),
# --- AUTO table
(r"""    peak window relative to the sidebands. The largest value,
    $R(m_{a}=30~\GeV)=1.30$, is not a concern: $R$ is evaluated on the simulated
    background, whereas the final background is the data-driven sideband fit, and a""",
 r"""    peak window relative to the sidebands. The largest value,
    $R(m_{a}=30~\GeV)=1.25$, is not a concern: $R$ is evaluated on the simulated
    background, whereas the final background is the data-driven fit to the data, and a"""),
(r"""    (Fig.~\ref{fig:sculpt_closure}), whose observed limit at $m_{a}=30~\GeV$ lies within
    the $\pm1\sigma$ expected band with a best-fit signal strength compatible with zero,
    showing that the $R>1$ there does not bias the result.}""",
 r"""    (Fig.~\ref{fig:sculpt_closure}), whose observed limit at $m_{a}=30~\GeV$ lies within
    the $\pm2\sigma$ expected band, with a local significance of $1.5\sigma$ and a best-fit
    signal strength of $+0.033^{+0.024}_{-0.022}$, showing that the $R>1$ there does not
    produce a significant fake signal.}"""),
(r"""      1 & 0.87 & 4  & 0.44 \\
      2 & 0.96 & 5  & 0.34 \\
      3 & 1.24 & 6  & 0.53 \\
        &      & 7  & 0.82 \\
        &      & 8  & 1.01 \\
        &      & 9  & 0.99 \\
        &      & 10 & 0.77 \\
        &      & 15 & 1.11 \\
        &      & 20 & 0.98 \\
        &      & 25 & 1.12 \\
        &      & 30 & 1.30 \\ \hline""",
 r"""      1 & 1.00 & 4  & 0.64 \\
      2 & 1.18 & 5  & 0.38 \\
      3 & 1.16 & 6  & 0.77 \\
        &      & 7  & 0.69 \\
        &      & 8  & 1.00 \\
        &      & 9  & 0.88 \\
        &      & 10 & 0.84 \\
        &      & 15 & 0.91 \\
        &      & 20 & 0.93 \\
        &      & 25 & 1.04 \\
        &      & 30 & 1.25 \\ \hline"""),
(r"""$m_{a}=3~\GeV$, whose working point (0.992) is chosen from the pseudo-data closure scan""",
 r"""$m_{a}=3~\GeV$, whose working point (0.988) is chosen from the pseudo-data closure scan"""),
(r"""    resulting values are $R=0.87$, $0.96$, and $1.24$ at $m_{a}=1$, $2$, and $3~\GeV$.
    The $m_{a}=3~\GeV$ working point (0.992) is set by the background-only pseudo-data
    closure scan of Section~\ref{sec:sculpt_closure} (Table~\ref{tab:closure_mA3_scan}),
    because no cut brings $R$ close to unity while keeping a background model that
    describes the data.""",
 r"""    resulting values are $R=1.00$, $1.18$, and $1.16$ at $m_{a}=1$, $2$, and $3~\GeV$.
    The $m_{a}=3~\GeV$ working point (0.988) is set by the background-only pseudo-data
    closure scan of Section~\ref{sec:sculpt_closure} (Table~\ref{tab:closure_mA3_scan}),
    because no cut brings $R$ close to unity."""),
(r"""    working point. At $m_{a}=4$--$6~\GeV$ fewer than five effective simulated background events
    remain in the peak window after the working point, so $R$ is not constrained there (bootstrap
    uncertainty 0.2--0.3). At $m_{a}=25$ and $30~\GeV$ $R$ exceeds unity over the whole scanned range.}""",
 r"""    working point. At $m_{a}=4$--$7~\GeV$ only a few effective simulated background events
    remain in the peak window after the tight working points (0.990--0.995), so $R$ is not
    constrained there. At $m_{a}=25$ and $30~\GeV$ $R$ exceeds unity over the whole scanned range.}"""),
# --- SR-adjacent slices
(r"""signal-region slice that is absent from the adjacent one: the fraction of events in
$120<m_{\ell\ell\gamma\gamma}<130~\GeV$ agrees between the two slices within two standard
deviations at 24 of the 30 points. The largest deviations are a nearly empty 125--130~\GeV bin
of the adjacent slice at $m_{a}=8~\GeV$ (negative-weight events), broad excesses without a local
peak at $m_{a}=19$ and $22~\GeV$ ($\chi^{2}$ $p$-values 0.18 and 0.49 for the two shapes), and a
depletion of the peak window in the signal-region slice at $m_{a}=1$ and $2~\GeV$. At
$m_{a}=25$--$30~\GeV$ both slices have a peak-window fraction of 0.18--0.22, above the value of
0.160 without a BDT requirement, so the $R>1$ of Table~\ref{tab:sculpt_R} is a smooth feature of
the whole high-score region, already present below the working point (for example 0.178 and
0.179 in the slices 0.960--0.980 and 0.980--1.000 at $m_{a}=25~\GeV$).""",
r"""signal-region slice that is absent from the adjacent one: the fraction of events in
$120<m_{\ell\ell\gamma\gamma}<130~\GeV$ agrees between the two slices within two standard
deviations at 26 of the 30 points. All four larger deviations, at $m_{a}=7$, 24, 28 and
$29~\GeV$ ($-2.9$, $-3.2$, $-2.0$ and $-2.6$ standard deviations), are a depletion of the peak
window in the signal-region slice, i.e.\ the opposite of a fake peak; no mass point shows an
enhancement of the signal-region slice by more than two standard deviations (the largest is
$+1.7$ at $m_{a}=16~\GeV$). At $m_{a}=25$--$30~\GeV$ both slices have a peak-window fraction of
0.17--0.23, above the value of about 0.16 without a BDT requirement, so the $R>1$ of
Table~\ref{tab:sculpt_R} is a smooth feature of the whole high-score region, already present
below the working point (for example 0.180 and 0.168 in the slices 0.960--0.980 and
0.980--1.000 at $m_{a}=25~\GeV$)."""),
# --- 30-point closure
(r"""the Asimov expectation: all 30 mass points lie within the $\pm2\sigma$ band and 25 within the
$\pm1\sigma$ band. Of the five points outside the $\pm1\sigma$ band, two ($m_{a}=6$ and
$23~\GeV$) have a stronger-than-expected limit, i.e.\ a downward fluctuation, and three
($m_{a}=26$, 28 and $30~\GeV$) a weaker one. No mass point shows a local significance above
$2\sigma$: the largest values are $1.77\sigma$ at $m_{a}=28~\GeV$ ($\hat{\mu}=+0.042$) and
$1.63\sigma$ at $m_{a}=26~\GeV$ ($\hat{\mu}=+0.039$), and $|\hat{\mu}|\le0.042$ for all
$m_{a}\ge3~\GeV$. At $m_{a}=5$ and $6~\GeV$, where the tight working points leave only 23 and
33 pseudo-data events, the signal-plus-background fit runs into the lower physical boundary
of $\mu$, i.e.\ the pseudo-data show a deficit in the signal window; the local significance
there is zero and $\hat{\mu}$ is not shown. The low-mass regime most susceptible to sculpting
is clean: the local significances are $0.00\sigma$, $0.65\sigma$ and $0.80\sigma$ at
$m_{a}=1$, 2 and $3~\GeV$, the $m_{a}=3~\GeV$ working point being the one selected with this
test (see below).""",
r"""the Asimov expectation: all 30 mass points lie within the $\pm2\sigma$ band and 20 within the
$\pm1\sigma$ band. Of the ten points outside the $\pm1\sigma$ band, six ($m_{a}=5$, 7, 9, 12, 22
and $23~\GeV$) have a stronger-than-expected limit, i.e.\ a downward fluctuation, and four
($m_{a}=25$, 28, 29 and $30~\GeV$) a weaker one. No mass point shows a local significance above
$2\sigma$: the largest values are $1.92\sigma$ at $m_{a}=28~\GeV$ ($\hat{\mu}=+0.045$),
$1.50\sigma$ at $m_{a}=30~\GeV$ ($\hat{\mu}=+0.033$), $1.39\sigma$ at $m_{a}=25~\GeV$
($\hat{\mu}=+0.030$) and $1.37\sigma$ at $m_{a}=29~\GeV$ ($\hat{\mu}=+0.034$), and
$|\hat{\mu}|\le0.045$ for all $m_{a}\ge3~\GeV$ except $m_{a}=5~\GeV$. At $m_{a}=5~\GeV$, where the
tight working point leaves pseudo-data with a sum of weights of only 34 events, the
signal-plus-background fit runs into the lower physical boundary of $\mu$, i.e.\ the
pseudo-data show a deficit in the signal window; the local significance there is zero. The
low-mass regime most susceptible to sculpting is clean: the local significances are
$0.00\sigma$, $1.13\sigma$ and $0.96\sigma$ at $m_{a}=1$, 2 and $3~\GeV$, the
$m_{a}=3~\GeV$ working point being the one selected with this test (see below)."""),
# --- mA3 working point
(r"""fit), so that the working point satisfies both requirements at once: the background model
describes the data, and the residual sculpting does not fake a signal. The results are
listed in Table~\ref{tab:closure_mA3_scan}. At the previous working point (0.988) the
pseudo-data yield $\hat{\mu}=+0.063$ with a local significance of $2.3\sigma$. Looser cuts
are worse: for 0.984 and 0.986 the fake signal grows to 3.3 and $4.1\sigma$, and no
background function describes the data. Among the cuts at which every data-envelope member
has a goodness-of-fit $p$-value above 0.01, 0.992 is the only one whose closure significance
is below $1.5\sigma$ ($0.8\sigma$, $\hat{\mu}=+0.018$), and its expected Asimov significance
is within 2\% of that at 0.988. The working point of $m_{a}=3~\GeV$ is therefore set to
0.992. The expected limit at $m_{a}=3~\GeV$ improves by 4\% relative to 0.988.

The number of simulated background events in the $120<m_{\ell\ell\gamma\gamma}<130~\GeV$
window is small at these tight cuts (effective size 36 at 0.988 and 20 at 0.992), so the
closure result itself fluctuates. This is quantified with 40 bootstrap replicas of the
simulation, in which every simulated event is weighted by a Poisson(1) random number and
the full closure procedure is repeated. The replicas give
$\hat{\mu}=+0.063\pm0.034$ and a significance of $2.2\pm1.1\sigma$ at 0.988, compared with
$\hat{\mu}=+0.027\pm0.027$ and $1.2\pm1.1\sigma$ at 0.992. The fake signal at 0.988 thus
differs from zero by about two standard deviations of the simulation statistics, whereas at
0.992 it is compatible with zero. The same limited simulation statistics also explain why
$R$ does not vary monotonically with the cut in Table~\ref{tab:closure_mA3_scan}.""",
r"""fit), so that the working point satisfies both requirements at once: the background model
describes the data, and the residual sculpting does not fake a signal. The scan was
repeated after the retraining of the BDT that followed the re-production of the 2024 DY+jets
samples with the overlap removal of Section~\ref{sec:overlap_removal}; the results are
listed in Table~\ref{tab:closure_mA3_scan}. At every scanned cut all data-envelope members
describe the data, with goodness-of-fit $p$-values of at least 0.34. The single-sample
closure significance fluctuates between $0.5\sigma$ and $2.2\sigma$ without a monotonic
trend: it is lowest at 0.984 ($0.5\sigma$, $\hat{\mu}=+0.019$) and 0.988 ($0.9\sigma$,
$\hat{\mu}=+0.032$), and it is $1.6\sigma$ ($\hat{\mu}=+0.042$) at 0.992, the working point
adopted before the retraining.

The number of simulated background events in the $120<m_{\ell\ell\gamma\gamma}<130~\GeV$
window is small at these tight cuts, so the closure result itself fluctuates. This is
quantified with 40 bootstrap replicas of the simulation at the three best candidate cuts, in
which every simulated event is weighted by a Poisson(1) random number and the full closure
procedure is repeated (38--39 replicas per cut converge). The replicas give a closure
significance of $1.06\pm1.23\sigma$ at 0.984, $1.03\pm1.07\sigma$ at 0.988 and
$1.62\pm1.17\sigma$ at 0.992, with 23\%, 15\% and 42\% of the replicas above $2\sigma$, and
$\hat{\mu}=+0.035\pm0.066$, $+0.032\pm0.050$ and $+0.048\pm0.039$, respectively. The cut 0.988
gives the most stable closure, with the smallest spread of the significance and the smallest
fraction of replicas with a fake signal above $2\sigma$; it also has the best description of
the data (goodness-of-fit $p$-values of at least 0.92) and the lowest $R$ of the scan (1.16).
The working point of $m_{a}=3~\GeV$ is therefore set to 0.988. The expected limit at
$m_{a}=3~\GeV$ is 5\% weaker than at 0.992. The same limited simulation statistics also
explain why $R$ does not vary monotonically with the cut in Table~\ref{tab:closure_mA3_scan}."""),
(r"""    the data, and the background-only pseudo-data closure result (best-fit $\hat{\mu}$ in
    units of the signal yield at 0.988, and local significance). Entries marked ``fail''
    have no background function with a $p$-value above 0.01. The adopted working point is
    0.992.}""",
 r"""    the data, and the background-only pseudo-data closure result (best-fit $\hat{\mu}$ in
    units of the signal yield at 0.992, and local significance), obtained after the retraining
    of the BDT. The adopted working point is 0.988.}"""),
(r"""      0.980 & 1.17 & 29.7 & 719 & 0.02--0.06 & $+0.082$ & 2.3 \\
      0.984 & 1.22 & 29.2 & 571 & fail       & $+0.104$ & 3.3 \\
      0.986 & 1.35 & 28.6 & 483 & fail       & $+0.121$ & 4.1 \\
      0.988 & 1.21 & 29.3 & 394 & 0.47--0.78 & $+0.063$ & 2.3 \\
      0.990 & 1.38 & 28.6 & 302 & 0.36--0.62 & $+0.056$ & 2.2 \\
      \textbf{0.992} & \textbf{1.24} & \textbf{28.7} & \textbf{228} & \textbf{0.21--0.62} & $\mathbf{+0.018}$ & \textbf{0.8} \\
      0.994 & 1.32 & 27.7 & 155 & 0.28--0.66 & $+0.026$ & 1.4 \\
      0.996 & 1.77 & 22.1 &  87 & 0.18--0.55 & $+0.034$ & 2.1 \\ \hline""",
 r"""      0.980 & 1.25 & 30.3 & 703 & 0.76--0.91 & $+0.080$ & 1.8 \\
      0.984 & 1.18 & 31.1 & 550 & 0.45--0.64 & $+0.019$ & 0.5 \\
      0.986 & 1.21 & 30.8 & 463 & 0.76--0.88 & $+0.066$ & 1.7 \\
      \textbf{0.988} & \textbf{1.16} & \textbf{30.5} & \textbf{370} & \textbf{0.92--0.99} & $\mathbf{+0.032}$ & \textbf{0.9} \\
      0.990 & 1.38 & 30.2 & 286 & 0.83--0.95 & $+0.072$ & 2.2 \\
      0.992 & 1.36 & 30.2 & 206 & 0.51--0.69 & $+0.042$ & 1.6 \\
      0.994 & 1.44 & 30.2 & 145 & 0.66--0.71 & $+0.045$ & 2.1 \\
      0.996 & 1.57 & 30.2 &  75 & 0.34--0.52 & $+0.039$ & 2.2 \\ \hline"""),
(r"""    compatible with zero across the range; no mass point reaches a local significance of
    $2\sigma$ (the largest is $1.77\sigma$, at $m_{a}=28~\GeV$). At $m_{a}=5$ and $6~\GeV$ the
    signal-plus-background fit reaches the lower physical boundary of $\mu$ (a deficit in the
    pseudo-data) and $\hat{\mu}$ is not shown.}""",
 r"""    compatible with zero across the range; no mass point reaches a local significance of
    $2\sigma$ (the largest is $1.92\sigma$, at $m_{a}=28~\GeV$). At $m_{a}=5~\GeV$ the
    signal-plus-background fit reaches the lower physical boundary of $\mu$ (a deficit in the
    pseudo-data), so the $\hat{\mu}$ marker there sits at the boundary.}"""),
]
