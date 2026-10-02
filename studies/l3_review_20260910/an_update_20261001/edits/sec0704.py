PAIRS = [
(r"""and at $m_{a}=20~\GeV$ the Bernstein order is capped at five, as motivated by the bias study (see text).""",
 r"""at $m_{a}=20$ and $25~\GeV$ the Bernstein order is capped at five, and the Laurent family is removed at $m_{a}=24$ and $25~\GeV$, as motivated by the bias study (see text)."""),
(r"""      & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp3  & Exp1  & Exp1  & Exp3  & Exp1  & Exp1  & Exp1  \\
      & Pow1  & Pow3  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  \\
      & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  \\
      &       &       & Bern3 &       &       &       &       &       & Bern6 & Bern5 & Bern2 & Bern5 & Bern4 & Bern5 & Bern5 \\\hline""",
 r"""      & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp3  & Exp3  & Exp1  & Exp3  & Exp1  & Exp1  & Exp3  \\
      & Pow1  & Pow3  & Pow3  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  \\
      & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  \\
      &       &       & Bern2 &       &       &       &       &       & Bern6 & Bern5 & Bern2 & Bern5 & Bern4 & Bern5 & Bern5 \\\hline"""),
(r"""      & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp3  & Exp3  & Exp3  & Exp3  & Exp3  & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  \\
      & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  \\
      & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  \\
      & Bern6 & Bern5 & Bern5 & Bern5 & Bern5 &       & Bern6 & Bern6 & Bern5 & Bern6 & Bern6 &       &       & Bern6 & Bern6 \\ \hline""",
 r"""      & Exp1  & Exp3  & Exp3  & Exp1  & Exp3  & Exp1  & Exp3  & Exp3  & Exp1  & Exp1  & Exp1  & Exp1  & Exp1  & Exp3  & Exp1  \\
      & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  & Pow1  \\
      & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  &       &       & Lau1  & Lau1  & Lau1  & Lau1  & Lau1  \\
      & Bern6 & Bern5 & Bern5 & Bern5 & Bern5 &       & Bern6 & Bern6 & Bern5 & Bern5 & Bern6 &       &       & Bern5 & Bern6 \\ \hline"""),
(r"""Of the 110 fits that enter the envelopes of the 30 mass points, 109 converged
with status 0; the remaining one (Exp1 at $m_{a}=26~\GeV$) returned status 1 and
has a goodness-of-fit probability of 0.31. Fifteen of the 110 functions, at
$m_{a}=2$, 24, 25, 26, 27, 28 and 29~GeV, are retained by the fallback rule
described above, i.e.\ they did not pass the combined goodness-of-fit and F-test
selection (all of them have a goodness-of-fit probability below 0.01) but are
kept so that every family remains represented in the envelope; these are flagged
in the tables. At $m_{a}=2$, 24, 25 and 28~GeV this applies to all envelope
members except, at 24 and 25~GeV, the exponential member (Exp3,
with probabilities of 0.076 and 0.068); at $m_{a}=26$, 27 and 29~GeV only the
Laurent member is affected. At $m_{a}=2~\GeV$ the low goodness of fit was
traced to the residual sculpting discussed in Section~\ref{sec:sculpt_R}. The
origin of the poor goodness of fit in the $24 \le m_{a} \le 29~\GeV$ range,
which appears with the final-state-radiation-corrected data, has not yet been
established. At all of these mass points the signal-injection bias study,
performed with each fallback function in turn as the truth model, gives a mean
pull within the 0.2 threshold (the largest being $-0.147$ for Bern5 at
$m_{a}=24~\GeV$), so the fallback members do not degrade the signal-strength
measurement. The bias values are discussed in Section~\ref{sec:bias_studies}.""",
r"""All 108 fits that enter the envelopes of the 30 mass points converged with
status 0, and the goodness-of-fit probabilities of the envelope members lie
between 0.06 and 0.99. The poor goodness of fit found in an earlier iteration in
the $24 \le m_{a} \le 29~\GeV$ range is no longer present after the parameter
ranges of the fit functions were widened there and new starting values were used
for the exponential family at $m_{a}=28$ and $29~\GeV$. A single function, Exp3 at
$m_{a}=29~\GeV$, is retained by the fallback rule described above, i.e.\ it did not
pass the combined goodness-of-fit and F-test selection but is kept so that every
family remains represented in the envelope; it is flagged in the tables. Its
goodness-of-fit probability is 0.98 and the signal-injection bias study performed
with it as the truth model gives a mean pull of $+0.006$, so it does not degrade
the signal-strength measurement. The bias values are discussed in
Section~\ref{sec:bias_studies}."""),
(r"""The envelope composition described above is the outcome of the automatic
selection at all but five mass points, where the bias study of""",
 r"""The envelope composition described above is the outcome of the automatic
selection at all but seven mass points, where the bias study of"""),
(r"""for a sample of about one hundred data events (111 with the current selection),
and the higher-order""",
 r"""for a sample of about one hundred data events,
and the higher-order"""),
(r"""Bernstein order is further limited to two: with Bern2 the largest bias of the
mass point is $-0.139$ and the goodness-of-fit probabilities of the four members
lie between 0.20 and 0.70. At $m_{a}=20~\GeV$ the sixth-order Bernstein member was
the only function above the threshold ($+0.208$); following the same rule, the
next lower order is used, and with Bern5 the largest bias of the mass point is
$-0.097$ (goodness-of-fit probabilities 0.28--0.49). At $m_{a}=2$, 4 and 27~\GeV the""",
 r"""Bernstein order is further limited to two: with Bern2 the largest bias of the
mass point is $-0.153$ and the goodness-of-fit probabilities of the four members
lie between 0.18 and 0.88. At $m_{a}=20~\GeV$ the sixth-order Bernstein member was
the only function above the threshold ($+0.208$); following the same rule, the
next lower order is used, and with Bern5 the largest bias of the mass point is
$-0.154$ (goodness-of-fit probabilities 0.49--0.95). After the retraining of the BDT
that followed the re-production of the 2024 DY+jets samples, three functions were
found above the threshold: Lau1 at $m_{a}=24~\GeV$ ($-0.232$), and Bern6 ($+0.212$)
and Lau1 ($-0.215$) at $m_{a}=25~\GeV$. Following the same rule, the Bernstein order
at $m_{a}=25~\GeV$ is capped at five, and the Laurent family, whose lowest order was
already in use, is removed at both mass points. With these changes the largest bias is
$-0.129$ (Exp1) at $m_{a}=24~\GeV$ and $-0.113$ (Exp1) at $m_{a}=25~\GeV$, with
goodness-of-fit probabilities of 0.58--0.89 and 0.52--0.88; the function preferred by
the penalized likelihood, and hence the expected limit, is unchanged at both mass
points. At $m_{a}=2$, 4 and 27~\GeV the"""),
]
