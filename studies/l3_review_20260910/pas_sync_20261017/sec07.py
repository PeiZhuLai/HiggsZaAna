#!/usr/bin/env python3
"""NPS-25-014 Sec 7 bias paragraph synced to AN-25-172 Sec 7.4/7.5: 108 truth functions, 107 within 0.2;
exception Lau1 at m_a = 1 (-0.261, envelope already lowest order, removing Laurent -> Exp1 -0.267, higher orders
-> Pow3 -0.317), negative pull = expected limit possibly too strong by up to 0.26 of the signal-strength
uncertainty; elsewhere largest -0.182 (m_a = 10, Exp3) -> 2% inflation; Bernstein excluded at 1, 2, 4-8, 21, 27, 28;
Bernstein order capped at 3 (Bern1), 11 (Bern2, lowest orders), 20 (Bern5), 25 (Bern5); Laurent removed at 24, 25."""
PAS = "/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/"
fn = "Sec-07-00-BackgroundExtraction.tex"
s = open(PAS + fn).read()
a = r"""The mean of the resulting pull distribution stays within 0.2 for all 110 truth
functions of the scan, the largest value being $-0.192$, which inflates the
total uncertainty on the signal strength by only 2\%. Reaching this required
excluding the Bernstein family at ten mass points, $\mass{a}=1$, 2, 4--8, 21,
27, and 28\GeV, and restricting the scan to the lowest orders at
$\mass{a}=11~\GeV$: where the data cannot constrain the free parameters of a
flexible function, that function absorbs part of an injected signal during the
discrete profiling and biases the fitted signal strength low."""
b = r"""The mean of the resulting pull distribution stays within 0.2 for 107 of the
108 truth functions of the scan. Reaching this required excluding the Bernstein
family at ten mass points, $\mass{a}=1$, 2, 4--8, 21, 27, and 28\GeV, and the
Laurent family at $\mass{a}=24$ and 25\GeV, limiting the Bernstein order at
$\mass{a}=3$, 20, and 25\GeV, and restricting the scan to the lowest orders at
$\mass{a}=11~\GeV$: where the data cannot constrain the free parameters of a
flexible function, that function absorbs part of an injected signal during the
discrete profiling and biases the fitted signal strength low. Apart from
$\mass{a}=1~\GeV$, the largest mean pull is $-0.182$, which inflates the total
uncertainty on the signal strength by only 2\%. At $\mass{a}=1~\GeV$ the
first-order Laurent function gives $-0.261$; the three functions of this
envelope are already at the lowest order, and removing the Laurent family or
allowing higher orders moves the excess to another function instead of removing
it, so the envelope is kept. Since the pull is negative, the expected limit at
this mass point may be too strong by up to about a quarter of the uncertainty
on the signal strength."""
assert s.count(a) == 1
open(PAS + fn, "w").write(s.replace(a, b)); print("edited", fn)
