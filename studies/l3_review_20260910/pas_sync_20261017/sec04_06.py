#!/usr/bin/env python3
"""NPS-25-014 Sec 4 and Sec 6 synced to AN-25-172 (current master).
Sec 4: FSR photon separation from every signal photon (AN Sec 4.4), not only the leading one.
Sec 6: low-mass inputs = low_swapPt_relE2 (AN Sec 6.1, Table tab:bdt_input_lowmass), m_a = 1 signal window,
sideband reweighting also applied to the signal (AN Sec 4.5), working points from the background-only closure
test (AN Sec 6.4: x0.5 data-sideband rule), R = 1.19/1.23/1.30, category gain 1.9% at m_a = 1 (AN Sec 6.5).
Backups: *.bak_ANsync_20261017"""
PAS = "/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/"

def edit(fn, pairs):
    s = open(PAS + fn).read()
    for a, b in pairs:
        assert s.count(a) == 1, (fn, a[:80])
        s = s.replace(a, b)
    open(PAS + fn, "w").write(s)
    print("edited", fn)

edit("Sec-04-00-ObjectDefinition.tex", [
    (r"""excluded from the muon isolation, and is required to be separated by
$\Delta R>0.2$ from the leading photon entering the H boson candidate.""",
     r"""excluded from the muon isolation, and is required to be separated by
$\Delta R>0.2$ from both photons entering the H boson candidate, so that the
subleading ALP photon cannot be reclassified as FSR and added to a muon."""),
])

edit("Sec-06-00-BDTOptimisation.tex", [
    (r"""mass spectrum. A dedicated low-mass model is therefore trained for
$\mass{a}=1,2,3~\GeV$ using a reduced set of five robust inputs:
$\pt(\gamma_{1})/\mllgg$, $R_{9}(\gamma_{1})$, the diphoton ALP isolation
$I^{\gamma,a}$, $\Delta R(\gamma_{1},\gamma_{2})$, and the mass-hypothesis
variable. Both models share the same mass-hypothesis parameterization and
training configuration.""",
     r"""mass spectrum. A dedicated low-mass model is therefore trained for
$\mass{a}=1,2,3~\GeV$ using a reduced set of five inputs: $R_{9}(\gamma_{1})$,
the diphoton ALP isolation $I^{\gamma,a}$, $\Delta R(\gamma_{1},\gamma_{2})$,
the relative energy uncertainty of the subleading photon
$\sigma_{E}/E(\gamma_{2})$, and the mass-hypothesis variable. This set was
selected in a dedicated search over input sets, requiring no excess of the
simulated background in the H boson peak window on a test sample not used in
the training and a single-peaked signal at $\mass{a}=1~\GeV$; in the training of
this model the $\mass{a}=1~\GeV$ signal enters only within
$115<\mllgg<133~\GeV$. Both models share the same mass-hypothesis
parameterization and training configuration."""),
    (r"""($95<\mllgg<115$ and $135<\mllgg<180~\GeV$) and applied to the background
template used for the training, so that the BDT input shapes match the data. The""",
     r"""($95<\mllgg<115$ and $135<\mllgg<180~\GeV$) and applied to the simulated
background and, with the same correction factors, to the simulated signal, so
that the BDT input shapes match the data; for the signal the weights are
normalized so that the yield after the preselection is unchanged. The"""),
    (r"""stable. The selected thresholds are generally above 0.9. For
$\mass{a}=1$, 2, and 3\GeV\ the threshold is instead chosen so that the residual
sculpting of the \mllgg\ spectrum induced by the BDT selection, quantified by
the ratio $R(\mass{a})$ of the fraction of background events in the peak region
after and before the selection, stays as close to unity as possible; this yields
$R=0.99$, 1.10, and 1.04 at these three mass points, at a small cost in raw
significance. For the \mass{a}\ hypotheses that are not directly simulated, the
threshold of the closest generated mass point is used. A single category per
\mass{a}\ is used, since increasing the number of categories yields a gain in
expected significance of less than 2\%.""",
     r"""stable. The selected thresholds are generally above 0.9. For
$\mass{a}=1$, 2, and 3\GeV\ the threshold is instead chosen with a
background-only closure test: the background model is fitted to pseudo-data
built from the simulated background passing the BDT selection, normalized to
the data sidebands, and a threshold is accepted only if the fitted signal is
compatible with zero within one standard deviation and every function of the
data envelope has a goodness-of-fit probability above 0.01. Among the accepted
thresholds, scanned in steps of the number of selected data sideband events,
the one with the best expected sensitivity is adopted. The residual sculpting of
the \mllgg\ spectrum, quantified by the ratio $R(\mass{a})$ of the fraction of
background events in the peak region after and before the selection, is
$R=1.19$, 1.23, and 1.30 at these three mass points; the closure test shows that
it does not produce a spurious signal. For the \mass{a}\ hypotheses that are
not directly simulated, the threshold of the closest generated mass point is
used. A single category per \mass{a}\ is used, since increasing the number of
categories yields a gain in expected significance of at most 1.9\% (at
$\mass{a}=1~\GeV$) and at most 0.5\% elsewhere."""),
])
