#!/usr/bin/env python3
"""PAS NPS-25-014 for the low-mass BDT low_relE2_w123 (relE2 inputs, signal of m_a = 1, 2, 3 trained inside 115-133 GeV;
production since 2026-10-17). Same numbers as the AN edits (deploy_w123_20261017/an_numbers/, an_edits/):
  R 1.30 / 1.27 / 1.26; category gain 2.0 % at m_a = 1, at most 0.4 % elsewhere; m_a = 1 Lau1 bias -0.246 (the two
  alternatives were measured on the previous version, Lau1 -0.261); overlapping-shower population after the BDT 0.2 %;
  in-window 95.9-98.8 %; effSigma muon 1.33-1.76, electron 1.65-2.87 GeV over all m_a and eras; expected limits
  30.5 / 9.7 / 5.5 fb at m_a 1-3 (Run 2 17.8 / 4.7 / 2.8 fb); efficiency 12 / 36 / 46 % of the preselected signal.
Figures (limits, effSigma) are copied from the AN by pas_w123_figures.sh. Backups: <file>.bak_w123_20261017"""
import shutil

PAS = "/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/"
TAG = ".bak_w123_20261017"


def edit(fn, pairs):
    p = PAS + fn
    s = open(p).read()
    for a, b in pairs:
        assert s.count(a) == 1, (fn, s.count(a), a[:100])
        s = s.replace(a, b)
    shutil.copy2(p, p + TAG)
    open(p, "w").write(s)
    print("edited", fn, len(pairs))


edit("Sec-06-00-BDTOptimisation.tex", [
    ("""the training and a single-peaked signal at $\\mass{a}=1~\\GeV$; in the training of
this model the $\\mass{a}=1~\\GeV$ signal enters only within
$115<\\mllgg<133~\\GeV$.""",
     """the training and a single-peaked signal at $\\mass{a}=1~\\GeV$; in the training of
this model the signal of all three mass points enters only within
$115<\\mllgg<133~\\GeV$, so that the three points follow the same training rule."""),
    ("""$R=1.19$, 1.23, and 1.30 at these three mass points;""",
     """$R=1.30$, 1.27, and 1.26 at these three mass points;"""),
    ("""categories yields a gain in expected significance of at most 1.9\\% (at
$\\mass{a}=1~\\GeV$) and at most 0.5\\% elsewhere.""",
     """categories yields a gain in expected significance of at most 2.0\\% (at
$\\mass{a}=1~\\GeV$) and at most 0.4\\% elsewhere."""),
])

edit("Sec-07-00-BackgroundExtraction.tex", [
    ("""first-order Laurent function gives $-0.261$; the three functions of this
envelope are already at the lowest order, and removing the Laurent family or
allowing higher orders moves the excess to another function instead of removing
it, so the envelope is kept.""",
     """first-order Laurent function gives $-0.246$; the three functions of this
envelope are already at the lowest order, and removing the Laurent family or
allowing higher orders, tested with the previous version of the low-mass
classifier, moved the excess to another function instead of removing it, so the
envelope is kept."""),
])

edit("Sec-08-00-SignalExtraction.tex", [
    ("""low-mass BDT selection they amount to 0.3\\% of the selected signal, and between
94.5 and 97.2\\% of the simulated yield, depending on era and channel, falls""",
     """low-mass BDT selection they amount to 0.2\\% of the selected signal, and between
95.9 and 98.8\\% of the simulated yield, depending on era and channel, falls"""),
    ("""across the entire \\mass{a}\\ range, $\\sigma_{\\mathrm{eff}}$ ranging from about
1.4 to 1.8\\GeV\\ against 1.7 to 2.9\\GeV\\ depending on \\mass{a}\\ and era, with""",
     """across the entire \\mass{a}\\ range, $\\sigma_{\\mathrm{eff}}$ ranging from about
1.3 to 1.8\\GeV\\ against 1.7 to 2.9\\GeV\\ depending on \\mass{a}\\ and era, with"""),
])

edit("Sec-10-00-ResultAndInterpretation.tex", [
    ("""at $\\mass{a}=5~\\GeV$, and about 28\\unit{fb} at $\\mass{a}=1~\\GeV$, where the two""",
     """at $\\mass{a}=5~\\GeV$, and about 30\\unit{fb} at $\\mass{a}=1~\\GeV$, where the two"""),
    ("""weaker: 28.3, 9.5, and 5.7\\unit{fb} at $\\mass{a}=1$, 2, and 3\\GeV, compared with""",
     """weaker: 30.5, 9.7, and 5.5\\unit{fb} at $\\mass{a}=1$, 2, and 3\\GeV, compared with"""),
    ("""quantify: the low-mass classifier uses neither shower-shape variables nor photon
transverse momenta, and its thresholds are chosen with the background-only""",
     """quantify: the low-mass classifier uses neither shower-shape variables nor photon
transverse momenta, its signal is trained only inside the window
$115<\\mllgg<133~\\GeV$, and its thresholds are chosen with the background-only"""),
    ("""11, 34, and 45\\% of the preselected signal at $\\mass{a}=1$, 2, and 3\\GeV,""",
     """12, 36, and 46\\% of the preselected signal at $\\mass{a}=1$, 2, and 3\\GeV,"""),
])

edit("Sec-11-00-Summary.tex", [
    ("""$2\\leq\\mass{a}\\leq30~\\GeV$ and about 28\\unit{fb} at $\\mass{a}=1~\\GeV$. Above""",
     """$2\\leq\\mass{a}\\leq30~\\GeV$ and about 30\\unit{fb} at $\\mass{a}=1~\\GeV$. Above"""),
])

edit("NPS-25-014.tex", [
    ("""reach about 28~fb at $m_{a} = 1$ GeV""", """reach about 30~fb at $m_{a} = 1$ GeV"""),
])
