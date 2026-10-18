#!/usr/bin/env python3
"""NPS-25-014 Sec 6: closure-test wording aligned with AN Sec 6.4 (sec:sculpt_closure, working points at m_a = 1-3):
pseudo-data = simulated Drell-Yan (with a prompt photon and with jets) selected at the same BDT threshold, which
reproduces the selected data yield; five thresholds defined by the number of selected data sideband events."""
PAS = "/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/"
fn = "Sec-06-00-BDTOptimisation.tex"
s = open(PAS + fn).read()
a = r"""background-only closure test: the background model is fitted to pseudo-data
built from the simulated background passing the BDT selection, normalized to
the data sidebands, and a threshold is accepted only if the fitted signal is
compatible with zero within one standard deviation and every function of the
data envelope has a goodness-of-fit probability above 0.01. Among the accepted
thresholds, scanned in steps of the number of selected data sideband events,
the one with the best expected sensitivity is adopted."""
b = r"""background-only closure test, in which the simulated Drell--Yan background
(with a prompt photon and with jets) selected with the same threshold, which
reproduces the selected data yield, is fitted as pseudo-data with exactly the
background-modeling and statistical procedure used for the data. Five
thresholds, defined by the number of selected data events in the \mllgg\
sidebands, are tested; a threshold is accepted only if the signal fitted to the
pseudo-data is compatible with zero within one standard deviation and every
function of the data envelope has a goodness-of-fit probability above 0.01, and
among the accepted thresholds the one with the best expected sensitivity is
adopted."""
assert s.count(a) == 1
open(PAS + fn, "w").write(s.replace(a, b)); print("edited", fn)
