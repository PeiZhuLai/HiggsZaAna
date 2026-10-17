#!/usr/bin/env python3
"""AN-25-172 Sec 6.3: low-mass K-S p-values from train/test subsets built as the training sample
(plot_overtrain_lowmass.py after the fix: m_a = 1 in 115-133 GeV, signal factor x sideband reweight) -> 0.625 / 0.304
(log: logs/overtrain_lowmass_KSrecipe.log; equal to the trainer self-check in logs/retrain_check.log).
AUC sentence: state the definition used by plot_mva_diagnostics.py (nominal weights, m_a = 1 not windowed).
Backup: Sec-06-03-BDTevaluation.tex.bak_KSrecipe"""
import shutil
AN = "/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172/"
fn = "Sec-06-03-BDTevaluation.tex"
shutil.copy2(AN + fn, AN + fn + ".bak_KSrecipe")
s = open(AN + fn).read()
pairs = [
    (r"""The high-mass model reaches a test-sample area under the ROC curve (AUC) of $0.983$, while the low-mass model reaches $0.940$;""",
     r"""The high-mass model reaches a test-sample area under the ROC curve (AUC) of $0.983$, while the low-mass model reaches $0.940$ (both evaluated with the nominal simulation weights, without the sideband reweighting, and for the low-mass model with the $m_{a}=1~\GeV$ signal not restricted to the training window);"""),
    (r"""The agreement between the training and test subsets is quantified with the Kolmogorov--Smirnov (K--S) test. The resulting p-values are $0.46$ (signal) and $0.49$ (background) for the high-mass model, and $0.35$ (signal) and $0.30$ (background) for the low-mass model.""",
     r"""The agreement between the training and test subsets is quantified with the Kolmogorov--Smirnov (K--S) test, with both subsets built exactly as the training sample of the corresponding model: the event weights include the sideband reweighting of Section~\ref{ssec:mcweights}, and for the low-mass model the $m_{a}=1~\GeV$ signal enters only within $115<m_{\ell\ell\gamma\gamma}<133~\GeV$. The resulting p-values are $0.46$ (signal) and $0.49$ (background) for the high-mass model, and $0.62$ (signal) and $0.30$ (background) for the low-mass model."""),
]
for a, b in pairs:
    assert s.count(a) == 1, a[:90]
    s = s.replace(a, b)
open(AN + fn, "w").write(s)
print("Sec 6.3 edited")
