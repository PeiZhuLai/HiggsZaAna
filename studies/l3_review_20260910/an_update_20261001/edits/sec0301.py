OL = "figure_ALP/run3/HiggsZaAna/Plot/plots/trigEffCompareVlepPt_trigdenomV2/denomOL/2024/mA_M1/"
OLD = "figure_ALP/run3/HiggsZaAna/Plot/plots/trigEffCompareVlepPt/2024/mA_M1/"
PAIRS = [
(r"""lead and sublead lepton \pt in the electron and muon channels. The full
per-era set of trigger-efficiency plots for the same mass point is collected in
Appendix~\ref{app:trigger-eff-pt}. Simultaneous use of single- and double-lepton triggers allows the analysis to achieve higher trigger efficiency-- the electron (muon) channel achieves a trigger efficiency of about 90\% (95\%) for the 2024 signal sample at $m_{a}=1~\GeV$ after applying the logical OR of single- and double-lepton triggers, compared to about 65\% (80\%) if only single-lepton triggers were used.""",
 r"""lead and sublead lepton \pt in the electron and muon channels. The efficiency of
each leg is measured in signal events with at least two selected same-flavor leptons in
which the other lepton passes its own offline double-lepton threshold: 15 (10)~\GeV for
the subleading electron (muon) when the leading leg is measured, and 25 (20)~\GeV for the
leading electron (muon) when the subleading leg is measured. Without this requirement on
the other leg, the double-lepton efficiency of the leading leg is diluted by events in which
the subleading lepton cannot satisfy the double-lepton trigger; on the plateau of the
leading leg it would then be 64--69\% (78--79\%) for electrons (muons) instead of
81--86\% (93--94\%) across the five eras, while the subleading-leg efficiencies are
unchanged. The full per-era set of trigger-efficiency plots for the same mass point is collected in
Appendix~\ref{app:trigger-eff-pt}. Simultaneous use of single- and double-lepton triggers allows the analysis to achieve higher trigger efficiency-- on the plateau of the leading lepton ($\pt>40~\GeV$), the electron (muon) channel achieves a trigger efficiency of about 95\% (99\%) for the 2024 signal sample at $m_{a}=1~\GeV$ after applying the logical OR of single- and double-lepton triggers, compared to about 82\% (94\%) if only the double-lepton triggers were used."""),
] + [(OLD + f"trigeffCompare_{x}_mA_M1_2024}}", OL + f"trigeffCompare_{x}_mA_M1_2024}}") for x in ("Electron_lead", "Electron_sublead", "Muon_lead", "Muon_sublead")] + [
(r"""shown separately for the leading and subleading leptons in the electron and muon channels. The dots with error bars represent the trigger efficiency, as given by the left y-axis:""",
 r"""shown separately for the leading and subleading leptons in the electron and muon channels. The efficiency of each leg is measured in events with at least two selected same-flavor leptons in which the other lepton passes its offline double-lepton threshold. The dots with error bars represent the trigger efficiency, as given by the left y-axis:"""),
]
