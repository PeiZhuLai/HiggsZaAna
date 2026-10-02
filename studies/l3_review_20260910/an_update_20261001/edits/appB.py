PAIRS = [
(r"""This appendix collects the per-era trigger-efficiency plots used in
Section~\ref{sec:datamc}. The individual Run 3 eras trigger-efficiency plots are shown. The figures below cover the full set of ALP mass hypotheses used in the analysis.

\newcommand{\trigpanel}[4]{%
  \subfigure[#1]{\includegraphics[width=0.155\textwidth]{figure_ALP/run3/HiggsZaAna/Plot/plots/trigEffCompareVlepPt/#2/mA_M#3/trigeffCompare_#4_mA_M#3_#2.pdf}}%
}

\newcommand{\trigmassfig}[1]{%""",
 r"""This appendix collects the per-era trigger-efficiency plots used in
Section~\ref{sec:datamc}. The individual Run 3 eras trigger-efficiency plots are shown. The figures below cover the full set of ALP mass hypotheses used in the analysis.

For $m_{a}=1~\GeV$ (Figure~\ref{fig:trigger_eff_pt_m1}) the efficiency of each leg is
measured, as in Figure~\ref{fig:triggersummary_run3_2024}, in events with at least two
selected same-flavor leptons in which the other lepton passes its own offline double-lepton
threshold (15 (10)~\GeV for the subleading electron (muon) when the leading leg is
measured, 25 (20)~\GeV for the leading electron (muon) when the subleading leg is
measured). On the plateau ($\pt>40~\GeV$ for the leading and $\pt>30~\GeV$ for the
subleading lepton) the leading-leg efficiency is 81--86\% (93--94\%) for the double-electron
(double-muon) triggers and 95--96\% (99\%) for the logical OR of single- and double-lepton
triggers across the five eras; for the subleading leg it is 83--88\% (93--94\%) and
96--98\% (99\%), respectively. This requirement on the other leg is not applied in the
figures for $m_{a}\geq2~\GeV$, which are kept for comparison: there the double-lepton
efficiency of the leading leg is diluted by events in which the subleading lepton cannot
satisfy the double-lepton trigger (64--69\% for electrons and 78--79\% for muons on the
plateau at $m_{a}=1~\GeV$ with this definition), while the subleading-leg curves are
unaffected. The trigger efficiencies are shown for illustration only; the analysis uses
the trigger scale factors of Appendix~\ref{app:mccorrections}, which do not depend on
this choice.

\newcommand{\trigpanel}[4]{%
  \subfigure[#1]{\includegraphics[width=0.155\textwidth]{figure_ALP/run3/HiggsZaAna/Plot/plots/trigEffCompareVlepPt/#2/mA_M#3/trigeffCompare_#4_mA_M#3_#2.pdf}}%
}

\newcommand{\trigpanelOL}[4]{%
  \subfigure[#1]{\includegraphics[width=0.155\textwidth]{figure_ALP/run3/HiggsZaAna/Plot/plots/trigEffCompareVlepPt_trigdenomV2/denomOL/#2/mA_M#3/trigeffCompare_#4_mA_M#3_#2.pdf}}%
}

\begin{figure}[p]
  \centering
  \trigpanelOL{22pre, e1}{2022preEE}{1}{Electron_lead}
  \trigpanelOL{22pre, e2}{2022preEE}{1}{Electron_sublead}
  \trigpanelOL{22pre, $\mu_{1}$}{2022preEE}{1}{Muon_lead}
  \trigpanelOL{22pre, $\mu_{2}$}{2022preEE}{1}{Muon_sublead}
  \trigpanelOL{22post, e1}{2022postEE}{1}{Electron_lead}
  \trigpanelOL{22post, e2}{2022postEE}{1}{Electron_sublead}\\
  \trigpanelOL{22post, $\mu_{1}$}{2022postEE}{1}{Muon_lead}
  \trigpanelOL{22post, $\mu_{2}$}{2022postEE}{1}{Muon_sublead}
  \trigpanelOL{23pre, e1}{2023preBPix}{1}{Electron_lead}
  \trigpanelOL{23pre, e2}{2023preBPix}{1}{Electron_sublead}
  \trigpanelOL{23pre, $\mu_{1}$}{2023preBPix}{1}{Muon_lead}
  \trigpanelOL{23pre, $\mu_{2}$}{2023preBPix}{1}{Muon_sublead}\\
  \trigpanelOL{23post, e1}{2023postBPix}{1}{Electron_lead}
  \trigpanelOL{23post, e2}{2023postBPix}{1}{Electron_sublead}
  \trigpanelOL{23post, $\mu_{1}$}{2023postBPix}{1}{Muon_lead}
  \trigpanelOL{23post, $\mu_{2}$}{2023postBPix}{1}{Muon_sublead}
  \trigpanelOL{24, e1}{2024}{1}{Electron_lead}
  \trigpanelOL{24, e2}{2024}{1}{Electron_sublead}\\
  \trigpanelOL{24, $\mu_{1}$}{2024}{1}{Muon_lead}
  \trigpanelOL{24, $\mu_{2}$}{2024}{1}{Muon_sublead}
  \caption{Trigger efficiency as a function of lepton $p_T$ for the Run 3
  signal sample with $m_{a}=1~\GeV$. The five eras are shown separately for
  the leading and subleading electrons and muons. The efficiency of each leg is
  measured in events with at least two selected same-flavor leptons in which the other
  lepton passes its offline double-lepton threshold.}
  \label{fig:trigger_eff_pt_m1}
\end{figure}

\newcommand{\trigmassfig}[1]{%"""),
(r"""  \caption{Trigger efficiency as a function of lepton $p_T$ for the Run 3
  signal sample with $m_{a}=#1~\GeV$. The five eras are shown separately for
  the leading and subleading electrons and muons.}""",
 r"""  \caption{Trigger efficiency as a function of lepton $p_T$ for the Run 3
  signal sample with $m_{a}=#1~\GeV$. The five eras are shown separately for
  the leading and subleading electrons and muons. Here the other lepton is not required
  to pass its double-lepton threshold, which dilutes the double-lepton efficiency of the
  leading leg (see text).}"""),
(r"""\trigmassfig{1}
\trigmassfig{2}""", r"""\trigmassfig{2}"""),
]
