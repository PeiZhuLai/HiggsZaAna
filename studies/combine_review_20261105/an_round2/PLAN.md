# NPS-25-014 Combine review: remaining items

Representative mass points: m_a = 1, 10, 20, 30 GeV (the points of the AN impact figures).
Everything is blind: Asimov datasets (-t -1) or toys generated from the background model; no observed data.
Cards and workspaces: current production (Combine/root_t2w_{xsec,wilson}, output_combine_results_{xsec,wilson}).
The twiki is not touched; reviewer questions are not answered through it.

## 1. Statistical checks

| Step | What | Where | Jobs |
|---|---|---|---|
| 1a | Asimov b-only dataset per mass (input of HybridNew) | local | 4 |
| 1b | Hesse covariance (correlation matrix): s+b and b-only fits on the r=0 Asimov, both cards | local | 8 |
| 1c | 1D likelihood scan of r on the r=0 Asimov, xsec card, r in [0, 3 x q50] | local | 4 |
| 1d | stat/syst split on the r=q50 Asimov: full fit and fit with all constrained nuisances frozen, both cards | local | 16 fits |
| 1e | HybridNew expected CLs limits (LHC-limits, 16 r points x 6 jobs x 500 toys), xsec card | condor | 384 |
| 1f | Bias with injection 0, 2 and 3 x q50 (2 chunks x 1000 toys per truth function), production bias fit options | condor | 24 |
| 1g | B-only (r=0 Asimov) impacts, fixed function, production ranges, both cards | condor | 8 |

## 2. Plots and tables

Correlation matrices, likelihood scans, HybridNew/asymptotic table, stat/syst table, B-only impact figures,
bias-vs-injection table, signal shape variations (nominal and +-1 sigma of the shape nuisances),
pre-/post-fit distributions on the Asimov dataset with the background uncertainty band.

## 3. AN

New subsection with the checks above (Sec. 10 or an appendix), compile with compile.sh, commit and push
(only the files the tex uses; no observed limits).

## 4. Datacard repository

Prepare input/mA10_xsec.txt, wilson/, shapes (LFS), input/systematics.yml, README, .gitlab-ci.yml with the
current q50 of mA10; run the CI commands locally in the container. Pushing to the review GitLab repository is
visible to the reviewers: done after the user's OK.

## 5. Bookkeeping

Update doc/HZa/review_comments/Combine_review/combine_review_summary.md and combine_review_checklist.md
(docdir: commit and push). Code commits of flashggFinalFit and HiggsZaAna: ask first.
