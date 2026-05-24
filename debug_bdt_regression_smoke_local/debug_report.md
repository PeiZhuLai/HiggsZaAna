# BDT Regression Debug Report

## Executive Summary

The strongest code-level regression candidates are ranked below.  When event-level inputs were available, the runtime findings are appended afterwards.

### [P0] Event-weight definition changed

Evidence:
- Old training weight expression: `factor*pho1SFs*pho2SFs`.
- New training weight expression: `factor`.
- Parque2Root_BDT.py currently sets `factor = weight_central`, while hzgml/skim_ntuples.py defines `weight_corr = weight_central × SFs × PU × trigger × ...`.
- run3_Za_BDT.py also converts negative training weights to abs(weight) and rebalances classes before fitting.
Recommendation: Dump per-sample weight stats and compare factor, factor_nominal, and weight_corr-style compositions before retraining.

### [P0] Observed performance regression is concentrated in low and high mass regions

Evidence:
- mA=6: threshold 0.990->0.985, eff ratio=0.870, bkg ratio=10.456
- mA=4: threshold 0.980->0.975, eff ratio=0.800, bkg ratio=4.659
- mA=30: threshold 0.980->0.980, eff ratio=0.746, bkg ratio=4.388
- mA=20: threshold 0.990->0.985, eff ratio=0.923, bkg ratio=4.345
- mA=1: threshold 0.955->0.920, eff ratio=0.683, bkg ratio=3.700
- mA=3: threshold 0.985->0.975, eff ratio=0.775, bkg ratio=3.554
- mA=7: threshold 0.985->0.990, eff ratio=0.750, bkg ratio=3.551
- mA=5: threshold 0.985->0.985, eff ratio=0.766, bkg ratio=3.502
- mA=2: threshold 0.980->0.965, eff ratio=0.764, bkg ratio=3.430
- mA=25: threshold 0.985->0.985, eff ratio=0.759, bkg ratio=3.107
Recommendation: Prioritize diagnostics for mA=1,2,3,20,25,30 first.

### [P0] Selection and control-region definition changed substantially

Evidence:
- Old tree2array selection is hardcoded to `passChaHadIso && passNeuHadIso && passdR_gl && passHOverE && H_m>110 && H_m<180`.
- New convert() mode is `forwarded_selection_argument`.
- New convert() selection is `None`.
- New background selection is `H_m>95 && H_m<180`.
- New signal selection is `H_m>115 && H_m<135`.
- The old script applies passChaHadIso/passNeuHadIso/passdR_gl/passHOverE inside convert(); the new script now forwards externally supplied selections instead.
Recommendation: Validate preselection yields before and after each selection component, especially in low-mass and high-mass regions.

### [P1] Background composition changed

Evidence:
- Old training background samples: ['DYJets'].
- New top-level background sample: ['All_Bkg'], built from ['DYJetsToLL', 'DYGto2LG'].
- Signal preparation years in run3 inputs: ['2022preEE', '2022postEE', '2023preBPix', '2023postBPix', '2024'].
Recommendation: Compare preselection yields separately for DYJets and DYGto2LG-like components, especially in mA=1,2,3 and mA=20,25,30.

### [P1] Current model selection is no longer pure Optuna best-trial

Evidence:
- KS_PVALUE_TARGET = 0.05.
- run3_Za_BDT.py defines `regularized_xgb_candidates()` and selects the first candidate passing a weighted-KS target.
- Legacy pickle hyperparameters: {'path': '/Users/laipeizhu/Documents/2026 HEP/HZa/HiggsZaAna/HZaMVA/keep/model_Za_BDT_run3.pkl', 'exists': True, 'objective': 'binary:logistic', 'eval_metric': 'logloss', 'num_trees': 700, 'max_depth': 3, 'learning_rate': 0.069853574, 'min_child_weight': 35, 'gamma': 0.257119179, 'reg_alpha': 3.36629653, 'reg_lambda': 5.62203979, 'subsample': 1, 'colsample_bytree': 1, 'scale_pos_weight': 2.2622757, 'best_score': 0.12962762791870672, 'best_iteration': 699, 'num_feature': 15}.
- Current pickle hyperparameters: {'path': '/Users/laipeizhu/Documents/2026 HEP/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_run3.pkl', 'exists': True, 'objective': 'binary:logistic', 'eval_metric': 'logloss', 'num_trees': 700, 'max_depth': 2, 'learning_rate': 0.0399999991, 'min_child_weight': 60, 'gamma': 2, 'reg_alpha': 1, 'reg_lambda': 30, 'subsample': 0.699999988, 'colsample_bytree': 0.699999988, 'scale_pos_weight': 1, 'best_score': 0.13367451742438063, 'best_iteration': 699, 'num_feature': 15}.
Recommendation: Check whether the KS-safe regularized model shifted the score tail enough to move thresholds down at low mA.

### [P1] Sideband reweighting is injected directly into current background training weights

Evidence:
- run3_Za_BDT.py loads sideband_run3_iterative.json and replaces `factor` with sideband-reweighted `factor` for background MC.
- The JSON performs 5 iterative reweighting rounds over all 15 BDT variables, including `param`.
- The reweight mass-hypothesis definition uses event-hash assignment over the 14 generated mass points, not the 1-30 scan.
Recommendation: Compare score tails and feature distributions with sideband reweight disabled versus enabled.

### [P1] There is no application split in the current run3 signal pipeline

Evidence:
- Parque2Root_BDT.py writes only `train` and `test` trees when `--split` is used.
- Background and data are stored only as `inclusive`.
- run3_Za_BDT.py then performs an additional 50/50 train_test_split inside Python.
Recommendation: Explicitly report tree-level and in-script splits before comparing any threshold optimization to older studies.

### [P1] Training mass hypotheses do not match the application scan

Evidence:
- Current background param assignment uses masses [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30].
- Application computes scores for masses [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30].
- run3_Za_BDT_v1.py previously assigned background masses from 1-30 for the parameterized training.
Recommendation: Quantify how much mA interpolation/extrapolation degrades score ranking for 1/2/3 and 20/25/30.

## Summary Comparison

| mA | old_threshold | new_threshold | threshold_diff | old_signal_eff | new_signal_eff | signal_eff_ratio | old_bkg_yield | new_bkg_yield | bkg_yield_ratio | warning_flag |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.955 | 0.92 | -0.03499999999999992 | 0.492 | 0.336 | 0.6829268292682927 | 83.195 | 307.845 | 3.7002824688983718 | background_yield_gt_2x_reference;signal_eff_lt_0p8_reference |
| 2 | 0.98 | 0.965 | -0.015000000000000013 | 0.669 | 0.511 | 0.7638266068759342 | 26.373 | 90.45 | 3.4296439540439083 | background_yield_gt_2x_reference;signal_eff_lt_0p8_reference |
| 3 | 0.985 | 0.975 | -0.010000000000000009 | 0.763 | 0.591 | 0.7745740498034076 | 7.888 | 28.033 | 3.553879310344828 | background_yield_gt_2x_reference;signal_eff_lt_0p8_reference |
| 4 | 0.98 | 0.975 | -0.0050000000000000044 | 0.842 | 0.674 | 0.8004750593824229 | 5.062 | 23.586 | 4.65942315290399 | background_yield_gt_2x_reference |
| 5 | 0.985 | 0.985 | 0.0 | 0.849 | 0.65 | 0.7656065959952886 | 5.123 | 17.943 | 3.502439976576225 | background_yield_gt_2x_reference;signal_eff_lt_0p8_reference |
| 6 | 0.99 | 0.985 | -0.0050000000000000044 | 0.817 | 0.711 | 0.8702570379436965 | 2.494 | 26.077 | 10.45589414595028 | background_yield_gt_2x_reference |
| 7 | 0.985 | 0.99 | 0.0050000000000000044 | 0.86 | 0.645 | 0.75 | 5.264 | 18.693 | 3.5511018237082066 | background_yield_gt_2x_reference;signal_eff_lt_0p8_reference |
| 8 | 0.99 | 0.99 | 0.0 | 0.799 | 0.648 | 0.8110137672090112 | 11.402 | 18.61 | 1.6321697947728468 |  |
| 9 | 0.99 | 0.99 | 0.0 | 0.784 | 0.643 | 0.8201530612244898 | 15.526 | 20.952 | 1.349478294473786 |  |
| 10 | 0.99 | 0.99 | 0.0 | 0.771 | 0.642 | 0.8326848249027238 | 10.795 | 21.75 | 2.014821676702177 | background_yield_gt_2x_reference |
| 15 | 0.99 | 0.99 | 0.0 | 0.701 | 0.572 | 0.8159771754636234 | 13.445 | 37.356 | 2.7784306433618444 | background_yield_gt_2x_reference |
| 20 | 0.99 | 0.985 | -0.0050000000000000044 | 0.626 | 0.578 | 0.9233226837060702 | 18.443 | 80.128 | 4.344629398687848 | background_yield_gt_2x_reference |
| 25 | 0.985 | 0.985 | 0.0 | 0.644 | 0.489 | 0.7593167701863354 | 36.615 | 113.778 | 3.1074149938549773 | background_yield_gt_2x_reference;signal_eff_lt_0p8_reference |
| 30 | 0.98 | 0.98 | 0.0 | 0.674 | 0.503 | 0.7462908011869436 | 43.64 | 191.503 | 4.388244729605866 | background_yield_gt_2x_reference;signal_eff_lt_0p8_reference |

## Runtime Notes

- Runtime diagnostics were skipped by command-line option.
