# Run3 limit improvement — progress and next steps

最後更新：2026-05-26

## 1. 起點問題

`flashggFinalFit/Plots/plot_limits/compare_xs_run2.pdf` 顯示 Run3 (172 fb⁻¹) limit
在低 ma 區域比 Run2 (138 fb⁻¹) 差很多：
- ma=1: Run3 ≈ 33 fb vs Run2 17.8 fb （差 ~85%）
- ma=2: Run3 ≈ 5.95 fb vs Run2 4.74 fb
- ma=3-30: 大致跟 Run2 ±10–20%

理論上 √(138/172) ≈ 0.89，加上 13.6 TeV ggH 截面 ~1.04 × 13 TeV，
Run3 應該比 Run2 好 ~15%，但實際在低 ma 端反而退化。

## 2. 已完成的診斷（全 6 項）

| # | 診斷 | 結論 |
|---|---|---|
| 1 | Signal σ_eff vs ma | ma=1 μ σ_eff=**3.76** GeV，e σ_eff=**5.77** GeV — 異常寬，σ_eff/ma 比值 ~4×ma |
| 2 | Signal yields | ma=1 yield ~198（vs ma=4 的 491），低 ma yield 塌陷 |
| 3 | BDT AUC | Run3 AUC = **0.960**，Run2 = **0.99** — Run3 BDT 變差 |
| 4 | bkg envelope dof | ma=1 用 **Bern5** (6 par)，ma=2 Bern5、ma=3 Bern4 — 低 ma 自由度過高，spurious-signal 風險高 |
| 5 | Systematics lnN | Run3 系統與 Run2 接近，不是主因 |
| 6 | Impact (ma=1) | `shapeBkg_norm` impact_r=**0.134** — bkg norm 在低 ma 大幅吃 sensitivity |

最終發現：低 ma σ_eff 變寬源自 **low-side tail**（photon-merging + electron bremsstrahlung
在 mγγ 小開角下，整個 m_llγγ 譜被拉低）。ee channel 比 μμ tail 更嚴重。

## 3. 已實作的修正

### (A) Signal fit — `Signal/scripts/signalFit.py`

| 修正 | 內容 |
|---|---|
| Mass-dependent fit range | ma≤1: `[118,135]`；ma=2: `[115,135]`；ma=3-4: `[110,140]`；其他維持 `[100,180]` |
| Low-ma 用 DCB | ma≤4 強制 `useDCB=True`（DCB+Gaussian，power-law tail 描述 brem） |
| `xvar.setRange()` | 確保 dataset reduce 後不含 tail |

副作用：DCB 的 pdf component key 非數字，舊 `plotPdfComponents` 的 `int(k[-1])` 會炸 →
已修 `Signal/tools/plottingTools.py` 用 try-except 容錯。

**驗證結果**（ma=1-4 σ_eff 重 fit）：

| ma | μμγγ σ_eff 改善 | eeγγ σ_eff 改善 |
|---|---|---|
| 1 | 3.76 → **2.27** GeV (**−40%**) | 5.77 → **3.50** GeV (**−39%**) |
| 2 | 1.91 → 1.70 GeV (−11%) | 2.84 → 2.62 GeV (−8%) |
| 3 | 1.79 → 1.69 GeV (−6%) | 2.78 → 2.67 GeV (−4%) |
| 4 | 1.74 → 1.63 GeV (−6%) | 2.72 → 2.73 GeV (~0%) |

### (A 下游) 跑了 ma=1-4 的 calcPhotonSyst + datacard + combine + plot

| ma | OLD Run3 | **NEW Run3** | Run2 | NEW/Run2 |
|---|---|---|---|---|
| 1 | 33.0 fb | **19.8 fb** | 17.8 fb | 1.11× |
| 2 | 5.95 fb | 5.62 fb | 4.74 fb | 1.19× |
| 3 | 3.20 fb | **2.29 fb** | 2.78 fb | **0.82×（打贏 Run2）** |
| 4 | 2.62 fb | 2.59 fb | 2.69 fb | 0.96× |

### (C) BDT — 加 4 個 derived features

| feature | 公式 | 物理意義 |
|---|---|---|
| `pho_pt_asym` | (pT₁−pT₂)/(pT₁+pT₂) | merged-γ 時 reco pT 不平衡 |
| `pho_dR_over_ma` | dR(γγ)/m_a | quasi-invariant；merged case 偏差大 |
| `min_pho_pt_over_ma` | min(pT)/m_a | 低 ma 的 soft photon 區分 |
| `ma_resid_norm` | (m_a − m_a,hyp)/m_a,hyp | normalized mass residual，跨 ma 一致 |

**修改的檔案**：
- `HZaMVA/scripts/run3_Za_BDT.py`：加 variables list、`add_derived_features()`、plot metadata、`bdt_input_branches` 過濾 derived
- `HZaMVA/scripts/1_make_sideband_reweight.py`：`REWEIGHT_VARS` 加 4 個 derived；`load_frame` 內計算 derived
- `HZaMVA/scripts/sideband_reweight.py`：`_dataframe_values`/`_object_value` 支援 derived，新增 `_derived_*` helpers
- `Parquet2Rootfile/Parque2Root_BDT.py`：`MVA_FEATURE_COLUMNS` 加 4 個 derived；inference 時計算 derived
- `Plot/scripts/1_prepare_dataVmc.py`：histogram booking、fill、systematics var_map、`pdfName_map` 編號 15-18
- `Plot/scripts/data_vmc_fast.py`：RDataFrame `Define()` 4 個 `rdf_*` 欄位；`base_specs` 加 4 個 HistSpec
- `Plot/lib/Plot_Configs.py`：`var_title_map` 加 4 個 LaTeX title

**BDT train 結果**（新 BDT，含 4 derived var + 新 sideband JSON）：
- Best AUC: **0.9725**（vs 舊 0.960，提升 1.3%）
- Top 5 trials AUC ≥ 0.970
- Selected: `ks_regularized`，max_depth=3, lr=0.06, n_estimators=900
- Saved to `HZaMVA/using/model_Za_BDT_run3.pkl`

### (D) Background fit envelope cap

`Background/test/fTest_ALP_turnOn.cpp` patched：限制 Bernstein max order ≤ 4（其他 family 維持 ≤ 6）。
- Patch 1: `while (prob<0.05 && order < (((*funcType)=="Bernstein") ? 5 : 7) ){`
- Patch 2: `while (prob<upperEnvThreshold && order < (((*funcType)=="Bernstein") ? 5 : 7) ){`
- Binary rebuilt: `Background/bin/fTest_ALP_turnOn`（2026-05-25 15:18）

預期效果：ma=1,2 不再用 Bern5 → 降低 spurious signal，shapeBkg_norm impact 應下降。

## 4. 目前正在跑的 pipeline

由兩個背景 orchestrator 串接：

### orchestrator 1: `/tmp/run_reweight_bdt_p2root_dataVmc.sh` (PID 2454233)

執行步驟：
1. ✅ Rebuild sideband reweight JSON（含 4 derived var）
2. ✅ Train BDT — AUC=0.9725
3. ⏳ **進行中**：P2Root MVA score + MVAcut + dataVmc plot

目前狀態：在 dataVmc resubmit attempt 2/5，**81 個 outputs missing**。
- Resubmit jobs 都 normal-terminate (exit 0) 但 output ROOT 不可讀 `(no_root_keys_or_unreadable)`
- 受影響 jobs 都是 `SR sideband_rwgt` 區的 data 和 DY samples
- 需要進一步查 jobs stderr/stdout 看為什麼 ROOT 沒寫出來

### orchestrator 2: `/tmp/run_flashgg_chain.sh` (PID 2604414)

在 wait `run_reweight_bdt_p2root_dataVmc.sh` 完成後接 flashgg chain：
4. flashgg MVA cut
5. Tree2WS
6. **bkg fit (Condor，Bern≤4 patched binary)**
7. signal fit (含 low-ma DCB) + photonSyst
8. datacard + combine
9. **final limit plot**（含 `compare_xs_run2.pdf`）

## 5. 結論

- 已找到並修正 Run3 ma=1-4 limit 變差的兩個直接根因（low-side tail、低 ma bkg envelope 太彈）
- 已重 train BDT 並 build 新 sideband reweight，AUC 從 0.960 → **0.9725**
- 進行中：用新 BDT 重 score、重做下游 flashgg chain，預期所有 ma 點 limit 都會改善
- ma=1 期望最終 limit ~15-20 fb（vs Run2 17.8 fb），ma=3 期望可以 < 2.5 fb

## 6. 接下來要做的事

### 短期（current pipeline 跑完之前）

- [ ] **修 dataVmc 81 個 missing outputs**：查 SR sideband_rwgt jobs 為什麼 ROOT 寫不出來。
  - 命令：看 `Plot/Condor/logs/dataVmc.8305007.*.err` 與 condor stdout
  - 可能原因：memory not enough（4000 MB），或新加的 derived feature columns 在 ntuples 沒被解析正確
- [ ] 等 P2Root score + MVAcut + dataVmc 完成後，**檢視 4 個 derived var 的 data/MC plots**
  - `Plot/plots/variables_dataVmc/15_pho_pt_asym.pdf` ... `18_ma_resid_norm.pdf`

### 中期（flashgg chain）

- [ ] 等 orchestrator 2 完成，看新 limit `Plots/plot_limits/compare_xs_run2.pdf`
- [ ] 驗證 bkg fit 的 envelope（`background_fit_summary_*.csv`）確認低 ma 不再選 Bern5
- [ ] 確認 ma=1 的 impact 中 `shapeBkg_norm` 不再 dominant

### 完成後

- [ ] 跑 Impact plot 確認最終 dominant nuisance（`RUN_FLASHGG_IMPACT=1`）
- [ ] 跑 Bias study（`RUN_FLASHGG_BIAS=1`）
- [ ] 更新 AN：`RUN_UPDATE_AN=1` 或手動 git push `AN-25-172` 與新圖
- [ ] 寫 Run2 vs Run3 比較段落：σ_eff 改善、BDT AUC、low-ma sensitivity

### 風險與待確認

- DCB low-ma signal fit 在 calcPhotonSyst 階段是否會引入新的 systematics shape 問題（要看 photon scale/smear 對 DCB α/n 參數的響應是否合理）
- Bernstein cap 至 4 是否會在某些 ma 造成 truth-model 找不到收斂（fTest 失敗）→ 看新 `background_fit_summary_*.csv`
- 新 BDT 對低 ma signal-eff 是否 trade off bkg-rej（看新 ROC、yields）

## 7. 重要檔案路徑

- 完整流程腳本：`HiggsZaAna/1_grand.sh`
- 環境：`source ~/setup_hza.sh`（含 conda activate）
- Reweight JSON：`HZaMVA/reweights/sideband_run3_iterative.json`
- BDT model：`HZaMVA/using/model_Za_BDT_run3.pkl`
- fTest binary：`flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/bin/fTest_ALP_turnOn`
- 進行中的 log：
  - `/tmp/run_reweight_bdt_p2root_dataVmc.log`
  - `/tmp/run_flashgg_chain.log`
- 最終 limit plot 將在：`flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Plots/plot_limits/`
