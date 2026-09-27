# RESUME — MLPhoton 自訓（HZa_merged）

**接手先讀這份，再讀 `README.md`（細節與所有踩過的坑都在那）。**
最後更新 2026-08-30。此刻沒有執行中的 job（queue 裡的是 HZgamma 的）。

---

## ✅ 背景 friend 生產已完成，dataVmc 三張圖已出（2026-08-30）

四個 2024 背景 tag 全部生產完並通過對帳，dataVmc 的 all / SR / CR 都畫好了。

圖：`/eos/home-p/pelai/HZa/output_plots/merged_dataVmc/region{0_all,1_SR,2_CR}/`（各 29 張）
驅動：`Plot/scripts/run_merged_datavmc.sh`（**必須 LCG_104**，見下）

| tag | chunks | rows | pass_merged |
|---|---|---|---|
| DYGto2LG_10to100 | 382 | 8,870,919 | 41.20% |
| DYJetsTo2E | 748 | 28,800,210 | 91.83% |
| DYJetsTo2Mu | 744 | 1,998,611 | 26.57% |
| DYJetsTo2Tau | 460 | 14,140 | 82.46% |

join 後（`add_merged_flag.py`）：Data 8,113 / DYGto2LG 52,688 / 2E 24,746 / 2Mu 11,246 / 2Tau 35。

### 🔴 圖能看形狀，不能看 normalization

`ScaleBkgToData` 三個 region 給幾乎相同的因子：all 0.5627、SR 0.5712、CR 0.5591。
物理效應不會讓 SR 和 CR 的 data/MC 比一致到這種程度 —— 這是 **data 端 MLPhoton friend
覆蓋不足**：Data 只有 **47.5%** 的事件有 ML 資訊（MC 是 100%）。`8113/0.475 = 17,080`，
與 bkg 的 14,418 同量級。

而且缺得**不均勻**，不能用 lumi 縮放修：

```
run 379416-383811   60-75%
run 383811-385695   7.6-9.5%   <- 幾乎沒有
run 386323-386951   68.1%
```

要修就是把 Data 的 MLPhoton friend 補齊（第 1 項的 Data 重產會一併解決）。

### 這一輪踩到的坑（全部會偽裝成正常結果）

1. **兩個 friend base 有相同的 run:lumi:event** → join 的 searchsorted 取到舊的那份
   （無 ML 欄位）→ 補 0 → `DYJetsTo2E` 報 pass_merged **0**，而它真實通過率 91.83%。
   已改為**跳過**沒有 `pass_allcuts_merged_ML` 的 friend 檔，不是補 0。
   （⚠️ skip 要放在 `keys.append()` **之前**，否則 keys/vals 長度不一致→IndexError。）
2. **LCG_104 的 pyarrow 太舊**，讀 HiggsDNA 寫的 chunk 會回 `ArrowNotImplementedError`，
   看起來像「460 個 chunk 全損壞」，其實一個都沒壞。真損壞是 `OSError`。
   `reconcile_bkg_ml_friend.py` 已加版本守衛（<16 直接拒跑）。**用 hza_ana 的 pyarrow 21。**
3. **parquet metadata 正常不代表檔案好**：DYJetsTo2Tau 的 merged 檔 metadata 報 14,140 rows
   一切正常，實際讀資料 460 個 row group 有 284 個炸掉。對帳已加實際讀取檢查。
   chunk 是好的，重 merge 即可（壞檔改名保留為 `.corrupt_*`，沒刪）。
4. **繪圖環境**：`hza_ana` 是 ROOT 6.34，拒絕 `from ROOT import *`（`Analyzer_Configs` 就是這樣寫）；
   system ROOT 6.40 會讓 `AsNumpy` 靜默失敗。**用 LCG_104 (6.28)**，見 skill `env-setup`。
5. **CC 重啟會殺掉 `setsid` 的 supervisor**（我試過，擋不住）。已在 condor queue 的 job 不受影響。
6. **`/tmp` 被清會帶走 grid proxy** → `sample_manager.py:136` 丟一個**沒有訊息**的
   `RuntimeError`。AFS 有 `x509up_condor` 備份，supervisor 現在會自動恢復。
7. HiggsDNA 的 **`condor_submit` 沒有 timeout**（`managers.py:617`）→ driver 卡死 9.5 小時，
   期間程序活著、queue 空、held=0，只有 log 的 mtime 不動。watchdog/supervisor 已加 `log_idle`。

### 次要待辦

* DYJetsTo2Tau 的 merged 是我用 python 合併的，少了 HiggsDNA 在 merge 階段加的
  `process_id` / `weight_central_no_lumi` / `year` 三欄（190 vs 193）。下游沒用到，
  但若要一致可讓 HiggsDNA 重 merge 一次。
* DYto2Mu 的 DAS dataset 已長到 2975 檔，parent_map 是照 2959 建的（少 16 檔的對應）。
* **DYJetsTo2E 的 friend rows (28.8M) 比 2Mu (2.0M) 多 14 倍，但 resolved 端 2Mu 反而更多** ——
  方向相反，尚未解釋。join 以 resolved 為準所以不影響這批圖，但值得查。

---

## 一句話狀態

merged photon 的 classifier 與 mass regressor 已經**自己訓練並部署**（regressor=v4、classifier=v3），
signal 端完整驗證過；**但 Data 還沒用新模型重產，所以整條分析鏈還不能跑**。

---

## 現在裝了什麼

| 檔案 | 版本 |
|---|---|
| `RecoEgamma/EgammaMLPhotonProducers/data/regressor.onnx` | **v4** |
| `RecoEgamma/EgammaMLPhotonProducers/data/classifier.onnx` | **v3** |

同一份也已同步到 `CMSSW_15_0_14/src/RecoEgamma/.../data/`（**兩處都要換，只換一處 CMSSW 仍載入舊的**）。

備份（可直接 cp 回去回滾）：
* `*.onnx.bak_20260811_110829` — 原始 EXO-22-022 模型
* `regressor.onnx.bak_20260818_031438` — v3 regressor

產物（原始的都沒被動過）：
* `MLNanoAOD_v3/`、`MLNanoAOD_v4/`（各 529 檔，signal）
* `parquet_merged_DNA_v3/`、`parquet_merged_DNA_v4/`（各 9 個 sub-GeV 質量點）
* `train_runs/`（regressor_v1..v4、classifier_v1..v3、export_v3、export_v4、unified_eval*）
* 全部在 `/eos/cms/store/group/phys_susy/pelai/HZa_merged/`

`MergedAna/merged_p2root.py` 的 `ROI_WINDOWS` 已更新為 v4 值（舊值在 `merged_p2root.py.bak_preV4_*`）。

---

## 🔴 現在不能直接跑分析

`ROI_WINDOWS` 是 v4 推導的，但 `DATA_ML_GLOB` 指向的 data friend parquet **仍是原始模型產的**。
重訓後整個質量尺度移了（M0p1 median 0.334 → 0.201），把 v4 的窗套到舊模型的 data 上會選到
完全錯的譜段，**而且不會有任何錯誤訊息**。

所以 `merged_p2root.py` 的 `SIG_ML` 目前用 `MERGED_ML_VERSION` 控制，**預設 `old`**（維持自洽）。
要切到 v4 必須先完成下面第 1 項。

---

## 未完成的事（依相依性排序）

### 1. Data 重產（必要，才能真正上線）
189,600 個檔案。`MLNanoAOD/` 下的 68 個 Data 目錄要用新模型重跑，再跑 HiggsDNA friend。
工具現成：
```bash
# MLNanoAOD（改 OUT_EOS_BASE 指向新版本）
OUT_EOS_BASE=/eos/.../HZa_merged/MLNanoAOD_v4 \
  bash RegressMergedPhoton/condor/submit_mass_point.sh <tag> <DAS_dataset> data <year>
# friend parquet
bash HiggsDNA/scripts/run_merged_data_ml_friend.sh
```
完成後把 `merged_p2root.py` 的 `DATA_ML_GLOB` 指向新產物，並設 `MERGED_ML_VERSION=v4`。

### 2. 背景 merged ML friend（dataVmc 的前提）— 🔄 2026-08-29 已提交，跑在 condor 上

**不用重產 MLNanoAOD。** 背景的 MLNanoAOD 早就在（2026-05-31、舊模型），4 個 2024 樣本共
39,960 檔。缺的只是 friend 那一步，而它之所以缺是因為 **join 從來沒被打開過**：

* `metadata/za_merged_bkgmc_run3.json` 少了 Data config 有的四個 key
  （`mlphoton_friend` / `mlphoton_parent_map` / `mlphoton_dir` / `mlphoton_tag`）。
* 結果：friend parquet 263 欄、**零個 `MLPhoton_*`**，只有舊的 `pass_allcuts_merged_AN2020`。
* ⚠️ 這個失敗會偽裝成物理結果：`add_merged_flag.py` 把缺欄位補 0，log 印
  `pass_merged 0`，讀起來就是「背景沒有事件通過 merged selection」。
  實測 smoke test 單一檔案就有 **1855/5445 通過**。

修法（已做）：`scripts/run_merged_bkg_ml_friend.sh`，per-tag patch config，仿 data 版。
2,330 個 condor job（DYGto2LG 382 / 2E 748 / 2Mu 740 / 2Tau 460，fpo=4，四個 tag 依序）。
輸出 `parquet_friend_ML/`（**不覆蓋**舊的 `parquet_friend/`）。

#### 2b. 附帶挖到的坑：DAS 對某些 dataset 沒有 file-level parentage

`Bkg_DYGto2LG_10to100_2024` 的 parent map **1495 個 entry 全是 `[]`**。
`dasgoclient parent file=<lfn>` 對那個 dataset 回 **rc=0 但無輸出**（`DYJetsTo2E` 正常），
所以 `build_parent_map.py` 忠實地寫出一張全空的 map，沒有任何錯誤。

dataset-level parent 是有的 → 改用 **lumi 交集**重建：
`scripts/build_parent_map_lumi.py`（兩次 `file,lumi dataset=` 查詢，不做 per-file DAS 呼叫）。

closure 驗過才用：對 `DYJetsTo2E`（有真 DAS parentage）重建 →
**2990/2990 exact match、16323 links 與 DAS 版一字不差**。
重建後的 DYGto2LG：1527 parents / 8496 links / 0 empty（5.56 mini/nano，2E 是 5.46）。

`run_merged_bkg_ml_friend.sh` 內建 **0-link 守衛**：parent map 沒有任何 link 就拒跑，
不讓它再產出一批沒有 MLPhoton 的 parquet。

### 3. dataVmc（工具已備妥，等第 2 項）
```bash
# Data 的 friend 在 parquet_friend/，背景在 parquet_friend_ML/ -- 兩個都要給
python Plot/scripts/add_merged_flag.py \
  --friend /eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend,/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend_ML \
  --samples Data,DYGto2LG_10to100,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau --eras 2024
python Plot/scripts/plot_fast_variable_dataVmc.py --mergedOnly --region 0 \
  --input-dir /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_mergedflag   # 0=all 1=SR 2=CR
```
Data 2024 的 join 已驗證：匹配 94.6%、565 個通過 merged ML selection。
region: 0=all (95-180)、1=SR (115-135)、2=CR (sideband)。

`add_merged_flag.py` 現在會在 friend 沒有 `pass_allcuts_merged_ML` 時**大聲警告**，
而不是安靜地補 0 —— 就是上面那個偽裝成物理結果的失敗。

模型自洽性：這一版 data 與背景**都是舊模型**（背景 MLNanoAOD 2026-05-31、Data per-run
friend 同一版），所以比較有效。v4 版要等第 1 項 + 背景 MLNanoAOD 也重產。

### 4. ~~下一版模型：把 R1/R2/R3 餵進 regressor~~ ❌ 已驗證無效（2026-08-23）
**不要再試。** 所有形狀變數對 m/E 都沒有判別力：
* C++ 的 R1/R2/R3 不是 shower shape moments（`compute_En` 純幾何、無能量加權、量到影像原點的距離）
  → corr(log R1, log m/E) = **-0.015**
* 正確算的能量加權二階矩也一樣：λ₁ **0.015**、σ_major **0.015**、npix -0.171
* 原因：分離 = 115×(m/E) 個 crystal，m/E=0.002 → **0.23 顆**、中位數 0.00775 → **0.89 顆**。
  兩光子落在同一顆或相鄰 crystal，影像裡沒有雙峰可以量。

**v4 已接近此任務的資訊上限**（詳見 README §5h）。attenuation 是輸入資訊不足時 MSE 最佳解的
必然結果，不是訓練問題。要突破得換更細的偵測器資訊。

### 5. M0p1 質量分布（已完成，供詮釋用）
`/eos/home-p/pelai/HZa/output_plots/MLPhoton_M0p1/m0p1_mass_versions.{pdf,png}`
窗寬縮到 0.28 倍、中心到真值距離縮到 0.42 倍 → **不只是先驗收縮**。
⚠️ 但峰仍在 0.15–0.2（真值 0.1），改善的是「先驗正確、尾巴收斂」而非「質量量準」。
描述 M0p1 的質量重建能力時措辭要小心。

---

## 成果摘要（判斷值不值得繼續投入時看這個）

* **classifier**：AUC 0.5665 → **0.7261**（舊模型幾乎等於隨機 0.5）。固定 hadronic 污染下
  diphoton 效率 **+74%**。通過選擇的事件數增加 4–69 倍。
* **regressor v4**：ROI 窗 vs 舊模型 — M0p1 **0.28**、M0p2 **0.24**、M0p3 **0.27**、M0p4 0.49、
  M0p5 0.72、M0p6 0.94；但 M0p7–M0p9 是 1.31/1.60/1.83（**比舊模型寬**）。
  交叉點在 M0p6 —— 贏在 merged 不可取代的 sub-GeV，輸在 resolved 本來就更強的區間。
* signal-only 靈敏度：保守假設下 M0p1 11.3×、M0p9 1.5×（區間太寬，不能當結論）。

---

## 會靜默失敗的陷阱（每一個都真的發生過）

1. **模型有兩份副本** — 只換 `RegressMergedPhoton/` 的話 CMSSW 仍載入舊模型，一切看起來正常。
2. **friend join 靜默跳過** — `mlphoton_parent_map` 給相對路徑時 loader 找不到檔案就跳過，
   parquet 照寫，只是少了所有 `MLPhoton_*`（263 欄 vs 290 欄），全程無錯誤。
   ⚠️ `mlphoton_parent_map` 要**絕對**路徑、`samples.catalog` 要**相對**路徑（慣例相反）。
3. **catalog 指向不存在的路徑** → 0 jobs → `ZeroDivisionError` in progress bar，訊息完全不提缺輸入。
4. **MLNanoAOD 是 18-branch slim friend**，不是完整 nano。catalog 要指 NanoAODv15 dataset。
5. **cmsRun 寫 EOS fuse 會 segfault** 在 `TTree::Fill`，看起來像自己的 code 壞了。寫本地再 cp。
6. **等待條件別用 `pgrep -f <script>`** — 會匹配到查進度的互動 shell，迴圈永不結束（白等一小時）。
   改成等產物（輪詢檔案大小，連續兩次相同才算寫完）。
7. **DAS 暫時性失敗 + `set -e`** 會讓整批中斷（M0p2 中招，dataset 其實好好的）。要重試。
8. **比較一定要用同一個 val / 同一個 selection** — 換了 val 集就會得到相反結論（v1/v2 那次）。

---

## 檔案地圖

```
RegressMergedPhoton/
  training/
    RESUME.md                 <- 這份
    README.md                 <- 完整記錄（§5c 訓練 / §5e v3 / §5f v4 / §5g dataVmc）
    cluster_py.py             Cluster.cc + DoPairings 的 numpy 移植（怪癖照抄，別「修好」）
    closure_test.py           preprocessing closure（改模型後必跑）
    run_closure.sh            一鍵 closure（已改成寫本地再 cp）
    models.py                 PyTorch 重建 + ONNX 載入/匯出
    verify_models.py          模型層 closure
    make_training_set.py      dump ROOT -> npz（含 e_frac 判準）
    build_dataset.py          打包（能量平衡 / log-flat / 分層切分 / m/E 視窗）
    summarize_training_set.py class 統計
    train_classifier.py       選模用 macro recall
    train_regressor.py        log-space head，選模用 log loss
    compare_classifiers.py    ROC 比較（比固定 cut 公平）
    bias_scan.py              完整分析鏈的 bias/解析度
    derive_roi_windows.py     重推 ROI_WINDOWS（換模型後必跑）
    estimate_sensitivity.py   signal-only 靈敏度（兩個極端假設）
    export_onnx.py            .pt -> ONNX + 三重檢查 + --install（自動備份、兩處都換）
    deploy_chain_v4.sh        完整驗證鏈（closure -> MLNanoAOD -> HiggsDNA -> bias/ROI）
    condor/                   訓練與 dump 的提交/監控腳本
  RecoEgamma/.../data/        部署中的 .onnx + 所有備份
Plot/scripts/
  add_merged_flag.py          merged 旗標 join 進 resolved ntuple
  plot_fast_variable_dataVmc.py  加了 --mergedOnly
HiggsDNA/scripts/
  run_merged_signal_v3.sh     per-mass-point config 生成 + 重試容錯
MergedAna/merged_p2root.py    ROI_WINDOWS(v4) + MERGED_ML_VERSION 開關
```

相關 memory：`project_hza_merged_mlphoton_retrain`、`ref_cern_condor_gpu_compute_capability`、
`ref_cmsrun_eos_fuse_write_segfault`。
