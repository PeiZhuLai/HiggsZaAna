# MLPhoton 重訓 framework（HZa_merged）

自行重新訓練 merged-photon 的 **classifier** 與 **mass regressor**（即
`RecoEgamma/EgammaMLPhotonProducers/data/{classifier,regressor}.onnx`）。

## 1. 現有模型的來源

現用的兩個 ONNX 來自 **EXO-22-022**（AN2020-159，Rutgers/Notre Dame），上游 repo
[atownse2/MLPhotons](https://github.com/atownse2/MLPhotons) **只有 CMSSW inference plugin，不含任何訓練程式碼**。
訓練樣本是該組私有的 particle gun MC，AN 未給 DAS 路徑，Rucio 要不到。

模型結構（由 ONNX 反解，pytorch 1.7 匯出、opset 9）：

| | classifier | regressor |
|---|---|---|
| 輸入 | `img[1,1,30,30]`（**原始能量影像，未正規化**） | `img[1,1,30,30]`（**除以總和正規化**）+ `eta[1,1]` |
| 骨架 | Conv(1→64,7×7)-ReLU-Conv(64,5×5)-ReLU-MaxPool-Conv(64,4×4)-ReLU-Conv(64,3×3)-MaxPool | 同結構但 16 filters |
| 頭部 | FC(4096→3) + LogSoftmax | FC(784→64)-LeakyReLU-concat(eta)-FC(65→16)-LeakyReLU-FC(16→1)-LeakyReLU |
| 輸出 | mono / di / hadron 三類分數 | `m/E` |

AN 的訓練設定：Adam 預設（lr=1e-3, β=0.9/0.999），classifier 350 epochs、regressor 500 epochs。
regressor 訓練樣本 500k（+50k 驗證），classifier 600k（三類各 200k）。

## 2. 為什麼要重訓：現有模型在 HZa 相位空間已經失效

AN Table 7 明載 regressor 的 **m/E 有效訓練範圍是 [0.005, 0.07]**，且 gun 的能量範圍是
**75 < E < 2500 GeV**。而 HZa merged 訊號的 cluster 能量中位數只有 **37–42 GeV**。

對 `Sig_MC_MLNANO_all` 的 truth-matched 事件實測（2026-08-09）：

| tag | m_gen | N | E_cluster | 真 m/E | 預測 m/E | m_rec 中位數 | m_rec/m_gen | IQR/med | 真 m/E < 0.005 |
|---|---|---|---|---|---|---|---|---|---|
| M0p1 | 0.100 | 485 | 36.6 | 0.00273 | 0.00852 | 0.334 | **3.34** | 2.25 | **93.4%** |
| M0p2 | 0.200 | 612 | 39.2 | 0.00510 | 0.00862 | 0.369 | **1.84** | 1.77 | 48.5% |
| M0p3 | 0.300 | 1078 | 37.2 | 0.00806 | 0.01050 | 0.395 | **1.32** | 1.55 | 17.1% |
| M0p4 | 0.400 | 1931 | 36.6 | 0.01091 | 0.01266 | 0.450 | 1.13 | 0.94 | 3.3% |
| M0p5 | 0.500 | 3667 | 38.1 | 0.01311 | 0.01446 | 0.534 | 1.07 | 0.66 | 1.3% |
| M0p6 | 0.600 | 5223 | 39.6 | 0.01515 | 0.01597 | 0.626 | 1.04 | 0.52 | 0.5% |
| M0p7 | 0.700 | 6545 | 40.7 | 0.01719 | 0.01751 | 0.723 | 1.03 | 0.41 | 0.4% |
| M0p8 | 0.800 | 7719 | 42.0 | 0.01906 | 0.01919 | 0.821 | 1.03 | 0.35 | 0.2% |
| M0p9 | 0.900 | 8511 | 42.1 | 0.02137 | 0.02120 | 0.917 | 1.02 | 0.31 | 0.1% |

結論：**M0p5 以上模型正常（ratio ≤ 1.07），M0p3 以下崩壞**。原因是真 m/E 掉出訓練域下界 0.005，
模型被「地板夾住」—— 不論真值多低，預測 m/E 都卡在 ~0.0085。M0p1 因此被重建成 0.334 GeV（高估 3.3 倍），
解析度（IQR/median）也從 0.31 惡化到 2.25。這正是 M0p1 的 per-mA ROI 窗會爛成 (0.154, 0.947) 的原因。

重製指令：`python bias_scan.py`（見 §6）。

## 3. 訓練樣本策略：零 particle gun

不複製論文的 private gun，改用**現有樣本**：

| 用途 | 來源 |
|---|---|
| regression target + diphoton class | HZa signal MiniAOD（M0p1…M10, RunIII2024Summer24MiniAODv6） |
| monophoton class | GJets / DY MC 中 truth-match 到 gen prompt γ 的 cluster |
| hadron class | GJets / QCD 中 match 不到 gen γ 的 cluster |

雖然 m_gen 只有離散幾個質量點，但 **E_cluster 連續分布在 28–74 GeV**，
所以 `m/E = m_gen / E_cluster` 天然連續覆蓋約 0.001–0.05，正好把現有模型死掉的低端補滿。
代價：模型只在 HZa 的相位空間可信，不通用。

## 4. 檔案

| 檔案 | 說明 |
|---|---|
| `cluster_py.py` | `Cluster.cc` + `DoPairings` + `calculateLorentzVector` 的 numpy 逐行移植，加上 ONNX 推論封裝 |
| `closure_test.py` | 比對「offline python 重現」vs「線上 C++ MLPhoton」，逐 cluster 檢查 11 個量 |
| `run_closure.sh` | 一鍵：DAS 找檔 → cmsRun → 比對 |
| `bias_scan.py` | §2 的 bias 表，重訓前後都要跑 |
| `make_training_set.py` | dump ROOT → `(image, eta, energy, label, m/E)` 的 npz |
| `condor/run_one_dump.sh` | condor worker：跑 `closure_cfg.py`，輸出到 EOS |
| `condor/submit_training_dump.sh` | 單一 dataset 的 stage + 提交 |
| `condor/submit_all_training.sh` | 全 campaign 提交（`DRY_RUN=1` 可先看不送） |
| `condor/monitor_and_convert.sh` | 背景監控 + drain 後自動轉檔 |
| `summarize_training_set.py` | class 統計、能量分布、m/E 覆蓋、現有 classifier 效率 |
| `models.py` | PyTorch 重建 classifier/regressor（與部署圖逐層一致）+ ONNX 權重載入/匯出 |
| `verify_models.py` | **模型層 closure**：torch == onnxruntime == C++ 存的值，含 export round-trip |
| `build_dataset.py` | 打包：classifier 能量平衡、regressor log(m/E) 平坦化，按檔案切 train/val |
| `train_classifier.py` | 分類器訓練（監控 diphoton recall，非 accuracy） |
| `train_regressor.py` | 迴歸器訓練（相對誤差 loss，逐 m/E band 報 bias/resolution） |
| `export_onnx.py` | `.pt` → ONNX + 三重檢查 + `--install`（自動備份舊檔） |
| `../../CMSSW_15_0_14/src/RegressMergedPhoton/RecHitDumper/plugins/EBRecHitDumper.cc` | RecHit + gen truth + 線上 MLPhoton mirror 的 dumper |
| `../../CMSSW_15_0_14/src/RegressMergedPhoton/RecHitDumper/test/closure_cfg.py` | 同一個 job 裡跑 producer + dumper |

## 5. 流程

```
[0] closure                                          <- 已 PASS（2026-08-09）
    run_closure.sh M1 200
    必須 PASS 才能往下走。offline 前處理若跟線上有任何差異，
    訓練出來的模型就是在學一個不存在的偵測器。

[1] 訓練資料生產（condor，CERN）   <- 進行中，2026-08-09 提交 156 jobs
    condor/submit_all_training.sh          （DRY_RUN=1 可先看不送）
    condor/monitor_and_convert.sh          （背景監控，queue drain 後自動轉 npz）
    註：cmsRun 讀遠端 MiniAOD 必須在 CERN 跑 —— IHEP batch node 沒有 grid CA。

[2] 訓練（IHEP GPU）
    PyTorch 重建上表結構，可載入現有 ONNX 權重當 warm start。
    classifier 與 regressor 分開訓。

[3] 匯出 + 部署 + 驗證
    export_onnx.py --install        （自動備份舊 .onnx）
    詳細步驟見 §7。
```

## 5d. 部署（[3]）—— 別低估這一步的成本

`export_onnx.py --install` 只是把 .onnx 換掉。**換完之後 `/eos/.../MLNanoAOD/` 底下所有既有產物
都是舊模型算的**，包括下游的 `parquet_merged_DNA_tmp`、`parquet_friend`、per-mA ROI 窗、
以及 `merged_p2root.py` 產出的所有 flashgg tree。要讓分析真正吃到新模型，必須：

1. `export_onnx.py --install`（會備份成 `*.onnx.bak_<timestamp>`）
2. `bash run_closure.sh M1 200` —— 換模型後重跑 preprocessing closure。
   這一步驗的不是「新模型比較好」，而是**新的 ONNX 在 CMSSW 裡真的被載入且離線可重現**。
3. 重跑 MLNanoAOD 生產（signal + data + bkg，約 500+ condor jobs）
4. 重跑 HiggsDNA merged tagger → 新的 `MLPhoton_*` / `MergedML_*`
5. `bias_scan.py` 重新量測 —— 這時才是在**完整分析鏈**上的數字，
   而不是 val 集上的數字。ROI 窗（`merged_p2root.py` 的 `ROI_WINDOWS`）也必須依新的
   signal 分布重新推導，舊的 q16/q84 是舊模型的。

也就是說：驗證集上的 bias/resolution 只是決定「值不值得部署」的依據；
真正的收益要到 (5) 才知道。在模型還沒明顯勝過現況前不要啟動 (3)。

目前完成：**[0] closure 已通過**。

### closure 結果（2026-08-09, M1, 200 事件）

輸入 `/eos/home-p/pelai/HZa/root_rechit/closure_M1_200ev_numEvent200.root`

```
events with matching cluster count : 200/200
cluster pairs compared             : 359
offline clusters with no online match: 0

     variable  max reldiff
         mass    2.163e-05
          moe    2.162e-05
           pt    5.341e-07
          eta    1.388e-06
          phi    6.710e-07
diphotonScore    8.837e-08
    r1/r2/r3   ~3e-07
RESULT: PASS
```

最大偏差 2×10⁻⁵ 出現在 mass/moe，屬 float32 往返誤差量級（C++ 用 float、python 用 float64）。
clustering、影像建構、兩個模型的輸入配置全部逐 cluster 一致。

## 5b. 訓練資料生產（[1]）

### Campaign（2026-08-09 提交）

| block | dataset | 檔/dataset | jobs | ~事件 |
|---|---|---|---|---|
| signal | `HZa-Zto2L-ato2G_Par-M-{0p1..0p9,1..10}`（19 點） | 6 | 114 | 2.2M |
| gjets | `GJ-4Jets_Bin-HT-{40to100,100to200,200to400}-PTG-10to100` | 8 | 24 | 1.2M |
| qcd | `QCD_Bin-PT-{30to50,50to80,80to120}_Fil-EMEnriched` | 6 | 18 | 0.9M |
| | | | **156** | **4.3M** |

輸出 `/eos/cms/store/group/phys_susy/pelai/HZa_merged/train_dump/<tag>/`，
npz 在 `.../train_npz/<tag>/`，log 在 `condor/stage/<tag>/logs/`、`condor/stage/monitor.log`。

GJets/QCD 用低 PTG / 低 PT 分箱，因為 HZa 的 cluster 能量只有 ~40 GeV；
高 HT 分箱的 photon 太硬，訓練起來對不上我們的相位空間。

### 小批量驗證（提交前，M1 signal 200 事件 + GJets 200 事件）

signal：264 clusters（diphoton 78 / mono 5 / hadron 181，另 36 個 dropped）。
GJets：122 clusters（diphoton **0**，正確；mono 68 e_frac 中位數 0.967；hadron 54）。

**依現有 classifier 分數拆解 M1 的 diphoton：**

| diphotonScore | N | m_rec 中位數 | m_rec/m_gen |
|---|---|---|---|
| >0.9 | 14（17.9%） | 0.9865 | **0.986** |
| 0.5–0.9 | 3 | 0.9520 | 0.952 |
| <0.5 | **61（78%）** | 0.2654 | **0.265** |

也就是說：**現有 classifier 只認出 17.9% 的真 a→γγ cluster**；被它認出的那些，
regressor 幾乎完美（ratio 0.986）。§2 的 bias 表看起來 M0p5 以上「正常」，是因為那張表只統計
通過完整分析選擇（含 diphotonScore 要求）的事件 —— 它看不到分類器整批漏掉的 78%。
所以重訓的主要收益很可能在 **classifier**，不在 regressor。

### 兩個在小批量測試中抓到的坑

1. **只用 dR 做 label 會毀掉 regression target**。ALP 的能量會被 clustering 切出碎片，
   碎片落在 dR<0.1 內卻只帶部分能量 → 全額 gen mass 配上部分 cluster 能量，target m/E 被灌大。
   實測 M1：純 dR 標記給出 target m/E 中位數 0.033，加上能量含量準則後是 0.030，
   而現有模型預測 0.012 —— 差別足以讓整個訓練走偏。
   修法：`--e-frac-lo/--e-frac-hi`（預設 0.7–1.3），不符者標為 ambiguous 直接丟棄，
   **不要丟進 hadron class**（否則等於教分類器「EM 碎片是強子」）。
2. **`ml_*` 不能用 eta 匹配**：`ml_eta` 是經過 PV 修正的 `etaprime`，不是 cluster eta。
   用 eta 比對會 96% 失敗。線上 collection 與 offline cluster 同序（closure 已證），直接用 index。

### 產出統計（2026-08-09 完成，156/156 job，4.17M clusters）

`summarize_training_set.py` 的結果。npz 在 `/eos/cms/store/group/phys_susy/pelai/HZa_merged/train_npz/`。

| class | N | 佔比 | AN2020-159 用量 |
|---|---|---|---|
| monophoton | 387,553 | 9.3% | 200k |
| diphoton | 702,800 | 16.8% | 200k（分類）/ 500k（迴歸） |
| hadronic | 3,082,531 | 73.9% | 200k |
| dropped | 284,191 | (6.4% of matched) | — |

三類統計都超過論文用量，**不需要補產**。

**能量分布（GeV）** — 三類差很多，訓練前必須做能量平衡：

| class | q05 | q16 | med | q84 | q95 |
|---|---|---|---|---|---|
| mono | 12.09 | 15.26 | 25.89 | 54.79 | 88.44 |
| di | 16.73 | 22.52 | 34.63 | 58.42 | 91.95 |
| had | **1.24** | 2.32 | 16.74 | 50.30 | 83.92 |

**迴歸目標 m/E 覆蓋** — 全體 signal 落在 [0.000108, 0.963]，
相對舊訓練域 [0.005, 0.07]：**14.6% 在下界外、17.8% 在上界外**。

| tag | N_di | m/E q05 | m/E med | m/E q95 | <0.005 | >0.07 |
|---|---|---|---|---|---|---|
| M0p1 | 82382 | 0.00117 | 0.00295 | 0.00613 | **88.8%** | 0.0% |
| M0p2 | 43434 | 0.00234 | 0.00588 | 0.01194 | 36.6% | 0.0% |
| M0p5 | 85691 | 0.00574 | 0.01457 | 0.03036 | 3.4% | 0.0% |
| M1 | 64399 | 0.01090 | 0.02917 | 0.06173 | 0.5% | 2.7% |
| M2 | 17420 | 0.02085 | 0.05876 | 0.12538 | 0.0% | **35.8%** |
| M5 | 14503 | 0.04841 | 0.14443 | 0.26738 | 0.0% | **88.7%** |
| M10 | 7584 | 0.08635 | 0.23516 | 0.45729 | 0.0% | **96.9%** |

⚠️ **新發現：外插不只發生在低端。** §2 只掃了 M0p1–M0p9，看不到這件事 ——
**mA ≥ 2 GeV 有 36%–97% 超出舊訓練域上界 0.07**。也就是說現有模型在
HZa merged 用到的 mA 1–10 GeV 大部分區域同樣是外插，不只 sub-GeV 端有問題。

**現有 classifier 對真 diphoton 的效率**（要補的洞）：

| tag | M0p1 | M0p3 | M0p5 | M1 | M2 | M5 | M10 |
|---|---|---|---|---|---|---|---|
| eff(score>0.9) | 6.9% | 15.5% | 18.8% | 23.9% | **27.2%** | 11.3% | 9.3% |
| eff(score>0.5) | 11.7% | 20.2% | 23.5% | 28.4% | 31.5% | 16.9% | 16.6% |

呈駝峰狀，峰值只有 27%：低質量太 collimated（看起來像單光子）、
高質量兩個光子分太開（看起來不像 merged）。整條曲線都遠低於可用水準。

### 訓練階段要注意

hadron class 的能量中位數只有 ~7 GeV，signal/mono 是 27–45 GeV。
AN §5.5.1 明確做了能量平衡（逐 bin 取最少類的事件數），重訓時必須照做，
否則分類器會退化成能量判別器。

## 5c. 訓練（[2]）

### 三層 closure 都必須 PASS

| 層 | 腳本 | 驗的是什麼 | 狀態 |
|---|---|---|---|
| preprocessing | `closure_test.py` | python clustering/影像 == 線上 C++ | ✅ 2.2e-5 |
| 模型 | `verify_models.py` | torch == onnxruntime == C++ 存的值 | ✅ 6.8e-5 |
| export | `verify_models.py`（round-trip） | torch → ONNX opset9 → onnxruntime | ✅ 6.1e-5 |

網路結構是從 ONNX 的 node attributes 精確還原的，不是猜的。注意
**classifier 的 conv4 是 pad=2 而 regressor 是 pad=1**，這個不對稱在部署圖裡就有，照抄。

### baseline：現有 regressor 在完整 diphoton 母體上的表現

`train_regressor.py --eval-only --head deployed --warm-start <舊 regressor.onnx>`

| m/E band | N | bias | res |
|---|---|---|---|
| [0.000, 0.005) | 256 | **1.282** | 0.657 |
| [0.005, 0.010) | 887 | 1.093 | 0.414 |
| [0.010, 0.020) | 2875 | 0.962 | 0.489 |
| [0.020, 0.050) | 10705 | **0.366** | 0.457 |
| [0.050, 0.100) | 7470 | **0.102** | 0.147 |
| [0.100, 0.200) | 4584 | **0.043** | 0.055 |
| [0.200, 1.000) | 6675 | **0.019** | 0.019 |

只有訓練域 [0.005, 0.07] 內 bias≈1；域外輸出被鎖死，高質量端把質量壓到真值的 2%。
重訓後用同一個表比對。

### ⚠️ 兩個踩過的坑（都已修，別再犯）

**1. relative MSE 會 collapse，不能當訓練目標。**
`L = ((pred-target)/target)²` 看似是跨數量級迴歸的正解，實際上預測 0 的代價恰好是 1.0（有界），
過衝的代價卻大得多 → 最佳化器把輸出停在 0。實測從零訓練收斂到 bias=-0.001、val loss 1.00，
**比未訓練的現有模型（0.63）還差**。這是 loss 的性質不是資料的問題。
改用 `MLPhotonRegressorLog`：網路輸出 `exp(raw·σ+μ)`（恆正），用 MSE 訓練標準化的 `raw`。
匯出的圖仍只回傳一個 m/E 純量，CMSSW 不用改。`relative_mse` 只留作報告用。

**2. 按檔案切 train/val 必須「依 tag 分層」。**
每個 npz 只來自單一 dataset，所以不分層的檔案切分會把整個質量點整批分到某一邊。
第一次打包的結果：**val 完全沒有 M0p1–M0p6**，m/E 中位數 0.058 vs train 0.0077 ——
驗證集描述的是另一個族群，early stopping 會挑錯模型。
修正後 val 涵蓋全部 19 個質量點（moe_med 0.00866 vs train 0.01063）。
按檔案切本身是對的（同事件的 cluster 在同一檔，row-wise 切會洩漏），錯的是沒分層。

### 算力：CERN condor GPU

CPU 實測 **63 分鐘/epoch**（466k 樣本、16 threads 已用滿），200 epochs = 8.7 天，不可行。
模型很小，GPU 上估計 30–60 秒/epoch。

**⚠️ CERN GPU 節點的 compute capability 陷阱。** GPU 節點的 driver 是 580.173.02
（DriverVersion 13.0），所以 `hzgdna` 的 torch 2.11.0+cu130 版本相容 —— 但 **V100 不行**：

```
Found GPU0 Tesla V100-PCIE-32GB which is of compute capability (CC) 7.0.
This version of PyTorch was built for: 7.5, 8.0, 8.6, 9.0
```

torch 2.11 官方 wheel 已放棄 Volta。`torch.cuda.is_available()` 仍回 **True**，
但實際跑 kernel 會失敗 —— probe 的 matmul benchmark 直接沒有輸出。
**這正是會靜默浪費一整天訓練的那種坑：能拿到 GPU、能 import、`is_available()` 是 True，然後什麼都算不出來。**

CERN condor 的 GPU 分布（`condor_status -constraint 'TotalGpus > 0' -long | grep DeviceName`）：

| GPU | CC | 數量 | torch 2.11 |
|---|---|---|---|
| A100-PCIE-40GB | 8.0 | 100 | ✅ |
| H100 NVL | 9.0 | 78 | ✅ |
| V100S-PCIE-32GB | 7.0 | 54 | ❌ |
| V100-PCIE-32GB | 7.0 | 41 | ❌ |
| T4 | 7.5 | 26 | ✅ |
| H100L MIG | 9.0 | 38 | ✅ |
| H200 | 9.0 | 18 | ✅ |

所以 JDL 一定要加：

```
request_GPUs = 1
require_gpus = Capability >= 7.5
```

備案（若 CC≥7.5 的機器排隊太久）：改用 env **`hza_NN`**（torch 2.6.0+cu124，仍支援 Volta CC 7.0），
V100 就能用。代價是 `hza_NN` 沒有 onnx 套件 → 不能用 `--warm-start <onnx>`（改從零訓練或用 `.pt`），
匯出一律回到 CERN 前台用 `hzgdna` 做（export 不需要 GPU）。

### v1 訓練結果（2026-08-09）

CERN condor GPU，**10–15 秒/epoch**。regressor 在 ep50 early stop、classifier 在 ep37。

**regressor（同一個 val 的公平對照）**：

| m/E band | N | 舊 bias | v1 bias | 舊 res | v1 res |
|---|---|---|---|---|---|
| [0.000, 0.005) | 30642 | **6.397** | **1.806** | 5.656 | 2.093 |
| [0.005, 0.010) | 6227 | 1.210 | 1.237 | 0.966 | 1.622 |
| [0.010, 0.020) | 6200 | 0.959 | 1.031 | 0.585 | 0.868 |
| [0.020, 0.050) | 7089 | 0.370 | **1.202** | 0.452 | 1.235 |
| [0.050, 0.100) | 4445 | 0.095 | **0.372** | 0.152 | 0.737 |
| [0.100, 0.200) | 4627 | 0.041 | **0.217** | 0.051 | 0.298 |
| [0.200, 1.000) | 11718 | 0.014 | **0.168** | 0.024 | 0.110 |

每個 band 的 bias 都更接近 1，但高質量端仍低估約 6 倍 —— 方向對，尚未可用。

⚠️ **比對一定要用同一個 val**。§5c 那張 baseline 是在舊的、未分層的 val 上算的，
換成分層 val 後舊模型最低 band 的 bias 是 6.397 而不是 1.282。拿兩個不同 val 的數字相減會得到相反結論。

**classifier**：

| 指標 | 舊模型 | v1 |
|---|---|---|
| accuracy | 0.3879 | **0.5455** |
| argmax recall (di) | 0.192 | **0.363** |
| di_recall(>0.5) | 0.1888 | **0.2689** |
| di_recall(>0.9) | 0.1455 | 0.1131 |
| had_fake(>0.9) | 0.04481 | **0.00640** |

舊模型在能量平衡後的 val 上 accuracy 只有 0.388（三類均衡的隨機值是 0.333），
且把 76% 的 cluster 判成 hadronic —— §5b 那張「7–27%」是在未平衡的自然樣本上算的，掩蓋了這件事。
confusion matrix 也顯示 44% 的真 diphoton 被判成 monophoton，這是低質量 a 的本質困難。

### classifier：唯一公平的比較是 ROC，不是固定 cut

`di_recall(>0.9)` 對兩個模型意義不同 —— 舊模型把分數推到兩極（76% 判成 hadronic），
0.9 這個閾值在兩個模型上切到的不是同一件事。正確的問法是：
**在相同的 hadronic 污染下，哪個模型留下比較多真 diphoton？**（`compare_classifiers.py`）

```
AUC (di vs rest):  deployed 0.5747   v1 0.7072   v2 0.7056
```

**舊模型的 AUC 0.575 幾乎等於隨機（0.5）** —— 它在 HZa 的相位空間基本上沒有判別能力。

固定 hadronic fake rate 下的 diphoton 效率：

| hadronic fake rate | deployed | v1 | v2 |
|---|---|---|---|
| 0.0010 | 0.0000 | **0.0376** | 0.0285 |
| 0.0050 | 0.0609 | **0.1031** | 0.0923 |
| 0.0100 | 0.0808 | **0.1358** | 0.1286 |
| 0.0200 | 0.1043 | 0.1775 | 0.1770 |
| 0.0448 | 0.1457 | **0.2365** | 0.2349 |
| 0.1000 | 0.2142 | **0.3222** | 0.3201 |

在同一污染水準下效率是舊模型的 **1.5–1.7 倍**（0.0448 那列 +62%）。
v1 與 v2 實質相同（AUC 0.707 vs 0.706），v2 的 argmax di-recall 略高（0.388 vs 0.363）。

### 🔴 第四個、也是最嚴重的坑：訓練樣本裡有物理上不可能的東西

v1/v2 的 regressor 在高 m/E 端怎麼訓都停在 bias ~0.2，原因不是模型容量，是**樣本本身錯了**。

開角 θ ≈ 2·m/E，而 `DoPairings` 的 `MATCH_DR = 0.15` 決定了兩束能量沉積會不會被併成一個
cluster。所以 **m/E > 0.075 時，兩個光子根本不可能在同一個 cluster 裡**：

| m/E band | m_gen 中位數 | θ=2m/E | crystal 分離 | cluster 像素數 |
|---|---|---|---|---|
| [0.000, 0.005) | 0.10 | 0.0017 | 0.1 | 30 |
| [0.020, 0.050) | 0.80 | 0.0607 | 3.5 | 20 |
| [0.100, 0.200) | 5.00 | 0.2798 | 16.1 | 19 |
| [0.200, 1.000) | **9.00** | **0.9972** | **57.3** | **17** |

最高 band 的兩個光子相隔 57 個 crystal，cluster 卻只有 17 個像素。
**那是「單一光子的影像」配上「兩個光子的不變質量」當 target** —— 網路被要求做不可能的事。
它們之所以通過 `e_frac ∈ [0.7,1.3]`，是因為高質量 a 的能量分配夠不對稱時，
單一光子的能量仍接近 a 的總能量。

用 gen photon 直接量測驗證（不是幾何推論）：

| 樣本 | gen dR(γγ) 中位數 | cluster 含兩光子 | 只含一個 | 都不含 |
|---|---|---|---|---|
| sig_M0p5 | 0.0482 | 16 | 15 | 10 |
| sig_M9 | 0.8134 | **0** | 32 | 7 |

M9 沒有任何一個 cluster 同時含兩個 gen photon。

**修正**：`build_dataset.py --moe-max 0.075`（預設開啟），對 **兩個網路都生效** ——
未合併的 pair 是單光子影像，留著當 diphoton 會同時教壞分類器和迴歸器。
實測丟掉 120,099 列（佔 diphoton class 的 17.1%），剩 582,701。

這也跟分析本身一致：`merged_p2root.py` 的 `ROI_WINDOWS` 只有 M0p1–M0p9，
**merged 分析本來就只用 sub-GeV**，mA ≥ 1 走 resolved。把 M1–M10 放進訓練從一開始就不該做。

⚠️ 這個切法仍不完美：即使 M0p5 也有 37%（15/41）的 cluster 只含一個光子。
若 v3 的低 band 仍不夠好，下一步是改用嚴格的 gen 判準（要求兩個 gen photon 都落在 cluster 影像內），
而不是繼續調訓練超參數。

### 三個版本的最終對照（全部在同一個 merged-only val）

`unified_eval/`。舊模型與 v1/v2/v3 都在 `train_packed_merged` 的 val（60,556 列 regressor /
65,113 列 classifier）上重評 —— **不同 val 的數字不能並排**，v1/v2 原本是在含 m/E>0.075 的 val 上訓練與評估的。

**regressor bias**：

| m/E band | N | 舊 | v1 | v2 | **v3** |
|---|---|---|---|---|---|
| [0.000, 0.005) | 29016 | **6.816** | 1.025 | 1.164 | 1.649 |
| [0.005, 0.010) | 9400 | 1.223 | 1.236 | 1.422 | **0.994** |
| [0.010, 0.020) | 9872 | 0.968 | 1.031 | 1.189 | 0.798 |
| [0.020, 0.050) | 9576 | 0.367 | 1.284 | 1.460 | 0.561 |
| [0.050, 0.100) | 2692 | 0.116 | 0.708 | 0.766 | 0.266 |
| 整體 | | 1.636 | 1.068 | 1.220 | **1.045** |

**regressor resolution**（不可校正，比 bias 重要）：

| m/E band | 舊 | v1 | v2 | **v3** |
|---|---|---|---|---|
| [0.000, 0.005) | **9.499** | 1.316 | 1.792 | **0.954** |
| [0.005, 0.010) | 0.945 | 1.629 | 2.070 | **0.586** |
| [0.010, 0.020) | 0.592 | 0.828 | 1.072 | **0.383** |
| [0.020, 0.050) | 0.451 | 1.348 | 1.554 | **0.320** |
| [0.050, 0.100) | 0.262 | 1.139 | 1.238 | 0.266 |
| 整體 | 4.555 | 1.240 | 1.581 | **0.892** |

**classifier**：

| | 舊 | v1 | v2 | **v3** |
|---|---|---|---|---|
| AUC (di vs rest) | 0.5665 | 0.7210 | 0.7224 | **0.7261** |
| argmax di recall | 0.2082 | 0.4101 | **0.4432** | 0.4145 |

固定 hadronic fake rate 下的 diphoton 效率：

| fake rate | 舊 | v1 | v2 | v3 |
|---|---|---|---|---|
| 0.0010 | 0.0000 | **0.0578** | 0.0377 | 0.0466 |
| 0.0050 | 0.0813 | **0.1252** | 0.1174 | 0.1182 |
| 0.0100 | 0.1020 | **0.1696** | 0.1638 | 0.1628 |
| 0.0200 | 0.1271 | **0.2218** | 0.2187 | 0.2200 |
| 0.0448 | 0.1671 | 0.2914 | **0.2919** | 0.2873 |
| 0.1000 | 0.2288 | 0.3852 | **0.3853** | 0.3778 |

**結論：部署 v3。**
regressor 的 resolution 在每個 band 都最好（最低 band 從舊模型的 9.499 降到 0.954，改善 10 倍）；
classifier 三版實質相同（AUC 差 0.005），但 v3 是唯一沒有把「未合併的單光子」當 diphoton 訓練的，
物理定義最乾淨。

取捨要講清楚：**v3 的 bias 有 m/E 相依性**（1.65→0.99→0.80→0.56→0.27），v1 的 bias 較平但
resolution 差 2–4 倍。選 v3 的理由是 per-mA ROI 窗本來就會依實際重建分布重新推導，
**bias 可以被 ROI 窗吸收，resolution 不能**。

### ⚠️ 第三個坑：選模準則要跟訓練目標一致，且對應物理需求

v1 的兩個模型都踩到：

* **regressor** 訓練用 log MSE，但 early stopping 用 relative MSE。後者被最小的 target 主導，
  一路從 19.5 爬到 37.5，而實際質量尺度（bias）始終在 1.0–1.3 —— 等於隨機挑 epoch。
  已改成選 log-space loss（`evaluate(..., head=)`）。
* **classifier** 用 val NLL 選到 ep12，但 diphoton recall 一路升到 ep35 才停
  （NLL 被網路已經學會的類別主導）。已改成選 **macro-averaged recall**（三類 recall 平均），
  且 `ReduceLROnPlateau` 必須跟著改 `mode="max"`，否則 scheduler 會往反方向調。

## 5e. v3 部署結果（2026-08-14）—— 一個明確的取捨，不是全面勝出

完整分析鏈上的量測（`bias_scan_v3.txt` / `bias_scan_old.txt` / `roi_windows_v3.txt`，都在
`/eos/cms/store/group/phys_susy/pelai/HZa_merged/`）。

| tag | N（舊→新） | ratio（舊→新） | IQR/med（舊→新） | ROI 寬度比 |
|---|---|---|---|---|
| M0p1 | 485 → **20186**（42×） | 3.34 → **2.14** | 2.25 → **0.96** | **0.29** |
| M0p2 | 612 → **42389**（69×） | 1.84 → **1.20** | 1.77 → **0.80** | **0.29** |
| M0p3 | 1078 → **52685**（49×） | 1.32 → **1.01** | 1.55 → **0.76** | **0.34** |
| M0p4 | 1931 → 54014（28×） | 1.13 → 0.93 | 0.94 → **0.73** | 0.60 |
| M0p5 | 3667 → 51116（14×） | 1.07 → 0.88 | 0.66 → 0.74 | 0.82 |
| M0p6 | 5223 → 46680（8.9×） | 1.04 → **0.83** | 0.52 → **0.77** | 1.02 |
| M0p7 | 6545 → 42148（6.4×） | 1.03 → **0.80** | 0.41 → **0.81** | **1.37** |
| M0p8 | 7719 → 38830（5.0×） | 1.03 → **0.78** | 0.35 → **0.85** | **1.66** |
| M0p9 | 8511 → 36008（4.2×） | 1.02 → **0.76** | 0.31 → **0.88** | **1.91** |

**贏**：M0p1–M0p3 是舊模型完全不可用的區域，v3 把 bias 從 3.34/1.84/1.32 降到 2.14/1.20/1.01，
解析度改善 2–2.3 倍，ROI 窗收窄到 **29%**。統計量增加 42–69 倍（classifier 效率的直接效果）。

**輸**：M0p6–M0p9 舊模型本來就好（bias 1.02–1.04、IQR/med 0.31–0.52，正好在它的訓練域
[0.005, 0.07] 內），v3 在那裡系統性低估 20–24%、解析度退步 2–3 倍、**ROI 窗變寬到 1.9 倍**。

**原因**：log-flat 重採樣涵蓋 m/E 0.0001–0.075，而 M0p1 佔了訓練集 38%（162373/429052），
模型被拉去照顧極低質量端，犧牲了 0.015–0.022 那段（= M0p6–M0p9）。

**下一步建議**：merged 分析實際只用 mA 0.1–0.9，對應 m/E 約 **0.002–0.03**。目前訓練的
0.0001–0.075 兩端都有模型不需要學的區域在搶容量。把重採樣改成覆蓋實際需要的區間，
應該能同時保住兩端。這只需重訓（~1 小時 GPU + 打包），不必重產任何 MLNanoAOD。

**目前部署狀態**：`data/*.onnx` 已經是 v3（`export_onnx.py --install`，兩處都換），
舊模型備份在 `*.onnx.bak_20260811_110829`，回滾隨時可以。
`MLNanoAOD_v3/` 與 `parquet_merged_DNA_v3/` 是 v3 產物；
原始的 `MLNanoAOD/` 與 `parquet_merged_DNA_tmp/` 未被動過。

### 部署途中踩到的五個坑（都已在腳本裡修掉）

1. **模型有兩份副本** —— `RegressMergedPhoton/` 與 CMSSW area。只換前者的話 CMSSW 仍載入舊模型，
   而且看起來一切正常。`--install` 現在會自動找出所有 CMSSW area 一起換。
2. **catalog 指向已不存在的 `/eos/project/...`** → 檔案列表空 → 0 jobs →
   `ZeroDivisionError: division by zero`（progress bar 算 `completed/all`）。錯誤訊息完全沒提到缺輸入。
3. **MLNanoAOD 是 18-branch slim friend，不是完整 NanoAOD**。把 catalog 指向它會得到
   `no field 'Photon' in record with 5 fields`。catalog 要指 **NanoAODv15 dataset**，
   MLPhoton 靠 friend join 進來。dataset 名稱從 parent map 的 key 推導，確保與 friend map 一致。
4. **`mlphoton_parent_map` 要絕對路徑，`samples.catalog` 要相對路徑** —— 兩者慣例相反
   （前者直接 `open()`，後者走 `expand_path()` 會加 repo 前綴）。parent map 給相對路徑時
   **loader 靜默跳過 friend join**：parquet 照寫，只是少了所有 MLPhoton 欄位（263 欄 vs 290 欄），
   全程無錯誤。
5. **DAS 暫時性失敗會讓整批中斷** —— M0p2 第一次失敗（dataset 其實好好的，DAS 上有 24 檔），
   `set -e` 讓後續 7 個質量點全部沒跑。現在每個質量點重試 3 次，失敗也只跳過該點。

另外：**等待條件不要用 `pgrep -f <script>`** —— 它會匹配到任何命令列含該字串的 shell，
包括用來查進度的互動指令，導致等待迴圈永遠不結束（實測白等一小時）。改成等產物
（輪詢檔案大小，連續兩次相同視為寫完）。

## 6. 執行環境

| 工作 | 環境 |
|---|---|
| cmsRun（dumper / producer） | `CMSSW_15_0_14`, `SCRAM_ARCH=el9_amd64_gcc12` |
| python（onnxruntime + uproot） | `/eos/home-p/pelai/App/Conda/.conda/envs/hzgdna/bin/python` — 是唯一同時有 onnx 1.21 / onnxruntime 1.26 / torch 2.11 的 env |
| parquet 讀取（bias 表） | `hza_ana` |

## 7. 已知陷阱

* **`makeImage()` 的怪癖必須照抄，不要「修好」它**。它做 `ieta += ieta_min + isize/2`（加最小值而非減），
  只因為下一步減掉 `cxval` 才剛好抵銷；影像是以**成員清單的第一顆晶體**為中心，不是最高能量那顆，
  所以成員順序有意義（`combine()` 把 C2 的晶體接在 C1 之後）；x 直接 clamp，y 則先平移再 clamp。
  這些全都寫進 `cluster_py.py` 了，改動任何一條都會讓訓練與線上脫鉤。
* **classifier 吃原始影像，regressor 吃正規化影像**。在 C++ 裡這點被 `img` 是 900×900 vector
  但只用第 0 列的寫法掩蓋掉了。
* **`std::sort` 不穩定**：`|p|²` 完全相等的 cluster，線上與 offline 的排序可能不同。float 能量下罕見，
  `closure_test.py` 會把中招的事件報出來。
* **dumper 的 `minHitEnergy` 要維持 0.0**：producer 的 zero suppression 是 `E > 0`，
  dumper 若設更高的門檻會靜默破壞 closure。
* **pfIsolation 沒有移植**（需要 packedPFCandidates，dumper 未存），因此不在 closure 比對項目內。
* **cmsRun 不要用 VarParsing 的 `outputFile`**：它會把檔名悄悄改成 `<stem>_numEvent<N>.root`，
  呼叫端無法預測路徑（腳本會找不到自己剛產的檔）。`closure_cfg.py` 因此加了 `outFile=`，原樣使用。
