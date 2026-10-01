# TLS Project

碩士論文: An Ensemble-Based Multi-Omics Feature Selection Framework for Cross-Cancer Identification of Tertiary Lymphoid Structure Signatures

整合 **RNA-seq、DNA Methylation、Somatic Mutations、Copy Number Variation (CNV)** 四種 omics 資料，透過各 omics 獨立前處理與特徵篩選（mRMR）→ concatenation → Ensemble 重要性排序（Random Forest / XGBoost / CatBoost）→ Consensus 篩選出 TLS signature → 建立最終預測模型。

支援的 `cancer_type`：`all`、`PAAD`、`CHOL`、`STAD`

---

## 目錄結構

```
TLS_project/
├── code/                          # 主要程式碼
│   ├── preprocessing/              # Step 1：各 omics 前處理（產生 gene count / beta value 等 matrix）
│   │   ├── preprocess_RNA.R
│   │   ├── preprocess_DNAmeth.R
│   │   ├── preprocess_DNAmut.R
│   │   └── preprocess_CNV/                 # CNV 前處理流程（CNVkit + GISTIC2.0）
│   │       ├── GDC_download.sh              # 下載 BAM 檔
│   │       ├── rename_bam.sh                # BAM 檔命名整理
│   │       ├── cnvkit_gistic.sh             # 單一樣本 CNVkit → GISTIC2.0 流程
│   │       └── cnvkit_gistic_batch.sh       # 批次執行版本
│   │
│   ├── align_sample/                # Step 2：跨 omics 樣本比對，產生 complete cases 及各 omics table
│   │   └── align_sample.ipynb
│   │
│   ├── feature_selection/           # Step 3：各 omics mRMR 特徵篩選 + Ensemble importance ranking + Consensus selection
│   │   ├── RNA_mrmr.py
│   │   ├── DNAmeth_mrmr.py
│   │   ├── DNAmut_mrmr.py
│   │   ├── CNV_mrmr.py
│   │   ├── ensemble_importance.py    # concat 四個 omics matrix，跑 RF/XGB/CatBoost 取得 feature importance
│   │   └── consensus_selection.py    # 選出 ≥2 個模型都選中的 feature，作為 TLS signature
│   │
│   ├── modeling/                     # Step 4：最終預測模型
│   │   ├── select_final_model.ipynb  # 觀察並評估模型表現，決定最終模型 
│   │   └── final_model.py            # 使用 TLS signature 建立最終模型
│   │
│   ├── utils/
│   │   └── model_utils.py            # 共用函式（讀寫、evaluation、模型相關 helper 等）
│   │
│   ├── visualization/                # 圖表繪製
│   │   ├── visualization.ipynb
│   │   ├── covariation.ipynb
│   │   ├── venn_plot.ipynb
│   │   └── plot_TMB.R
│   │
│   └── enrichment/                   # TLS signature 的 enrichment 分析
│       └── enrichment.ipynb
│
├── CNV/                           # CNV 資料與輸出檔案 (配合 CNV 的 preprocess)
│
├── R/                             # R 相關檔案與輸出結果 / 套件環境設定 (配合 RNA、DNAmeth、DNAmut 的 preprocess)
│
├── config/                        # 各步驟的參數設定檔（json）
│   ├── config_rna.json
│   ├── config_meth.json
│   ├── config_mut.json
│   ├── config_cnv.json
│   ├── config_ensemble.json
│   └── config_final_model.json
│
├── data/                          # 包含各個癌症的 raw matrix / metadata (not uploaded)
│
├── results/                       # Pipeline 輸出結果 (not uploaded)
│   └── {cancer_type}/               # e.g. PAAD, CHOL, STAD, all
│       ├── clustering/                # 各步驟執行 log
│       ├── fig/                       # 儲存圖檔
│       ├── importance_df/             # ensemble model 所計算的 feature importance
│       ├── logs/                      # 各步驟執行 log
│       ├── selected_features/         # 各步驟所選擇的特徵
│       ├── table/                     # 各 omics 及最終結果 table（如 rna_table_full.sqlite）
│       └── common_samples.npy         # 各 omics 的 complete cases
│
├── test/                          # 測試用 script (.ipynb)
├── tiff_figure/                   # 最終輸出圖檔 (not uploaded)
│
├── run_seed.sh                    # 可重複執行多個 random seed 的自動化腳本
├── script.sh                      # 主要 pipeline 執行入口（依序呼叫上述所有步驟）
├── workflow.tiff                  # 主要流程圖
└── TLS.yaml                       # Conda 環境設定檔
```

---

## Pipeline 流程說明

整體流程對應 workflow 圖片所示: 四大 omics 各自前處理 → 篩選 → mRMR → concat → ensemble importance → consensus → TLS signature → 最終模型：

| Step | 說明 | 對應程式 |
|---|---|---|
| 1 | 由 R 產生 sample sheet 並進行前處理，得到 gene counts、beta value 以及 mutation binary matrix | `code/preprocessing/preprocess_RNA.R`、`preprocess_DNAmeth.R`、`preprocess_DNAmut.R` |
| 1' | CNV 另外用 CNVkit 產生 log2 ratio by segments，再用 GISTIC2.0 得到 gene-level log2 ratio | `code/preprocessing/preprocess_CNV/*.sh` |
| 2 | 依各 omics 樣本取交集（complete cases），並輸出各 omics 對應 table | `code/align_sample/align_sample.ipynb` |
| 3 | 各 omics 分別做前處理後過濾（variance / MAD / mutation prevalence / clustering）+ mRMR 特徵篩選 | `code/feature_selection/{RNA,DNAmeth,DNAmut,CNV}_mrmr.py` |
| 4 | 合併四個 omics 篩選後的特徵矩陣，訓練 100 個 RF / 100 個 XGBoost / 100 個 CatBoost，取得各自 feature importance | `code/feature_selection/ensemble_importance.py` |
| 5 | 取至少被 2 種模型選中的 feature，作為最終 **TLS signature** | `code/feature_selection/consensus_selection.py` |
| 6 | 使用 TLS signature 特徵建立最終預測模型 | `code/modeling/final_model.py` |
| 7 | 視覺化與 enrichment 分析 | `code/visualization/*`、`code/enrichment/enrichment.ipynb` |

---

## 如何執行

### 1. 環境設定

```bash
conda env create -f TLS.yaml
conda activate TLS
```

### 2. 準備 config 檔

`config/` 底下各 json 對應各步驟參數（資料路徑、mRMR 選取的 feature 數量、模型超參數等）。
執行前請先確認 / 修改對應的 `config_*.json`（實際路徑欄位請填入自己環境的路徑，勿沿用他人絕對路徑，詳見下方注意事項）。

### 3. 執行完整 pipeline

依照 `script.sh` 的順序執行（各 omics mRMR → ensemble importance → consensus → 最終模型）：

```bash
conda activate TLS

python ./code/feature_selection/RNA_mrmr.py       --config ./config_rna.json     > ./results/PAAD/logs/RNA_mrmr_42.log
python ./code/feature_selection/DNAmeth_mrmr.py   --config ./config_meth.json    > ./results/PAAD/logs/DNAmeth_mrmr_42.log
python ./code/feature_selection/DNAmut_mrmr.py    --config ./config_mut.json     > ./results/PAAD/logs/DNAmut_mrmr_42.log
python ./code/feature_selection/CNV_mrmr.py       --config ./config_cnv.json     > ./results/PAAD/logs/CNV_mrmr_42.log
python ./code/feature_selection/ensemble_importance.py --config ./config_ensemble.json > ./results/PAAD/logs/ensemble_importance_mrmr_42.log
python ./code/feature_selection/consensus_selection.py --config ./config_ensemble.json > ./results/PAAD/logs/consensus_selection_mrmr_42.log

python ./code/modeling/final_model.py --config ./config_final_model.json > ./results/PAAD/logs/final_model.log
#python ./code/modeling/final_model.py --config ./config_final_model.json > ./results/all/logs/final_model_consensus.log
```

> 上述範例以 `PAAD` 為例，`cancer_type` 可替換為 `all` / `PAAD` / `CHOL` / `STAD`。

### 4. 多 random seed 重複執行（穩健性驗證）

`run_seed.sh` 用於自動化跑多個 seed，避免結果受單一 seed 影響：

```bash
# 用法
sh run_seed.sh [config_path] [start_seed] [end_seed] [model] [cancer_type]

# 範例：對四個 omics 及 ensemble 分別跑 seed 0~9
sh run_seed.sh ./config_rna.json      0 9 rna      PAAD
sh run_seed.sh ./config_meth.json     0 9 meth     PAAD
sh run_seed.sh ./config_mut.json      0 9 mut      PAAD
sh run_seed.sh ./config_cnv.json      0 9 cnv      PAAD
sh run_seed.sh ./config_ensemble.json 0 9 ensemble PAAD
```

腳本特性：
- `set -e`：任一指令失敗即中止，避免後續覆蓋錯誤結果
- 自動建立 `./results/{cancer_type}/logs` 及暫存資料夾 `./tmp/run_seed_{RUN_ID}`
- 若執行中被中斷（`Ctrl+C`），會自動清理暫存檔（`trap cleanup EXIT INT TERM`）

### 5. CNV 前處理（獨立流程，需先於主 pipeline 執行）

```bash
cd ./code/preprocessing/preprocess_CNV
sh GDC_download.sh          # 下載 BAM 檔
sh rename_bam.sh            # 整理 BAM 檔命名
sh cnvkit_gistic_batch.sh   # 批次跑 CNVkit → GISTIC2.0，產出 log2 ratio
```

---

## 各資料夾 / 檔案細節

### `code/preprocessing/`
由 R 產生 sample sheet 並進行 QC / normalization：
- `preprocess_RNA.R`：gene counts → log2(TPM+1)
- `preprocess_DNAmeth.R`：IDAT files（Grn/Red）→ preprocessFunnorm & QC → beta value
- `preprocess_DNAmut.R`：MAF files → non-synonymous mutation filtering → mutation binary matrix

### `code/align_sample/align_sample.ipynb`
比對四種 omics 的樣本 ID，取出 complete cases（四種 omics 都有的樣本），
並輸出每種 omics 對齊後的 table，供後續 feature selection 使用。

### `code/feature_selection/`
- `{RNA,DNAmeth,DNAmut,CNV}_mrmr.py`：各 omics 各自的前處理過濾（variance / MAD / mutation prevalence / clustering）+ mRMR 篩選，並依 80/20 切分 train/test
- `ensemble_importance.py`：合併四個 omics 篩選後特徵，訓練 100 個 Random Forest / 100 個 XGBoost / 100 個 CatBoost，輸出各模型的 feature importance
- `consensus_selection.py`：從三種模型的 top features 中，選出至少被 2 種模型選中的 feature，作為最終 TLS signature

### `code/modeling/final_model.py`
使用 consensus 選出的 TLS signature 特徵，建立最終預測模型並輸出評估結果。

### `code/utils/model_utils.py`
共用工具函式（模型訓練 / 評估 等，供各步驟 import 使用）。

### `code/visualization/`
- `visualization.ipynb`：主要結果視覺化
- `covariation.ipynb`：特徵間 covariation 分析圖
- `venn_plot.ipynb`：consensus feature 的 Venn diagram
- `plot_TMB.R`：Tumor Mutation Burden 相關圖表

### `code/enrichment/enrichment.ipynb`
TLS signature 的 pathway / gene set enrichment 分析。

### `preprocess_CNV/`
獨立於主 pipeline 之外的 CNV 前處理流程：GDC 下載 BAM → 命名整理 → CNVkit 產生 log2 ratio by segments → GISTIC2.0 得到 gene-level log2 ratio。輸出結果會被 `align_sample` 讀取整合。

---

## Config 檔說明

各步驟參數統一透過 `config/*.json` 傳入。常見欄位包括：

```json
{
   "input_table": "meth_table_full.sqlite",
   "fs_method": "mrmr",
   "cancer_type": "PAAD",
   "random_seed": 42,
   "model": ["RF", "XGB", "CatBoost"],
   "scaler_dict": {},
   "suffix": ""
}
```

> 使用者接手後請自行建立 / 修改對應 `config_*.json`

