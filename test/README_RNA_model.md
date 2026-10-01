# feature_selection/RNA_model.py 使用說明

本文件說明 `RNA_model.py` 的用途、設定檔參數與輸出結果。

## 腳本位置

- `feature_selection/RNA_model.py`

## 功能概述

`RNA_model.py` 用於建立 RNA-based TLS 分類模型，主要流程如下：

1. 讀取 JSON config。
2. 載入對應 cancer type 的 RNA `PivotTable`。
3. 進行特徵前處理（protein coding、移除 PAR_Y、移除 chrX/chrY）。
4. 建立/重用 `split_<random_seed>` 的 train/test 分割。
5. 先跑全特徵基準評估。
6. 依 `fs_method` 執行特徵選擇（`nvar` / `mrmr` / `rf`）。
7. 使用選後特徵重新評估模型。

## 執行方式

在 `TLS_model` 根目錄執行：

```bash
python feature_selection/RNA_model.py --config config.json
```

## 必要 config 欄位

以下鍵值為腳本會檢查的必要欄位：

- `input_table`
- `selected_sample`
- `random_seed`
- `cancer_type`
- `save_feature_importance`
- `save_top_features`

## 常用可選欄位

- `fs_method`: 特徵選擇方法，預設 `rf`
  - 可用值：`nvar`、`mrmr`、`rf`
- `fs_cutoff`: 前篩 variance cutoff，預設 `0.25`

## 輸入資料來源

- RNA table（預設路徑）：
  - `results/{cancer_type}/table/rna_table.sqlite`
- sample 子集（npy）：
  - 由 `selected_sample` 指定
- 基因註解：
  - `data/metadata/gencode_anno_v36.csv`

## 輸出結果

### 1. 特徵重要度（`fs_method=rf`）

- `results/{cancer_type}/importance_df/{save_feature_importance}_rf_{random_seed}.csv`

### 2. 選後特徵（npy）

- `nvar`:
  - `results/{cancer_type}/selected_features/{save_top_features}_nvar_{random_seed}.npy`
- `mrmr`:
  - 前篩特徵：
    - `results/{cancer_type}/selected_features/{save_top_features}_preselect_{random_seed}.npy`
  - mRMR 選後特徵：
    - `results/{cancer_type}/selected_features/{save_top_features}_mrmr_{random_seed}.npy`
- `rf`:
  - `results/{cancer_type}/selected_features/{save_top_features}_rf_{random_seed}.npy`

### 3. split 寫回

若 `sample_metadata` 尚無 `split_<random_seed>`，腳本會新增並寫回：

- `results/{cancer_type}/table/rna_table.sqlite`

## 模型與評估

- Base model 由 `utils/model_utils.py` 的 `build_base_model("RandomForest", ...)` 建立。
- 評估由 `evaluate_model(...)` 進行，會輸出 train/test 預測與報表。

## 相依套件

- `pymaftools`
- `numpy`
- `pandas`
- `scikit-learn`
- `mrmr-selection`

## 注意事項

1. `cancer_type="all"` 會嘗試 merge CHOL/PAAD/STAD 的 RNA table，並存成 `results/all/table/rna_table.sqlite`。
2. `mrmr` 分支目前固定 `mrmr_k = 150`。
3. 路徑多為絕對路徑（`/home/data/data_YM313112004/TLS_model/...`），若搬移專案需同步調整。
