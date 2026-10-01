# feature_selection/CNV_model.py 使用說明

本文件說明 `CNV_model.py` 的用途、設定檔參數、主要流程與輸出結果。

## 腳本位置

- `feature_selection/CNV_model.py`

## 功能概述

`CNV_model.py` 用於建立 CNV-based TLS 分類模型，核心流程如下：

1. 讀取 JSON config。
2. 載入對應 cancer type 的 CNV `PivotTable`。
3. 進行特徵前處理（移除 chrX/chrY、保留 protein-coding 基因）。
4. 建立/重用 `split_<random_seed>` 的 train/test 分割。
5. （可選）執行 clustering 評估與 best-k 選擇。
6. 以 best-k 對基因做 re-clustering，並聚合成 cluster-level 特徵。
7. 先跑 cluster-level 全特徵基準評估。
8. 依 `fs_method` 做特徵選擇（`nvar` / `mrmr` / `rf`）。
9. 用選後特徵進行最終模型評估。

## 執行方式

在 `TLS_model` 根目錄執行：

```bash
python feature_selection/CNV_model.py --config config_cnv.json
```

若你有其他設定檔，也可改為：

```bash
python feature_selection/CNV_model.py --config <your_config>.json
```

## 必要 config 欄位

腳本會檢查下列必要鍵值：

- `input_table`
- `selected_sample`
- `min_clusters`
- `max_clusters`
- `random_seed`
- `cancer_type`
- `save_feature_importance`
- `save_top_features`

## 常用可選欄位

- `clustering`: 是否重跑 clustering，預設 `False`
  - `True`: 重新跑 `k_fold_clustering_evaluation`
  - `False`: 直接讀既有 `metric_df_{random_seed}.csv`
- `fs_method`: 特徵選擇方法，預設 `rf`
  - 可用值：`nvar`、`mrmr`、`rf`
- `fs_cutoff`: 程式中有保留欄位，預設 `0`（目前分支中未實際用於篩選）

## 輸入資料來源

- CNV table（預設路徑）：
  - `results/{cancer_type}/table/cnv_table.sqlite`
- 若 `cancer_type != "all"`：
  - 會改讀 `results/{cancer_type}/table/{input_table}`，再用 `selected_sample` 做 sample 子集
- sample 子集（npy）：
  - 由 `selected_sample` 指定
- 基因註解（protein coding）：
  - `data/metadata/gencode_anno_v36.csv`

## 輸出結果

### 1. split 寫回

若 `sample_metadata` 尚無 `split_<random_seed>`，會新增並寫回：

- `results/{cancer_type}/table/cnv_table.sqlite`

### 2. Clustering 相關

- k-fold cluster labels：
  - `results/{cancer_type}/clustering/cluster_labels_{min_clusters}_{max_clusters}_{random_seed}.json`
- 指標表（silhouette/ARI/coherence）：
  - `results/{cancer_type}/clustering/metric_df_{random_seed}.csv`
- majority vote labels：
  - `results/{cancer_type}/clustering/majority_vote_dict_{random_seed}.json`
- penalty sweep 圖：
  - `results/{cancer_type}/clustering/penalty_sweep_{random_seed}.png`
- best-k 指標圖：
  - `results/{cancer_type}/clustering/find_best_k_plot_{random_seed}.tiff`
- best-k 重新分群標籤：
  - `results/{cancer_type}/clustering/best_k_cluster_labels_{random_seed}.csv`

### 3. Cluster table

- 聚合後的 cluster-level table：
  - `results/{cancer_type}/table/cnv_cluster_table_{random_seed}.sqlite`

### 4. 特徵重要度（`fs_method=rf`）

- `results/{cancer_type}/importance_df/{save_feature_importance}_rf_{random_seed}.csv`

### 5. 選後特徵（npy）

- `nvar`:
  - `results/{cancer_type}/selected_features/{save_top_features}_nvar_{random_seed}.npy`
- `mrmr`:
  - `results/{cancer_type}/selected_features/{save_top_features}_mrmr_{random_seed}.npy`
- `rf`:
  - `results/{cancer_type}/selected_features/{save_top_features}_rf_{random_seed}.npy`

## 模型與評估

- Base model 由 `utils/model_utils.py` 的 `build_base_model("RandomForest", ...)` 建立。
- 評估由 `evaluate_model(...)` 進行。
- `best_k` 由 `sweep_penalty_select_best_k` 與 `plot_clustering_metrics_and_find_best_k` 共同決定。

## 目前程式中的固定值

- `mrmr` 分支固定 `mrmr_k = 20`
- `nvar` 分支測試特徵數：`[2, 4, 6, ..., 20]`

## 相依套件

- `pymaftools`
- `numpy`
- `pandas`
- `scikit-learn`
- `mrmr-selection`
- `matplotlib`
- `tqdm`

## 注意事項

1. `clustering=True` 計算時間很長（程式註解顯示約數小時）。
2. 若 `clustering=False`，必須先有既存的 `metric_df_{random_seed}.csv`，否則讀檔會失敗。
3. 程式內使用多個絕對路徑（`/home/data/data_YM313112004/TLS_model/...`），若搬移專案請同步調整。
4. `cancer_type="all"` 時會 merge CHOL/PAAD/STAD 的 CNV table 並寫回 `results/all/table/cnv_table.sqlite`。
