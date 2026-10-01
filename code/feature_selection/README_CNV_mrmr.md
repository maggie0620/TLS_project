# `CNV_mrmr.py` 使用說明

本文件說明 `CNV_mrmr.py` 的 CNV 特徵處理、CNV gene clustering、mRMR 特徵選擇與 TLS 分類流程。

## 腳本位置

- `code/feature_selection/CNV_mrmr.py`

## 功能概述

腳本以 CNV gene-level `PivotTable` 建立 TLS 分類模型，流程如下：

1. 讀取 JSON config。
2. 載入 CNV table，並取出共同樣本。
3. 移除 chrX/chrY，保留 protein-coding genes。
4. 建立或重用 `split_<random_seed>` 的 80/20 分層 train/test split。
5. 依 `clustering` 設定執行或讀取 gene clustering 結果。
6. 以 cosine similarity 和 agglomerative clustering 產生 cluster-level CNV table。
7. 先評估 cluster-level 全特徵模型，再依 `fs_method` 選擇特徵。
8. 使用選出的特徵重新評估 Random Forest TLS 分類器。

## 執行方式

在 `TLS_model` 根目錄執行：

```bash
python code/feature_selection/CNV_mrmr.py --config config_cnv.json
```

也可以指定其他設定檔：

```bash
python code/feature_selection/CNV_mrmr.py --config <config>.json
```

## Config

必要欄位：`input_table`、`clustering`、`min_clusters`、`max_clusters`、`random_seed`、`cancer_type`。

可選欄位：

- `fs_method`: `mrmr`、`nvar` 或 `rf`，預設為 `mrmr`。
- `fs_cutoff`: 保留欄位，預設為 `0`；CNV 的 mRMR 分支目前未用此值做前篩選。
- `suffix`: 加在 selected feature 檔名中的字串，預設為空字串。
- `clustering`: `true` 時重跑 k-fold clustering；`false` 時讀取既有的 clustering metric。
- `min_clusters` / `max_clusters`: clustering 搜尋範圍，預設為 2/350。

預設設定檔為 `config_cnv.json`。

## 輸入資料

當 `cancer_type` 不是 `all` 時，腳本使用：

- CNV table：`results/{cancer_type}/table/{input_table}`
- 共同樣本：`results/{cancer_type}/common_samples.npy`
- gene annotation：`data/gencode_anno_v36.csv`

當 `cancer_type` 為 `all` 時，會合併：

- `results/CHOL/table/cnv_table.sqlite`
- `results/PAAD/table/cnv_table.sqlite`
- `results/STAD/table/cnv_table.sqlite`

## 輸出結果

所有 `{cancer_type}` 路徑均相對於 `TLS_model` 根目錄。

### Clustering

- `results/{cancer_type}/clustering/cluster_labels_{min_clusters}_{max_clusters}_{random_seed}.json`
- `results/{cancer_type}/clustering/metric_df_{random_seed}.csv`
- `results/{cancer_type}/clustering/majority_vote_dict_{random_seed}.json`
- `results/{cancer_type}/clustering/penalty_sweep_{random_seed}.png`
- `results/{cancer_type}/clustering/find_best_k_plot_{random_seed}.tiff`
- `results/{cancer_type}/clustering/best_k_cluster_labels_{random_seed}.csv`

### CNV table

- `results/{cancer_type}/table/cnv_table.sqlite`
- `results/{cancer_type}/table/cnv_cluster_table_{random_seed}.sqlite`

### Selected features

- `mrmr`: `results/{cancer_type}/selected_features/top_features_cnv_mrmr{suffix}_{random_seed}.npy`
- `nvar`: `results/{cancer_type}/selected_features/top_features_cnv_nvar{suffix}_{random_seed}.npy`
- `rf`: `results/{cancer_type}/selected_features/top_features_cnv_rf{suffix}_{random_seed}.npy`

`mrmr` 固定選取 `20` 個 cluster-level features；`nvar` 會測試 2 至 20 個特徵。

### Random Forest feature importance

使用 `fs_method=rf` 時輸出：

- `results/{cancer_type}/importance_df/cnv_importance_{fs_method}{suffix}_{random_seed}.csv`

## 模型與評估

- 模型：`build_base_model("RandomForest", ...)`
- mRMR：`mrmr.mrmr_classif`
- 評估：`evaluate_model`
- clustering metric：silhouette、ARI 與 chromosome coherence ratio

## 注意事項

1. `clustering=false` 時，必須先存在 `metric_df_{random_seed}.csv` 與 `majority_vote_dict_{random_seed}.json` 等 clustering 結果。
2. `clustering=true` 可能需要數小時，實際時間取決於樣本數、基因數與 clustering 範圍。
3. 腳本以目前工作目錄解析 `./results`、`./data` 與 `./code/utils`，請從 `TLS_model` 根目錄執行。
4. `cancer_type` 需要與既有結果目錄及 `PivotTable` 的 `TLS` sample metadata 相符。
