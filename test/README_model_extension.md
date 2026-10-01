# feature_selection/model_extension.py 使用說明

本文件說明 `model_extension.py` 的用途、設定檔參數、資料需求與輸出結果。

## 腳本位置

- `feature_selection/model_extension.py`

## 功能概述

`model_extension.py` 用於建立整合型 TLS 分類模型，將 RNA、Meth、Mutation 與 CNV 四種資料整合後訓練主模型，並評估模型在目標 cancer 與其他 cancer 的表現。

主要流程如下：

1. 讀取 JSON config。
2. 依 cancer type 載入 RNA / Meth / Mutation / CNV 的 PivotTable。
3. 將每個 omics 資料依 `split_<random_seed>` 切成 train / test。
4. 讀取各模態已選好的特徵檔，組成整合特徵矩陣。
5. 建立 base model，並先做全特徵評估。
6. 以 Random Forest 做 feature importance，決定整合模型的 top features。
7. 用選出的 top features 再做一次最終評估。
8. 對其他 cancer 的 independent dataset 做外部驗證。

## 執行方式

在 `TLS_model` 根目錄執行：

```bash
python feature_selection/model_extension.py --config config_integrated.json
```

## 必要 config 欄位

腳本會檢查以下必要鍵值：

- `random_seed`
- `cancer_type`
- `scaler_dict`
- `save_feature_importance`
- `save_top_features`

## 常用可選欄位

- `model`: base model 名稱或可被 `eval` 的模型字串
- `fs_prefix`: 前綴字串，會插入特徵檔名稱中，預設空字串

## 輸入資料來源

### 1. 四種 omics table

依 cancer type 讀取：

- RNA：`results/{cancer_type}/table/rna_table.sqlite`
- Meth：`results/{cancer_type}/table/meth_table.sqlite`
- Mutation：`results/{cancer_type}/table/mutation_table.sqlite`
- CNV：
  - 目標 cancer：`results/{cancer_type}/table/cnv_cluster_table_{random_seed}.sqlite`
  - 非目標 cancer 的 independent set：`results/{cancer_type}/table/cnv_table.sqlite`

### 2. 目標樣本清單

- `selected_sample`

### 3. 各模態已選特徵檔

腳本會從下列路徑讀取 top features：

- `results/{cancer_type}/selected_features/top_features_rna_{fs_prefix}{random_seed}.npy`
- `results/{cancer_type}/selected_features/top_features_meth_{fs_prefix}{random_seed}.npy`
- `results/{cancer_type}/selected_features/top_features_mut_{fs_prefix}{random_seed}.npy`
- `results/{cancer_type}/selected_features/top_features_cnv_{fs_prefix}{random_seed}.npy`

## CNV clustering 依賴

`model_extension.py` 會直接使用已完成的 CNV clustering 結果：

- `results/{cancer_type}/clustering/best_k_cluster_labels_{random_seed}.csv`

這些 cluster labels 會被 `process_clustering()` 套用到 independent dataset 的 CNV table 上，並產生 cluster-level CNV table。

## 輸出結果

### 1. 特徵重要度

- `results/{cancer_type}/importance_df/{save_feature_importance}_{random_seed}.csv`

### 2. 選後整合特徵

- `results/{cancer_type}/selected_features/{save_top_features}_{random_seed}.npy`

### 3. 其他中間結果

腳本會在執行過程中印出：

- train/test 的全特徵報告
- top features
- target cancer 與 independent cancer 的評估報告

## 模型與評估

- Base model 由 `utils/model_utils.py` 的 `build_base_model(model, ...)` 建立。
- 評估由 `evaluate_model(...)` 完成。
- feature importance 由 `run_feature_importance(...)` 產生。
- 最佳截點由 `two_line_elbow(...)` 決定。

## 主要注意事項

1. `cancer_type="PAAD"`、`CHOL`、`STAD` 時，會把其他兩個 cancer 當作 independent datasets 做驗證。
2. `cancer_type="all"` 時，主資料會使用 CHOL / PAAD / STAD merge 後的 table。
3. `process_clustering()` 需要先有 `best_k_cluster_labels_{random_seed}.csv`，否則 independent dataset 的 CNV clustering 無法進行。
4. 本腳本使用多個絕對路徑（`/home/data/data_YM313112004/TLS_model/...`），若搬移專案需同步調整。
5. `model` 若傳入字串，程式會用 `eval` 轉成模型物件，請確保字串安全且可解析。
