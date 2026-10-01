# `ensemble_importance.py` 使用說明

本文件說明 `ensemble_importance.py` 如何整合 RNA、DNA methylation、mutation 與 CNV 的已選特徵，並計算多個分類模型的 feature importance。

## 腳本位置

- `code/feature_selection/ensemble_importance.py`

## 功能概述

此腳本位於各 omics 個別特徵選擇之後、consensus selection 之前。主要流程如下：

1. 讀取 JSON config。
2. 載入指定 cancer type 的 RNA、methylation、mutation 與 CNV `PivotTable`。
3. CNV 固定使用已完成 clustering 的 cluster-level table。
4. 依 `split_<random_seed>` 取出四種資料的 train/test samples。
5. 讀取四種 omics 個別產生的 `.npy` selected feature 檔。
6. 將四種 omics 的 selected features 串接成整合 train/test matrix。
7. 對每個指定模型先評估全特徵，再執行 100 次 feature importance。
8. 以 two-line elbow 決定整合模型的 top features，並進行最終評估。

## 執行方式

在 `TLS_model` 根目錄執行：

```bash
python code/feature_selection/ensemble_importance.py --config config_ensemble.json
```

也可以指定其他設定檔：

```bash
python code/feature_selection/ensemble_importance.py --config <config>.json
```

## Config

必要欄位：

- `random_seed`
- `cancer_type`
- `scaler_dict`

可選欄位：

- `model`: 模型名稱字串或模型名稱清單，例如 `["RF", "XGB", "CatBoost"]`。
- `fs_method`: 個別 omics selected feature 檔名中的方法字串，預設為空字串。
- `suffix`: 個別與整合 selected feature 檔名使用的 suffix，預設為空字串。
- `scaler_dict`: 各 omics 的 scaling 設定；空字典 `{}` 代表不額外指定 scaler。

目前預設設定檔 `config_ensemble.json` 使用：

```json
{
    "random_seed": 42,
    "cancer_type": "PAAD",
    "fs_method": "mrmr",
    "model": ["RF", "XGB", "CatBoost"],
    "scaler_dict": {},
    "suffix": ""
}
```

## 前置條件

執行前必須先完成四種 omics 的個別特徵選擇，並準備下列檔案：

- `results/{cancer_type}/table/rna_table.sqlite`
- `results/{cancer_type}/table/meth_table.sqlite`
- `results/{cancer_type}/table/mutation_table.sqlite`
- `results/{cancer_type}/table/cnv_cluster_table_{random_seed}.sqlite`
- `results/{cancer_type}/selected_features/top_features_rna_{fs_method}{suffix}_{random_seed}.npy`
- `results/{cancer_type}/selected_features/top_features_meth_{fs_method}{suffix}_{random_seed}.npy`
- `results/{cancer_type}/selected_features/top_features_mut_{fs_method}{suffix}_{random_seed}.npy`
- `results/{cancer_type}/selected_features/top_features_cnv_{fs_method}{suffix}_{random_seed}.npy`

四個 SQLite table 的 `sample_metadata` 都必須包含同一個 `split_{random_seed}` 欄位，且至少有 `train` 與 `test` 標籤。

## 資料整合方式

四個 selected feature table 會依下列順序串接：

```text
RNA + Meth + Mutation + CNV
```

`concat_omics` 使用的 omics keys 為：`RNA`、`Meth`、`Mutation`、`CNV`。整合後的 `feature_metadata` 會保留 omics 資訊，供 `scaler_dict` 依 omics 套用 scaling。

## 模型與 feature importance

- 模型由 `build_base_model(model_entry, ...)` 建立。
- `model` 可以是單一字串，也可以是模型清單。
- 每個模型固定執行 `100` 次 `run_feature_importance`。
- 使用 `two_line_elbow` 從平均 feature importance 排序結果決定 top features。
- 最終評估由 `evaluate_model` 完成，並輸出 train/test 分類指標與最佳 cutoff。

## 輸出結果

所有 `{cancer_type}` 路徑均相對於 `TLS_model` 根目錄。對每個模型名稱 `{model_name}`，其中 `{model_name}` 會轉成小寫：

### Feature importance

```text
results/{cancer_type}/importance_df/concat_importance_{fs_method}{suffix}_{model_name}_{random_seed}.csv
```

CSV 內容包含各次 feature selection 結果、selected count 與 selection frequency 等資訊。

### 整合後 selected features

```text
results/{cancer_type}/selected_features/top_features_concat_{fs_method}{suffix}_{model_name}_{random_seed}.npy
```

此 `.npy` 檔會被後續 `consensus_selection.py` 讀取，用於比較不同模型選出的 features。

### 終端輸出

執行時會印出：

- 每個模型的全特徵評估結果
- feature importance 計算進度
- 選出的 top features
- 每個模型使用 top features 後的 train/test 評估結果

## 與後續 consensus selection 的銜接

當 RF、XGB 與 CatBoost 都完成後，可將本腳本輸出的三組 `top_features_concat_*.npy` 提供給 `consensus_selection.py`，再選出至少被兩種模型共同選中的 TLS signatures。

## 注意事項

1. 本腳本目前在 `load_omics_data` 呼叫時固定使用 `cluster=True`，因此一定需要 `cnv_cluster_table_{random_seed}.sqlite`。
2. 腳本不會重新建立各 omics 的 train/test split；它只讀取既有的 `split_{random_seed}` 欄位。
3. `fs_method` 與 `suffix` 必須和四種個別 mRMR 腳本實際輸出的檔名一致，否則會找不到 `.npy` 檔。
4. 腳本以目前工作目錄解析 `./results` 與 `./code/utils`，請從 `TLS_model` 根目錄執行。
5. 若 `model` 為清單，所有模型都會依序執行，計算時間會隨模型數量與特徵數增加。
