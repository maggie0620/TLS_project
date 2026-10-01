# `DNAmut_mrmr.py` 使用說明

本文件說明 `DNAmut_mrmr.py` 的 mutation gene 前處理、mRMR 特徵選擇與 TLS 分類流程。

## 腳本位置

- `code/feature_selection/DNAmut_mrmr.py`

## 功能概述

腳本以 mutation gene `PivotTable` 建立 TLS 分類模型，流程如下：

1. 讀取 JSON config。
2. 載入 mutation table，並取出共同樣本。
3. 保留 protein-coding genes、補上 GENCODE annotation，並移除 chrX/chrY。
4. 建立或重用 `split_<random_seed>` 的 80/20 分層 train/test split。
5. 建立 Random Forest，先評估全特徵基準模型。
6. 以 mutation prevalence、mRMR 或 Random Forest 進行特徵選擇。
7. 使用選出的特徵重新評估 TLS 分類模型。

## 執行方式

在 `TLS_model` 根目錄執行：

```bash
python code/feature_selection/DNAmut_mrmr.py --config config_mut.json
```

也可以指定其他設定檔：

```bash
python code/feature_selection/DNAmut_mrmr.py --config <config>.json
```

## Config

必要欄位：`input_table`、`random_seed`、`cancer_type`。

可選欄位：

- `fs_method`: `mrmr`、`nvar` 或 `rf`，預設為 `mrmr`。
- `fs_cutoff`: mutation prevalence 前篩門檻，預設為 `0.01`。
- `suffix`: 加在 selected feature 檔名中的字串，預設為空字串。

預設設定檔為 `config_mut.json`。

## 輸入資料

當 `cancer_type` 不是 `all` 時，腳本使用：

- mutation table：`results/{cancer_type}/table/{input_table}`
- 共同樣本：`results/{cancer_type}/common_samples.npy`
- gene annotation：`data/gencode_anno_v36.csv`

mutation feature metadata 需要包含 `BIOTYPE` 與 `Gene` 欄位；補充 annotation 後會使用 `chr` 欄位移除性染色體基因。

當 `cancer_type` 為 `all` 時，會合併 CHOL、PAAD、STAD 的 `mutation_table.sqlite`，缺失 mutation value 以 0 補值。

## 特徵選擇

- `mrmr`: 先保留 training data 中 prevalence `>= fs_cutoff` 的 genes，再固定使用 `mrmr_k=10`。
- `nvar`: 實際使用 mutation prevalence，在 `[2, 4, ..., 20]` 中以 5-fold AUC 選擇最佳數量。
- `rf`: 先以 prevalence 過濾，再以 100 次 Random Forest feature importance 與 two-line elbow 選擇特徵。

## 輸出結果

所有 `{cancer_type}` 路徑均相對於 `TLS_model` 根目錄。

### Table 與 split

- `results/{cancer_type}/table/mutation_table.sqlite`
- 若不存在 `split_<random_seed>`，會將 train/test 標籤寫回上述 SQLite table。

### Selected features

- mRMR 前篩：`results/{cancer_type}/selected_features/top_features_mut_preselect{suffix}_{random_seed}.npy`
- mRMR 選後：`results/{cancer_type}/selected_features/top_features_mut_mrmr{suffix}_{random_seed}.npy`
- nvar：`results/{cancer_type}/selected_features/top_features_mut_nvar{suffix}_{random_seed}.npy`
- rf：`results/{cancer_type}/selected_features/top_features_mut_rf{suffix}_{random_seed}.npy`

### Random Forest feature importance

使用 `fs_method=rf` 時輸出：

- `results/{cancer_type}/importance_df/mut_importance_rf_{random_seed}.csv`

## 模型與評估

- 模型：`build_base_model("RandomForest", ...)`
- mRMR：`mrmr.mrmr_classif`
- 評估：`evaluate_model`，會輸出 train/test 的分類指標與最佳 cutoff。

## 注意事項

1. 腳本以目前工作目錄解析 `./results`、`./data` 與 `./code/utils`，請從 `TLS_model` 根目錄執行。
2. mutation prevalence 是 training samples 中非零值的比例；test table 不參與特徵選擇。
3. 若指定 `cancer_type=all`，需先準備 CHOL、PAAD、STAD 三個 cancer 的 mutation SQLite table。
