# final_model.py 使用說明

本文件說明 `final_model.py` 的用途、設定檔參數、輸入資料與輸出結果。

## 腳本位置

- `final_model.py`

## 功能概述

`final_model.py` 用於建立最終整合模型，流程如下：

1. 讀取 `config_integrated.json`。
2. 載入 RNA、Meth、Mutation 與 CNV 四種 omics 資料。
3. 依 `split_<random_seed>` 切 train / test。
4. 載入已選好的 TLS signatures，作為最終模型使用的特徵。
5. 依 `model` 參數進行模型 tuning。
6. 以 selected TLS signatures 訓練與評估最終模型。
7. 對其他 cancer 的 independent dataset 做外部驗證。
8. 輸出 ROC、PR、Confusion Matrix 圖。

## 執行方式

在 `TLS_model` 根目錄執行：

```bash
python final_model.py --config config_integrated.json
```

如果要改成其他設定檔，也可以指定：

```bash
python final_model.py --config <your_config>.json
```

## 必要 config 欄位

腳本會檢查以下欄位：

- `random_seed`
- `cancer_type`
- `scaler_dict`

另外，`model` 不能是空值，因為程式會用它來做 tuning。

## 常用 config 欄位

- `model`
  - 可指定要 tuning 的模型名稱或模型物件
  - 程式內支援的模型名稱包含：`RandomForest`、`XGBoost`、`CatBoost`
- `random_seed`
  - 預設 `42`
- `cancer_type`
  - 預設 `PAAD`
  - 常見值：`PAAD`、`CHOL`、`STAD`、`all`
- `scaler_dict`
  - 傳給 `evaluate_model(...)` 與 `run_feature_importance(...)` 用的標準化/縮放設定

## 輸入資料來源

### 1. 四種 omics table

依 cancer type 載入：

- RNA：`results/{cancer_type}/table/rna_table.sqlite`
- Meth：`results/{cancer_type}/table/meth_table.sqlite`
- Mutation：`results/{cancer_type}/table/mutation_table.sqlite`
- CNV：
  - 目標 cancer：`results/{cancer_type}/table/cnv_cluster_table_{random_seed}.sqlite`
  - independent dataset：`results/{cancer_type}/table/cnv_table.sqlite`

### 2. TLS signatures

最終模型實際使用的特徵會從下列檔案載入：

- `results/{cancer_type}/selected_features/TLS_signatures_{random_seed}.npy`

### 3. CNV clustering 標籤

外部驗證時，CNV 需要先套用 clustering label：

- `results/{cancer_type}/clustering/best_k_cluster_labels_{random_seed}.csv`

## 模型流程細節

### 1. train / test split

腳本依 `split_{random_seed}` 將樣本切成 train / test，並保留原本順序。

### 2. 特徵整合

先將 RNA / Meth / Mutation / CNV concat 成整合矩陣，再只保留 `TLS_signatures_{random_seed}.npy` 中的特徵。

### 3. 模型 tuning

程式會呼叫 `utils/model_utils.py` 裡的 `tune_model(...)` 自動搜尋最佳參數，回傳：

- `best_model`
- `best_params`
- `best_score`

### 4. 最終評估

- 目標 cancer：使用 `concat_train` / `concat_test` 評估
- independent datasets：逐一評估 CHOL / PAAD / STAD 的外部資料

### 5. 圖表輸出

會輸出：

- ROC curve
- PR curve
- Confusion matrix

## 輸出結果

### 1. 特徵重要度

- `results/{cancer_type}/importance_df/final_model_importance_{random_seed}.csv`

### 2. 圖檔

會輸出到：

- `results/{cancer_type}/fig/ROC_curve_{cancer_type}_{random_seed}.tiff`
- `results/{cancer_type}/fig/PR_curve_{cancer_type}_{random_seed}.tiff`
- `results/{cancer_type}/fig/confusion_matrix_{cancer_type}_{random_seed}.tiff`

## 相依套件

- `pymaftools`
- `numpy`
- `pandas`
- `scikit-learn`
- `matplotlib`
- `utils/model_utils.py` 內的輔助函式

## 注意事項

1. `final_model.py` 依賴已完成的前處理與特徵選擇結果，包含 TLS signatures 與 CNV cluster table。
2. 若 `results/{cancer_type}/selected_features/TLS_signatures_{random_seed}.npy` 不存在，腳本會無法往下執行。
3. `process_clustering()` 會讀取 `best_k_cluster_labels_{random_seed}.csv`，因此 CNV clustering 必須先完成。
4. 腳本使用多個絕對路徑（`/home/data/data_YM313112004/TLS_model/...`），若搬移專案需要同步調整。
5. `cancer_type="all"` 時，主資料會使用 CHOL / PAAD / STAD merge 後的 table。
