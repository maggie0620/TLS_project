import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, RandomizedSearchCV, StratifiedKFold, ParameterGrid, GridSearchCV
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from catboost import CatBoostClassifier
from tqdm import tqdm
from kneefinder import KneeFinder
from sklearn.metrics import roc_curve, f1_score, accuracy_score, roc_auc_score, average_precision_score, precision_recall_curve, auc
from sklearn.base import clone
from sklearn.inspection import permutation_importance
import warnings
from collections import Counter
import shap
from sklearn.linear_model import LinearRegression, LogisticRegressionCV, LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix

def sweep_penalty_select_best_k(metric_df, penalties_range=None, points=1000, method='multiplicative'):
    """
    1) 先做 coarse sweep（預設 logspace(-6,0,100)）找出 penalty 使得 selected k 變化的區間。
    2) 在該區間做 fine sweep（fine_points 個點），回傳細掃結果與建議的區間/k。

    Parameters
    ----------
    metric_df : DataFrame (包含 fold1..fold5_silhouette 欄位)
    penalties_range : array-like or None (若 None 使用 np.logspace(-6,0,100))
    points : int 精掃點數, default=1000
    method : 決定懲罰形式, default='multiplicative'
        - "additive": mean_silhouette - p * k
        - "multiplicative": mean_silhouette * (1 - p * k)

    Returns
    ----------
    df : 更新後的 metric_df
    fine_pens : array
        精掃的 penalty array
    fine_best_k : dict
        key=penalty, value=best_k
    most_common_k : int
        最穩定/常見的 k
    freq_dict : dict
        key=k, value=出現次數
    chosen_penalty : float
        在 most_common_k 的區間內，計算 penalized_silhouette 的 penalty
    """

    df = metric_df.copy()
    fold_cols = [c for c in df.columns if c.startswith("fold")]
    if len(fold_cols) == 0:
        raise ValueError("metric_df must contain fold* columns (e.g. fold1_silhouette)")
    df["mean_silhouette"] = df[fold_cols].mean(axis=1)

    # 在區間進行 penalty sweep
    if penalties_range is None:
        penalties_range = np.logspace(-6, 0, 100)  # 在 penalty 數值小的範圍有較多細節的測試點，在 penalty 大則不用浪費太多點
    coarse_best_k = {}
    ks_list = []
    pens_sorted = np.array(sorted(penalties_range))
    for p in pens_sorted:
        if method == 'additive':
            scores = df["mean_silhouette"] - p * df.index
        else:
            scores = df["mean_silhouette"] * (1 - p * df.index)
        best_k = int(scores.idxmax())
        coarse_best_k[p] = best_k
        ks_list.append(best_k)

    ks = np.array(ks_list)

    # 找出發生變化的區段
    diff_idx = np.where(np.diff(ks) != 0)[0]
    if len(diff_idx) > 0:
        i0 = diff_idx[0]
        i1 = diff_idx[-1] + 1  # inclusive
        p_min = pens_sorted[max(i0 - 1, 0)]
        p_max = pens_sorted[min(i1 + 1, len(pens_sorted)-1)]
    else:
        # 若沒有變化，擴展到整個 coarse 範圍的內側  (避免用到 extremes)
        p_min, p_max = pens_sorted[0], pens_sorted[-1]

    # 若 p_min == p_max（幾乎沒變化），擴展小範圍避免 logspace error
    if np.isclose(p_min, p_max):
        p_min = max(p_min / 10.0, pens_sorted[0])
        p_max = min(p_max * 10.0, pens_sorted[-1])
    
    print(f"Trimmed penalty range: [{p_min:.2e}, {p_max:.2e}]")

    # 再次進行 penalty sweep in trimmed range
    fine_pens = np.logspace(np.log10(p_min), np.log10(p_max), points)
    fine_best_k = {}
    for p in fine_pens:
        if method == 'additive':
            scores = df["mean_silhouette"] - p * df.index
        else:
            scores = df["mean_silhouette"] * (1 - p * df.index)
        fine_best_k[p] = int(scores.idxmax())

    # 選擇最穩定/常見的 k
    freq = Counter(fine_best_k.values())
    most_common_k, count = freq.most_common(1)[0]
    print(f"Most stable best_k = {most_common_k}, appearing {count}/{len(fine_best_k)} times")

    # 找出 most_common_k 的區域
    ranges = []
    start = None
    prev_p = None
    for p, k in zip(fine_pens, list(fine_best_k.values())):
        if k == most_common_k and start is None:
            start = p
        if k != most_common_k and start is not None:
            ranges.append((start, prev_p))
            start = None
        prev_p = p
    if start is not None:
        ranges.append((start, prev_p))

    # 利用 range 中的 chosen_penalty 計算 penalized_silhouette
    longest = max(ranges, key=lambda x: np.log10(x[1]) - np.log10(x[0]))   # 找出第一段或最長的穩定區段
    chosen_penalty = 10 ** np.mean([np.log10(longest[0]), np.log10(longest[1])])
    # if method == 'additive':
    #     df["penalized_silhouette"] = df["mean_silhouette"] - chosen_penalty * df.index
    # else:
    #     df["penalized_silhouette"] = df["mean_silhouette"] * (1 - chosen_penalty * df.index)

    return df, fine_pens, fine_best_k, most_common_k, dict(freq), ranges, chosen_penalty


def plot_clustering_metrics_and_find_best_k(
        metric_df, 
        filename, 
        title=None,
        dpi=300, 
        target_col="penalized_silhouette",
        penalty=None,
        method='multiplicative',
        plot_coherence=False,
        bbox_inches='tight', 
        transparent=True, 
        format=None,
        **kwargs):

    # === 計算 5-fold mean ± std ===
    fold_columns = [c for c in metric_df.columns if c.startswith("fold")]
    fold_df = metric_df.loc[:, fold_columns]
    fold_means = fold_df.mean(axis=1)
    fold_stds = fold_df.std(axis=1)

    if method == 'additive':
        metric_df["penalized_silhouette"] = metric_df["mean_silhouette"] - penalty * metric_df.index
    elif method == 'multiplicative':
        metric_df["penalized_silhouette"] = metric_df["mean_silhouette"] * (1 - penalty * metric_df.index)
    else:
        metric_df["penalized_silhouette"] = metric_df["mean_silhouette"] * (1 - penalty * metric_df.index) * (metric_df["chr_coherence_ratio"] ** 0.5)

    best_k = metric_df[target_col].idxmax()

    # === Plot ===
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(metric_df.index, fold_means, label="Mean silhouette (5-fold)", color="#1f77b4", linestyle='-')
    ax.errorbar(metric_df.index, fold_means, yerr=fold_stds, fmt='none', capsize=2, color='#bbbbbb', alpha=0.6)
    ax.plot(metric_df.index, metric_df["mean_ari_5_fold"], label="Mean ARI (5-fold)", color="#2ca02c", linestyle='-.')
    ax.plot(metric_df.index, metric_df["penalized_silhouette"], label=f"Penalized silhouette, penalty={penalty:.2e}", color="#7e40b8", linestyle='-.')
    if plot_coherence:
        ax.plot(metric_df.index, metric_df["chr_coherence_ratio"], label="Coherence ratio", color="#cb46a1", linestyle='-.')

    ax.axvline(best_k, color='red', linestyle=':', alpha=0.7, label=f'Selected k={best_k}')
    ax.set_xlabel("Number of Clusters (k)", fontsize=12)
    ax.set_ylabel("Score", fontsize=12)
    if title:
        ax.set_title(title, fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    # === Save figure ===
    format = filename.split('.')[-1].lower()
    pil_kwargs = {"compression": "tiff_lzw"} if format == "tiff" else {}
    fig.savefig(
        filename,
        dpi=dpi,
        bbox_inches=bbox_inches,
        transparent=transparent,
        format=format,
        pil_kwargs=pil_kwargs,
        **kwargs
    )
    plt.close(fig)

    return best_k


def clone_and_set_seed(model, seed):
    """
    Clone estimator and set random_state or random_seed if supported.
    Returns cloned estimator.
    """
    if model is None:
        return None
    m = clone(model)
    try:
        params = m.get_params()
        # 優先用 set_params 修改
        if 'random_state' in params:
            m.set_params(random_state=int(seed))
        elif 'random_seed' in params:
            m.set_params(random_seed=int(seed))
        else:
            # fallback: 嘗試直接設定屬性
            if hasattr(m, 'random_seed'):
                setattr(m, 'random_seed', int(seed))
            else:
                setattr(m, 'random_state', int(seed))
    except Exception:
        # 若無法設定 random_state，直接回傳 clone 的 model（不影響原行為）
        pass
    return m


def build_base_model(model, concat_train, random_state=42):
    model = str(model).strip()
    n_samples = concat_train.shape[1]
    class_weight = sum(concat_train.sample_metadata["TLS"] == 0) / sum(concat_train.sample_metadata["TLS"] == 1)

    # ---- RandomForest ----
    if model in ["RandomForest", "RandomForestClassifier", "RF"]:
        if n_samples < 100:
            return RandomForestClassifier(
                n_estimators=100,
                max_depth=3,
                min_samples_split=8,            # 分裂內部節點所需的最小樣本數 (default=2)
                min_samples_leaf=5,              # 葉節點所需的最小樣本數 (default=1)
                class_weight="balanced",
                random_state=random_state,
                n_jobs=10
            )
        else:
            return RandomForestClassifier(
                n_estimators=300,
                max_depth=5,
                min_samples_split=10,            # 分裂內部節點所需的最小樣本數 (default=2)
                min_samples_leaf=5,              # 葉節點所需的最小樣本數 (default=1)
                class_weight="balanced",
                random_state=random_state,
                n_jobs=10
            )
    
    # ---- XGBoost ----
    if model in ["XGBoost", "XGB", "XGBClassifier"]:
        if n_samples < 100:
            return XGBClassifier(
                n_estimators=150,
                max_depth=3,
                learning_rate=0.05,
                subsample=0.8,
                colsample_bytree=0.8,
                reg_lambda=8,
                min_child_weight=3,
                random_state=random_state,
                scale_pos_weight=class_weight,
                objective="binary:logistic",     # 預設 binary:logistic 就是二分類的 sigmoid，輸出機率
                eval_metric='auc', 
                importance_type="gain",
                n_jobs=10
            )
        else:
            return XGBClassifier(
                n_estimators=400,
                max_depth=5,
                learning_rate=0.01,
                subsample=0.8,
                colsample_bytree=0.8,
                reg_lambda=3,
                min_child_weight=5,
                random_state=random_state,
                scale_pos_weight=class_weight,
                objective="binary:logistic",    # 預設 binary:logistic 就是二分類的 sigmoid，輸出機率
                eval_metric='auc', 
                importance_type="gain",
                n_jobs=10
            )
  
    # ---- CatBoost ----
    if model in ["CatBoost", "CatBoostClassifier"]:
        if n_samples < 100:
            return CatBoostClassifier(
                iterations=150,
                depth=3,
                learning_rate=0.01,
                l2_leaf_reg=10,
                min_data_in_leaf=8,
                random_seed=random_state,
                class_weights=[1, class_weight],
                eval_metric='AUC',
                verbose=0,
                thread_count=10
            )   
        else:
            return CatBoostClassifier(
                iterations=200,
                depth=5,
                learning_rate=0.01,
                l2_leaf_reg=10,
                min_data_in_leaf=10,
                random_seed=random_state,
                class_weights=[1, class_weight],
                eval_metric='AUC',
                verbose=0,
                thread_count=10
            )   
    raise ValueError(f"Unsupported model: {model}")
 

def tune_model(concat_train, model_name, random_state=42, n_splits=5):
    X = concat_train.T.values
    y = concat_train.sample_metadata.loc[concat_train.columns]["TLS"].values
    class_weight = sum(concat_train.sample_metadata["TLS"] == 0) / sum(concat_train.sample_metadata["TLS"] == 1)
    n_samples = concat_train.shape[1]
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    model_name = str(model_name).strip()

    # ---- RandomForest ----
    if model_name in ["RandomForest", "RandomForestClassifier", "RF"]:
        base_model = RandomForestClassifier(random_state=random_state, class_weight="balanced", n_jobs=10)
        if n_samples < 100:
            param_grid_rf = {
                "n_estimators": [100, 150, 200],
                "max_depth": [3, 4, 5],
                "min_samples_split": [8, 10, 12],  # 節點要有足夠樣本才能分裂
                "min_samples_leaf": [5, 8, 10],   # 葉節點最少樣本
                "max_features": ["sqrt", "log2"],  # 減少每次分裂看的 feature，增加多樣性
            }
        else:
            param_grid_rf = {
                "n_estimators": [200, 300, 400],
                "max_depth": [4, 5, 6, 7],
                "min_samples_split": [8, 10, 12],
                "min_samples_leaf": [5, 8, 10],
                "max_features": ["sqrt", "log2"]
            }

        search = GridSearchCV(
            estimator=base_model,
            param_grid=param_grid_rf,
            scoring="roc_auc",
            cv=cv,
            n_jobs=5,
            verbose=0
        )
        search.fit(X, y)
        return search.best_estimator_, search.best_params_, float(search.best_score_)

    # ---- XGBoost ----
    if model_name in ["XGBoost", "XGB", "XGBClassifier"]:
        base_model = XGBClassifier(
            random_state=random_state, 
            n_jobs=10,
            eval_metric="auc",
            objective="binary:logistic", # 預設 binary:logistic 就是二分類的 sigmoid，輸出機率
            scale_pos_weight=class_weight,
            importance_type="gain"
        )
        param_grid_xgb = {
            "n_estimators": [150, 200, 250, 300, 350, 400],
            "max_depth": [3, 4, 5],
            "learning_rate": [0.01, 0.03, 0.05],
            "subsample": [0.6, 0.7, 0.8],  # colsample_bytree=0.8, subsample=0.8
            "colsample_bytree": [0.6, 0.8],
            "reg_lambda": [3, 5, 8, 10],   # 強制 L2 正則化，壓縮權重
            "min_child_weight": [3, 5, 8, 10]   # 防止模型為極少數樣本建立節點
        }

        search = GridSearchCV(
            estimator=base_model,
            param_grid=param_grid_xgb,
            scoring="roc_auc",
            cv=cv,
            n_jobs=10,
            verbose=0
        )
        search.fit(X, y)
        return search.best_estimator_, search.best_params_, float(search.best_score_)

    # ---- CatBoost ----
    if model_name in ["CatBoost", "CatBoostClassifier"]:
        param_grid_cat = {
            "depth": [3, 4, 5, 6],
            "iterations": [100, 200, 300],
            "learning_rate": [0.01, 0.03, 0.05],
            "l2_leaf_reg": [10, 12, 15],
            "min_data_in_leaf": [8, 10, 12]
        }

        best_score = -np.inf
        best_params = None

        for params in ParameterGrid(param_grid_cat):
            fold_scores = []
            for tr_idx, va_idx in cv.split(X, y):
                clf = CatBoostClassifier(
                    **params,
                    class_weights=[1, class_weight],
                    eval_metric="AUC",
                    verbose=0,
                    random_seed=random_state,
                    thread_count=10
                )
                clf.fit(
                    X[tr_idx], y[tr_idx], 
                    eval_set=(X[va_idx], y[va_idx])  # 提供驗證集讓 CatBoost 自動使用 early stopping
                )
                va_prob = clf.predict_proba(X[va_idx])[:, 1]
                fold_scores.append(roc_auc_score(y[va_idx], va_prob))

            mean_auc = float(np.mean(fold_scores))
            if mean_auc > best_score:
                best_score = mean_auc
                best_params = params

        best_model = CatBoostClassifier(
            **best_params,
            class_weights=[1, class_weight],
            eval_metric="AUC",
            verbose=0,
            thread_count=10
        )
        return best_model, best_params, float(best_score)

    raise ValueError(f"Unsupported model for tuning: {model_name}")


def run_feature_importance(train_table, run_seeds=100, concat=False, scaler_dict=None, save_path=None, plot_top=20, model=None):
    """
    計算 model 的基因重要性
    
    Parameters
    ----------
    train_table : PivotTable
        單一或合併過後的多組學資料
    run_seeds : int, default=100
        外層 random split 次數
    concat : bool, default=False
        輸入的多組學資料是否為 concat 的資料
    scaler_dict : dict or None, default=None
        每個 omics 的 scaling 方法字典，例如：
        {
            "RNA": StandardScaler(),
            "METH": MinMaxScaler(),
            "CNV": None
        }
    save_path : str or None, default=None
        如果指定，會將 importance DataFrame 存成 CSV 檔
    plot_top : int or None, default=20
        畫出前 N 個基因的 heatmap；設為 None 則不畫圖
    model : sklearn estimator instance or None, default=None
        如果提供，會使用此模型來計算 feature importance；否則使用 RandomForestClassifier

    Returns
    -------
    global_importance : Series
        平均 feature importance，index=基因，value=importance
    """
    
    # 轉成 samples × genes，並取得 TLS labels
    X_train = train_table.T
    y_train = train_table.sample_metadata.loc[train_table.columns]["TLS"].values

    if concat:
        # 對每個 omics 資料進行特定的 scaling
        scaled_train_omics = {}

        for omics in train_table.feature_metadata["omics"].unique():
            omics_features = train_table.feature_metadata[train_table.feature_metadata["omics"] == omics].index
            scaler = scaler_dict.get(omics, None) if scaler_dict else None
            
            if scaler is not None:
                # Fit-transform on training data, transform on validation data
                scaled_train_omics[omics] = pd.DataFrame(
                    scaler.fit_transform(X_train[omics_features]),
                    index=X_train.index,
                    columns=omics_features
                )
            else:
                # No scaling applied
                scaled_train_omics[omics] = X_train[omics_features].copy()

        # 合併所有 omics 資料
        X_train_scaled = pd.concat(scaled_train_omics.values(), axis=1)

    else:
        X_train_scaled = X_train.copy()

    # 準備 DataFrame 儲存 importance
    importance_df = pd.DataFrame(index=train_table.index, columns=np.arange(run_seeds))

    for seed in tqdm(range(run_seeds), desc="Feature importance"):
        # 每次都 clone 並設定不同 seed（若使用者提供 model）
        if model is None:
            local_model = RandomForestClassifier(random_state=seed, class_weight="balanced")
        else:
            # clone user model 並嘗試設定 seed (使用 seed 保證多樣性)
            local_model = clone_and_set_seed(model, seed)
        local_model.fit(X_train_scaled, y_train)

        # 取得 Feature Importance
        if hasattr(local_model, "feature_importances_"):
            # Tree-based
            importance_df.iloc[:, seed] = local_model.feature_importances_
        elif hasattr(local_model, "coef_"):
            # Linear models (Logistic, Linear SVC, etc.)
            # coef_ shape 可能是 (1, n_features) 或 (n_classes, n_features)
            imp = np.abs(local_model.coef_)
            if imp.ndim > 1:
                imp = np.mean(imp, axis=0) # 若是多類別，取平均
            importance_df.iloc[:, seed] = imp.flatten()
        else:
            # Fallback: Permutation Importance (for KNN, MLP, RBF SVM)
            # 注意: 在 100 次迴圈中執行 permutation importance 會非常耗時
            # 這裡為了效能，設定 n_repeats=5
            if seed == 42:
                print(f"[Info] Using Permutation Importance for {type(local_model).__name__}. "
                      f"Computing once and replicating to save time.")
                try:
                    r = permutation_importance(local_model, X_train_scaled, y_train, scoring="roc_auc",
                                            n_repeats=10, random_state=42, n_jobs=5)
                    importance_df.iloc[:, seed] = r.importances_mean
                except Exception as e:
                    print(f"Warning: Failed to calculate importance for seed {seed}: {e}")
                    importance_df.iloc[:, seed] = 0

    # 對所有 100 seed 平均
    global_importance = importance_df.mean(axis=1).sort_values(ascending=False)

    # Export CSV & heatmap
    importance_df = importance_df.astype(np.float64)
    importance_df["mean_importance"] = importance_df.mean(axis=1)
    importance_df = importance_df.sort_values(by="mean_importance", ascending=False)
    importance_df = importance_df.drop(columns="mean_importance")

    if save_path is not None:
        importance_df.to_csv(save_path)

    # Plot heatmap
    if plot_top is not None:
        plt.figure(figsize=(20, 8))
        sns.heatmap(importance_df.head(plot_top), cmap="Blues")
        plt.title(f"Top {plot_top} Feature Importance (scaler={scaler_dict})")
        plt.show()

    return global_importance


def run_feature_importance_freq(train_table, run_seeds=100, concat=False, scaler_dict=None, save_path=None, plot_top=20, model=None):
    """
    計算 model 的基因重要性
    
    Parameters
    ----------
    train_table : PivotTable
        單一或合併過後的多組學資料
    run_seeds : int, default=100
        外層 random split 次數
    concat : bool, default=False
        輸入的多組學資料是否為 concat 的資料
    scaler_dict : dict or None, default=None
        每個 omics 的 scaling 方法字典，例如：
        {
            "RNA": StandardScaler(),
            "METH": MinMaxScaler(),
            "CNV": None
        }
    save_path : str or None, default=None
        如果指定，會將 importance DataFrame 存成 CSV 檔
    plot_top : int or None, default=20
        畫出前 N 個基因的 heatmap；設為 None 則不畫圖
    model : sklearn estimator instance or None, default=None
        如果提供，會使用此模型來計算 feature importance；否則使用 RandomForestClassifier

    Returns
    -------
    global_importance : Series
        每個基因被選取到的頻率（被選次數 / run_seeds），index=基因，value=frequency
    """
    
    # 轉成 samples × genes，並取得 TLS labels
    X_train = train_table.T
    y_train = train_table.sample_metadata.loc[train_table.columns]["TLS"].values

    if concat:
        # 對每個 omics 資料進行特定的 scaling
        scaled_train_omics = {}

        for omics in train_table.feature_metadata["omics"].unique():
            omics_features = train_table.feature_metadata[train_table.feature_metadata["omics"] == omics].index
            scaler = scaler_dict.get(omics, None) if scaler_dict else None
            
            if scaler is not None:
                # Fit-transform on training data, transform on validation data
                scaled_train_omics[omics] = pd.DataFrame(
                    scaler.fit_transform(X_train[omics_features]),
                    index=X_train.index,
                    columns=omics_features
                )
            else:
                # No scaling applied
                scaled_train_omics[omics] = X_train[omics_features].copy()

        # 合併所有 omics 資料
        X_train_scaled = pd.concat(scaled_train_omics.values(), axis=1)

    else:
        X_train_scaled = X_train.copy()

    # 準備 DataFrame 儲存每個 seed 是否被選中（0/1）
    importance_df = pd.DataFrame(0.0, index=train_table.index, columns=np.arange(run_seeds))

    for seed in tqdm(range(run_seeds), desc="Feature importance"):
        # 每次都 clone 並設定不同 seed（若使用者提供 model）
        if model is None:
            local_model = RandomForestClassifier(random_state=seed, class_weight="balanced")
        else:
            # clone user model 並嘗試設定 seed (使用 seed 保證多樣性)
            local_model = clone_and_set_seed(model, seed)
        local_model.fit(X_train_scaled, y_train)

        # 取得 Feature Importance，並用 two-line elbow 決定該 seed 的 top features
        if hasattr(local_model, "feature_importances_"):
            feat_imp = pd.Series(local_model.feature_importances_, index=X_train_scaled.columns).sort_values(ascending=False)
        elif hasattr(local_model, "coef_"):
            imp = np.abs(local_model.coef_)
            if imp.ndim > 1:
                imp = np.mean(imp, axis=0)
            feat_imp = pd.Series(imp.flatten(), index=X_train_scaled.columns).sort_values(ascending=False)
        else:
            warnings.warn(
                f"Model {type(local_model).__name__} does not provide feature importance/coefs. Skipped.",
                UserWarning,
            )
            continue

        if feat_imp.sum() == 0:
            continue

        knee_x, rss = two_line_elbow(feat_imp, plot=False)
        selected_features = feat_imp.head(knee_x).index
        importance_df.loc[selected_features, seed] = 1.0

    # 被選次數與頻率
    selection_count = importance_df.sum(axis=1)
    global_importance = (selection_count / float(run_seeds)).sort_values(ascending=False)

    # Export CSV & heatmap
    importance_df = importance_df.astype(np.float64)
    importance_df["selected_count"] = selection_count
    importance_df["selection_frequency"] = selection_count / float(run_seeds)
    importance_df = importance_df.sort_values(by="selection_frequency", ascending=False)

    if save_path is not None:
        importance_df.to_csv(save_path)

    # Plot heatmap
    if plot_top is not None:
        plt.figure(figsize=(20, 8))
        heatmap_df = importance_df.drop(columns=["selected_count", "selection_frequency"]).head(plot_top)
        sns.heatmap(heatmap_df, cmap="Blues")
        plt.title(f"Top {plot_top} Feature Selection Frequency (scaler={scaler_dict})")
        plt.show()

    return global_importance


def run_feature_importance_l1(train_table, run_seeds=1, concat=False, scaler_dict=None, save_path=None, plot_top=20):
    """
    計算 model 的基因重要性
    
    Parameters
    ----------
    train_table : PivotTable
        單一或合併過後的多組學資料
    run_seeds : int, default=100
        外層 random split 次數
    concat : bool, default=False
        輸入的多組學資料是否為 concat 的資料
    scaler_dict : dict or None, default=None
        每個 omics 的 scaling 方法字典，例如：
        {
            "RNA": StandardScaler(),
            "METH": MinMaxScaler(),
            "CNV": None
        }
    save_path : str or None, default=None
        如果指定，會將 importance DataFrame 存成 CSV 檔
    plot_top : int or None, default=20
        畫出前 N 個基因的 heatmap；設為 None 則不畫圖

    Returns
    -------
    global_importance : Series
        平均 feature importance，index=基因，value=importance
    """
    
    # 轉成 samples × genes，並取得 TLS labels
    X_train = train_table.T
    y_train = train_table.sample_metadata.loc[train_table.columns]["TLS"].values

    if concat:
        # 對每個 omics 資料進行特定的 scaling
        scaled_train_omics = {}

        for omics in train_table.feature_metadata["omics"].unique():
            omics_features = train_table.feature_metadata[train_table.feature_metadata["omics"] == omics].index
            scaler = scaler_dict.get(omics, None) if scaler_dict else None
            
            if scaler is not None:
                # Fit-transform on training data, transform on validation data
                scaled_train_omics[omics] = pd.DataFrame(
                    scaler.fit_transform(X_train[omics_features]),
                    index=X_train.index,
                    columns=omics_features
                )
            else:
                # No scaling applied
                scaled_train_omics[omics] = X_train[omics_features].copy()

        # 合併所有 omics 資料
        X_train_scaled = pd.concat(scaled_train_omics.values(), axis=1)

    else:
        X_train_scaled = X_train.copy()

    # 準備 DataFrame 儲存 importance
    importance_df = pd.DataFrame(index=train_table.index, columns=np.arange(run_seeds))

    # 每次都 clone 並設定不同 seed（若使用者提供 model）
    for seed in tqdm(range(run_seeds), desc="Feature importance"):
        from sklearn.linear_model import LogisticRegression
        local_model = LogisticRegression(l1_ratio=0.5, solver="saga", C=0.5, max_iter=5000, random_state=seed, class_weight="balanced")
        local_model.fit(X_train_scaled, y_train)
        importance_df.iloc[:, seed] = np.abs(local_model.coef_).flatten()
        #print(f"Best C for seed {seed}: {local_model.C_[0]:.4f}, best l1_ratio: {local_model.l1_ratio_[0]:.2f}")

    # 保留非零 feature
    global_importance = importance_df.mean(axis=1).sort_values(ascending=False)

    # Export CSV & heatmap
    importance_df = importance_df.astype(np.float64)
    importance_df["mean_importance"] = importance_df.mean(axis=1)
    importance_df = importance_df.sort_values(by="mean_importance", ascending=False)
    importance_df = importance_df.drop(columns="mean_importance")

    if save_path is not None:
        importance_df.to_csv(save_path)

    # Plot heatmap
    if plot_top is not None:
        plt.figure(figsize=(20, 8))
        sns.heatmap(importance_df.head(plot_top), cmap="Blues")
        plt.title(f"Top {plot_top} Feature Importance (scaler={scaler_dict})")
        plt.show()

    return global_importance


def two_line_elbow(global_importance, plot=True):
    """
    利用兩段線性回歸 (piecewise linear) 找 elbow point。
    
    Parameters
    ----------
    global_importance : Series
        基因平均重要性，index=基因，value=importance

    Returns
    -------
    best_k : int
        Elbow point (Cutoff feature number)
    rss_list : list
        Residual sum of squares for each candidate breakpoint
    """
    global_importance_filter = global_importance[global_importance != 0]
    y = np.array(global_importance_filter.values)
    x = np.arange(1, len(y) + 1).reshape(-1, 1)

    # possible breakpoints (avoid first/last point)
    # 這裡確保至少有 3 個點才能進行切分
    if len(x) < 3:
        best_k = len(x)
        rss_list = []
        print(f"Too few features ({len(x)}) for two-line elbow. Selecting all non-zero features: {best_k}")
        return best_k, rss_list

    # possible breakpoints (avoid first/last point)
    candidate_ks = range(2, len(x)-1)

    rss_list = []
    for k in candidate_ks:
        # left segment
        model_left = LinearRegression().fit(x[:k], y[:k])
        pred_left = model_left.predict(x[:k])
        rss_left = np.sum((y[:k] - pred_left)**2)

        # right segment
        model_right = LinearRegression().fit(x[k:], y[k:])
        pred_right = model_right.predict(x[k:])
        rss_right = np.sum((y[k:] - pred_right)**2)

        rss_list.append(rss_left + rss_right)

    best_k = candidate_ks[np.argmin(rss_list)]
    print("Elbow point (Feature top n):", best_k)

    # Plot
    if plot:
        plt.figure(figsize=(12,8))
        plt.plot(x, y, '.-', label="Feature Importance")
        plt.axvline(best_k, color='red', linestyle='--', label=f"Cutpoint = {best_k}")
        plt.legend()
        plt.title("Two-line (piecewise) regression for elbow detection")
        plt.xlabel("Feature rank")
        plt.ylabel("Mean abs importance")
        plt.show()

    return best_k, rss_list


def three_line_elbow(global_importance, plot=True):
    """
    利用三段線性回歸 (piecewise linear) 找兩個 elbow points。

    Parameters
    ----------
    global_importance : Series
        基因平均重要性，index=基因，value=importance
    plot : bool, default=True
        是否繪製 importance 曲線與兩個 cutpoints

    Returns
    -------
    best_k1 : int
        第一個 cutpoint (Feature top n 的第一段結束點)
    best_k2 : int
        第二個 cutpoint (Feature top n 的第二段結束點)
    rss_list : list
        每組候選 (k1, k2) 的 total RSS，格式為 [((k1, k2), rss), ...]
    """
    global_importance_filter = global_importance[global_importance != 0]
    y = np.array(global_importance_filter.values)
    x = np.arange(1, len(y) + 1).reshape(-1, 1)

    # 三段切分時，為避免不穩定擬合，每段至少保留 2 個點
    # x[:k1], x[k1:k2], x[k2:]
    if len(x) < 6:
        best_k1 = len(x)
        best_k2 = len(x)
        rss_list = []
        print(
            f"Too few features ({len(x)}) for three-line elbow. "
            f"Selecting all non-zero features: {len(x)}"
        )
        return best_k1, best_k2, rss_list

    rss_list = []
    best_k1 = None
    best_k2 = None
    best_rss = np.inf

    # k1: 第一段結束位置；k2: 第二段結束位置
    for k1 in range(2, len(x) - 3):
        for k2 in range(k1 + 2, len(x) - 1):
            # left segment: [0, k1)
            model_left = LinearRegression().fit(x[:k1], y[:k1])
            pred_left = model_left.predict(x[:k1])
            rss_left = np.sum((y[:k1] - pred_left) ** 2)

            # middle segment: [k1, k2)
            model_mid = LinearRegression().fit(x[k1:k2], y[k1:k2])
            pred_mid = model_mid.predict(x[k1:k2])
            rss_mid = np.sum((y[k1:k2] - pred_mid) ** 2)

            # right segment: [k2, end)
            model_right = LinearRegression().fit(x[k2:], y[k2:])
            pred_right = model_right.predict(x[k2:])
            rss_right = np.sum((y[k2:] - pred_right) ** 2)

            total_rss = rss_left + rss_mid + rss_right
            rss_list.append(((k1, k2), total_rss))

            if total_rss < best_rss:
                best_rss = total_rss
                best_k1 = k1
                best_k2 = k2

    print("Elbow points (Feature top n):", best_k1, best_k2)

    if plot:
        plt.figure(figsize=(12, 8))
        plt.plot(x, y, '.-', label="Feature Importance")
        plt.axvline(best_k1, color='red', linestyle='--', label=f"Cutpoint 1 = {best_k1}")
        plt.axvline(best_k2, color='orange', linestyle='--', label=f"Cutpoint 2 = {best_k2}")
        plt.legend()
        plt.title("Three-line (piecewise) regression for elbow detection")
        plt.xlabel("Feature rank")
        plt.ylabel("Mean abs importance")
        plt.show()

    return best_k1, best_k2, rss_list

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def compute_stability_stats(importance_df, top_k=50):
    """
    計算 stability-based statistics
    
    Parameters
    ----------
    importance_df : DataFrame
        rows = seeds, columns = features, values = importance
    top_k : int
        每個 seed 取前 top_k 當作「被選到」
    
    Returns
    -------
    stats_df : DataFrame
        index = feature
        columns = ['mean_importance', 'frequency']
    """
    # mean importance
    mean_imp = importance_df.mean(axis=0)

    # selection frequency
    selected = []
    for _, row in importance_df.iterrows():
        top_features = row.sort_values(ascending=False).head(top_k).index
        selected.append(top_features)

    selected = pd.Series(np.concatenate(selected))
    freq = selected.value_counts() / importance_df.shape[0]

    stats_df = pd.DataFrame({
        "mean_importance": mean_imp,
        "frequency": freq
    }).fillna(0)

    # sort by mean importance
    stats_df = stats_df.sort_values("mean_importance", ascending=False)

    return stats_df

def stability_elbow(stats_df, min_freq=0.6, plot=True):
    """
    利用 stability (selection frequency) 做 elbow-based cutpoint
    
    Parameters
    ----------
    stats_df : DataFrame
        index = feature
        columns = ['mean_importance', 'frequency']
        必須已依 mean_importance 排序
    min_freq : float
        最低穩定性門檻（例如 0.6 或 0.7）
    
    Returns
    -------
    best_k : int
        Stability-based cutpoint
    """
    # 只看穩定 feature
    stable_df = stats_df[stats_df["frequency"] >= min_freq]

    if stable_df.shape[0] < 3:
        print("Too few stable features, selecting all.")
        return stable_df.shape[0]

    y = stable_df["frequency"].values
    x = np.arange(1, len(y) + 1).reshape(-1, 1)

    candidate_ks = range(2, len(x)-1)
    rss_list = []

    for k in candidate_ks:
        # left
        left_model = LinearRegression().fit(x[:k], y[:k])
        left_pred = left_model.predict(x[:k])
        rss_left = np.sum((y[:k] - left_pred) ** 2)

        # right
        right_model = LinearRegression().fit(x[k:], y[k:])
        right_pred = right_model.predict(x[k:])
        rss_right = np.sum((y[k:] - right_pred) ** 2)

        rss_list.append(rss_left + rss_right)

    best_k = candidate_ks[np.argmin(rss_list)]

    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, ".-", label="Selection frequency")
        plt.axvline(best_k, color="red", linestyle="--", label=f"Cutpoint = {best_k}")
        plt.xlabel("Feature rank (by mean importance)")
        plt.ylabel("Selection frequency")
        plt.title("Stability-based elbow detection")
        plt.legend()
        plt.show()

    return best_k



def plot_kneefinder(global_importance):
    """
    利用 KneeFinder 找到 elbow point
    
    Parameters
    ----------
    global_importance : Series
        基因平均重要性，index=基因，value=importance

    Returns
    -------
    knee_x : int
        Elbow point (Cutoff feature number)
    knee_y : float
        Elbow point importance value
    """

    # 找 elbow 
    data_y = global_importance.values
    data_x = range(1, len(data_y) + 1)
    kf = KneeFinder(data_x, data_y)
    knee_x, knee_y = kf.find_knee()
    print("Elbow point (Feature top n):", knee_x)

    # Plot
    plt.figure(figsize=(12,8))
    plt.plot(data_x, data_y, marker='o', markersize=3, label="Feature importance")
    plt.axvline(knee_x, color='r', linestyle='--', label=f"Elbow = {knee_x}")
    plt.xlabel("Feature by rank")
    plt.ylabel("Mean abs importance")
    plt.title("Elbow method for feature selection")
    plt.legend()

    return knee_x, knee_y


def plot_metric_boxplot(metric_list, run_seeds, num_features, metric_name):
    """
    繪製 Metric 的 boxplot 並顯示統計資訊。
    """
    mean_metric = np.mean(metric_list)
    std_metric = np.std(metric_list)
    min_metric = np.min(metric_list)
    max_metric = np.max(metric_list)

    # Boxplot
    plt.boxplot(metric_list)
    plt.ylabel(f"{metric_name}")
    plt.xticks([])
    plt.title(f"{metric_name} over {run_seeds} runs (selected top {num_features} features)")
    plt.text(0.78, 0.15, f"Mean={mean_metric:.3f}, Std={std_metric:.3f}\n\nMin={min_metric:.3f}, Max={max_metric:.3f}",
             horizontalalignment='center', transform=plt.gca().transAxes)
    plt.show()

    return mean_metric, std_metric, min_metric, max_metric

def evaluate_model(train_table, test_table, top_features, concat=False, scaler_dict=None, random_state=42, model=None, cutoff=None):
    """
    使用 top_features 訓練模型，計算 F1 score、Accuracy、AUC 和 AUPRC。

    Parameters
    ----------
    train_table : PivotTable
        單一或合併過後的多組學資料，包含所有 omics 的特徵和樣本。
    test_table : PivotTable
        單一或合併過後的多組學資料，包含所有 omics 的特徵和樣本。
    top_features : list
        選擇的前 N 個重要特徵。
    concat : bool, default=False
        是否合併多組學資料。
    scaler_dict : dict or None, default=None
        每個 omics 的 scaling 方法字典，例如：
        {
            "RNA": StandardScaler(),
            "METH": MinMaxScaler(),
            "CNV": None
        }
    random_state : int, default=42
        隨機狀態。
    model : sklearn estimator instance or None, default=None
        如果提供，會使用此模型來訓練；否則使用 RandomForestClassifier。
    cutoff : float or None, default=None
        如果提供，會使用此 cutoff 來分類；否則根據 train set 的 AUROC 計算最佳 cutoff。

    Returns
    -------
    f1_train : float
        訓練集的 F1 score。
    acc_train : float
        訓練集的 Accuracy。
    auc_train : float
        訓練集的 AUC。
    auprc_train : float
        訓練集的 AUPRC。
    f1_test : float
        測試集的 F1 score。
    acc_test : float
        測試集的 Accuracy。
    auc_test : float
        測試集的 AUC。
    auprc_test : float
        測試集的 AUPRC。
    """

    # 轉成 samples × genes，並取得 TLS labels
    X_train = train_table.T
    y_train = train_table.sample_metadata.loc[train_table.columns]["TLS"].values

    X_test = test_table.T
    y_test = test_table.sample_metadata.loc[test_table.columns]["TLS"].values   

    if concat:
        # 對每個 omics 資料進行特定的 scaling
        scaled_train_omics = {}
        scaled_test_omics = {}

        for omics in train_table.feature_metadata["omics"].unique():
            omics_features = train_table.feature_metadata[train_table.feature_metadata["omics"] == omics].index

            # 只保留 test set 真的存在的 features
            test_omics_features = [f for f in omics_features if f in X_test.columns]
            if len(test_omics_features) == 0:
                warnings.warn(f"[WARNING] Omics '{omics}' has no overlapping features between train and test. Skipped.", UserWarning)
                continue

            scaler = scaler_dict.get(omics, None) if scaler_dict else None

            if scaler is not None:
                # Fit-transform on training data, transform on validation data
                scaled_train_omics[omics] = pd.DataFrame(
                    scaler.fit_transform(X_train[omics_features]),
                    index=X_train.index,
                    columns=omics_features
                )
                scaled_test_omics[omics] = pd.DataFrame(
                    scaler.transform(X_test[test_omics_features]),
                    index=X_test.index,
                    columns=test_omics_features
                )
            else:
                # No scaling applied
                scaled_train_omics[omics] = X_train[omics_features].copy()
                scaled_test_omics[omics] = X_test[test_omics_features].copy()

        # 合併所有 omics 資料
        X_train_scaled = pd.concat(scaled_train_omics.values(), axis=1)
        X_test_scaled = pd.concat(scaled_test_omics.values(), axis=1)

    else:
        X_train_scaled = X_train.copy()
        X_test_scaled = X_test.copy()

    # 再次檢查並過濾 top_features
    top_features = sorted(top_features)  # 強制排序
    missing_in_train = [f for f in top_features if f not in X_train_scaled.columns]
    if missing_in_train:
        warnings.warn(f"下列 top_features 在 train_table 中找不到，將被忽略：{missing_in_train}", UserWarning)

    selected_features_train = [f for f in top_features if f in X_train_scaled.columns]
    missing_in_test = [f for f in selected_features_train if f not in X_test_scaled.columns]
    if missing_in_test:
        warnings.warn(f"下列 top_features 在 test_table 中找不到，將被忽略：{missing_in_test}", UserWarning)

    common_features = [f for f in selected_features_train if f in X_test_scaled.columns]
    if not common_features:
        raise ValueError("沒有可用的 top_features（train 或 test 中皆未找到相符的 features）。")

    # 標註缺失的 features
    missing_features = [feature for feature in common_features if feature not in X_train_scaled.columns]
    if missing_features:
        print(f"Missing features: {missing_features}")

    # 準備訓練與測試資料 (samples × features)
    X_train_top = X_train_scaled[common_features].copy()
    X_test_top = X_test_scaled[common_features].copy()

    # 訓練 model
    if model is None:
        model = RandomForestClassifier(random_state=random_state, class_weight="balanced")
    model.fit(X_train_top, y_train)

    # 預測類別和機率
    #y_pred = model.predict(X_train_top)
    y_pred_prob = model.predict_proba(X_train_top)[:, 1]
    #y_test_pred = model.predict(X_test_top)
    y_test_pred_prob = model.predict_proba(X_test_top)[:, 1]

    # AUROC-based cutoff
    if cutoff is not None:
        best_cutoff = cutoff
    else:
        fpr, tpr, thresholds = roc_curve(y_train, y_pred_prob, drop_intermediate=False)  # default 會刪掉部分 ROC 上中間點，只保留較精簡的點
        youden = tpr - fpr
        valid = np.isfinite(thresholds)  # 排除 inf
        if valid.any():
            # 先用 Youden 最大值決定初始 cutoff，如果有多個最佳 cutoff，取最佳區間中心值： (min + max) / 2
            max_youden_idx = np.argmax(youden[valid])
            init_cutoff = float(thresholds[valid][max_youden_idx])

            # 依初始 cutoff 先產生標籤，再用「1 的最小分數」與「0 的最大分數」中點更新 cutoff
            y_pred_init = (y_pred_prob >= init_cutoff).astype(int)
            pos_scores = y_pred_prob[y_pred_init == 1]
            neg_scores = y_pred_prob[y_pred_init == 0]

            if len(pos_scores) > 0 and len(neg_scores) > 0:
                min_pos_score = float(np.min(pos_scores))
                max_neg_score = float(np.max(neg_scores))
                best_cutoff = (min_pos_score + max_neg_score) / 2.0
            elif len(pos_scores) > 0:
                best_cutoff = float(np.min(pos_scores))
            elif len(neg_scores) > 0:
                best_cutoff = float(np.max(neg_scores))
            else:
                best_cutoff = init_cutoff
        else:
            best_cutoff = 0.5

    y_pred = (y_pred_prob >= best_cutoff).astype(int)
    y_test_pred = (y_test_pred_prob >= best_cutoff).astype(int)

    # Model Performance
    print(f"Selected top {len(common_features)} features")
    print(f"Best cutoff in AUROC: {best_cutoff}")

    f1_train = f1_score(y_train, y_pred)
    acc_train = accuracy_score(y_train, y_pred)
    auc_train = roc_auc_score(y_train, y_pred_prob)
    auprc_train = average_precision_score(y_train, y_pred_prob)
    print(f"Train | F1-score: {f1_train:.3f} | Accuracy: {acc_train:.3f} | AUC: {auc_train:.3f} | AUPRC: {auprc_train:.3f}")

    f1_test = f1_score(y_test, y_test_pred)
    acc_test = accuracy_score(y_test, y_test_pred)
    auc_test = roc_auc_score(y_test, y_test_pred_prob)
    auprc_test = average_precision_score(y_test, y_test_pred_prob)
    print(f"Test  | F1-score: {f1_test:.3f} | Accuracy: {acc_test:.3f} | AUC: {auc_test:.3f} | AUPRC: {auprc_test:.3f}")

    report_train = classification_report(
    y_train,
    y_pred,
    target_names=['TLS-negative', 'TLS-positive'],
    digits=3,
    zero_division=0,
    )
    
    report = classification_report(
        y_test,
        y_test_pred,
        target_names=['TLS-negative', 'TLS-positive'],
        digits=3,
        zero_division=0
    )

    return report_train, report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff


from pymaftools import PivotTable 
from typing import List, Optional, Any, Literal

def concat_omics(
        tables: List["PivotTable"],
        keys: Optional[List[str]] = None,
        axis: int = 0,
        join: Literal["inner", "outer"] = "outer",
        fill_value: Any = np.nan,
        feature_fill_value: Any = np.nan,
        sample_fill_value: Any = np.nan,
        sep: str = "_"
    ) -> "PivotTable":
    """
    Concatenate PivotTables along axis 0 (features) or axis 1 (samples).

    Parameters
    ----------
    tables : List[PivotTable]
        List of PivotTable objects to concatenate.
    keys : List[str], optional
        Keys to use for constructing a MultiIndex or prefixed index.
        Only used when axis=0 to prefix feature names (e.g. RNA_Gene1).
    axis : {0, 1}, default 0
        The axis to concatenate along.
        0: Stack features (vertical). Useful for multi-omics integration.
        1: Stack samples (horizontal). Same as PivotTable.merge.
    join : {'inner', 'outer'}, default 'outer'
        How to handle indexes on other axis.
    fill_value : scalar, default np.nan
        Value for missing data.
    feature_fill_value : scalar, default np.nan
        Value to fill missing values in feature_metadata.
    sample_fill_value : scalar, default np.nan
        Value to fill missing values in sample_metadata.
    sep : str, default "_"
        Separator for joining key and index when axis=0 and keys are provided.

    Returns
    -------
    PivotTable
        Concatenated PivotTable with synchronized metadata.
    """
    # 如果是 axis=1 (水平合併樣本)，直接呼叫現有的 merge
    if axis == 1:
        return PivotTable.merge(
            tables, 
            fill_value=fill_value, 
            feature_fill_value=feature_fill_value,
            sample_fill_value=sample_fill_value,
            join=join
        )
        
    # Axis 0: Concatenate features (Vertical stacking)
    ## Concatenate Data
    if keys:
        if len(keys) != len(tables):
            raise ValueError("Length of keys must match length of tables")
        # Use keys to create MultiIndex then flatten
        concat_df = pd.concat(tables, axis=0, keys=keys, join=join)
        # Flatten MultiIndex: (Key, Feature) -> "Key_Feature"
        concat_df.index = [f"{k}{sep}{i}" for k, i in concat_df.index]
    else:
        concat_df = pd.concat(tables, axis=0, join=join)
        
    ## Feature Metadata
    feat_meta_list = []
    for i, tbl in enumerate(tables):
        meta = tbl.feature_metadata.copy()
        if keys:
            # Prefix index to match data
            meta.index = [f"{keys[i]}{sep}{idx}" for idx in meta.index]
        # Add a column to record the source omics
        meta["omics"] = keys[i] if keys else f"omics_{i}"
        feat_meta_list.append(meta)
        
    # Concatenate feature metadata
    merged_feature_meta = pd.concat(feat_meta_list, axis=0, join="outer")
    
    ## Sample Metadata
    sample_meta_list = [t.sample_metadata for t in tables]
    merged_sample_meta = pd.concat(sample_meta_list, axis=0)
    # Remove duplicates (keep first occurrence for each sample)
    merged_sample_meta = merged_sample_meta[~merged_sample_meta.index.duplicated(keep='first')]
    
    ## Output PivotTable
    result = PivotTable(concat_df.fillna(fill_value))
    # Align metadata to the new data shape
    result.feature_metadata = merged_feature_meta.reindex(result.index, fill_value=feature_fill_value)
    result.sample_metadata = merged_sample_meta.reindex(result.columns, fill_value=sample_fill_value)
    
    return result


def tune_final_model(X_train, top_features=None, model=None, param_grid=None, n_splits=5, scoring='roc_auc', n_jobs=10, random_state=42):
    """
    針對 Final Model (Concat 後) 進行超參數調優。
    這是在 Feature Selection 後，針對選定的特徵組合，尋找最佳的模型參數。
    
    Parameters
    ----------
    X_train : DataFrame (features x samples)
        包含 sample_metadata 的訓練資料。
    top_features : list or None
        若有指定，只使用這些特徵進行調參。
    model : estimator
        要調參的基底模型 (e.g. XGBClassifier, RandomForestClassifier)。
        若為 None，預設使用 RandomForestClassifier。
    param_grid : dict
        GridSearchCV 的參數網格。
    n_splits : int
        Cross-validation 的 fold 數。
    scoring : str
        評估指標 (e.g. 'roc_auc', 'f1', 'accuracy')。
    
    Returns
    -------
    best_model : estimator
        調參後的最佳模型 (已 fit 在整個 X_train 上)。
    best_params : dict
        最佳參數組合。
    """
    # 1. 準備資料
    X = X_train.T
    y = X_train.sample_metadata.loc[X_train.columns]["TLS"].values

    # 2. 特徵篩選
    if top_features is not None:
        # 過濾不存在的特徵
        valid_features = [f for f in top_features if f in X.columns]
        if len(valid_features) < len(top_features):
            print(f"Warning: {len(top_features) - len(valid_features)} features dropped (not in X_train).")
        X = X[valid_features]
    
    # 3. 設定模型與參數
    if model is None:
        # 預設使用 Random Forest
        model = RandomForestClassifier(random_state=random_state, class_weight="balanced")
    
    if param_grid is None:
        # 預設參數網格 (針對 RF)
        param_grid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [None, 10, 20],
            'min_samples_leaf': [1, 2, 4],
            'class_weight': ['balanced', 'balanced_subsample']
        }

    # 4. 執行 GridSearchCV
    print(f"Starting GridSearchCV (n_splits={n_splits}, scoring={scoring})...")
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    
    grid_search = GridSearchCV(estimator=model, param_grid=param_grid, 
                               cv=cv, scoring=scoring, n_jobs=n_jobs, verbose=1)
    
    grid_search.fit(X, y)
    
    print(f"Best Parameters: {grid_search.best_params_}")
    print(f"Best CV Score ({scoring}): {grid_search.best_score_:.4f}")
    
    return grid_search.best_estimator_, grid_search.best_params_

from sklearn.model_selection import KFold
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
def k_fold_clustering_evaluation_step(table: PivotTable,
                                 min_clusters: int = 2,
                                 max_clusters: int = 50,
                                 step: int = 1,
                                 metric: Literal['cosine', 'hamming', 'jaccard'] = "cosine",
                                 random_state: int = 42, 
                                 group_col="subtype"):
    """
    使用K-fold交叉驗證評估基因聚類的最佳聚類數

    Parameters
    ----------
    table : PivotTable
        基因表達或CNV數據表
    max_clusters : int
        最大聚類數
    metric : str
        相似性度量方法
    random_state : int
        隨機種子

    Returns
    -------
    pd.DataFrame
        包含每個fold和聚類數的silhouette分數
    """
    group = table.sample_metadata[group_col].values
    kf = KFold(n_splits=5, shuffle=True, random_state=random_state)
    if not isinstance(table, PivotTable):
        raise ValueError("Input table must be a PivotTable instance.")
    if min_clusters < 2 or max_clusters < min_clusters:
        raise ValueError(
            "min_clusters must be at least 2 and max_clusters must be greater than or equal to min_clusters.")
    all_results = []
    cluster_label_dict = {}
    # only use (k-1)/k train part for clustering
    for fold, (train_idx, _) in enumerate(kf.split(table.T.values, np.array(group))):
        print(f"Processing fold {fold + 1}/5")

        # 獲取訓練樣本
        sample_train = table.sample_metadata.iloc[train_idx].index
        table_train = table.subset(samples=sample_train)
        similarity_matrix = table_train.T.compute_similarity(method=metric)
        distance_matrix_train = 1 - similarity_matrix # dist = 1 - similarity
        # fill diagonal to 0
        for i in range(len(distance_matrix_train)):
            distance_matrix_train.iloc[i, i] = 0.0

        # test different number of clusters
        for k in tqdm(range(min_clusters, max_clusters + 1, step), desc=f"Fold {fold + 1}"):
            model = AgglomerativeClustering(
                n_clusters=k,
                metric='precomputed',
                linkage='average'
            )
            labels = model.fit_predict(distance_matrix_train)
            if k not in cluster_label_dict:
                cluster_label_dict[k] = {}
            cluster_label_dict[k][fold+1] = labels

            # 計算silhouette分數
            score = silhouette_score(table_train, labels, metric=metric)

            all_results.append({
                'fold': fold + 1,
                'n_clusters': k,
                'silhouette_score': score
            })

    return pd.DataFrame(all_results), cluster_label_dict


def compute_chr_coherence_ratio(labels, feature_cytoband):
    """
    Dominant chromosome coherence:
    For each cluster, compute the proportion of features
    belonging to the most frequent chromosome.
    Then average across clusters.
    """
    df = pd.DataFrame({"cluster": labels,
                       "cytoband": feature_cytoband.values}, index=feature_cytoband.index)

    cluster_scores = []

    for _, sub_df in df.groupby("cluster"):
        # extract chromosome (e.g. 1p36.33 -> 1, 17q21 -> 17)
        chrs = (sub_df["cytoband"].str.replace(r"(p|q).*", "", regex=True))

        if len(chrs) == 0:
            continue

        chr_counts = chrs.value_counts()
        dominant_ratio = chr_counts.iloc[0] / chr_counts.sum()
        cluster_scores.append(dominant_ratio)

    return float(np.mean(cluster_scores))
