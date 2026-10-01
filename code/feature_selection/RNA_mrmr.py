#!/usr/bin/env python
# coding: utf-8

from pymaftools import *
import json
import sys
import os
sys.path.append("./code/utils/")

def load_config(config_path="config_rna.json"):
    with open(config_path, "r") as f:
        return json.load(f)

def run_rna_model(config):
    # 檢查必要參數
    required_keys = ["input_table", "random_seed", "cancer_type"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")

    # 讀取參數
    input_table = config["input_table"]
    random_seed = config.get("random_seed", 42)
    cancer_type = config.get("cancer_type", "PAAD")
    fs_method = config.get("fs_method", "mrmr")
    fs_cutoff = config.get("fs_cutoff", 0.25)
    suffix = config.get("suffix", "")
    os.makedirs(f"./results/{cancer_type}/table", exist_ok=True)
    os.makedirs(f"./results/{cancer_type}/selected_features", exist_ok=True)
    print(f"\n##### Running model with random seed: {random_seed} and cancer type: {cancer_type} #####")

    # Load Data
    from pymaftools import PivotTable
    table_path = f"./results/{cancer_type}/table/rna_table.sqlite"
    if os.path.exists(table_path):
        print(f"Loading existing table.")
        table = PivotTable.read_sqlite(table_path)
        
    elif cancer_type == "all":
        chol_table = PivotTable.read_sqlite(f"./results/CHOL/table/rna_table.sqlite")
        paad_table = PivotTable.read_sqlite(f"./results/PAAD/table/rna_table.sqlite")
        stad_table = PivotTable.read_sqlite(f"./results/STAD/table/rna_table.sqlite")
        table = PivotTable.merge([chol_table, paad_table, stad_table])
        table.to_sqlite(db_path=table_path)

    else:
        all_table = PivotTable.read_sqlite(f"./results/{cancer_type}/table/{input_table}")
        common_samples = np.load(f"./results/{cancer_type}/common_samples.npy", allow_pickle=True)
        table = all_table.subset(samples=common_samples)

        # ---- Feature Filtering ----
        # Protein coding genes
        table = table.subset(features=table.feature_metadata["gene_type"] == "protein_coding")

        # Filter pseudoautosomal region genes
        non_par_mask = ~table.feature_metadata["ensg_name"].str.contains("PAR_Y", na=False)
        table = table.subset(features=non_par_mask)

        # Add feature metadata
        gencode = pd.read_csv("./data/gencode_anno_v36.csv", index_col=False)
        gencode = gencode[~gencode["gene_id"].str.contains("PAR_Y", na=False)]  # 移除包含 "PAR_Y" 的行
        gencode.index = gencode["gene_name"]
        merged_gencode = pd.merge(table.feature_metadata, gencode, 
                                left_on="ensg_name", right_on="gene_id", how="left") # 使用左合併，保留 table.feature_metadata 的所有行
        merged_gencode.index = merged_gencode["gene_name_x"]
        table = table.add_feature_metadata(merged_gencode, fill_value=np.nan, force=True)

        # Filter sex chromosome genes
        clean_table = pd.Series(table.feature_metadata["chr"]).dropna()
        mask = clean_table.apply(lambda x: x not in ["chrX", "chrY"])
        non_sex_chr = mask[mask].index.tolist()
        table = table.subset(features=non_sex_chr)

    # ---- Train/Test Split ---- 
    from sklearn.model_selection import train_test_split
    X = table.T  # shape: samples × genes
    y = table.sample_metadata["TLS"]  # TLS label (0 or 1)

    split_col = f"split_{random_seed}"
    if split_col not in table.sample_metadata.columns:
        # Add train/test labels to sample_metadata and sort
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_seed, stratify=y)
        table.sample_metadata.loc[y_train.index, split_col] = "train"
        table.sample_metadata.loc[y_test.index, split_col] = "test"
        # 只在第一次存檔時 reindex（以第一次的 seed 順序為主）
        if not os.path.exists(table_path):
            table = table.reindex(columns=y_train.index.tolist() + y_test.index.tolist())
        table.to_sqlite(db_path=table_path)
        print(f"Created '{split_col}' and saved PivotTable.")
    else:
        print(f"Using existing {split_col}")

    # Keep sample order as train then test for the current seed
    train_idx = table.sample_metadata.index[table.sample_metadata[split_col] == "train"].tolist()
    test_idx = table.sample_metadata.index[table.sample_metadata[split_col] == "test"].tolist()
    table = table.reindex(columns=train_idx + test_idx)
    print(f"Shape of PivotTable after feature filtering: {table.shape}")

    train_samples = table.sample_metadata[split_col] == "train" 
    train_table = table.subset(samples=train_samples)
    test_samples = table.sample_metadata[split_col] == "test" 
    test_table = table.subset(samples=test_samples)

    # Define model
    from model_utils import build_base_model
    clf = build_base_model("RandomForest", train_table, random_state=42)

    # All features performance
    all_features = train_table.index
    from model_utils import evaluate_model
    report_train, report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff = evaluate_model(train_table, test_table, top_features=all_features, model=clf)
    print(report)

    # ---- Feature Selection Method ----
    # Top N Variance
    from sklearn.model_selection import StratifiedKFold
    from mrmr import mrmr_classif
    from sklearn.metrics import roc_auc_score

    if fs_method == "nvar":
        print(f"Feature selection method: Top N Variance / MAD / Prevalence")
        X = pd.DataFrame(train_table.T)
        y = pd.Series(train_table.sample_metadata.loc[train_table.columns]["TLS"])
        top_n = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]  # 你要選的特徵數
        outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        auc_scores = []
        std_auc_scores = []
        best_auc = 0
        best_k = None

        for K in top_n:
            auc_scores_K = []
            for train_idx, val_idx in outer_cv.split(X, y):
                X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
                y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]
                
                # 只用 training fold 做 top n variance
                row_var = X_train.T.var(axis=1)
                threshold = row_var.nlargest(K).min()
                selected_features = row_var[row_var >= threshold].index
                X_train_sel = X_train[selected_features]
                X_val_sel = X_val[selected_features]
                
                clf.fit(X_train_sel, y_train)
                auc = roc_auc_score(y_val, clf.predict_proba(X_val_sel)[:,1])
                auc_scores_K.append(auc)
            
            mean_auc = np.mean(auc_scores_K)
            std_auc = np.std(auc_scores_K)
            auc_scores.append(mean_auc)
            std_auc_scores.append(std_auc)
            print(f"K={K}, threshold={threshold:.3f}, mean AUC={mean_auc:.3f}, std={std_auc:.3f}")
            
            if mean_auc > best_auc:
                best_auc = mean_auc
                best_k = K

        print(f"Best K: {best_k}, AUC: {best_auc:.3f}")

        row_var = train_table.var(axis=1)
        threshold = row_var.nlargest(best_k).min()
        top_features = row_var[row_var >= threshold].index # 取前百分比的 features
        top_features_path = f"./results/{cancer_type}/selected_features/top_features_rna_{fs_method}{suffix}_{random_seed}.npy"
        np.save(top_features_path, top_features)
        print(f"Top features: {top_features}")
        print(f"Top features saved to: {top_features_path}")
        train_table = train_table.subset(features=top_features)
        test_table = test_table.subset(features=top_features)
    
    # Minimum Redundancy Maximum Relevance (mRMR)
    elif fs_method == "mrmr":
        print(f"Feature selection method: Minimum Redundancy Maximum Relevance (mRMR) to select K features (fs_cutoff={fs_cutoff})")
        # 保留在 training validation set 中 high variance 的 features (fs_cutoff default = 0.25)
        row_var = train_table.var(axis=1)
        selected_features = row_var[row_var >= fs_cutoff].index 
        train_table = train_table.subset(features=selected_features)
        test_table = test_table.subset(features=selected_features)
        selected_features_path = f"./results/{cancer_type}/selected_features/top_features_rna_preselect{suffix}_{random_seed}.npy"
        np.save(selected_features_path, selected_features)
        print(f"Selected features after variance filtering: {len(selected_features)}")
        print(f"Selected features saved to: {selected_features_path}")

        X = pd.DataFrame(train_table.T)
        y = pd.Series(train_table.sample_metadata.loc[train_table.columns]["TLS"])

        mrmr_k = 150
        top_features = mrmr_classif(X=X, y=y, K=mrmr_k)  # 取前百分比的 features
        top_features_path = f"./results/{cancer_type}/selected_features/top_features_rna_{fs_method}{suffix}_{random_seed}.npy"
        np.save(top_features_path, top_features)
        print(f"Top features: {top_features}")
        print(f"Top features saved to: {top_features_path}")
        train_table = train_table.subset(features=top_features)
        test_table = test_table.subset(features=top_features)
    
    # Random Forest
    elif fs_method == "rf":
        print(f"Feature selection method: Random Forest (fs_cutoff={fs_cutoff})")
        # 保留在 training validation set 中 high variance 的 features (fs_cutoff default = 0.25)
        row_var = train_table.var(axis=1)
        top_features = row_var[row_var >= fs_cutoff].index
        train_table = train_table.subset(features=top_features)
        test_table = test_table.subset(features=top_features)
        
        # Feature Importance (Run 100 times)
        from model_utils import run_feature_importance
        os.makedirs(f"./results/{cancer_type}/importance_df", exist_ok=True)
        feature_importance_path = f"./results/{cancer_type}/importance_df/rna_importance_{fs_method}{suffix}_{random_seed}.csv"
        global_importance = run_feature_importance(train_table, run_seeds=100, model=clf, save_path=feature_importance_path, plot_top=20)
        print(f"Feature importance saved to: {feature_importance_path}")

        # Choose cutoff for selected features
        from model_utils import two_line_elbow
        knee_x, rss = two_line_elbow(global_importance)
        top_features = global_importance.head(knee_x).index
        top_features_path = f"./results/{cancer_type}/selected_features/top_features_rna_{fs_method}{suffix}_{random_seed}.npy"
        np.save(top_features_path, top_features)
        print(f"Top features: {top_features}")
        print(f"Top features saved to: {top_features_path}")

    else:
        raise ValueError(f"Unknown feature selection method: {fs_method}")

    # ---- Model Evaluation ----
    from model_utils import evaluate_model
    report_train, report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff = evaluate_model(train_table, test_table, top_features=top_features, model=clf)
    print(report)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_rna.json", help="Path to config JSON")
    args = parser.parse_args()
    config = load_config(args.config)
    run_rna_model(config)

