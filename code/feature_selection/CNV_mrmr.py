#!/usr/bin/env python
# coding: utf-8

from pymaftools import *
import json
import os
import sys
import os
sys.path.append("./code/utils/")

def load_config(config_path="config_cnv.json"):
    with open(config_path, "r") as f:
        return json.load(f)

def run_clustering(train_table, min_clusters=2, max_clusters=350, random_seed=42, cancer_type="PAAD"):
    # run k fold clustering evaluation (8.5HR)
    results_df, cluster_label_dict = k_fold_clustering_evaluation(
                                        train_table, 
                                        min_clusters=min_clusters,
                                        max_clusters=max_clusters, 
                                        group_col="TLS"
                                        )
    converted_data = convert_ndarray_to_list(cluster_label_dict)
    with open(f"./results/{cancer_type}/clustering/cluster_labels_{min_clusters}_{max_clusters}_{random_seed}.json", 'w', encoding='utf-8') as f:
        json.dump(converted_data, f)

    # 將 clustering 結果進行對齊（以 fold 1 為參考）
    with open(f"./results/{cancer_type}/clustering/cluster_labels_{min_clusters}_{max_clusters}_{random_seed}.json", 'r', encoding='utf-8') as f:
        cluster_label_dict = json.load(f)
    aligned_all_cluster_label_dict = align_cluster_label_dict(cluster_label_dict)

    # Calculate metrics (4.5HR)
    import pandas as pd
    from sklearn.metrics import silhouette_score, adjusted_rand_score
    from model_utils import compute_chr_coherence_ratio
    k_range = range(min_clusters, max_clusters + 1)
    feature_cytoband = (train_table.feature_metadata.loc[train_table.index, "Cytoband"])
    metric_df = pd.DataFrame(index=k_range, columns=["majority_vote_silhouette", 
                                                     "mean_ari_5_fold",
                                                     "chr_coherence_ratio"])
    majority_vote_dict = {}

    for k in tqdm(k_range): 
        # Load aligned cluster labels for this k
        label_df = aligned_all_cluster_label_dict[str(k)]   # shape: (n_samples, n_folds)
        majority_vote_labels = label_df.mode(axis=1)[0].values    # Use majority vote across folds (眾數)
        majority_vote_dict[str(k)] = majority_vote_labels
        for fold in label_df.columns:
            fold_labels = label_df[fold].values
            metric_df.loc[k, f"fold{fold}_silhouette"] = silhouette_score(train_table, fold_labels, metric='cosine')

        majority_vote_silhouette = silhouette_score(train_table, majority_vote_labels, metric='cosine')
        mean_ari_5_fold = calculate_ari_matrix(aligned_all_cluster_label_dict, str(k)).mean().mean()  # fold 之間分群結果的一致程度
        chr_coherence_ratio =  compute_chr_coherence_ratio(majority_vote_labels, feature_cytoband)
    
        metric_df.loc[k, "majority_vote_silhouette"] = majority_vote_silhouette
        metric_df.loc[k, "mean_ari_5_fold"] = mean_ari_5_fold
        metric_df.loc[k, "chr_coherence_ratio"] = chr_coherence_ratio

    metric_df = metric_df.astype(float)
    metric_df.to_csv(f"./results/{cancer_type}/clustering/metric_df_{random_seed}.csv")
    with open(f"./results/{cancer_type}/clustering/majority_vote_dict_{random_seed}.json", 'w', encoding='utf-8') as f:
        json.dump(convert_ndarray_to_list(majority_vote_dict), f)

    return metric_df


def run_cnv_model(config):
    # 檢查必要參數
    required_keys = ["input_table", "clustering", "min_clusters", "max_clusters", "random_seed", "cancer_type"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")

    # 讀取參數
    input_table = config["input_table"]
    clustering = config.get("clustering", False)
    min_clusters = config.get("min_clusters", 2)
    max_clusters = config.get("max_clusters", 350)
    random_seed = config.get("random_seed", 42)
    cancer_type = config.get("cancer_type", "PAAD")
    fs_method = config.get("fs_method", "mrmr")
    fs_cutoff = config.get("fs_cutoff", 0)
    suffix = config.get("suffix", "")
    os.makedirs(f"./results/{cancer_type}/table", exist_ok=True)
    os.makedirs(f"./results/{cancer_type}/selected_features", exist_ok=True)
    os.makedirs(f"./results/{cancer_type}/clustering", exist_ok=True)
    print(f"\n##### Running model with random seed: {random_seed} and cancer type: {cancer_type} #####")

    # Load Data
    from pymaftools import PivotTable
    table_path = f"./results/{cancer_type}/table/cnv_table.sqlite"
    if os.path.exists(table_path):
        print(f"Loading existing table.")
        table = PivotTable.read_sqlite(table_path)

    elif cancer_type == "all":
        chol_table = PivotTable.read_sqlite(f"./results/CHOL/table/cnv_table.sqlite")
        paad_table = PivotTable.read_sqlite(f"./results/PAAD/table/cnv_table.sqlite")
        stad_table = PivotTable.read_sqlite(f"./results/STAD/table/cnv_table.sqlite")
        table = PivotTable.merge([chol_table, paad_table, stad_table])
        table.to_sqlite(db_path=table_path)

    else:
        all_table = PivotTable.read_sqlite(f"./results/{cancer_type}/table/{input_table}")
        common_samples = np.load(f"./results/{cancer_type}/common_samples.npy", allow_pickle=True)
        table = all_table.subset(samples=common_samples)

        # ---- Feature Filtering ----
        # Filter sex chromosome genes
        clean_table = pd.Series(table.feature_metadata["chromosome"]).dropna()
        mask = clean_table.apply(lambda x: x not in ["chrX", "chrY"])
        non_sex_chr = mask[mask].index.tolist()
        table = table.subset(features=non_sex_chr)

        # Protein coding genes
        gencode = pd.read_csv("./data/gencode_anno_v36.csv", index_col=False)
        pc_genes = gencode[gencode["gene_type"] == "protein_coding"]
        pc_gene_list = pc_genes["gene_name"].unique().tolist()
        clean_gene = pd.Series(table.feature_metadata["gene_name"])
        mask = clean_gene.isin(pc_gene_list)
        pc_gene_cnv = mask[mask].index.tolist()
        table = table.subset(features=pc_gene_cnv)
        print("Number of Protein coding genes:", len(pc_gene_cnv))

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

    # ---- Clustering ----
    # run k fold clustering evaluation (8.5HR)
    if clustering == True:
        print("Running Clustering...")
        metric_df = run_clustering(train_table, min_clusters=min_clusters, max_clusters=max_clusters, random_seed=random_seed, cancer_type=cancer_type)
    else:
        print("Skipped Clustering...")
        metric_df = pd.read_csv(f"./results/{cancer_type}/clustering/metric_df_{random_seed}.csv", index_col=0)
    
    # 找出 best_k 最穩定的結果與 penalty 建議區間
    from model_utils import sweep_penalty_select_best_k
    new_metric_df, penalties, best_k_for_penalty, most_common_k, freq_dict, ranges, chosen_penalty = sweep_penalty_select_best_k(metric_df, penalties_range=np.logspace(-6, np.log10(1/max_clusters), 100))

    # Plot
    plt.figure(figsize=(8,5))
    plt.plot(penalties, list(best_k_for_penalty.values()), marker='o', markersize=4, linestyle='-')
    plt.xscale('log')
    plt.xlabel('penalty (log scale)')
    plt.ylabel('clusters (k)')
    plt.title(f'Different penalty to select best_k')
    plt.grid(True)
    plt.axhline(most_common_k, color='red', linestyle='-.', alpha=0.8, label=f'best_k={most_common_k}')
    for a,b in ranges:
        plt.axvspan(a, b, color='red', alpha=0.12)
    plt.legend()
    plt.savefig(f"./results/{cancer_type}/clustering/penalty_sweep_{random_seed}.png", dpi=300, bbox_inches="tight")
    plt.show()

    # Find the best_k using chosen penalty
    from model_utils import plot_clustering_metrics_and_find_best_k
    tiff = f"./results/{cancer_type}/clustering/find_best_k_plot_{random_seed}.tiff"
    best_k = plot_clustering_metrics_and_find_best_k(new_metric_df, tiff, penalty=chosen_penalty, transparent=False)  # 選擇 penalized_silhouette 最大值為 best_k
    print(f"Best k: {best_k}")

    # 計算 by gene 相似度，類別會變成 SimilarityMatrix
    similarity_matrix = train_table.T.compute_similarity(method="cosine")
    distance_matrix = 1 - similarity_matrix.values
    np.fill_diagonal(distance_matrix, 0) # 用相減的對角線會是很小的數值 為了精準會直接設為零

    # 重新用 best_k 分群
    from sklearn.cluster import AgglomerativeClustering
    clustering_model = AgglomerativeClustering(
        n_clusters=best_k,         # 要分幾群
        metric='precomputed',       # 表示我們餵的是距離矩陣
        linkage='average'             # 適合非 Euclidean 距離（不能用 ward）
    )
    re_clustering_labels = clustering_model.fit_predict(distance_matrix)

    # 比較 majority_vote 和 re_clustering 的分群結果
    from sklearn.metrics import silhouette_score
    with open(f"./results/{cancer_type}/clustering/majority_vote_dict_{random_seed}.json", 'r', encoding='utf-8') as f:
        majority_vote_dict = json.load(f)
    majority_vote_labels = majority_vote_dict[str(best_k)]
    majority_vote_silhouette = silhouette_score(train_table, majority_vote_labels, metric='cosine')
    re_clustering_silhouette = silhouette_score(train_table, re_clustering_labels, metric='cosine')
    print(f"Majority Vote Silhouette: {majority_vote_silhouette}")
    print(f"Re-clustering Silhouette: {re_clustering_silhouette}")

    # Save clustering results
    table.feature_metadata["cluster"] = re_clustering_labels
    table.feature_metadata.cluster.to_csv(f"./results/{cancer_type}/clustering/best_k_cluster_labels_{random_seed}.csv", index=True)
    
    # summarize clusters
    cnv_tmp = pd.DataFrame(table)
    cnv_tmp["cluster"] = re_clustering_labels

    # to cluster table
    cluster_table = PivotTable(cnv_tmp.groupby("cluster").mean())
    cluster_table.sample_metadata = table.sample_metadata
    cluster_table.feature_metadata["unique_cytoband"] = (
        table.feature_metadata.groupby("cluster")["Cytoband"]
        .agg(lambda s: ",".join(pd.Series(s).dropna().astype(str).unique()))
    )
    cluster_table.feature_metadata["ids_in_cluster"] = (
        table.feature_metadata.groupby("cluster")["Gene ID"]
        .agg(lambda s: ",".join(pd.Series(s).dropna().astype(str).unique()))
    )
    cluster_table.feature_metadata["genes_in_cluster"]  = table.feature_metadata.groupby("cluster").apply(lambda df: ",".join(df.index), include_groups=False)
    cluster_table.feature_metadata["gene_count"] = table.feature_metadata.groupby("cluster").size()

    cluster_table.to_sqlite(db_path=f"./results/{cancer_type}/table/cnv_cluster_table_{random_seed}.sqlite")

    # 取出 cluster 過後的 training testing set
    train_samples = cluster_table.sample_metadata[split_col] == "train" 
    cluster_train = cluster_table.subset(samples=train_samples)
    test_samples = cluster_table.sample_metadata[split_col] == "test" 
    cluster_test = cluster_table.subset(samples=test_samples)

    # Define model
    from model_utils import build_base_model
    clf = build_base_model("RandomForest", cluster_train, random_state=42)

    # All features performance
    all_cluster_features = cluster_train.index
    from model_utils import evaluate_model
    report_train, report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff = evaluate_model(cluster_train, cluster_test, top_features=all_cluster_features, model=clf)
    print(report)

    # ---- Feature Selection Method ----
    # Top N Variance
    from sklearn.model_selection import StratifiedKFold
    from mrmr import mrmr_classif
    from sklearn.metrics import roc_auc_score

    if fs_method == "nvar":
        print(f"Feature selection method: Top N Variance / MAD / Prevalence")
        X = pd.DataFrame(cluster_train.T)
        y = pd.Series(cluster_train.sample_metadata.loc[cluster_train.columns]["TLS"])
        top_n = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20] # 你要選的特徵數
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

        row_var = cluster_train.var(axis=1)
        threshold = row_var.nlargest(best_k).min()
        top_features = row_var[row_var >= threshold].index
        top_features_path = f"./results/{cancer_type}/selected_features/top_features_cnv_{fs_method}{suffix}_{random_seed}.npy"
        np.save(top_features_path, top_features)
        print(f"Top features: {top_features}")
        print(f"Top features saved to: {top_features_path}")
        cluster_train = cluster_train.subset(features=top_features)
        cluster_test = cluster_test.subset(features=top_features)
    
    # Minimum Redundancy Maximum Relevance (mRMR)
    elif fs_method == "mrmr":
        print(f"Feature selection method: Minimum Redundancy Maximum Relevance (mRMR) to select K features (fs_cutoff={fs_cutoff})")
        X = pd.DataFrame(cluster_train.T)
        y = pd.Series(cluster_train.sample_metadata.loc[cluster_train.columns]["TLS"])

        mrmr_k = 20
        top_features = mrmr_classif(X=X, y=y, K=mrmr_k)  # 取前百分比的 features
        top_features_path = f"./results/{cancer_type}/selected_features/top_features_cnv_{fs_method}{suffix}_{random_seed}.npy"
        np.save(top_features_path, top_features)
        print(f"Top features: {top_features}")
        print(f"Top features saved to: {top_features_path}")
        cluster_train = cluster_train.subset(features=top_features)
        cluster_test = cluster_test.subset(features=top_features)
    
    # Random Forest
    elif fs_method == "rf":
        print(f"Feature selection method: Random Forest (fs_cutoff={fs_cutoff})")
        
        ### Feature Importance (Run 100 times)
        from model_utils import run_feature_importance
        os.makedirs(f"./results/{cancer_type}/importance_df", exist_ok=True)
        feature_importance_path = f"./results/{cancer_type}/importance_df/cnv_importance_{fs_method}{suffix}_{random_seed}.csv"
        global_importance = run_feature_importance(cluster_train, run_seeds=100, model=clf, save_path=feature_importance_path, plot_top=20)
        print(f"Feature importance saved to: {feature_importance_path}")

        ### Choose cutoff for selected features
        from model_utils import two_line_elbow
        knee_x, rss = two_line_elbow(global_importance)
        top_features = global_importance.head(knee_x).index
        top_features_path = f"./results/{cancer_type}/selected_features/top_features_cnv_{fs_method}{suffix}_{random_seed}.npy"
        np.save(top_features_path, top_features)
        print(f"Top features: {top_features}")
        print(f"Top features saved to: {top_features_path}")

    else:
        raise ValueError(f"Unknown feature selection method: {fs_method}")

    ### Model Evaluation
    from model_utils import evaluate_model
    report_train, report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff = evaluate_model(cluster_train, cluster_test, top_features=top_features, model=clf)
    print(report)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_cnv.json", help="Path to config JSON")
    args = parser.parse_args()
    config = load_config(args.config)
    run_cnv_model(config)

