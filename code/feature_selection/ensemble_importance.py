#!/usr/bin/env python
# coding: utf-8

from pymaftools import *
import json
import sys
import os
sys.path.append("./code/utils/")

def load_config(config_path="config_ensemble.json"):
    with open(config_path, "r") as f:
        return json.load(f)

def load_omics_data(cancer_type, random_seed, cluster=False):
    base_path = f"./results/{cancer_type}"

    from pymaftools import PivotTable
    if cluster:
        rna = PivotTable.read_sqlite(f"{base_path}/table/rna_table.sqlite")
        meth = PivotTable.read_sqlite(f"{base_path}/table/meth_table.sqlite")
        mut = PivotTable.read_sqlite(f"{base_path}/table/mutation_table.sqlite")
        cnv = PivotTable.read_sqlite(f"{base_path}/table/cnv_cluster_table_{random_seed}.sqlite")
    else:
        rna = PivotTable.read_sqlite(f"{base_path}/table/rna_table.sqlite")
        meth = PivotTable.read_sqlite(f"{base_path}/table/meth_table.sqlite")
        mut = PivotTable.read_sqlite(f"{base_path}/table/mutation_table.sqlite")
        cnv = PivotTable.read_sqlite(f"{base_path}/table/cnv_table.sqlite")  # 尚未 cluster

    return rna, meth, mut, cnv

def split_train_test(table, random_seed):
    split_col = f"split_{random_seed}"

    # Keep sample order as train then test for the current seed
    train_idx = table.sample_metadata.index[table.sample_metadata[split_col] == "train"].tolist()
    test_idx = table.sample_metadata.index[table.sample_metadata[split_col] == "test"].tolist()
    table = table.reindex(columns=train_idx + test_idx)

    train_samples = table.sample_metadata[split_col] == "train"
    test_samples = table.sample_metadata[split_col] == "test"
    train_table = table.subset(samples=train_samples)
    test_table = table.subset(samples=test_samples)
    return train_table, test_table


def run_ensemble_model(config):
    # 檢查必要參數
    required_keys = ["random_seed", "cancer_type", "scaler_dict"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")

    # 讀取參數
    random_seed = config.get("random_seed", 42)
    cancer_type = config.get("cancer_type", "PAAD")
    model = config.get("model", None)
    fs_method = config.get("fs_method", "")
    scaler_dict = config.get("scaler_dict", None)
    suffix = config.get("suffix", "")
    if isinstance(model, str):
        models_list = [model]
    elif isinstance(model, list):
        models_list = model
    else:
        raise ValueError("config['model'] must be a string or list")

    os.makedirs(f"./results/{cancer_type}/selected_features", exist_ok=True)
    os.makedirs(f"./results/{cancer_type}/importance_df", exist_ok=True)
    print(f"\n##### Running model with random seed: {random_seed} and cancer type: {cancer_type} (scaler_dict = {scaler_dict}) #####")

    # Load Data and assign training and independent datasets
    rna_table, meth_table, mut_table, cnv_table = load_omics_data(cancer_type, random_seed, cluster=True)

    # 取出該 cohort 的 train test dataset
    rna_train, rna_test = split_train_test(rna_table, random_seed)
    meth_train, meth_test = split_train_test(meth_table, random_seed)
    mut_train, mut_test = split_train_test(mut_table, random_seed)
    cnv_train, cnv_test = split_train_test(cnv_table, random_seed)

    top_features_rna = np.load(f"./results/{cancer_type}/selected_features/top_features_rna_{fs_method}{suffix}_{random_seed}.npy", allow_pickle=True)
    rna_train_top = rna_train.subset(features=top_features_rna)
    rna_test_top = rna_test.subset(features=top_features_rna)
    top_features_meth = np.load(f"./results/{cancer_type}/selected_features/top_features_meth_{fs_method}{suffix}_{random_seed}.npy", allow_pickle=True)
    meth_train_top = meth_train.subset(features=top_features_meth)
    meth_test_top = meth_test.subset(features=top_features_meth)
    top_features_mut = np.load(f"./results/{cancer_type}/selected_features/top_features_mut_{fs_method}{suffix}_{random_seed}.npy", allow_pickle=True)
    mut_train_top = mut_train.subset(features=top_features_mut)
    mut_test_top = mut_test.subset(features=top_features_mut)
    top_features_cnv = np.load(f"./results/{cancer_type}/selected_features/top_features_cnv_{fs_method}{suffix}_{random_seed}.npy", allow_pickle=True)
    cnv_train_top = cnv_train.subset(features=top_features_cnv)
    cnv_test_top = cnv_test.subset(features=top_features_cnv)

    # Concat matrix
    from model_utils import concat_omics
    concat_train = concat_omics([rna_train_top, meth_train_top, mut_train_top, cnv_train_top], keys=["RNA", "Meth", "Mutation", "CNV"], axis=0)
    concat_test = concat_omics([rna_test_top, meth_test_top, mut_test_top, cnv_test_top], keys=["RNA", "Meth", "Mutation", "CNV"], axis=0)

    for model_entry in models_list:
        print(f"\n==== Running model: {model_entry} ====")

        # 從頭到尾固定，不調參 (保守參數確保 importance 可信即可)
        from model_utils import build_base_model
        base_model = build_base_model(model_entry, concat_train, random_state=42)

        model_name = model_entry.strip().lower()

        # All features performance
        all_features = concat_train.index
        from model_utils import evaluate_model
        report_train, report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff = evaluate_model(concat_train, concat_test, top_features=all_features, model=base_model, scaler_dict=scaler_dict)
        print(report)

        ### Feature Importance (Run 100 times)
        from model_utils import run_feature_importance
        feature_importance_path = f"./results/{cancer_type}/importance_df/concat_importance_{fs_method}{suffix}_{model_name}_{random_seed}.csv"
        global_importance = run_feature_importance(concat_train, run_seeds=100, model=base_model, concat=True, scaler_dict=scaler_dict,
                                                   save_path=feature_importance_path, plot_top=20)

        # Choose cutoff for selected features
        from model_utils import two_line_elbow
        knee_x, rss = two_line_elbow(global_importance)
        top_features = global_importance.head(knee_x).index
        top_features_path = f"./results/{cancer_type}/selected_features/top_features_concat_{fs_method}{suffix}_{model_name}_{random_seed}.npy"
        np.save(top_features_path, top_features)
        print(f"Top features: {top_features}")
        print(f"Top features saved to: {top_features_path}")

        # Evaluate on test dataset
        from model_utils import evaluate_model
        print(f"### Performance on {cancer_type} dataset ###")
        report_train, report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff = evaluate_model(concat_train, concat_test, top_features=top_features, model=base_model, concat=True, scaler_dict=scaler_dict)
        print(report)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_ensemble.json", help="Path to config JSON")
    args = parser.parse_args()
    config = load_config(args.config)
    run_ensemble_model(config)
