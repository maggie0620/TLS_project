#!/usr/bin/env python
# coding: utf-8

from pymaftools import *
import json
import sys
import os
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score, confusion_matrix
sys.path.append("./code/utils/")

def load_config(config_path="config_final_model.json"):
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

def process_clustering(cnv_table, cancer_type, random_seed):
    """
    Process clustering for a given CNV table.
    """
    cluster_labels_df = pd.read_csv(f"./results/{cancer_type}/clustering/best_k_cluster_labels_{random_seed}.csv", index_col=0)
    if "cluster" in cluster_labels_df.columns:
        gene_cluster_series = cluster_labels_df["cluster"].astype(int)
    else:
        gene_cluster_series = cluster_labels_df.iloc[:, 0].astype(int)

    # Map clustering results
    mapping = gene_cluster_series.reindex(cnv_table.index)
    cnv_table.feature_metadata["cluster"] = mapping

    # Summarize clusters
    cnv_tmp = pd.DataFrame(cnv_table)
    cnv_tmp["cluster"] = mapping.loc[cnv_tmp.index].values

    # Create cluster table
    cnv_cluster = PivotTable(cnv_tmp.groupby("cluster").mean())
    cnv_cluster.sample_metadata = cnv_table.sample_metadata
    cnv_cluster.feature_metadata["unique_cytoband"] = (
        cnv_table.feature_metadata.groupby("cluster")["Cytoband"]
        .agg(lambda s: ",".join(pd.Series(s).dropna().astype(str).unique()))
    )
    cnv_cluster.feature_metadata["ids_in_cluster"] = (
        cnv_table.feature_metadata.groupby("cluster")["Gene ID"]
        .agg(lambda s: ",".join(pd.Series(s).dropna().astype(str).unique()))
    )
    cnv_cluster.feature_metadata["genes_in_cluster"] = cnv_table.feature_metadata.groupby("cluster").apply(lambda df: ",".join(df.index), include_groups=False)
    cnv_cluster.feature_metadata["gene_count"] = cnv_table.feature_metadata.groupby("cluster").size()

    return cnv_cluster


def plot_metrics(train_curve, curve_data, cancer_type, random_seed, suffix):
    """Plot and save training and evaluation ROC/PR curves and confusion matrices."""
    fig_dir = f"./results/{cancer_type}/fig"
    os.makedirs(fig_dir, exist_ok=True)

    # Keep colors consistent with cancer types across ROC/PR figures.
    cancer_colors = {
        "PAAD": "#7E5A92",  # purple
        "CHOL": "#4FA6A6",  # teal
        "STAD": "#ffcf33",  # yellow
    }

    def _unpack_curve(item, default_name=None):
        if item is None:
            return None
        if isinstance(item, dict):
            return (
                item.get("name", default_name),
                np.asarray(item["y_true"]),
                np.asarray(item["y_score"]),
                float(item.get("threshold", 0.5)),
            )
        if len(item) == 4:
            name, y_true, y_score, threshold = item
        else:
            name, y_true, y_score = item
            threshold = 0.5
        return name, np.asarray(y_true), np.asarray(y_score), float(threshold)

    def _plot_threshold_point(ax, x_values, y_values, thresholds, threshold, color):
        if thresholds is None or len(thresholds) == 0:
            return
        finite_mask = np.isfinite(thresholds)
        if not finite_mask.any():
            return
        finite_thresholds = thresholds[finite_mask]
        finite_x = np.asarray(x_values)[finite_mask]
        finite_y = np.asarray(y_values)[finite_mask]
        idx = np.abs(finite_thresholds - threshold).argmin()
        ax.scatter(
            finite_x[idx],
            finite_y[idx],
            color=color,
            marker="o",
            s=80,
            zorder=7,
            edgecolors="white",
            linewidths=1
        )
    train_item = _unpack_curve(train_curve, default_name=cancer_type)

    # ROC curve
    plt.figure(figsize=(7, 6), dpi=300)
    valid_roc = 0
    if train_item is not None:
        name, y_true, y_score, threshold = train_item
        if len(np.unique(y_true)) >= 2:
            fpr, tpr, roc_thresholds = roc_curve(y_true, y_score, drop_intermediate=False)
            roc_auc = auc(fpr, tpr)
            line_color = cancer_colors.get(name, "#4C72B0")
            plt.plot(fpr, tpr, lw=2, color=line_color, linestyle="--", drawstyle="steps-post",
                     label=f"{name} train (AUC={roc_auc:.3f})")  # , thr={threshold:.3f}
            _plot_threshold_point(plt.gca(), fpr, tpr, roc_thresholds, threshold, line_color)
            valid_roc += 1
        else:
            print(f"[Warning] Skip train ROC for {name}: only one class in y_true.")

    for name, y_true, y_score, threshold in curve_data:
        y_true = np.asarray(y_true)
        y_score = np.asarray(y_score)

        if len(np.unique(y_true)) < 2:
            print(f"[Warning] Skip ROC for {name}: only one class in y_true.")
            continue

        fpr, tpr, roc_thresholds = roc_curve(y_true, y_score, drop_intermediate=False)
        roc_auc = auc(fpr, tpr)
        line_color = cancer_colors.get(name, "#4C72B0")  # default blue if cancer type not in dict
        plt.plot(fpr, tpr, lw=2, color=line_color, drawstyle="steps-post", label=f"{name} test (AUC={roc_auc:.3f})")
        _plot_threshold_point(plt.gca(), fpr, tpr, roc_thresholds, threshold, line_color)
        valid_roc += 1

    plt.plot([0, 1], [0, 1], linestyle="--", color="gray", lw=1)
    plt.xlabel("False Positive Rate", fontsize=16)
    plt.ylabel("True Positive Rate", fontsize=16)
    plt.title(f"ROC Curve Across Cancers ({cancer_type}-trained model)", fontsize=18)
    if valid_roc > 0:
        plt.legend(loc="lower right", fontsize=14)
    plt.tight_layout()
    roc_path = f"{fig_dir}/ROC_curve{suffix}_{random_seed}.tiff"
    plt.savefig(roc_path, dpi=300, format='tiff', bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()

    # PR curve
    plt.figure(figsize=(7, 6), dpi=300)
    valid_pr = 0
    if train_item is not None:
        name, y_true, y_score, threshold = train_item
        if len(np.unique(y_true)) >= 2:
            precision, recall, pr_thresholds = precision_recall_curve(y_true, y_score)
            pr_auc = average_precision_score(y_true, y_score)
            line_color = cancer_colors.get(name, "#4C72B0")
            plt.plot(recall, precision, lw=2, color=line_color, linestyle="--", drawstyle="steps-post",
                     label=f"{name} train (AUPRC={pr_auc:.3f})")
            finite_mask = np.isfinite(pr_thresholds)
            if finite_mask.any():
                finite_thresholds = pr_thresholds[finite_mask]
                idx = np.abs(finite_thresholds - threshold).argmin()
                # precision_recall_curve returns thresholds of length n-1
                point_idx = np.flatnonzero(finite_mask)[idx]
                plt.scatter(
                    recall[point_idx + 1],
                    precision[point_idx + 1],
                    color=line_color,
                    marker="o",
                    s=80,
                    zorder=7,
                    edgecolors="white",
                    linewidths=1
                )
            valid_pr += 1
        else:
            print(f"[Warning] Skip train PR for {name}: only one class in y_true.")

    for name, y_true, y_score, threshold in curve_data:
        y_true = np.asarray(y_true)
        y_score = np.asarray(y_score)

        if len(np.unique(y_true)) < 2:
            print(f"[Warning] Skip PR for {name}: only one class in y_true.")
            continue

        precision, recall, pr_thresholds = precision_recall_curve(y_true, y_score)
        pr_auc = average_precision_score(y_true, y_score)
        line_color = cancer_colors.get(name, "#4C72B0")  # default blue if cancer type not in dict
        plt.plot(recall, precision, lw=2, color=line_color, drawstyle="steps-post", label=f"{name} test (AUPRC={pr_auc:.3f})")
        finite_mask = np.isfinite(pr_thresholds)
        if finite_mask.any():
            finite_thresholds = pr_thresholds[finite_mask]
            idx = np.abs(finite_thresholds - threshold).argmin()
            point_idx = np.flatnonzero(finite_mask)[idx]
            plt.scatter(
                recall[point_idx + 1],
                precision[point_idx + 1],
                color=line_color,
                marker="o",
                s=80,
                zorder=7,
                edgecolors="white",
                linewidths=1
            )
        valid_pr += 1

    plt.xlabel("Recall", fontsize=16)
    plt.ylabel("Precision", fontsize=16)
    plt.title(f"PR Curve Across Cancers ({cancer_type}-trained model)", fontsize=18)
    if valid_pr > 0:
        plt.legend(loc="lower left", fontsize=14)
    plt.tight_layout()
    pr_path = f"{fig_dir}/PR_curve{suffix}_{random_seed}.tiff"
    plt.savefig(pr_path, dpi=300, format='tiff', bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()

    # Confusion matrix
    total_panels = len(curve_data) + (1 if train_item is not None else 0)
    fig, axes = plt.subplots(
        1,
        total_panels,
        figsize=(5 * total_panels, 4.5),
        dpi=300,
        constrained_layout=True,
    )
    if total_panels == 1:
        axes = [axes]
    elif train_item is not None and total_panels == 2:
        axes = list(np.atleast_1d(axes))
    else:
        axes = list(np.atleast_1d(axes))

    cm_image = None
    panel_items = []
    if train_item is not None:
        panel_items.append(("train", train_item))
    panel_items.extend(("test", _unpack_curve(item)) for item in curve_data)

    for ax, (split_name, item) in zip(axes, panel_items):
        if len(item) == 4:
            name, y_true, y_score, threshold = item
        else:
            name, y_true, y_score = item
            threshold = 0.5

        y_true = np.asarray(y_true)
        y_score = np.asarray(y_score)
        y_pred = (y_score >= threshold).astype(int)
        cm = confusion_matrix(y_true, y_pred)

        cm_image = ax.imshow(cm, cmap="Blues")
        text_color_threshold = cm.max() * 0.5 if cm.size > 0 else 0
        for (row, col), value in np.ndenumerate(cm):
            text_color = "white" if value > text_color_threshold else "black"
            ax.text(col, row, int(value), ha="center", va="center", fontsize=16, color=text_color)

        ax.set_title(f"{name} {split_name} (cutoff={threshold:.3f})", fontsize=18)
        ax.set_xlabel("Predicted", fontsize=16)
        ax.set_ylabel("Actual", fontsize=16)
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(["0", "1"], fontsize=14)
        ax.set_yticklabels(["0", "1"], fontsize=14)

    if cm_image is not None:
        fig.colorbar(cm_image, ax=axes, fraction=0.02, pad=0.04)
    fig.suptitle(f"Confusion Matrix Across Cancers ({cancer_type}-trained model)", y=1.08, fontsize=18)
    cm_path = f"{fig_dir}/confusion_matrix{suffix}_{random_seed}.tiff"
    fig.savefig(cm_path, dpi=300, format='tiff', bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)

    print(f"Saved ROC figure: {roc_path}")
    print(f"Saved PR figure: {pr_path}")
    print(f"Saved CM figure: {cm_path}")


def train_final_model(config):
    # 檢查必要參數
    required_keys = ["random_seed", "cancer_type", "suffix"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")

    # 讀取參數
    random_seed = config.get("random_seed", 42)
    cancer_type = config.get("cancer_type", "PAAD")
    model = config.get("model", None)
    scaler_dict = config.get("scaler_dict", None)
    suffix = config["suffix"]  # consensus or pooled
    print(f"\n##### Running model with random seed: {random_seed} and cancer type: {cancer_type} (scaler_dict = {scaler_dict}) #####")
    
    # Load Data and assign training and independent datasets
    if cancer_type == "all":
        rna_table, meth_table, mut_table, cnv_table = load_omics_data("all", random_seed, cluster=True)
        independent_datasets = []
    elif cancer_type == "PAAD":
        PAAD_rna, PAAD_meth, PAAD_mut, PAAD_cnv = load_omics_data("PAAD", random_seed, cluster=True)
        CHOL_rna, CHOL_meth, CHOL_mut, CHOL_cnv = load_omics_data("CHOL", random_seed)
        STAD_rna, STAD_meth, STAD_mut, STAD_cnv = load_omics_data("STAD", random_seed)
        rna_table, meth_table, mut_table, cnv_table = PAAD_rna, PAAD_meth, PAAD_mut, PAAD_cnv
        independent_datasets = [("CHOL", CHOL_rna, CHOL_meth, CHOL_mut, CHOL_cnv),
                                ("STAD", STAD_rna, STAD_meth, STAD_mut, STAD_cnv)]
    elif cancer_type == "CHOL":
        PAAD_rna, PAAD_meth, PAAD_mut, PAAD_cnv = load_omics_data("PAAD", random_seed)
        CHOL_rna, CHOL_meth, CHOL_mut, CHOL_cnv = load_omics_data("CHOL", random_seed, cluster=True)
        STAD_rna, STAD_meth, STAD_mut, STAD_cnv = load_omics_data("STAD", random_seed)
        rna_table, meth_table, mut_table, cnv_table = CHOL_rna, CHOL_meth, CHOL_mut, CHOL_cnv
        independent_datasets = [("PAAD", PAAD_rna, PAAD_meth, PAAD_mut, PAAD_cnv),
                                ("STAD", STAD_rna, STAD_meth, STAD_mut, STAD_cnv)]
    else:
        PAAD_rna, PAAD_meth, PAAD_mut, PAAD_cnv = load_omics_data("PAAD", random_seed)
        CHOL_rna, CHOL_meth, CHOL_mut, CHOL_cnv = load_omics_data("CHOL", random_seed)
        STAD_rna, STAD_meth, STAD_mut, STAD_cnv = load_omics_data("STAD", random_seed, cluster=True)
        rna_table, meth_table, mut_table, cnv_table = STAD_rna, STAD_meth, STAD_mut, STAD_cnv
        independent_datasets = [("CHOL", CHOL_rna, CHOL_meth, CHOL_mut, CHOL_cnv),
                                ("PAAD", PAAD_rna, PAAD_meth, PAAD_mut, PAAD_cnv)]

    # 取出該 cohort 的 train test dataset
    rna_train, rna_test = split_train_test(rna_table, random_seed)
    meth_train, meth_test = split_train_test(meth_table, random_seed)
    mut_train, mut_test = split_train_test(mut_table, random_seed)
    cnv_train, cnv_test = split_train_test(cnv_table, random_seed)
    
    tls_sig = np.load(f"/home/data/data_YM313112004/TLS_model/results/{cancer_type}/selected_features/TLS_signatures{suffix}_{random_seed}.npy", allow_pickle=True)

    # Concat matrix
    from model_utils import concat_omics
    all_concat_train = concat_omics([rna_train, meth_train, mut_train, cnv_train], keys=["RNA", "Meth", "Mutation", "CNV"], axis=0)
    concat_train = all_concat_train[all_concat_train.index.isin(tls_sig)]
    
    all_concat_test = concat_omics([rna_test, meth_test, mut_test, cnv_test], keys=["RNA", "Meth", "Mutation", "CNV"], axis=0)
    concat_test = all_concat_test[all_concat_test.index.isin(tls_sig)]
    print(f"TLS Signatures on {cancer_type}: {tls_sig.tolist()}")
    
    # # 從頭到尾固定，不調參 (保守參數確保 importance 可信即可)
    # from model_utils import build_base_model
    # base_model = build_base_model(model, all_concat_train, random_state=42)

    # TLS signatures performance
    # from model_utils import evaluate_model
    # print(f"### Performance on {cancer_type} dataset ###")
    # report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff = evaluate_model(all_concat_train, all_concat_test, top_features=tls_sig, model=base_model, concat=True, scaler_dict=scaler_dict)
    # print(report)

    # Define & tune model after determine TLS signatures
    if isinstance(model, str):
        model_name = model
    elif model is not None:
        model_name = model.__class__.__name__
    else:
        raise ValueError("Use: RandomForest | XGBoost | CatBoost")

    print(f"Tuning {model_name} with selected features...")
    from model_utils import tune_model
    best_model, best_params, best_score = tune_model(concat_train, model_name=model_name, random_state=42, n_splits=5)
    print(f"Usage model: {best_model}")
    print("Best params:", best_params)
    print("Best CV AUC:", best_score)

    # Evaluate on test dataset
    from model_utils import evaluate_model
    print(f"### Performance on {cancer_type} dataset ###")
    report_train, report, y_train, y_pred_prob, y_test, y_test_pred_prob, best_cutoff = evaluate_model(concat_train, concat_test, top_features=tls_sig, model=best_model, concat=True, scaler_dict=scaler_dict)
    print(report_train)
    print(report)

    train_curve = (cancer_type, y_train, y_pred_prob, best_cutoff)
    curve_data = [(cancer_type, y_test, y_test_pred_prob, best_cutoff)]

    # Evaluate on independent datasets
    from model_utils import evaluate_model
    for ind_cancer, ind_rna, ind_meth, ind_mut, ind_cnv in independent_datasets:
        ind_cnv_cluster = process_clustering(ind_cnv, cancer_type, random_seed)
        ind_concat = concat_omics([ind_rna, ind_meth, ind_mut, ind_cnv_cluster], keys=["RNA", "Meth", "Mutation", "CNV"], axis=0)
        print(f"### Performance on {ind_cancer} independent dataset ###")
        report_train, report, _, _, y_test, y_test_pred_prob, best_cutoff = evaluate_model(concat_train, ind_concat, top_features=tls_sig, model=best_model, concat=True, scaler_dict=scaler_dict)
        print(report_train)
        print(report)
        curve_data.append((ind_cancer, y_test, y_test_pred_prob, best_cutoff))

    # Plot ROC curves, PR curves and confusion matrix for the three cancers
    plot_metrics(train_curve=train_curve, curve_data=curve_data, cancer_type=cancer_type, random_seed=random_seed, suffix=suffix)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_final_model.json", help="Path to config JSON")
    args = parser.parse_args()
    config = load_config()
    train_final_model(config)
