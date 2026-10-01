#!/usr/bin/env python
# coding: utf-8

from pymaftools import *
import json
import sys
sys.path.append("./code/utils/")

def load_config(config_path="config_ensemble.json"):
    with open(config_path, "r") as f:
        return json.load(f)


def run_ensemble_model(config):
    # 檢查必要參數
    required_keys = ["random_seed", "cancer_type", "model"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")

    # 讀取參數
    random_seed = config.get("random_seed", 42)
    cancer_type = config.get("cancer_type", "PAAD")
    model = config.get("model", None)
    fs_method = config.get("fs_method", "")
    suffix = config.get("suffix", "")
    if len(model) == 3 and isinstance(model, list):
        models_list = model
        print(f"Models list: {models_list}")
    else:
        raise ValueError("config['model'] must be a list with 3 models")

    # 讀取每個模型的 top features
    top_features_list = []
    for model_entry in models_list:
        model_name = model_entry.strip().lower()
        top_features = np.load(f"./results/{cancer_type}/selected_features/top_features_concat_{fs_method}{suffix}_{model_name}_{random_seed}.npy", allow_pickle=True)
        top_features_list.append(top_features)

    top_features_1 = top_features_list[0]
    top_features_2 = top_features_list[1]
    top_features_3 = top_features_list[2]

    #c3_3 = list(set(top_features_rf) & set(top_features_xgb) & set(top_features_cat))
    c3_2 = list(set(top_features_1) & set(top_features_2) | set(top_features_1) & set(top_features_3) | set(top_features_2) & set(top_features_3))
    #c3_1 = list(set(top_features_rf) | set(top_features_xgb) | set(top_features_cat))
    print(f"Selected by at least 2 models: {len(c3_2)} TLS signatures")
    print(f"TLS signatures: {c3_2}")

    arr = np.array(c3_2, dtype=object)
    signatures_path = f"./results/{cancer_type}/selected_features/TLS_signatures_{random_seed}.npy"
    np.save(signatures_path, arr)
    print(f"TLS signatures saved to: {signatures_path}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_ensemble.json", help="Path to config JSON")
    args = parser.parse_args()
    config = load_config(args.config)
    run_ensemble_model(config)
