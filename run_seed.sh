#!/bin/sh

# 任何指令失敗就停止，避免後續覆蓋錯誤結果
set -e

# 用法:
# sh run_seed.sh [config_path] [start_seed] [end_seed] [model] [cancer_type]
# 範例:
# sh run_seed.sh ./config.json 0 9 all PAAD
# sh run_seed.sh ./config.json 0 9 cnv PAAD

CONFIG_PATH="${1:-./config.json}"
START_SEED="${2:-0}"
END_SEED="${3:-9}"
MODEL="${4:-all}"
CANCER_TYPE="${5:-PAAD}"
RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"
TMP_DIR="./tmp/run_seed_${RUN_ID}"

if [ ! -f "$CONFIG_PATH" ]; then
	echo "Config not found: $CONFIG_PATH"
	exit 1
fi

# 建立資料夾
mkdir -p ./results/${CANCER_TYPE}/logs
mkdir -p "$TMP_DIR"

# 若中途中斷（例如 Ctrl+C），會清理暫存檔
cleanup() {
	rm -rf "$TMP_DIR"
}
trap cleanup EXIT INT TERM

# 依序跑 random seed
for seed in $(seq "$START_SEED" "$END_SEED"); do
	echo "Running seed=${seed}"
	SEED_CONFIG="$TMP_DIR/config_seed_${seed}.json"

	# 由指定 config 複製出本 seed 專屬設定檔，覆寫 random_seed 與 cancer_type
	python -c "import json; src='${CONFIG_PATH}'; dst='${SEED_CONFIG}'; c=json.load(open(src)); c['random_seed']=${seed}; c['cancer_type']='${CANCER_TYPE}'; json.dump(c, open(dst, 'w'), indent=4)"

	# 依指定 model 執行
	if [ "$MODEL" = "all" ] || [ "$MODEL" = "rna" ]; then
		python ./code//feature_selection/RNA_mrmr.py --config "$SEED_CONFIG" > "./results/${CANCER_TYPE}/logs/RNA_mrmr_${seed}.log"
	fi
	if [ "$MODEL" = "all" ] || [ "$MODEL" = "meth" ]; then
		python ./code//feature_selection/DNAmeth_mrmr.py --config "$SEED_CONFIG" > "./results/${CANCER_TYPE}/logs/DNAmeth_mrmr_${seed}.log"
	fi
	if [ "$MODEL" = "all" ] || [ "$MODEL" = "mut" ]; then
		python ./code/feature_selection/DNAmut_mrmr.py --config "$SEED_CONFIG" > "./results/${CANCER_TYPE}/logs/DNAmut_mrmr_${seed}.log"
	fi
	if [ "$MODEL" = "all" ] || [ "$MODEL" = "cnv" ]; then
		python ./code/feature_selection/CNV_mrmr.py --config "$SEED_CONFIG" > "./results/${CANCER_TYPE}/logs/CNV_mrmr_${seed}.log"
	fi
    if [ "$MODEL" = "all" ] || [ "$MODEL" = "ensemble" ]; then
		python ./code/feature_selection/ensemble_importance.py --config "$SEED_CONFIG" > "./results/${CANCER_TYPE}/logs/ensemble_importance_mrmr_${seed}.log"
		python ./code/feature_selection/consensus_selection.py --config "$SEED_CONFIG" > "./results/${CANCER_TYPE}/logs/consensus_selection_mrmr_${seed}.log"
	fi
done

# 正常結束時取消 trap 並清理暫存檔
trap - EXIT INT TERM
cleanup