#!/bin/bash

# Load packages in conda
source "/staging/biology/$USER/miniconda3/etc/profile.d/conda.sh"

# 需要 token
GDC_token="/staging/biology/$USER/gdc-user-token.2025-10-07T10_43_58.108Z.txt"
download_dir="/staging/biology/$USER/CNV/PAAD/normal_bam_file/"
sample_sheet="/staging/biology/$USER/CNV/PAAD/data/gdc_sample_sheet_PAAD_normal.tsv"
log_dir="/staging/biology/$USER/CNV/PAAD/logs/"

[[ ! -d $download_dir ]] && mkdir -p "$download_dir"
[[ ! -d $log_dir ]] && mkdir -p "$log_dir"

download_log="${log_dir}/gdc_download.log"
success_log="${log_dir}/gdc_download_success.log"
failed_ids_file="${log_dir}/gdc_failed_ids.txt"

# 清除先前的 failed list（若希望累積可拿掉）
> "$download_log"
> "$failed_ids_file"

echo "Download started at $(date -u +'%Y-%m-%dT%H:%M:%SZ')" | tee -a "$download_log"

# 讀取 sample_sheet 的第一欄 (跳過 header)，忽略空白列
while IFS=$'\t' read -r first_col _; do
    GDC_ID="$first_col"
    # skip header or empty
    if [[ "$GDC_ID" == "File ID" || -z "$GDC_ID" ]]; then
        continue
    fi

    echo "----------------------------------------" | tee -a "$download_log"
    echo "Downloading $GDC_ID ..." | tee -a "$download_log"

    timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

    # single attempt: run download, capture stdout/stderr to log
    gdc-client download -t "$GDC_token" -d "$download_dir" "$GDC_ID" >> "$download_log" 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then
        echo "$timestamp $GDC_ID SUCCESS" >> "$success_log"
    else
        echo "$timestamp $GDC_ID FAIL rc=$rc" >> "$download_log"
        echo "$GDC_ID" >> "$failed_ids_file"
        echo "$(date -u +"%Y-%m-%dT%H:%M:%SZ") $GDC_ID FAILED" >> "$download_log"
    fi

done < <(tail -n +2 "$sample_sheet" | awk -F"\t" '{print $1}')

echo "Download finished at $(date -u +'%Y-%m-%dT%H:%M:%SZ')" | tee -a "$download_log"

# summary
total=$(tail -n +2 "$sample_sheet" | awk -F"\t" 'NF{print $1}' | wc -l)
failed=$(wc -l < "$failed_ids_file" || true)
success_count=$(wc -l < "$success_log" || true)

echo "Total IDs in sheet: $total" | tee -a "$download_log"
echo "Successful downloads (logged): $success_count" | tee -a "$download_log"
echo "Failed downloads (record file: $failed_ids_file): $failed" | tee -a "$download_log"

echo "All done. See $download_log for details."
