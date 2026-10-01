# 設定來源資料夾和目標資料夾
source_dir="/home/maggie0620/TLS_meth/TCGA-CHOL_RNA/"  # 注意最後的 "/"
target_dir="/home/maggie0620/TLS_meth/DATA_CHOL_T"

# 創建目標資料夾
mkdir -p "$target_dir"

# 找到來源資料夾中的所有檔案並複製到目標資料夾
find "$source_dir" -type f -exec cp {} "$target_dir" \;

echo "All the files are copied successfully to: $target_dir"