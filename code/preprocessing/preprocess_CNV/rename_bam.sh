# Define variables
BASE_DIR="/staging/biology/$USER/CNV/PAAD"
sample_sheet="$BASE_DIR/data/gdc_sample_sheet_PAAD_tumor_filtered.tsv"
BAM_DIR="$BASE_DIR/tumor_bam_file"

# Read sample sheet from file
while IFS=$'\t' read -r File_ID File_Name Sample_ID; do

    echo "Rename $Sample_ID ..."
    Front_Name="${File_Name%.*}"

    # Rename BAM file and its index file
    mv $BAM_DIR/$File_ID/${Front_Name}.bam $BAM_DIR/${Sample_ID}.bam
    mv $BAM_DIR/$File_ID/${Front_Name}.bai $BAM_DIR/${Sample_ID}.bai

done < <(tail -n +2 "$sample_sheet" | awk -F"\t" '{print $1"\t"$2"\t"$7}')  # 跳過 header

# mkdir bad_bam
# samtools quickcheck -v /staging/biology/$USER/CNV/PAAD/normal_bam_done/*.bam | cut -d' ' -f1 | sort -u \
#   | xargs -I{} mv {} bad_bam/

