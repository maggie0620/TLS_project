#!/usr/bin/sh
#SBATCH -A MST109178     # 計畫系統代號
#SBATCH -J cnvkit  # Job name
#SBATCH -p ngs186G        # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 28              # 使用的core數
#SBATCH --mem=186G        # 使用的記憶體量
#SBATCH --mail-user=maggie910620.md13@nycu.edu.tw   # Email
#SBATCH --mail-type=ALL    # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL

# Load packages in conda
source "/staging/biology/$USER/miniconda3/etc/profile.d/conda.sh"
conda activate CNVkit  # conda install pandas=1.5.3

# 下載 refFlat.txt 和 access-10kb.hg38.bed
BAM_DIR="/staging/biology/$USER/CNV/PAAD"
refFLAT="/staging/biology/$USER/CNV/refFlat.txt"
TARGET_BED="/staging/biology/$USER/CNV/SureSelectHumanAllExonV5_50Mbp_hg38.bed"
ACCESS_BED="/staging/biology/$USER/CNV/access-10kb.hg38.bed"
FASTA="/staging/biology/$USER/CNV/Homo_sapiens_assembly38.fasta"
THREADS=20

CNVKIT_DIR="${BAM_DIR}/CNVkit_results"
[[ ! -d $CNVKIT_DIR ]] && mkdir -p "$CNVKIT_DIR"
REF_CNN="${CNVKIT_DIR}/my_reference.cnn"

###################### Run CNVkit ############################

echo --- Running CNVkit reference build ---

# Create normal sample reference (--normal 表示這批都是正常樣本，用於建立 baseline reference)
cnvkit.py batch \
    --normal ${BAM_DIR}/normal_bam_file/*.bam \
    --targets $TARGET_BED \
    --annotate $refFLAT \
    --fasta $FASTA \
    --access $ACCESS_BED \
    -p $THREADS \
    --output-reference $REF_CNN \
    --output-dir $CNVKIT_DIR

rm ${CNVKIT_DIR}/*.targetcoverage.cnn
rm ${CNVKIT_DIR}/*.antitargetcoverage.cnn


echo --- Running CNVkit for tumor samples ---

# For tumor samples only
cnvkit.py batch ${BAM_DIR}/tumor_bam_file/*.bam \
    --reference $REF_CNN \
    --output-dir $CNVKIT_DIR \
    -p $THREADS \
    --diagram --scatter
    # --diagram: 根據 gain 或 loss 進行結果可視化; --scatter: 根據 log2 值進行結果可視化

rm ${CNVKIT_DIR}/*.cnr
rm ${CNVKIT_DIR}/*.bintest.cns
rm ${CNVKIT_DIR}/*.call.cns
rm ${CNVKIT_DIR}/*.targetcoverage.cnn
rm ${CNVKIT_DIR}/*.antitargetcoverage.cnn

echo --- Export segment file ---

# Export segment file for GISTIC
cnvkit.py export seg \
    ${CNVKIT_DIR}/*.cns \
    -o ${CNVKIT_DIR}/merged_segments.seg

# ###################### Run GISTIC ############################

conda activate gistic

GISTIC_RESULTS_DIR="${BAM_DIR}/gistic_results"
[[ ! -d $GISTIC_RESULTS_DIR ]] && mkdir -p "$GISTIC_RESULTS_DIR"
SEGFILE="${CNVKIT_DIR}/merged_segments.seg"
REFGENE="/staging/biology/$USER/CNV/hg38.UCSC.add_miR.160920.refgene.mat"

echo --- running GISTIC ---
gistic2 -b $GISTIC_RESULTS_DIR -seg $SEGFILE -refgene $REFGENE \
    -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme

echo "GISTIC completed."

