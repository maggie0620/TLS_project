cancer_type = "PAAD"

# TLS label
library(readxl)
library(tidyverse)
library(maftools)
TLS_sheet <- read_excel(paste0("./data/", cancer_type, "/TLS vs. model_filtered (TCGA-", cancer_type, ").xlsx"))
if (cancer_type == "CHOL"){
  colnames(TLS_sheet)[8] <- "Sample.ID"
  colnames(TLS_sheet)[12] <- "model"
} else{
  colnames(TLS_sheet)[7] <- "Sample.ID"
  colnames(TLS_sheet)[11] <- "model"
}
TLS_label <- select(TLS_sheet, Sample.ID, TLS, model)

# 只要有一個TLS = 1，就是TLS-positive
TLS_new_label <- TLS_label %>%
  group_by(Sample.ID) %>%
  summarise(TLS = as.integer(any(TLS == 1)), model = as.integer(any(model == 1)), .groups = "drop")

# Change sample sheet.tsv into .csv
gdc_data <- read.delim(paste0("./data/", cancer_type, "/metadata/gdc_sample_sheet_mut.tsv"), header = TRUE, sep = "\t")

# Get file path
files_dir <- gdc_data %>%
  filter(grepl("aliquot_ensemble_masked.maf", File.Name)) %>%
  rowwise() %>%  # 使用 rowwise 操作每一列
  mutate(
    sample_list = list(str_split(Sample.ID, ",\\s*")[[1]]),
    Tumor.Sample.ID = paste(sample_list[str_detect(sample_list, "01A|06A")], collapse = ", "),
    Normal.Sample.ID = paste(sample_list[str_detect(sample_list, "10A|11A|10B|10E")], collapse = ", ")
  ) %>%
  ungroup() %>%
  mutate(Basename = paste0("~/TLS_meth/DATA_", cancer_type, "_mut/", sub("\\.gz$", "", File.Name))) %>%
  select(-sample_list) %>%  # 移除中間欄位
  distinct(Tumor.Sample.ID, Normal.Sample.ID, Basename)

# Merge with TLS labels
sample_sheet <- merge(x = files_dir, y = TLS_new_label, by.x = "Tumor.Sample.ID", by.y = "Sample.ID")
write.csv(sample_sheet, paste0("./data/", cancer_type, "/metadata/sample_sheet_mut.csv"), row.names = FALSE)
table(sample_sheet$TLS)


# 用 reticulate 呼叫 Python 的 numpy 來讀取
library(reticulate)
np <- import("numpy")
common_samples <- np$load(path.expand(paste0("~/TLS_meth/common_samples_", cancer_type)), allow_pickle = TRUE)
common_samples <- as.character(common_samples) # 轉成 R 的字串向量

# 篩選出 Tumor.Sample.ID 不在 complete cases 的 samples
sample_sheet <- sample_sheet[sample_sheet$Tumor.Sample.ID %in% common_samples, ] # only complete cases retained

############################################# MAFtools #################################################

maf_list <- lapply(sample_sheet$Basename, function(file) {
  tryCatch({
    maf <- read.maf(maf = file)
    
    # Check if non-synonymous mutation exists
    if (nrow(maf@data) > 0) {
      print(paste('Clinical Data Barcode :', maf@clinical.data$Tumor_Sample_Barcode))
      return(maf)
    } else { # 成功讀取 MAF 檔案，但沒有非同義突變
      warning(paste("Skipping file:", file, "No non-synonymous mutations found."))
      return(NULL)
    }
  }, error = function(e) { # 1. MAF檔案格式錯誤 2.檔案遺失或無法讀取 3.Variant_Classification欄位缺失
    warning(paste("Error reading file:", file, "Skipping..."))
    message("Reason: ", conditionMessage(e))
    return(NULL)
  })
})

# 看有多少檔案失敗
failed <- sum(sapply(maf_list, is.null))
cat("Failed or empty files:", failed, "/", length(maf_list), "\n")

# error_reasons <- sapply(sample_sheet$Basename, function(file) {
#   tryCatch({
#     maf <- read.maf(maf = file)
#     return("OK")
#   }, error = function(e) {
#     return(conditionMessage(e))
#   })
# })
# table(grepl("No non-synonymous mutations found", error_reasons))

# 移除 NULL，確保只有成功讀取的 MAF 檔案
maf_list <- Filter(Negate(is.null), maf_list) # Negate反轉TRUE/FALSE

# Merge MAF files
filtermaf = merge_mafs(mafs = maf_list, verbose = TRUE)
write.mafSummary(maf = filtermaf, basename = paste0("./data/", cancer_type, "/maf/filtermaf"))


## MAFtools
plotmafSummary(maf = filtermaf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# 準備 TLS 分組資訊
filtermaf@clinical.data$ShortID <- substr(filtermaf@clinical.data$Tumor_Sample_Barcode, 1, 16)
tls_df <- data.frame(ShortID = sample_sheet$Tumor.Sample.ID,
                     group = ifelse(sample_sheet$TLS == 1, "TLS-positive", "TLS-negative"), stringsAsFactors = FALSE)
filtermaf@clinical.data <- merge(filtermaf@clinical.data, tls_df, by = "ShortID", all.x = TRUE)

# 準備 topBarData（sample name + TMB）
tmb_result <- tmb(maf = filtermaf, captureSize = 35.8, plotType = NA)
tmb_data <- data.frame(tmb_result$Tumor_Sample_Barcode, tmb_result$total_perMB)
colnames(tmb_data) <- c("Tumor_Sample_Barcode", "TMB")

## Oncoplot (Add label TLS+/TLS-) 
tiff(filename = paste0("./tiff_figure/oncoplot_", cancer_type, ".tiff"), width = 8, height = 6, res = 300, 
     compression = "lzw", units = "in", bg = "white")
oncoplot(maf = filtermaf, top = 20, clinicalFeatures = "group", sortByAnnotation = TRUE,
         annotationColor = list(group = c("TLS-positive" = "#CC252D", "TLS-negative" = "#2065B4")),
         topBarData = tmb_data)
dev.off()

############################################# 匯出 metadata #################################################

sample_summary <- merge(x = filtermaf@variant.type.summary, y = filtermaf@variant.classification.summary, by = "Tumor_Sample_Barcode")
sample_summary$Tumor_Sample_Barcode <- sub("^([^-]+-[^-]+-[^-]+-[^-]+).*", "\\1", sample_summary$Tumor_Sample_Barcode)
write.csv(sample_summary, paste0("./data/", cancer_type, "/metadata/mutation_sample_metadata.csv"), row.names = FALSE)

gene_summary <- as.data.frame(filtermaf@gene.summary)
gene_silent <- as.data.frame(filtermaf@maf.silent)
write.csv(gene_summary, paste0("./data/", cancer_type, "/metadata/mutation_feature_metadata.csv"), row.names = FALSE)
write.csv(gene_silent, paste0("./data/", cancer_type, "/metadata/silentmutation_feature_metadata.csv"), row.names = FALSE)


library(dplyr)
gene_info_unique <- unique(as.data.table(filtermaf@data)[, .(Hugo_Symbol, BIOTYPE, Gene, Entrez_Gene_Id)])
gene_info <- gene_info_unique %>%
  group_by(Hugo_Symbol, Gene, Entrez_Gene_Id) %>%
  summarise(BIOTYPE = paste(unique(BIOTYPE), collapse = ";"), .groups = "drop")
write.csv(gene_info, paste0("./data/", cancer_type, "/metadata/mutation_geneinfo.csv"), row.names = FALSE)


######################################## 計算兩組的 TMB 是否有顯著差異 ##########################################

# 製作 TLS+ 和 TLS- 的 MAF files
maf_list_pos <- lapply(sample_sheet$Basename[sample_sheet$TLS == 1], function(file) {
  tryCatch({
    maf <- read.maf(maf = file)
    
    # Check if non-synonymous mutation exists
    if (nrow(maf@data) > 0) {
      print(paste('Clinical Data Barcode :', maf@clinical.data$Tumor_Sample_Barcode))
      return(maf)
    } else { # 成功讀取 MAF 檔案，但沒有非同義突變
      warning(paste("Skipping file:", file, "No non-synonymous mutations found."))
      return(NULL)
    }
  }, error = function(e) { # 1. MAF檔案格式錯誤 2.檔案遺失或無法讀取 3.Variant_Classification欄位缺失
    warning(paste("Error reading file:", file, "Skipping..."))
    return(NULL)
  })
})
# 移除 NULL，確保只有成功讀取的 MAF 檔案 (剩下56個)
maf_list_pos <- Filter(Negate(is.null), maf_list_pos) # Negate反轉TRUE/FALSE

# Merge MAF files
posmaf = merge_mafs(mafs = maf_list_pos, verbose = TRUE)
write.mafSummary(maf = posmaf, basename = paste0("./data/", cancer_type, "/maf/posmaf"))

maf_list_neg <- lapply(sample_sheet$Basename[sample_sheet$TLS == 0], function(file) {
  tryCatch({
    maf <- read.maf(maf = file)
    
    # Check if non-synonymous mutation exists
    if (nrow(maf@data) > 0) {
      print(paste('Clinical Data Barcode :', maf@clinical.data$Tumor_Sample_Barcode))
      return(maf)
    } else { # 成功讀取 MAF 檔案，但沒有非同義突變
      warning(paste("Skipping file:", file, "No non-synonymous mutations found."))
      return(NULL)
    }
  }, error = function(e) { # 1. MAF檔案格式錯誤 2.檔案遺失或無法讀取 3.Variant_Classification欄位缺失
    warning(paste("Error reading file:", file, "Skipping..."))
    return(NULL)
  })
})
# 移除 NULL，確保只有成功讀取的 MAF 檔案 (剩下106個)
maf_list_neg <- Filter(Negate(is.null), maf_list_neg) # Negate反轉TRUE/FALSE

# Merge MAF files
negmaf = merge_mafs(mafs = maf_list_neg, verbose = TRUE)
write.mafSummary(maf = negmaf, basename = paste0("./data/", cancer_type, "/maf/negmaf"))

## TMB 沒有顯著區別
TMB_positive <- tmb(maf = posmaf, captureSize = 35.8, plotType = NA)$total_perMB
TMB_negative <- tmb(maf = negmaf, captureSize = 35.8, plotType = NA)$total_perMB
TMB_pvalue <- wilcox.test(TMB_positive, TMB_negative, alternative = "two.sided", exact = FALSE)$p.value

# 準備資料
tmb_df <- data.frame(TMB = c(TMB_positive, TMB_negative),
                     Group = c(rep("TLS-positive", length(TMB_positive)), rep("TLS-negative", length(TMB_negative))))

# 決定 y 軸上限（例如排除最極端的 outlier）
library(ggplot2)
y_limit <- 3  # 調整這個數字

tiff(filename = paste0("./tiff_figure/TMB_", cancer_type, ".tiff"), width = 6, height = 6, res = 300, 
     compression = "lzw", units = "in", bg = "white")
ggplot(tmb_df, aes(x = Group, y = TMB, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.1) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 1.5, width = 0.5) +
  annotate("text", 
           x = 1.5,
           y = y_limit * 1.01,             # p value 位置跟著調整
           label = paste0("p = ", round(TMB_pvalue, 3)), 
           size = 4.5) +
  annotate("segment",
           x = 1, xend = 2,
           y = y_limit * 0.97, yend = y_limit* 0.97) +  # 畫橫線連接兩組
  scale_fill_manual(values = c("TLS-positive" = "#CC252D", "TLS-negative" = "#2065B4")) +
  coord_cartesian(ylim = c(0, y_limit)) +  # 限制顯示範圍但保留 outlier
  labs(x = NULL, y = "TMB (mutations/Mb)") +
  theme_bw() +                            
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


# 準備 TLS+ 分組資訊
posmaf@clinical.data$ShortID <- substr(posmaf@clinical.data$Tumor_Sample_Barcode, 1, 16)
tls_df <- data.frame(ShortID = sample_sheet$Tumor.Sample.ID[sample_sheet$TLS == 1],
                     group = "TLS-positive", stringsAsFactors = FALSE)
posmaf@clinical.data <- merge(posmaf@clinical.data, tls_df, by = "ShortID", all.x = TRUE)

## Oncoplot (TLS+) 
oncoplot(maf = posmaf, top = 10, clinicalFeatures = "group", sortByAnnotation = TRUE, 
         annotationColor = list(group = c("TLS-positive" = "#E41A1C","TLS-negative" = "#377EB8")))

# 準備 TLS- 分組資訊
negmaf@clinical.data$ShortID <- substr(negmaf@clinical.data$Tumor_Sample_Barcode, 1, 16)
tls_df <- data.frame(ShortID = sample_sheet$Tumor.Sample.ID[sample_sheet$TLS == 0],
                     group = "TLS-negative", stringsAsFactors = FALSE)
negmaf@clinical.data <- merge(negmaf@clinical.data, tls_df, by = "ShortID", all.x = TRUE)
## Oncoplot (TLS-) 
oncoplot(maf = negmaf, top = 10, clinicalFeatures = "group", sortByAnnotation = TRUE, 
         annotationColor = list(group = c("TLS-positive" = "#E41A1C","TLS-negative" = "#377EB8")))

#getGeneSummary(posmaf)
#getGeneSummary(negmaf)

######################################## Create mutation count matrix ##########################################

# 建立 mutation count matrix
all_genes <- unique(filtermaf@data$Hugo_Symbol)
mat <- maftools:::createOncoMatrix(m = filtermaf, g = all_genes)
mut_num_matrix <- mat$numericMatrix  # 以 filtermaf@clinical.data$Tumor_Sample_Barcode 對欄位命名
colnames(mut_num_matrix) <- sub("^([^-]+-[^-]+-[^-]+-[^-]+).*", "\\1", colnames(mut_num_matrix))

# createOncoMatrix() 會把完全沒有突變的 sample 移除，因此需補回被移除的 samples
all_samples <- sample_sheet$Tumor.Sample.ID
rm_samples <- setdiff(all_samples, colnames(mut_num_matrix)) 
# 填入突變資料
filled_num_matrix <- matrix(0, nrow = length(all_genes), ncol = length(all_samples),
                            dimnames = list(all_genes, all_samples))
filled_num_matrix[rownames(mut_num_matrix), colnames(mut_num_matrix)] <- mut_num_matrix
write.csv(filled_num_matrix, paste0("./data/", cancer_type, "/raw_matrix/mut_num_matrix_tomodel.csv"), row.names = TRUE)

######################################## Create mutation binary matrix ##########################################

# 建立 binary matrix
all_genes <- unique(filtermaf@data$Hugo_Symbol)
mat <- maftools:::createOncoMatrix(m = filtermaf, g = all_genes)
mut_matrix <- mat$numericMatrix  # 以 filtermaf@clinical.data$Tumor_Sample_Barcode 對欄位命名
colnames(mut_matrix) <- sub("^([^-]+-[^-]+-[^-]+-[^-]+).*", "\\1", colnames(mut_matrix))
mut_matrix[mut_matrix != 0] <- 1

# createOncoMatrix() 會把完全沒有突變的 sample 移除，因此需補回被移除的 samples
all_samples <- sample_sheet$Tumor.Sample.ID
rm_samples <- setdiff(all_samples, colnames(mut_matrix))
# 填入突變資料
filled_matrix <- matrix(0, nrow = length(all_genes), ncol = length(all_samples),
                        dimnames = list(all_genes, all_samples))
filled_matrix[rownames(mut_matrix), colnames(mut_matrix)] <- mut_matrix
write.csv(filled_matrix, paste0("./data/", cancer_type, "/raw_matrix/mut_matrix_tomodel.csv"), row.names = TRUE)
