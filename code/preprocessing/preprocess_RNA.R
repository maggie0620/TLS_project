cancer_type = "PAAD"

# TLS label
library(readxl)
library(tidyverse)
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
gdc_data <- read.delim(paste0("./data/", cancer_type, "/metadata/gdc_sample_sheet_RNA.tsv"), header = TRUE, sep = "\t")

# Get file path
if (cancer_type == "CHOL"){
  files_dir <- gdc_data %>%
    filter(grepl("gene_counts.tsv", File.Name)) %>%
    mutate(Basename = paste0("~/TLS_meth/DATA_", cancer_type, "_RNA/", File.Name)
    ) %>%
    distinct(Case.ID, Sample.ID, Sample.Type, Basename)
} else {
  files_dir <- gdc_data %>%
    filter(grepl("gene_counts.tsv", File.Name)) %>%
    mutate(Basename = paste0("~/TLS_meth/DATA_", cancer_type, "_RNA/", File.Name)
    ) %>%
    distinct(Case.ID, Sample.ID, Tissue.Type, Tumor.Descriptor, Basename)
}

# Merge with TLS labels
sample_sheet <- merge(x = files_dir, y = TLS_new_label, by = "Sample.ID")
write.csv(sample_sheet, paste0("./data/", cancer_type, "/metadata/sample_sheet_RNA.csv"), row.names = FALSE)
table(sample_sheet$TLS)

# Create matrix with gene type (V4: unstranded, V7: tpm_unstranded)
file <- read.csv(sample_sheet$Basename[1], header = F, sep='\t', skip=6)
unique(file$V3)
gene_ids <- as.character(file$V1)
data_tpm <- file$V7

for (i in 2:nrow(sample_sheet)) {
  file <- read.csv(sample_sheet$Basename[i], header = F, sep='\t', skip=6)
  data_tpm <- cbind(data_tpm, file$V7)
}

data_tpm <- as.data.frame(data_tpm)
data_tpm <- log2(data_tpm + 1)
new_data_tpm <- cbind(data_tpm, file$V2, file$V3) # 直接使用 gencode v36的 gene name
colnames(new_data_tpm) <- c(sample_sheet$Sample.ID, "gene_name","gene_type")
rownames(new_data_tpm) <- as.character(gene_ids)
write.csv(new_data_tpm, paste0("./data/", cancer_type, "/raw_matix/RNA_tomodel.csv"), row.names = TRUE)

###############################################################################################

# Protein-coding genes for DESeq2
file <- read.csv(sample_sheet$Basename[1], header = F, sep='\t', skip=6)
unique(file$V3)
new_file <- file %>% filter(grepl("protein_coding", V3))   # 只剩19962個gene
gene_ids <- as.character(new_file$V1)
data_raw_pc <- new_file$V4

for (i in 2:nrow(sample_sheet)) {
  file <- read.csv(sample_sheet$Basename[i], header = F, sep='\t', skip=6)
  new_file <- file %>% filter(grepl("protein_coding", V3))   # 只剩19962個gene
  data_raw_pc <- cbind(data_raw_pc, new_file$V4)
}

data_raw_pc <- as.data.frame(data_raw_pc)
colnames(data_raw_pc) <- sample_sheet$Sample.ID
rownames(data_raw_pc) <- gene_ids
#saveRDS(data_raw_pc, file = paste0("./RDS/", cancer_type, "/RNA/data_raw_pc.rds"))


## Perform DESeq2
library("DESeq2")
library(ggplot2)
sample_sheet$TLS <- factor(sample_sheet$TLS)  # 若 TLS 是數值（不是 factor），那 DESeq2 會認為 TLS 從 0 到 1 是一種「強度遞增」的數值效應
dds <- DESeqDataSetFromMatrix(countData = data_raw_pc, 
                              colData = sample_sheet, 
                              design = ~TLS)
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))

# Summary of differential gene expression
summary(res)
# Sort summary list by p-value
res <- res[order(res$padj),]
head(res)

# plotCounts
par(mfrow=c(2,3))
for (i in 1:6){
  plotCounts(dds, gene=rownames(res)[i], intgroup="TLS")
}

## Volcano plot (用log2FC = 1作為標準)
par(mfrow=c(1,1))
res_df <- as.data.frame(res)
res_df$change <- ifelse(res_df$padj < 10^-2 & abs(res_df$log2FoldChange) > 1,
                        ifelse(res_df$log2FoldChange > 1,"UP","DOWN"),"NOT")
table(res_df$change)
ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(alpha=0.4, size=1) +
  theme_set(theme_set(theme_gray(base_size = 10))) +
  xlab("log2(FC)") +
  ylab("-log10(p.adj)") +
  theme(plot.title = element_text(size=10, hjust=0.5)) +
  scale_color_manual(values = c("blue","black","red"))

write.csv(res_df, paste0("./data/", cancer_type, "/metadata/res_df_pc.csv"), row.names = TRUE)
