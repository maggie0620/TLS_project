library(minfi)
library(IlluminaHumanMethylation450kmanifest)
# library(minfiData)
# library("ChAMP")
# #library("ChAMPdata")
library(devtools)
library(maxprobes)
# library(shinyMethyl)
# library(doParallel)
# library(ENmix)
library(RColorBrewer)
# library(openxlsx)
library(VennDiagram)
# library(tidyverse)
# library(qqman)
library(gtools)
# library("ggplot2")
# library(DMRcate)
# library(Gviz)
library("missMethyl")

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

# Change sample sheet.tsv into .csv (450K)
gdc_data <- read.delim(paste0("./data/", cancer_type, "/metadata/gdc_sample_sheet_meth.tsv"), header = TRUE, sep = "\t")

# Keep only IDAT files
idat_files <- gdc_data %>%
  filter(grepl("idat", File.Name)) %>%
  mutate(Sentrix_ID = sub("_.*", "", File.Name),                   # 抽取 Sentrix_ID
         Sentrix_Position = sub(".*_", "", sub("_(Grn|Red)\\.idat", "", File.Name)),  # 抽取 Sentrix_Position
         Basename = paste0("~/TLS_meth/DATA_", cancer_type, "/", Sentrix_ID, "_", Sentrix_Position)
  ) %>%
  distinct(Case.ID, Sample.ID, Sample.Type, Sentrix_ID, Sentrix_Position, Basename)

# Merge with TLS labels
sample_sheet <- merge(x = idat_files, y = TLS_new_label, by = "Sample.ID")
write.csv(sample_sheet, paste0("./data/", cancer_type, "/metadata/sample_sheet_meth.csv"), row.names = FALSE)
table(sample_sheet$TLS)

#############################################################################################

## Information of the sample sheet
RGset_extend <- read.metharray.exp(targets = sample_sheet, force = TRUE, extended = TRUE)
saveRDS(RGset_extend, file = paste0("./RDS/", cancer_type, "/RGset_extend.rds"))
#RGset_extend <- readRDS(paste0("./RDS/", cancer_type, "/RGset_extend.rds"))
phenodata <- pData(RGset_extend)
mani <- getManifest(RGset_extend)
#head(getProbeInfo(mani))
#annotation(RGset_extend)  #使用的是Illumina Human Methylation 450K

## Filter probes with <3 beads in more than 5% of samples per probe (Remove necessary control probes?)
few_beads <- getNBeads(RGset_extend) < 3
keep_probe <- rowMeans(few_beads) <= 0.05
table(keep_probe)  # 排除FALSE
#RGset <- RGset_extend[keep_probe,]

## Plot Control QC
#plotCtrl(RGset)
qcReport(RGset, pdf= "qcReport.pdf")

## Preprocessing by preprocessFunnorm
GMset.funnorm <- preprocessFunnorm(RGset_extend, ratioConvert = F) #Background and dye bias correction with noob

## 過濾品質較差的sample (badSampleCutoff = 10.5)
GMset_normQC <- minfiQC(GMset.funnorm, fixOutliers = T, verbose = T)
tiff(filename = paste0("./tiff_figure/sampleQC_", cancer_type, ".tiff"), width = 6, height = 5, res = 300, 
     compression = "lzw", units = "in", bg = "transparent")
plotQC(GMset_normQC$qc) # No bad samples
dev.off()

# 取得 RGSet 的 annotation 資訊 (Name 會儲存標準 Illumina CpG ID)
anno <- getAnnotation(RGset_extend)
length(anno$Name)
# Apply <3 beads filtering to GMset.funnorm
valid_probes <- intersect(anno$Name[keep_probe], rownames(GMset.funnorm))
length(valid_probes)
GMset.funnorm.filtered <- GMset.funnorm[valid_probes, ]  # 過濾<3 beads的probes

## 計算detection p-values (需要使用完整原始的原始訊號強度資料!!!!!)
detP <- detectionP(RGset_extend)
head(detP)
saveRDS(detP, file = paste0("./RDS/", cancer_type, "/detP.rds"))
## 過濾可信度較低的probe (detection p-value > 0.01，若>10% samples中failed會被移除)
keep_probe <- rowMeans(detP > 0.01) <= 0.1  # 同一個probe在所有sample中failed的比例
table(keep_probe)  # 排除FALSE

# Apply detP > 0.01 filtering to GMset.funnorm.filtered
valid_probes <- intersect(anno$Name[keep_probe], rownames(GMset.funnorm.filtered))
length(valid_probes)  
GMset.funnorm.filtered <- GMset.funnorm.filtered[valid_probes, ]  # 過濾detP > 0.01的probes
saveRDS(GMset.funnorm.filtered, file = paste0("./RDS/", cancer_type, "/GMset_funnorm_filtered.rds"))  #GenomicMethylSet

## Filter all SNP-related probes
GRset <- ratioConvert(GMset.funnorm.filtered, what = "both", keepCN = TRUE)  #GenomicRatioSet
#snp <- getSnpInfo(GRset)
GRset <- addSnpInfo(GRset)
GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0) # maf為minor allele frequencies
saveRDS(GRset, file = paste0("./RDS/", cancer_type, "/GRset_pre.rds"))

## Filter cross-reactive probes (搬到筆電 RStudio 4.3.2 進行)
GRset <- dropXreactiveLoci(GRset)
GRset <- readRDS(paste0("./RDS/", cancer_type, "/GRset.rds"))

## Probes located in chromosome X and Y
anno = getAnnotation(GRset)
XY_probe = rownames(anno)[anno$chr %in% c("chrX","chrY")]
keep = !(featureNames(GRset) %in% XY_probe)
GRset = GRset[keep,]
saveRDS(GRset, file = paste0("./RDS/", cancer_type, "/GRset.rds"))
anno_new <- getAnnotation(GRset)
write.csv(anno_new, paste0(cancer_type, "_meth_anno_tomodel.csv"), row.names = TRUE)

GRset <- readRDS(paste0("./RDS/", cancer_type, "/GRset.rds"))
## Get beta-value
beta <- getBeta(GRset)
head(beta)
saveRDS(beta, file = paste0("./RDS/", cancer_type, "/beta.rds"))
colnames(beta) <- sample_sheet$Sample.ID
beta_matrix <- as.data.frame(beta)
write.csv(beta_matrix, paste0("./data/", cancer_type, "/raw_matrix/meth_tomodel.csv"), row.names = TRUE)

# Check if beta value has any 0
if (any(beta == 0)) {
  cat("Beta has at least one zero value.\n")
} else {
  cat("Beta does not have any zero values.\n")
}
zero_positions <- which(beta == 0, arr.ind = TRUE)
length(zero_positions)
na_positions <- which(is.na(beta), arr.ind = TRUE)
length(na_positions)

# Density Plot
tiff(filename = paste0("./tiff_figure/beta_densityplot_", cancer_type, ".tiff"), width = 8, height = 6, res = 300, 
     compression = "lzw", units = "in", bg = "transparent")
densityPlot(beta, sampGroups = sample_sheet$TLS, main = cancer_type, legend = FALSE, pal = c("#2065B4", "#CC252D"))
legend("topright", legend = c("TLS−","TLS+"), text.col=c("#2065B4", "#CC252D"), fill=c("#2065B4", "#CC252D"), cex=1)
dev.off()

# Get M-value
m_value <- getM(GRset)
head(m_value)
saveRDS(m_value, file = paste0("./RDS/", cancer_type, "/m_value.rds"))

###########################################################################################

# Raw Data對應methylation locus (未經normalization)
pre_Mset <- preprocessRaw(RGset_extend)  #MethylSet

## 比較normalization前後 beta value的差異 (by ProbeType)
par(mfrow=c(1,2))
plotBetasByType(pre_Mset[,1], main="Raw", legendPos = "top", cex=0.55)
# GMset不含ProbeInfo 須從Mset取出
typeI <- getProbeInfo(pre_Mset, type = "I")[, c("Name","nCpG")]
typeII <- getProbeInfo(pre_Mset, type = "II")[, c("Name","nCpG")]
probeTypes <- rbind(typeI, typeII)
probeTypes$Type <- rep(x = c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
plotBetasByType(beta[,1], probeTypes = probeTypes, main="Normalization", legendPos = "top", cex=0.55)

## 比較normalization前後 beta value的差異 (by SampleGroup)
densityPlot(RGset_extend, sampGroups=sample_sheet$TLS, main="Raw", legend=FALSE)
legend("top", legend = c("TLS-","TLS+"), text.col=brewer.pal(8,"Dark2"), fill=brewer.pal(8,"Dark2"), cex=0.6)
densityPlot(beta, sampGroups=sample_sheet$TLS, main="Normalization", legend=FALSE)
legend("top", legend = c("TLS-","TLS+"), text.col=brewer.pal(8,"Dark2"), fill=brewer.pal(8,"Dark2"), cex=0.6)

## MDS Plot
# Raw
mdsPlot(getBeta(RGset_extend), numPositions = 1000, sampGroups = phenodata$TLS, pch=16, legendPos = NA)
legend("bottomright", legend = c("TLS-","TLS+"), text.col=brewer.pal(8,"Dark2"), fill=brewer.pal(8,"Dark2"), cex=0.6)
# Norm & filter
mdsPlot(beta, numPositions = 1000, sampGroups = phenodata$TLS, pch=16, legendPos = NA)
legend("bottomright", legend = c("TLS-","TLS+"), text.col=brewer.pal(8,"Dark2"), fill=brewer.pal(8,"Dark2"), cex=0.6)

###########################################################################################

## Differentially Methylated Positions (DMPs)
# 使用Beta value：在預測方面可能會更好 (m-value會取log造成部分值特別大)
dmp <- dmpFinder(beta, pheno = phenodata$TLS, type = "categorical") # sorted by p-value
head(dmp)
length(dmp[dmp$qval < 0.05,1])
#saveRDS(dmp, file = paste0("./RDS/", cancer_type, "/dmp.rds"))

dmp_m <- dmpFinder(m_value, pheno = phenodata$TLS, type = "categorical")
head(dmp_m)
length(dmp_m[dmp_m$qval < 0.05,1])
#saveRDS(dmp_m, file = paste0("./RDS/", cancer_type, "/dmp_m.rds"))

#重疊率(79635)
dmpB <- rownames(dmp[dmp$qval < 0.05,])
dmpM <- rownames(dmp_m[dmp_m$qval < 0.05,])
ven <- venn.diagram(x=list(dmpM=dmpM, dmpB=dmpB), filename = NULL, fill = c("#F8756D","#01bec4"))
cowplot::plot_grid(ven)

## Plot DMPs
cpgs <- rownames(dmp)[1:6]
par(mfrow = c(2,3))
plotCpg(beta, cpg=cpgs, pheno = phenodata$TLS, ylab = "Beta")

# TLS+與TLS-的差異分析
TLS_pos_ind <- which(sample_sheet$TLS == 1, arr.ind = TRUE)
TLS_neg_ind <- which(sample_sheet$TLS == 0, arr.ind = TRUE)
# Volcano plot (Fold Change版本)
beta_pos <- beta[which(row.names(beta) %in% rownames(dmp)), TLS_pos_ind]
beta_neg <- beta[which(row.names(beta) %in% rownames(dmp)), TLS_neg_ind]
FC <- foldchange(rowMeans(beta_pos), rowMeans(beta_neg))
logFC <- foldchange2logratio(FC, base = 2)  # base為log底數
logFC_order <- logFC[order(names(logFC))]
dmp_order <- dmp[order(row.names(dmp)),]
dmp_order$logFC <- logFC_order

dmp_order$change <- ifelse(dmp_order$qval < 10^-2 & abs(dmp_order$logFC) > 1,
                           ifelse(dmp_order$logFC > 1,"Hyper","Hypo"),"NOT")
table(dmp_order$change)
ggplot(data = dmp_order, aes(x = logFC, y = -log10(qval), color = change)) +
  geom_point(alpha=0.4, size=1) +
  theme_set(theme_set(theme_gray(base_size = 10))) +
  xlab("log2(FC)") +
  ylab("-log10(P-value)") +
  theme(plot.title = element_text(size=10, hjust=0.5)) +
  scale_color_manual(values = c("red","blue","black"))


# Volcano plot (delta beta版本: 用 0.1 作為標準)
delbeta <- rowMeans(beta_pos) - rowMeans(beta_neg)
delbeta_order <- delbeta[order(names(delbeta))]
dmp_order$delta_beta <- delbeta_order

dmp_order$changeBeta <- ifelse(dmp_order$qval < 10^-2 & abs(dmp_order$delta_beta) > 0.1,
                               ifelse(dmp_order$delta_beta > 0.1,"Hyper","Hypo"),"NOT")
table(dmp_order$changeBeta)
ggplot(data = dmp_order, aes(x = delta_beta, y = -log10(qval), color = changeBeta)) +
  geom_point(alpha=0.4, size=1) +
  theme_set(theme_set(theme_gray(base_size = 10))) +
  xlab("delta beta") +
  ylab("-log10(P-value)") +
  theme(plot.title = element_text(size=10, hjust=0.5)) +
  scale_color_manual(values = c("red","blue","black"))

write.csv(dmp_order , paste0("./data/", cancer_type, "/metadata/meth_dmp_tomodel.csv"), row.names = TRUE)
saveRDS(dmp_order, file = paste0("./RDS/", cancer_type, "/dmp_order.rds"))
table(dmp_order$changeBeta)


## RefGene in DMPs
dmp_sig <- row.names(dmp_order[dmp_order$changeBeta %in% c("Hyper","Hypo"),])
length(dmp_sig)  # Total number of significant CpGs
anno <- getAnnotation(GRset)
anno_pick <- anno[,c("chr","pos","Relation_to_Island","UCSC_RefGene_Name","UCSC_RefGene_Group","UCSC_RefGene_Accession")]
dmp_anno <- as.data.frame(subset(anno_pick, row.names(anno_pick) %in% c(dmp_sig)))
saveRDS(dmp_anno, file = paste0("./RDS/", cancer_type, "/dmp_anno.rds"))
