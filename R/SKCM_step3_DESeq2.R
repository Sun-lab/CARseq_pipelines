
library(DESeq2)

# ------------------------------------------------------------------------
# read in bulk RNA-seq data
# ------------------------------------------------------------------------

counts = readRDS("../data/TCGA_SKCM_raw_counts.rds")
dim(counts)
counts[1:2,1:5]

# remove genes with low read counts: 
# genes with 75% percentile of read counts rd75 < 20 are removed
rd75   = apply(counts, 1, function(x) quantile(x, 0.75))
counts = counts[rd75 >= 20, ]
counts = round(counts)

dim(counts)
counts[1:2,1:5]

# ------------------------------------------------------------------------
# read in covariate data
# ------------------------------------------------------------------------

col_data = readRDS("../data/SKCM_cavariates.rds")
class(col_data)

dim(col_data)
col_data[1:2,]

is.factor(col_data$tss)

colData = col_data
for(i in 1:ncol(colData)){
  if(is.character(colData[[i]])){
    colData[[i]] = as.factor(colData[[i]])
  }
}
dim(colData)
colData[1:2,]
summary(colData)

# ------------------------------------------------------------------------
# Run DESeq2
# ------------------------------------------------------------------------

countData = counts
colnames(countData) = gsub(".", "-", substr(colnames(counts), 1, 12), fixed=T)
dim(countData)
countData[1:2,1:3]
countData = countData[,match(colData$bcr_patient_barcode, colnames(countData))]
dim(countData)

dds = DESeqDataSetFromMatrix(countData = countData, 
                             colData = colData,
                             design = ~ gender + scaled_age + 
                               scaled_log_depth + stage + tss + Survival)
dds = DESeq(dds)

res = results(dds)
dim(res)
head(res)
summary(res)

pdf("SKCM_pval_hist.pdf", width=3, height=3)
par(mar=c(5,4,1,1), bty="n")
hist(res$pvalue, main="", xlab="p-value")
dev.off()

saveRDS(res, file="../results/SKCM_DESeq2.rds")

gc()

sessionInfo()
q(save="no")
