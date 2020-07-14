
# compare the results of DESeq2 with or without cell type proportion 
# (log ratios) as covaraites

library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(CARseq)
library(DESeq2)
theme_set(theme_classic2())

# -------------------------------------------------------------------------
# compare DESeq2 analysis results
# -------------------------------------------------------------------------

res1 = fread("../results/SCZ_step1_DESeq2_bulk_adj_covariates_SV2.txt")

dim(res1)
res1[1:2,]

res2 = fread("../results/SCZ_step1_DESeq2_bulk_adj_covariates_SV2_log_ct_fractions.txt")

dim(res2)
res2[1:2,]

res1$qvalue = get_qvalues_one_inflated(res1$pvalue)
res2$qvalue = get_qvalues_one_inflated(res2$pvalue)

table(res1$V1 == res2$V1)

sort(res1$pvalue)[1:5]
sort(res2$pvalue)[1:5]

sort(res1$qvalue)[1:5]
sort(res2$qvalue)[1:5]

summary(res1$pvalue)
summary(res2$pvalue)

summary(res1$qvalue)
summary(res2$qvalue)

res.df = merge(res1, res2, by="V1", suffixes = c(".noCT", ".wCT"))
dim(res.df)
res.df[1:2,]

table(res.df$qvalue.noCT < 0.1)
table(res.df$qvalue.wCT  < 0.1)

table(res.df$qvalue.noCT < 0.05, res.df$qvalue.wCT < 0.05)
table(res.df$qvalue.noCT < 0.1, res.df$qvalue.wCT  < 0.1)

cut1 = 1e-7
res.df$qvalue.noCT[which(res.df$qvalue.noCT < cut1)] = cut1
res.df$qvalue.wCT[which(res.df$qvalue.wCT < cut1)]   = cut1

g1 = ggplot(res.df, aes(x=-log10(qvalue.noCT),y=-log10(qvalue.wCT))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(qvalue no CT adj.)", y= "-log10(qvalue with CT adj.)", 
       title="SCZ vs. control") + xlim(0, 7) + ylim(0, 7) + 
  geom_abline(intercept = 0, slope=1)

summary(res.df$log2FoldChange.noCT)
summary(res.df$log2FoldChange.wCT)

g2 = ggplot(res.df, aes(x=log2FoldChange.noCT,
                        y=log2FoldChange.wCT)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "log2 fc, no CT adj.", y= "log2 fc, with CT adj.)", 
       title="SCZ vs. control") + 
  geom_abline(intercept = 0, slope=1)

pdf("../figures/SCZ_DESeq2_comparison_log10pval.pdf", width=4, height=3.2)
print(g1)
dev.off()

pdf("../figures/SCZ_DESeq2_comparison_log2fc.pdf", width=4, height=3.2)
print(g2)
dev.off()

# read files copied from "CARseq_pipelines" repository
trec = readRDS("../data/trec_filtered_scz_control.rds")
clinical_variables_raw = readRDS("../data/dat_cavariates_scz_control_with_svs.rds")

# for code compatibility:
col_data = as.data.frame(clinical_variables_raw)
col_data$Dx = ifelse(col_data$DxSCZ == 1, "SCZ", "Control")
col_data$RNAseq.Sample_RNA_ID = colnames(trec)
col_data$Diagnosis = col_data$Dx

# rescale age_death, PMI, RIN, RIN^2, log_depth
col_data$scaled_age_death = scale(col_data$age_death)
col_data$scaled_log_depth = scale(col_data$log_depth)
col_data$scaled_PMI = scale(col_data$PMI)
col_data$scaled_RIN = scale(col_data$RIN)
col_data$scaled_RIN2 = scale(col_data$RIN2)

# rescale PCs and SVs -- the names are still sv1, sv2, etc, but they have been scaled:
covariates_to_scale = grep("sv|PC", colnames(col_data))
for (m in covariates_to_scale) {
  col_data[, m] = scale(col_data[, m])
}

# permuted disease label
set.seed(1234)
col_data$Diagnosis = col_data$Dx
col_data$DiagnosisP = CARseq:::permute_case_and_controls(col_data$Diagnosis)
table(col_data$Diagnosis, col_data$DiagnosisP)

# load gene annotation:
rse_filtered = readRDS("../data/rse_filtered_SV.rds")

# cellular proportions
data_folder = "../data"
prop_output_file = "../data/SCZ_prop.rds"
prop_list = readRDS(prop_output_file)

prop_list$ICeDT = prop_list$ICeDT[
  match(col_data$RNAseq.Sample_RNA_ID, rownames(prop_list$ICeDT)), ]
prop_list$CIBERSORT = prop_list$CIBERSORT[
  match(col_data$RNAseq.Sample_RNA_ID, rownames(prop_list$CIBERSORT)), ]
#saveRDS(prop_list, file="../data/prop_resized.rds")

cell_sizes = readRDS("../MTG/cell_sizes_MTG.rds")
cell_sizes = cell_sizes[c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")]
icedt_rho = prop_list$ICeDT
cibersort_rho = prop_list$CIBERSORT

# -------------------------------------------------------------------------
# read in data
# -------------------------------------------------------------------------

# total read count
trec = readRDS("../data/trec_filtered_scz_control.rds")  # 20788x527 matrix
clinical_variables_raw = readRDS("../data/dat_cavariates_scz_control_with_svs.rds")

dim(clinical_variables_raw)
clinical_variables_raw[1:2,]

# for code compatibility:
col_data = as.data.frame(clinical_variables_raw)
col_data$Dx = ifelse(col_data$DxSCZ == 1, "SCZ", "Control")
col_data$RNAseq.Sample_RNA_ID = colnames(trec)

# rescale age_death, PMI, RIN, RIN^2, log_depth
col_data$scaled_age_death = scale(col_data$age_death)
col_data$scaled_log_depth = scale(col_data$log_depth)
col_data$scaled_PMI = scale(col_data$PMI)
col_data$scaled_RIN = scale(col_data$RIN)
col_data$scaled_RIN2 = scale(col_data$RIN2)

# rescale PCs and SVs -- the names are still sv1, sv2, etc, but they have been scaled:
covariates_to_scale = grep("sv|PC", colnames(col_data))
for (m in covariates_to_scale) {
  col_data[, m] = scale(col_data[, m])
}

# cellular proportions
prop_output_file = "../data/SCZ_prop.rds"
prop_list = readRDS(prop_output_file)

prop_list$ICeDT = prop_list$ICeDT[
  match(col_data$RNAseq.Sample_RNA_ID, rownames(prop_list$ICeDT)), ]
prop_list$CIBERSORT = prop_list$CIBERSORT[
  match(col_data$RNAseq.Sample_RNA_ID, rownames(prop_list$CIBERSORT)), ]

icedt_rho = prop_list$ICeDT
cibersort_rho = prop_list$CIBERSORT

dim(icedt_rho)
icedt_rho[1:2,]

rownames(col_data) = col_data$RNAseq.Sample_RNA_ID
table(rownames(col_data) == rownames(icedt_rho))

for(i in 1:ncol(col_data)){
  if(is.character(col_data[[i]])){
    col_data[[i]] = as.factor(col_data[[i]])
  }
}

col_data = cbind(col_data, log(icedt_rho + 0.01))
dim(col_data)
col_data[1:2,]

# -------------------------------------------------------------------------
# check whether those genes with smaller p-values are associated with
# cell type compositions
# -------------------------------------------------------------------------

crs = abs(cor(t(trec), icedt_rho, method = "spearman"))
dim(crs)
crs[1:2,1:3]

crs.max  = apply(crs, 1, max)
crs.wmax = apply(crs, 1, which.max)

table(crs.wmax)

table(res1$qvalue < 0.1)
table(res1$qvalue < 0.05)

summary(crs.max[which(res1$qvalue < 0.05)])
summary(crs.max[which(res1$qvalue < 0.1)])
summary(crs.max[which(res1$qvalue >= 0.1)])

# -------------------------------------------------------------------------
# Re-run DESeq2 given log transformed cell type composition of Exc
# -------------------------------------------------------------------------

dim(trec)
trec[1:2,1:3]

dim(col_data)
col_data[1:3,]

col_data$Diagnosis = col_data$Dx
table(col_data$Diagnosis)

dd3 = DESeqDataSetFromMatrix(countData = trec, 
                             colData = col_data,
                             design = ~ scaled_log_depth + InstitutionPenn + InstitutionPitt +
                               genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 +
                               genoPC1 + genoPC2 + libclustB + libclustbase + libclustC + libclustD + 
                               libclustE + libclustF + libclustG + Exc + sv1 + sv2 + Diagnosis)
date()
dd3 = DESeq(dd3, parallel=TRUE)
date()

res3 = results(dd3)
dim(res3)
res3[1:2,]

res3 = as.data.frame(res3)

res3$qvalue = get_qvalues_one_inflated(res3$pvalue)
table(rownames(res3) == res2$V1)

# -------------------------------------------------------------------------
# Re-run DESeq2 given log transformed cell type composition of Inh
# -------------------------------------------------------------------------

dd4 = DESeqDataSetFromMatrix(countData = trec, 
                             colData = col_data,
                             design = ~ scaled_log_depth + InstitutionPenn + InstitutionPitt +
                               genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 +
                               genoPC1 + genoPC2 + libclustB + libclustbase + libclustC + libclustD + 
                               libclustE + libclustF + libclustG + Inh + sv1 + sv2 + Diagnosis)
date()
dd4 = DESeq(dd4, parallel=TRUE)
date()

res4 = results(dd4)
dim(res4)
res4[1:2,]

res4 = as.data.frame(res4)

res4$qvalue = get_qvalues_one_inflated(res4$pvalue)
table(rownames(res4) == res2$V1)

table(res1$qvalue < 0.1)
table(res2$qvalue < 0.1)
table(res3$qvalue < 0.1)
table(res4$qvalue < 0.1)

table(res1$qvalue < 0.05, res2$qvalue < 0.05)
table(res1$qvalue < 0.05, res3$qvalue < 0.05)
table(res1$qvalue < 0.05, res4$qvalue < 0.05)
table(res2$qvalue < 0.05, res3$qvalue < 0.05)
table(res2$qvalue < 0.05, res4$qvalue < 0.05)
table(res3$qvalue < 0.05, res4$qvalue < 0.05)

write.table(res3, file="../results/SCZ_step1_DESeq2_bulk_adj_covariates_SV2_log_Exc.txt",  
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(res4, file="../results/SCZ_step1_DESeq2_bulk_adj_covariates_SV2_log_Inh.txt",  
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

gc()

sessionInfo()
q(save="no")


