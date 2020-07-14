
library(CARseq)
library(SummarizedExperiment)
library(matrixStats)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(tidyr)
theme_set(theme_bw())
library(tidyverse)
library(stargazer)

rse_filtered = readRDS("../data/ASD_rse_filtered_with_SVs.rds")
cellular_proportions = readRDS("../data/ASD_prop.rds")

set.seed(1234)
assay(rse_filtered, "counts") = round(assay(rse_filtered, "counts"))

trec = assays(rse_filtered)$counts
col_data = colData(rse_filtered)

dim(col_data)
col_data[1:2,]

dim(trec)
trec[1:2,1:3]

table(col_data$BrainBank)
table(col_data$CellType)
table(col_data$Sequencing.Batch)

col_data$scaled_AgeDeath = scale(col_data$AgeDeath)
col_data$scaled_log_depth = scale(col_data$log_depth)
col_data$scaled_RIN = scale(col_data$RIN)
for (i in 1:8) {
  col_data[[paste0("scaled_sv", i)]] = scale(col_data[[paste0("sv", i)]])
}

cor(model.matrix(~ BrainBank + Sequencing.Batch + scaled_AgeDeath + 
                   scaled_log_depth + scaled_RIN + scaled_sv1 + scaled_sv2 +  
                   scaled_sv3 + scaled_sv4 + Diagnosis, col_data)[,-1])
col_data$DiagnosisP = CARseq:::permute_case_and_controls(col_data$Diagnosis)
table(col_data$Diagnosis, col_data$DiagnosisP)


# ------------------------------------------------------------------------
# plot CIBERSORT vs. ICeDT cellular frequency estimates
# ------------------------------------------------------------------------

sapply(cellular_proportions, dim)

icedt_rho = cellular_proportions$ICeDT
cibersort_rho = cellular_proportions$CIBERSORT

dim(icedt_rho)
dim(cibersort_rho)
icedt_rho[1:2,]
cibersort_rho[1:2,]
table(rownames(icedt_rho) == rownames(cibersort_rho))

icedt_rho_df = data.frame(icedt_rho)
cibersort_rho_df = data.frame(cibersort_rho)

icedt_rho_df$sample = rownames(icedt_rho)
cibersort_rho_df$sample = rownames(cibersort_rho)

rho_df = merge(icedt_rho_df, cibersort_rho_df, by="sample", 
               suffixes=c(".ICeDT", ".CIBERSORT"))
dim(rho_df)
rho_df[1:2,]

rhos = rho_df %>% pivot_longer(-sample, names_to=c("cell_type", "method"), 
                               names_sep="\\.") 

dim(rhos)
rhos[1:2,]

summary(rho_df[,-1])
stargazer(rho_df[,-1])

g1 =  ggplot(rhos, aes(x=cell_type, y=value, fill=factor(method))) +
  geom_boxplot(outlier.alpha=0.5, outlier.size=1) + theme_bw()

pdf("../figures/rho_boxplot_UCLA_ASD.pdf", width=6, height=3)
g1
dev.off()

g2 = list()
g2[[1]] = ggplot(rho_df, aes(x=Astro.ICeDT, y=Astro.CIBERSORT)) + 
  geom_point(colour = "#009E73", alpha=0.5) + geom_abline(intercept = 0, slope =1) + 
  labs(x = "ICeDT", y="CIBERSORT", title="Astro")

g2[[2]] = ggplot(rho_df, aes(x=Exc.ICeDT, y=Exc.CIBERSORT)) + 
  geom_point(colour = "#009E73", alpha=0.5) + geom_abline(intercept = 0, slope =1) + 
  labs(x = "ICeDT", y="CIBERSORT", title="Exc")

g2[[3]] = ggplot(rho_df, aes(x=Inh.ICeDT, y=Inh.CIBERSORT)) + 
  geom_point(colour = "#009E73", alpha=0.5) + geom_abline(intercept = 0, slope =1) + 
  labs(x = "ICeDT", y="CIBERSORT", title="Inh")

g2[[4]] = ggplot(rho_df, aes(x=Micro.ICeDT, y=Micro.CIBERSORT)) + 
  geom_point(colour = "#009E73", alpha=0.5) + geom_abline(intercept = 0, slope =1) + 
  labs(x = "ICeDT", y="CIBERSORT", title="Micro")

g2[[5]] = ggplot(rho_df, aes(x=Oligo.ICeDT, y=Oligo.CIBERSORT)) + 
  geom_point(colour = "#009E73", alpha=0.5) + geom_abline(intercept = 0, slope =1) + 
  labs(x = "ICeDT", y="CIBERSORT", title="Oligo")

g2[[6]] = ggplot(rho_df, aes(x=OPC.ICeDT, y=OPC.CIBERSORT)) + 
  geom_point(colour = "#009E73", alpha=0.5) + geom_abline(intercept = 0, slope =1) + 
  labs(x = "ICeDT", y="CIBERSORT", title="OPC")

pdf("../figures/rho_scatterplot_UCLA_ASD.pdf", width=7.5, height=5)
ggarrange(g2[[1]], g2[[2]], g2[[3]], g2[[4]], g2[[5]], g2[[6]], ncol = 3, nrow = 2)
dev.off()

# ------------------------------------------------------------------------
# run DESeq2
# ------------------------------------------------------------------------

ctypes2use = which(colnames(icedt_rho) != "Exc")
plog = log(icedt_rho[,ctypes2use] + 0.01)
plog = plog - log(icedt_rho[,which(colnames(icedt_rho) == "Exc")])
colnames(plog) = paste("log", colnames(plog), sep="_")

grp = "bulk"
table(col_data$BrainBank, col_data$Sequencing.Batch)

for(i in 1:ncol(col_data)){
  if(is.character(col_data[[i]])){
    col_data[[i]] = as.factor(col_data[[i]])
  }
}
dim(col_data)
col_data[1:2,]

table(rownames(col_data) == rownames(icedt_rho))
col_data = cbind(col_data, plog)
dim(col_data)
col_data[1:2,c(1:6,(ncol(col_data)-6):ncol(col_data))]

col_data$Diagnosis = factor(col_data$Diagnosis, levels=c("Control", "ASD"))
col_data$DiagnosisP = factor(col_data$DiagnosisP, levels=c("Control", "ASD"))

if (!file.exists(sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4.txt", grp))) {
  ddseqSV4 = DESeqDataSetFromMatrix(countData = trec, 
                               colData = col_data,
                               design = ~ BrainBank + Sequencing.Batch + scaled_AgeDeath + 
                                 scaled_log_depth + scaled_RIN + seqSV1 + seqSV2 + seqSV3 + seqSV4 + Diagnosis)
  ddseqSV4  = DESeq(ddseqSV4)
  #saveRDS(ddseqSV4, file=sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4.rds", grp))
  
  resseqSV4 = results(ddseqSV4)
  dim(resseqSV4)
  resseqSV4[1:2,]
  resseqSV4 = as.data.frame(resseqSV4)
  write.table(resseqSV4, file=sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4.txt", grp), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  resseqSV4_lfcShrink = lfcShrink(ddseqSV4, coef="Diagnosis_ASD_vs_Control")
  resseqSV4_lfcShrink = as.data.frame(resseqSV4_lfcShrink)
  write.table(resseqSV4_lfcShrink, file=sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4_lfcShrink.txt", grp), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  ddseq.CT = DESeqDataSetFromMatrix(countData = trec, 
                                    colData = col_data,
                                    design = ~ BrainBank + Sequencing.Batch + scaled_AgeDeath + 
                                      scaled_log_depth + scaled_RIN + seqSV1 + seqSV2 + seqSV3 + seqSV4 + 
                                      log_Astro + log_Inh + log_Micro + log_Oligo + log_OPC + Diagnosis)
  ddseq.CT  = DESeq(ddseq.CT)
  #saveRDS(ddseqSV4, file=sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4.rds", grp))
  
  ddseq.CT = results(ddseq.CT)
  dim(ddseq.CT)
  ddseq.CT[1:2,]
  ddseq.CT = as.data.frame(ddseq.CT)
  write.table(ddseq.CT, file=sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4_log_ct_fractions.txt", grp), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}


if (!file.exists(sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4_permuted.txt", grp))) {
  ddseqSV4 = DESeqDataSetFromMatrix(countData = trec, 
                               colData = col_data,
                               design = ~ BrainBank + Sequencing.Batch + scaled_AgeDeath + scaled_log_depth + scaled_RIN + seqSV1 + seqSV2 + seqSV3 + seqSV4 + DiagnosisP)
  ddseqSV4  = DESeq(ddseqSV4)
  resseqSV4 = results(ddseqSV4)
  dim(resseqSV4)
  resseqSV4[1:2,]
  resseqSV4 = as.data.frame(resseqSV4)
  
  #saveRDS(ddseqSV4, file=sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4_permuted.rds", grp))
  write.table(resseqSV4, file=sprintf("../results/ASD_step1_DESeq2_%s_adj_covariates_seqSV4_permuted.txt", grp), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}



# TOAST (TPM non-log scale, ICeDT cellular frequencies)

library(TOAST)

cellular_proportions_from_TPM = readRDS("../data/ASD_prop_from_TPM.rds")
rho_from_TPM = cellular_proportions_from_TPM$ICeDT
with_covariates = without_covariates = list()
empty_pval_matrix = matrix(nrow=nrow(rse_filtered), ncol=ncol(rho_from_TPM))
rownames(empty_pval_matrix) = rownames(rse_filtered)
colnames(empty_pval_matrix) = colnames(rho_from_TPM)

##################################################
# permuted
##################################################
with_covariates$pval_matrix = without_covariates$pval_matrix = empty_pval_matrix

# w/ covariates
cell_types = colnames(rho_from_TPM)
design = model.matrix(~ DiagnosisP + BrainBank + Sequencing.Batch + scaled_AgeDeath + scaled_log_depth + scaled_RIN + scaled_sv1 + scaled_sv2 + scaled_sv3 + scaled_sv4,
                      col_data)[, -1]
design_out = makeDesign(as.data.frame(design), rho_from_TPM)
fitted_model = fitModel(design_out, assays(rse_filtered)$TPM)

res_table_list = list()
for (cell_type in cell_types) {
  res_table_list[[cell_type]] = csTest(fitted_model, 
                     coef = "DiagnosisP", 
                     cell_type = cell_type)
  with_covariates$pval_matrix[, cell_type] = res_table_list[[cell_type]][rownames(rse_filtered), "p_value"]
}

opar = par(mfrow=c(2,3))
for (cell_type in cell_types) {
  hist(with_covariates$pval_matrix[, cell_type], breaks=20, main=cell_type, xlab="TOAST with 4 seqSVs (permuted)")
}
par(opar)

write.table(with_covariates$pval_matrix, file=file.path("../results/ASD_step1_TOAST_with_covariates_seqSV4_permuted.txt"), sep="\t", quote=FALSE)


##################################################
# unpermuted
##################################################
with_covariates$pval_matrix = without_covariates$pval_matrix = empty_pval_matrix

# w/ covariates
cell_types = colnames(rho_from_TPM)
design = model.matrix(~ Diagnosis + BrainBank + Sequencing.Batch + scaled_AgeDeath + scaled_log_depth + scaled_RIN + scaled_sv1 + scaled_sv2 + scaled_sv3 + scaled_sv4,
                      col_data)[, -1]
design_out = makeDesign(as.data.frame(design), rho_from_TPM)
fitted_model = fitModel(design_out, assays(rse_filtered)$TPM)

res_table_list = list()
for (cell_type in cell_types) {
  res_table_list[[cell_type]] = csTest(fitted_model, 
                     coef = "Diagnosis", 
                     cell_type = cell_type)
  with_covariates$pval_matrix[, cell_type] = res_table_list[[cell_type]][rownames(rse_filtered), "p_value"]
}

opar = par(mfrow=c(2,3))
for (cell_type in cell_types) {
  hist(with_covariates$pval_matrix[, cell_type], breaks=20, main=cell_type, xlab="TOAST with 4 seqSVs")
}
par(opar)

write.table(with_covariates$pval_matrix, file=file.path("../results/ASD_step1_TOAST_with_covariates_seqSV4.txt"), sep="\t", quote=FALSE)

