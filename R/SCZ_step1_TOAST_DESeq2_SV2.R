
library(SummarizedExperiment)
library(CARseq)
library(ggplot2)
library(ggpubr)
library(tidyr)
theme_set(theme_bw())
library(tidyverse)
library(stargazer)

# total read count
trec = readRDS("../data/trec_filtered_scz_control.rds")  # 20788x527 matrix
clinical_variables_raw = readRDS("../data/dat_cavariates_scz_control_with_svs.rds")

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
icedt_rho = prop_list$ICeDT
cibersort_rho = prop_list$CIBERSORT

# plot CIBERSORT vs. ICeDT cellular frequency estimates

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

pdf("../figures/rho_boxplot_CMC_SCZ.pdf", width=6, height=3)
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

pdf("../figures/rho_scatterplot_CMC_SCZ.pdf", width=7.5, height=5)
ggarrange(g2[[1]], g2[[2]], g2[[3]], g2[[4]], g2[[5]], g2[[6]], ncol = 3, nrow = 2)
dev.off()

# TOAST (TPM non-log scale, ICeDT cellular frequencies)

library(TOAST)

# the cellular frequencies directly obtained from the deconvolution of TPM, without adjustment of cell sizes
cellular_proportions_from_TPM = readRDS("../data/SCZ_prop_from_TPM.rds")  
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
design = model.matrix(~ DiagnosisP + scaled_log_depth + InstitutionPenn + InstitutionPitt +
                        genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 +
                        genoPC1 + genoPC2 +
                        libclustB + libclustbase + libclustC + libclustD + libclustE + libclustF + libclustG +
                        sv1 + sv2, data = col_data)[, -1]
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
  hist(with_covariates$pval_matrix[, cell_type], breaks=20, main=cell_type, xlab="TOAST with 2 SVs (permuted)")
}
par(opar)

write.table(with_covariates$pval_matrix, file=file.path("../results/step1_TOAST_with_covariates_SV2_permuted.txt"), sep="\t", quote=FALSE)


##################################################
# unpermuted
##################################################
with_covariates$pval_matrix = without_covariates$pval_matrix = empty_pval_matrix

# w/ covariates
cell_types = colnames(rho_from_TPM)
design = model.matrix(~ Diagnosis + scaled_log_depth + InstitutionPenn + InstitutionPitt +
                        genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 +
                        genoPC1 + genoPC2 +
                        libclustB + libclustbase + libclustC + libclustD + libclustE + libclustF + libclustG +
                        sv1 + sv2, data = col_data)[, -1]
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
  hist(with_covariates$pval_matrix[, cell_type], breaks=20, main=cell_type, xlab="TOAST with 2 SVs")
}
par(opar)

write.table(with_covariates$pval_matrix, file=file.path("../results/step1_TOAST_with_covariates_SV2.txt"), sep="\t", quote=FALSE)


library(DESeq2)
grp = "bulk"
# table(col_data$BrainBank, col_data$Sequencing.Batch)
# ------------------------------------------------------------------------
# run DESeq2
# ------------------------------------------------------------------------

ctypes2use = which(colnames(icedt_rho) != "Exc")
plog = log(icedt_rho[,ctypes2use] + 0.01)
plog = plog - log(icedt_rho[,which(colnames(icedt_rho) == "Exc")])
colnames(plog) = paste("log", colnames(plog), sep="_")

dim(plog)
plog[1:2,]
summary(plog)
apply(plog, 2, sd)

for(i in 1:ncol(col_data)){
  if(is.character(col_data[[i]])){
    col_data[[i]] = as.factor(col_data[[i]])
  }
}
dim(col_data)
col_data[1:2,]
summary(col_data)
rownames(col_data) = col_data$RNAseq.Sample_RNA_ID

table(rownames(col_data) == rownames(icedt_rho))
col_data = cbind(col_data, plog)
dim(col_data)
col_data[1:2,]

if (!file.exists(sprintf("../results/SCZ_step1_DESeq2_%s_adj_covariates_SV2.txt", grp))) {

  dd2 = DESeqDataSetFromMatrix(countData = trec, 
                               colData = col_data,
                               design = ~ scaled_log_depth + InstitutionPenn + InstitutionPitt + 
                                 genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 + 
                                 genoPC1 + genoPC2 + libclustB + libclustbase + libclustC + libclustD + 
                                 libclustE + libclustF + libclustG + sv1 + sv2 + Diagnosis)
  dd2 = DESeq(dd2)
  #saveRDS(dd2, file=sprintf("results/SCZ_step1_DESeq2_%s_adj_covariates_SV2.rds", grp))
  
  res2 = results(dd2)
  dim(res2)
  res2[1:2,]
  res2 = as.data.frame(res2)
  write.table(res2, file=sprintf("../results/SCZ_step1_DESeq2_%s_adj_covariates_SV2.txt", grp), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  res2_lfcShrink = lfcShrink(dd2, coef="Diagnosis_SCZ_vs_Control")
  res2_lfcShrink = as.data.frame(res2_lfcShrink)
  write.table(res2_lfcShrink, file=sprintf("../results/SCZ_step1_DESeq2_%s_adj_covariates_SV2_lfcShrink.txt", grp), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  dd3 = DESeqDataSetFromMatrix(countData = trec, 
                               colData = col_data,
                               design = ~ scaled_log_depth + InstitutionPenn + InstitutionPitt +
                                 genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 +
                                 genoPC1 + genoPC2 + libclustB + libclustbase + libclustC + libclustD + 
                                 libclustE + libclustF + libclustG + log_Astro + log_Inh + log_Micro + 
                                 log_Oligo + log_OPC + sv1 + sv2 + Diagnosis)
  dd3 = DESeq(dd3)
  res3 = results(dd3)
  dim(res3)
  res3[1:2,]
  res3 = as.data.frame(res3)
  write.table(res3, file=sprintf("../results/SCZ_step1_DESeq2_%s_adj_covariates_SV2_log_ct_fractions.txt", grp), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
}

if (!file.exists(sprintf("../results/SCZ_step1_DESeq2_%s_adj_covariates_SV2_permuted.txt", grp))) {

  dd2 = DESeqDataSetFromMatrix(countData = trec, 
                               colData = col_data,
                               design = ~ scaled_log_depth + InstitutionPenn + InstitutionPitt +
                          genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 +
                          genoPC1 + genoPC2 +
                          libclustB + libclustbase + libclustC + libclustD + libclustE + libclustF + libclustG +
                          sv1 + sv2 + DiagnosisP)
  dd2 = DESeq(dd2)
  
  res2 = results(dd2)
  dim(res2)
  res2[1:2,]
  
  #saveRDS(dd2, file=sprintf("../results/SCZ_step1_DESeq2_%s_adj_covariates_SV2_permuted.rds", grp))
  
  res2 = as.data.frame(res2)
  
  write.table(res2, file=sprintf("../results/SCZ_step1_DESeq2_%s_adj_covariates_SV2_permuted.txt", grp), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}




