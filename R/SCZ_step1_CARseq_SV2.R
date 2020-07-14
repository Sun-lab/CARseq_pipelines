library(CARseq)
library(SummarizedExperiment)

# read files from "CARseq_pipelines" repository
trec = readRDS("../data/trec_filtered_scz_control.rds")  # 20788x527 matrix
clinical_variables_raw = readRDS("../data/dat_cavariates_scz_control_with_svs.rds")    # 527x33 matrix; there are 527 instead of 537 samples
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

cellular_proportions = readRDS("../data/SCZ_prop.rds")

cellular_proportions$ICeDT = cellular_proportions$ICeDT[
    match(col_data$RNAseq.Sample_RNA_ID, rownames(cellular_proportions$ICeDT)), ]
cellular_proportions$CIBERSORT = cellular_proportions$CIBERSORT[
    match(col_data$RNAseq.Sample_RNA_ID, rownames(cellular_proportions$CIBERSORT)), ]





res_CARseq = run_CARseq(count_matrix = trec,
                                 cellular_proportions = cellular_proportions$ICeDT,
                                 groups = col_data$Diagnosis,
                                 formula = ~ scaled_log_depth + InstitutionPenn + InstitutionPitt +
                                   genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 +
                                   genoPC1 + genoPC2 +
                                   libclustB + libclustbase + libclustC + libclustD + libclustE + libclustF + libclustG +
                                   sv1 + sv2,
                                 data = col_data,
                                 read_depth = 1,
                                 shrunken_lfc = TRUE,
                                 cores = 32,
                                 fix_overdispersion = FALSE
)

pdf("../figures/SCZ_CARseq_ICeDT_SV2.pdf",
    width=6, height=6)
par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE), " genes q value < 0.1)"))
}
dev.off()
saveRDS(res_CARseq, "../results/SCZ_CARseq_ICeDT_SV2.rds")



res_CARseq = run_CARseq(count_matrix = trec,
                                 cellular_proportions = cellular_proportions$CIBERSORT,
                                 groups = col_data$Diagnosis,
                                 formula = ~ scaled_log_depth + InstitutionPenn + InstitutionPitt +
                                   genderMale + scaled_age_death + scaled_PMI + scaled_RIN + scaled_RIN2 +
                                   genoPC1 + genoPC2 +
                                   libclustB + libclustbase + libclustC + libclustD + libclustE + libclustF + libclustG +
                                   sv1 + sv2,
                                 data = col_data,
                                 read_depth = 1,
                                 shrunken_lfc = TRUE,
                                 cores = 32,
                                 fix_overdispersion = FALSE
)

pdf("../figures/SCZ_CARseq_CIBERSORT_SV2.pdf",
    width=6, height=6)
par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE), " genes q value < 0.1)"))
}
dev.off()
saveRDS(res_CARseq, "../results/SCZ_CARseq_CIBERSORT_SV2.rds")

