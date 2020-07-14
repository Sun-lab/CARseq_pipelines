library(CARseq)
library(SummarizedExperiment)
rse_filtered = readRDS("../data/ASD_rse_filtered_with_SVs.rds")
cellular_proportions = readRDS("../data/ASD_prop.rds")
set.seed(1234)
assay(rse_filtered, "counts") = round(assay(rse_filtered, "counts"))

trec = assays(rse_filtered)$counts
col_data = colData(rse_filtered)

col_data$scaled_AgeDeath = scale(col_data$AgeDeath)
col_data$scaled_log_depth = scale(col_data$log_depth)
col_data$scaled_RIN = scale(col_data$RIN)
for (i in 1:8) {
  col_data[[paste0("scaled_sv", i)]] = scale(col_data[[paste0("sv", i)]])
}

cor(model.matrix(~ BrainBank + Sequencing.Batch + scaled_AgeDeath + scaled_log_depth + scaled_RIN + seqSV1 + seqSV2 + seqSV3 + seqSV4 + Diagnosis, col_data)[,-1])

# permuted disease label
set.seed(1234)
# col_data$DiagnosisP = rep("Control", nrow(col_data))
# wControl = which(col_data$Diagnosis == "Control")
# wCase    = which(col_data$Diagnosis == "ASD")
# 
# ww1 = sample(wControl, size=round(length(wControl)/2))
# ww2 = sample(wCase, size=round(length(wCase)/2))
# 
# col_data$DiagnosisP[c(ww1, ww2)] = "ASD"
col_data$DiagnosisP = CARseq:::permute_case_and_controls(col_data$Diagnosis)
col_data$Diagnosis = factor(col_data$Diagnosis)
col_data$Diagnosis = relevel(col_data$Diagnosis, ref="Control")
col_data$DiagnosisP = relevel(col_data$DiagnosisP, ref="Control")
table(col_data$Diagnosis, col_data$DiagnosisP)

# CARseq, CIBERSORT
# There are too few samples with nonzero estimates for OPC (3 / 85), so OPC is collapsed to Oligo:
table(cellular_proportions$CIBERSORT[, "OPC"] == 0)
cellular_proportions_CIBERSORT_collapsed = cellular_proportions$CIBERSORT
cellular_proportions_CIBERSORT_collapsed[, "Oligo"] = 
  cellular_proportions$CIBERSORT[, "Oligo"] + cellular_proportions$CIBERSORT[, "OPC"]
cellular_proportions_CIBERSORT_collapsed = cellular_proportions_CIBERSORT_collapsed[, -6]  # remove OPC

# permuted
res_CARseq_permuted = run_CARseq(count_matrix = assay(rse_filtered, "counts"),
                 cellular_proportions =  cellular_proportions$CIBERSORT,
                 groups = col_data$DiagnosisP,
                 formula = ~ BrainBank + Sequencing.Batch + scaled_AgeDeath + scaled_log_depth + scaled_RIN + seqSV1 + seqSV2 + seqSV3 + seqSV4,
                 data = col_data,
                 read_depth = 1,
                 shrunken_lfc = TRUE,
                 cores = 16,
                 fix_overdispersion = FALSE
)


pdf("../figures/ASD_CARseq_CIBERSORT_permuted_seqSV4.pdf",
    width=6, height=6)
par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq_permuted$p))) {
  hist(res_CARseq_permuted$p[,k], main=colnames(res_CARseq_permuted$p)[k], breaks=20,
       xlab = paste0("p-value (", sum(CARseq:::get_qvalues_one_inflated(res_CARseq_permuted$p[, k]) < 0.1, na.rm=TRUE), " genes q value < 0.1)"))
}
dev.off()
saveRDS(res_CARseq_permuted, "../results/ASD_CARseq_CIBERSORT_permuted_seqSV4.rds")

# permuted
res_CARseq_permuted = run_CARseq(count_matrix = assay(rse_filtered, "counts"),
                 cellular_proportions =  cellular_proportions$ICeDT,
                 groups = col_data$DiagnosisP,
                 formula = ~ BrainBank + Sequencing.Batch + scaled_AgeDeath + scaled_log_depth + scaled_RIN + seqSV1 + seqSV2 + seqSV3 + seqSV4,
                 data = col_data,
                 read_depth = 1,
                 shrunken_lfc = TRUE,
                 cores = 16,
                 fix_overdispersion = FALSE
)


pdf("../figures/ASD_CARseq_ICeDT_permuted_seqSV4.pdf",
    width=6, height=6)
par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq_permuted$p))) {
  hist(res_CARseq_permuted$p[,k], main=colnames(res_CARseq_permuted$p)[k], breaks=20,
       xlab = paste0("p-value (", sum(CARseq:::get_qvalues_one_inflated(res_CARseq_permuted$p[, k]) < 0.1, na.rm=TRUE), " genes q value < 0.1)"))
}
dev.off()
saveRDS(res_CARseq_permuted, "../results/ASD_CARseq_ICeDT_permuted_seqSV4.rds")
