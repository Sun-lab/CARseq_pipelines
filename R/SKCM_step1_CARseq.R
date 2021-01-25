
library(CARseq)
library(matrixStats)
library(sva)
library(readxl)
library(parallel)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
source("base_functions_plittle.R")

# samples are from 8 Tissue Source Sites, after filtering.
# we can group them into 5 clusters or just use the TSS with 8 levels
cluster_tss = TRUE

# cell size factor
# here we use the cell size factors in github.com/GfellerLab/EPIC:
# data/mRNA_cell_default.rda:
# Bcells Macrophages   Monocytes Neutrophils     NKcells      Tcells 
# 0.4016      1.4196      1.4196      0.1300      0.4396      0.3952 
# CD4_Tcells  CD8_Tcells  Thelper        Treg  otherCells     default 
# 0.3952      0.3952      0.3952      0.3952      0.4000      0.4000 

cell_sizes = c(0.4016, 0.4, 0.3952, 0.3952, 0.4, 1.4196, 0.4396)
names(cell_sizes) = c("Bcells", "CAFs", "CD4_Tcells", "CD8_Tcells",
                      "Endothelial", "Macrophages", "NKcells")

trec_file_name = "../data/SKCM_trec.rds"
covariate_file_name = "../data/SKCM_cavariates.rds"
cell_fractions_file_name = "../data/SKCM_cell_fraction.rds"

# ------------------------------------------------------------------------
# read in gene expression and cell type proportion estimates
# ------------------------------------------------------------------------

# read counts
counts = readRDS("../data/TCGA_SKCM_raw_counts.rds")
dim(counts)
counts[1:5,1:5]

# remove genes with low read counts: 
# genes with 75% percentile of read counts rd75 < 20 are removed
rd75   = apply(counts, 1, function(x) quantile(x, 0.75))
counts = counts[rd75 >= 20, ]
counts = round(counts)
dim(counts)

# read cell fractions
ct_prop = read.table("../CIBERSORT/CIBERSORTx_SKCM_EPIC_Adjusted.txt",
                                  header = TRUE, row.names = 1)
ct_prop = ct_prop[,1:7]

dim(ct_prop)
ct_prop[1:2,]

stopifnot(all (rownames(ct_prop) %in% names(counts)))
counts = counts[,match(rownames(ct_prop), names(counts))]
table(rownames(ct_prop) == names(counts))

cell_sizes = cell_sizes[names(ct_prop)]
stopifnot(all(names(ct_prop) == names(cell_sizes)))

ct_prop_sized = t(apply(ct_prop, 1, 
                        function(xx){yy = xx / cell_sizes; yy / sum(yy)}))

dim(ct_prop_sized)
ct_prop_sized[1:2,]

# ------------------------------------------------------------------------
# read ABSOLUTE cancer cell fraction estimates
# ------------------------------------------------------------------------

abs_tbl = read.delim("../data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", 
                     as.is = TRUE)
dim(abs_tbl)
abs_tbl[1:2,]

table(substr(rownames(ct_prop_sized), 15, 16))

rownames(abs_tbl) = paste0(gsub("-", ".", abs_tbl$array), "A")
table(rownames(ct_prop_sized) %in% rownames(abs_tbl))

abs_tbl_match = match(rownames(ct_prop_sized), rownames(abs_tbl))
print(table(is.na(abs_tbl_match)))

tumor_purity = abs_tbl$purity[abs_tbl_match]
summary(tumor_purity)
table(tumor_purity < 0.10)
table(tumor_purity < 0.20)
table(tumor_purity < 0.30)

table(tumor_purity > 0.90)
table(tumor_purity > 0.95)
table(tumor_purity > 0.99)

# ------------------------------------------------------------------------
# renormalize cell type proportion based on tumor purity
# ------------------------------------------------------------------------

ct_prop_sized = data.frame(as.matrix(ct_prop_sized) * (1 - tumor_purity), 
                           Cancer_cells = tumor_purity)
ww1 = which(tumor_purity > 0.95)
tumor_purity[ww1]
ct_prop_sized[ww1[1:4],]

ct_prop_long = pivot_longer(ct_prop_sized, everything(), names_to = "cell_type",
                            values_to = "proportion")
dim(ct_prop_long)
ct_prop_long[1:2,]

p = ggplot(ct_prop_long, aes(x=cell_type, y=proportion, color=cell_type)) +
  geom_boxplot() + coord_flip()
ggsave("../figures/SKCM_ct_prop.pdf", p, width=6, height=6)

dim(ct_prop_sized)
ct_prop_sized[1:2,]

# ------------------------------------------------------------------------
# read clinical data
# ------------------------------------------------------------------------

cdat = read_excel("../data/TCGA-CDR-SupplementalTableS1.xlsx", 
                  sheet = "TCGA-CDR", na="#NA", guess_max=2000)
dim(cdat)
cdat[1:2,1:5]
names(cdat)

patient = substr(colnames(counts), 1, 12)
patient = gsub(".", "-", patient, fixed = TRUE)
table(patient %in% cdat$bcr_patient_barcode)

cdat = cdat[match(patient, cdat$bcr_patient_barcode),]
dim(cdat)
cdat[1:2,]

# ------------------------------------------------------------------------
# aggregate patient group into 4 main stages
# ------------------------------------------------------------------------

summary(cdat$DSS.time)
cdat$stage = cdat$ajcc_pathologic_tumor_stage
table(cdat$stage, useNA="ifany")

cdat$stage[grep("Stage I", cdat$ajcc_pathologic_tumor_stage)] = "I"
cdat$stage[grep("Stage II", cdat$ajcc_pathologic_tumor_stage)] = "II"
cdat$stage[grep("Stage III", cdat$ajcc_pathologic_tumor_stage)] = "III"
cdat$stage[grep("Stage IV", cdat$ajcc_pathologic_tumor_stage)] = "IV"
cdat$stage[!(cdat$stage %in% c("I","II","III","IV"))] = NA
table(cdat$stage, useNA="ifany")

# ------------------------------------------------------------------------
# generate an indicator of 5 year survival 
# ------------------------------------------------------------------------

table(cdat$vital_status, cdat$OS)
table(cdat$vital_status, cdat$DSS)
summary(cdat$DSS.time)

cdat$five_year_DSS = rep(NA, nrow(cdat))
cdat$five_year_DSS[which(cdat$DSS.time > 365 * 5)] = 1
cdat$five_year_DSS[which(cdat$DSS.time <= 365 * 5 & cdat$DSS == 1)] = 0
table(cdat$five_year_DSS, useNA="ifany")

# ------------------------------------------------------------------------
# select samples without NA and with purity < 0.90
# ------------------------------------------------------------------------

notNA = !(is.na(cdat$stage)) & !(is.na(cdat$DSS.time)) & !(is.na(tumor_purity))
table(notNA)

notNA = notNA & !(is.na(cdat$five_year_DSS))
table(notNA)

notNA = notNA & (tumor_purity < 0.90)
table(notNA)

trec = as.matrix(counts[notNA])
cell_fractions = ct_prop_sized[notNA, ]
cdat = cdat[notNA,]

rm(counts)
rm(ct_prop)
rm(ct_prop_sized)

# ----------------------------------------------------------------------
# remove lonely samples that from the TSS with few samples
# ----------------------------------------------------------------------

trec[1:2,1:4]
tss = substr(colnames(trec), 6, 7)

tb_tss = table(tss)
tb_tss

lonely_ones = which(tss %in% names(tb_tss)[which(tb_tss <= 4)])
colnames(trec)[lonely_ones]

trec = trec[,-lonely_ones]
cdat = cdat[-lonely_ones,]
cell_fractions = cell_fractions[-lonely_ones,]

patient = substr(colnames(trec), 1, 12)
patient = gsub(".", "-", patient, fixed = TRUE)
table(patient == cdat$bcr_patient_barcode)

saveRDS(trec, file = trec_file_name)
saveRDS(cell_fractions, file = cell_fractions_file_name)

# ----------------------------------------------------------------------
# TSS clustering
# ----------------------------------------------------------------------

tss = substr(colnames(trec), 6, 7)

log_depth = log(apply(trec, 2, function(x) quantile(x, 0.75)))
cdat$log_depth = log_depth

log_trec  = t(log(t(trec)) -log_depth)
log_trec[1:2,1:4]

ge_tss = numeric(0)
u_tss = unique(tss)
for(tss1 in u_tss){
  w1 = which(tss == tss1)
  ge_tss = cbind(ge_tss, rowMeans(log_trec[,w1]))
}
colnames(ge_tss) = u_tss
dim(ge_tss)
ge_tss[1:2,1:3]

ge_tss = ge_tss - rowMeans(ge_tss)
dim(ge_tss)
ge_tss[1:2,1:3]

h1 = hclust(dist(t(ge_tss)))
pdf("../figures/SKCM_tss_cluster.pdf", width=4, height=4)
plot(h1)
dev.off()
c1 = cutree(h1, k=5)
c1

if(cluster_tss){
  cdat$tss = rep(NA, nrow(cdat))
  for(k in 1:5){
    cdat$tss[which(tss %in% names(c1)[c1==k])] = k
  }
}else{
  cdat$tss = tss
}
cdat$tss = as.factor(cdat$tss)
table(cdat$tss, tss)

# ------------------------------------------------------------------------
# generate covariate data
# ------------------------------------------------------------------------

col_data = cdat[, c("age_at_initial_pathologic_diagnosis",
                    "gender", "stage", "five_year_DSS", "tss")]

col_data$scaled_age = scale(col_data$age_at_initial_pathologic_diagnosis)
col_data$scaled_log_depth = scale(cdat$log_depth)

dim(col_data)
col_data[1:2,]

cor(model.matrix(~ gender + scaled_age + scaled_log_depth + stage, col_data)[,-1])

# permuted disease label
set.seed(1234)
col_data = as.data.frame(col_data)
col_data$SurvivalP = CARseq:::permute_case_and_controls(col_data$five_year_DSS)
col_data$Survival  = factor(col_data$five_year_DSS)
table(col_data$Survival, col_data$SurvivalP)

saveRDS(col_data, file = covariate_file_name)

# CARseq using permuted response variable
date()

f1 = ~ gender + scaled_age + scaled_log_depth + stage + tss
f1

res_CARseq = run_CARseq(count_matrix = trec,
                        cellular_proportions = cell_fractions,
                        groups = col_data$SurvivalP,
                        formula = f1, 
                        data = col_data,
                        read_depth = 1,
                        shrunken_lfc = TRUE,
                        cores = 15, 
                        fix_overdispersion = FALSE, 
                        useSocket = FALSE
)
date()

saveRDS(res_CARseq, "../results/SKCM_CARseq_permuted.rds")
pdf("../figures/SKCM_CARseq_permuted.pdf", width=9, height=5)
par(mfrow=c(2,4), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  sk = sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE)
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sk, " genes w/ q < 0.1)"))
}
dev.off()



# CARseq, not permuted
date()

res_CARseq = run_CARseq(count_matrix = trec,
                        cellular_proportions = cell_fractions,
                        groups = col_data$Survival,
                        formula = f1, 
                        data = col_data,
                        read_depth = 1,
                        shrunken_lfc = TRUE,
                        cores = 15, 
                        fix_overdispersion = FALSE,
                        useSocket = FALSE
)
date()

saveRDS(res_CARseq, "../results/SKCM_CARseq.rds")

pdf("../figures/SKCM_CARseq.pdf", width=9, height=5)
par(mfrow=c(2,4), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  sk = sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE)
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sk, " genes w/ q < 0.1)"))
}
dev.off()


gc()

sessionInfo()
q(save="no")
