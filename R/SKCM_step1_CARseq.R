
library(CARseq)
library(matrixStats)
library(sva)
library(readxl)
library(parallel)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
source("base_functions_plittle.R")

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
covariate_file_name = "../data/SKCM_dat_cavariates_with_svs.rds"
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

table(tss)

cdat$tss = rep(NA, nrow(cdat))
for(k in 1:5){
  cdat$tss[which(tss %in% names(c1)[c1==k])] = k
}
table(cdat$tss)
table(cdat$tss, tss)
cdat$tss = as.factor(cdat$tss)
# if (!file.exists(covariate_file_name)) {
# ------------------------------------------------------------------------
# generate covariate data
# ------------------------------------------------------------------------

col_data = cdat[, c("age_at_initial_pathologic_diagnosis",
                    "gender", "stage", "five_year_DSS", "tss")]

col_data$scaled_age = scale(col_data$age_at_initial_pathologic_diagnosis)
col_data$scaled_log_depth = scale(cdat$log_depth)

cor(model.matrix(~ gender + scaled_age + scaled_log_depth + stage, col_data)[,-1])

# ----------------------------------------------------------------------
# first run PCA on residuals to help determin the number of PCs to use
# ----------------------------------------------------------------------

prop = cell_fractions
dat  = col_data

dim(prop)
prop[1:2,]

dim(dat)
dat[1:2,]

plog = log(prop + 1e-5)
plog = plog - plog[,"Cancer_cells"]
colnames(plog) = paste("log", colnames(plog), sep="_")
plog = plog[, !grepl("log_Cancer_cells", colnames(plog))]

dat = cbind(dat, plog)
dim(dat)
dat[1:2,]

mod0.terms = c("gender", "scaled_age", "scaled_log_depth", "stage", 
               "tss", colnames(plog))
length(mod0.terms)

mod0.str = paste(mod0.terms, collapse = " + ")

# Looping over genes for lm() and residuals
GG = nrow(trec); NN = ncol(trec)
rr = matrix(NA,GG,NN) # residual matrix

for(gg in seq(GG)){
  # gg = 1
  if(gg %% 1e2 == 0) cat(".")
  if(gg %% 2e3 == 0 || gg == GG) cat(sprintf("%s out of %s\n",gg,GG))
  
  log.trec.gg = log(trec[gg,] + 1)
  
  lm_out = lm(formula(sprintf("log.trec.gg ~ %s", mod0.str)), data = dat)
  rr[gg,] = as.numeric(lm_out$residuals)
  aa = drop1(lm_out,.~.,test="F")
  
  if(gg == 1){
    pp = matrix(NA,GG,nrow(aa)-1)
    colnames(pp) = rownames(aa)[-1]
  }
  
  pp[gg,] = aa[-1,6]
  rm(lm_out)
}

pdf(file.path("..", "figures", "SKCM_expression_PCA_all_log_Prop.pdf"),
    height = 8,width = 8)

show_pvalue_hist(mat_pvalues = pp, test_type0 = 3)

rr2 = rr - rowMeans(rr,na.rm = TRUE)
cov_rr2 = t(rr2) %*% rr2 / nrow(rr2); pca_rr2 = eigen(cov_rr2)

show_screeplot(pca_rr2,main = "")
par(mfrow = c(2,1),mar = c(4,4,1,1),oma = c(0,0,2,0))
barplot(pca_rr2$values[3:21], names.arg =3:21, main = "", 
        xlab = "Index", ylab = "Eigen-value")
barplot(diff(-pca_rr2$values[3:22]), names.arg =3:21, main = "", 
        xlab = "Index", ylab = "Eigen-value[i] - Eigen-value[i+1]")

num_pcs = 5; pcs = smart_df(pca_rr2$vectors[,1:num_pcs])
names(pcs) = paste0("PC",seq(num_pcs))
show_pc_color(pcs,submain = "")

dev.off()

pcs = pca_rr2$vectors[,1:20]
cor(pcs, dat$five_year_DSS)

# ----------------------------------------------------------------------
# SVA
# ----------------------------------------------------------------------

mod.terms = c("five_year_DSS", mod0.terms)
length(mod.terms)

mod0 = model.matrix(as.formula(paste("~", paste(mod0.terms, collapse=" + "))), 
                    data=dat)
mod  = model.matrix(as.formula(paste("~", paste(mod.terms, collapse=" + "))), 
                    data=dat)

dim(mod0)
mod0[1:2,]

dim(mod)
mod[1:2,]

n.sv = num.sv(log(trec + 1), mod, method="leek")
print(n.sv)

n.sv = num.sv(log(trec + 1), mod, method="be")
print(n.sv)

sv4 = sva(log(trec + 1), mod, mod0, n.sv=4)
sv8 = sva(log(trec + 1), mod, mod0, n.sv=8)


# ----------------------------------------------------------------------
# check associatoin between PCs and svs
# ----------------------------------------------------------------------

pcs = pca_rr2$vectors[,1:8]
npc = 1:8
round(cor(sv4$sv, pcs), 2)
round(cor(sv8$sv, pcs), 2)
round(cor(sv4$sv, sv8$sv), 2)

r2 = matrix(NA, nrow=8, ncol=length(npc))
for(k in 1:8){
  svk = as.numeric(sv8$sv[,k])
  for(m in 1:length(npc)){
    lmk = summary(lm(svk ~ pcs[,1:npc[m]]))
    r2[k,m] = lmk$r.squared
  }
}

round(r2, 2)

pdf(file.path("..","figures", "SKCM_svs_vs_pcs.pdf"),
    height = 6,width = 6)
par(mar=c(5,4,1,1), bty="n")
for(i in 1:nrow(r2)){
  if(i==1){
    plot(npc, r2[1,], ylim=c(0,1), type="l", xlim=c(0,9), 
         xlab="# of PCs", ylab="R2 of surrogate variables explained by PCs")
  }else{
    lines(npc, r2[i,])
  }
}
abline(v=4, lty=2, col="darkred")
abline(v=8, lty=2, col="darkred")
text(rep(9,8), r2[,8], labels=1:8)
dev.off()

# ----------------------------------------------------------------------
# R squared vs. the number of surrogate variables to use
# ----------------------------------------------------------------------

colnames(sv8$sv) = paste0("sv", 1:8)
dat = cbind(dat, sv8$sv)
dim(dat)
dat[1:2,]

log_trec = log(trec + 1)
final.mod0.terms = mod0.terms
final.mod0.str = paste(final.mod0.terms, collapse = " + ")

max_number_of_SV = 8
final.mod0.str.vector = rep(NA, max_number_of_SV + 1)
final.mod0.str.vector[1] = final.mod0.str

for (number_of_SV in (seq_len(max_number_of_SV))) {
  final.mod0.str.vector[1 + number_of_SV] = 
    sprintf("%s + sv%d", final.mod0.str.vector[number_of_SV], number_of_SV)
}

GG = nrow(log_trec); NN = ncol(log_trec)
rr = matrix(NA,GG,NN) # residual matrix

rsquared = matrix(NA, nrow = GG, ncol = 1 + max_number_of_SV)
for(gg in seq_len(GG)){
  # gg = 1
  if(gg %% 1e2 == 0) cat(".")
  if(gg %% 2e3 == 0 || gg == GG) cat(sprintf("%s out of %s\n",gg,GG))
  
  log_trec.gg = log_trec[gg,]
  
  # starts from 0
  for (number_of_SV in (-1 + seq_len(max_number_of_SV + 1))) {
    lm_out = lm(formula(sprintf("log_trec.gg ~ %s", 
                                final.mod0.str.vector[1 + number_of_SV])), 
                data = dat)
    rsquared[gg, 1 + number_of_SV] = summary(lm_out)$r.squared
    rm(lm_out)
  }
}

colnames(rsquared) = paste0("SV",  -1 + seq_len(max_number_of_SV + 1))

saveRDS(rsquared, file.path("..", "results", "SKCM_rsquared.rds"))

colMedians(rsquared, na.rm=TRUE)
colnames(rsquared) = sprintf("%s\n%.2f", 
                             paste0("SV",  -1 + seq_len(max_number_of_SV + 1)), 
                             colMedians(rsquared, na.rm=TRUE))

png(file.path("..", "figures", "SKCM_rsquared.png"),
    width=9, height=5, units="in", res=400)
par(mar = c(6.5, 3.0, 0.5, 1.0))
boxplot(rsquared, sub = paste(strwrap(paste("log_trec ~", final.mod0.str), 110), 
                              collapse = "\n"))
dev.off()

round(cor(sv8$sv, cell_fractions),2)

# ----------------------------------------------------------------------
# save results
# ----------------------------------------------------------------------

mod[1:2,]
colnames(mod)

modDat = cbind(mod[,2:12], sv8$sv)
colnames(modDat)[1:2] = c("five_year_DSS", "gender") 
colnames(modDat)[12:19] = paste0("sv", 1:8)
dim(modDat)
modDat[1:2,]

# rescale PCs and SVs -- the names are still sv1, sv2, etc, 
# but they have been scaled:
covariates_to_scale = grep("sv|PC", colnames(modDat))
for (m in covariates_to_scale) {
  modDat[, m] = scale(modDat[, m])
}

col_data = modDat
saveRDS(col_data, file=covariate_file_name)


# permuted disease label
set.seed(1234)
col_data = as.data.frame(col_data)
col_data$SurvivalP = CARseq:::permute_case_and_controls(col_data$five_year_DSS)
col_data$Survival  = factor(col_data$five_year_DSS)
table(col_data$Survival, col_data$SurvivalP)


# CARseq, CIBERSORT, 0 SV
# permuted
date()
res_CARseq = run_CARseq(count_matrix = trec,
                        cellular_proportions = cell_fractions,
                        groups = col_data$SurvivalP,
                        formula = ~ gender + scaled_age + scaled_log_depth + 
                          stageII + stageIII + stageIV + tss,
                        data = col_data,
                        read_depth = 1,
                        shrunken_lfc = TRUE,
                        cores = 15, 
                        fix_overdispersion = FALSE, 
                        useSocket = FALSE
)
date()

saveRDS(res_CARseq, "../results/SKCM_CARseq_CIBERSORT_permuted_SV0.rds")
pdf("../figures/SKCM_CARseq_CIBERSORT_permuted_SV0.pdf",
    width=6, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  sk = sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE)
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sk, " genes w/ q < 0.1)"))
}
dev.off()



# CARseq, CIBERSORT, 0 SV
# not permuted
date()

res_CARseq = run_CARseq(count_matrix = trec,
                        cellular_proportions = cell_fractions,
                        groups = col_data$Survival,
                        formula = ~ gender + scaled_age + scaled_log_depth + 
                          stageII + stageIII + stageIV + tss,
                        data = col_data,
                        read_depth = 1,
                        shrunken_lfc = TRUE,
                        cores = 15, 
                        fix_overdispersion = FALSE,
                        useSocket = FALSE
)
date()

saveRDS(res_CARseq, "../results/SKCM_CARseq_CIBERSORT_SV0.rds")

pdf("../figures/SKCM_CARseq_CIBERSORT_SV0.pdf", width=6, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  sk = sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE)
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sk, " genes w/ q < 0.1)"))
}
dev.off()


# CARseq, CIBERSORT, 4 SV
# permuted
date()

res_CARseq = run_CARseq(count_matrix = trec,
                        cellular_proportions = cell_fractions,
                        groups = col_data$SurvivalP,
                        formula = ~ gender + scaled_age + scaled_log_depth + 
                          stageII + stageIII + stageIV + tss + 
                          sv1 + sv2 + sv3 + sv4,
                        data = col_data,
                        read_depth = 1,
                        shrunken_lfc = TRUE,
                        cores = 15, 
                        fix_overdispersion = FALSE,
                        useSocket = FALSE
)
date()

saveRDS(res_CARseq, "../results/SKCM_CARseq_CIBERSORT_permuted_SV4.rds")

pdf("../figures/SKCM_CARseq_CIBERSORT_permuted_SV4.pdf", width=6, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  sk = sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE)
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sk, " genes w/ q < 0.1)"))
}
dev.off()


# CARseq, CIBERSORT, 4 SV
# not permuted
date()

res_CARseq = run_CARseq(count_matrix = trec,
                        cellular_proportions = cell_fractions,
                        groups = col_data$Survival,
                        formula = ~ gender + scaled_age + scaled_log_depth + 
                          stageII + stageIII + stageIV + tss + 
                          sv1 + sv2 + sv3 + sv4,
                        data = col_data,
                        read_depth = 1,
                        shrunken_lfc = TRUE,
                        cores = 15, 
                        fix_overdispersion = FALSE,
                        useSocket = FALSE
)
date()

saveRDS(res_CARseq, "../results/SKCM_CARseq_CIBERSORT_SV4.rds")

pdf("../figures/SKCM_CARseq_CIBERSORT_SV4.pdf", width=6, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  sk = sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE)
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sk, " genes w/ q < 0.1)"))
}
dev.off()




# CARseq, CIBERSORT, 1 SV
# permuted
date()

res_CARseq = run_CARseq(count_matrix = trec,
                        cellular_proportions = cell_fractions,
                        groups = col_data$SurvivalP, 
                        formula = ~ gender + scaled_age + scaled_log_depth + 
                          stageII + stageIII + stageIV + tss + sv1, 
                        data = col_data, 
                        read_depth = 1, 
                        shrunken_lfc = TRUE, 
                        cores = 15, 
                        fix_overdispersion = FALSE,
                        useSocket = FALSE
)
date()

saveRDS(res_CARseq, "../results/SKCM_CARseq_CIBERSORT_permuted_sv1.rds")

pdf("../figures/SKCM_CARseq_CIBERSORT_permuted_sv1.pdf", width=6, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  sk = sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE)
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sk, " genes w/ q<0.1)"))
}
dev.off()


# CARseq, CIBERSORT, 1 SV
# not permuted
date()

res_CARseq = run_CARseq(count_matrix = trec,
                        cellular_proportions = cell_fractions,
                        groups = col_data$Survival,
                        formula = ~ gender + scaled_age + scaled_log_depth + 
                          stageII + stageIII + stageIV + tss + sv1,
                        data = col_data,
                        read_depth = 1,
                        shrunken_lfc = TRUE,
                        cores = 15, 
                        fix_overdispersion = FALSE,
                        useSocket = FALSE
)
date()

saveRDS(res_CARseq, "../results/SKCM_CARseq_CIBERSORT_sv1.rds")

pdf("../figures/SKCM_CARseq_CIBERSORT_sv1.pdf", width=6, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(res_CARseq$p))) {
  sk = sum(CARseq:::get_qvalues_one_inflated(res_CARseq$p[, k]) < 0.1, na.rm=TRUE)
  hist(res_CARseq$p[,k], main=colnames(res_CARseq$p)[k], breaks=20,
       xlab = paste0("p-value (", sk, " genes w/ q<0.1)"))
}
dev.off()

gc()

sessionInfo()
q(save="no")
