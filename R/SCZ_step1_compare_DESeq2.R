
# compare the results of DESeq2 with or without cell type proportion 
# (log ratios) as covaraites

library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(CARseq)
library(DESeq2)
theme_set(theme_classic())
library(SummarizedExperiment)
library(foreach)
library(doParallel)

# -------------------------------------------------------------------------
# compare DESeq2 analysis results
# -------------------------------------------------------------------------

res1 = fread("../results/SCZ_step1_DESeq2_bulk_adj_covariates_SV2.txt")

dim(res1)
res1[1:2,]

fnm2 = "SCZ_step1_DESeq2_bulk_adj_covariates_SV2_log_ct_fractions.txt"
res2 = fread(file.path("../results", fnm2))

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

pdf("../figures/SCZ_DESeq2_compare_log10pval.pdf", width=4, height=3.2)
print(g1)
dev.off()

pdf("../figures/SCZ_DESeq2_compare_log2fc.pdf", width=4, height=3.2)
print(g2)
dev.off()

# -------------------------------------------------------------------------
# read in gene expressio and clinical data
# -------------------------------------------------------------------------

private_data = "../../CARseq_pipelines_private/data"
trec = readRDS(file.path(private_data, "trec_filtered_scz_control.rds"))

fnm = "dat_cavariates_scz_control_with_svs.rds"
clinical_variables_raw = readRDS(file.path(private_data, fnm))

dim(clinical_variables_raw)
clinical_variables_raw[1:2,]

# for code compatibility:
col_data    = as.data.frame(clinical_variables_raw)
col_data$Dx = ifelse(col_data$DxSCZ == 1, "SCZ", "Control")
col_data$RNAseq.Sample_RNA_ID = colnames(trec)

# rescale age_death, PMI, RIN, RIN^2, log_depth
col_data$scaled_age_death = scale(col_data$age_death)
col_data$scaled_log_depth = scale(col_data$log_depth)
col_data$scaled_PMI  = scale(col_data$PMI)
col_data$scaled_RIN  = scale(col_data$RIN)
col_data$scaled_RIN2 = scale(col_data$RIN2)

# rescale PCs and SVs -- the names are still sv1, sv2, etc, 
# but they have been scaled:
covariates_to_scale = grep("sv|PC", colnames(col_data))
for (m in covariates_to_scale) {
  col_data[, m] = scale(col_data[, m])
}

# cellular proportions
prop_output_file = "../data/SCZ_prop.rds"
prop_list = readRDS(prop_output_file)
names(prop_list)

icedt_rho    = prop_list$ICeDT
match_sample = match(col_data$RNAseq.Sample_RNA_ID, rownames(icedt_rho))
icedt_rho    = icedt_rho[match_sample,]

dim(icedt_rho)
icedt_rho[1:2,]

rownames(col_data) = col_data$RNAseq.Sample_RNA_ID
table(rownames(col_data) == rownames(icedt_rho))

for(i in 1:ncol(col_data)){
  if(is.character(col_data[[i]])){
    col_data[[i]] = as.factor(col_data[[i]])
  }
}

ctypes2use = which(colnames(icedt_rho) != "Exc")
plog = log(icedt_rho[,ctypes2use] + 0.01)
plog = plog - log(icedt_rho[,which(colnames(icedt_rho) == "Exc")])
colnames(plog) = paste("log", colnames(plog), sep="_")

col_data = cbind(col_data, plog)
dim(col_data)
col_data[1:2,c(1:6,(ncol(col_data)-6):ncol(col_data))]

col_data$Diagnosis = col_data$Dx

# ------------------------------------------------------------------------
# run linear regression, estimate the R2 by all the covariates 
# with or without cell type proportions
# ------------------------------------------------------------------------

rd = colSums(trec)
summary(rd)
rd = rd/median(rd)
summary(rd)

log_trec = log(t(t(trec + 1)/rd))
dim(log_trec)
log_trec[1:2,1:5]

R2s = matrix(NA, nrow=nrow(trec), ncol=2)

nCore = 15

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

R2s = foreach(i = 1:nrow(trec), .combine = "rbind") %dopar% {
  y   = log_trec[i,]
  
  lm0 = lm(y ~ scaled_log_depth + InstitutionPenn + InstitutionPitt + 
             genderMale + scaled_age_death + scaled_PMI + scaled_RIN + 
             scaled_RIN2 + genoPC1 + genoPC2 + libclustB + 
             libclustbase + libclustC + libclustD + libclustE + 
             libclustF + libclustG + sv1 + sv2 + Diagnosis, 
           data=col_data)

  lm1 = lm(y ~ scaled_log_depth + InstitutionPenn + InstitutionPitt + 
             genderMale + scaled_age_death + scaled_PMI + scaled_RIN + 
             scaled_RIN2 + genoPC1 + genoPC2 + libclustB + 
             libclustbase + libclustC + libclustD + libclustE + 
             libclustF + libclustG + sv1 + sv2 + Diagnosis + 
             log_Astro + log_Inh + log_Micro + log_Oligo + log_OPC, 
           data=col_data)
  
  a1 = anova(lm0, lm1)
  s0 = summary(lm0)
  s1 = summary(lm1)
  
  c(s0$r.squared, s1$r.squared, a1$`Pr(>F)`[2])
}

dim(R2s)
R2s[1:2,]

R2s = as.data.frame(R2s)
rownames(R2s) = rownames(trec)
names(R2s) = c("R2_no_ct", "R2_with_ct", "pval")
table(rownames(R2s) == res.df$V1)

g2 = ggplot(R2s, aes(x=R2_no_ct, y=R2_with_ct)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "R2, no CT adj.", y= "R2, with CT adj.") + 
  geom_abline(intercept = 0, slope=1)

pdf("../figures/SCZ_DESeq2_compare_R2_with_or_without_log_ct_prop.pdf", 
    width=4, height=3.2)
print(g2)
dev.off()

# ------------------------------------------------------------------------
# check the association between ct-R2 and wether the gene expression
# is assocaited with case/control status with/without ct correction
# ------------------------------------------------------------------------

table(res.df$qvalue.noCT < 0.1)
table(res.df$qvalue.wCT  < 0.1)

table(res.df$qvalue.noCT < 0.2)
table(res.df$qvalue.wCT  < 0.2)

summary(R2s)
table(R2s$pval < 0.05/nrow(R2s))

table(res.df$qvalue.noCT < 0.1, res.df$qvalue.wCT < 0.1)
table(res.df$qvalue.noCT < 0.2, res.df$qvalue.wCT < 0.2)

stats = list()
qcuts = c(0.1, 0.2)

for(k in 1:length(qcuts)){
  qcut = qcuts[k]
  w00 = res.df$qvalue.noCT >= qcut & res.df$qvalue.wCT >= qcut
  w10 = res.df$qvalue.noCT < qcut  & res.df$qvalue.wCT >= qcut
  w01 = res.df$qvalue.noCT >= qcut & res.df$qvalue.wCT < qcut
  w11 = res.df$qvalue.noCT < qcut  & res.df$qvalue.wCT < qcut
  
  f00 = fisher.test(R2s$pval < 0.05/nrow(R2s), w00)
  f10 = fisher.test(R2s$pval < 0.05/nrow(R2s), w10)
  f01 = fisher.test(R2s$pval < 0.05/nrow(R2s), w01)
  f11 = fisher.test(R2s$pval < 0.05/nrow(R2s), w11)
  
  stat1 = c(f00$estimate, f00$conf.int, f00$p.value)
  stat1 = rbind(stat1, c(f10$estimate, f10$conf.int, f10$p.value))
  stat1 = rbind(stat1, c(f01$estimate, f01$conf.int, f01$p.value))
  stat1 = rbind(stat1, c(f11$estimate, f11$conf.int, f11$p.value))
  
  stats[[k]] = stat1 
}

stats

table(R2s$pval < 0.05/nrow(R2s))

qcut = 0.1

w00 = res.df$qvalue.noCT >= qcut & res.df$qvalue.wCT >= qcut
c1 = chisq.test(R2s$pval < 0.05/nrow(R2s), w00)
c1$expected/rowSums(c1$observed)
c1$observed/rowSums(c1$observed)

w10 = res.df$qvalue.noCT < qcut & res.df$qvalue.wCT >= qcut
c1 = chisq.test(R2s$pval < 0.05/nrow(R2s), w10)
c1$expected/rowSums(c1$observed)
c1$observed/rowSums(c1$observed)

w01 = res.df$qvalue.noCT >= qcut & res.df$qvalue.wCT < qcut
c1 = chisq.test(R2s$pval < 0.05/nrow(R2s), w01)
c1$expected/rowSums(c1$observed)
c1$observed/rowSums(c1$observed)

w11 = res.df$qvalue.noCT < qcut & res.df$qvalue.wCT < qcut
c1 = chisq.test(R2s$pval < 0.05/nrow(R2s), w11)
c1$expected/rowSums(c1$observed)
c1$observed/rowSums(c1$observed)

# ------------------------------------------------------------------------
# plot the odds ratio
# ------------------------------------------------------------------------

grp = c("q < 0.1 with or w/o accounting for\n cell type (ct) proportion (prop.)", 
        "q < 0.1 with ct prop. and \nq >= 0.1 w/o ct prop.", 
        "q >= 0.1 with ct prop. and \nq < 0.1 w/o ct prop. ", 
        "q < 0.1 with or w/o ct prop.")
yAxis = length(grp):1

df_or = data.frame(group = grp, yAxis = yAxis, OR = stats[[1]][,1], 
                   CI_low = stats[[1]][,2], CI_high = stats[[1]][,3], 
                   pval = paste0("p = ", signif(stats[[1]][,4],2)))

p1 = ggplot(df_or, aes(x = OR, y = yAxis, label= pval)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = CI_high, xmin = CI_low), size = .5, 
                 height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = yAxis, labels = grp) + ylab("") + 
  xlab("Odds Ratio") + theme_bw() + 
  ggtitle("SCZ association vs. cell type prop. association") +
  geom_text(nudge_x = 0.15, nudge_y = - 0.2)

pdf("../figures/SCZ_DESeq2_qval_vs_log_ct_association.pdf", 
    width=6, height=2.5)
print(p1)
dev.off()

# ------------------------------------------------------------------------
# Compare the significance level for different groups of genes
# ------------------------------------------------------------------------

R2s$R2_diff = R2s$R2_with_ct - R2s$R2_no_ct
R2s$grp = R2s$group = rep(NA, nrow(R2s))
R2s$grp[w01] = "01"
R2s$grp[w10] = "10"
R2s$grp[w11] = "11"

R2s$group[(w00 | w01)] = "q >= 0.1"
R2s$group[(w10 | w11)] = "q < 0.1"

res.df$grp   = R2s$grp
res.df$group = R2s$group

ttl = "while accounting for cell type prop."
g1 = ggplot(subset(res.df, grp %in% c("01", "11")), 
            aes(x=-log10(qvalue.wCT), color=group)) + 
  geom_density() + theme(legend.position=c(0.55,0.8)) + xlim(1,5) + 
  xlab("-log10(qvalue for SCZ association \n given cell type prop.)") + 
  ggtitle(paste0("SCZ-associated genes (q < 0.1) \n ", ttl)) + 
  labs(color = "q-value without accounting \n for cell type prop.")

pdf("../figures/SCZ_DESeq2_qval_density.pdf", 
    width=4, height=3.5)
print(g1)
dev.off()

# -------------------------------------------------------------------------
# Re-run DESeq2 given log transformed cell type composition of Exc
# -------------------------------------------------------------------------

dim(trec)
trec[1:2,1:3]

col_data = cbind(col_data, log(icedt_rho + 0.01))

dim(col_data)
col_data[1:2,]

dd3 = DESeqDataSetFromMatrix(countData = trec, 
                             colData = col_data,
                             design = ~ scaled_log_depth + InstitutionPenn + 
                               InstitutionPitt + genderMale + scaled_age_death + 
                               scaled_PMI + scaled_RIN + scaled_RIN2 +
                               genoPC1 + genoPC2 + libclustB + libclustbase + 
                               libclustC + libclustD + libclustE + libclustF + 
                               libclustG + Exc + sv1 + sv2 + Diagnosis)
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
                             design = ~ scaled_log_depth + InstitutionPenn + 
                               InstitutionPitt + genderMale + scaled_age_death + 
                               scaled_PMI + scaled_RIN + scaled_RIN2 +
                               genoPC1 + genoPC2 + libclustB + libclustbase + 
                               libclustC + libclustD + libclustE + libclustF + 
                               libclustG + Inh + sv1 + sv2 + Diagnosis)
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

fnm = "../results/SCZ_step1_DESeq2_bulk_adj_covariates_SV2_log_Exc.txt"
write.table(res3, file=fnm, quote = FALSE, sep = "\t", 
            row.names = TRUE, col.names = TRUE)

fnm = "../results/SCZ_step1_DESeq2_bulk_adj_covariates_SV2_log_Inh.txt"
write.table(res4, file=fnm, quote = FALSE, sep = "\t", 
            row.names = TRUE, col.names = TRUE)

gc()

sessionInfo()
q(save="no")


