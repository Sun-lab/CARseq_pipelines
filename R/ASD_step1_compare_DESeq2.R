
# compare the results of DESeq2 with or without cell type proportion 
# (log ratios) as covaraites

library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
theme_set(theme_classic2())
library(CARseq)

res1 = fread("../results/ASD_step1_DESeq2_bulk_adj_covariates_seqSV4.txt")

dim(res1)
res1[1:2,]

res2 = fread("../results/ASD_step1_DESeq2_bulk_adj_covariates_seqSV4_log_ct_fractions.txt")

dim(res2)
res2[1:2,]

table(res1$V1 == res2$V1)

sort(res1$pvalue)[1:5]
sort(res2$pvalue)[1:5]

summary(res1$pvalue)
summary(res2$pvalue)

res1$qvalue = get_qvalues_one_inflated(res1$pvalue)
res2$qvalue = get_qvalues_one_inflated(res2$pvalue)

summary(res1$qvalue)
summary(res2$qvalue)

res.df = merge(res1, res2, by="V1", suffixes = c(".noCT", ".wCT"))
dim(res.df)
res.df[1:2,]

res.df$pvalue.noCT[which(res.df$pvalue.noCT < 1e-8)] = 1e-8
res.df$pvalue.wCT[which(res.df$pvalue.wCT < 1e-8)]   = 1e-8

g1 = ggplot(res.df, aes(x=-log10(pvalue.noCT),y=-log10(pvalue.wCT))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval no CT adj.)", y= "-log10(pval with CT adj.)", 
       title="ASD vs. control") + xlim(0, 8) + ylim(0, 8) + 
  geom_abline(intercept = 0, slope=1)

summary(res.df$log2FoldChange.noCT)
summary(res.df$log2FoldChange.wCT)

g2 = ggplot(res.df, aes(x=log2FoldChange.noCT,
                        y=log2FoldChange.wCT)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "log2 fc, no CT adj.", y= "log2 fc, with CT adj.)", 
       title="ASD vs. control") + 
  geom_abline(intercept = 0, slope=1)

pdf("../figures/ASD_DESeq2_comparison_log10pval.pdf", width=4, height=3.2)
print(g1)
dev.off()

pdf("../figures/ASD_DESeq2_comparison_log2fc.pdf", width=4, height=3.2)
print(g2)
dev.off()

table(res.df$qvalue.noCT < 0.1)
table(res.df$qvalue.wCT  < 0.1)

table(res.df$qvalue.noCT < 0.05)
table(res.df$qvalue.wCT  < 0.05)

table(res.df$qvalue.noCT < 0.05, res.df$qvalue.wCT < 0.05)
table(res.df$qvalue.noCT < 0.1, res.df$qvalue.wCT  < 0.1)

gc()

sessionInfo()
q(save="no")
