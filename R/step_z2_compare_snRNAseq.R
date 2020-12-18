
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
theme_set(theme_classic2())
library(ggcorrplot)
library(qvalue)
library(ggpointdensity)

# ----------------------------------------------------------------------
# read in p-values for ASD
# ----------------------------------------------------------------------

pval.path = "../results/_pvalues"

pval_files = list.files(path=pval.path, pattern="ASD", full.names=TRUE)
pval_files

f1 = pval_files[1]
pval_carseq = fread(f1)
lb = str_extract(basename(f1), "(?<=ASD_)\\w+_\\w+(?=.txt)")
names(pval_carseq)[3] = lb

dim(pval_carseq)
pval_carseq[1:2,]

length(unique(pval_carseq$gene_id))
length(unique(pval_carseq$gene_name))

for(f1 in pval_files[-1]){
  d1 = fread(f1)
  if(any(d1$gene_id != pval_carseq$gene_id)){
    stop("gene_id do not match\n")
  }
  
  bn = basename(f1)
  lb = str_extract(bn, "(?<=ASD_)\\w+_\\w+(?=.txt)")
  
  pval_carseq[[lb]] = d1$pvalue
}

dim(pval_carseq)
pval_carseq[1:2,]

# ----------------------------------------------------------------------
# read in p-values from snRNA-seq analysis
# ----------------------------------------------------------------------

cts = c("Astro", "Inh", "Exc", "Micro", "Oligo", "OPC")
cts

grp = cts[1]
fnm = sprintf("../results/stepz1_DESeq2_%s.txt", grp)

deseq2 = read.table(file=fnm, sep="\t", header=TRUE)
dim(deseq2)
deseq2[1:2,]

pval_snseq = data.frame(gene_name = rownames(deseq2), 
                        stringsAsFactors=FALSE)
percent_zeros = pval_snseq

pval_snseq[[grp]] = deseq2$pvalue
percent_zeros[[grp]] = deseq2$percent_zeros

for(grp in cts[-1]){
  fnm = sprintf("../results/stepz1_DESeq2_%s.txt", grp)
  
  deseq2 = read.table(file=fnm, sep="\t", header=TRUE)
  
  if(nrow(deseq2) != nrow(pval_snseq)){
    stop("the number of rows of deseq2 results do not match\n")
  }
  
  if(any(rownames(deseq2) != pval_snseq$gene_name)){
    stop("rownames of deseq2 results do not match\n")
  }
  
  pval_snseq[[grp]] = deseq2$pvalue
  percent_zeros[[grp]] = deseq2$percent_zeros
}

dim(pval_snseq)
pval_snseq[1:2,]

dim(percent_zeros)
percent_zeros[1:2,]

colSums(pval_snseq[,-1] < 0.001)
colSums(percent_zeros[,-1] > 0.5)

# ----------------------------------------------------------------------
# take the intersection of genes from the two studies
# ----------------------------------------------------------------------

nrow(pval_snseq)
nrow(pval_carseq)

pval_snseq$gene_name[1:2]
pval_carseq$gene_name[1:2]

gene2use = intersect(pval_carseq$gene_name, pval_snseq$gene_name)
length(gene2use)

pval_snseq  = pval_snseq[match(gene2use,  pval_snseq$gene_name),]
pval_carseq = pval_carseq[match(gene2use, pval_carseq$gene_name),]

dim(pval_snseq)
pval_snseq[1:2,]

pval_carseq = as.data.frame(pval_carseq)
dim(pval_carseq)
pval_carseq[1:2,]

table(pval_snseq$gene_name == pval_carseq$gene_name)

saveRDS(pval_carseq, "../results/step_z2_pval_carseq.rds")
saveRDS(pval_snseq,  "../results/step_z2_pval_snseq.rds")

# ----------------------------------------------------------------------
# check the correlation matrix within each method
# ----------------------------------------------------------------------

cor_carseq = cor(pval_carseq[,-(1:2)], method = "spearman", use="pair")
summary(c(cor_carseq))

cor_snseq = cor(pval_snseq[,-(1:1)], method = "spearman", use="pair")
summary(c(cor_snseq))

gc1 = ggcorrplot(cor_carseq, tl.cex = 6)  

pdf("../figures/step_z2_cor_CARseq.pdf", width=3.6, height=2.8)
print(gc1)
dev.off()

gc2 = ggcorrplot(cor_snseq, tl.cex = 6)  

pdf("../figures/step_z2_cor_snseq.pdf", width=2.7, height=2.1)
print(gc2)
dev.off()

# ----------------------------------------------------------------------
# check the correlation matrix between methods
# ----------------------------------------------------------------------

cormat = cor(pval_carseq[,-(1:2)], pval_snseq[,-1], 
             method = "spearman", use="pair")
summary(c(cormat))

gc1 = ggcorrplot(t(cormat), tl.cex = 6, title="Spearman correlation") + 
  scale_fill_gradient2(limit = c(-0.08,0.08), low = "blue", high =  "red", 
                       mid = "white", midpoint = 0) 

pdf("../figures/step_z2_cor_CARseq_vs_snRNAseq.pdf", width=3, height=3.5)
print(gc1)
dev.off()

cor_pval = matrix(nrow=nrow(cormat), ncol=ncol(cormat))
for(i in 1:nrow(cor_pval)){
  for(j in 1:ncol(cor_pval)){
    xi = pval_carseq[,2+i]
    yj = pval_snseq[,1+j]
    
    cor_pval[i,j] = cor.test(xi, yj, alternative = "greater", 
                             method = "spearman", use="pair")$p.value
  }
}

rownames(cor_pval) = names(pval_carseq)[-(1:2)]
colnames(cor_pval) = names(pval_snseq)[-(1)]

summary(c(cor_pval))
sort(c(cor_pval))[1:10]
cor_pval[which(cor_pval < 1e-5)] = 1e-5

gc2 = ggcorrplot(t(-log10(cor_pval)), tl.cex = 6, 
                 title="-log10(p-value)") + 
  scale_fill_gradient2(limit = c(0,5.01), low = "blue", high =  "red", 
                       mid = "white", midpoint = 1) 

pdf("../figures/step_z2_cor_pval_CARseq_vs_snRNAseq.pdf", width=3, height=3.5)
print(gc2)
dev.off()

# ----------------------------------------------------------------------
# check the association using fisher exact test
# ----------------------------------------------------------------------

pcut = 0.05
fisher_pval = matrix(nrow=nrow(cormat), ncol=ncol(cormat))

for(i in 1:nrow(cor_pval)){
  for(j in 1:ncol(cor_pval)){
    xi = pval_carseq[,2+i]
    yj = pval_snseq[,1+j]
    f1 = fisher.test(xi < pcut, yj < pcut, alternative = "greater")
    fisher_pval[i,j] = f1$p.value
  }
}

rownames(fisher_pval) = names(pval_carseq)[-(1:2)]
colnames(fisher_pval) = names(pval_snseq)[-(1)]

summary(c(fisher_pval))
sort(c(fisher_pval))[1:10]
fisher_pval[which(fisher_pval < 1e-4)] = 1e-4

gc2 = ggcorrplot(t(-log10(fisher_pval)), tl.cex = 6) + 
  scale_fill_gradient2(limit = c(0,4.01), low = "blue", high =  "red", 
                       mid = "white", midpoint = 1) 

pdf("../figures/step_z2_fisher_pval_CARseq_vs_snRNAseq.pdf", 
    width=3, height=3.5)
print(gc2)
dev.off()

# ----------------------------------------------------------------------
# read in the correlation between bulk RNA-seq and pseudo bulk
# ----------------------------------------------------------------------

sn_cr = fread("../data/sn_bulk_corr.txt")
dim(sn_cr)
sn_cr[1:5,]

table(sn_cr$r > 0.5)
table(sn_cr$r > 0.3)

table(pval_snseq$gene_name == gene2use)
table(pval_carseq$gene_name == gene2use)

sn_cr_v = rep(NA, length(gene2use))
match_gene = match(gene2use, sn_cr$gene)
table(is.na(match_gene))
sn_cr_v[!is.na(match_gene)] = sn_cr$r[match_gene[!is.na(match_gene)]]

table(sn_cr_v >= 0.31)
w2kp = which(sn_cr_v >= 0.31)

# ----------------------------------------------------------------------
# check the overall correlation matrix
# ----------------------------------------------------------------------

cormat = cor(pval_carseq[w2kp,-(1:2)], pval_snseq[w2kp,-1], 
             method = "spearman", use="pair")
summary(c(cormat))

gc1 = ggcorrplot(t(cormat), tl.cex = 6) + 
  scale_fill_gradient2(limit = c(-0.08,0.08), low = "blue", high =  "red", 
                       mid = "white", midpoint = 0) 

pdf("../figures/step_z2_cor_CARseq_vs_snRNAseq_selected_genes.pdf", 
    width=3, height=3.5)
print(gc1)
dev.off()

cor_pval = matrix(nrow=nrow(cormat), ncol=ncol(cormat))
for(i in 1:nrow(cor_pval)){
  for(j in 1:ncol(cor_pval)){
    xi = pval_carseq[w2kp,2+i]
    yj = pval_snseq[w2kp,1+j]
    
    cor_pval[i,j] = cor.test(xi, yj, alternative = "greater", 
                             method = "spearman", use="pair")$p.value
  }
}

rownames(cor_pval) = names(pval_carseq)[-(1:2)]
colnames(cor_pval) = names(pval_snseq)[-(1)]

summary(c(cor_pval))
sort(c(cor_pval))[1:10]
cor_pval[which(cor_pval < 1e-3)] = 1e-3

gc2 = ggcorrplot(t(-log10(cor_pval)), tl.cex = 6) + 
  scale_fill_gradient2(limit = c(0,3.01), low = "blue", high =  "red", 
                       mid = "white", midpoint = 1) 

pdf("../figures/step_z2_cor_pval_CARseq_vs_snRNAseq_selected_genes.pdf", 
    width=3, height=3.5)
print(gc2)
dev.off()

# ----------------------------------------------------------------------
# check the association using fisher exact test
# ----------------------------------------------------------------------

pcut = 0.05
fisher_pval = matrix(nrow=nrow(cormat), ncol=ncol(cormat))

for(i in 1:nrow(fisher_pval)){
  for(j in 1:ncol(fisher_pval)){
    xi = pval_carseq[w2kp,2+i]
    yj = pval_snseq[w2kp,1+j]
    f1 = fisher.test(xi < pcut, yj < pcut, alternative = "greater")
    fisher_pval[i,j] = f1$p.value
  }
}

rownames(fisher_pval) = names(pval_carseq)[-(1:2)]
colnames(fisher_pval) = names(pval_snseq)[-(1)]

summary(c(fisher_pval))
sort(c(fisher_pval))[1:10]
fisher_pval[which(fisher_pval < 1e-3)] = 1e-3

gc2 = ggcorrplot(t(-log10(fisher_pval)), tl.cex = 6) + 
  scale_fill_gradient2(limit = c(0,3.01), low = "blue", high =  "red", 
                       mid = "white", midpoint = 1) 

pdf("../figures/step_z2_fisher_pval_CARseq_vs_snRNAseq_selected_genes.pdf", 
    width=3, height=3.5)
print(gc2)
dev.off()

# ----------------------------------------------------------------------
# check the association within carseq using fisher exact test
# ----------------------------------------------------------------------

pcut = 0.05
fisher_carseq = matrix(nrow=nrow(cor_carseq), ncol=ncol(cor_carseq))

for(i in 1:nrow(fisher_carseq)){
  for(j in 1:ncol(fisher_carseq)){
    xi = pval_carseq[,2+i]
    yj = pval_carseq[,2+j]
    f1 = fisher.test(xi < pcut, yj < pcut, alternative = "greater")
    fisher_carseq[i,j] = f1$p.value
  }
}

rownames(fisher_carseq) = names(pval_carseq)[-(1:2)]
colnames(fisher_carseq) = names(pval_carseq)[-(1:2)]

summary(c(fisher_carseq))
sort(c(fisher_carseq[upper.tri(fisher_carseq)]))[1:20]
fisher_carseq[which(fisher_carseq < 1e-10)] = 1e-10

gc2 = ggcorrplot(t(-log10(fisher_carseq)), tl.cex = 6) + 
  scale_fill_gradient2(limit = c(0,10.01), low = "blue", high =  "red", 
                       mid = "white", midpoint = 2) 

pdf("../figures/step_z2_fisher_pval_CARseq.pdf", 
    width=4.5, height=3.5)
print(gc2)
dev.off()


# ----------------------------------------------------------------------
# check the association within snSeq using fisher exact test
# ----------------------------------------------------------------------

pcut = 0.05
fisher_snseq = matrix(nrow=nrow(cor_snseq), ncol=ncol(cor_snseq))

for(i in 1:nrow(fisher_snseq)){
  for(j in 1:ncol(fisher_snseq)){
    xi = pval_snseq[,1+i]
    yj = pval_snseq[,1+j]
    f1 = fisher.test(xi < pcut, yj < pcut, alternative = "greater")
    fisher_snseq[i,j] = f1$p.value
  }
}

rownames(fisher_snseq) = names(pval_snseq)[-(1)]
colnames(fisher_snseq) = names(pval_snseq)[-(1)]

summary(c(fisher_snseq))
sort(c(fisher_snseq[upper.tri(fisher_snseq)]))[1:10]

fisher_snseq_bounded = fisher_snseq
fisher_snseq_bounded[which(fisher_snseq < 1e-80)] = 1e-80

gc2 = ggcorrplot(t(-log10(fisher_snseq_bounded)), tl.cex = 6) + 
  scale_fill_gradient2(limit = c(0,80.01), low = "blue", high =  "red", 
                       mid = "white", midpoint = 2) 

pdf("../figures/step_z2_fisher_pval_snseq.pdf", 
    width=3, height=2.4)
print(gc2)
dev.off()

# ----------------------------------------------------------------------
# check a few cases
# ----------------------------------------------------------------------

w1 = which(fisher_snseq < 1e-100 & fisher_snseq > 0, arr.ind = TRUE)
cts = rownames(fisher_snseq)
df1 = data.frame(ct1=cts[w1[,1]], ct2=cts[w1[,2]], 
                 fisher_pval=fisher_snseq[w1], stringsAsFactors=FALSE)
df1[order(df1$fisher_pval),]

summary(pval_snseq$Inh)
summary(pval_snseq$Exc)

gs1 = ggplot(pval_snseq,aes(x=-log10(Exc),y=-log10(Inh))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c()

pdf("../figures/step_z2_snseq_Exc_Inh.pdf", 
    width=4.5, height=3.5)
print(gs1)
dev.off()

c1 = chisq.test(pval_snseq$Inh < pcut, pval_snseq$Exc < pcut)
c1$p.value
c1$expected
c1$observed

gc()

sessionInfo()
q(save="no")


