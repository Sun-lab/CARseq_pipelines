
# ========================================================================
# libraries and path
# ========================================================================

library(MASS)
library(data.table)
library(foreach)
library(doParallel)
library(qvalue)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(fgsea)
library(stringr)
library(ggcorrplot)

# ------------------------------------------------------------------------
# read in pathway information
# ------------------------------------------------------------------------

gmtfile_reactome  = "../data/c2.cp.reactome.v7.1.symbols.gmt"
pathways_reactome = gmtPathways(gmtfile_reactome)

reactome_genes = unique(unlist(pathways_reactome))
length(reactome_genes)

# ------------------------------------------------------------------------
# read in p-values
# ------------------------------------------------------------------------

pval_carseq = readRDS("../results/step_z2_pval_carseq.rds")
pval_snseq  = readRDS("../results/step_z2_pval_snseq.rds")

dim(pval_carseq)
dim(pval_snseq)
pval_carseq[1:2,]
pval_snseq[1:2,]

table(pval_carseq$gene_name == pval_snseq$gene_name)
genes = pval_carseq$gene_name
length(genes)
length(unique(genes))

# ------------------------------------------------------------------------
# filter genes in reactome annoations
# ------------------------------------------------------------------------

length(pathways_reactome)
summary(sapply(pathways_reactome, length))

for(p1 in names(pathways_reactome)){
  genes_p1 = intersect(pathways_reactome[[p1]], genes)
  if(length(genes_p1) < 10 || length(genes_p1) > 1000){
    pathways_reactome[[p1]] = NULL
  }else{
    pathways_reactome[[p1]] = genes_p1
  }
}

length(pathways_reactome)
summary(sapply(pathways_reactome, length))

# ------------------------------------------------------------------------
# GESA
# ------------------------------------------------------------------------

pvals = cbind(pval_carseq[,-(1:2)], pval_snseq[,-1])
rownames(pvals) = pval_carseq$gene_name

dim(pvals)
pvals[1:2,1:4]
names(pvals)

gsea = list()
for(i in 1:ncol(pvals)){
  cat(i, date(), "\n")
  
  stats = -log10(pvals[, i])
  names(stats) =  rownames(pvals)
  length(stats)
  
  stats = na.omit(stats)
  length(stats)

  set.seed(1234)
  fgseaRes = fgseaMultilevel(pathways_reactome, stats, minSize=10, 
                             maxSize=1000)
  od1 = order(fgseaRes[,"padj"], -fgseaRes[,"NES"])
  fgseaRes = fgseaRes[od1,]
  gsea[[names(pvals)[i]]] = fgseaRes
}

lapply(gsea, dim)
gsea$CARseq_Micro[1:2,]

# ------------------------------------------------------------------------
# check the overlap of top pathways
# ------------------------------------------------------------------------

pathways_nms = names(pathways_reactome)

fisher_pval = matrix(NA, nrow=ncol(pvals), ncol=ncol(pvals))
rownames(fisher_pval) = names(pvals)
colnames(fisher_pval) = names(pvals)

for(i in 1:nrow(fisher_pval)){
  label_i = names(pvals)[i]
  gsea_i  = gsea[[label_i]]
  gsea_i  = gsea_i[which(gsea_i$NES > 0 & gsea_i$padj < 0.2),]
  xi = gsea_i$pathway
  if(length(xi) <= 1) { next }
  
  for(j in 1:ncol(fisher_pval)){
    label_j = names(pvals)[j]
    gsea_j  = gsea[[label_j]]
    gsea_j  = gsea_j[which(gsea_j$NES > 0 & gsea_j$padj < 0.2),]
    yj = gsea_j$pathway
    if(length(yj) <= 1) { next }
    
    f1 = fisher.test(pathways_nms %in% xi, 
                     pathways_nms %in% yj, alternative = "greater")
    fisher_pval[i,j] = f1$p.value
  }
}

dim(fisher_pval)
summary(c(fisher_pval))
rowSums(is.na(fisher_pval[1:13,14:19]))
colSums(is.na(fisher_pval[1:13,14:19]))

sort(c(fisher_pval[upper.tri(fisher_pval)]))[1:30]
fisher_pval[which(fisher_pval < 1e-10)] = 1e-10

gc2 = ggcorrplot(t(-log10(fisher_pval[1:13,14:19])), tl.cex = 6) + 
  scale_fill_gradient2(limit = c(0,10.01), low = "blue", high =  "red", 
                       mid = "white", midpoint = 2) 

pdf("../figures/step_z3_gsea_top10_pathways_overlap_fisher_pval.pdf", 
    width=5, height=3.5)
print(gc2)
dev.off()

sessionInfo()
q(save="no")
