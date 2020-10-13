
library(data.table)
library(stringr)

library(ggplot2)
library(ggpubr)
library(ggpointdensity)
theme_set(theme_classic2())
library(ggcorrplot)

path_ideas = "../../ideas/Autism"

# ----------------------------------------------------------------------
# read in p-values for ASD
# ----------------------------------------------------------------------

pval.path = "../results/_pvalues"

pval.files = list.files(path=pval.path, pattern="ASD", full.names=TRUE)
pval.files

f1 = pval.files[1]
pval.ASD = fread(f1)
lb = str_extract(basename(f1), "(?<=ASD_)\\w+_\\w+(?=.txt)")
names(pval.ASD)[3] = lb

dim(pval.ASD)
pval.ASD[1:2,]

length(unique(pval.ASD$gene_id))
length(unique(pval.ASD$gene_name))

for(f1 in pval.files[-1]){
  d1 = fread(f1)
  if(any(d1$gene_id != pval.ASD$gene_id)){
    stop("gene_id do not match\n")
  }
  
  bn = basename(f1)
  lb = str_extract(bn, "(?<=ASD_)\\w+_\\w+(?=.txt)")
  
  pval.ASD[[lb]] = d1$pvalue
}

dim(pval.ASD)
pval.ASD[1:2,]

# ----------------------------------------------------------------------
# read in p-values from snRNA-seq analysis
# ----------------------------------------------------------------------

pval_mtrx = data.frame(pval.ASD[,-(1:2)])
rownames(pval_mtrx) = pval.ASD$gene_id

dim(pval_mtrx)
pval_mtrx[1:2,]


grp = "PFC_L2_3"
fnm = sprintf("%s/res/step1_DESeq2_%s_adj_covariates.txt", path_ideas, grp)

deseq2 = read.table(file=fnm, sep="\t", header=TRUE)
dim(deseq2)
deseq2[1:2,]

gene2use = intersect(pval.ASD$gene_name, rownames(deseq2))
length(gene2use)

mat1 = match(gene2use, pval.ASD$gene_name)
mat2 = match(gene2use, rownames(deseq2))

pval_mtrx[[grp]] = rep(NA, nrow(pval.ASD))
pval_mtrx[[grp]][mat1] = deseq2$pvalue[mat2]

min(pval_mtrx, na.rm = TRUE)
pval_mtrx = -log10(pval_mtrx)

dim(pval_mtrx)
pval_mtrx[1:2,]

cor(pval_mtrx[,1:6],  pval_mtrx[,"DESeq2_bulk"], use="pairwise")
cor(pval_mtrx[,8:14], pval_mtrx[,"DESeq2_bulk"], use="pairwise")

cor(pval_mtrx[,1:6],  pval_mtrx[,grp], use="pairwise")
cor(pval_mtrx[,8:13], pval_mtrx[,grp], use="pairwise")

cor(pval_mtrx[,1:6],  pval_mtrx[,grp], use="pairwise", method="spearman")
cor(pval_mtrx[,8:13], pval_mtrx[,grp], use="pairwise", method="spearman")

gc()

sessionInfo()
q(save="no")


