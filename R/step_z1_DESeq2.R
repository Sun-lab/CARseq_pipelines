
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
theme_set(theme_classic2())
library(ggcorrplot)
library(qvalue)
library(ggpointdensity)
library(DESeq2)
library(dplyr)

path_ideas = "../../ideas/Autism"


qqp <- function(pvals, main, confidence=.95, cutoff=1){
  
  alpha = 1-confidence
  n     = length(pvals)
  
  pvals[is.na(pvals)]=1
  pvals=sort(pvals)
  
  k=c(1:n)
  
  lower = cutoff*qbeta(alpha/2, k, n+1-k)
  upper = cutoff*qbeta((1-alpha/2), k, n+1-k)
  
  expected = cutoff*k/(n+1)
  n0 = length(which(pvals ==0))
  
  if(n0 > 0){
    warning(sprintf("there are %d p-values being 0\n", n0))
  }
  
  biggest= max(-log10(pvals[which(pvals > 0)]), -log10(expected))
  
  plot(-log10(expected), -log10(pvals), xlim=c(0,biggest),
       ylim=c(0,biggest), pch=20, xlab="-log10(expected p-value)",
       ylab="-log10(observed p-value)", cex=0.6, bty="n", main=main)
  
  lines(-log10(expected), -log10(lower), lty=2)
  lines(-log10(expected), -log10(upper), lty=2)
  
}


# ----------------------------------------------------------------------
# check the number of cells and gene expression levels
# ----------------------------------------------------------------------

cell_info = fread(file.path(path_ideas, "data/meta.tsv"))
dim(cell_info)
cell_info[1:2,]

table(cell_info$region)

cell_info = cell_info[which(cell_info$region=="PFC"),]

sort(table(paste(cell_info$diagnosis, cell_info$sample, sep=":")))

p1 = ggplot(cell_info, aes(x=cluster, y=log10(UMIs), color=diagnosis)) +
  geom_boxplot(outlier.size = 0.3) + coord_flip() 

pdf("../figures/step_z1_n_UMI_per_cluster.pdf", 
    width=5, height=6)
print(p1)
dev.off()

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
# ------------------------------------------------------------------------

cts = list.files(file.path(path_ideas, "data/ct_mtx/"))
cts = sort(gsub(".rds", "", cts))
cts

n_zeros = list()
cell_rd = list()
n_cells = NULL

for(grp in cts){
  dat1 = readRDS(file.path(path_ideas, sprintf("data/ct_mtx/%s.rds", grp)))
  dim(dat1)
  class(dat1)
  dat1[1:5,1:4]
  
  n_zeros[[grp]] = rowSums(dat1 == 0)
  cell_rd[[grp]] = colSums(dat1)
  n_cells = c(n_cells, ncol(dat1))
}

table(sapply(n_zeros, length))
sapply(cell_rd, median)
sapply(cell_rd, mean)

n_zeros = as.data.frame(n_zeros)
dim(n_zeros)
n_zeros[1:2,1:5]

names(n_cells) = cts
sort(n_cells)

percent_zeros = t(t(n_zeros)/n_cells)
dim(percent_zeros)
percent_zeros[1:2,1:5]

sort(colSums(percent_zeros < 0.8))
sort(round(colSums(percent_zeros < 0.8)/nrow(percent_zeros),2))

for(i in 1:5){ gc() }
gc()

# ----------------------------------------------------------------------
# conduct DESeq2 analysis
# ----------------------------------------------------------------------

n_zeros = q75_trec = list()
 
n_cells = NULL

cts

ct_grps = list()
ct_grps[["Astro"]] = cts[1:2]
ct_grps[["Inh"]]   = cts[4:7]
ct_grps[["Exc"]]   = cts[8:11]
ct_grps[["Micro"]] = cts[12]
ct_grps[["Oligo"]] = cts[16]
ct_grps[["OPC"]]   = cts[17]

ct_grps

for(ct1 in names(ct_grps)){
  grps  = ct_grps[[ct1]]
  
  cat(ct1, date(), "\n")
  
  dat1 = NULL
  for(grp1 in grps){
    d1 = readRDS(file.path(path_ideas, sprintf("data/ct_mtx/%s.rds", grp1)))
    dat1 = cbind(dat1, d1)
  }
  
  dim(dat1)
  dat1[1:2,1:4]
  
  n_zeros[[ct1]] = rowSums(dat1 == 0)
  n_cells = c(n_cells, ncol(dat1))
  
  cat(sprintf("there are %d cells\n", ncol(dat1)))
  
  stopifnot(all(colnames(dat1) %in% cell_info$cell))
  meta = cell_info[match(colnames(dat1), cell_info$cell),]
  dim(meta)
  meta[1:2,]
  
  meta_ind = distinct(meta[,3:12])
  dim(meta_ind)
  meta_ind[1:2,]
  names(meta_ind)[9:10] = c("PMI", "RIN")
  
  if(nrow(meta_ind) != length(unique(meta$individual))){
    stop("there is non-unique information\n")
  }
  
  table(meta_ind$Seqbatch, meta_ind$Capbatch)
  
  # ------------------------------------------------------------------------
  # collect count data
  # ------------------------------------------------------------------------
  
  trec1 = matrix(NA, nrow=nrow(dat1), ncol=nrow(meta_ind))
  colnames(trec1) = meta_ind$sample
  rownames(trec1) = rownames(dat1)
  dim(trec1)
  trec1[1:2,1:3]
  
  for(i in 1:ncol(trec1)){
    wi = which(meta$sample == meta_ind$sample[i])
    trec1[,i] = rowSums(dat1[,wi])
  }
  
  dim(trec1)
  trec1[1:2,1:3]
  
  summary(apply(trec1, 1, median))
  q75 = apply(trec1, 1, quantile, probs=0.75)
  summary(q75)
  table(q75 >= 20)
  
  q75_trec[[ct1]] = q75
  
  # ------------------------------------------------------------------------
  # run DESeq2
  # ------------------------------------------------------------------------
  
  colData = meta_ind
  for(i in 1:ncol(colData)){
    if(is.character(colData[[i]])){
      colData[[i]] = as.factor(colData[[i]])
    }
  }
  dim(colData)
  colData[1:2,]
  summary(colData)
  
  colData$diagnosis = factor(colData$diagnosis, levels=c("Control", "ASD"))
  
  
  dds = DESeqDataSetFromMatrix(countData = trec1, 
                               colData = colData,
                               design = ~ age + sex + Seqbatch + RIN + diagnosis)
  dds = DESeq(dds)
  
  res = results(dds)
  dim(res)
  res[1:2,]
  summary(res)
  
  res  = as.data.frame(res)
  n_zeros_ct1 = rowSums(trec1 == 0)
  
  res[["n_zeros"]] = n_zeros_ct1
  res[["percent_zeros"]] = n_zeros_ct1/ncol(trec1)
  dim(res)
  res[1:2,]
  
  nms = resultsNames(dds)
  nms
  nms = nms[-1]
  
  pvals2 = matrix(NA, nrow=nrow(trec1), ncol=length(nms))
  
  for(k in 1:length(nms)){
    rk = results(dds, name=nms[k])
    pvals2[,k] = rk$pvalue
  }
  
  colnames(pvals2) = nms
  dim(pvals2)
  cat("summary of p-values\n")
  print(summary(pvals2))
  
  
  png(sprintf("../figures/stepz1_DESeq2/%s_pval_hist.png", ct1), 
      width=7.5, height=5, units="in", res=400)
  par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
  for(k in 1:length(nms)){
    hist(pvals2[,k], main=nms[k], xlab="p-value", breaks=50)
  }
  
  mm = "diagnosis_ASD_vs_Control"
  qqp(pvals2[,mm], main=mm)
      
  dev.off()
  
  write.table(res, file=sprintf("../results/stepz1_DESeq2_%s.txt", ct1), 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  gc()
}


# ----------------------------------------------------------------------
# sumarize number of zeros at cell level and the 75th percentile after 
# collpasing cell level counts as individual level counts
# ----------------------------------------------------------------------

table(sapply(n_zeros, length))
table(sapply(q75_trec, length))

n_zeros = as.data.frame(n_zeros)
dim(n_zeros)
n_zeros[1:2,]

q75_trec = as.data.frame(q75_trec)
dim(q75_trec)
q75_trec[1:2,]

names(n_cells) = names(ct_grps)
sort(n_cells)

percent_zeros = t(t(n_zeros)/n_cells)
dim(percent_zeros)
percent_zeros[1:2,1:5]

sort(colSums(percent_zeros < 0.8))
sort(round(colSums(percent_zeros < 0.8)/nrow(percent_zeros),2))

pdf("../figures/step_z1_snseq_percent_zero_across_cells.pdf", width=7.5, height=5)
par(mfrow=c(2,3), mar=c(5,4,2,1))
for(k in 1:6){
  hist(1 - percent_zeros[,k], main=names(ct_grps)[k], breaks=20, 
       xlab="proportion of cells with any expression")
}
dev.off()


pdf("../figures/step_z1_snseq_percent_q75_across_individuals.pdf", width=7.5, height=5)
par(mfrow=c(2,3), mar=c(5,4,2,1))
for(k in 1:6){
  hist(log10(q75_trec[,k]+1), main=names(ct_grps)[k], breaks=seq(0,6,by=0.2), 
       xlab="log10(75-th percentile of read counts +1)")
}
dev.off()

gc()

sessionInfo()
q(save="no")


