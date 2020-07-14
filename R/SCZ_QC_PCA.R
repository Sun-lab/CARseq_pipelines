
# ----------------------------------------------------------------------
# QC_PCA: Quality Control and PCA on Total Read Counts
# ----------------------------------------------------------------------

rm(list=ls())
# repo_dir = "~/github/CSeQTL"
# cmc_dir  = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/CMC"
# repo_dir = "~/github/CARseq_pipelines"
# s3_dir = "/datastore/hayeslab/Current_members/Paul_Little/CS_eQTL/s3_Real"
# cmc_dir = file.path(s3_dir,"CMC")
# scRNA_dir = file.path(s3_dir,"scRNAseq_pipelines")

repo_dir = "~/research/GitHub/CARseq_pipelines"
cmc_dir  = "~/research/data/CMC"

# ----------------------------------------------------------------------
# Library/Source Basic Functions
# ----------------------------------------------------------------------

source(file.path(repo_dir,"R","base_functions_plittle.R"))
library(data.table)
library(sva)
# library(biomaRt)
# library(Rcpp); library(RcppArmadillo)
# sourceCpp(file.path(repo_dir,"src","QC_PCA.cpp"),showOutput = TRUE)

CMC_clinical_prep = function(cmc_dir){

	# Get DATA
	out = readRDS(file.path(cmc_dir,"output","DATA_TREC.rds"))
	# lapply(out, dim)
	dat = out$DATA; 
	cat("initial dim of dat from TReC data:", dim(dat), "\n")
	dat[1:3,]
	# trec = out$TREC; trec[1:5,1:4]

	# Key file of subject ids
	cat("Import mapping between sample RNAseq, Individual, and genotype ids...\n")
	key_fn = file.path(cmc_dir,"CMC_Human_sampleIDkey_metadata")
	key = smart_RT(key_fn,sep=",",header=TRUE)
	key = name_change(key,"RNAseq.Sample_RNA_ID","RNAseq_sample_id")
	key = name_change(key,"Genotypes.Genotyping_Sample_ID","Genotype_sample_id")
	dat = smart_merge(dat,key,all.x = TRUE)
	dim(key); key[1:2,]
	dim(dat); dat[1:2,]

	# Clinical variables
	cat("Import and check clinical variables...\n")
	# Source: https://www.synapse.org/#!Synapse:syn3354385 
	clin_fn = file.path(cmc_dir,"CMC_Human_clinical_metadata.csv")
	clin = smart_RT(clin_fn,sep = ",",header = TRUE)
	dim(clin); clin[1:3,]
	clin = name_change(clin,"Individual.ID","Individual_ID")
	clin = name_change(clin,"Reported.Gender","gender")
	clin = name_change(clin,"Brain.Weight..in.grams.","brain_weight")
	clin = name_change(clin,"PMI..in.hours.","PMI") # Post Mortem Interval
	clin = name_change(clin,"Year.of.Autopsy","autopsy_year")
	clin = name_change(clin,"Age.of.Death","age_death")
	clin = name_change(clin,"Cause.of.Death","cause_death")
	clin = name_change(clin,"Height..Inches.","height")
	clin = name_change(clin,"Weight..pounds.","weight")
	clin$age_death[clin$age_death == "90+"] = 90
	clin$age_death = as.numeric(clin$age_death)
	dat = smart_merge(dat,clin,all.x = TRUE)
  table(clin$age_death==90)
  
	# RNAseq QC metrics
	cat("Import RNAseq QC file ...\n")
	# Source: https://www.synapse.org/#!Synapse:syn18358379
	qc_fn = file.path(cmc_dir,"CMC_Human_rnaSeq_metadata.csv")
	qc = smart_RT(qc_fn,sep=",",header=TRUE)
	dim(qc); qc[1:2,]

	qc = name_change(qc,"SampleID","RNAseq_sample_id")
	dat = smart_merge(dat,qc,all.x=TRUE)
	smart_table(!is.na(dat$RIN))
	smart_table(!is.na(dat$Study))

	cat("dim of dat after merging clinc and qc info:", dim(dat), "\n")
	
	# Check variables
	# sort(names(dat))
	print(sapply(c("Institution","Dx","gender","Sex", "Study", 
								 "Ethnicity","autopsy_year","cause_death","Hemisphere"),
							 function(xx) smart_table(dat[,xx])))
	options(width=100)
	print(t(sapply(c("brain_weight","height","weight",
									 "PMI","pH","age_death",
									 "Total_RNA_Yield","RIN","X28S.18S",
									 "Total_Reads","Mapped_Reads","Genes_Detected",
									 "Percent_Aligned", "rRNA_Rate"),
								 function(xx) summary(dat[,xx]))))

	# Before subsetting, try to replicate sample size from paper
	dim(dat)
	dat = dat[which(dat$Study == "CMC"),]
	dim(dat)
	
	cat("dim of dat after filter 'dat$Study == \"CMC\"':", dim(dat), "\n")
	
	# remove duplicaets based on paper supplements
	# "Of the 10 Pitt control samples that were sequenced twice, only 
	#  the first sequencing run was included in our analysis"
	tab = table(dat$Individual_ID)
	tab = tab[tab > 1]; tab

	tmp_df = dat[which(dat$Individual_ID %in% names(tab)),
							 c("Individual_ID","RNAseq_sample_id","RIN",
								 "Total_Reads","Mapped_Reads","Genes_Detected",
								 "RNA_Prep_Date")]
	tmp_df = tmp_df[order(tmp_df$Individual_ID,tmp_df$RNA_Prep_Date),]
	tmp_df$order = rep(c(1,2),10)
	tmp_df = tmp_df[which(tmp_df$order == 2),] # keep RNAseq_sample_id's to exclude
	dat = dat[which(!(dat$RNAseq_sample_id %in% tmp_df$RNAseq_sample_id)),]
	dim(dat)
	cat("dim of dat after filter duplicates:", dim(dat), "\n")
	
	smart_table(dat[,c("gender","Sex")])
	dat = dat[which(dat$Sex %in% c("XX","XY")),]
	dim(dat)
	cat("dim of dat after filter Sex as XX or XY:", dim(dat), "\n")
	
	smart_table(dat$Dx)
	smart_table(!is.na(dat$RIN))
	smart_table(dat[,c("Dx","Institution")])
	smart_table(dat[,c("Dx","gender")])
	
	## Create Library cluster based on CMC's Supplement
	dat$libclust = NA
	dat$libclust[which(dat$Library_Batch %in% c("1_11/11/13","1_11/26/13"))] = "base"
	dat$libclust[which(dat$Library_Batch %in% c("1_17","1_18","1_10","1_8/28/13"))] = "A"
	dat$libclust[which(dat$Library_Batch %in% c("1_14","1_7","1_11","1_6","1_16","1_25"))] = "B"
	dat$libclust[which(dat$Library_Batch %in% c("1_5","1_1","1_3"))] = "C"
	dat$libclust[which(dat$Library_Batch %in% c("1_2","1_28","1_8","1_12","1_4"))] = "D"
	dat$libclust[which(dat$Library_Batch %in% c("1_20","1_24","1_26","1_27","1_9"))] = "E"
	dat$libclust[which(dat$Library_Batch %in% c("1_13","1_15","1_22","1_21","1_23"))] = "F"
	dat$libclust[which(dat$Library_Batch %in% c("1_10/15/13","1_19"))] = "G"
	dat$libclust[which(dat$Library_Batch %in% c("1_10/9/13"))] = "H"
	smart_table(dat$libclust)
	dat$libclust[dat$libclust == "H"] = "G"
	smart_table(dat$libclust)
	
	dat
}


TREC_QC = function(TREC){
	cat("\nExploring gene statistics ...\n")
	rMin = apply(TREC,1,min)
	rMed = apply(TREC,1,median)
	r75  = apply(TREC,1,quantile,probs=0.75)
	r90  = apply(TREC,1,quantile,probs=0.90)
	
	par(mfrow=c(3,3), mar=c(4,4,1,1), bty="n")
	smart_hist(log10(1+rMin),xlab="log10(min + 1)",main="")
	smart_hist(log10(1+rMed),xlab="log10(median + 1)",main="")
	smart_hist(log10(1+r75),xlab="log10(75 percentile + 1)",main="")
	smart_hist(log10(1+r90),xlab="log10(90 percentile + 1)",main="")
	plot(log10(1+rMin),log10(1+rMed),xlab="log10(1+rMin)",ylab="log10(1+rMed)")
	plot(log10(1+r75),log10(1+rMed),xlab="log10(1+r75)",ylab="log10(1+rMed)")
	plot(log10(1+r90),log10(1+rMed),xlab="log10(1+r90)",ylab="log10(1+rMed)")
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,bty="o")

	print("smart_table(rMed >= 10)")
	print(smart_table(rMed >= 10))
	print("smart_table(r75 >= 20)")
	print(smart_table(r75 >= 20))
	
	print("summary(rMin[r75 >= 20])")
	print(summary(rMin[r75 >= 20]))
	print("summary(rMed[r75 >= 20])")
	print(summary(rMed[r75 >= 20]))
}


TREC_PCA = function(DATA,TREC,MODEL,PROP,SIG,GTF,show = TRUE){
	if(FALSE){
		DATA = dat2; TREC = trec; MODEL = MODEL_LOG_PROP;
		PROP = prop; SIG = sig; GTF = gtf; show = TRUE
	}
	
  int_subj = intersect(DATA$RNAseq_sample_id,colnames(TREC)); length(int_subj)
  TREC = TREC[,int_subj]
  DATA = DATA[which(DATA$RNAseq_sample_id %in% int_subj),]
  TREC = TREC[,match(DATA$RNAseq_sample_id,colnames(TREC))]
  
  # Log normalized counts
  TREC = t(log10(t(TREC + 1)/exp(DATA$log_depth))) 
  
  ## Sort subjects in PROP
  PROP = PROP[int_subj,]
  PROP = PROP[match(DATA$RNAseq_sample_id,rownames(PROP)),]
  
	# Make sure subjects and genes are correctly ordered using DATA
	## Get genes, use GTF to order genes
	dim(GTF)
	GTF = GTF[which(GTF$ensg %in% rownames(TREC) & GTF$gene %in% rownames(SIG)),]
	dim(GTF)
	
	int_gene_ensg = intersect(GTF$ensg,rownames(TREC)); length(int_gene_ensg)
	int_gene_hugo = intersect(GTF$gene,rownames(SIG)); length(int_gene_hugo)
	SIG = SIG[int_gene_hugo,]
	SIG = SIG[match(GTF$gene,rownames(SIG)),]
	
	TREC = TREC[int_gene_ensg,]
	TREC = TREC[match(GTF$ensg,rownames(TREC)),]

	stopifnot(all(DATA$RNAseq_sample_id == colnames(TREC)))
	stopifnot(all(DATA$RNAseq_sample_id == rownames(PROP)))
	stopifnot(all(GTF$gene == rownames(SIG)))
	stopifnot(all(GTF$ensg == rownames(TREC)))
	
	# Run raw PCA on normalized TREC
	cat("\nRun PCA on row-centered outcomes ...\n")
	TREC = TREC - rowMeans(TREC,na.rm = TRUE); TREC[is.na(TREC)] = 0
	cov_TREC = t(TREC) %*% TREC / nrow(TREC); dim(cov_TREC)
	pca_TREC = eigen(cov_TREC)
	
	if( show ) show_screeplot(pca_TREC,main = "trec PCA")
	pca_TREC$values[1:22]
	cumsum(pca_TREC$values)[1:22]/sum(pca_TREC$values)
	
	num_pcs = 5; pcs = smart_df(pca_TREC$vectors[,1:num_pcs])
	names(pcs) = paste0("PC",seq(num_pcs))
	if( show ) show_pc_color(pcs,submain = "trec PCA")
	
	# Check if model includes "CS"
	MODEL = paste(MODEL,collapse = " + ")
	include_CS = grepl("CS",MODEL)
	include_PROP = grepl("prop_Astro",MODEL)
	include_LOG_PROP = grepl("log_Astro",MODEL)
	
	# Looping over genes for lm() and residuals
	GG = nrow(TREC); NN = ncol(TREC)
	rr = matrix(NA,GG,NN) # residual matrix
	
	for(gg in seq(GG)){
	 	# gg = 1
	 	if(gg %% 1e2 == 0) cat(".")
	 	if(gg %% 2e3 == 0 || gg == GG) cat(sprintf("%s out of %s\n",gg,GG))
	 	# variables to adjust for: gender + institution + libclust + 
		#		age_death + PMI + RIN + genoPC + log_depth + log(rho^T mu_sig)
		DATA2 = DATA
		if( include_CS) DATA2$CS = log(c(PROP %*% SIG[gg,]))
		ctypes2use = which(colnames(PROP) != "Exc")
		
		if( include_PROP ) {
		  prop.gg = PROP[,ctypes2use]
		  colnames(prop.gg) = paste("prop", colnames(prop.gg), sep="_")
		  DATA2 = cbind(DATA2, prop.gg)
		}
		
		if( include_LOG_PROP ) {
		  prop.gg = log(PROP[,ctypes2use] + 1e-5)
		  prop.gg = prop.gg - log(PROP[,which(colnames(PROP) == "Exc")])
		  colnames(prop.gg) = paste("log", colnames(prop.gg), sep="_")
		  DATA2 = cbind(DATA2, prop.gg)
		}
		
		DATA2$TREC = TREC[gg,]
		lm_out = lm(formula(sprintf("TREC ~ %s",MODEL)),data = DATA2)
	 	rr[gg,] = as.numeric(lm_out$residuals)
		aa = drop1(lm_out,.~.,test="F")
		
		if(gg == 1){
	 		pp = matrix(NA,GG,nrow(aa)-1)
			colnames(pp) = rownames(aa)[-1]
	 	}
		
		pp[gg,] = aa[-1,6]
		rm(lm_out)
	}
	
	if( show ) show_pvalue_hist(mat_pvalues = pp,test_type0 = 3)
	rr2 = rr - rowMeans(rr,na.rm = TRUE)
	cov_rr2 = t(rr2) %*% rr2 / nrow(rr2); pca_rr2 = eigen(cov_rr2)
	if( show ) show_screeplot(pca_rr2,main = "")
	num_pcs = 5; pcs = smart_df(pca_rr2$vectors[,1:num_pcs])
	names(pcs) = paste0("PC",seq(num_pcs))
	if( show ) show_pc_color(pcs,submain = "")
	nADD_PCs = 20; add_pcs = smart_df(pca_rr2$vectors[,1:nADD_PCs])
	
	if( include_CS ){
		names(add_pcs) = paste0("CSadjustPC",seq(nADD_PCs))
	} else if(include_PROP) {
	  names(add_pcs) = paste0("PROPadjustPC",seq(nADD_PCs))
	}else if(include_LOG_PROP) {
	  names(add_pcs) = paste0("logPROPadjustPC",seq(nADD_PCs))
	}else{
		names(add_pcs) = paste0("noCSadjustPC",seq(nADD_PCs))
	}
	DATA = cbind(DATA,add_pcs)
	# print(head(DATA))
	
	DATA
}

## Check correlation between PCs and prop
plotTwo = function(MAT1,MAT2,...){
  if(FALSE){
    MAT1 = noCS_PCs; MAT2 = prop
  }
  
  num_rows = ncol(MAT1)
  num_cols = ncol(MAT2)
  par(mfrow=c(num_rows,num_cols),mar=c(5,4.2,1,0.5))
  for(ii in seq(num_rows)){
    for(jj in seq(num_cols)){
      xlab = ""; ylab = ""
      if(jj == 1){
        ylab = colnames(MAT1)[ii]
      }
      if(ii == num_rows){
        xlab = colnames(MAT2)[jj]
      }
      pears = cor(MAT2[,jj],MAT1[,ii],method="pears")
      spear = cor(MAT2[,jj],MAT1[,ii],method="spear")
      plot(x = MAT2[,jj],y = MAT1[,ii],xlab = xlab,
           ylab = ylab,bty = "n",
           main = sprintf("Pear=%s;Spear=%s",round(pears,2),round(spear,2)),
           ...)
    }}
  par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
}

# ----------------------------------------------------------------------
# Prepare clinical data: Merge key, clinical, and other metadata files
# ----------------------------------------------------------------------

dat = CMC_clinical_prep(cmc_dir = cmc_dir)
dat$RIN2 = dat$RIN*dat$RIN
dim(dat)
head(dat)

# log_depth will be recalculated later
dat = dat[,names(dat) != "log_depth"]

# ----------------------------------------------------------------------
# get outllier samples
# ----------------------------------------------------------------------

outlier.file = "CMC_MSSM-Penn-Pitt_DLPFC_mRNA_outlierSamples.txt"
outliers = read.table(file.path(cmc_dir, outlier.file), sep="\t", 
                      header=TRUE, as.is=TRUE)
outliers

table(outliers$Individual.ID %in% dat$Individual_ID)
dim(dat)

table(dat$Individual_ID  %in% outliers$Individual.ID)

dat = dat[which(! dat$Individual_ID %in% outliers$Individual.ID),]
dim(dat)

# ----------------------------------------------------------------------
# Import geno: a set of snps from chr1-chr22, and do PCA
# ----------------------------------------------------------------------

geno = readRDS(file.path(cmc_dir,"output","prune_geno.rds"))
dim(geno); geno[1:10,1:5]

## Run genotype PCA on all subjects
table(dat$Genotype_sample_id %in% colnames(geno))
int_subj = intersect(colnames(geno),dat$Genotype_sample_id)
length(int_subj)

dat 	= dat[which(dat$Genotype_sample_id %in% int_subj),]
geno 	= geno[,colnames(geno) %in% int_subj]
geno 	= geno[,match(dat$Genotype_sample_id,colnames(geno))]
rNA = rowSums(is.na(geno))
sum(rNA)/(nrow(geno)*ncol(geno))
summary(rNA/ncol(geno))
table(rNA/ncol(geno) > 0.05)

## Code 0/1/2/3 to 0/1/2 genotype
geno[geno==2] = 1; geno[geno==3] = 2

geno_c = geno - rowMeans(geno,na.rm = TRUE)
geno_c[is.na(geno_c)] = 0 # b/c there are missing genotypes
date()
cov_geno = t(geno_c) %*% geno_c / nrow(geno_c)
date()
pca_geno = eigen(cov_geno)

pdf(file.path(cmc_dir,"output","genotype_batches.pdf"),height=8,width=8)
show_screeplot(pca_geno,main = "Genotype PCA")
num_pcs = 5
pcs = smart_df(pca_geno$vectors[,1:num_pcs])
names(pcs) = paste0("PC",seq(num_pcs))
show_pc_color(pcs,submain = "Genotype PCA")
show_pc_color(pcs,DATA = dat,VAR = "Institution",CAT = TRUE,submain = "Genotype PCA")
show_pc_color(pcs,DATA = dat,VAR = "Ethnicity",CAT = TRUE,submain = "Genotype PCA")
show_pc_color(pcs,DATA = dat,VAR = "gender",CAT = TRUE,submain = "Genotype PCA")
dev.off()

## Retain top genotype PCs
num_geno_pcs = 5 # number of genotype pcs
pcs = smart_df(pca_geno$vectors[,1:num_geno_pcs])
names(pcs) = paste0("genoPC",seq(num_geno_pcs))
dat = cbind(dat,pcs)
rm(geno,geno_c,cov_geno,pca_geno,pcs)
dat[1:3,]

# ----------------------------------------------------------------------
# Import gene expressio data
# ----------------------------------------------------------------------

# Import trec: Results from "SO1", corresponds to summerizeOverlaps() with
# 	mode = "Union",singleEnd = FALSE,ignore.strand = TRUE,fragments = TRUE
rds = readRDS(file.path(cmc_dir,"output","DATA_TREC.rds"))
str(rds); names(rds); trec = rds$TREC; rm(rds)
dim(trec); trec[1:5,1:4]

int_subj  = intersect(dat$RNAseq_sample_id, colnames(trec)); length(int_subj)
trec  	  = trec[,which(colnames(trec) %in% int_subj)]
dat  		  = dat[which(dat$RNAseq_sample_id %in% int_subj),]
trec  	  = trec[,match(dat$RNAseq_sample_id,colnames(trec))]

dim(trec); trec[1:5,1:4]

# Run TReC QC on all genes
for(condition in c("Control","SCZ")){
  trec2 	  = trec[,which(dat$Dx == condition)]
  png(file.path(cmc_dir,"output",sprintf("expression_QC_%s.png",condition)),
      units='px',height=2500,width=2500,res=250,type='cairo',pointsize=20)
  TREC_QC(TREC = trec2)
  dev.off()
}

# ----------------------------------------------------------------------
# subset gene expressio data
# ----------------------------------------------------------------------

cat("\nSubset genes with gene75 >= 20 ...\n")
cat(sprintf("Number of genes before filtering = %s\n",nrow(trec)))

r75 = apply(trec,1,quantile,probs=0.75)
trec = trec[which(r75 >= 20),]

cat(sprintf("Number of genes after filtering = %s\n",nrow(trec)))

s75 = apply(trec,2,function(xx) quantile(xx,0.75))
dat$log_depth = as.numeric(log(s75)) # calculate log_depth
print(head(dat))

pdf(file.path(cmc_dir,"output","log_read_depth.pdf"), height=3,width=4)
par(mar=c(5,4,1,1))
smart_hist(dat$log_depth, breaks=40, 
           xlab="log(Subj75th Percentile depth)", main="")
dev.off()

# ----------------------------------------------------------------------
# Import cell type-specific gene expression data
# ----------------------------------------------------------------------

# Import gene/cell type signature matrix and cell type proportions
prop = readRDS(file.path(repo_dir,"MTG","prop_MTG.rds"))
names(prop)
prop = prop$ICeDT; dim(prop); head(prop)

int_subj = intersect(rownames(prop),dat$RNAseq_sample_id); length(int_subj)
prop = prop[which(rownames(prop) %in% int_subj),]
prop = prop[match(dat$RNAseq_sample_id,rownames(prop)),]
dim(prop); head(prop)

sig = readRDS(file.path(repo_dir,"MTG","MTG_sce_full.rds"))
dim(sig); head(sig)

# Need to map ENSG ids (in trec) with Hugo ids (or entrez id in sig)
# Also tried to use ensembl v75 to search the matching between ensembl 
# id and entrez id or gene name. v70 is used to annoate expression data
# but it is not available anymore.
#
# ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",
# host="feb2014.archive.ensembl.org")
# 
# however the number of genes that we can find is even smaller. 
# here since we only need the matching of genes to perform PCA
# it is ok that we miss some genes. 

gtf = data.table::fread(file.path(cmc_dir,"gencode.v29.annotation.gtf"),
	header = FALSE,data.table = FALSE)
dim(gtf); head(gtf); gtf = gtf[which(gtf$V3 == "gene"),]; dim(gtf); head(gtf)
gtf$ensg = sapply(gtf$V9,function(xx) gsub("\"","",strsplit(xx,"[ .;]")[[1]][2]),USE.NAMES=FALSE)
gtf$gene = sapply(gtf$V9,function(xx) gsub("\"","",strsplit(xx,"[ .;]")[[1]][9]),USE.NAMES=FALSE)
gtf = unique(gtf[,c("ensg","gene")]); rownames(gtf) = NULL
dim(gtf)
table(rownames(trec) %in% gtf$ensg)

tab = table(gtf$gene); tab = tab[tab > 1]; length(tab)
gtf = gtf[which(!(gtf$gene %in% names(tab))),]
dim(gtf)
gtf[1:10,]
table(rownames(trec) %in% gtf$ensg)

# ----------------------------------------------------------------------
# perform PCA
# ----------------------------------------------------------------------

## Specify model variables
ctypes = setdiff(colnames(prop), "Exc")

MODEL_noCS = c("gender","Institution","libclust","age_death","PMI",
               "RIN","RIN2", paste0("genoPC",1:5),"log_depth")
MODEL_CS   = c(MODEL_noCS,"CS")
MODEL_PROP = c(MODEL_noCS, paste("prop", ctypes, sep="_"))
MODEL_LOG_PROP = c(MODEL_noCS, paste("log", ctypes, sep="_"))

modelList = list(MODEL_noCS, MODEL_CS, MODEL_PROP, MODEL_LOG_PROP)
names(modelList) = c("noCS", "CS", "PROP", "logPROP")

for(condition in c("Control","SCZ")){
  cat("\n------------------------------------------------\n")
  cat(condition)
  cat("\n------------------------------------------------\n")
  
  dat2  = dat[which(dat$Dx == condition),]
  prop2 = prop[which(dat$Dx == condition),]
  
  for(mk in names(modelList)){
    cat("\n------------------------------------------------\n")
    cat(mk)
    cat("\n------------------------------------------------\n")
    
    pdf(file.path(cmc_dir,"output",sprintf("expression_PCA_%s_%s.pdf", condition, mk)),
        height = 8,width = 8)
    out1 = TREC_PCA(DATA = dat2, TREC = trec, MODEL = modelList[[mk]],
                        PROP = prop, SIG = sig, GTF = gtf, show = TRUE)
    dev.off()
    
    datPCs = out1[,grep(paste0("^", mk, "adjustPC"), colnames(out1))]
    dat2   = cbind(dat2, datPCs)

    cr.pearson  = round(cor(datPCs,  prop2, method = "pearson"),2)
    cr.spearman = round(cor(datPCs,  prop2, method = "spearman"),2)
    max.pearson  = apply(abs(cr.pearson), 1, max)
    max.spearman = apply(abs(cr.spearman), 1, max)
    
    pv.pearson = pv.spearman = rep(NA, ncol(datPCs))
    for(i in 1:ncol(datPCs)){
      pv.pearson[i] = cor.test(datPCs[,i], prop2[,which.max(abs(cr.pearson[i,]))], 
                               method="pearson")$p.value
      pv.spearman[i] = cor.test(datPCs[,i], prop2[,which.max(abs(cr.spearman[i,]))], 
                                method="spearman")$p.value
    }
    
    print(cbind(cr.pearson, max.pearson, signif(pv.pearson,2)))
    print(cbind(cr.spearman, max.spearman, signif(pv.spearman,2)))
    
    pdf(file.path(cmc_dir,"output",sprintf("compare_PCs_%s_Prop_%s.pdf", mk, condition)),
        height=12,width=15)
    
    plotTwo(MAT1 = datPCs[,1:5], MAT2 = prop2, col = rgb(0,0,0,0.5), pch = 16,
            cex.axis = 1.3,cex.lab = 1.5)
    plotTwo(MAT1 = datPCs[,6:10], MAT2 = prop2, col = rgb(0,0,0,0.5), pch = 16,
            cex.axis = 1.3,cex.lab = 1.5)
    
    dev.off()
    
  }
  
  ## Retain PCs
  
  write.table(dat2, file.path(cmc_dir,"output",sprintf("dat_post_QC_%s.tsv",condition)),
              sep = "\t",row.names = FALSE,quote = FALSE)
  
}

# ----------------------------------------------------------------------
# write out gene expression data
# ----------------------------------------------------------------------

dim(trec)
trec[1:2,1:5]

dim(dat)
dat[1:2,]

table(colnames(trec) == dat$RNAseq_sample_id)

saveRDS(trec, file=file.path(cmc_dir, "output/trec_filtered.rds"))

saveRDS(dat, file=file.path(cmc_dir, "output/dat_cavariates.rds"))

log.trec = t(log10(t(trec + 1)/exp(dat$log_depth))) 
dim(log.trec)
log.trec[1:2,1:5]

saveRDS(log.trec, 
        file=file.path(cmc_dir, "output/trec_filtered_log_rd_corrected.rds"))

# ----------------------------------------------------------------------
# compare with gene expression from CMC
# https://www.synapse.org/#!Synapse:syn3346749
# ----------------------------------------------------------------------

cmcFnm = "CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv"
cmcDat = fread(file.path(cmc_dir, cmcFnm))
dim(cmcDat)
cmcDat[1:2,1:5]

table(colnames(trec) %in% names(cmcDat))
cmcDat1 = cmcDat[,colnames(trec), with = FALSE]
dim(cmcDat1)

table(rownames(trec) %in% cmcDat$V1)
cmcDat1 = cmcDat1[match(rownames(trec), cmcDat$V1),]
dim(cmcDat1)
cmcDat1[1:2,1:5]

cmcDat1 = data.matrix(cmcDat1)
dim(cmcDat1)
cmcDat1[1:2,1:5]

mean.trec = rowMeans(trec)
mean.cmc  = rowMeans(cmcDat1)

sd.trec = apply(trec, 1, sd)
sd.cmc  = apply(cmcDat1, 1, sd)

summary(sd.trec)
summary(sd.cmc)
table(sd.cmc==0)

figNm = file.path(cmc_dir, "output/compare_expression_trec_cmc.png")

png(figNm, height=8, width=8, res=250, units="in")
par(mfrow=c(2,2), mar=c(5,6,1,1), bty="n", cex=0.5, cex.axis=2, 
    cex.lab=2, col=rgb(1,0.2,0.2,0.8))

plot(log10(mean.trec+1), log10(sd.trec+1))
abline(0,1, col="blue")

plot(log10(mean.cmc+1), log10(sd.cmc+1))
abline(0,1, col="blue")

plot(log10(mean.trec+1), log10(mean.cmc+1))
abline(0,1, col="blue")

plot(log10(sd.trec+1), log10(sd.cmc+1))
abline(0,1, col="blue")

dev.off()

cr1 = rep(NA, nrow(trec))
w2check = which(sd.trec > 0 & sd.cmc > 0)
for(i in w2check){
  cr1[i] = cor(cmcDat1[i,], trec[i,])
}
summary(cr1)

figNm = file.path(cmc_dir, "output/compare_expression_trec_cmc_cr.pdf")
pdf(figNm, width=4, height=3)
par(mar=c(5,4,1,1))
hist(cr1, main="", xlab="correlation", breaks=50)
dev.off()

# mem_used()
gc()

sessionInfo()
q(save="no")


