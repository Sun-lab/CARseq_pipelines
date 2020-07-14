# Assuming we're on longleaf
rm(list=ls())
cs_dir 		= "/pine/scr/p/l/pllittle/CS_eQTL"
cmc_dir 	= file.path(cs_dir,"s3_Real/CMC")
repo_dir 	= "~/github/CARseq_pipelines"

if(FALSE){ # Assuming mycompH
rm(list=ls()); desk_dir = "C:/Users/plittle/Desktop"
cmc_dir = file.path(desk_dir,"CS_eQTL/s3_Real/CMC")
repo_dir = file.path(desk_dir,"github","CARseq_pipelines")
}

# Source functions
source(file.path(repo_dir,"R","base_functions_plittle.R"))
run_hapQC = function(hap1,hap2,hapN){
	if(FALSE){
		hap1 = rds$hap1; hap2 = rds$hap2; hapN = rds$hapN
	}
	
	# Double check rownames and column names are matched
	check1 = any(colnames(hap1) != colnames(hap2))
	check2 = any(colnames(hap2) != colnames(hapN))
	check3 = any(rownames(hap1) != rownames(hap2))
	check4 = any(rownames(hap2) != rownames(hapN))
	if( check1 || check2 || check3 || check4 ){
		stop("Row or column names aren't matching up!")
	}
	
	tmp_df = smart_df(RNAseq_ID = colnames(hap1),
		subj_hapN = colSums(hapN),
		subj_hapTot = colSums(hap1 + hap2 + hapN))
	tmp_df$prop_hapN = tmp_df$subj_hapN / tmp_df$subj_hapTot
	rownames(tmp_df) = NULL
	# head(tmp_df)
	
	par(mfrow=c(1,2),mar=c(5,4,1,0.5))
	smart_hist(tmp_df$prop_hapN,breaks = 40,main = "",
		xlab = "Proportion of unassigned haplotype counts \n(hapN / (hap1 + hap2 + hapN))")
	plot(tmp_df$subj_hapN/1e3,tmp_df$subj_hapTot/1e6,bty = "n",
		xlab = "# hapN Counts (thousand)",
		ylab = "# hapTotal Counts (millon)",
		col = rgb(0,0,0,0.5),pch = 16)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
	
	tmp_df
	
}
run_TReC_ASReC = function(hap1,hap2,trec){
	if(FALSE){
		hap1 = rds$hap1; hap2 = rds$hap2; trec = rds$trec
	}
	
	tmp_df = smart_df(subj_asrec = colSums(hap1 + hap2),
		subj_trec = colSums(trec))
	tmp_df$prop_asrec = tmp_df$subj_asrec / tmp_df$subj_trec
	tmp_df[1:5,]
	par(mfrow=c(1,2),mar=c(4,4,1,0.5))
	smart_hist(tmp_df$prop_asrec,breaks = 50,main = "",
		xlab = "# ASReC / # TReC")
	plot(tmp_df$subj_trec/1e6,tmp_df$subj_asrec/1e6,
		bty = "n",xlab = "# TReC (million)",
		ylab = "# ASReC (million)",pch = 16,col = rgb(0,0,0,0.5))
	lm_out = lm(subj_asrec ~ -1 + subj_trec,data = tmp_df); lm_out
	tmp_slope = as.numeric(lm_out$coefficients)
	abline(a=0,b=tmp_slope,lwd=2,lty=2,col="red")
	text(min(tmp_df$subj_trec/1e6),
		max(tmp_df$subj_asrec/1e6),adj = 0,
		labels = sprintf("y = %s x",round(tmp_slope,3)))
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
	
}

## Write functions modeled off Vasyl's code
## Want to filter out subjects, then filter out genes, then filter out subject/gene pairs (aka set ASREC set to 0)


# ----------
# CMC-related analyses
# ----------

# Import list of 10 outliers
outlier_fn = file.path(cmc_dir,"CMC_MSSM-Penn-Pitt_DLPFC_mRNA_outlierSamples.txt")
outlier = smart_RT(outlier_fn,header = TRUE,sep = "\t")
outlier

# Import Release3 sampleID key
key_fn = file.path(cmc_dir,"Release3_SampleID_key.csv")
key = smart_RT(key_fn,sep = ",",header = TRUE)
outliers = key[which(key$Individual_ID %in% outlier$Individual.ID),]
outliers

# Import total and haploytype counts
rds_fn = file.path(cmc_dir,"trec_hap1_hap2_hapN.rds")
if( !file.exists(rds_fn) ){
	all_bam_fns = list.files(file.path(cmc_dir,"BAM"))
	length(all_bam_fns)
	bad_samps = c()
	for(ii in seq(length(all_bam_fns))){
		# ii = 1
		if(ii %% 5 == 0) cat(".")
		if(ii %% 1e2 == 0 || ii == length(all_bam_fns)) cat(sprintf("%s\n",ii))
		
		samp = all_bam_fns[ii]
		samp_dir = file.path(cmc_dir,"BAM",samp)
		
		trec_fn = file.path(samp_dir,"gene_level_counts_filterIt_total.txt")
		hap1_fn = file.path(samp_dir,"gene_level_counts_filterIt_hap1.txt")
		hap2_fn = file.path(samp_dir,"gene_level_counts_filterIt_hap2.txt")
		hapN_fn = file.path(samp_dir,"gene_level_counts_filterIt_hapN.txt")
		
		if( !file.exists(hap1_fn) ){
			warning(sprintf("%s files have an issue!\n",samp))
			bad_samps = c(bad_samps,samp)
			next
		}
		
		tmp_trec = smart_RT(trec_fn,sep = "\t",header = TRUE)
		tmp_hap1 = smart_RT(hap1_fn,sep = "\t",header = TRUE)
		tmp_hap2 = smart_RT(hap2_fn,sep = "\t",header = TRUE)
		tmp_hapN = smart_RT(hapN_fn,sep = "\t",header = TRUE)
		
		if(ii == 1){
			trec = matrix(NA,nrow(tmp_trec),length(all_bam_fns))
			rownames(trec) = tmp_trec$genes
			colnames(trec) = all_bam_fns
			hap1 = trec; hap2 = trec; hapN = trec
		}
		
		trec[,samp] = as.numeric(tmp_trec[,3]) # SO_1
		hap1[,samp] = tmp_hap1$counts
		hap2[,samp] = tmp_hap2$counts
		hapN[,samp] = tmp_hapN$counts
		
		rm(tmp_trec,tmp_hap1,tmp_hap2,tmp_hapN)
	}
	if( length(bad_samps) > 0 ){
		trec = trec[,!(colnames(trec) %in% bad_samps)]
		hap1 = hap1[,!(colnames(hap1) %in% bad_samps)]
		hap2 = hap2[,!(colnames(hap2) %in% bad_samps)]
		hapN = hapN[,!(colnames(hapN) %in% bad_samps)]
		print(bad_samps)
	}
	# MSSM_RNA_PFC_331 seems to be a bad quality sample b/c no heterozygous SNPs!
	saveRDS(list(outliers = outliers,trec = trec,
		hap1 = hap1,hap2 = hap2,hapN = hapN),rds_fn)
}

# Analysis/Plotting
rds = readRDS(rds_fn)
names(rds)
str(rds)
rds$outliers
rds$trec[1:10,1:5]

## The ten outliers labeled by CMC's analysis
rds$outliers[,c("Individual_ID","RNAseq.Sample_RNA_ID")]

## Remove ten outliers
int_subj = intersect(colnames(rds$trec),key$RNAseq.Sample_RNA_ID)
int_subj = int_subj[!(int_subj %in% rds$outliers$RNAseq.Sample_RNA_ID)]
length(int_subj)
rds$trec = rds$trec[,int_subj]
rds$hap1 = rds$hap1[,int_subj]
rds$hap2 = rds$hap2[,int_subj]
rds$hapN = rds$hapN[,int_subj]

## Isolating four potential additional CMC outliers
png(file.path(repo_dir,"results","hapQC_CMC.png"),
	units = 'px',height = 1500,width = 3000,res = 250,
	type = 'cairo',pointsize = 20)
out_hapQC = run_hapQC(hap1 = rds$hap1,hap2 = rds$hap2,hapN = rds$hapN)
dev.off()
smart_WT(out_hapQC,file.path(repo_dir,"results","hapQC_CMC.tsv"),sep = "\t")

outlier2 = out_hapQC[which(out_hapQC$subj_hapN/1e3 > 8
	& out_hapQC$subj_hapTot/1e6 < 1.5),]; outlier2

## Collect outlier samples
outlie_df = smart_df(RNAseq_ID = rds$outliers$RNAseq.Sample_RNA_ID,
	outlier_reason = "CMC")
outlie_df = rbind(outlie_df,smart_df(RNAseq_ID = rownames(outlier2),
	outlier_reason = "hapQC"))
# outlie_df = rbind(outlie_df,smart_df(RNAseq_ID = "MSSM_RNA_PFC_331",outlier_reason = "no_hetSNPs"))
outlie_df = outlie_df[order(outlie_df$RNAseq_ID),]
rownames(outlie_df) = seq(nrow(outlie_df))
outlie_df

## It seems MSSM_RNA_PFC_260 is an outlier from both CMC and our hapQC plot!

## TReC and ASReC
png(file.path(repo_dir,"results","TReC_ASReC_CMC.png"),
	units = 'px',height = 1500,width = 3000,res = 250,
	type = 'cairo',pointsize = 20)
run_TReC_ASReC(hap1 = rds$hap1,hap2 = rds$hap2,trec = rds$trec)
dev.off()

###


