# Phasing and Imputation from bed/bim/fam file

# ----------
# Notes
# ----------
## work_dir = working directory to place files
## input_fn: E.g. "CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_QCed"
## plink_dir = plink directory
## shapeit_dir = shapeit directory
## impute2_dir = impute2 directory
## src_dir = CARseq_pipelines/src directory


# ----------
# Source base functions and Libraries
# ----------
# source("base_functions_plittle.R")
library(Rcpp); library(RcppArmadillo)
library(data.table)


# ----------
# Phasing/Imputation Functions
# ----------
setdirs_SPI = function(work_dir){
	# Make/Set directories
	chr_dir 						= file.path(work_dir,"chr"); smart_mkdir(chr_dir)
	precheck_dir 				= file.path(chr_dir,"precheck_log"); smart_mkdir(precheck_dir)
	impute_dir 					= file.path(work_dir,"impute"); smart_mkdir(impute_dir)
	thsdGP3_dir 				= file.path(impute_dir,"1000GP_Phase3"); smart_mkdir(thsdGP3_dir)
	chr_flip_dir 				= file.path(work_dir,"chr_flipped"); smart_mkdir(chr_flip_dir)
	precheck_flip_dir 	= file.path(chr_flip_dir,"precheck_log"); smart_mkdir(precheck_flip_dir)
	phased_dir 					= file.path(work_dir,"phased"); smart_mkdir(phased_dir)
	res_dir 						= file.path(work_dir,"shapeit"); smart_mkdir(res_dir)
	imputation_dir 			= file.path(work_dir,"imputation"); smart_mkdir(imputation_dir)
	imputed_dir 				= file.path(work_dir,"imputed"); smart_mkdir(imputed_dir)
	data_genotype_dir 	= file.path(work_dir,"data_genotype"); smart_mkdir(data_genotype_dir)
	data_snp_dir				= file.path(work_dir,"data_snp"); smart_mkdir(data_snp_dir)
	dg38_dir 						= file.path(work_dir,"dg38"); smart_mkdir(dg38_dir)
	
	list(spi_dir = work_dir,chr_dir = chr_dir,precheck_dir = precheck_dir,
		impute_dir = impute_dir,thsdGP3_dir = thsdGP3_dir,
		chr_flip_dir = chr_flip_dir,
		precheck_flip_dir = precheck_flip_dir,
		phased_dir = phased_dir,res_dir = res_dir,
		imputation_dir = imputation_dir,imputed_dir = imputed_dir,
		data_genotype_dir = data_genotype_dir,data_snp_dir = data_snp_dir,
		dg38_dir = dg38_dir)
}
step1a_bed2ped = function(plink_dir,work_dir,input_fn,excludesnps=FALSE){
	if(FALSE){
		plink_dir = plink_dir
		work_dir = cmc_dir
		input_fn = "CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_QCed"
	}
	
	# Convert bim/bed/fam files to ped/map and by chromosome
	my_dirs = setdirs_SPI(work_dir = work_dir)
	
	system(smart_sprintf("cd %s; echo '.' > excludesnps.txt",work_dir))
	
	cmd = smart_sprintf("
		plink_dir=%s; 
		cd %s; 
		$plink_dir/plink --bfile %s \
			--exclude excludesnps.txt \
			--recode A-transpose \
			--out geno.call; ",
		plink_dir,work_dir,input_fn)
	system(cmd)
	
	# smart_mkdir(file.path(work_dir,"chr"))
	if( excludesnps ){
		cmd = smart_sprintf("
			plink_dir=%s; 
			cd %s; 
			for chr in `seq 1 22`; 
				do 
				$plink_dir/plink --bfile %s \
					--chr $chr \
					--exclude excludesnps.txt \
					--recode \
					--out chr/geno.chr$chr; 
				$plink_dir/plink --file chr/geno.chr$chr \
					--recode A-transpose \
					--out chr/geno.call.chr$chr; 
			done",
			plink_dir,work_dir,input_fn)
	} else {
		cmd = smart_sprintf("
			plink_dir=%s; 
			cd %s; 
			for chr in `seq 1 22`; 
				do 
				$plink_dir/plink --bfile %s \
					--chr $chr \
					--recode \
					--out chr/geno.chr$chr; 
				$plink_dir/plink --file chr/geno.chr$chr \
					--recode A-transpose \
					--out chr/geno.call.chr$chr; 
			done",
			plink_dir,work_dir,input_fn)
	}
	
	cat(paste0("\n#########\n",cmd,"\n#########\n"))
	system(cmd)
}
step1b_shapeit = function(shapeit_dir,work_dir,chr,run_shapeit=FALSE,threads=1){
	# IMPUTE files: http://mathgen.stats.ox.ac.uk/impute/
	my_dirs = setdirs_SPI(work_dir = work_dir)
	precheck_dir 	= my_dirs$precheck_dir
	impute_dir 		= my_dirs$impute_dir
	thsdGP3_dir 	= my_dirs$thsdGP3_dir
	
	log_fn 				= file.path(precheck_dir,paste0("gwas.alignments_chr",chr))
	input_ped 		= file.path("chr",paste0("geno.chr",chr))
	input_map 		= file.path(thsdGP3_dir,paste0("genetic_map_chr",chr,"_combined_b37.txt"))
	input_ref 		= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".hap.gz"))
	input_legend 	= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".legend.gz"))
	input_sample 	= file.path(thsdGP3_dir,"1000GP_Phase3.sample")
	
	# Download reference files
	if( !file.exists(input_map) ){
		map_link 	= paste0("http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/genetic_map_chr",chr,"_combined_b37.txt")
		cmd 			= smart_sprintf("cd %s; wget %s",thsdGP3_dir,map_link)
		system(cmd)
	}
	if( !file.exists(input_ref) ){
		map_link 	= paste0("http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3_chr",chr,".hap.gz")
		cmd 			= smart_sprintf("cd %s; wget %s",thsdGP3_dir,map_link)
		system(cmd)
	}
	if( !file.exists(input_legend) ){
		map_link 	= paste0("http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3_chr",chr,".legend.gz")
		cmd 			= smart_sprintf("cd %s; wget %s",thsdGP3_dir,map_link)
		system(cmd)
	}
	if( !file.exists(input_sample) ){
		map_link 	= "http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3.sample"
		cmd 			= smart_sprintf("cd %s; wget %s",thsdGP3_dir,map_link)
		system(cmd)
	}
	
	# Refer to this to interpret output https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#reference
	# Basically this runs a check on inputted snps against the 1000Genomes reference panel of haplotypes
	# If there are matches, good, if not the output files will contains snps to be excluded and that have strand issues
	cmd = smart_sprintf(
		"cd %s; 
		%s/shapeit -check \
			--input-ped %s \
			--input-map %s \
			--input-ref %s \
			%s \
			%s \
			--thread %s \
			--output-log %s \
			> /dev/null",
		work_dir,shapeit_dir,input_ped,
		input_map,input_ref,input_legend,
		input_sample,threads,log_fn)
	if( run_shapeit ){
		cat(paste0("\n#########\n",cmd,"\n#########\n"))
		system(cmd)
	}
}
step2_flipcheck = function(plink_dir,shapeit_dir,work_dir){
	if(FALSE){
		plink_dir = "/nas/longleaf/apps/plink/1.90b3/bin"
		shapeit_dir = "/nas/longleaf/apps/shapeit/2.837/bin"
		work_dir = spi_dir
	}
	
	my_dirs = setdirs_SPI(work_dir = work_dir)
	chr_dir				= my_dirs$chr_dir
	precheck_dir 	= my_dirs$precheck_dir
	# impute_dir 		= my_dirs$impute_dir
	thsdGP3_dir 	= my_dirs$thsdGP3_dir
	
	input_sample 	= file.path(thsdGP3_dir,"1000GP_Phase3.sample")
	checki = matrix(NA, nrow = 22, ncol = 2)
	rownames(checki) = paste0("chr",1:22)
	
	for(chr in seq(22)){
		# chr = 1
		cat(paste0("Checking chr",chr,"...\n"))
		
		# flip the snps of interest
		con = file(file.path(precheck_dir,paste0("gwas.alignments_chr",chr,".snp.strand")))
		lns = readLines(con)
		lns = lns[substr(lns, 1, 6) == "Strand"] # extract SNP with strand issue 
		strand = matrix(unlist(strsplit(lns, split = "\t")), nrow = 11)[4,]
		close(con)
		tmp_strand_fn = file.path(precheck_dir,"tmp_strand.txt")
		write.table(x = unique(strand),
					file = tmp_strand_fn,
					row.names = FALSE,
					col.names = FALSE,
					quote = FALSE) # write misstrand SNP to a tmp file
		
		## copy chr map/ped files to strand_chr
		geno_chr = file.path(chr_dir,paste0("geno.chr",chr))
		tmp_geno_chr = file.path(precheck_dir,"tmp")
		
		old_map = paste0(geno_chr,".map")
		tmp_map = paste0(tmp_geno_chr,".map")
		smart_remove(tmp_map)
		file.copy(from = old_map,to = tmp_map)
		
		old_ped = paste0(geno_chr,".ped")
		tmp_ped = paste0(tmp_geno_chr,".ped")
		smart_remove(tmp_ped)
		file.copy(from = old_ped,to = tmp_ped)
		
		## run plink to flip subset of strand issue snps
		cmd = smart_sprintf("
			cd %s; 
			%s/plink --file %s --flip %s --recode \
			> /dev/null",
			precheck_dir,plink_dir,
			tmp_geno_chr,tmp_strand_fn)
		cat(paste0("\n#########\n",cmd,"\n#########\n"))
		system(cmd)
		
		## check that the plink really flipped those snps
		## we see that plink does flip those strands
		con = file(tmp_ped)
		tmp = readLines(con)
		close(con)

		con = file(file.path(precheck_dir,"plink.ped"))
		tmp0 = readLines(con)
		close(con)
		message(!all(unlist(tmp) == unlist(tmp0)))
		
		## recheck using shapeit on flipped data
		input_map 		= file.path(thsdGP3_dir,paste0("genetic_map_chr",chr,"_combined_b37.txt"))
		input_ref 		= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".hap.gz"))
		input_legend 	= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".legend.gz"))
		cmd = smart_sprintf("
			cd %s; 
			%s/shapeit -check \
				--input-ped %s \
				--input-map %s \
				--input-ref %s %s %s \
				--output-log tmp \
				> /dev/null",
			precheck_dir,shapeit_dir,
			file.path(precheck_dir,"plink"),
			input_map,input_ref,input_legend,
			input_sample)
		cat(paste0("\n#########\n",cmd,"\n#########\n"))
		system(cmd)
		
		con = file(file.path(precheck_dir,"tmp.snp.strand"))
		lns2 = readLines(con)
		lns2 = lns2[substr(lns2, 1, 6) == "Strand"]
		strand2 = matrix(unlist(strsplit(lns2, split = "\t")), nrow = 11)[4,]
		close(con)

		# count how many snps are in the one exclusion set, but not in the other
		checki[chr, 1] = length(setdiff(unique(strand), unique(strand2)))
		checki[chr, 2] = length(setdiff(unique(strand2), unique(strand)))
		
		# Clean up
		system(smart_sprintf("cd %s; rm tmp*; rm plink*; ",precheck_dir))
	}

	print(checki)
	smart_WT(checki,file.path(precheck_dir,"step2_check_chr.log"),sep="\t")
	
	return(NULL)
}
step3_chr_summary = function(work_dir){
	if(FALSE){
		work_dir = cmc_dir
	}
	
	my_dirs = setdirs_SPI(work_dir = work_dir)
	chr_dir 			= my_dirs$chr_dir
	precheck_dir 	= my_dirs$precheck_dir
	
	chrs = seq(22)
	summ = matrix(NA,nrow = 22,ncol = 4)
	colnames(summ) = c("total","not.comp","strand","missing")
	rownames(summ) = paste0("chr",1:22)
	
	for(chr in chrs){
		# chr = 1
		cmd = smart_sprintf("wc -l %s/geno.chr%s.map | cut -d ' ' -f1",chr_dir,chr)
		summ[chr, 1] = as.numeric(system(cmd,intern=TRUE))
		gwas_align_fn = file.path(precheck_dir,paste0("gwas.alignments_chr",chr,".snp.strand"))
		dat = read.table(gwas_align_fn,header = FALSE,as.is = TRUE,sep = "\t",skip = 1)
		names(dat)[-1] = strsplit(readLines(gwas_align_fn,n=1),"\t")[[1]]
		summ[chr,2] = nrow(dat)
		summ[chr,3] = length(which(dat$V1 == "Strand"))
		summ[chr,4] = length(which(dat$V1 == "Missing"))
	}
	
	print(summ)
	smart_WT(summ,file.path(precheck_dir,"step3_summ.log"),sep="\t")
	
	return(NULL)
}
step4_flip_snp = function(plink_dir,shapeit_dir,work_dir){
	if(FALSE){
		work_dir = cmc_dir
	}
	
	my_dirs = setdirs_SPI(work_dir = work_dir)
	chr_dir						= my_dirs$chr_dir
	precheck_dir 			= my_dirs$precheck_dir
	impute_dir 				= my_dirs$impute_dir
	thsdGP3_dir 			= my_dirs$thsdGP3_dir
	chr_flip_dir			= my_dirs$chr_flip_dir
	precheck_flip_dir	= my_dirs$precheck_flip_dir
	
	input_sample 		= file.path(thsdGP3_dir,"1000GP_Phase3.sample")
	for(chr in seq(22)){
		# chr = 1
		cat(paste0("Flip snps for chr",chr,"...\n"))
		
		# flip the snps of interest
		gwas_align_fn = file.path(precheck_dir,paste0("gwas.alignments_chr",chr,".snp.strand"))
		dat = smart_RT(gwas_align_fn,header = FALSE,sep = "\t",skip = 1)
		names(dat)[-1] = strsplit(readLines(gwas_align_fn,n=1),"\t")[[1]]
		dat = dat[which(dat$type == "Strand"),]
		tmp_strand_fn = file.path(precheck_dir,"tmp_strand.txt")
		write.table(x = unique(dat$main_id),file = tmp_strand_fn,
			row.names = FALSE,col.names = FALSE,quote = FALSE)
		
		geno_chr = file.path(chr_dir,paste0("geno.chr",chr))
		tmp_geno_chr = file.path(precheck_dir,"tmp")
		
		old_map = paste0(geno_chr,".map")
		tmp_map = paste0(tmp_geno_chr,".map")
		smart_remove(tmp_map)
		file.copy(from = old_map,to = tmp_map)
		
		old_ped = paste0(geno_chr,".ped")
		tmp_ped = paste0(tmp_geno_chr,".ped")
		smart_remove(tmp_ped)
		file.copy(from = old_ped,to = tmp_ped)
		
		flipped_fn = file.path(chr_flip_dir,paste0("geno.chr",chr,"_flip"))
		cmd = smart_sprintf("
			cd %s; 
			%s/plink --file tmp \
				--flip tmp_strand.txt \
				--recode --out %s \
				> /dev/null",
			precheck_dir,plink_dir,flipped_fn)
		cat(paste0("\n#########\n",cmd,"\n#########\n"))
		system(cmd)
		
		# recheck using shapeit on flipped data
		input_map 		= file.path(thsdGP3_dir,paste0("genetic_map_chr",chr,"_combined_b37.txt"))
		input_ref 		= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".hap.gz"))
		input_legend 	= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".legend.gz"))
		out_fn 				= file.path(precheck_flip_dir,paste0("geno.chr",chr,"_flip_check"))
		
		cmd = smart_sprintf("
			%s/shapeit -check \
				--input-ped %s \
				--input-map %s \
				--input-ref %s \
				%s \
				%s \
				--output-log %s > /dev/null",
			shapeit_dir,flipped_fn,input_map,
			input_ref,input_legend,input_sample,
			out_fn)
		cat(paste0("\n#########\n",cmd,"\n#########\n"))
		system(cmd)
		
		# Clean up
		system(smart_sprintf("cd %s; rm tmp*;",precheck_dir))
		
	}

	return(NULL)
}
step5_shapeit = function(shapeit_dir,work_dir,chr,threads=1){
	if(FALSE){
		work_dir = cmc_dir
		chr = 1
		threads = 1
	}
	
	my_dirs = setdirs_SPI(work_dir = work_dir)
	impute_dir 					= my_dirs$impute_dir
	thsdGP3_dir 				= my_dirs$thsdGP3_dir
	chr_flip_dir 				= my_dirs$chr_flip_dir
	precheck_flip_dir		= my_dirs$precheck_flip_dir
	phased_dir 					= my_dirs$phased_dir
	res_dir 						= my_dirs$res_dir
	setwd(res_dir)
	
	flip_ped 				= file.path(chr_flip_dir,paste0("geno.chr",chr,"_flip"))
	input_map 			= file.path(thsdGP3_dir,paste0("genetic_map_chr",chr,"_combined_b37.txt"))
	input_ref 			= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".hap.gz"))
	input_legend 		= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".legend.gz"))
	input_sample 		= file.path(thsdGP3_dir,"1000GP_Phase3.sample")
	exclu_snp_fn 		= file.path(precheck_flip_dir,paste0("geno.chr",chr,"_flip_check.snp.strand.exclude"))
	phased_chr_fn 	= file.path(phased_dir,paste0("phasing_chr",chr))
	stdout_fn 			= file.path(res_dir,paste0("shapeit_chr",chr,".o"))
	stderr_fn 			= file.path(res_dir,paste0("shapeit_chr",chr,".e"))
	
	cmd = smart_sprintf("%s/shapeit \
		--input-ped %s \
		--input-map %s \
		--input-ref %s %s %s \
		--exclude-snp %s \
		-O %s \
		--effective-size 20000 \
		--seed 1234567 \
		--thread %s \
		> %s \
		2> %s",
		shapeit_dir,flip_ped,input_map,
		input_ref,input_legend,input_sample,
		exclu_snp_fn,phased_chr_fn,threads,
		stdout_fn,stderr_fn)
	
	cat(paste0("\n#########\n",cmd,"\n#########\n"))
	system(cmd)
	
	return(NULL)
}
step6_imputation = function(impute2_dir,work_dir,chr,mem=25000){
	if(FALSE){
		# impute2_dir = impute2_dir
		work_dir = cmc_dir
		chr = 1
		mem = 4000
	}
	
	# Input: work_dir/phased/{*.haps,*.sample}
	# Output: 
	
	# ensure that when we write positions to a file R doesn't switch them to scientific
	options("scipen"=999,"digits"=4)
	
	my_dirs = setdirs_SPI(work_dir = work_dir)
	impute_dir 				= my_dirs$impute_dir
	thsdGP3_dir 			= my_dirs$thsdGP3_dir
	phased_dir 				= my_dirs$phased_dir
	imputation_dir 		= my_dirs$imputation_dir
	imputed_dir 			= my_dirs$imputed_dir
	
	# since we prephased using shapeit we need to use -use_prephased_g and -known_haps_g
	input_map 		= file.path(thsdGP3_dir,paste0("genetic_map_chr",chr,"_combined_b37.txt"))
	input_ref 		= file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".hap.gz"))
	
	#######################
	##### Editting block ##
	#######################
	get_block = function(str,split=" ",block=2){
		unlist(strsplit(str,split=split))[block]
	}
	
	# using the appropriate positions from both data and 1000G we find the range in which we will do imputation
	phased_chr_fn = file.path(phased_dir,paste0("phasing_chr",chr,".haps"))
	cat(paste0("Reading in ",phased_chr_fn,"...\n"))
	dati = read.table(phased_chr_fn,as.is = TRUE)
	rng = range(dati[,3])
	
	input_legend = file.path(thsdGP3_dir,paste0("1000GP_Phase3_chr",chr,".legend.gz"))
	cat(paste0("Reading in ",input_legend,"...\n"))
	con = gzfile(input_legend)
	refi = readLines(con)
	close(con)
	
	cat(paste0("Making rng2...\n"))
	rng2 = range(as.numeric(sapply(refi[-1],get_block)))
	rng[1] = min(rng[1],rng2[1])
	rng[2] = max(rng[2],rng2[2])
	message(paste(rng2-rng,collapse=" "))
	
	mb = 1e6
	cat(paste0("Making indj...\n"))
	indj = sprintf("%se6",c(seq(floor(rng[1]/1e6),floor(rng2[2]/mb),by=5),ceiling(rng[2]/mb)))
	write.table(data.frame(chr,paste(indj,collapse=";")),
		file.path(imputation_dir,"indicies.txt"),
		append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
	
	# here we will run by the 5MB increment from the rng[1] until we exhaust the chromosome
	cat(paste0("Submitting jobs...\n"))
	for(jj in 2:length(indj)){
		# jj = length(indj)
		# jj = 2
		output_fn = file.path(imputed_dir,paste0("phased_imputed_chr",chr,"_",indj[jj-1],"_",indj[jj]))
		impute_chr_prefix_fn = file.path(imputation_dir,paste0("impute_chr",chr,"_",indj[jj-1],"_",indj[jj]))
		cmd = smart_sprintf("
			%s/impute2 -use_prephased_g \
				-m %s -h %s -l %s \
				-known_haps_g %s \
				-align_by_maf_g -Ne 20000 -seed 12345 \
				-phase -o %s -int %s %s \
				> %s 2> %s",
			impute2_dir,input_map,input_ref,
			input_legend,phased_chr_fn,output_fn,
			indj[jj-1],indj[jj],
			paste0(impute_chr_prefix_fn,".o"),
			paste0(impute_chr_prefix_fn,".e"))
		
		# Show command, inputs, options, outputs
		cat(paste0("\n###############\n",
			gsub("-","\\\\\n\t-",cmd),
			"\n###############\n\n"))
		
		# Write job command
		job_cmd = smart_sprintf("
			curr_dir=%s; 
			cd $curr_dir; 
			chr=%s; 
			start=%s; 
			end=%s; 
			mem=%s; 
			sbatch -n 1 -c 1 \
				--job-name=impute2_chr${chr}_${start}_${end} \
				--mem=${mem} \
				-o chr${chr}_${start}_${end}.log \
				-t 0-5 \
				--wrap='%s'",
			imputation_dir,chr,indj[jj-1],indj[jj],mem,cmd)
		
		# Run job
		system(job_cmd)
	}

}
get_anno_geno = function(work_dir,sample_id,chr){
	if(FALSE){
		work_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/CMC"
		sample_id = "0_MSSM_13"  # usually the format is "TCGA-B0-5099-01A"
		chr = 1
		
	}
	
	# Load the SNP 6 csv annotation file
	snp6_anno_fn = file.path(work_dir,"GenomeWideSNP_6.na35.annot.csv")
	if( !file.exists(snp6_anno_fn) ){
		# Source: http://www.affymetrix.com/support/technical/byproduct.affx?product=genomewidesnp_6 
		src_link = "http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip"
		down_fn = strsplit(src_link,"/")[[1]]
		down_fn = down_fn[length(down_fn)]
		cmd = sprintf("cd %s; wget %s; unzip %s",work_dir,src_link,down_fn)
		system(cmd)
	}
	
	cat(paste0("Reading in ",snp6_anno_fn,"...\n"))
	num_lines = as.numeric(system(sprintf("wc -l %s | cut -d ' ' -f1",snp6_anno_fn),intern = TRUE)) - 1
	snpcsv = smart_RT(snp6_anno_fn,sep = ',',header = TRUE,as.is = TRUE,nrows = num_lines)
	snpcsv = snpcsv[,c("Probe.Set.ID","dbSNP.RS.ID","Chromosome","Physical.Position","Allele.A","Allele.B")]
	# print(head(snpcsv))
	
	# Import Genotype calls for a chromosome and subset a sample
	geno_chr_fn = file.path(work_dir,"chr",paste0("geno.call.chr",chr,".traw"))
	cat(paste0("Reading in ",geno_chr_fn,"...\n"))
	num_lines = as.numeric(system(sprintf("wc -l %s | cut -d ' ' -f1",geno_chr_fn),intern = TRUE))
	geno_chr = smart_RT(geno_chr_fn,sep = "\t",header = TRUE,as.is = TRUE, check.names = FALSE,nrows = num_lines)
	# print(dim(geno_chr))
	# print(geno_chr[1:5,1:10])
	
	# 2.1. Get the column corresponding to sample_id
	sample_geno_chr = geno_chr[,c("CHR","SNP","POS","COUNTED","ALT",sample_id)]
	sample_geno_chr = name_change(sample_geno_chr,"POS","Physical.Position")
	sample_geno_chr = name_change(sample_geno_chr,"SNP","dbSNP.RS.ID")
	sample_geno_chr = name_change(sample_geno_chr,"CHR","Chromosome")
	sample_geno_chr$Chromosome = as.character(sample_geno_chr$Chromosome)
	# print(head(sample_geno_chr))
	
	# 2.2 Merge birdseed calls with SNP 6 reference
	snpcsv = snpcsv[which(snpcsv$Chromosome %in% unique(sample_geno_chr$Chromosome)),]
	snpcsv$Physical.Position = as.numeric(snpcsv$Physical.Position)
	geno_pre_phase = smart_merge(snpcsv,sample_geno_chr)
	# print(nrow(geno_pre_phase))
	geno_pre_phase = geno_pre_phase[which(geno_pre_phase$Chromosome != "---"),]
	# print(nrow(geno_pre_phase))
	# print(head(geno_pre_phase))
	
	# 2.3 Sort
	geno_pre_phase = geno_pre_phase[order(geno_pre_phase$Chromosome,
		geno_pre_phase$Physical.Position),]
	# print(head(geno_pre_phase))
	
	geno_pre_phase
}
step7_check_inputGeno_outputPhase = function(work_dir,sample_id,chr,geno_pre_phase){
	if(FALSE){
		work_dir = spi_dir
		sample_id = "0_MSSM_BP_10"
		chr = 1
		geno_pre_phase = geno_pre_phase
	}
	
	my_dirs = setdirs_SPI(work_dir = work_dir)
	phased_dir 				= my_dirs$phased_dir
	
	# Checking the pipeline input genotype calls vs. the phased output
	
	# 3. Load the genotype after phasing using shapeit
	sample_table_fn = file.path(phased_dir,paste0("phasing_chr",chr,".sample"))
	cat(paste0("Reading in ",sample_table_fn,"...\n"))
	sample_table = smart_RT(sample_table_fn,sep = " ",header = TRUE,as.is = TRUE)
	print(head(sample_table))
	sample_id_1 = strsplit(sample_id,"_")[[1]][1]
	sample_id_2 = paste(strsplit(sample_id,"_")[[1]][-1],collapse="_")
	sample_index = which(sample_table$ID_1 == sample_id_1
		& sample_table$ID_2 == sample_id_2)
	
	haps_table_fn = file.path(phased_dir,paste0("phasing_chr",chr,".haps"))
	num_lines = as.numeric(system(sprintf("wc -l %s | cut -d ' ' -f1",haps_table_fn),intern=TRUE))
	haps_table = smart_RT(haps_table_fn,sep = " ",header = FALSE,nrows = num_lines)
	hap_cols = 4 + c(0,1) + sample_index
	haps_table2 = haps_table[,c(1:5,hap_cols)]
	# head(haps_table2)
	names(haps_table2) = c("Chromosome","dbSNP.RS.ID","Physical.Position","Allele.A","Allele.B","hap1","hap2")
	haps_table2$Chromosome = as.character(haps_table2$Chromosome)
	print(head(haps_table2))
	
	# 4. Merge the two lists
	genotype_merged = merge(x = geno_pre_phase,y = haps_table2,
		by = c("dbSNP.RS.ID","Chromosome","Physical.Position"))
	genotype_merged = genotype_merged[order(genotype_merged$Chromosome,
		genotype_merged$Physical.Position),]
	# nrow(genotype_merged)
	print(head(genotype_merged))
	print(smart_table(genotype_merged[,sample_id]))
	print(smart_table(genotype_merged$Allele.A.x, genotype_merged$Allele.A.y))
	print(smart_table(genotype_merged$Allele.B.x, genotype_merged$Allele.B.y))
	print(smart_table(genotype_merged[,sample_id],genotype_merged$hap1+genotype_merged$hap2))

	#write.table(genotype_merged, file="step7_genotype_calls_vs_shapeit.txt", row.names=FALSE, quote=FALSE, sep="\t")

}
step8_check_inputGeno_outputFlip = function(work_dir,sample_id,chr,geno_pre_phase){
	if(FALSE){
		work_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/CMC"
		sample_id = "0_MSSM_BP_10"  # usually the format is "TCGA-B0-5099-01A"
		chr = 1
		geno_pre_phase = geno_pre_phase
	}

	# Checking the pipeline input genotype calls vs. the flipped output
	my_dirs = setdirs_SPI(work_dir = work_dir)
	chr_flip_dir = my_dirs$chr_flip_dir
	
	# 3. Load the genotype after phasing using shapeit
	sample_table_fn = file.path(chr_flip_dir,paste0("geno.chr",chr,"_flip.ped"))
	sample_table = system(sprintf("cut -d ' ' -f1-2 %s | sed 's/ /_/g'",sample_table_fn),intern=TRUE)
	sample_index = which(sample_table == sample_id)
	rm(sample_table)
	
	ped_fn = file.path(chr_flip_dir,paste0("geno.chr",chr,"_flip.ped"))
	pedfile = readLines(ped_fn)[sample_index]
	pedfile = strsplit(pedfile," ")[[1]]
	geno1 = pedfile[seq(7,length(pedfile)-1,by=2)]
	geno2 = pedfile[seq(8,length(pedfile),by=2)]
	
	map_fn = file.path(chr_flip_dir,paste0("geno.chr",chr,"_flip.map"))
	num_lines = as.numeric(system(sprintf("wc -l %s | cut -d ' ' -f1",map_fn),intern = TRUE))
	cat(paste0("Reading in ",map_fn,"...\n"))
	mapfile = smart_RT(map_fn,sep = "\t",header = FALSE,nrows = num_lines)
	mapfile = mapfile[,c(1,2,4)]
	
	haps_table = smart_df(mapfile,geno1,geno2)
	names(haps_table) = c("Chromosome", "dbSNP.RS.ID", "Physical.Position", "genotype1", "genotype2")
	# print(head(haps_table))
	# write.table(haps_table, file="genotype_flipped.txt", row.names=FALSE, quote=FALSE, sep="\t")
	
	# 4. Merge the two lists
	genotype_merged = merge(x = geno_pre_phase,y = haps_table,
		by = c("dbSNP.RS.ID","Chromosome","Physical.Position"))
	genotype_merged = genotype_merged[order(genotype_merged$Chromosome,
		genotype_merged$Physical.Position),]
	# nrow(genotype_merged)
	# head(genotype_merged)
	print(smart_table(genotype_merged[,sample_id]))
	print(smart_table(genotype_merged[,sample_id],genotype_merged$genotype1 == genotype_merged$genotype2))

	# write.table(genotype_merged, file="genotype_calls_vs_flipped.txt", row.names=FALSE, quote=FALSE, sep="\t")
}
step9_check_inputGeno_outputUnflip = function(work_dir,sample_id,chr,geno_pre_phase){
	if(FALSE){
		work_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/CMC"
		sample_id = "0_MSSM_13"  # usually the format is "TCGA-B0-5099-01A"
		chr = 1
		geno_pre_phase = geno_pre_phase
	}
	
	# Checking the pipeline input genotype_calls.txt vs. the unflipped input
	my_dirs = setdirs_SPI(work_dir = work_dir)
	chr_dir	= my_dirs$chr_dir
	
	# 3. Load the genotype after phasing using shapeit
	sample_table_fn = file.path(chr_dir,paste0("geno.chr",chr,".ped"))
	sample_table = system(sprintf("cut -d ' ' -f1-2 %s | sed 's/ /_/g'",sample_table_fn),intern = TRUE)
	sample_index = which(sample_table == sample_id)
	rm(sample_table)
	
	ped_fn = file.path(chr_dir,paste0("geno.chr",chr,".ped"))
	pedfile = readLines(ped_fn)[sample_index]
	pedfile = strsplit(pedfile," ")[[1]]
	geno1 = pedfile[seq(7,length(pedfile)-1,by=2)]
	geno2 = pedfile[seq(8,length(pedfile),by=2)]
	
	map_fn = file.path(chr_dir,paste0("geno.chr",chr,".map"))
	num_lines = as.numeric(system(sprintf("wc -l %s | cut -d ' ' -f1",map_fn),intern = TRUE))
	cat(paste0("Reading in ",map_fn,"...\n"))
	mapfile = smart_RT(map_fn,sep = "\t",header = FALSE,nrows = num_lines)
	mapfile = mapfile[,c(1,2,4)]
	
	haps_table = smart_df(mapfile,geno1,geno2)
	names(haps_table) = c("Chromosome", "dbSNP.RS.ID", "Physical.Position", "genotype1", "genotype2")
	# head(haps_table)
	# write.table(haps_table, file="genotype_unflipped.txt", row.names=FALSE, quote=FALSE, sep="\t")
	
	# 4. Merge the two lists
	genotype_merged = merge(x = geno_pre_phase,y = haps_table,
		by = c("dbSNP.RS.ID","Chromosome","Physical.Position"))
	genotype_merged = genotype_merged[order(genotype_merged$Chromosome, genotype_merged$Physical.Position),]
	# nrow(genotype_merged)
	# head(genotype_merged)
	print(smart_table(genotype_merged[,sample_id]))
	print(smart_table(genotype_merged[,sample_id],genotype_merged$genotype1 == genotype_merged$genotype2))

	# write.table(genotype_merged, file="genotype_calls_vs_unflipped.txt", row.names=FALSE, quote=FALSE, sep="\t")
	
}
step10_most_like_geno = function(work_dir,chr,geno_thres=0.8,src_dir){
	if(FALSE){
		work_dir = cmc_dir
		chr = 1; geno_thres = 0.8
		src_dir = "~/github/CARseq_pipelines/src"
	}
	
	# ensure that when we write positions to a file R doesn't switch them to scientific
	options("scipen"=999,"digits"=4)
	my_dirs = setdirs_SPI(work_dir = work_dir)
	phased_dir 	 			= my_dirs$phased_dir
	imputed_dir				= my_dirs$imputed_dir
	data_genotype_dir	= my_dirs$data_genotype_dir
	
	# Load cpp function
	Rcpp::sourceCpp(file.path(src_dir,"test_asSeq.cpp"),showOutput = TRUE)
	
	# get list of samples
	phased_chr_fn = file.path(phased_dir,paste0("phasing_chr",chr,".sample"))
	samples = smart_RT(phased_chr_fn,sep = " ",header = TRUE)
	samples = samples[-1,]
	# samples = paste0(samples$ID_1,"_",samples$ID_2)
	samples = samples$ID_2
	# define snps (I'll exclude indels with matching to the base
	bases = c("A","C","G","T")
	
	# get all the appropriate chunks for the given chromosome
	fls = list.files(path = imputed_dir,pattern = paste0("chr",chr,"_"))
	fls = fls[-grep("warnings",fls)]
	fls = fls[-grep("summary",fls)]
	fls = fls[-grep("info",fls)]
	if( length(grep(".txt",fls)) > 0 ) fls = fls[-grep(".txt",fls)]
	hps = fls[grep("_haps",fls)]
	alp = fls[grep("_allele_probs",fls)]
	fls = setdiff(fls,union(hps,alp))
	
	# Ensure genomic order
	fls = smart_df(label = fls)
	fls$start = as.numeric(sapply(fls$label,function(xx) strsplit(xx,"_")[[1]][4],USE.NAMES=FALSE))
	fls$end = as.numeric(sapply(fls$label,function(xx) strsplit(xx,"_")[[1]][5],USE.NAMES=FALSE))
	fls = fls[order(fls$start),]
	fls = fls$label
	
	# append = FALSE
	off = 5; info_thres = 0.3; min_af = 0.05; # geno_thres = 0.8
	sample_count = length(samples)
	out_fn = file.path(data_genotype_dir,paste0("chr",chr,".txt"))
	smart_remove(out_fn) # Reset the appending of output
	
	# phe
	for(jj in 1:length(fls)){
		# jj = 1
		cat(paste0("Progress: ",jj," out of ",length(fls)," - Reading in ",fls[jj],"...\n"))
		snpj_fn 	= file.path(imputed_dir,fls[jj])
		infoj_fn 	= file.path(imputed_dir,paste0(fls[jj],"_info"))
		hpsj_fn 	= file.path(imputed_dir,paste0(fls[jj],"_haps"))
		
		# Import dosages, haplotypes, QC metrics
		cat("\tReading in dosages, haplotypes, QC metrics...\n")
		snpj = data.table::fread(snpj_fn,sep = " ",
			header = FALSE,stringsAsFactors = FALSE,
			quote = "",data.table = FALSE,showProgress = FALSE)
		infoj = data.table::fread(infoj_fn,sep = " ",
			header = TRUE,stringsAsFactors = FALSE,
			quote = "",data.table = FALSE,showProgress = FALSE)
		hpsj = data.table::fread(hpsj_fn,sep = " ",
			header = FALSE,stringsAsFactors = FALSE,
			quote = "",data.table = FALSE,showProgress = FALSE)
		
		# Filters
		snps_filter 	= snpj$V4 %in% bases & snpj$V5 %in% bases
		info_filter 	= infoj$info > info_thres
		maf_filter 		= infoj$exp_freq_a1 > min_af & infoj$exp_freq_a1 < 1-min_af
		final_filter 	= which(snps_filter & info_filter & maf_filter)
		length(final_filter)
		if( length(final_filter) == 0 ){
			cat(paste0("\tSkipping ",fls[jj],"!!!\n"))
			next
		}
		snpj 	= snpj[final_filter,]
		infoj = infoj[final_filter,]
		hpsj 	= hpsj[final_filter,]
		rm(snps_filter,info_filter,maf_filter,final_filter)
		
		pg = matrix(NA,nrow(snpj),sample_count) # phased genotype matrix
		colnames(pg) = samples
		rownames(pg) = paste0("chr",chr,":",snpj$V3,":",snpj$V4,":",snpj$V5)
		
		# Loop thru each sample to infer the genotype
		thres = 0
		cat("\t")
		for(ind in 1:sample_count){
			# ind = 1
			if(ind/sample_count*100 >= thres){
				cat(paste0(thres,"% "))
				thres = thres + 10
			}
			snp_index = off + (ind-1)*3 + 1 + c(0,1,2)
			hps_index = off + (ind-1)*2 + 1 + c(0,1)
			tmp_hps = hpsj[,hps_index]
			
			# choose the most likely genotype if its probability/dosage is larger than 0.8
			tmp_geno = Rcpp_dose2geno(
				doses = as.matrix(snpj[,snp_index]),
				geno_thres = geno_thres)[,1]
			tmp_geno[which(tmp_geno == 5)] = NA
			tmp_geno[which(tmp_geno == 2)] = 3
			geno_hetind = which(tmp_geno == 1)
			# head(tmp_hps)
			flip = tmp_hps[geno_hetind,1] == 1
			tmp_geno[geno_hetind[flip]] = 2
			# smart_table(tmp_geno)
			pg[,ind] = tmp_geno
			rm(tmp_hps,tmp_geno,geno_hetind,flip)
		}
		cat("\n")
		
		# Output and append to file
		app = file.exists(out_fn)
		pg2 = data.table::data.table(pg)
		rownames(pg2) = rownames(pg)
		data.table::fwrite(x = pg2,file = out_fn,sep = " ",
			na = "NA",row.names = TRUE,col.names = !app,
			quote = FALSE,append = app)
		
		# Clean up
		rm(pg,pg2,snpj,infoj,hpsj)
	}
	
	return(NULL)
}
step11_get_snp = function(work_dir,chr){
	if(FALSE){
		work_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/CMC"
		chr = 1
	}
	
	#ensure that when we write positions to a file R doesn't switch them to scientific
	options("scipen"=999,"digits"=4)
	
	#note that we also need to convert these snps to build 38
	#I have converted by liftover positions of those snps
	#also these snps build37 (under R302) build38 (under R312)
	#source("http://bioconductor.org/biocLite.R");biocLite("SNPlocs.Hsapiens.dbSNP.20120608")
	#library(SNPlocs.Hsapiens.dbSNP.20120608)
	#library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
	my_dirs = setdirs_SPI(work_dir = work_dir)
	phased_dir 					= my_dirs$phased_dir
	imputed_dir					= my_dirs$imputed_dir
	data_snp_dir				= my_dirs$data_snp_dir
	setwd(imputed_dir)
	
	# the list of samples for which we will get snps
	tmp_fn = file.path(phased_dir,paste0("phasing_chr",chr,".sample"))
	samples = data.table::fread(tmp_fn,header = TRUE,sep = " ",data.table = FALSE)
	samples = samples[-1,2]
	# define snps (I'll exclude indels with matching to the base)
	bases = c("A","C","G","T")
	
	# get all the appropriate chunks for the given chromosome
	fls = list.files(pattern = paste0("chr",chr,"_"))
	to_rm = union(union(union(grep(".txt",fls),grep("info_by_sample",fls)),
		grep("warnings",fls)),grep("summary",fls))
	fls = fls[-to_rm]
	hps = fls[grep("_haps",fls)]
	alp = fls[grep("_allele_probs",fls)]
	info = fls[grep("info$",fls)]
	fls = setdiff(fls,union(union(hps,alp), info))
	# fls; info; hps; alp
	
	append = FALSE
	j = 1; off = 5
	probs = c(.9,.95,.99,.995,.999)
	
	sample_count = length(samples)
	summ_het = matrix(0,nrow=sample_count,ncol=length(probs))
	rownames(summ_het) = samples
	colnames(summ_het) = paste0("p",probs)
	tot = rep(0,sample_count)
	
	# phe
	for(jj in 1:length(fls)){
		# jj = match("phased_imputed_chr2_85e6_90e6",fls)+1
		# jj = 1
		
		cat(paste0("Progress: ",jj," out of ",length(fls)," - Reading in ",fls[jj],"...\n"))
		
		# flj = 
		# snpj = smart_RT(fls[jj],as.is = TRUE)
		# hpj = smart_RT(hps[jj],as.is = TRUE)
		# alj = smart_RT(alp[jj],as.is = TRUE)
		# infoj = smart_RT(info[jj],as.is = TRUE,header = TRUE)
		snpj_fn 	= file.path(imputed_dir,fls[jj])
		hpsj_fn 	= file.path(imputed_dir,paste0(fls[jj],"_haps"))
		alj_fn 		= file.path(imputed_dir,paste0(fls[jj],"_allele_probs"))
		infoj_fn 	= file.path(imputed_dir,paste0(fls[jj],"_info"))
		
		snpj = data.table::fread(snpj_fn,sep = " ",
			header = FALSE,stringsAsFactors = FALSE,
			quote = "",data.table = FALSE,showProgress = FALSE)
		hpj = data.table::fread(hpsj_fn,sep = " ",
			header = FALSE,stringsAsFactors = FALSE,
			quote = "",data.table = FALSE,showProgress = FALSE)
		alj = data.table::fread(alj_fn,sep = " ",
			header = FALSE,stringsAsFactors = FALSE,
			quote = "",data.table = FALSE,showProgress = FALSE)
		infoj = data.table::fread(infoj_fn,sep = " ",
			header = FALSE,stringsAsFactors = FALSE,
			quote = "",data.table = FALSE,showProgress = FALSE)
		
		snpj$V4 = as.character(snpj$V4)
		snpj$V5 = as.character(snpj$V5)

		rm_indels = which(snpj$V4%in%bases & snpj$V5%in%bases)
		snpj = snpj[rm_indels,]
		alj = alj[rm_indels,]
		hpj = hpj[rm_indels,]
		infoj = infoj[rm_indels,]
		# dim(hpj)
		# dim(infoj)
		# dim(alj)
		# dim(snpj)

		# remove R2 < 0.3  "'info' is similar to the r-squared 
		# metrics reported by other programs like MaCH and Beagle"

		R2kp = which(infoj$info > 0.3)
		snpj = snpj[R2kp,]
		alj = alj[R2kp,]
		hpj = hpj[R2kp,]
		infoj = infoj[R2kp,]
		# dim(hpj)
		# dim(infoj)
		# dim(alj)
		# dim(snpj)
		
		thres = 0
		cat("\t")
		for(ind in 1:sample_count){
			# ind = 1
			# ind = ind + 1
			if(ind/sample_count*100 >= thres){
				cat(paste0(thres,"% "))
				thres = thres + 10
			}
			samp_geno_dir = file.path(data_snp_dir,samples[ind])
			smart_mkdir(samp_geno_dir)
			hetloc = off + (ind-1)*3+2
			indloc = off + (ind-1)*2+1
			for(probk in 1:length(probs)){
				summ_het[ind,probk] = summ_het[ind,probk] + sum(snpj[,hetloc] >= probs[probk])
			}
			tot[ind] = tot[ind] + nrow(snpj)
			# I'll use only the snps with probability>0.999
			flag = snpj[,hetloc] > probs[5]
			if( sum(flag) > 0 ){
				snpind = snpj[flag,]
				alind = alj[flag,c(1:5,indloc,indloc+1)]
				hpind = hpj[flag,c(1:5,indloc,indloc+1)]
				infoind = infoj[flag, ]
				flip = (hpind[,6]==1)
				hpind[,1] = paste0("chr",chr)
				hpind[flip,4:5] = hpind[flip,5:4]     
				outfile = file.path(samp_geno_dir,paste0("chr",chr,".txt"))
				app = file.exists(outfile)
				write.table(hpind[,c(1,3,4,5,2)],outfile,
					row.names = FALSE,col.names = FALSE,
					quote = FALSE,append = app)
			}
		}
		cat("\n")
		
		# message(flj)
	}

	# summary of
	outfile = file.path(data_snp_dir,paste0("chr_",chr,"_num_snps_1.txt"))
	write.table(cbind(summ_het,tot),outfile,
		row.names = TRUE,col.names = TRUE,
		quote = FALSE)
	
}
step12_checkHet_probesGenoCalls_vs_ImputedSnps = function(work_dir,sample_id){
	if(FALSE){
		work_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/CMC"
		sample_id = "0_MSSM_13"
	}
	
	my_dirs = setdirs_SPI(work_dir = work_dir)
	data_snp_dir = my_dirs$data_snp_dir
	
	# Checking the pipeline input genotype calls vs. the imputed output
	genotype_fn = file.path(work_dir,"geno.call.traw")
	snp6_anno_fn = file.path(work_dir,"GenomeWideSNP_6.na35.annot.csv")
	
	# 1.1 Load the SNP 6 csv file
	# num_lines = as.numeric(system(sprintf("wc -l %s | cut -d ' ' -f1",snp6_anno_fn),intern = TRUE)) - 1
	# snpcsv = smart_RT(snp6_anno_fn,sep = ',',header = TRUE,as.is = TRUE,nrows = num_lines)
	# snpcsv = snpcsv[,c(1:4,9:10)]
	# head(snpcsv)
	
	# Get original observed genotypes
	row_idx = data.table::fread(file.path(work_dir,"chr","geno.chr1.ped"),sep = " ",header = FALSE)
	row_idx = paste0(row_idx$V1,"_",row_idx$V2); row_idx = which(row_idx == sample_id)
	geno_pre_phase = c()
	for(chr in seq(22)){
		# chr = 1
		cat(sprintf("Chr = %s\n",chr))
		cmd = sprintf("sed -n %sp %s",row_idx,file.path(work_dir,"chr",sprintf("geno.chr%s.ped",chr)))
		genos = system(cmd,intern = TRUE); genos = strsplit(genos," ")[[1]]
		genos = genos[-c(1:6)]; genos = matrix(genos,ncol=2,byrow=TRUE)
		genos = paste0(genos[,1],genos[,2])
		genos[genos == "00"] = NA
		map = data.table::fread(file.path(work_dir,"chr",sprintf("geno.chr%s.map",chr)),
			sep = "\t",header = FALSE,data.table = FALSE)
		names(map) = c("Chromosome","dbSNP.RS.ID","V3","Physical.Position")
		map$pre_phased_geno = genos; rm(genos)
		geno_pre_phase = rbind(geno_pre_phase,map)
	}
	
	# 3. Load the combined SNPs list
	sample_id_2 = paste(strsplit(sample_id,"_")[[1]][-1],collapse="_")
	# snps_dir = dir(path=data_snp_dir,pattern=sample_id_2,full.names=TRUE)
	snps_dir = file.path(data_snp_dir,sample_id_2)
	setwd(snps_dir)
	system(smart_sprintf("
		cd %s; 
		[ -f combined.txt ] && rm combined.txt; 
		for chr in `seq 1 22`; do 
			cut -d' ' -f1-4 chr${chr}.txt | sed -n 's/^chr//p' >> combined.txt; 
		done",snps_dir)) # rbind snp files for all chr
	genotype_imputed = data.table::fread(file.path(snps_dir,"combined.txt"),
		header = FALSE,sep = " ",data.table = FALSE)
	names(genotype_imputed) = c("Chromosome","Physical.Position","hap1","hap2")
	print(dim(genotype_imputed)); print(head(genotype_imputed))
	genotype_imputed$imputed_phased_geno = paste0(genotype_imputed$hap1,genotype_imputed$hap2)

	# 4. Merge the two lists
	genotype_merged = merge(x = geno_pre_phase,y = genotype_imputed,
		by=c("Chromosome","Physical.Position"))
	print(nrow(genotype_merged)); print(head(genotype_merged))
	print(smart_table(genotype_merged$pre_phased_geno)); print(smart_table(genotype_merged$imputed_phased_geno))
	print(smart_table(genotype_merged[,c("pre_phased_geno","imputed_phased_geno")]))
	
}
step12a_check_ourImpute_studyImpute = function(work_dir,sample_id){
	if(FALSE){
		work_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s4_CMC/geno"
		sample_id = "0_MSSM_13"
	}
	
	options(width=120)
	cmc_impute_dir = file.path(work_dir,"cmc_impute")
	sample_id_2 = paste(strsplit(sample_id,"_")[[1]][-1],collapse="_"); sample_id_2
	bases = c("A","C","G","T")
	info_thres = 0.3; min_af = 0.05
	get_phased_imputed_geno = function(hap1,hap2,geno){
		if( is.na(geno) ){
			NA
		} else if(geno == 0){
			paste0(hap1,hap1)
		} else if(geno == 1){
			paste0(hap1,hap2)
		} else if(geno == 2){
			paste0(hap2,hap1)
		} else if(geno == 3){
			paste0(hap2,hap2)
		} else {
			print(hap1); print(hap2); print(geno)
			stop("Something's wrong")
		}
	}
	
	for(chr in seq(22)){
		# chr = 22
		cat(sprintf("Chr = %s\n",chr))
		ci_chr_fn = file.path(cmc_impute_dir,sprintf("chr%s",chr))
		fls = list.files(ci_chr_fn); fls = fls[grep(".gen$",fls)]
		samps_fn = file.path(ci_chr_fn,sprintf("CMC_MSSM-Penn-Pitt_DLPFC_DNA_imputed_chr%s.sample",chr))
		samps = data.table::fread(samps_fn,sep = " ",header = TRUE,data.table = FALSE)
		samp_idx = which(samps$ID_2 == sample_id_2)
		col_idx = 4 + 0:2 + samp_idx; col_idx
		
		geno = c()
		for(jj in seq(length(fls))){
			# jj = 1
			cat(sprintf("\tProgress: %s out of %s\n",jj,length(fls)))
			gen_fn = file.path(ci_chr_fn,fls[jj])
			info_fn = sprintf("%s_info",gen_fn)
			gen = data.table::fread(gen_fn,sep = " ",header = FALSE,data.table = FALSE)
			gen = gen[,c(1:5,col_idx)]
			names(gen) = c("V1","dbSNP_ID","Position","Ref","Alt",paste0("post",0:2))
			info = data.table::fread(info_fn,sep = " ",header = TRUE,data.table = FALSE)
			info = info[,c("exp_freq_a1","info","certainty")]
			
			# Filter snps: base subs, high info, high dosage
			snps_filter 	= gen$Ref %in% bases & gen$Alt %in% bases
			info_filter 	= info$info > info_thres
			maf_filter 		= info$exp_freq_a1 > min_af & info$exp_freq_a1 < 1-min_af
			final_filter 	= which(snps_filter & info_filter & maf_filter)
			
			# info = info[final_filter,]
			gen = gen[final_filter,]
			infer_geno = t(apply(gen[,paste0("post",0:2)],1,function(xx) c(which.max(xx)-1,xx[which.max(xx)])))
			infer_geno = smart_df(infer_geno)
			names(infer_geno) = c("geno","max_dose")
			gen = cbind(gen,infer_geno)
			gen$cmc_geno = apply(gen[,c("Ref","Alt","geno")],1,function(xx) ifelse(xx[3]==0,paste0(xx[1],xx[1]),
				ifelse(xx[3]==1,paste0(xx[1],xx[2]),paste0(xx[2],xx[2]))))
			gen$cmc_geno[which(gen$max_dose < 0.8)] = NA
			gen$ChrPos = paste0("chr",chr,":",gen$Position)
			geno = rbind(geno,gen)
		}
		
		# Get our pipeline's imputed results
		data_geno_dir = file.path(work_dir,"data_genotype")
		chr_fn = file.path(data_geno_dir,sprintf("chr%s.txt",chr))
		res = data.table::fread(chr_fn,sep = " ",header = TRUE,data.table = FALSE)
		tmp_names = readLines(chr_fn,n=1); tmp_names = strsplit(tmp_names," ")[[1]][-1]
		names(res) = c("ChrPosHap1Hap2",tmp_names)
		res = res[,c("ChrPosHap1Hap2",sample_id_2)]
		res[,sample_id_2] = gsub(" ","",res[,sample_id_2])
		res$ChrPos = sapply(res$ChrPosHap1Hap2,function(xx) 
			paste(strsplit(xx,":")[[1]][1:2],collapse=":"),USE.NAMES = FALSE)
		res$hap1 = sapply(res$ChrPosHap1Hap2,function(xx) strsplit(xx,":")[[1]][3],USE.NAMES = FALSE)
		res$hap2 = sapply(res$ChrPosHap1Hap2,function(xx) strsplit(xx,":")[[1]][4],USE.NAMES = FALSE)
		res$our_geno = apply(res[,c("hap1","hap2",sample_id_2)],1,function(xx)
			get_phased_imputed_geno(hap1 = xx[1],hap2 = xx[2],geno = xx[3]))
		
		# Merge
		mm = smart_merge(res,geno)
		cat(sprintf("\tNumber of our workflow's loci = %s\n",nrow(res)))
		cat(sprintf("\tNumber of CMC workflow loci = %s\n",nrow(geno)))
		cat(sprintf("\tNumber of shared loci = %s\n",nrow(mm)))
		print(smart_table(mm[,c("our_geno","cmc_geno")]))
		
	}
	
	
}
step13_liftOver_hg19_hg38 = function(liftover_dir,work_dir){
	if(FALSE){
		liftover_dir = liftover_dir
		work_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s4_CMC/geno"
	}
	
	my_dirs 						= setdirs_SPI(work_dir = work_dir)
	data_genotype_dir		= my_dirs$data_genotype_dir
	dg38_dir						= my_dirs$dg38_dir
	# setwd(dg38_dir)
	
	chain_fn = file.path(liftover_dir,"hg19ToHg38.over.chain")
	chain_url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
	smart_DL(down_dir = liftover_dir,url_link = chain_url)
	
	for(chr in seq(22)){
		# chr = 1
		cat(sprintf("Chr = %s\n",chr))
		
		chrpos_fn = file.path(dg38_dir,sprintf("map_hg19_hg38_chr%s.tsv",chr))
		if( !file.exists(chrpos_fn) ){
			geno_chr_fn = file.path(data_genotype_dir,sprintf("chr%s.txt",chr))
			chrpos = data.table::fread(geno_chr_fn,sep = " ",header = FALSE,data.table = FALSE,select = 1)
			names(chrpos) = "ChrPosHap1Hap2_hg19"
			chrpos$Chr_hg19 = sapply(chrpos$ChrPosHap1Hap2_hg19,function(xx)
				strsplit(xx,":")[[1]][1],USE.NAMES = FALSE)
			chrpos$Position_hg19 = as.integer(sapply(chrpos$ChrPosHap1Hap2_hg19,function(xx) 
				strsplit(xx,":")[[1]][2],USE.NAMES = FALSE))
			chrpos = chrpos[order(chrpos$Position_hg19),]
			chrpos$Start_hg19 = as.integer(chrpos$Position_hg19 - 1)
			chrpos$End_hg19 = as.integer(chrpos$Position_hg19)
			chrpos$ChrPos_hg19 = paste0(chrpos$Chr_hg19,":",chrpos$End_hg19)
			chrpos$hap1 = sapply(chrpos$ChrPosHap1Hap2_hg19,function(xx)
				strsplit(xx,":")[[1]][3],USE.NAMES = FALSE)
			chrpos$hap2 = sapply(chrpos$ChrPosHap1Hap2_hg19,function(xx)
				strsplit(xx,":")[[1]][4],USE.NAMES = FALSE)
			# head(chrpos)
			chrpos_input = chrpos[,c("Chr_hg19","Start_hg19","End_hg19")]
			names(chrpos_input) = c("Chr","Start","End")
			# dim(chrpos_input); head(chrpos_input); str(chrpos_input)
			chr_hg19_fn = file.path(dg38_dir,sprintf("input_chr%s_hg19.txt",chr))
			smart_WT(chrpos_input,chr_hg19_fn,sep = " ",col.names = FALSE)
			
			# Run liftOver
			chr_hg38_fn = file.path(dg38_dir,sprintf("output_chr%s_hg38.txt",chr))
			chr_unlift_fn = file.path(dg38_dir,sprintf("output_chr%s_unlifted.txt",chr))
			cmd = sprintf("%s/liftOver %s %s %s %s",liftover_dir,chr_hg19_fn,chain_fn,chr_hg38_fn,chr_unlift_fn)
			# print(cmd); 
			system(cmd)
			
			# Import and map coordinates
			chrpos_hg38 = data.table::fread(chr_hg38_fn,sep = "\t",header = FALSE,data.table = FALSE)
			names(chrpos_hg38) = c("Chr_hg38","Start_hg38","End_hg38")
			chrpos_unlift = smart_RT(chr_unlift_fn,sep = "\t",header = FALSE)
			names(chrpos_unlift) = c("Chr_hg19","Start_hg19","End_hg19")
			chrpos_unlift$ChrPos_hg19 = paste0(chrpos_unlift$Chr_hg19,":",chrpos_unlift$End_hg19)
			chrpos$Chr_hg38 = NA; chrpos$Start_hg38 = NA; chrpos$End_hg38 = NA
			unlift_idx = which(chrpos$ChrPos_hg19 %in% chrpos_unlift$ChrPos_hg19)
			hg38_cols = c("Chr_hg38","Start_hg38","End_hg38")
			chrpos[-unlift_idx,hg38_cols] = chrpos_hg38[,hg38_cols]
			chrpos = chrpos[which(!is.na(chrpos$End_hg38)),]
			chrpos = chrpos[which(chrpos$Chr_hg38 %in% paste0("chr",1:22)),]
			chrpos$ChrPosHap1Hap2_hg38 = paste0(chrpos$Chr_hg38,":",chrpos$End_hg38,
				":",chrpos$hap1,":",chrpos$hap2)
			chrpos = chrpos[order(chrpos$End_hg38),]
			# smart_WT(chrpos,file.path(dg38_dir,sprintf("map_chr%s.txt",chr)),sep = "\t")
			smart_remove(chr_hg19_fn); smart_remove(chr_hg38_fn); smart_remove(chr_unlift_fn)
			
			smart_WT(chrpos,chrpos_fn,sep = "\t")
		}
		
	}
	
}
step14_prepChrGeno_hg38 = function(work_dir,chr){
	if(FALSE){
		work_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s4_CMC/geno"
		chr = 10
	}
	
	my_dirs 						= setdirs_SPI(work_dir = work_dir)
	data_genotype_dir		= my_dirs$data_genotype_dir
	dg38_dir						= my_dirs$dg38_dir
	
	# chr = 22
	
	cat(sprintf("Chr = %s\n",chr))
	geno38_chr = c()
	
	for(chr2 in seq(22)){
		# chr2 = 2
		cat(sprintf("\t%s\n",chr2))
		# Import hg19 chr file
		geno19_fn = file.path(data_genotype_dir,sprintf("chr%s.txt",chr2))
		geno19 = data.table::fread(geno19_fn,sep = " ",header = FALSE,data.table = FALSE)
		tmp_names = readLines(geno19_fn,n = 1); tmp_names = strsplit(tmp_names," ")[[1]]
		tmp_names[1] = "ChrPosHap1Hap2_hg19"
		names(geno19) = tmp_names
		
		# Import hg19 to hg38 map file
		map_fn = file.path(dg38_dir,sprintf("map_hg19_hg38_chr%s.tsv",chr2))
		map = data.table::fread(map_fn,sep = "\t",header = TRUE,data.table = FALSE)
		smart_table(map$Chr_hg38)
		
		# Subset specific chr#
		map = map[which(map$Chr_hg38 == sprintf("chr%s",chr)),]
		if(nrow(map) == 0){
			rm(geno19,map)
			next
		}
		dim(map)
		
		int_loci = intersect(geno19$ChrPosHap1Hap2_hg19,map$ChrPosHap1Hap2_hg19)
		geno19 = geno19[which(geno19$ChrPosHap1Hap2_hg19 %in% int_loci),]
		map = map[which(map$ChrPosHap1Hap2_hg19 %in% int_loci),]
		
		geno19 = geno19[match(map$ChrPosHap1Hap2_hg19,geno19$ChrPosHap1Hap2_hg19),]
		all(geno19$ChrPosHap1Hap2_hg19 == map$ChrPosHap1Hap2_hg19)
		geno19[,1] = map$ChrPosHap1Hap2_hg38
		names(geno19)[1] = "ChrPosHap1Hap2_hg38"
		
		geno38_chr = rbind(geno38_chr,geno19); rm(geno19,map)
	}
	
	# Order chromosome's positions: Double check for multiple alleles at same locus
	tmp_df = smart_df(t(sapply(geno38_chr$ChrPosHap1Hap2_hg38,function(xx) 
		strsplit(xx,":")[[1]],USE.NAMES = FALSE)))
	names(tmp_df) = c("Chr","Pos","hap1","hap2")
	tmp_df$Pos = as.integer(tmp_df$Pos)
	tmp_df$ChrPosHap1Hap2_hg38 = geno38_chr$ChrPosHap1Hap2_hg38
	tmp_df$PropMiss = apply(geno38_chr[,-1],1,function(xx) mean(is.na(xx)))
	tmp_df = tmp_df[which(tmp_df$PropMiss < 0.05),]
	
	tmp_df = tmp_df[!duplicated(tmp_df$ChrPosHap1Hap2_hg38),]
	geno38_chr = geno38_chr[which(geno38_chr$ChrPosHap1Hap2_hg38 
		%in% tmp_df$ChrPosHap1Hap2_hg38),]
	
	geno38_chr = geno38_chr[!duplicated(geno38_chr$ChrPosHap1Hap2_hg38),]
	
	tt = table(tmp_df$ChrPosHap1Hap2_hg38)
	if( length(which(tt > 1)) > 0 ){
		stop("Multiple ChrPosHap1Hap2_hg38!")
	}
	
	tmp_df = tmp_df[order(tmp_df$Pos,tmp_df$hap1,tmp_df$hap2),]
	geno38_chr = geno38_chr[which(geno38_chr$ChrPosHap1Hap2_hg38 %in% tmp_df$ChrPosHap1Hap2_hg38),]
	geno38_chr = geno38_chr[match(tmp_df$ChrPosHap1Hap2_hg38,geno38_chr$ChrPosHap1Hap2_hg38),]
	rm(tmp_df)
	rownames(geno38_chr) = geno38_chr$ChrPosHap1Hap2_hg38
	geno38_chr = geno38_chr[,-1]
	data.table::fwrite(geno38_chr,file.path(dg38_dir,sprintf("chr%s.txt",chr)),
		sep = " ",quote = FALSE,row.names = TRUE,col.names = TRUE,na = "NA")
	cat("\n")
	
}

if(FALSE){

rm(list=ls()); git_dir = "~/github"
source(file.path(git_dir,"CARseq_pipelines/R/base_functions_plittle.R"))
source(file.path(git_dir,"CARseq_pipelines/R/snp_phase_impute.R"))
cs_dir = "/pine/scr/p/l/pllittle/CS_eQTL"
cmc_dir = file.path(cs_dir,"s3_Real/CMC")
my_dirs = setdirs_SPI(work_dir = cmc_dir)




}


