There are the codes to analyze gene expression data of SCZ vs. control from CMC, as well as gene expression data of ASD vs. control from UCLA study. Here is the list of all the files, and we provide some explanation for each file. 

├── ASD_gene_annotations.R
├── ASD_step1_CARseq_seqSV4.R
├── ASD_step1_CARseq_seqSV4_permuted.R
├── ASD_step1_TOAST_and_DESeq2_SeqSV4.R
├── SCZ_QC_PCA.R
├── SCZ_QC_PCA.Rout
├── SCZ_RNA_QC.R
├── SCZ_RNA_QC.Rout
├── SCZ_SVA.R
├── SCZ_SVA.Rout
├── SCZ_gene_annotations.R
├── SCZ_snp_phase_impute.R
├── SCZ_step1_CARseq_SV2.R
├── SCZ_step1_CARseq_SV2_permuted.R
├── SCZ_step1_TOAST_DESeq2_SV2.R
├── base_functions_plittle.R
├── simulation_step1_get_distribution_of_parameters_from_real_data.R
├── simulation_step2_simulate_data.R
├── simulation_step3_test_methods.R
├── simulation_step4_compare_methods_multiple_replicates.R
├── simulation_step4_compare_methods_one_replicate.R
├── simulation_submit_step1.sh
├── simulation_submit_step3.sh
├── simulation_test_CARseq.R
├── simulation_test_CIBERSORTx.R
├── simulation_test_TOAST_count.R
├── simulation_test_TOAST_TPM.R
├── simulation_test_csSAM.R
├── compute_MTG_size_factor.R
└── summarize_pval.R

simulation_step1_get_distribution_of_parameters_from_real_data.R
	Get distribution of parameters from (unfiltered) CMC schizophrenia data.

simulation_step2_simulate_data.R
	Simulate data. The simulated data are not included in the repository
	 because they were very large in file sizes.

simulation_step3_test_methods.R
	Conduct tests on simulation data using CARseq, TOAST, csSAM, etc.

simulation_step4_compare_methods_multiple_replicates.R
	Collect the DE test results from 10 replicates and make plots.

simulation_step4_compare_methods_one_replicate.R
	Collect the DE test results from 1 replicate and make plots.
	This covers CIBERSORTx and TOAST using read counts.

simulation_submit_step1.sh
	Submit "simulation_step1_get_distribution_of_parameters_from_real_data.R" in batch using SLURM.

simulation_submit_step3.sh
	Submit "simulation_step3_test_methods.R" in batch using SLURM.

simulation_test_CARseq.R
	A wrapper of CARseq that need to be sourced when conducting DE tests on simulation data.

simulation_test_CIBERSORTx.R
	A script that processes the CIBERSORTx high resolution model output
	that need to be sourced when conducting DE tests on simulation data.

simulation_test_TOAST_count.R
	A wrapper of TOAST that takes expression in read counts
	that need to be sourced when conducting DE tests on simulation data.

simulation_test_TOAST_TPM.R
	A wrapper of TOAST that takes expression in TPM
	that need to be sourced when conducting DE tests on simulation data.

simulation_test_csSAM.R
	A wrapper of csSAM that need to be sourced when conducting DE tests on simulation data.

ASD_gene_annotations.R
	Add gene annotations to create a SummarizedExperiment object for the ASD dataset.

ASD_step1_CARseq_seqSV4.R
	Run CARseq on the ASD dataset.

ASD_step1_CARseq_seqSV4_permuted.R
	Run CARseq on the ASD dataset when the case-control label is permuted.
	This is an approximation for the null distribution of p-values.
	There might be problems because the cell fraction is not permuted.

ASD_step1_TOAST_and_DESeq2_SeqSV4.R
	Run TOAST on expression in TPM and DESeq2 on the ASD dataset.

SCZ_QC_PCA.R
	Quality Control and PCA on gene expression (Total Read Counts per gene 
	per sample) and genotype data.

SCZ_RNA_QC.R
	Quality Control for RNA-seq data based on Total Read Counts and 
	Allele-specific read counts, this is only revelant for eQTL mapping

SCZ_snp_phase_impute.R
	phasing and imputation of SNP genotype data. This is only revelant for 
	eQTL mapping.

SCZ_SVA.R
	Surrogate variable analysis.

SCZ_gene_annotations.R
	Add gene annotations to create a SummarizedExperiment object for the SCZ dataset.

SCZ_step1_CARseq_SV2.R
	Run CARseq on the SCZ dataset.

SCZ_step1_CARseq_SV2_permuted.R
	Run CARseq on the SCZ dataset when the case-control label is permuted.
	This is an approximation for the null distribution of p-values.
	There might be problems because the cell fraction is not permuted.

SCZ_step1_TOAST_DESeq2_SV2.R
	Run TOAST on expression in TPM and DESeq2 on the SCZ dataset.

base_functions_plittle.R
	A suite of practical functions for plotting, data manipulation, etc.
	maintained by Dr. Paul Little <https://github.com/pllittle>.

compute_MTG_size_factor.R
	Calculate size factor from genomewide TPM in MTG single cell data.

summarize_pval.R
	Compare differential expression analysis results in ASD vs. SCZ.
