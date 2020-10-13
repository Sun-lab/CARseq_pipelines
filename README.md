# CARseq_pipelines

## Data availabilty

The following data files used in this analysis (in data folder) has restricted access and thus cannot be shared in this GitHub. 

* UCLA-ASD processed data
    
    * data/ASD_rse_filtered_with_SVs.rds
    * data/ASD_rse_filtered.rds
    
    * produced using code in data/UCLA_ASD/step2_check_RSEM_ucla.R
    
        * data/ucla_eDat.rds
        * data/ucla_cData.ds 

* CMC-SCZ raw data

    * data/CMC-CMC_HBCC_clinical_.csv
    * data/Release3_SampleID_key.csv
    * data/CMC_MSSM-Penn-Pitt_Paul_geneExpressionRaw.rds

 * CMC-SCZ processed data
 
    * data/dat_cavariates_scz_control_with_svs.rds
    * data/dat_cavariates.rds
    * data/rse.rds
    * data/rse_filtered_SV.rds
    * data/trec_dat_post_QC_PCA_SCZ.rds
    * data/trec_dat_post_QC_PCA_Control.rds
    * data/trec_filtered_scz_control.rds
    * data/trec_filtered.rds
    * data/trec_filtered_log_rd_corrected_scz_control.rds
    * data/trec_filtered_log_rd_corrected.rds



Here are the URLs to download the relevant data

* CMC gene expression data: https://www.synapse.org/#!Synapse:syn3346749
  * CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv.gz

* CMC gene expression meta data: https://www.synapse.org/#!Synapse:syn18103174
  * CMC_Human_rnaSeq_metadata.csv (the same as CMC_Human_DLPFC_rnaSeq.csv, except with an additional column: rRNA_Rate)
  
* CMC clinical data: https://www.synapse.org/#!Synapse:syn12080241
  * CMC_Human_clinical_metadata.csv (the same as CMC-CMC_HBCC_clinical_.csv)
  * CMC_Human_sampleIDkey_metadata (the same as Release3_SampleID_key.csv)

* UCLA-ASD gene expression data: https://www.synapse.org/#!Synapse:syn8365527
    * RSEM_Quant.genes.results for each sample, e.g., UCLA_ASD/318/14855318/AN14757_vermis_39.RSEM_Quant.genes.results

* UCLA-ASD gene expression meta data: https://www.synapse.org/#!Synapse:syn5602933
    * UCLA_R01MH094714_ASD_Metadata_RNAseq.csv
    
* UCLA-ASD clinical data: https://www.synapse.org/#!Synapse:syn5602932
    * UCLA_R01MH094714_ASD_Metadata_Clinical_August2016Release.csv


## Reproducible figures

Most of the reproducible figures can be generated using
`reproducible_figures.Rmd`.

## _p_-values from real data analysis

Please check `results/_pvalues`.

## MTG

Information of cell type-specific gene expression or cell type proportions derived from [Allen Brain Atlas MTG (middle temporal gyrus) dataset](http://celltypes.brain-map.org/rnaseq). Some details of their analysis is provided [here](http://help.brain-map.org/download/attachments/8323525/CellTypes_Transcriptomics_Overview.pdf). Briefly, these data were generated using SMART-Seq v4 Ultra Low Input RNA Kit, which is an improved version of SMART-seq2 protocol. 

We have explored this dataset in the other repo: https://github.com/Sun-lab/scRNAseq_pipelines/tree/master/MTG

## Schizophrenia dataset

1. Annotation and cell fractions

Surrogate variables are added to the dataset in `R/SVA.R`.

Gene annotations are added in `R/SCZ_step0_gene_annotations.R`.

The cell fractions are estimated using ICeDT and CIBERSORT using
MTG signature genes as reference,
and they are stored in `data/prop_from_TPM.rds`.
After adjustment using cell sizes calculated using the sum of
TPMs across the signature genes for each cell type,
the cell fractions are stored in `data/prop_resized.rds`.

2. Run CARseq

    - `R/SCZ_step1_CARseq_SV2.R`

    - `R/SCZ_step1_CARseq_SV2_permuted.R`

3. Run TOAST and DESeq2
    
    - `R/SCZ_step1_TOAST_and_DESeq2_SV2.R`

## Autism spectrum disorder (ASD) dataset

1. Annotation and cell fractions

Surrogate variables are added to the dataset and gene annotations
 are added in `R/ASD_step0_gene_annotations.R`.

The cell fractions are estimated using ICeDT and CIBERSORT using
MTG signature genes as reference,
and they are stored in `data/ASD_prop_from_TPM.rds`.
After adjustment using cell sizes calculated using the sum of
TPMs across the signature genes for each cell type,
the cell fractions are stored in `data/ASD_prop.rds`.

2. Run CARseq

    - `R/ASD_step1_CARseq_seqSV4.R`

    - `R/ASD_step1_CARseq_seqSV4_permuted.R`

3. Run TOAST and DESeq2
    
    - `R/ASD_step1_TOAST_and_DESeq2_seqSV4.R`

## Simulation

1. Get distribution of parameters from CMC schizophrenia data

    - `R/simulation_submit_step1.sh`

    - `R/simulation_step1_get_distribution_of_parameters_from_real_data.R`

2. Simulate data

    - `R/simulation_step2_simulate_data.R`

3. Conduct tests using CARseq, TOAST and csSAM, as well as checking the validity of methods in
   terms of replicability and when noise is present:
    
    - `R/simulation_submit_step3.sh`
    
    - `R/simulation_step3_test_method.R`

       The script sources the following scripts:

        + `R/simulation_test_CARseq.R`

        + `R/simulation_test_TOAST_TPM.R`

        + `R/simulation_test_csSAM.R`

        + `R/simulation_test_CIBERSORTx.R`
      
    - `R/simulation_step3_test_methods_reproducibility.R`
    
    - `R/simulation_step3_test_methods_with_noise.R`
    
4. Summarize the simulation results using tables and plots
    
    - `R/simulation_step4_compare_methods_multiple_replicates.R`
    
    - `R/simulation_step4_compare_methods_one_replicate.R`
  
    - `R/simulation_step4_compare_methods_with_noise.R`

    - Simulationa differential expression (DE) testing results, in tables:

      + `results/compare_method_metrics_multiple_replicates.csv`

    - Simulation DE test results, in figures:

      + `figures/compare_methods_fdr_multiple_replicates_draft_version_p1.pdf`

      + `figures/compare_methods_fdr_multiple_replicates_draft_version_p2.pdf`

      + `figures/compare_methods_fdr_multiple_replicates_draft_version_p3.pdf`

      + `figures/compare_methods_fdr_multiple_replicates_draft_version_p4.pdf`
      
    - Simulation DE test ratio of sensitivity between TOAST and CARseq, in figures:
      
      + `figures/ratio_of_sensitivity_TOAST_over_CARseq_p1.pdf`
      
      + `figures/ratio_of_sensitivity_TOAST_over_CARseq_p2.pdf`
      
      + `figures/ratio_of_sensitivity_TOAST_over_CARseq_p3.pdf`
      
      + `figures/ratio_of_sensitivity_TOAST_over_CARseq_p4.pdf`

    - Simulation reproducibility measured by effect size estimation, in figures:
      
      + `figures/reproducibility_n_50_DE_pattern_2_1_1_replicate_1.pdf`
      
      + `figures/reproducibility_n_100_DE_pattern_2_1_1_replicate_1.pdf`
      
      + `figures/reproducibility_n_200_DE_pattern_2_1_1_replicate_1.pdf`
      
    - Simulation DE test results when size factor is misspecified, in figures:
      
      + `figures/compare_methods_fdr_size_factor_1.2_draft_version_p1.pdf`
      
      + `figures/compare_methods_fdr_size_factor_2_draft_version_p1.pdf`
      
    - Simulation DE test results when noise has been added to cell fractions, in figures:
    
      + `figures/compare_methods_fdr_with_noise_draft_version_p1.pdf`

