├── dat_cavariates_scz_control_with_svs.rds
├── trec_filtered_log_rd_corrected_scz_control.rds
└── trec_filtered_scz_control.rds

------------------------------------------------------------------
aat8127_Table_S1.xlsx
------------------------------------------------------------------

Supplementary Table S1 from Gandal et al (2018) Transcriptome-wide 
isoform-level dysregulation in ASD, schizophrenia, and bipolar disorder. 
Science 362(6420) eaat8127

------------------------------------------------------------------
DER-23_Cell_fractions_Raw.xlsx
DER-24_Cell_fractions_Normalized.xlsx
------------------------------------------------------------------

Results from Wang et al (2018) Comprehensive functional genomic 
resource and integrative model for the human brain. 
Science 362(6420), eaat8464

Downloaded from http://resource.psychencode.org/

------------------------------------------------------------------
dat_cavariates_scz_control_with_svs.rds
------------------------------------------------------------------

An R data file of the covaraites per sample for data anlaysis. This is 
already in the format of a design matrix (i.e., factors are encoded 
as indicators for different levels). Here are the frist two lines of 
this data matrix


> dim(modDat)
[1] 527  33
> modDat[1:2,]
   DxSCZ genderMale InstitutionPenn InstitutionPitt libclustB libclustbase
21     1          1               0               0         0            0
22     0          1               0               0         1            0
   libclustC libclustD libclustE libclustF libclustG age_death  PMI RIN  RIN2
21         0         1         0         0         0        68  8.9 7.5 56.25
22         0         0         0         0         0        58 12.3 8.8 77.44
         genoPC1     genoPC2    genoPC3      genoPC4      genoPC5 log_depth
21 -0.0205963704 -0.02205074 -0.0105337 -0.007387133 -0.004576622  6.410175
22 -0.0004623463 -0.01721824  0.3105687 -0.085514534  0.100494270  7.060691
           sv1        sv2        sv3         sv4        sv5         sv6
21  0.02169023 0.01800533 0.03972740 -0.06987432 0.01924304 -0.01205863
22 -0.05525339 0.01258435 0.01896169 -0.04361832 0.07727248 -0.02826308
          sv7         sv8         sv9        sv10          sv11       sv12
21 0.07962704 -0.02158386 -0.04022296  0.05123833 -0.0004505038 0.06564886
22 0.03670519  0.02245130 -0.07717329 -0.02117659 -0.0226474310 0.01271174

DxSCZ is an indicator whether the diagnosis is SCZ (1) or control (0)

genderMale is an indicator whether the gender is Male (1) or femaile (0)

Institution: a factor of three levels are coded as two indicators. Indicate 
the source of sample. 

libclust: cluster of RNA processing libraries. a batch effect. 

age, PMI

RIN, RNA integrity number that quantifies RNA quality, RIN2 = RIN*RIN

genoPCs, PCs from genotype data

log_depth, log transformed read-depth measurement, the 75 percentile of RNAseq
fragment count across all the genes of one sample. Note that the genes are 
first filtered so that each gene is expressed at least 20 in at least 25% of 
the samples in a slightly larger dataset that inlcuding these and ~50 other 
samples. 

sv1-sv12, surroage variables that aim to capture unmeasured confoudning effect. 


------------------------------------------------------------------
trec_filtered_scz_control.rds
------------------------------------------------------------------

a data matirx of gene expression, each row is one gene and each column 
is a sample. Here are the first two rows of this matrix. 

trec stands for Toteal Read Count. Though more precisely, this is 
total RNA-seq fragment count, where one fragment has two reads. 

The count were generated using R function summarizeOverlap. 

> trec0 = trec0[,match(dat$RNAseq_sample_id, colnames(trec0))]
> dim(trec0)
[1] 20788   527
> trec0[1:2,1:5]
                MSSM_RNA_PFC_1 MSSM_RNA_PFC_10 MSSM_RNA_PFC_100
ENSG00000000003            137             320              251
ENSG00000000419            416            1242              625
                MSSM_RNA_PFC_101 MSSM_RNA_PFC_102
ENSG00000000003              101              209
ENSG00000000419              366              701

------------------------------------------------------------------
trec_filtered_log_rd_corrected_scz_control.rds
------------------------------------------------------------------

log tansformed and read-depth corrected read count. 

its j-th row (gene) and i-th column (sample) is log(TReC[j,i]/q_i), 
where q_i, the read-depth measurement, is the 75 percentile of RNAseq
fragment count across all the genes of one sample. Note that the genes are 
first filtered so that each gene is expressed at least 20 in at least 25% of 
the samples in a slightly larger dataset that inlcuding these and ~50 other 
samples. 
