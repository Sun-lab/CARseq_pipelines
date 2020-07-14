args = commandArgs(trailingOnly=TRUE)

if (length(args) == 1) {
  batch = as.numeric(args[1])
}
batchsize = 100
# input

library(CARseq)
library(MASS)
library(zetadiv)
library(multcomp)

#Read counts:
# CMC data in counts -- CARseq needs raw counts
CMC_count = readRDS("../data/CMC_MSSM-Penn-Pitt_Paul_geneExpressionRaw.rds")$so1
dim(CMC_count)
CMC_count[1:5,1:5]
# CMC label
CMC_clinical = read.csv("../data/CMC-CMC_HBCC_clinical_.csv") 
dim(CMC_clinical)
names(CMC_clinical)
CMC_key = read.csv("../data/Release3_SampleID_key.csv")
dim(CMC_key)
names(CMC_key)
CMC_key["Individual_ID","RNAseq.Sample_RNA_ID"]
CMC_clinical_merged = merge(CMC_clinical, CMC_key, by.x="Individual.ID", by.y="Individual_ID")
CMC_clinical_merged = CMC_clinical_merged[match(colnames(CMC_count), CMC_clinical_merged$RNAseq.Sample_RNA_ID), ]
CMC_clinical_merged = CMC_clinical_merged[CMC_clinical_merged$Dx %in% c("Control","SCZ"), ]
# read clinical variables
clinical_Control = readRDS("../data/trec_dat_post_QC_PCA_Control.rds")
clinical_SCZ = readRDS("../data/trec_dat_post_QC_PCA_SCZ.rds")
clinical_variables_all = rbind(clinical_Control$DATA, clinical_SCZ$DATA)
CMC_clinical_merged2 = merge(CMC_clinical_merged, clinical_variables_all, by.x="RNAseq.Sample_RNA_ID", by.y="RNAseq_sample_id",
                             all = FALSE, suffixes = c(".x", ""))
clinical_variables = model.matrix( ~ 1 + RIN,
                                   data = CMC_clinical_merged2
)[, -1, drop=FALSE]


disease = CMC_clinical_merged2$Dx
table(disease, useNA="ifany")
CMC_count = CMC_count[, match(CMC_clinical_merged2$RNAseq.Sample_RNA_ID, colnames(CMC_count))]

stopifnot(as.numeric(R.version$minor) >= 6.0)  # R < 3.6.0 has a different RNG
set.seed(1234)

# proportion estimates from ICeDT
prop = as.data.frame(readRDS("../MTG/prop_MTG.rds")$ICeDT)
prop$Other = rowSums(prop[, c("Astro", "Micro", "Oligo", "OPC")])
prop = as.matrix(prop[, c("Other", "Exc", "Inh")])

# differential expression
d = exp(CMC_clinical_merged2$log_depth)
x = rep(0, length(disease))
x[disease == "SCZ"] = 1
rho = prop[match(CMC_clinical_merged2$RNAseq.Sample_RNA_ID, rownames(prop)), ]
# rho = na.omit(rho)
H = ncol(rho)
n_B = nrow(rho)
K = ncol(clinical_variables)

# only use about 15000 ~ 20000 highly expressed genes:
rd75 = apply(CMC_count, 1, function(x) quantile(x, 0.75))

# specify design matrix
cell_type_specific_variables_simulation = array(1, dim=c(n_B, H, 1))

# H cell types, 1 clinical variable, 1 overdispersion parameter,
number_of_parameters = 3
estimates_mat = matrix(nrow=batchsize, ncol=number_of_parameters)
# Start hypothesis testing
for (j in (1+(batch-1)*batchsize):(batch*batchsize)) {
  
  # simulate data (simulate t_{ji}):
  if (j > nrow(CMC_count)) break
  
  if (rd75[j] <= 20) next
  
  estimates_mat[j - (batch-1)*batchsize, ] = tryCatch({
    # filter samples by Cook's distance
    glmmodel = MASS::glm.nb(CMC_count[j, ]~offset(log(d))+clinical_variables)
    c(glmmodel$coefficients, glmmodel$theta)
    
  }, error = function(e) {
    NA
  }, finally = {
  })
  
  ##################################################
  # end testing one cell type using nloptr
  ##################################################
}

colnames(estimates_mat) = c("Intercept", "RIN", "overdispersion")
saveRDS(estimates_mat, sprintf("../results/get_distribution_of_parameters_from_real_data_results/estimates_mat%s.rds", batch))
