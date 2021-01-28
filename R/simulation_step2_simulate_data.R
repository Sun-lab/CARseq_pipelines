library(CARseq)
library(MASS)
library(zetadiv)
library(multcomp)

# use mclust to fit multivariate normal distribution of
# cell type-specific expression & effect size of RIN
library(mclust)
library(Rfast)
library(Compositional)
library(openssl)

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

# six cell types in proportion estimates from ICeDT
prop = readRDS("../MTG/prop_MTG.rds")$ICeDT

# differential expression
d = exp(CMC_clinical_merged2$log_depth)
x = rep(0, length(disease))
x[disease == "SCZ"] = 1
rho = prop[match(CMC_clinical_merged2$RNAseq.Sample_RNA_ID, rownames(prop)), ]
# rho = na.omit(rho)
H = ncol(rho)
n_B = nrow(rho)
K = ncol(clinical_variables)

# Read parameter estimates from "get_distribution_of_parameters_from_real_data.R"
# obtain cooks distance matrix
estimates_list = list()
for (batch in 1:567) {
  estimates_list[[1+length(estimates_list)]] = readRDS(sprintf("../results/get_distribution_of_parameters_from_real_data_results/estimates_mat%s.rds", batch))
}
estimates = do.call(rbind, estimates_list)
estimates = estimates[seq_len(nrow(CMC_count)), ]
rownames(estimates) = rownames(CMC_count)
estimates = na.omit(estimates)
estimates[, "overdispersion"] = log(estimates[, "overdispersion"])
pdf("../figures/log_scale_expression_RIN_overdispersion.pdf")
pairs(~Intercept+RIN+overdispersion,
      as.data.frame(estimates),
      pch=19, cex=0.05, col=rgb(0,0,0,0.1))
dev.off()
table(estimates[, "overdispersion"] > 5)
# Note there is a small cluster of genes with overdispersion > 5.
# This is totally fine, but we can no longer claim that log(overdispersion) is normal.
# We will just remove the corresponding genes.
estimates = estimates[estimates[, "overdispersion"] <= 5, ]
cor(estimates)
dim(estimates)
# Notice there is almost no correlation between overdispersion and other parameters.

# Some statistics about the reference matrix derived from MTG single cell data
# We have estimated the cor_between_cell_types from single cell reference.
# Generated using _brain_cell_type/step3_DE.R
MTG_statistics = readRDS("../MTG/all_genes_MTG.rds")

multivariate_gaussian_model = mclust::mvn(modelName = "Ellipsoidal", estimates[,1:3])
multivariate_gaussian_model$parameters$mean
multivariate_gaussian_model$parameters$variance$Sigma
# Six cell types, RIN effect, overdispersion, and gene length
param_mean = c(rep(multivariate_gaussian_model$parameters$mean[1], H),
               multivariate_gaussian_model$parameters$mean[2:3],
               MTG_statistics$mean_log_gene_length)

# TODO: may need to change the order of cell types (or just 1 is minor, 2 is major)
# 9 dims for six cell types, RIN effect, overdispersion, and gene length
param_sigma = matrix(0, nrow=H+3, ncol=H+3)
param_sigma[H:(H+2), H:(H+2)] = multivariate_gaussian_model$parameters$variance$Sigma  # Cov of last cell type, RIN effect, and overdispersion
param_sigma[1:(H-1), H+2] = param_sigma[H+2, 1:(H-1)] = param_sigma[H, H+2]  # Cov of overdispersion and cell type expr
param_sigma[1:(H-1), H+1] = param_sigma[H+1, 1:(H-1)] = param_sigma[H, H+1]  # RIN effect and cell type expr
param_sigma[1:H, 1:H][upper.tri(param_sigma[1:H, 1:H])] = MTG_statistics$cor_between_cell_types * param_sigma[H, H]
param_sigma[1:H, 1:H][lower.tri(param_sigma[1:H, 1:H])] = MTG_statistics$cor_between_cell_types * param_sigma[H, H]

diag(param_sigma)[1:(H-1)] = diag(param_sigma)[H]
param_sigma[1:H, H+3] = param_sigma[H+3, 1:H] = MTG_statistics$cor_between_cell_type_and_gene_length
param_sigma[H+3, H+3] = (MTG_statistics$sd_log_gene_length)^2

# We will simulate cell type-specific expression using:
# Rfast::rmvnorm(1, param_mean, param_sigma)

# The distribution of read depth
log_read_depth_model = mclust::mvn(modelName = "X", clinical_variables_all$log_depth)
log_read_depth_mean = log_read_depth_model$parameters$mean
log_read_depth_sd = sqrt(log_read_depth_model$parameters$variance$sigmasq)
# rnorm(1, log_read_depth_mean, log_read_depth_sd)

# The distribution of cellular frequencies
estimates_diri = Rfast::diri.nr2(rho, type = 1, tol = 1e-07)
# Compositional::rdiri(1, estimates_diri$param)

# The distribution of RIN
RIN_model = mclust::mvn(modelName = "X", clinical_variables_all$RIN)
RIN_mean = RIN_model$parameters$mean
RIN_sd = sqrt(RIN_model$parameters$variance$sigmasq)
# rnorm(1, RIN_mean, RIN_sd)

# Function to generate data. This is not a standalone function as it needs the help from
# parameters that we have just fitted.
simulate_data = function(n = 100, DE_pattern = c(2, 1, 1), replicate = 1) {
  # Recast the triple "DE_pattern" (Exc, Inh, Other) into 
  # a six dimension vector of (Astro, Exc, Inh, Micro, Oligo, OPC):
  DE_pattern_full = DE_pattern[c(3, 1, 2, 3, 3, 3)]
  # Generate the RData file name
  RDatafolder = sprintf("n_%s_DE_pattern_%s_replicate_%s",
                        n,
                        paste(DE_pattern, collapse="_"),
                        replicate)
  dir.create(file.path("../simulation", RDatafolder), showWarnings = FALSE, recursive = TRUE)
  # Compute a deterministic seed from the RData file name
  # md5 is overkill but it works
  seed_from_config = strtoi(substr(openssl::md5(RDatafolder), 1, 7), 16L)
  set.seed(seed_from_config)
  
  # Among all the genes, without loss of generality, suppose DE genes come at the top:
  prop_of_DE_genes = 0.2  # proportion of differentially expressed genes
  G = 10000  # total number of genes
  G_signature = 100  # total number of signature genes per cell type
  
  # generate cell type-specific expression
  # Three cell types, RIN effect, overdispersion, and gene length
  cell_type_specific_expression = Rfast::rmvnorm(G, param_mean, param_sigma)
  cell_type_specific_expression[, c(1:H, (H+2):(H+3))] = exp(cell_type_specific_expression[, c(1:H, (H+2):(H+3))])
  colnames(cell_type_specific_expression) = c(colnames(rho), "RIN", "overdispersion", "gene_lengths")
  clinical_variables_effect_size = cell_type_specific_expression[,"RIN", drop=FALSE]
  
  
  # generate clinical variables (RIN)
  clinical_variables = matrix(rnorm(n, RIN_mean, RIN_sd), ncol=1)
  colnames(clinical_variables) = "RIN"
  
  # pick the signature genes by computing log fold change from 
  # cell type-specific mean among the genes without differential expression
  fold_change_among_non_DE_genes = matrix(0, nrow=G, ncol=H)
  for (h in 1:H) {
    fold_change_among_non_DE_genes[, h] = cell_type_specific_expression[, h] / rowSums(cell_type_specific_expression[, 1:H][, -h])
  }
  signature_gene_list = list()
  for (h in seq_len(ncol(fold_change_among_non_DE_genes))) {
    # need to exclude the DE genes at the top of the table
    signature_gene_list[[h]] = head(setdiff(order(fold_change_among_non_DE_genes[, h], decreasing = TRUE),
                                            1:round(prop_of_DE_genes * G)),
                                    n = G_signature)
  }
  signature_gene_unique_indices = unlist(signature_gene_list)
  signature_gene_unique_indices = signature_gene_unique_indices[!duplicated(signature_gene_unique_indices)]
  
  # generate gene lengths
  # To make the gene lengths more realistic, they should be rounded.
  # However, not rounding will not make much difference in the simulation since the gene lengths are quite large. 
  gene_lengths = round(cell_type_specific_expression[, "gene_lengths"])
  
  # generate read depth
  d = round(exp(rnorm(n, log_read_depth_mean, log_read_depth_sd)))
  
  # generate cellular frequencies
  rho = Compositional::rdiri(n, estimates_diri$param)
  
  # generate overdispersion
  overdispersion = cell_type_specific_expression[, "overdispersion"]
  
  # among the DE genes, assume half of them are higher in group 1 and half are higher in group 2:
  expected_mean_expression = array(0, dim=c(G, n, 2))
  gene_indices = 1:round(0.5 * prop_of_DE_genes * G)
  expected_mean_expression[gene_indices, , 1] = 
    cell_type_specific_expression[gene_indices, 1:H] %*%
    (t(rho) * DE_pattern_full) %*%
    diag(d) *
    exp(tcrossprod(clinical_variables_effect_size[gene_indices,], clinical_variables))
  expected_mean_expression[gene_indices, , 2] = 
    cell_type_specific_expression[gene_indices, 1:H] %*%
    t(rho) %*%
    diag(d) *
    exp(tcrossprod(clinical_variables_effect_size[gene_indices,], clinical_variables))
  gene_indices = (1+round(0.5 * prop_of_DE_genes * G)):round(prop_of_DE_genes * G)
  expected_mean_expression[gene_indices, , 1] = 
    cell_type_specific_expression[gene_indices, 1:H] %*%
    t(rho) %*%
    diag(d) *
    exp(tcrossprod(clinical_variables_effect_size[gene_indices,], clinical_variables))
  expected_mean_expression[gene_indices, , 2] =   
    cell_type_specific_expression[gene_indices, 1:H] %*%
    (t(rho) * DE_pattern_full) %*%
    diag(d) *
    exp(tcrossprod(clinical_variables_effect_size[gene_indices,], clinical_variables))
  gene_indices = (1+round(prop_of_DE_genes * G)):G
  expected_mean_expression[gene_indices, , 1] =  
      expected_mean_expression[gene_indices, , 2] = 
    cell_type_specific_expression[gene_indices, 1:H] %*%
    t(rho) %*%
    diag(d) *
    exp(tcrossprod(clinical_variables_effect_size[gene_indices,], clinical_variables))
  
  # random generation for the negative binomial distribution
  observed_read_count = matrix(NA, nrow=G, ncol=n)
  # group 1
  sample_indices = 1:round(0.5*n)
  observed_read_count[,sample_indices] = rnbinom(n = G*length(sample_indices),
                                                 mu = expected_mean_expression[,sample_indices,1],
                                                 size = overdispersion)
  # group 2
  sample_indices = (1+round(0.5*n)):n
  observed_read_count[,sample_indices] = rnbinom(n = G*length(sample_indices),
                                                     mu = expected_mean_expression[,sample_indices,2],
                                                     size = overdispersion)
  observed_read_count_in_signature_genes = observed_read_count[signature_gene_unique_indices, ]
  
  # # generate gene lengths
  # # suppose a high corr between log gene lengths and 75% percentile of log expression
  # # gene length 
  # log_expression = log(cell_type_specific_expression[, 1:H] %*% t(rho))
  # scaled_log_expression_third_quartile = scale(apply(log_expression, 1, function(x) quantile(x, 0.75)))
  # gene_lengths = as.numeric(round(1000*exp(scaled_log_expression_third_quartile * 0.7 + rnorm(G) * (1 - 0.7))))
  
  # cell size
  cell_sizes = colSums(cell_type_specific_expression[, 1:H] / gene_lengths)
  cell_sizes
  
  # cellular frequencies derived from TPM
  # adjust cellular frequencies by cell size
  rho_from_TPM = rho * rep(cell_sizes, times=n)
  rho_from_TPM = rho_from_TPM / rowSums(rho_from_TPM)
  
  # TPM of all genes
  # michael's version
  # https://support.bioconductor.org/p/91218/
  tpm3 = function(counts,len) {
    x = counts/len
    return(t(t(x)*1e6/colSums(x)))
  }
  observed_TPM = tpm3(observed_read_count, gene_lengths)
  
  # TPM of signature matrix
  signature_gene_TPM = tpm3(cell_type_specific_expression[, 1:H], gene_lengths)[signature_gene_unique_indices, ]
  
  # TPM of signature matrix after adjustment by considering the batch effect RIN.
  # The adjustment is needed because there is no "centering" after the batch effect RIN is added.
  # After the adjustment, the mixture TPM will be roughly weighted sum of reference cell type-specific TPM.
  # We multiply the reference cell type-specific expression by mean of effects stemming from batch effect RIN.
  # Note that we did not consider that there is differential expression between case/control samples.
  adjusted_cell_type_specific_expression = cell_type_specific_expression[,1:H] * colMeans(exp(tcrossprod(clinical_variables, clinical_variables_effect_size)))
  # The result is similar to:
  # adjusted_cell_type_specific_expression = cell_type_specific_expression[,1:H] * as.numeric(exp(mean(clinical_variables) * clinical_variables_effect_size))
  adjusted_signature_gene_TPM = tpm3(adjusted_cell_type_specific_expression, gene_lengths)[signature_gene_unique_indices, 1:H]
  # Well, we actually only wanted cell type-specific expression (the first H columns) in cell_type_specific_expression:
  cell_type_specific_expression = cell_type_specific_expression[, 1:H]
  
  # Add gene ID
  rownames(observed_read_count) = rownames(observed_TPM) = paste0("gene", 1:10000)
  rownames(signature_gene_TPM) =  rownames(adjusted_signature_gene_TPM) = paste0("gene", signature_gene_unique_indices)
  
  # prepare CIBERSORTx input
  geneSymbol_and_observed_TPM = cbind(rownames(observed_TPM), observed_TPM)
  colnames(geneSymbol_and_observed_TPM) = c("geneSymbol", paste0("sample", seq_len(n)))
  geneSymbol_and_signature_gene_TPM = cbind(rownames(adjusted_signature_gene_TPM), adjusted_signature_gene_TPM)
  colnames(geneSymbol_and_signature_gene_TPM) = c("geneSymbol", colnames(adjusted_signature_gene_TPM))
  # generate matrices as input for CIBERSORTx online
  write.table(geneSymbol_and_observed_TPM,
              file = file.path("../simulation", RDatafolder, sprintf("CIBERSORTx_input_observed_TPM_%s.txt", RDatafolder)),
              quote = FALSE,
              row.names = FALSE,
              sep="\t")
  write.table(geneSymbol_and_signature_gene_TPM,
              file = file.path("../simulation", RDatafolder, sprintf("CIBERSORTx_input_signature_gene_TPM_%s.txt", RDatafolder)),
              quote = FALSE,
              row.names = FALSE,
              sep="\t")
  write(paste0("gene", 1:10000),
        file = file.path("../simulation", RDatafolder, sprintf("CIBERSORTx_input_gene_subset_1_10000_%s.txt", RDatafolder)))
  
  
  # save them into RData
  fold_change = max(DE_pattern)  # probably not used
  save(list = c(
    "fold_change",
    "n",
    "rho",
    "d",
    "overdispersion",
    "gene_lengths",
    "signature_gene_unique_indices",
    "cell_type_specific_expression",
    "adjusted_cell_type_specific_expression",
    "observed_read_count",
    "observed_read_count_in_signature_genes",
    "clinical_variables",
    "clinical_variables_effect_size",
    "observed_TPM",
    "signature_gene_TPM",
    "adjusted_signature_gene_TPM",
    "rho_from_TPM"
    ),
    file = file.path("../simulation", RDatafolder, "simulation.RData")
  )
  # return the configuration
  RDatafolder
}

# NOTE:
# Only simulation of groups "simulate_data" is included in this script.
# The function relevant to "simulate_data_continuous_covariates" is not included.

# setup 1
n_list = c(50, 100, 200)  # total number of samples, 50, 100, 200
# add cell type-specific DE and compute the expected mean
fold_changes = c(2, 3, 4)  # 2,3,4
replicates = 1:10
for (n in n_list) {
  for (fold_change in fold_changes) {
    DE_pattern = c(fold_change, 1, 1)
    for (replicate in replicates) {
      simulate_data(n = n, DE_pattern = DE_pattern, replicate = replicate)
    }
  }
}


# setup 2
n_list = c(50, 100, 200)  # total number of samples, 50, 100, 200
# add cell type-specific DE and compute the expected mean
fold_changes = c(2, 3, 4)  # 2,3,4
replicates = 1:10
for (n in n_list) {
  for (fold_change in fold_changes) {
    DE_pattern = c(fold_change, fold_change, 1)
    for (replicate in replicates) {
      simulate_data(n = n, DE_pattern = DE_pattern, replicate = replicate)
    }
  }
}


# setup 3
n_list = c(50, 100, 200)  # total number of samples, 50, 100, 200
# add cell type-specific DE and compute the expected mean
fold_changes = c(2, 3, 4)  # 2,3,4
replicates = 1:10
for (n in n_list) {
  for (fold_change in fold_changes) {
    DE_pattern = c(fold_change, round(1000/fold_change)/1000, 1)
    for (replicate in replicates) {
      simulate_data(n = n, DE_pattern = DE_pattern, replicate = replicate)
    }
  }
}


# setup 4
n_list = c(50, 100, 200)  # total number of samples, 50, 100, 200
# add cell type-specific DE and compute the expected mean
fold_changes = c(2, 3, 4)  # 2,3,4
replicates = 1:10
for (n in n_list) {
  for (fold_change in fold_changes) {
    DE_pattern = c(1, fold_change, 1)
    for (replicate in replicates) {
      simulate_data(n = n, DE_pattern = DE_pattern, replicate = replicate)
    }
  }
}
