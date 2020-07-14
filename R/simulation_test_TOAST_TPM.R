library(TOAST)

# config = "n_100_DE_pattern_2_1_1_replicate_1"

test_TOAST_TPM = function(config) {
  RDatafile = sprintf(file.path(config, "simulation.RData"))
  load(RDatafile)
  
  # Will use: observed_TPM, rho_from_TPM, clinical_variables, adjusted_signature_gene_TPM (only for cell type names)
  colnames(rho_from_TPM) = colnames(adjusted_signature_gene_TPM)
  with_RIN = without_RIN = list()
  empty_pval_matrix = matrix(nrow=nrow(observed_TPM), ncol=ncol(rho_from_TPM))
  rownames(empty_pval_matrix) = rownames(observed_TPM)
  colnames(empty_pval_matrix) = colnames(rho_from_TPM)
  with_RIN$pval_matrix = without_RIN$pval_matrix = empty_pval_matrix
  
  # Without RIN
  cell_types = colnames(rho_from_TPM)
  design = data.frame(group=gl(2, round(n/2)))
  design_out = makeDesign(design, rho_from_TPM)
  fitted_model = fitModel(design_out, observed_TPM)
  
  for (cell_type in cell_types) {
    res_table = csTest(fitted_model, 
                       coef = "group", 
                       cell_type = cell_type)
    without_RIN$pval_matrix[, cell_type] = res_table[rownames(observed_TPM), "p_value"]
  }
  
  # With RIN
  cell_types = colnames(rho_from_TPM)
  design = data.frame(RIN=clinical_variables, group=gl(2, round(n/2)))
  design_out = makeDesign(design, rho_from_TPM)
  fitted_model = fitModel(design_out, observed_TPM)
  
  for (cell_type in cell_types) {
    res_table = csTest(fitted_model, 
                       coef = "group", 
                       cell_type = cell_type)
    with_RIN$pval_matrix[, cell_type] = res_table[rownames(observed_TPM), "p_value"]
  }
  
  pdf(file.path(config, "TOAST_TPM_pvalue_distribution.pdf"), height=8, width=10)
  par(mfrow=c(2,3))
  
  hist(without_RIN$pval_matrix[1:2000, 1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$pval_matrix[1:2000, 2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$pval_matrix[1:2000, 3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$pval_matrix[2001:10000, 1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$pval_matrix[2001:10000, 2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$pval_matrix[2001:10000, 3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  
  hist(with_RIN$pval_matrix[1:2000, 1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$pval_matrix[1:2000, 2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$pval_matrix[1:2000, 3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$pval_matrix[2001:10000, 1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$pval_matrix[2001:10000, 2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$pval_matrix[2001:10000, 3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  dev.off()
  
  save(with_RIN,
       without_RIN,
       file=file.path(config, "TOAST_TPM_res.RData"))
}
