library(TOAST)

# config = "n_100_DE_pattern_2_1_1_replicate_1"

test_TOAST_ICeDT = function(config) {
  RDatafile = sprintf(file.path(config, "simulation.RData"))
  load(RDatafile)
  
  # run ICeDT
  icedt_output_file =  file.path("../simulation", config, "ICeDT_output.rds")
  if (!file.exists(icedt_output_file)) {
    set.seed(1234)
    icedt_output = ICeDT::ICeDT(
      Y = observed_TPM[rownames(adjusted_signature_gene_TPM), ],
      Z = adjusted_signature_gene_TPM,
      tumorPurity = rep(0, ncol(observed_TPM)),
      refVar = NULL)
    # save to ICeDT cellular frequency RDS file
    saveRDS(icedt_output, icedt_output_file)
  } else {
    icedt_output = readRDS(icedt_output_file)
  }
  cell_sizes = colSums(cell_type_specific_expression / gene_lengths)
  icedt_rho_from_TPM = t(icedt_output$rho)[,c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")]
  icedt_rho = t(apply(t(icedt_output$rho)[,c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")],1,function(xx){yy = xx / cell_sizes; yy / sum(yy)}))
  
  # Will use: observed_read_count, rho, clinical_variables, signature_gene_TPM (only for cell type names)
  cell_types = colnames(icedt_rho) = colnames(signature_gene_TPM)
  with_RIN = without_RIN = list()
  empty_pval_matrix = matrix(nrow=nrow(observed_read_count), ncol=ncol(icedt_rho))
  rownames(empty_pval_matrix) = rownames(observed_read_count)
  colnames(empty_pval_matrix) = colnames(icedt_rho)
  with_RIN$pval_matrix = without_RIN$pval_matrix = empty_pval_matrix
  
  # Without RIN
  design = data.frame(group=gl(2, round(n/2)), depth=d)
  design_out = makeDesign(design, icedt_rho)
  fitted_model = fitModel(design_out, observed_read_count)
  
  for (cell_type in cell_types) {
    res_table = csTest(fitted_model, 
                       coef = "group", 
                       cell_type = cell_type)
    without_RIN$pval_matrix[, cell_type] = res_table[rownames(observed_read_count), "p_value"]
  }
  
  # With RIN
  design = data.frame(RIN=clinical_variables, group=gl(2, round(n/2)), depth=d)
  design_out = makeDesign(design, icedt_rho)
  fitted_model = fitModel(design_out, observed_read_count)
  
  for (cell_type in cell_types) {
    res_table = csTest(fitted_model, 
                       coef = "group", 
                       cell_type = cell_type)
    with_RIN$pval_matrix[, cell_type] = res_table[rownames(observed_read_count), "p_value"]
  }
  
  pdf(file.path(config, "TOAST_ICeDT_pvalue_distribution.pdf"), height=8, width=10)
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
       file=file.path(config, "TOAST_ICeDT_res.RData"))
}
