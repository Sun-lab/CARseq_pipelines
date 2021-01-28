library(csSAM)
library(peakRAM)

# config = "n_100_DE_pattern_2_1_1_replicate_1"

test_csSAM_ICeDT = function(config) {
  RDatafile = sprintf(file.path("../simulation", config, "simulation.RData"))
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
  RMSE_rho = apply(rho - icedt_rho, 2, function(x) sqrt(mean(x^2)))
  RMSE_rho_from_TPM = apply(rho_from_TPM - icedt_rho_from_TPM, 2, function(x) sqrt(mean(x^2)))
  rho_summary = list(rho, icedt_rho, RMSE_rho, rho_from_TPM, icedt_rho_from_TPM, RMSE_rho_from_TPM, cell_sizes)
  
  csSAM_bench = peakRAM({
    # From the demo
    fileName = file.path("../simulation", config, "csSAM_ICeDT_program_output.pdf")
    # Now run, either using the wrapper
    # NB: more permutations would be needed for real data
    deconvResults = csSamWrapper(t(observed_TPM), rho_from_TPM, gl(2, round(n/2)), nperms = 50, alternative = "two.sided"
                                 , nonNeg = TRUE
                                 , standardize = TRUE
                                 , medianCenter = TRUE
                                 , fileName = fileName)
    fdr_matrix = t(deconvResults$sigGene.csSAM)
    # str(fdr_matrix)
    fdr_matrix
  })
  Elapsed_Time_sec  = csSAM_bench$Elapsed_Time_sec
  Peak_RAM_Used_MiB = csSAM_bench$Peak_RAM_Used_MiB
  
  pdf(file.path("../simulation", config, "csSAM_ICeDT_fdr_distribution.pdf"), height=8, width=10)
  par(mfrow=c(2,3))
  hist(fdr_matrix[1:2000, 1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(fdr_matrix[1:2000, 2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(fdr_matrix[1:2000, 3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(fdr_matrix[2001:10000, 1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(fdr_matrix[2001:10000, 2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(fdr_matrix[2001:10000, 3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  dev.off()
  
  save(deconvResults,
       fdr_matrix,
       Elapsed_Time_sec,
       Peak_RAM_Used_MiB,
       Elapsed_Time_sec,
       rho_summary,
       file=file.path("../simulation", config, "csSAM_ICeDT_res.RData"))
}