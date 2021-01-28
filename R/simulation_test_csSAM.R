library(csSAM)
library(peakRAM)

# config = "n_100_DE_pattern_2_1_1_replicate_1"

test_csSAM = function(config) {
  RDatafile = sprintf(file.path("../simulation", config, "simulation.RData"))
  load(RDatafile)
  
  csSAM_bench = peakRAM({
    # From the demo
    fileName = file.path("../simulation", config, "csSAM_program_output.pdf")
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
  
  pdf(file.path("../simulation", config, "csSAM_fdr_distribution.pdf"), height=8, width=10)
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
       file=file.path("../simulation", config, "csSAM_res.RData"))
}