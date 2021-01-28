library(CARseq)
library(peakRAM)

################################################### 
# We test if "group" (there are two) is significant in each cell type.
# M = 2, and m = {1, 2} stands for the cell type-specific expression within a group.
# Data is simulated by "simulate_data()":
# config = "n_100_DE_pattern_2_1_1_replicate_1"
################################################### 
# Suppose we have a cell type-specific covariate called "age".
# We test if "age" is significant in each cell type.
# M = 2, and the cell type-specific intercept is m = 1 while age is m = 2 .
# Data is simulated by "simulate_data_continuous_covariates()":
# let age ~ U(30, 80) and let the fold change between cell type expression
# conditional on the two boundary ages being (2, 1, 1), etc.
# WARNING: this example is problematic; design matrix is likely to be singular
# config = "n_50_DE_pattern_4_1_1_replicate_101"
################################################### 
test_CARseq = function(config, ncpus) {
  
  RDatafile = sprintf(file.path("../simulation", config, "simulation.RData"))
  load(RDatafile)
  
  # CARseq (w/ and w/o clinical variables)
  # Input: read counts; cellular frequencies.
  
  H = 6
  n_B = n
  M = 2
  K = 1
  
  x = gl(2, round(n/2))
  
  without_RIN_bench_no_parallel = peakRAM({
    without_RIN = run_CARseq(count_matrix = observed_read_count,
                             cellular_proportions = rho,
                             groups = x,
                             formula = NULL,
                             data = clinical_variables,
                             read_depth = d,
                             shrunken_lfc = FALSE,
                             cores = 1,
                             fix_overdispersion = FALSE
    )
  })
  without_RIN_bench_parallel = peakRAM({
    without_RIN = run_CARseq(count_matrix = observed_read_count,
                             cellular_proportions = rho,
                             groups = x,
                             formula = NULL,
                             data = clinical_variables,
                             read_depth = d,
                             shrunken_lfc = FALSE,
                             cores = ncpus,
                             fix_overdispersion = FALSE
    )
  })
  without_RIN$Elapsed_Time_sec_no_parallel  = without_RIN_bench_no_parallel$Elapsed_Time_sec
  without_RIN$Peak_RAM_Used_MiB_no_parallel = without_RIN_bench_no_parallel$Peak_RAM_Used_MiB
  without_RIN$Elapsed_Time_sec_parallel  = without_RIN_bench_parallel$Elapsed_Time_sec
  # When using parallel, peakRAM is not able to get all the memory usage:
  without_RIN$Peak_RAM_Used_MiB_parallel = without_RIN_bench_parallel$Peak_RAM_Used_MiB
  without_RIN$Elapsed_Time_sec  = without_RIN_bench_parallel$Elapsed_Time_sec
  without_RIN$Peak_RAM_Used_MiB = without_RIN_bench_no_parallel$Peak_RAM_Used_MiB
  
  with_RIN_bench_no_parallel = peakRAM({
    with_RIN = run_CARseq(count_matrix = observed_read_count,
                             cellular_proportions = rho,
                             groups = x,
                             formula = ~ RIN,
                             data = clinical_variables,
                             read_depth = d,
                             shrunken_lfc = FALSE,
                             cores = 1,
                             fix_overdispersion = FALSE
    )
  })
  with_RIN_bench_parallel = peakRAM({
    with_RIN = run_CARseq(count_matrix = observed_read_count,
                             cellular_proportions = rho,
                             groups = x,
                             formula = ~ RIN,
                             data = clinical_variables,
                             read_depth = d,
                             shrunken_lfc = FALSE,
                             cores = ncpus,
                             fix_overdispersion = FALSE
    )
  })
  with_RIN$Elapsed_Time_sec_no_parallel  = with_RIN_bench_no_parallel$Elapsed_Time_sec
  with_RIN$Peak_RAM_Used_MiB_no_parallel = with_RIN_bench_no_parallel$Peak_RAM_Used_MiB
  with_RIN$Elapsed_Time_sec_parallel  = with_RIN_bench_parallel$Elapsed_Time_sec
  # When using parallel, peakRAM is not able to get all the memory usage:
  with_RIN$Peak_RAM_Used_MiB_parallel = with_RIN_bench_parallel$Peak_RAM_Used_MiB
  with_RIN$Elapsed_Time_sec  = with_RIN_bench_parallel$Elapsed_Time_sec
  with_RIN$Peak_RAM_Used_MiB = with_RIN_bench_no_parallel$Peak_RAM_Used_MiB

  save(without_RIN, with_RIN, file=file.path("../simulation", config, "CARseq_res.RData"))
  
  pdf(file.path("../simulation", config, "CARseq_pvalue_distribution.pdf"), height=8, width=10)
  par(mfrow=c(2,3))
  
  hist(without_RIN$p[1:2000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$p[1:2000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$p[1:2000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$p[2001:10000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$p[2001:10000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_RIN$p[2001:10000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  
  hist(with_RIN$p[1:2000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$p[1:2000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$p[1:2000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$p[2001:10000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$p[2001:10000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_RIN$p[2001:10000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  
  dev.off()
}
