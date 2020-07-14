args = commandArgs(trailingOnly=TRUE)

if (length(args) == 1) {
  i = as.numeric(args[1])
}

configs = list.files(pattern="../data/^n_.*_DE_pattern_.*_replicate_.*", full.names = TRUE) # all 10 replicates

config = configs[i]

# config = "n_100_DE_pattern_2_1_1_replicate_1"

source("simulation_test_csSAM.R")

source("simulation_test_CARseq.R")

source("simulation_test_CIBERSORTx.R")

source("simulation_test_TOAST_TPM.R")

source("simulation_test_TOAST_count.R")

for (config in configs) {

  if (!file.exists(file.path(config, "csSAM_res.RData"))) {
    test_csSAM(config)
  }

  # DE test using TOAST with expression in TPM as response (the setting shown in main figures)
  if (!file.exists(file.path(config, "TOAST_res_TPM.RData"))) {
    test_TOAST_TPM(config)
  }

  if (!file.exists(file.path(config, "CARseq_res.RData"))) {
    test_CARseq(config, ncpus = 12)
  }
  
  # These tests are experimental, only done for one (1) replicate, and are not covered in main text:
  if (grepl("config_1$", config)) {
    # DE test using TOAST with expression in counts as response
    if (!file.exists(file.path(config, "TOAST_res.RData"))) {
      test_TOAST(config)
    }
    
    # This requires additional CIBERSORT output files to run.
    # The computational work is done by https://cibersortx.stanford.edu/runcibersortx.php
    # We submit high-resolution w/o B mode correction as Job 1.
    # We submit high-resolution w/ B mode correction as Job 2.
    # Then unzip both of them under the $config/ folder.
    # Otherwise, there is no guarantee against the error messages you will see.
    if (!file.exists(file.path(config, "CIBERSORTx_res.RData"))) {
      test_CIBERSORTx(config)
    }
  }
}