args = commandArgs(trailingOnly=TRUE)

if (length(args) == 1) {
  i = as.numeric(args[1])
}

configs = list.files(path="../simulation", pattern="^n_.*_DE_pattern_.*_replicate_[0-9]?[0-9]?$", full.names = TRUE) # all 10 replicates

config = configs[i]

# config = "n_100_DE_pattern_2_1_1_replicate_1"

source("simulation_test_csSAM.R")

source("simulation_test_csSAM_ICeDT.R")

source("simulation_test_CARseq.R")

source("simulation_test_CARseq_ICeDT.R")

source("simulation_test_TOAST_TPM.R")

source("simulation_test_TOAST_TPM_ICeDT.R")

source("simulation_test_TOAST_count.R")

source("simulation_test_TOAST_count_ICeDT.R")


#for (config in configs) {

  if (!file.exists(file.path(config, "csSAM_res.RData"))) {
    test_csSAM(config)
  }

  if (!file.exists(file.path(config, "csSAM_ICeDT_res.RData"))) {
    test_csSAM_ICeDT(config)
  }

  # DE test using TOAST with expression in TPM as response (the setting shown in main figures)
  
  if (!file.exists(file.path(config, "TOAST_TPM_res.RData"))) {
    test_TOAST_TPM(config)
  }

  if (!file.exists(file.path(config, "TOAST_TPM_ICeDT_res.RData"))) {
    test_TOAST_TPM_ICeDT(config)
  }

  if (!file.exists(file.path(config, "CARseq_res.RData"))) {
    test_CARseq(config, ncpus = 12)
  }

  if (!file.exists(file.path(config, "CARseq_ICeDT_res.RData"))) {
    test_CARseq_ICeDT(config, ncpus = 12)
  }

  
  ## These tests are experimental, only done for one (1) replicate, and are not covered in main text:
  #if (grepl("config_1$", config)) {
    # DE test using TOAST with expression in counts as response
    if (!file.exists(file.path(config, "TOAST_res.RData"))) {
      test_TOAST(config)
    }
    
    if (!file.exists(file.path(config, "TOAST_ICeDT_res.RData"))) {
      test_TOAST_ICeDT(config)
    }

    # This requires additional CIBERSORT output files to run.
    # The computational work is done by https://cibersortx.stanford.edu/runcibersortx.php
    # We submit high-resolution w/o B mode correction as Job 1.
    # We submit high-resolution w/ B mode correction as Job 2.
    # Then unzip both of them under the $config/ folder.
    # Otherwise, there is no guarantee against the error messages you will see.
    # if (!file.exists(file.path(config, "CIBERSORTx_res.RData"))) {
    #   test_CIBERSORTx(config)
    # }
  #}
# }