args = commandArgs(trailingOnly=TRUE)

if (length(args) == 1) {
  i = as.numeric(args[1])
}

configs = list.files(path="../simulation", pattern="^n_.*_DE_pattern_.*_1_1_replicate_.?[0123456789]$", full.names=TRUE) # d_1_1 from replicate 1 to 99
# configs = list.files(pattern="^n_.*_DE_pattern_.*_replicate_.*")             # replicates 1_1_d
# configs = list.files(pattern="n_.*_DE_pattern_.*_replicate_.?[0123456789]$") # from replicate 1 to 99 
# configs = list.files(pattern="n_.*_DE_pattern_.*_replicate_.[023456789]$")   # from replicate 2 to 10
# configs = list.files(pattern="n_.*_DE_pattern_.*_replicate_.*101$")          # from replicate 100

config = configs[i]

# config = "n_100_DE_pattern_2_1_1_replicate_1"

source("simulation_test_csSAM.R")

source("simulation_test_CARseq.R")

source("simulation_test_TOAST_TPM.R")

to_logit = function(x) {log(x / (1 - x))}
to_inverse_logit = function(x) {1/(1 + exp(-x))}

# for (config in configs) {
  
  input_data_names = load(file.path(config, "simulation.RData"))
  rho_true = rho
  cell_sizes = colSums(cell_type_specific_expression / gene_lengths)
  cell_sizes
  config_strsplit_list = unlist(strsplit(config, split = "_replicate_"))
  
  # config_noise: we add noise to cell fractions. The noise follows N(mean=0, sd=0.01) on logit scale and
  # we re-scale cellular proportions within a sample so that the sum is 1.
  # The replicate number starts with 201:
  # config_noise = "n_100_DE_pattern_2_1_1_replicate_201"
  config_noise = sprintf("%s_replicate_%d", config_strsplit_list[1], as.integer(config_strsplit_list[2])+200L)
  dir.create(config_noise, recursive = TRUE)
  set.seed(1234)
  rho_with_noise = rho_true
  rho_with_noise[] = pmax(0, to_inverse_logit(to_logit(rho_true) + rnorm(n = length(rho_true), mean = 0, sd = 0.1)))
  rho_with_noise = t(t(rho_with_noise) / rowSums(rho_with_noise))
  rho = rho_with_noise
  rho_from_TPM = rho * rep(cell_sizes, each=nrow(rho))
  rho_from_TPM = rho_from_TPM / rowSums(rho_from_TPM)
  save(list = input_data_names, file = file.path(config_noise, "simulation.RData"))
  
  if (!file.exists(file.path(config_noise, "csSAM_res.RData"))) {
    test_csSAM(config_noise)
  }

  if (!file.exists(file.path(config_noise, "TOAST_TPM_res.RData"))) {
    test_TOAST_TPM(config_noise)
  }
  
  # This requires a lot of computing resource:
  # with 12 CPUs, (100 samples, 10000 genes, 3 cell types, 1 batch effect) takes 10 minutes.
  # Since we fit both models w/ and w/o batch effect, it takes roughly 20 minutes.
  set.seed(1234)
  if (!file.exists(file.path(config_noise, "CARseq_res.RData"))) {
    test_CARseq(config_noise, ncpus = 1)
  }
  
  input_data_names = load(file.path(config, "simulation.RData"))
  rho_true = rho
  cell_sizes_true = colSums(cell_type_specific_expression / gene_lengths)
  cell_sizes_true
  config_strsplit_list = unlist(strsplit(config, split = "_replicate_"))
  
  # config_size_factor: we add a cell size factor of c(2, 1, 1) to the cellular frequencies. 
  # (They should be around c(1, 1, 1) when the expression data was simulated.)
  # The replicate number starts with 301:
  # config_size_factor = "n_100_DE_pattern_2_1_1_replicate_301"
  config_size_factor = sprintf("%s_replicate_%d", config_strsplit_list[1], as.integer(config_strsplit_list[2])+300L)
  dir.create(config_size_factor, recursive = TRUE)
  cell_sizes = cell_sizes_true
  cell_sizes[2] = cell_sizes[2] * 2.0
  rho_incorrect_cell_size = rho_from_TPM / rep(cell_sizes, each=nrow(rho_from_TPM))
  rho_incorrect_cell_size = rho_incorrect_cell_size / rowSums(rho_incorrect_cell_size)
  rho = rho_incorrect_cell_size
  save(list = input_data_names, file = file.path(config_size_factor, "simulation.RData"))
  
  # TOAST and csSAM use rho_TPM, so cell size factor is not involved. The following results should be exactly same
  # as the replicate where cell size factor is not misspecified.
  if (!file.exists(file.path(config_size_factor, "csSAM_res.RData"))) {
    test_csSAM(config_size_factor)
  }

  if (!file.exists(file.path(config_size_factor, "TOAST_TPM_res.RData"))) {
    test_TOAST_TPM(config_size_factor)
  }

  # This requires a lot of computing resource:
  # with 12 CPUs, (100 samples, 10000 genes, 3 cell types, 1 batch effect) takes 10 minutes.
  # Since we fit both models w/ and w/o batch effect, it takes roughly 20 minutes.
  set.seed(1234)
  if (!file.exists(file.path(config_size_factor, "CARseq_res.RData"))) {
    test_CARseq(config_size_factor, ncpus = 1)
  }
  
  # config_size_factor: we add a cell size factor of c(1.2, 1, 1) to the cellular frequencies. 
  # (They should be c(1, 1, 1) when the expression data was simulated.)
  # The replicate number starts with 401:
  # config_size_factor = "n_100_DE_pattern_2_1_1_replicate_401"
  config_size_factor = sprintf("%s_replicate_%d", config_strsplit_list[1], as.integer(config_strsplit_list[2])+400L)
  dir.create(config_size_factor, recursive = TRUE)
  cell_sizes = cell_sizes_true
  cell_sizes[2] = cell_sizes[2] * 1.2
  rho_incorrect_cell_size = rho_from_TPM / rep(cell_sizes, each=nrow(rho_from_TPM))
  rho_incorrect_cell_size = rho_incorrect_cell_size / rowSums(rho_incorrect_cell_size)
  rho = rho_incorrect_cell_size
  save(list = input_data_names, file = file.path(config_size_factor, "simulation.RData"))
  
  # TOAST and csSAM use rho_TPM, so cell size factor is not involved. The following results should be exactly same
  # as the replicate where cell size factor is not misspecified.
  if (!file.exists(file.path(config_size_factor, "csSAM_res.RData"))) {
    test_csSAM(config_size_factor)
  }
  
  if (!file.exists(file.path(config_size_factor, "TOAST_TPM_res.RData"))) {
    test_TOAST_TPM(config_size_factor)
  }
  
  # This requires a lot of computing resource:
  # with 12 CPUs, (100 samples, 10000 genes, 3 cell types, 1 batch effect) takes 10 minutes.
  # Since we fit both models w/ and w/o batch effect, it takes roughly 20 minutes.
  set.seed(1234)
  if (!file.exists(file.path(config_size_factor, "CARseq_res.RData"))) {
    test_CARseq(config_size_factor, ncpus = 1)
  }
  
# }
