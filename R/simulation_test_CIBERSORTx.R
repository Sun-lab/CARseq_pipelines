# Submit jobs on CIBERSORTx.
# The truth in pattern_2_1_1: the first cell type is differentially expressed across two groups for genes {1..2000}.
# Submit jobs using "Impute Cell Expression" and "High-Resolution"
# Using no batch correction & B-mode batch correction
# read and analyze CIBERSORTx results

# config = "n_100_DE_pattern_2_1_1_replicate_1"

test_CIBERSORTx = function(config) {
  RDatafile = sprintf(file.path(config, "simulation.RData"))
  load(RDatafile)
  settings = grep("CIBERSORTx_Job.*_output", list.dirs(path = config, full.names = FALSE), value = TRUE)
  if (length(settings) != 2) {
    message(sprintf("CIBERSORTx output files not found under %s; skipping the DE test using CIBERSORTx", config))
    return()
  }
  jobnames = gsub("CIBERSORTx_|_output", "", settings)
  
  celltypes = c("Exc", "Inh", "Other")
  
  # read high resolution expression into a nested list
  without_B_mode = with_B_mode = list()
  without_B_mode_without_RIN = with_B_mode_without_RIN = list()
  without_B_mode$high_resolution = with_B_mode$high_resolution = list()
  for (j in seq_along(celltypes)) {
    celltype = celltypes[j]
    without_B_mode$high_resolution[[j]] = read.table(file.path(config, settings[[1]], sprintf("CIBERSORTxHiRes_%s_%s_Window12.txt", jobnames[1], celltype)),
                                         header=TRUE, row.names=1)
    with_B_mode$high_resolution[[j]] = read.table(file.path(config, settings[[2]], sprintf("CIBERSORTxHiRes_%s_%s_Window12.txt", jobnames[2], celltype)),
                                         header=TRUE, row.names=1)
  }
  
  # test
  group_labels = gl(2, round(n/2))
  without_B_mode$pval_matrix = with_B_mode$pval_matrix = 
          without_B_mode_without_RIN$pval_matrix = with_B_mode_without_RIN$pval_matrix = 
    matrix(NA, nrow=nrow(without_B_mode$high_resolution[[1]]), ncol=length(celltypes))
  colnames(without_B_mode$pval_matrix) = colnames(with_B_mode$pval_matrix) = 
          colnames(without_B_mode_without_RIN$pval_matrix) = colnames(with_B_mode_without_RIN$pval_matrix) = 
    celltypes
  
  for (gene_ind in seq_len(nrow(without_B_mode$high_resolution[[1]]))) {
    for (j in seq_along(celltypes)) {
      sd_gene_expression = sd(without_B_mode$high_resolution[[j]][gene_ind, ])
      if (is.na(sd_gene_expression)) {
        next
      } else if (sd_gene_expression == 0) {
        without_B_mode$pval_matrix[gene_ind, j] = 1
        without_B_mode_without_RIN$pval_matrix[gene_ind, j] = 1
      } else {
        without_B_mode$pval_matrix[gene_ind, j] = 
          summary(lm(as.numeric(without_B_mode$high_resolution[[j]][gene_ind, ]) ~ group_labels + clinical_variables))$coefficients["group_labels2","Pr(>|t|)"]
        without_B_mode_without_RIN$pval_matrix[gene_ind, j] = 
          summary(lm(as.numeric(without_B_mode$high_resolution[[j]][gene_ind, ]) ~ group_labels))$coefficients["group_labels2","Pr(>|t|)"]
      }
  
      sd_gene_expression = sd(with_B_mode$high_resolution[[j]][gene_ind, ])
      if (is.na(sd_gene_expression)) {
        next
      } else if (sd_gene_expression == 0) {
        with_B_mode$pval_matrix[gene_ind, j] = 1
        with_B_mode_without_RIN$pval_matrix[gene_ind, j] = 1
      } else {
        with_B_mode$pval_matrix[gene_ind, j] = 
          summary(lm(as.numeric(with_B_mode$high_resolution[[j]][gene_ind, ]) ~ group_labels + clinical_variables))$coefficients["group_labels2","Pr(>|t|)"]
        with_B_mode_without_RIN$pval_matrix[gene_ind, j] = 
          summary(lm(as.numeric(with_B_mode$high_resolution[[j]][gene_ind, ]) ~ group_labels))$coefficients["group_labels2","Pr(>|t|)"]
      }
    }
  }
  
  save(without_B_mode, with_B_mode,
       without_B_mode_without_RIN, with_B_mode_without_RIN,
       file=file.path(config, "CIBERSORTx_res.RData"))
  
  pdf(file.path(config, "CIBERSORTx_pvalue_distribution.pdf"), height=8, width=10)
  par(mfrow=c(2,3))
  
  hist(without_B_mode$pval_matrix[1:2000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode$pval_matrix[1:2000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode$pval_matrix[1:2000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode$pval_matrix[2001:10000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode$pval_matrix[2001:10000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode$pval_matrix[2001:10000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  
  hist(without_B_mode_without_RIN$pval_matrix[1:2000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode_without_RIN$pval_matrix[1:2000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode_without_RIN$pval_matrix[1:2000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode_without_RIN$pval_matrix[2001:10000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode_without_RIN$pval_matrix[2001:10000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(without_B_mode_without_RIN$pval_matrix[2001:10000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  
  hist(with_B_mode$pval_matrix[1:2000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode$pval_matrix[1:2000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode$pval_matrix[1:2000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode$pval_matrix[2001:10000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode$pval_matrix[2001:10000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode$pval_matrix[2001:10000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))

  hist(with_B_mode_without_RIN$pval_matrix[1:2000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode_without_RIN$pval_matrix[1:2000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode_without_RIN$pval_matrix[1:2000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode_without_RIN$pval_matrix[2001:10000,1], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode_without_RIN$pval_matrix[2001:10000,2], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  hist(with_B_mode_without_RIN$pval_matrix[2001:10000,3], breaks=seq(0, 1, by = 0.05), xlim=c(0,1))
  
  dev.off()
  
}
