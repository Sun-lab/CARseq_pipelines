# source("test_CARseq.R")
# source("test_TOAST_TPM.R")
library(TOAST)
library(CARseq)

n_list = c(50, 100, 200)

for (n in n_list) {
  config = sprintf("../data/n_%d_DE_pattern_2_1_1_replicate_1", n)
  RData_file = sprintf("../results/replicability_n_%d_DE_pattern_2_1_1_replicate_1.RData", n)
  
  if (!file.exists(RData_file)) {
    load(RData_file)
  
    # partition samples into two halves and compare them:
    first_half = seq(1, n-1, by = 2)
    second_half = seq(2, n, by = 2)
    
    # TOAST_TPM
    RDatafile = sprintf(file.path(config, "simulation.RData"))
    load(RDatafile)
    
    # Will use: observed_TPM, rho_from_TPM, clinical_variables, adjusted_signature_gene_TPM (only for cell type names)
    colnames(rho_from_TPM) = colnames(adjusted_signature_gene_TPM)
    TOAST_first_half = TOAST_second_half = list()
    empty_pval_matrix = matrix(nrow=nrow(observed_TPM), ncol=ncol(rho_from_TPM))
    rownames(empty_pval_matrix) = rownames(observed_TPM)
    colnames(empty_pval_matrix) = colnames(rho_from_TPM)
    TOAST_first_half$pval_matrix = TOAST_second_half$pval_matrix = empty_pval_matrix
    TOAST_first_half$effect_size = TOAST_second_half$effect_size = empty_pval_matrix
    TOAST_first_half$lfc = TOAST_second_half$lfc = empty_pval_matrix
    
    
    # With RIN -- first half
    cell_types = colnames(rho_from_TPM[first_half, ])
    design = data.frame(RIN=clinical_variables[first_half, ], group=gl(2, round(n/2))[first_half])
    design_out = makeDesign(design, rho_from_TPM[first_half, ])
    fitted_model = fitModel(design_out, observed_TPM[, first_half])
    
    for (cell_type in cell_types) {
      res_table = csTest(fitted_model, 
                         coef = "group", 
                         cell_type = cell_type)
      # ordered_res_table = res_table[order(as.numeric(substring(rownames(res_table), 5))),]
      TOAST_first_half$pval_matrix[, cell_type] = res_table[rownames(observed_TPM), "p_value"]
      TOAST_first_half$effect_size[, cell_type] = res_table[rownames(observed_TPM), "effect_size"]
      #  define a robust effect size in TOAST as \( \beta / \max (|\mu| + |\mu + \beta|) \)
      # TOAST_first_half$robust_effect_size[, cell_type] =
      #   res_table[rownames(observed_TPM), "beta"] / max(
      #     abs(res_table[rownames(observed_TPM), "mu"]),
      #     abs(res_table[rownames(observed_TPM), "mu"] + res_table[rownames(observed_TPM), "beta"]))
      # TOAST_first_half$robust_effect_size[, cell_type] = 2 * res_table[rownames(observed_TPM), "beta"] /
      #   (abs(res_table[rownames(observed_TPM), "mu"]) +
      #    abs(res_table[rownames(observed_TPM), "mu"] + res_table[rownames(observed_TPM), "beta"]))
      
      TOAST_first_half$lfc[, cell_type] =
        log(abs(res_table[rownames(observed_TPM), "mu"] + res_table[rownames(observed_TPM), "beta"]) /
          abs(res_table[rownames(observed_TPM), "mu"]))
    }
    
    # With RIN -- second half
    cell_types = colnames(rho_from_TPM[second_half, ])
    design = data.frame(RIN=clinical_variables[second_half, ], group=gl(2, round(n/2))[second_half])
    design_out = makeDesign(design, rho_from_TPM[second_half, ])
    fitted_model = fitModel(design_out, observed_TPM[, second_half])
    
    for (cell_type in cell_types) {
      res_table = csTest(fitted_model, 
                         coef = "group", 
                         cell_type = cell_type)
      TOAST_second_half$pval_matrix[, cell_type] = res_table[rownames(observed_TPM), "p_value"]
      TOAST_second_half$effect_size[, cell_type] = res_table[rownames(observed_TPM), "effect_size"]
      #  define a robust effect size in TOAST as \( \beta / \max (|\mu| + |\mu + \beta|) \)
      # TOAST_second_half$robust_effect_size[, cell_type] =
      #   res_table[rownames(observed_TPM), "beta"] / max(
      #     abs(res_table[rownames(observed_TPM), "mu"]),
      #     abs(res_table[rownames(observed_TPM), "mu"] + res_table[rownames(observed_TPM), "beta"]))
      # TOAST_second_half$robust_effect_size[, cell_type] = 2 * res_table[rownames(observed_TPM), "beta"] /
      #   (abs(res_table[rownames(observed_TPM), "mu"]) +
      #    abs(res_table[rownames(observed_TPM), "mu"] + res_table[rownames(observed_TPM), "beta"]))
      
      TOAST_second_half$lfc[, cell_type] =
        log(abs(res_table[rownames(observed_TPM), "mu"] + res_table[rownames(observed_TPM), "beta"]) /
              abs(res_table[rownames(observed_TPM), "mu"]))
    }
    
    # CARseq
    RDatafile = sprintf(file.path(config, "simulation.RData"))
    load(RDatafile)

    # CARseq (w/ and w/o clinical variables)
    # Input: read counts; cellular frequencies.

    H = 3
    n_B = n
    M = 2
    K = 1

    x = gl(2, round(n/2))

    CARseq_first_half = run_CARseq(count_matrix = observed_read_count[, first_half],
                          cellular_proportions = rho[first_half, ],
                          groups = x[first_half],
                          formula = ~ RIN,
                          data = clinical_variables[first_half, ,drop=FALSE],
                          read_depth = d[first_half],
                          shrunken_lfc = TRUE,
                          cores = 12,
                          fix_overdispersion = FALSE
    )

    CARseq_second_half = run_CARseq(count_matrix = observed_read_count[, second_half],
                          cellular_proportions = rho[second_half, ],
                          groups = x[second_half],
                          formula = ~ RIN,
                          data = clinical_variables[second_half, ,drop=FALSE],
                          read_depth = d[second_half],
                          shrunken_lfc = TRUE,
                          cores = 12,
                          fix_overdispersion = FALSE
    )

    save(list=c("TOAST_first_half","TOAST_second_half","CARseq_first_half","CARseq_second_half"),
         file=RData_file)
  } else {
    load(RData_file)
  }
    
  # plot
  pdf(sprintf("../figures/replicability_n_%d_DE_pattern_2_1_1_replicate_1.pdf", n), width=6.5, height=2.2)  # 5, 2
  opar = par(mfrow=c(1,4), mar=c(6,2,1,1))
  # opar = par(mfrow=c(3,4), mar=c(3,2,1,1))
  
  # create a color palette to use in smoothed scatterplot
  # https://rstudio-pubs-static.s3.amazonaws.com/151690_ac65a180e03641e2adc3cb2ecf6306c3.html
  # library(RColorBrewer)
  
  smoothScatter(CARseq_first_half$lfc[1:2000,1], CARseq_second_half$lfc[1:2000,1],
                xlab = "CARseq LFC", ylab = "", xlim=c(-5,5), ylim=c(-5,5), pch = 15, cex = .1)
  abline(h=0, v=0, col=rgb(0.7,0.7,0.7,0.5), lty="dotted", lwd = 0.75)
  abline(h=c(-log(2), log(2)), v=c(-log(2), log(2)), col=rgb(0.7,0.7,0.7,0.5), lty="dashed", lwd = 0.75)
  legend("topright", bty="n", legend=sprintf("cor=%.3f", cor(CARseq_first_half$lfc[1:2000,1], CARseq_second_half$lfc[1:2000,1], method="spearman", use = "complete.obs")))
  
  smoothScatter(CARseq_first_half$shrunken_lfc[1:2000,1], CARseq_second_half$shrunken_lfc[1:2000,1],
                xlab = "CARseq shrunken LFC", ylab = "", xlim=c(-5,5), ylim=c(-5,5), pch = 15, cex = .1)
  abline(h=0, v=0, col=rgb(0.7,0.7,0.7,0.5), lty="dotted", lwd = 0.75)
  abline(h=c(-log(2), log(2)), v=c(-log(2), log(2)), col=rgb(0.7,0.7,0.7,0.5), lty="dashed", lwd = 0.75)
  legend("topright", bty="n", legend=sprintf("cor=%.3f", cor(CARseq_first_half$shrunken_lfc[1:2000,1], CARseq_second_half$shrunken_lfc[1:2000,1], method="spearman", use = "complete.obs")))
  
  smoothScatter(TOAST_first_half$lfc[1:2000,1], TOAST_second_half$lfc[1:2000,1],
                xlab = "TOAST LFC", ylab = "", xlim=c(-5,5), ylim=c(-5,5), pch = 15, cex = .1)
  abline(h=0, v=0, col=rgb(0.7,0.7,0.7,0.5), lty="dotted", lwd = 0.75)
  abline(h=c(-log(2), log(2)), v=c(-log(2), log(2)), col=rgb(0.7,0.7,0.7,0.5), lty="dashed", lwd = 0.75)
  legend("topright", bty="n",  legend=sprintf("cor=%.3f", cor(TOAST_first_half$lfc[1:2000,1], TOAST_second_half$lfc[1:2000,1], method="spearman", use = "complete.obs")))
  
  smoothScatter(TOAST_first_half$effect_size[1:2000,1], TOAST_second_half$effect_size[1:2000,1],
                xlab = "TOAST effect size", ylab = "", xlim=c(-5,5), ylim=c(-5,5), pch = 15, cex = .1)
  abline(h=0, v=0, col=rgb(0.7,0.7,0.7,0.5), lty="dotted", lwd = 0.75)
  legend("topright", bty="n", legend=sprintf("cor=%.3f", cor(TOAST_first_half$effect_size[1:2000,1], TOAST_second_half$effect_size[1:2000,1], method="spearman", use = "complete.obs")))
  
  par(opar)
  dev.off()
}