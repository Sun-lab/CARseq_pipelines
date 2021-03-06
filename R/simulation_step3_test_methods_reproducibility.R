# source("test_CARseq.R")
# source("test_TOAST_TPM.R")
library(TOAST)
library(CARseq)
library(ggplot2)
library(ggpointdensity)
library(magrittr)

n_list = c(50, 100, 200)

label_data_frame = effect_sizes = NULL

for (n in n_list) {
  config = sprintf("../simulation/n_%d_DE_pattern_2_1_1_replicate_1", n)
  RData_file = sprintf("../results/reproducibility_n_%d_DE_pattern_2_1_1_replicate_1.RData", n)
  
  if (!file.exists(RData_file)) {
  
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

    H = 6
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
    
  # library(patchwork)
  # plot
  # opar = par(mfrow=c(1,4), mar=c(6,2,1,1))
  # opar = par(mfrow=c(3,4), mar=c(3,2,1,1))
  
  # create a color palette to use in smoothed scatterplot
  # https://rstudio-pubs-static.s3.amazonaws.com/151690_ac65a180e03641e2adc3cb2ecf6306c3.html
  # library(RColorBrewer)
  
  h = 2  # the cell type having DE
  
  frac_of_correct_direction_in_both_reps = rep(NA, 4)   # fraction of genes that are in the same direction of DE in both half of samples
  CARseq_lfc = data.frame(first_half = CARseq_first_half$lfc[1:2000,h],
                          second_half = CARseq_second_half$lfc[1:2000,h],
                          method = "CARseq_LFC",
                          n = sprintf("Total sample size in a half = %d", n * 0.5))
  frac_of_correct_direction_in_both_reps[1] = (sum(CARseq_first_half$lfc[1:1000,h] < 0 & CARseq_second_half$lfc[1:1000,h] < 0, na.rm = TRUE) +
                                               sum(CARseq_first_half$lfc[1001:2000,h] > 0 & CARseq_second_half$lfc[1001:2000,h] > 0, na.rm = TRUE)) / 2000
  
  CARseq_shrunken_lfc = data.frame(first_half = CARseq_first_half$shrunken_lfc[1:2000,h],
                                   second_half = CARseq_second_half$shrunken_lfc[1:2000,h],
                                   method = "CARseq_shrunken_LFC",
                                   n = sprintf("Total sample size in a half = %d", n * 0.5))
  frac_of_correct_direction_in_both_reps[2] = (sum(CARseq_first_half$shrunken_lfc[1:1000,h] < 0 & CARseq_second_half$shrunken_lfc[1:1000,h] < 0, na.rm = TRUE) +
                                               sum(CARseq_first_half$shrunken_lfc[1001:2000,h] > 0 & CARseq_second_half$shrunken_lfc[1001:2000,h] > 0, na.rm = TRUE)) / 2000
  
  TOAST_lfc = data.frame(first_half = TOAST_first_half$lfc[1:2000,h],
                         second_half = TOAST_second_half$lfc[1:2000,h],
                         method = "TOAST_LFC",
                         n = sprintf("Total sample size in a half = %d", n * 0.5))
  frac_of_correct_direction_in_both_reps[3] = (sum(TOAST_first_half$lfc[1:1000,h] < 0 & TOAST_second_half$lfc[1:1000,h] < 0, na.rm = TRUE) +
                                               sum(TOAST_first_half$lfc[1001:2000,h] > 0 & TOAST_second_half$lfc[1001:2000,h] > 0, na.rm = TRUE)) / 2000
  
  TOAST_effect_size = data.frame(first_half = TOAST_first_half$effect_size[1:2000,h],
                                 second_half = TOAST_second_half$effect_size[1:2000,h],
                                 method = "TOAST_effect_size",
                                 n = sprintf("Total sample size in a half = %d", n * 0.5))
  frac_of_correct_direction_in_both_reps[4] = (sum(TOAST_first_half$effect_size[1:1000,h] < 0 & TOAST_second_half$effect_size[1:1000,h] < 0, na.rm = TRUE) +
                                               sum(TOAST_first_half$effect_size[1001:2000,h] > 0 & TOAST_second_half$effect_size[1001:2000,h] > 0, na.rm = TRUE)) / 2000
  
  effect_sizes = effect_sizes %>%
    rbind(CARseq_lfc) %>%
    rbind(CARseq_shrunken_lfc) %>%
    rbind(TOAST_lfc) %>%
    rbind(TOAST_effect_size)
  label_data_frame = label_data_frame %>% rbind(data.frame(
      label = sprintf("%.1f%% correct", 100 * frac_of_correct_direction_in_both_reps),
      method = c("CARseq_LFC", "CARseq_shrunken_LFC", "TOAST_LFC", "TOAST_effect_size"),
      n = sprintf("Total sample size in a half = %d", n * 0.5)
    ))
  
  # h = 2
  # 
  # cor_between_estimates = rep(NA, 4)
  # CARseq_lfc = data.frame(first_half = CARseq_first_half$lfc[1:2000,h],
  #                         second_half = CARseq_second_half$lfc[1:2000,h],
  #                         method = "CARseq_LFC",
  #                         n = sprintf("Total sample size in a half = %d", n * 0.5))
  # cor_between_estimates[1] = cor(CARseq_first_half$lfc[1:2000,h], CARseq_second_half$lfc[1:2000,h], method="spearman", use = "complete.obs")
  # 
  # CARseq_shrunken_lfc = data.frame(first_half = CARseq_first_half$shrunken_lfc[1:2000,h],
  #                                  second_half = CARseq_second_half$shrunken_lfc[1:2000,h],
  #                                  method = "CARseq_shrunken_LFC",
  #                                  n = sprintf("Total sample size in a half = %d", n * 0.5))
  # cor_between_estimates[2] = cor(CARseq_first_half$shrunken_lfc[1:2000,h], CARseq_second_half$shrunken_lfc[1:2000,h], method="spearman", use = "complete.obs")
  # 
  # TOAST_lfc = data.frame(first_half = TOAST_first_half$lfc[1:2000,h],
  #                        second_half = TOAST_second_half$lfc[1:2000,h],
  #                        method = "TOAST_LFC",
  #                        n = sprintf("Total sample size in a half = %d", n * 0.5))
  # cor_between_estimates[3] = cor(TOAST_first_half$lfc[1:2000,h], TOAST_second_half$lfc[1:2000,h], method="spearman", use = "complete.obs")
  # 
  # TOAST_effect_size = data.frame(first_half = TOAST_first_half$effect_size[1:2000,h],
  #                                second_half = TOAST_second_half$effect_size[1:2000,h],
  #                                method = "TOAST_effect_size",
  #                                n = sprintf("Total sample size in a half = %d", n * 0.5))
  # cor_between_estimates[4] = cor(TOAST_first_half$effect_size[1:2000,h], TOAST_second_half$effect_size[1:2000,h], method="spearman", use = "complete.obs")
  # 
  # effect_sizes = effect_sizes %>%
  #   rbind(CARseq_lfc) %>%
  #   rbind(CARseq_shrunken_lfc) %>%
  #   rbind(TOAST_lfc) %>%
  #   rbind(TOAST_effect_size)
  # label_data_frame = label_data_frame %>% rbind(data.frame(
  #   label = sprintf("cor=%.3f", cor_between_estimates),
  #   method = c("CARseq_LFC", "CARseq_shrunken_LFC", "TOAST_LFC", "TOAST_effect_size"),
  #   n = sprintf("Total sample size in a half = %d", n * 0.5)
  # ))
  
  
  # # obsolete code of smoothScatter
  # smoothScatter(CARseq_first_half$lfc[1:2000,h], CARseq_second_half$lfc[1:2000,h],
  #               xlab = "CARseq LFC", ylab = "", xlim=c(-5,5), ylim=c(-5,5), pch = 15, cex = .1)
  # abline(h=0, v=0, col=rgb(0.7,0.7,0.7,0.5), lty="dotted", lwd = 0.75)
  # abline(h=c(-log(2), log(2)), v=c(-log(2), log(2)), col=rgb(0.7,0.7,0.7,0.5), lty="dashed", lwd = 0.75)
  # legend("topright", bty="n", legend=sprintf("cor=%.3f", cor(CARseq_first_half$lfc[1:2000,h], CARseq_second_half$lfc[1:2000,h], method="spearman", use = "complete.obs")))
  # 
  # smoothScatter(CARseq_first_half$shrunken_lfc[1:2000,h], CARseq_second_half$shrunken_lfc[1:2000,h],
  #               xlab = "CARseq shrunken LFC", ylab = "", xlim=c(-5,5), ylim=c(-5,5), pch = 15, cex = .1)
  # abline(h=0, v=0, col=rgb(0.7,0.7,0.7,0.5), lty="dotted", lwd = 0.75)
  # abline(h=c(-log(2), log(2)), v=c(-log(2), log(2)), col=rgb(0.7,0.7,0.7,0.5), lty="dashed", lwd = 0.75)
  # legend("topright", bty="n", legend=sprintf("cor=%.3f", cor(CARseq_first_half$shrunken_lfc[1:2000,h], CARseq_second_half$shrunken_lfc[1:2000,h], method="spearman", use = "complete.obs")))
  # 
  # smoothScatter(TOAST_first_half$lfc[1:2000,h], TOAST_second_half$lfc[1:2000,h],
  #               xlab = "TOAST LFC", ylab = "", xlim=c(-5,5), ylim=c(-5,5), pch = 15, cex = .1)
  # abline(h=0, v=0, col=rgb(0.7,0.7,0.7,0.5), lty="dotted", lwd = 0.75)
  # abline(h=c(-log(2), log(2)), v=c(-log(2), log(2)), col=rgb(0.7,0.7,0.7,0.5), lty="dashed", lwd = 0.75)
  # legend("topright", bty="n",  legend=sprintf("cor=%.3f", cor(TOAST_first_half$lfc[1:2000,h], TOAST_second_half$lfc[1:2000,h], method="spearman", use = "complete.obs")))
  # 
  # smoothScatter(TOAST_first_half$effect_size[1:2000,h], TOAST_second_half$effect_size[1:2000,h],
  #               xlab = "TOAST effect size", ylab = "", xlim=c(-5,5), ylim=c(-5,5), pch = 15, cex = .1)
  # abline(h=0, v=0, col=rgb(0.7,0.7,0.7,0.5), lty="dotted", lwd = 0.75)
  # legend("topright", bty="n", legend=sprintf("cor=%.3f", cor(TOAST_first_half$effect_size[1:2000,h], TOAST_second_half$effect_size[1:2000,h], method="spearman", use = "complete.obs")))

}

effect_sizes$method = gsub("_", " ", effect_sizes$method)
label_data_frame$method = gsub("_", " ", label_data_frame$method)

# store the data illustrated in Fig 3 of the main text
write.csv(effect_sizes, file = "../results/_figures_data/Fig3.csv")

pdf(sprintf("../figures/reproducibility_DE_pattern_2_1_1_replicate_1.pdf", n), width=8.5, height=7.2)  # 6.5, 2.2
ggplot(effect_sizes, aes(x=first_half, y=second_half)) +
  geom_pointdensity(size = .25) + scale_color_viridis_c() +
  facet_grid(n ~ method) +
  xlab("First half of samples") + ylab("Second half of samples") +
  xlim(c(-5,5)) + ylim(c(-5,5)) +
  geom_hline(yintercept = c(-log(2), log(2)), lty="dashed", lwd = 0.5, col=rgb(0.7,0.7,0.7,0.5)) +
  geom_hline(yintercept = 0, lty="dotted", lwd = 0.5, col=rgb(0.7,0.7,0.7,0.5)) +
  geom_vline(xintercept = c(-log(2), log(2)), lty="dashed", lwd = 0.5, col=rgb(0.7,0.7,0.7,0.5)) +
  geom_vline(xintercept = 0, lty="dotted", lwd = 0.5, col=rgb(0.7,0.7,0.7,0.5)) +
  # theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom") +
  geom_text(data = label_data_frame, mapping = aes(x = -2.6, y = 4.9, label = label))
dev.off()
