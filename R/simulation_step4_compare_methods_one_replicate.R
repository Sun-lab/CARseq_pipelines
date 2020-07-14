# This script compares csSAM, TOAST, CIBERSORTx, CARseq in one replicate.
# FDR:
#   sensitivity
#   1-precision (false discovery rate)
# p-values (using 0.05 as cutoff):
#   power (gene/cell type pairs being differentially expressed)
#   type I error (genes being differentially expressed in other cell types; genes not being differentially expressed)

# config = "n_100_DE_pattern_2_1_1_replicate_1"

# Calculate relevant metrics for a given config and method.
# The res need to have the $pval_matrix member.
compute_metrics_from_pval_matrix = function(res, method, config) {
  if (is.null(res)) return(NA)
  if (is.null(res$pval_matrix)) {
    if ("pval_mat" %in% names(res)) {
      res$pval_matrix = res$pval_mat
    } else if ("p" %in% names(res)) {
      res$pval_matrix = res$p
    }
  }
  pvalue_cutoff = 0.05
  fdr_cutoff = 0.1
  H = 3  # number of cell types
  # Genes that are differentially expressed in cell type 1
  gene_indices_DE = 1:2000
  # Genes that are not differentially expressed in any cell type
  gene_indices_non_DE = 2001:10000
  
  number_of_discoveries_list = list()
  number_of_true_discoveries_list = list()
  for (h in seq_len(H)) {
    padj = p.adjust(res$pval_matrix[, h], method="BH")
    number_of_true_discoveries_list[[h]] = sum(padj[gene_indices_DE] < fdr_cutoff, na.rm=TRUE)
    number_of_discoveries_list[[h]] = sum(padj < fdr_cutoff, na.rm=TRUE)
  }
  # Calculate total FDR and sensitivity across cell types
  fold_changes = as.numeric(unlist(strsplit(gsub(".*_pattern_|_replicate_.*", "", config), "_")))
  has_fold_changes = which(fold_changes != 1)
  fdr = 1 - sum(unlist(number_of_true_discoveries_list[has_fold_changes])) /
      sum(unlist(number_of_discoveries_list))
  sensitivity = sum(unlist(number_of_true_discoveries_list[has_fold_changes])) /
      (length(has_fold_changes) * length(gene_indices_DE))
  names(fdr) = "fdr"
  names(sensitivity) = "sensitivity"
  names(number_of_discoveries_list) = paste0("number_of_discoveries_ct", seq_len(H))
  names(number_of_true_discoveries_list) = paste0("number_of_true_discoveries_ct", seq_len(H))

  # power (DE genes)
  power_list = list()
  for (h in seq_len(H)) {
    pval = na.omit(res$pval_matrix[gene_indices_DE, h])
    power_list[[h]] = sum(pval < pvalue_cutoff) / length(gene_indices_DE)
  }
  names(power_list) = paste0("power_ct", seq_len(H))
  # type I error (non DE genes)
  type_I_error_list = list()
  for (h in seq_len(H)) {
    pval = na.omit(res$pval_matrix[gene_indices_non_DE, h])
    type_I_error_list[[h]] = sum(pval < pvalue_cutoff) / length(pval)
  }
  names(type_I_error_list) = paste0("type_I_error_ct", seq_len(H))

  data.frame(method = method,
             sensitivity,
             fdr,
             number_of_true_discoveries_list,
             number_of_discoveries_list,
             power_list,
             type_I_error_list)
}

# Calculate relevant metrics for a given config and method.
# The res need to have the $fdr_matrix member.
compute_metrics_from_fdr_matrix = function(res, method, config) {
  if (is.null(res)) return(NA)
  pvalue_cutoff = 0.05
  fdr_cutoff = 0.1
  H = 3  # number of cell types
  # Genes that are differentially expressed in cell type 1
  gene_indices_DE = 1:2000
  # Genes that are not differentially expressed in any cell type
  gene_indices_non_DE = 2001:10000
  
  number_of_discoveries_list = list()
  number_of_true_discoveries_list = list()
  for (h in seq_len(H)) {
    padj = res$fdr_matrix[, h]
    number_of_true_discoveries_list[[h]] = sum(padj[gene_indices_DE] < fdr_cutoff, na.rm=TRUE)
    number_of_discoveries_list[[h]] = sum(padj < fdr_cutoff, na.rm=TRUE)
  }
  # Calculate total FDR and sensitivity across cell types
  fold_changes = as.numeric(unlist(strsplit(gsub(".*_pattern_|_replicate_.*", "", config), "_")))
  has_fold_changes = which(fold_changes != 1)
  fdr = 1 - sum(unlist(number_of_true_discoveries_list[has_fold_changes])) /
    sum(unlist(number_of_discoveries_list))
  sensitivity = sum(unlist(number_of_true_discoveries_list[has_fold_changes])) /
    (length(has_fold_changes) * length(gene_indices_DE))
  names(fdr) = "fdr"
  names(sensitivity) = "sensitivity"
  names(number_of_discoveries_list) = paste0("number_of_discoveries_ct", seq_len(H))
  names(number_of_true_discoveries_list) = paste0("number_of_true_discoveries_ct", seq_len(H))
  
  # power (DE genes)
  power_list = as.list(rep(NA, H))
  names(power_list) = paste0("power_ct", seq_len(H))
  # type I error (non DE genes)
  type_I_error_list = as.list(rep(NA, H))
  names(type_I_error_list) = paste0("type_I_error_ct", seq_len(H))

  data.frame(method = method,
             sensitivity,
             fdr,
             number_of_true_discoveries_list,
             number_of_discoveries_list,
             power_list,
             type_I_error_list)
}

collect_evaluation_metrics = function(config) {
  csSAM_res = new.env()
  TOAST_res = new.env()
  TOAST_TPM_res = new.env()
  CIBERSORTx_res = new.env()
  CARseq_res = new.env()
  tryCatch(load(file=file.path(config, "csSAM_res.RData"), envir = csSAM_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path(config, "TOAST_res.RData"), envir = TOAST_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path(config, "TOAST_TPM_res.RData"), envir = TOAST_TPM_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path(config, "CIBERSORTx_res.RData"), envir = CIBERSORTx_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path(config, "CARseq_res.RData"), envir = CARseq_res), error=function(e) {print(e)})

  metrics_list = list()
  ############################################################
  #          CARseq with RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(CARseq_res$with_RIN, "CARseq_with_RIN", config)

  ############################################################
  #          CARseq without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(CARseq_res$without_RIN, "CARseq_without_RIN", config)
  
  ############################################################
  #          CIBERSORTx with B mode
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(CIBERSORTx_res$with_B_mode, "CIBERSORTx_with_B_mode", config)
  
  ############################################################
  #          CIBERSORTx with B mode without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(CIBERSORTx_res$with_B_mode_without_RIN, "CIBERSORTx_with_B_mode_without_RIN", config)
  
  ############################################################
  #          CIBERSORTx without B mode
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(CIBERSORTx_res$without_B_mode, "CIBERSORTx_without_B_mode", config)
  
  ############################################################
  #          CIBERSORTx without B mode without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(CIBERSORTx_res$without_B_mode_without_RIN, "CIBERSORTx_without_B_mode_without_RIN", config)
  
  ############################################################
  #          TOAST using counts with RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_res$with_RIN, "TOAST_with_RIN", config)
  
  ############################################################
  #          TOAST using counts without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_res$without_RIN, "TOAST_without_RIN", config)
  
  ############################################################
  #          TOAST using TPM with RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_TPM_res$with_RIN, "TOAST_TPM_with_RIN", config)

  ############################################################
  #          TOAST using TPM without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_TPM_res$without_RIN, "TOAST_TPM_without_RIN", config)

  ############################################################
  #          csSAM (fdr only, config)
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_fdr_matrix(csSAM_res, "csSAM", config)

  # return a data frame
  do.call(rbind, metrics_list)
}


collect_evaluation_metrics_with_params = function(n, DE_pattern, DE_pattern_index, replicate) {
  # Generate the RData file name
  RDatafolder = sprintf("n_%s_DE_pattern_%s_replicate_%s",
                        n,
                        paste(DE_pattern, collapse="_"),
                        replicate)
  metrics = collect_evaluation_metrics(RDatafolder)
  m = nrow(metrics)
  data.frame(n = rep(n, m),
             DE_pattern = rep(paste(DE_pattern, collapse="_"), m),
             DE_pattern_index = rep(DE_pattern_index, m),
             replicate = rep(replicate, m),
             metrics)
}

# setup
fdr_cutoff = 0.1
pvalue_cutoff = 0.05
n_list = c(50, 100, 200)  # total number of samples, 50, 100, 200
# add cell type-specific DE and compute the expected mean
fold_changes = c(2, 3, 4)  # 2,3,4
replicates = 1
metrics_list = list()
for (n in n_list) {
  for (fold_change in fold_changes) {
    DE_pattern_list = list(c(fold_change, 1, 1),
                           c(fold_change, fold_change, 1),
                           c(fold_change, round(1000/fold_change)/1000, 1),
                           c(1, fold_change, 1))
    for (DE_pattern_index in seq_along(DE_pattern_list)) {
      DE_pattern = DE_pattern_list[[DE_pattern_index]]
      for (replicate in replicates) {
        metrics_list[[1+length(metrics_list)]] = collect_evaluation_metrics_with_params(n = n,
                                                                                        DE_pattern = DE_pattern,
                                                                                        DE_pattern_index = DE_pattern_index,
                                                                                        replicate = replicate)
      }
    }
  }
}
# create a data frame
metrics = do.call(rbind, metrics_list)
write.csv(metrics, file = "../results/compare_methods_metrics_one_replicate.csv", quote = FALSE, row.names = FALSE)
# If there are no discoveries at all, then sensitivity is 0, but fdr is not defined.
# Since we are not actually making false discoveries,
# to make the result show on the plot, when fdr = 0/0 = NaN, we set it to 0.
# Likewise, if the # of total discoveries are no more than 20 for a cell type,
# (there are 2000 genes with DE and 8000 genes without DE), then they are all discarded,
# and we set fdr = 0 & sensitivity = 0.

library(tidyr)
library(ggplot2)
metrics = read.csv(file = "../results/compare_methods_metrics_one_replicate.csv", as.is = TRUE)
metrics$fdr[(metrics$number_of_discoveries_ct1 +
               metrics$number_of_discoveries_ct2 +
               metrics$number_of_discoveries_ct3) <= 10] = 0

names(metrics) = names(metrics) %>%
  gsub("power", "DEGenes__PositiveRate", x=.) %>%
  gsub("type_I_error", "NonDEGenes_PositiveRate", x=.)

metrics$method = metrics$method %>%
  gsub("TOAST_w", "TOAST_count_w", x=.) %>%
  gsub("_with_RIN", "", x=.) %>%
  gsub("without_RIN", "w/o_covariates", x=.) %>%
  gsub("_without_", "_w/o_", x=.) %>%
  gsub("_with_", "_w/_", x=.)

pvalue_cutoff = 0.05
fdr_cutoff = 0.1

# Figure 1:
# The cutoff being used is FDR < 0.1.

# only show methods with covariates
metrics_with_covariates = metrics[!grepl("csSAM|w/o_covariates", metrics$method), ]
metrics_with_covariates$method = metrics_with_covariates$method %>% gsub("_w/?o_covariates", "", x=.)
metrics_with_covariates$method = as.factor(metrics_with_covariates$method)

# only show methods without covariates
metrics_without_covariates = metrics[grepl("csSAM|w/o_covariates", metrics$method), ]
metrics_without_covariates$method = metrics_without_covariates$method %>% gsub("_w/?o_covariates", "", x=.)
metrics_without_covariates$method = factor(metrics_without_covariates$method,
                                           c(levels=levels(metrics_with_covariates$method), "csSAM"))

############################################################
# Figure 1 (combined draft version)
############################################################
library(ggpubr)
library(wesanderson)
g1_with_covariates = ggplot(metrics_with_covariates, aes(x = fdr, y = sensitivity)) +
  geom_point(aes(col = method, shape = method), stroke = 0.2) +
  facet_grid(DE_pattern ~ n) +
  geom_vline(xintercept=fdr_cutoff, linetype="solid", color = "black", size = 0.2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_shape_manual(values=1:nlevels(metrics_with_covariates$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1", n = 6, type = "continuous")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        # legend.position = "bottom",
        legend.text=element_text(size=8),
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

g1_without_covariates = ggplot(metrics_without_covariates, aes(x = fdr, y = sensitivity)) +
  geom_point(aes(col = method, shape = method), stroke = 0.2) +
  facet_grid(DE_pattern ~ n) +
  geom_vline(xintercept=fdr_cutoff, linetype="solid", color = "black", size = 0.2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_shape_manual(values=1:nlevels(metrics_without_covariates$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1", n = 6, type = "continuous")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        # legend.position = "bottom",
        legend.text=element_text(size=8),
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

pdf("../figures/compare_methods_fdr_one_replicate_draft_version_p1.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 1),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 1),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_one_replicate_draft_version_p2.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 2),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 2),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_one_replicate_draft_version_p3.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 3),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 3),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_one_replicate_draft_version_p4.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 4),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 4),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()
