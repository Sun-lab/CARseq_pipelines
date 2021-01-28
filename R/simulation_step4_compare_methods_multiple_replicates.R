# This script compares csSAM, TOAST, CARseq in 10 replicates.
# FDR:
#   sensitivity
#   1-precision (false discovery rate)
# p-values (using 0.05 as cutoff):
#   power (gene/cell type pairs being differentially expressed)
#   type I error (genes being differentially expressed in other cell types; genes not being differentially expressed)

# config = "n_100_DE_pattern_2_1_1_replicate_1"

# Calculate relevant metrics for a given config and method.
# The res need to have the $pval_matrix member.
compute_metrics_from_pval_matrix = function(res, method, config, gene_expression_quantile_range = NULL) {
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
  H = 6  # number of cell types
  # Genes that are differentially expressed in cell type 1
  gene_indices_DE = 1:2000
  # Genes that are not differentially expressed in any cell type
  gene_indices_non_DE = 2001:10000
  # If we are required to stratify by gene expression, we need to load the counts matrix and rank by geometric mean:
  if (!is.null(gene_expression_quantile_range)) {
    simulation_RDatafile = sprintf(file.path("..", "simulation", config, "simulation.RData"))
    simulation_env = new.env()
    load(simulation_RDatafile, envir = simulation_env)
    gene_expression_medians = rowMedians(simulation_env$observed_read_count)
    gene_expression_medians_quantile = quantile(gene_expression_medians, probs = gene_expression_quantile_range)
    gene_indices_within_quantile_range = which(gene_expression_medians >= gene_expression_medians_quantile[1] &
                                               gene_expression_medians <= gene_expression_medians_quantile[2])
    gene_indices_DE = intersect(gene_indices_within_quantile_range, gene_indices_DE)
    gene_indices_non_DE = intersect(gene_indices_within_quantile_range, gene_indices_non_DE)
  }
  
  number_of_discoveries_list = list()
  number_of_true_discoveries_list = list()
  for (h in seq_len(H)) {
    padj = p.adjust(res$pval_matrix[, h], method="BH")
    number_of_true_discoveries_list[[h]] = sum(padj[gene_indices_DE] < fdr_cutoff, na.rm=TRUE)
    number_of_discoveries_list[[h]] = sum(padj < fdr_cutoff, na.rm=TRUE)
  }
  # Calculate total FDR and sensitivity across cell types
  fold_changes = as.numeric(unlist(strsplit(gsub(".*_pattern_|_replicate_.*", "", config), "_")))[c(3, 1, 2, 3, 3, 3)]
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

  Elapsed_Time_sec = Peak_RAM_Used_MiB = NA
  # replicates based on "ICeDT" have the time and memory set to be NA:
  if (!is.null(res$Elapsed_Time_sec)) {
    if (!grepl("ICeDT", method)) {
      Elapsed_Time_sec = res$Elapsed_Time_sec
      if (grepl("CARseq", method)) {
        Peak_RAM_Used_MiB = res$Peak_RAM_Used_MiB_no_parallel * 12 + res$Peak_RAM_Used_MiB_parallel
      } else {
        Peak_RAM_Used_MiB = res$Peak_RAM_Used_MiB
      }
    }
  }
  
  data.frame(method = method,
             sensitivity,
             fdr,
             number_of_true_discoveries_list,
             number_of_discoveries_list,
             power_list,
             type_I_error_list,
             Elapsed_Time_sec = Elapsed_Time_sec,
             Peak_RAM_Used_MiB = Peak_RAM_Used_MiB)
}

# Calculate relevant metrics for a given config and method.
# The res need to have the $fdr_matrix member.
compute_metrics_from_fdr_matrix = function(res, method, config) {
  if (is.null(res)) return(NA)
  pvalue_cutoff = 0.05
  fdr_cutoff = 0.1
  H = 6  # number of cell types
  # Genes that are differentially expressed in cell type 1
  gene_indices_DE = 1:2000
  # Genes that are not differentially expressed in any cell type
  gene_indices_non_DE = 2001:10000
  
  number_of_discoveries_list = list()
  number_of_true_discoveries_list = list()
  if (is.null(res$fdr_matrix)) {
    res$fdr_matrix = t(res$deconvResults$sigGene.csSAM)
  }
  for (h in seq_len(H)) {
    padj = res$fdr_matrix[, h]
    number_of_true_discoveries_list[[h]] = sum(padj[gene_indices_DE] < fdr_cutoff, na.rm=TRUE)
    number_of_discoveries_list[[h]] = sum(padj < fdr_cutoff, na.rm=TRUE)
  }
  # Calculate total FDR and sensitivity across cell types
  fold_changes = as.numeric(unlist(strsplit(gsub(".*_pattern_|_replicate_.*", "", config), "_")))[c(3, 1, 2, 3, 3, 3)]
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

  Elapsed_Time_sec = Peak_RAM_Used_MiB = NA
  # replicates based on "ICeDT" have the time and memory set to be NA:
  if (!grepl("ICeDT", method)) {
    Elapsed_Time_sec = res$Elapsed_Time_sec
    if (grepl("CARseq", method)) {
      Peak_RAM_Used_MiB = res$Peak_RAM_Used_MiB_no_parallel * 12 + res$Peak_RAM_Used_MiB_parallel
    } else {
      Peak_RAM_Used_MiB = res$Peak_RAM_Used_MiB
    }
  }
  
  data.frame(method = method,
             sensitivity,
             fdr,
             number_of_true_discoveries_list,
             number_of_discoveries_list,
             power_list,
             type_I_error_list,
             Elapsed_Time_sec = Elapsed_Time_sec,
             Peak_RAM_Used_MiB = Peak_RAM_Used_MiB)
}

collect_evaluation_metrics = function(config) {
  csSAM_res = new.env()
  csSAM_ICeDT_res = new.env()
  TOAST_TPM_res = new.env()
  TOAST_TPM_ICeDT_res = new.env()
  TOAST_counts_res = new.env()
  TOAST_counts_ICeDT_res = new.env()
  CARseq_res = new.env()
  CARseq_ICeDT_res = new.env()
  tryCatch(load(file=file.path("..", "simulation", config, "csSAM_res.RData"), envir = csSAM_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path("..", "simulation", config, "csSAM_ICeDT_res.RData"), envir = csSAM_ICeDT_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path("..", "simulation", config, "TOAST_TPM_res.RData"), envir = TOAST_TPM_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path("..", "simulation", config, "TOAST_TPM_ICeDT_res.RData"), envir = TOAST_TPM_ICeDT_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path("..", "simulation", config, "TOAST_res.RData"), envir = TOAST_counts_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path("..", "simulation", config, "TOAST_ICeDT_res.RData"), envir = TOAST_counts_ICeDT_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path("..", "simulation", config, "CARseq_res.RData"), envir = CARseq_res), error=function(e) {print(e)})
  tryCatch(load(file=file.path("..", "simulation", config, "CARseq_ICeDT_res.RData"), envir = CARseq_ICeDT_res), error=function(e) {print(e)})

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
  #          CARseq (ICeDT) with RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(CARseq_ICeDT_res$with_RIN, "CARseq_ICeDT_with_RIN", config)
  
  ############################################################
  #          CARseq (ICeDT) without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(CARseq_ICeDT_res$without_RIN, "CARseq_ICeDT_without_RIN", config)
  
  ############################################################
  #          TOAST using TPM with RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_TPM_res$with_RIN, "TOAST_TPM_with_RIN", config)

  ############################################################
  #          TOAST using TPM without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_TPM_res$without_RIN, "TOAST_TPM_without_RIN", config)

  ############################################################
  #          TOAST using TPM (ICeDT) with RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_TPM_ICeDT_res$with_RIN, "TOAST_TPM_ICeDT_with_RIN", config)
  
  ############################################################
  #          TOAST using TPM (ICeDT) without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_TPM_ICeDT_res$without_RIN, "TOAST_TPM_ICeDT_without_RIN", config)
  
  ############################################################
  #          TOAST using counts with RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_counts_res$with_RIN, "TOAST_counts_with_RIN", config)
  
  ############################################################
  #          TOAST using counts without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_counts_res$without_RIN, "TOAST_counts_without_RIN", config)
  
  ############################################################
  #          TOAST using counts (ICeDT) with RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_counts_ICeDT_res$with_RIN, "TOAST_counts_ICeDT_with_RIN", config)
  
  ############################################################
  #          TOAST using counts (ICeDT) without RIN
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_pval_matrix(TOAST_counts_ICeDT_res$without_RIN, "TOAST_counts_ICeDT_without_RIN", config)
  
  ############################################################
  #          csSAM (fdr only, config)
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_fdr_matrix(csSAM_res, "csSAM", config)
  
  ############################################################
  #          csSAM (ICeDT) (fdr only, config)
  ############################################################
  metrics_list[[1+length(metrics_list)]] = compute_metrics_from_fdr_matrix(csSAM_ICeDT_res, "csSAM_ICeDT", config)
  
  ############################################################
  #          ICeDT cell fraction estimate RMSE
  ############################################################
  # # RMSE_rho
  # names(csSAM_ICeDT_res$rho_summary[[3]]) =
  #     paste0("RMSE_rho_", seq_along(csSAM_ICeDT_res$rho_summary[[3]]))
  # # RMSE_rho_from_TPM
  # names(csSAM_ICeDT_res$rho_summary[[6]]) =
  #     paste0("RMSE_rho_from_TPM_", seq_along(csSAM_ICeDT_res$rho_summary[[6]]))

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
replicates = 1:10
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
write.csv(metrics, file = "../results/compare_methods_metrics_multiple_replicates.csv", quote = FALSE, row.names = FALSE)
# If there are no discoveries at all, then sensitivity is 0, but fdr is not defined.
# Since we are not actually making false discoveries,
# to make the result show on the plot, when fdr = 0/0 = NaN, we set it to 0.
# Likewise, if the # of total discoveries are no more than 20 for a cell type,
# (there are 2000 genes with DE and 8000 genes without DE), then they are all discarded,
# and we set fdr = 0 & sensitivity = 0.


############################################################
# Read parsed summary
############################################################

library(tidyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)
library(matrixStats)
metrics = read.csv(file = "../results/compare_methods_metrics_multiple_replicates.csv", as.is = TRUE)
metrics$fdr[(metrics$number_of_discoveries_ct1 +
             metrics$number_of_discoveries_ct2 +
             metrics$number_of_discoveries_ct3 +
             metrics$number_of_discoveries_ct4 +
             metrics$number_of_discoveries_ct5 +
             metrics$number_of_discoveries_ct6) <= 10] = 0
metrics_all = metrics  # a copy of metrics that we will keep intact


############################################################
# plot time & memory usage
############################################################
metrics_without_covariates_n_200 = metrics[metrics$n == 200 & !grepl("with_RIN|ICeDT|TOAST_counts", metrics$method), ]
metrics_without_covariates_n_200$method = as.factor(as.character(metrics_without_covariates_n_200$method))
levels(metrics_without_covariates_n_200$method) = c("CARseq (12 threads)", "csSAM", "TOAST")
# metrics_without_covariates_n_200$Elapsed_Time_millisec = 1000 * metrics_without_covariates_n_200$Elapsed_Time_sec
pdf("../figures/simulation_time_and_memory.pdf", height=5, width=4.5)
ggplot(metrics_without_covariates_n_200, aes(x = method, y = Elapsed_Time_sec)) +
  geom_boxplot(aes(col = method, shape = method)) +
  ggtitle("Elapsed Time (sec)", subtitle = "Total sample size 200, 10K genes, no covariates") +
  scale_shape_manual(values=1:nlevels(metrics_without_covariates_n_200$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1")) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8))
ggplot(metrics_without_covariates_n_200, aes(x = method, y = Peak_RAM_Used_MiB)) +
  geom_boxplot(aes(col = method, shape = method)) +
  ggtitle("Peak RAM Used (MB)", subtitle = "Total sample size 200, 10K genes, no covariates") +
  scale_shape_manual(values=1:nlevels(metrics_without_covariates_n_200$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1")) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8))
dev.off()
by(metrics_without_covariates_n_200$Elapsed_Time_sec, metrics_without_covariates_n_200$method, summary)
by(metrics_without_covariates_n_200$Peak_RAM_Used_MiB, metrics_without_covariates_n_200$method, summary)



############################################################
# all methods
############################################################
metrics = metrics_all
# metrics = metrics_all[metrics_all$replicate == 1, ]

names(metrics) = names(metrics) %>%
  gsub("power", "DEGenes__PositiveRate", x=.) %>%
  gsub("type_I_error", "NonDEGenes_PositiveRate", x=.)

metrics$method = metrics$method %>%
  gsub("_with_RIN", "", x=.) %>%
  gsub("without_RIN", "w/o_covariates", x=.) %>%
  gsub("_without_", "_w/o_", x=.) %>%
  gsub("_with_", "_w/_", x=.)

pvalue_cutoff = 0.05
fdr_cutoff = 0.1

# Figure 1:
# The cutoff being used is FDR < 0.1.

# only show methods with covariates
metrics_with_covariates = metrics
# metrics_with_covariates = metrics_with_covariates[grepl("ICeDT", metrics_with_covariates$method), ]
# metrics_with_covariates$method = metrics_with_covariates$method %>% gsub("_ICeDT", "", x=.)
metrics_with_covariates = metrics_with_covariates[!grepl("csSAM|w/o_covariates", metrics_with_covariates$method), ]
metrics_with_covariates$method = metrics_with_covariates$method %>% gsub("_w/?o_covariates", "", x=.)
metrics_with_covariates$method = as.factor(metrics_with_covariates$method)

# only show methods without covariates
metrics_without_covariates = metrics
# metrics_without_covariates = metrics_without_covariates[grepl("ICeDT", metrics_without_covariates$method), ]
# metrics_without_covariates$method = metrics_without_covariates$method %>% gsub("_ICeDT", "", x=.)
metrics_without_covariates = metrics_without_covariates[grepl("csSAM|w/o_covariates", metrics_without_covariates$method), ]
metrics_without_covariates$method = metrics_without_covariates$method %>% gsub("_w/?o_covariates", "", x=.)
metrics_without_covariates$method = factor(metrics_without_covariates$method,
                                           levels=c(levels(metrics_with_covariates$method), "csSAM", "csSAM_ICeDT"))

g1_with_covariates = ggplot(metrics_with_covariates, aes(x = fdr, y = sensitivity)) +
  geom_point(aes(col = method, shape = method), stroke = 0.2) +
  facet_grid(DE_pattern ~ n) +
  geom_vline(xintercept=fdr_cutoff, linetype="solid", color = "black", size = 0.2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_shape_manual(values=1:nlevels(metrics_with_covariates$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1", n = 8, type = "continuous")) +
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
  scale_color_manual(values=wes_palette(name="Darjeeling1", n = 8, type = "continuous")) +
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

pdf("../figures/compare_methods_fdr_all_methods_draft_version_p1.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 1),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 1),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_all_methods_draft_version_p2.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 2),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 2),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_all_methods_draft_version_p3.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 3),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 3),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_all_methods_draft_version_p4.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 4),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 4),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

############################################################
# TOAST counts vs. TPM 
############################################################
metrics = metrics_all
# metrics = metrics_all[metrics_all$replicate == 1, ]

names(metrics) = names(metrics) %>%
  gsub("power", "DEGenes__PositiveRate", x=.) %>%
  gsub("type_I_error", "NonDEGenes_PositiveRate", x=.)

metrics$method = metrics$method %>%
  gsub("_with_RIN", "", x=.) %>%
  gsub("without_RIN", "w/o_covariates", x=.) %>%
  gsub("_without_", "_w/o_", x=.) %>%
  gsub("_with_", "_w/_", x=.)

pvalue_cutoff = 0.05
fdr_cutoff = 0.1

# Figure 1:
# The cutoff being used is FDR < 0.1.

# only show methods with covariates
metrics_with_covariates = metrics
metrics_with_covariates = metrics_with_covariates[grepl("ICeDT", metrics_with_covariates$method), ]
metrics_with_covariates$method = metrics_with_covariates$method %>% gsub("_ICeDT", "", x=.)
metrics_with_covariates = metrics_with_covariates[!grepl("csSAM|w/o_covariates", metrics_with_covariates$method), ]
metrics_with_covariates$method = metrics_with_covariates$method %>% gsub("_w/?o_covariates", "", x=.)
metrics_with_covariates$method = as.factor(metrics_with_covariates$method)

# only show methods without covariates
metrics_without_covariates = metrics
metrics_without_covariates = metrics_without_covariates[grepl("ICeDT", metrics_without_covariates$method), ]
metrics_without_covariates$method = metrics_without_covariates$method %>% gsub("_ICeDT", "", x=.)
metrics_without_covariates = metrics_without_covariates[grepl("csSAM|w/o_covariates", metrics_without_covariates$method), ]
metrics_without_covariates$method = metrics_without_covariates$method %>% gsub("_w/?o_covariates", "", x=.)
metrics_without_covariates$method = factor(metrics_without_covariates$method,
                                           levels=c(levels(metrics_with_covariates$method), "csSAM", "csSAM_ICeDT"))

g1_with_covariates = ggplot(metrics_with_covariates, aes(x = fdr, y = sensitivity)) +
  geom_point(aes(col = method, shape = method), stroke = 0.2) +
  facet_grid(DE_pattern ~ n) +
  geom_vline(xintercept=fdr_cutoff, linetype="solid", color = "black", size = 0.2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_shape_manual(values=1:nlevels(metrics_with_covariates$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1", n = 8, type = "continuous")) +
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
  scale_color_manual(values=wes_palette(name="Darjeeling1", n = 8, type = "continuous")) +
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

pdf("../figures/compare_methods_fdr_ICeDT_methods_draft_version_p1.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 1),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 1),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_ICeDT_methods_draft_version_p2.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 2),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 2),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_ICeDT_methods_draft_version_p3.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 3),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 3),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

pdf("../figures/compare_methods_fdr_ICeDT_methods_draft_version_p4.pdf", height=7, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 4),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 4),
  ncol = 1, nrow = 2, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL),
  align = "v"
)
dev.off()

############################################################
# Supplementary Figure (combined draft version, w/o ICeDT)
############################################################

metrics = metrics_all
# use true cell fractions only:
metrics = metrics[!grepl("_ICeDT", metrics$method), ]

names(metrics) = names(metrics) %>%
  gsub("power", "DEGenes__PositiveRate", x=.) %>%
  gsub("type_I_error", "NonDEGenes_PositiveRate", x=.)

metrics$method = metrics$method %>%
  gsub("_with_RIN", "", x=.) %>%
  gsub("without_RIN", "w/o_covariates", x=.)
metrics$method = as.factor(metrics$method)

pvalue_cutoff = 0.05
fdr_cutoff = 0.1

# Figure 1:
# The cutoff being used is FDR < 0.1.

# only show methods with covariates
metrics_with_covariates =  metrics[!grepl("csSAM|w/o_covariates", metrics$method), ]
metrics_with_covariates$method = as.factor(metrics_with_covariates$method)

# only show methods without covariates
metrics_without_covariates =  metrics[grepl("csSAM|w/o_covariates", metrics$method), ]
metrics_without_covariates$method = as.factor(metrics_without_covariates$method)


metrics_with_covariates$method = as.factor(sapply(metrics_with_covariates$method, function(x) {strsplit(as.character(x), split="_")[[1]][1]}))
metrics_without_covariates$method = as.factor(sapply(metrics_without_covariates$method, function(x) {strsplit(as.character(x), split="_")[[1]][1]}))
metrics_without_covariates$method = factor(metrics_without_covariates$method, levels=c("CARseq", "TOAST", "csSAM"))
library(ggpubr)
library(wesanderson)
g1_with_covariates = ggplot(metrics_with_covariates, aes(x = fdr, y = sensitivity)) +
  geom_point(aes(col = method, shape = method), stroke = 0.2) +
  facet_grid(DE_pattern ~ n) +
  geom_vline(xintercept=fdr_cutoff, linetype="solid", color = "black", size = 0.2) +
  scale_x_continuous(limits = c(0, 0.25)) +
  scale_shape_manual(values=1:nlevels(metrics_with_covariates$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1")) +
  # scale_shape_manual(values=c(15, 0, 7, 16, 1, 25, 6, 25, 6, 17, 2, 18)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8))
g1_without_covariates = ggplot(metrics_without_covariates, aes(x = fdr, y = sensitivity)) +
  geom_point(aes(col = method, shape = method), stroke = 0.2) +
  facet_grid(DE_pattern ~ n) +
  geom_vline(xintercept=fdr_cutoff, linetype="solid", color = "black", size = 0.2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_shape_manual(values=1:nlevels(metrics_without_covariates$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1")) +
  # scale_shape_manual(values=c(15, 0, 7, 16, 1, 25, 6, 25, 6, 17, 2, 18)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8))

pdf("../figures/compare_methods_fdr_true_cell_fractions_draft_version_p1.pdf", height=4, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 1),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 1),
  ncol = 2, nrow = 1, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL)
)
dev.off()

pdf("../figures/compare_methods_fdr_true_cell_fractions_draft_version_p2.pdf", height=4, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 2),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 2),
  ncol = 2, nrow = 1, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL)
)
dev.off()

pdf("../figures/compare_methods_fdr_true_cell_fractions_draft_version_p3.pdf", height=4, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 3),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 3),
  ncol = 2, nrow = 1, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL)
)
dev.off()

pdf("../figures/compare_methods_fdr_true_cell_fractions_draft_version_p4.pdf", height=4, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 4),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 4),
  ncol = 2, nrow = 1, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL)
)
dev.off()


############################################################
# Figure 1 (combined draft version, use ICeDT)
############################################################

metrics = metrics_all
# use the estimates derived from ICeDT cell fractions only:
metrics = metrics[grepl("_ICeDT", metrics$method), ]
metrics$method = sub("_ICeDT", "", metrics$method)

names(metrics) = names(metrics) %>%
  gsub("power", "DEGenes__PositiveRate", x=.) %>%
  gsub("type_I_error", "NonDEGenes_PositiveRate", x=.)

metrics$method = metrics$method %>%
  gsub("_with_RIN", "", x=.) %>%
  gsub("without_RIN", "w/o_covariates", x=.)
metrics$method = as.factor(metrics$method)

pvalue_cutoff = 0.05
fdr_cutoff = 0.1

# Figure 1:
# The cutoff being used is FDR < 0.1.

# only show methods with covariates
metrics_with_covariates =  metrics[!grepl("csSAM|w/o_covariates", metrics$method), ]
metrics_with_covariates$method = as.factor(metrics_with_covariates$method)

# only show methods without covariates
metrics_without_covariates =  metrics[grepl("csSAM|w/o_covariates", metrics$method), ]
metrics_without_covariates$method = as.factor(metrics_without_covariates$method)


metrics_with_covariates$method = as.factor(sapply(metrics_with_covariates$method, function(x) {strsplit(as.character(x), split="_")[[1]][1]}))
metrics_without_covariates$method = as.factor(sapply(metrics_without_covariates$method, function(x) {strsplit(as.character(x), split="_")[[1]][1]}))
metrics_without_covariates$method = factor(metrics_without_covariates$method, levels=c("CARseq", "TOAST", "csSAM"))
library(ggpubr)
library(wesanderson)
g1_with_covariates = ggplot(metrics_with_covariates, aes(x = fdr, y = sensitivity)) +
  geom_point(aes(col = method, shape = method), stroke = 0.2) +
  facet_grid(DE_pattern ~ n) +
  geom_vline(xintercept=fdr_cutoff, linetype="solid", color = "black", size = 0.2) +
  scale_x_continuous(limits = c(0, 0.25)) +
  scale_shape_manual(values=1:nlevels(metrics_with_covariates$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1")) +
  # scale_shape_manual(values=c(15, 0, 7, 16, 1, 25, 6, 25, 6, 17, 2, 18)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8))
g1_without_covariates = ggplot(metrics_without_covariates, aes(x = fdr, y = sensitivity)) +
  geom_point(aes(col = method, shape = method), stroke = 0.2) +
  facet_grid(DE_pattern ~ n) +
  geom_vline(xintercept=fdr_cutoff, linetype="solid", color = "black", size = 0.2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_shape_manual(values=1:nlevels(metrics_without_covariates$method)) +
  scale_color_manual(values=wes_palette(name="Darjeeling1")) +
  # scale_shape_manual(values=c(15, 0, 7, 16, 1, 25, 6, 25, 6, 17, 2, 18)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8))

pdf("../figures/compare_methods_fdr_multiple_replicates_draft_version_p1.pdf", height=4, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 1),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 1),
  ncol = 2, nrow = 1, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL)
)
dev.off()

pdf("../figures/compare_methods_fdr_multiple_replicates_draft_version_p2.pdf", height=4, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 2),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 2),
  ncol = 2, nrow = 1, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL)
)
dev.off()

pdf("../figures/compare_methods_fdr_multiple_replicates_draft_version_p3.pdf", height=4, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 3),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 3),
  ncol = 2, nrow = 1, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL)
)
dev.off()

pdf("../figures/compare_methods_fdr_multiple_replicates_draft_version_p4.pdf", height=4, width=7)
ggarrange(
  g1_with_covariates %+% subset(metrics_with_covariates, DE_pattern_index == 4),
  g1_without_covariates %+% subset(metrics_without_covariates, DE_pattern_index == 4),
  ncol = 2, nrow = 1, labels = c("with covariates", "without covariates"),
  font.label = list(size = 12, color = "black", face = "plain", family = NULL)
)
dev.off()

################################################################################
# show sensitivity of TOAST vs. CARseq with respect to sample size (figure 2)
################################################################################

# only show methods with covariates
metrics_TOAST = metrics[grepl("TOAST_TPM", metrics$method), ]
# metrics_TOAST_without_covariates = metrics[grepl("TOAST_TPM_w/o_covariates", metrics$method), ]
metrics_CARseq = metrics[grepl("CARseq", metrics$method), ]
# metrics_CARseq_without_covariates = metrics[grepl("CARseq_w/o_covariates", metrics$method), ]
# check if these metric tables are already ranked and are directly comparable:
stopifnot(metrics_TOAST$DE_pattern == metrics_CARseq$DE_pattern)
stopifnot(metrics_TOAST$DE_pattern_index == metrics_CARseq$DE_pattern_index)
stopifnot(metrics_TOAST$replicate == metrics_CARseq$replicate)
metrics_ratio_of_sensitivity = metrics_CARseq
metrics_ratio_of_sensitivity$ratio_of_sensitivity = metrics_TOAST$sensitivity / metrics_CARseq$sensitivity 

metrics_ratio_of_sensitivity$method = factor(metrics_ratio_of_sensitivity$method)  # drop unused levels
levels(metrics_ratio_of_sensitivity$method) = c("w/ covariates", "w/o covariates")
# only w/ covariates
metrics_ratio_of_sensitivity = metrics_ratio_of_sensitivity[metrics_ratio_of_sensitivity$method == "w/ covariates",]
metrics_ratio_of_sensitivity$n = factor(metrics_ratio_of_sensitivity$n)
library(ggpubr)
library(wesanderson)
g_ratio_of_sensitivity = ggplot(metrics_ratio_of_sensitivity, aes(x = n, y = ratio_of_sensitivity)) +
  geom_boxplot(aes(shape = method)) + 
  geom_jitter(size = 0.5, col = "red", position=position_jitter(0.2)) +
  facet_grid(rows = vars(DE_pattern)) +
  scale_y_continuous(limits = c(0, 1)) +
  # scale_color_manual(values=wes_palette(name="Darjeeling1")) +
  # geom_hline(yintercept=1, linetype="solid", color = "black", size = 0.2) +
  theme_bw() +
  ylab("ratio of TOAST's sensitivity to CARseq's") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8))

pdf("../figures/ratio_of_sensitivity_TOAST_over_CARseq_p1.pdf", height=4, width=2)
  g_ratio_of_sensitivity %+% subset(metrics_ratio_of_sensitivity, DE_pattern_index == 1)
dev.off()

pdf("../figures/ratio_of_sensitivity_TOAST_over_CARseq_p2.pdf", height=4, width=2)
  g_ratio_of_sensitivity %+% subset(metrics_ratio_of_sensitivity, DE_pattern_index == 2)
dev.off()

pdf("../figures/ratio_of_sensitivity_TOAST_over_CARseq_p3.pdf", height=4, width=2)
  g_ratio_of_sensitivity %+% subset(metrics_ratio_of_sensitivity, DE_pattern_index == 3)
dev.off()

pdf("../figures/ratio_of_sensitivity_TOAST_over_CARseq_p4.pdf", height=4, width=2)
  g_ratio_of_sensitivity %+% subset(metrics_ratio_of_sensitivity, DE_pattern_index == 4)
dev.off()

# sensitivity ratio stratified by quantiles of gene expression
metrics_list = list()
DE_pattern = "2_1_1"
for (n in c(50, 100, 200)) {
  for (method in c("CARseq_ICeDT", "TOAST_TPM_ICeDT")) {
    for (replicate in 1:10) {
      config = sprintf("n_%d_DE_pattern_%s_replicate_%d", n, DE_pattern, replicate)
      CARseq_ICeDT_res = new.env()
      TOAST_TPM_ICeDT_res = new.env()
      if (method == "CARseq_ICeDT") {
        tryCatch(load(file=file.path("..", "simulation", config, "CARseq_ICeDT_res.RData"), envir = CARseq_ICeDT_res), error=function(e) {print(e)})
        metrics_list[[1+length(metrics_list)]] = c(list(gene_expression = "low", 
                                                        n = n, DE_pattern = DE_pattern, DE_pattern_index = 1, replicate = replicate, method = method), 
                                                        compute_metrics_from_pval_matrix(CARseq_ICeDT_res$with_RIN, method, config, c(0.0, 1/3)))
        metrics_list[[1+length(metrics_list)]] = c(list(gene_expression = "medium", 
                                                        n = n, DE_pattern = DE_pattern, DE_pattern_index = 1, replicate = replicate, method = method), 
                                                        compute_metrics_from_pval_matrix(CARseq_ICeDT_res$with_RIN, method, config, c(1/3, 2/3)))
        metrics_list[[1+length(metrics_list)]] = c(list(gene_expression = "high", 
                                                        n = n, DE_pattern = DE_pattern, DE_pattern_index = 1, replicate = replicate, method = method),  
                                                        compute_metrics_from_pval_matrix(CARseq_ICeDT_res$with_RIN, method, config, c(2/3, 1.0)))
      } else if (method == "TOAST_TPM_ICeDT") {
        tryCatch(load(file=file.path("..", "simulation", config, "TOAST_TPM_ICeDT_res.RData"), envir = TOAST_TPM_ICeDT_res), error=function(e) {print(e)})
        metrics_list[[1+length(metrics_list)]] = c(list(gene_expression = "low", 
                                                        n = n, DE_pattern = DE_pattern, DE_pattern_index = 1, replicate = replicate, method = method),  
                                                        compute_metrics_from_pval_matrix(TOAST_TPM_ICeDT_res$with_RIN, method, config, c(0.0, 1/3)))
        metrics_list[[1+length(metrics_list)]] = c(list(gene_expression = "medium", 
                                                        n = n, DE_pattern = DE_pattern, DE_pattern_index = 1, replicate = replicate, method = method),  
                                                        compute_metrics_from_pval_matrix(TOAST_TPM_ICeDT_res$with_RIN, method, config, c(1/3, 2/3)))
        metrics_list[[1+length(metrics_list)]] = c(list(gene_expression = "high", 
                                                        n = n, DE_pattern = DE_pattern, DE_pattern_index = 1, replicate = replicate, method = method),  
                                                        compute_metrics_from_pval_matrix(TOAST_TPM_ICeDT_res$with_RIN, method, config, c(2/3, 1.0)))
      }
    }
  }
}
metrics = do.call(rbind, metrics_list)
write.csv(metrics, file = "../results/compare_methods_metrics_stratified_by_gene_expression.csv", quote = FALSE, row.names = FALSE)

metrics = read.csv("../results/compare_methods_metrics_stratified_by_gene_expression.csv", as.is = TRUE)
# only show methods with covariates
metrics_TOAST = metrics[grepl("TOAST_TPM_ICeDT", metrics$method), ]
# metrics_TOAST_without_covariates = metrics[grepl("TOAST_TPM_w/o_covariates", metrics$method), ]
metrics_CARseq = metrics[grepl("CARseq_ICeDT", metrics$method), ]
# metrics_CARseq_without_covariates = metrics[grepl("CARseq_w/o_covariates", metrics$method), ]
# check if these metric tables are already ranked and are directly comparable:
stopifnot(metrics_TOAST$DE_pattern == metrics_CARseq$DE_pattern)
stopifnot(metrics_TOAST$DE_pattern_index == metrics_CARseq$DE_pattern_index)
stopifnot(metrics_TOAST$replicate == metrics_CARseq$replicate)
metrics_ratio_of_sensitivity = metrics_CARseq
metrics_ratio_of_sensitivity$ratio_of_sensitivity = metrics_TOAST$sensitivity / metrics_CARseq$sensitivity 

metrics_ratio_of_sensitivity$method = factor(metrics_ratio_of_sensitivity$method)  # drop unused levels
levels(metrics_ratio_of_sensitivity$method) = c("w/ covariates", "w/o covariates")
# only w/ covariates
metrics_ratio_of_sensitivity = metrics_ratio_of_sensitivity[metrics_ratio_of_sensitivity$method == "w/ covariates",]
metrics_ratio_of_sensitivity$n = factor(metrics_ratio_of_sensitivity$n)
metrics_ratio_of_sensitivity$gene_expression = factor(metrics_ratio_of_sensitivity$gene_expression, 
                                                      levels = c("low", "medium", "high"))
library(ggpubr)
library(wesanderson)
g_ratio_of_sensitivity = ggplot(metrics_ratio_of_sensitivity, aes(x = gene_expression, y = ratio_of_sensitivity)) +
  geom_boxplot(aes(shape = method)) + 
  geom_jitter(size = 0.5, col = "red", position=position_jitter(0.2)) +
  facet_grid(DE_pattern ~ n) +
  scale_y_continuous(limits = c(0, 1)) +
  # scale_color_manual(values=wes_palette(name="Darjeeling1")) +
  # geom_hline(yintercept=1, linetype="solid", color = "black", size = 0.2) +
  theme_bw() +
  ylab("ratio of TOAST's sensitivity to CARseq's") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = .35, r = .1, b = 0, l = .1, unit = "in"),
        legend.position = "bottom",
        axis.text.x = element_text(size=8, angle=90),
        axis.text.y = element_text(size=8))

pdf("../figures/ratio_of_sensitivity_TOAST_over_CARseq_p1_stratified.pdf", height=3.5, width=6)
g_ratio_of_sensitivity %+% subset(metrics_ratio_of_sensitivity, DE_pattern_index == 1)
dev.off()    



# show how much noise was incorporated in the examples when using the ICeDT estimate
# # rho
# names(csSAM_ICeDT_res$rho_summary[[3]]) =
#     paste0("RMSE_rho_", seq_along(csSAM_ICeDT_res$rho_summary[[1]]))
# # icedt_rho
# names(csSAM_ICeDT_res$rho_summary[[6]]) =
#     paste0("RMSE_rho_from_TPM_", seq_along(csSAM_ICeDT_res$rho_summary[[2]]))
envir1 = new.env()
load("../simulation/n_50_DE_pattern_2_1_1_replicate_1/csSAM_ICeDT_res.RData", envir=envir1)
pdf("../figures/cell_fraction_ICeDT_in_simulation.pdf", width=5.5, height=4)
opar = par(mfrow=c(2, 3), mar=c(4.4, 4.2, 1, 0.8))
H = 6
# y axis: cell fraction from ICeD-T adjusted by "empirical" cell size factors
for (h in 1:H) {
  plot(envir1$rho_summary[[1]][,h], envir1$rho_summary[[2]][,h],
       cex=0.25, col=rgb(0,0,0,0.5),
       xlab="true cell fraction", ylab="cell fraction from ICeD-T (adj.)",
       main=paste("cell type", h), cex.main=1)
  abline(0, 1, lty="dashed", col=rgb(0.7,0.7,0.7,0.5))
}
dev.off()


# Show precision-recall curve
library(PRROC)
method_colors = wes_palette(name="Darjeeling1", 2)
pdf("../figures/precision_recall.pdf", width=3.5, height=6.3)
opar = par(mfrow=c(3, 2), mar=c(5.2, 4.2, 2.8, 1.2))
for (DE_pattern in c("2_1_1", "3_1_1", "4_1_1")) {
  PR_AUCs = as.table(matrix(NA, nrow=2, ncol=3, dimnames = list(rownames=c("CARseq", "TOAST_TPM"), colnames=c(50, 100, 200))))
  for (method in c("CARseq", "TOAST_TPM")) {
    for (n in c(50, 100, 200)) {
      config = sprintf("n_%d_DE_pattern_%s_replicate_1", n, DE_pattern)
      envir = new.env()
      load(file.path("..", "simulation", config, paste0(method, "_ICeDT_res.RData")), envir=envir)
      fold_changes = as.numeric(unlist(strsplit(gsub(".*_pattern_|_replicate_.*", "", config), "_")))[c(3, 1, 2, 3, 3, 3)]
      has_fold_changes = which(fold_changes != 1)
      gene_indices_DE = 1:2000
      gene_indices_non_DE = 2001:10000
      H = 6
      has_fold_changes_matrix = matrix(0, nrow = length(gene_indices_DE) + length(gene_indices_non_DE), ncol = H)
      has_fold_changes_matrix[gene_indices_DE, has_fold_changes] = 1
      
      scores = NULL 
      if (grepl("csSAM", method)) {
        scores = -envir$fdr_matrix
      } else {
        scores = -envir$with_RIN$p
      }
      scores[!is.finite(scores)] = Inf
      roc_pval = pr.curve(scores.class0  = scores,
                          weights.class0 = has_fold_changes_matrix,
                          curve = TRUE)
      PR_AUCs[method, as.character(n)] = roc_pval$auc.integral
      if (grepl("CARseq", method)) {
        roc_pval_CARseq = roc_pval
      } else if (grepl("TOAST", method)) {
        roc_pval_TOAST = roc_pval
      }
  
    }
  }
  # for each DE pattern:
  # left: barplot
  stopifnot(n == 200)  # make sure we are plotting the case when n is 200
  rownames(PR_AUCs) = c("CARseq", "TOAST")
  # The first row has some more legends:
  legend = main_barplot = main_curve = NULL
  if (DE_pattern == "2_1_1") {
    legend = rownames(PR_AUCs)
    main_barplot = "Precision-Recall AUC"
    main_curve = "Precision-Recall curve"
  }
  xx = barplot(PR_AUCs,
          xlab = "n", ylab = "AUC", 
          col = method_colors, ylim = c(0, 1),
          sub = sprintf("%s", DE_pattern),
          legend = legend, beside = TRUE,
          args.legend=list(x="topleft"),
          main = main_barplot, cex.main = 0.9)
  ## Add text at top of bars
  text(x = xx, y = as.numeric(PR_AUCs),
       label = sprintf("  %.2f", as.numeric(PR_AUCs)),
       pos = 3, cex = 0.7, offset = 0.1,
       col = rep(method_colors, 3))
  # right: two PR curves
  plot(NA, col = method_colors[1],
       sub = sprintf("%s, n = %d", DE_pattern, n),
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Recall", ylab = "Precision",
       main = main_curve, cex.main = 0.9)
  plot(roc_pval_CARseq, col = method_colors[1], add = TRUE, lwd = 1)
  plot(roc_pval_TOAST,  col = method_colors[2], add = TRUE, lwd = 1)
}
dev.off()

