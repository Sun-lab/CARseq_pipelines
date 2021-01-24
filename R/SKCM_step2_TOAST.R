
library(TOAST)
library(data.table)

# ------------------------------------------------------------------------
# read in bulk RNA-seq data
# ------------------------------------------------------------------------

counts = readRDS("../data/TCGA_SKCM_raw_counts.rds")
dim(counts)
counts[1:2,1:5]

# remove genes with low read counts: 
# genes with 75% percentile of read counts rd75 < 20 are removed
rd75   = apply(counts, 1, function(x) quantile(x, 0.75))
counts = counts[rd75 >= 20, ]
counts = round(counts)
dim(counts)
counts[1:2,1:5]

calculate_TPM = function(count, gene_length) {
  if (nrow(count) != length(gene_length)) {
    stop("Number of rows of the count matrix does not match gene lengths!")
  }
  TPM = count / gene_length
  t(t(TPM)*1e6/colSums(TPM))
}

load("../data/Gene_Lengths.RData")
dim(GeneLengths.Mat)
GeneLengths.Mat[1:2,]

table(rownames(counts) %in% GeneLengths.Mat$Gencode.ID)

mat1 = match(rownames(counts), GeneLengths.Mat$Gencode.ID)
geneLength = GeneLengths.Mat$Exonic[mat1]
SKCM_TPM   = calculate_TPM(counts, geneLength) 
dim(SKCM_TPM)
SKCM_TPM[1:2,1:5]

# ------------------------------------------------------------------------
# read in cell type proportion and covariate data
# ------------------------------------------------------------------------

rho_SKCM = readRDS("../data/SKCM_cell_fraction.rds")
dim(rho_SKCM)
rho_SKCM[1:2,]

col_data = readRDS("../data/SKCM_dat_cavariates_with_svs.rds")
dim(col_data)
col_data[1:2,]

table(rownames(rho_SKCM) == rownames(col_data))
table(rownames(rho_SKCM) %in% colnames(SKCM_TPM))
SKCM_TPM = SKCM_TPM[,match(rownames(rho_SKCM), colnames(SKCM_TPM))]
dim(SKCM_TPM)
SKCM_TPM[1:2,1:3]

class(col_data)
col_data = data.frame(col_data)
dim(col_data)
col_data[1:2,]

design = model.matrix(~ gender + scaled_age + scaled_log_depth + 
                        stageII + stageIII + stageIV + tss + sv1,
                      col_data)[, -1]
design_out = makeDesign(as.data.frame(design), rho_SKCM)
fitted_model = fitModel(design_out, assays(rse_filtered)$TPM)


gc()

sessionInfo()
q(save="no")
