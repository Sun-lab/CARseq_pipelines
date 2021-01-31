# estimate cell type proportion of SKCM samples

library(data.table)

# ------------------------------------------------------------------------
# read in bulk RNA-seq data
# ------------------------------------------------------------------------

counts = readRDS("../data/TCGA_SKCM_raw_counts.rds")
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

table(GeneLengths.Mat$Gencode.ID == rownames(counts))

geneLength = GeneLengths.Mat$Exonic
SKCM_TPM   = calculate_TPM(counts, geneLength) 
dim(SKCM_TPM)
SKCM_TPM[1:2,1:5]

geneinfo = fread("../data/gencode.v22.genes.txt", header = TRUE, 
                 sep = "\t", drop = 8, na.strings = c("NA", ""))
head(geneinfo)

# ------------------------------------------------------------------------
# filter genes based on gene annoation
# ------------------------------------------------------------------------

missing.id = which(is.na(geneinfo$hgnc_symbol))
length(missing.id)
geneinfo.nomissing = geneinfo[-missing.id,]

# Some ensembl IDs have duplicated hgnc symbols
duplicated.id = which(duplicated(geneinfo.nomissing$hgnc_symbol))
length(duplicated.id)
geneinfo.filtered = geneinfo.nomissing[-duplicated.id,]
rm(geneinfo.nomissing)

dim(geneinfo.filtered)
geneinfo.filtered[1:2.]
table(geneinfo.filtered$geneId %in% row.names(SKCM_TPM))

# Update SKCM_TPM to only genes with 1-1 mapping between 
# their ensembl ID and hgnc symbols, 35295 genes remained in the end. 
w2match = match(geneinfo.filtered$geneId, row.names(SKCM_TPM))
SKCM_TPM_filtered = data.frame(geneSymbol = geneinfo.filtered$hgnc_symbol) 
SKCM_TPM_filtered = cbind(SKCM_TPM_filtered, SKCM_TPM[w2match,])
dim(SKCM_TPM_filtered)
SKCM_TPM_filtered[1:5,1:5]

rm(SKCM_TPM)

# ------------------------------------------------------------------------
# read in reference
# ------------------------------------------------------------------------

load("../../EPIC/data/TRef.rda")
length(TRef)
lapply(TRef, function(x){if(is.vector(x)) length(x) else dim(x)})
TRef$refProfiles[1:2,1:4]
TRef$sigGenes[1:3]

table(TRef$sigGenes %in% rownames(TRef$refProfiles))
EPIC_ref = TRef$refProfiles[match(TRef$sigGenes, rownames(TRef$refProfiles)),]
dim(EPIC_ref)
EPIC_ref[1:2,]

EPIC_sigmat = data.frame(geneSymbol = rownames(EPIC_ref))
EPIC_sigmat = cbind(EPIC_sigmat, EPIC_ref)
dim(EPIC_sigmat)
EPIC_sigmat[1:2,]

# ------------------------------------------------------------------------
# take intersection
# ------------------------------------------------------------------------

genes2use  = intersect(EPIC_sigmat$geneSymbol, SKCM_TPM_filtered$geneSymbol)
length(genes2use)
match.TCGA = match(genes2use, SKCM_TPM_filtered$geneSymbol)
match.EPIC = match(genes2use, EPIC_sigmat$geneSymbol)

SKCM_TPM_EPIC = SKCM_TPM_filtered[match.TCGA,]
dim(SKCM_TPM_EPIC)
SKCM_TPM_EPIC[1:2,1:5]

EPIC_sigmat = EPIC_sigmat[match.EPIC,] 
dim(EPIC_sigmat)
EPIC_sigmat[1:2,1:5]

table(SKCM_TPM_EPIC$geneSymbol == EPIC_sigmat$geneSymbol)

# ------------------------------------------------------------------------
# write out input for CIBERSORT
# ------------------------------------------------------------------------

CBSTfolder = "../CIBERSORT"
sig_file_EPIC = file.path(CBSTfolder, "signature_gene_EPIC.txt")
mix_file_EPIC = file.path(CBSTfolder, "TCGA_SKCM_TPM_EPIC.txt")

# generate matrices as input for CIBERSORT online
write.table(SKCM_TPM_EPIC, file = mix_file_EPIC, quote = FALSE,
            row.names = FALSE, sep="\t")

write.table(EPIC_sigmat, file = sig_file_EPIC, quote = FALSE,
            row.names = FALSE, sep="\t")


gc()
sessionInfo()

