  
# ----------------------------------------------------------------------
# SVA: Calculate surrogate variables
# ----------------------------------------------------------------------

rm(list=ls())

repo_dir = ".."
cmc_dir  = ".."

# ----------------------------------------------------------------------
# Library/Source Basic Functions
# ----------------------------------------------------------------------

source(file.path(repo_dir,"R","base_functions_plittle.R"))
library(data.table)
library(sva)
library(matrixStats)

# ----------------------------------------------------------------------
# read in gene expression and covariate data
# ----------------------------------------------------------------------

trec = readRDS(file=file.path(cmc_dir, "data/trec_filtered_log_rd_corrected.rds"))
dim(trec)
trec[1:2,1:5]

dat = readRDS(file=file.path(cmc_dir, "data/dat_cavariates.rds"))
dim(dat)
dat[1:2,]

CMC_count = readRDS("../data/trec_filtered_scz_control.rds")  # 20788x527 matrix
dim(CMC_count)
CMC_count[1:2,1:5]

table(dat$Brain_Region)
table(dat$Dx)

table(dat$RNAseq_sample_id == colnames(trec))

w2kp = which(dat$Dx %in% c("Control", "SCZ"))
length(w2kp)

dat = dat[w2kp,]
trec = trec[,w2kp]

dim(dat)
dat[1:2,1:5]

dim(trec)
trec[1:2,1:5]

# ----------------------------------------------------------------------
# generate cell type proportions
# ----------------------------------------------------------------------

# fetch gene names
# obtain gene length from gtf (used in bulk expression quantification):
# NOTE: the actual genelength_file was generated using 56,632 genes, while we only need a subset of 20,788 genes now.
input_folder = "../data"
data_folder = "../data"
gencode_file = file.path(input_folder, "Homo_sapiens.GRCh37.70.processed.gtf.gz")
genelength_file = file.path(input_folder, "ExonicGeneLengths_GRCh37_70.RData")
txdb_file = file.path(input_folder, "Homo_sapiens_ensembl_70_GRCh37.sqlite")
gtf_link = "Homo_sapiens.GRCh37.70.processed.gtf"
library(GenomicFeatures)
if( !file.exists(genelength_file) ){
  if( !file.exists(txdb_file) ){
    # prepare txdb
    exdb = GenomicFeatures::makeTxDbFromGFF(file = gencode_file,
                                            format="gtf", dataSource = gtf_link)
  } else {
    exdb = AnnotationDbi::loadDb(txdb_file)
  }
  exons.list.per.gene = GenomicFeatures::exonsBy(exdb,by="gene")
  exonic_gene_sizes = lapply(exons.list.per.gene, function(x){sum(width(reduce(x)))})
  
  # We also need to extract other information from gtf files since they are not included in txdb:
  gencode_gtf = rtracklayer::import(gencode_file)
  gencode_gtf = gencode_gtf[!duplicated(gencode_gtf$gene_id), ]
  gencode_gtf = gencode_gtf[match(rownames(trec), gencode_gtf$gene_id),
                            c("gene_id", "source", "gene_biotype", "gene_name")]
  gencode_gtf$gene_length = as.numeric(exonic_gene_sizes)
  
  save(file = genelength_file, exonic_gene_sizes, gencode_gtf)
} else {
  load(genelength_file)
  gencode_gtf = gencode_gtf[match(rownames(trec), gencode_gtf$gene_id),
                            c("gene_id", "source", "gene_biotype", "gene_name", "gene_length")]
}

# Calculate TPM
# indices: gene indices (or names matching rownames(count)) to use to scale 
#          so that the sum of TPM among them is 1 million.
calculate_TPM = function(count, gene_length, indices) {
  if (nrow(count) != length(gene_length)) stop("Number of rows of the count matrix does not match gene lengths!")
  TPM = count / gene_length
  t(t(TPM)*1e6/colSums(TPM[indices, ]))
}

cell_sizes = readRDS("../MTG/cell_sizes_MTG.rds")

# read the signature matrix and prepare the CIBERSORT input
signature_matrix_file = "../MTG/signature_MTG.rds"
signature_matrix = readRDS(signature_matrix_file)$SIG
signature_gene_file = file.path(data_folder, "CIBERSORT_input_signature_gene_SCZ.txt")
mixture_file = file.path(data_folder, "CIBERSORT_input_observed_TPM_SCZ.txt")
# prepare CIBERSORT input
gene_name_match = pmatch(row.names(signature_matrix), gencode_gtf$gene_name)
stopifnot(length(unique(na.omit(gene_name_match))) == length(na.omit(gene_name_match)))
# use the subset of matched genes
signature_matrix = signature_matrix[!is.na(gene_name_match), ]

# requires some rescaling to work properly
# To get proper scaling, we need to obtain a matrix of reference expression
# across all the genes, normalize to TPM for the 15k genes, and then
# take a subset restricted to the signature genes.
# The mixture expression in TPM also needs to be prepared in the same fashion.
signature_gene_names = rownames(signature_matrix)
length(signature_gene_names)
signature_gene_names = intersect(signature_gene_names, gencode_gtf$gene_name)
reference_expression = readRDS("../MTG/all_genes_MTG.rds")
reference_TPM = calculate_TPM(reference_expression$SIG, reference_expression$anno$gene_length, cell_sizes$gene_names)
SCZ_TPM = calculate_TPM(CMC_count, gencode_gtf$gene_length, match(cell_sizes$gene_names, gencode_gtf$gene_name))
signature_matrix_scaled = reference_TPM[signature_gene_names, ]
mixture_TPM_scaled = SCZ_TPM[match(signature_gene_names, gencode_gtf$gene_name), ]
rownames(mixture_TPM_scaled) = signature_gene_names
geneSymbol_and_observed_TPM = cbind(signature_gene_names, mixture_TPM_scaled)
colnames(geneSymbol_and_observed_TPM)[1] = "geneSymbol"
geneSymbol_and_signature_gene_TPM = cbind(signature_gene_names, signature_matrix_scaled)
colnames(geneSymbol_and_signature_gene_TPM)[1] = "geneSymbol"
# generate matrices as input for CIBERSORT online
write.table(geneSymbol_and_observed_TPM,
            file = mixture_file,
            quote = FALSE,
            row.names = FALSE,
            sep="\t")
write.table(geneSymbol_and_signature_gene_TPM,
            file = signature_gene_file,
            quote = FALSE,
            row.names = FALSE,
            sep="\t")

# run CIBERSORT
# https://cibersort.stanford.edu/runcibersort.php
cibersort_output_file = file.path(data_folder, "SCZ_CIBERSORT.Output.csv")
cibersort_output = read.csv(cibersort_output_file)

# run ICeDT
signature_matrix = as.matrix(signature_matrix)
icedt_output_file =  file.path(data_folder, "SCZ_ICeDT_output.rds")
if (!file.exists(icedt_output_file)) {
  set.seed(1234)
  icedt_output = ICeDT::ICeDT(
    Y = mixture_TPM_scaled,
    Z = as.matrix(signature_matrix_scaled),
    tumorPurity = rep(0, ncol(mixture_TPM_scaled)),
    refVar = NULL)
  # save to ICeDT cellular frequency RDS file
  saveRDS(icedt_output, icedt_output_file)
} else {
  icedt_output = readRDS(icedt_output_file)
}

cell_sizes = readRDS("../MTG/cell_sizes_MTG.rds")
cell_sizes = cell_sizes$cell_sizes[c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")]
icedt_rho = t(apply(t(icedt_output$rho)[,c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")],1,function(xx){yy = xx / cell_sizes; yy / sum(yy)}))
cibersort_rho = t(apply(cibersort_output[,c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")],1,function(xx){yy = xx / cell_sizes; yy / sum(yy)}))
rownames(icedt_rho) = rownames(cibersort_rho) = rownames(cibersort_output) = rownames(t(icedt_output$rho))
prop_output_file = file.path(data_folder, "SCZ_prop.rds")
prop_list = list(ICeDT=icedt_rho, CIBERSORT=cibersort_rho)
saveRDS(prop_list, file=prop_output_file)
prop_from_TPM_output_file = file.path(data_folder, "SCZ_prop_from_TPM.rds")
prop_from_TPM_list = list(ICeDT=t(icedt_output$rho)[,c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")],
                          CIBERSORT=cibersort_output[,c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")])
saveRDS(prop_from_TPM_list, file=prop_from_TPM_output_file)

# Note that Chong's MTG SCZ cell fraction estimates "data/SCZ_prop.rds" differs from Paul's "MTG/prop_MTG.rds".
# compare new (by Chong) vs. old (by Paul) SCZ MTG cell fraction estimates
prop = readRDS("data/SCZ_prop.rds")
prop_old = readRDS("MTG/prop_MTG.rds")
prop_old$ICeDT = prop_old$ICeDT[rownames(prop$ICeDT), ]
prop_old$CIBERSORT = prop_old$CIBERSORT[rownames(prop$ICeDT), ]
pdf("figures/SCZ_cell_fraction_estimates_comparison.pdf", height = 6, width = 8)
opar = par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
for(k in seq_len(ncol(prop$ICeDT))) {
  plot(prop_old$ICeDT[,k], prop$ICeDT[,k], main=paste(colnames(prop$ICeDT)[k], "ICeDT"),
       xlab="Chong's est. in data/SCZ_prop.rds", ylab = "Paul's est. in MTG/prop_MTG.rds")
  abline(a=0, b=1)
}
for(k in seq_len(ncol(prop$ICeDT))) {
  plot(prop_old$CIBERSORT[,k], prop$CIBERSORT[,k], main=paste(colnames(prop$ICeDT)[k], "CIBERSORT"),
       xlab="Chong's est. in data/SCZ_prop.rds", ylab = "Paul's est. in MTG/prop_MTG.rds")
  abline(a=0, b=1)
}
par(opar)
dev.off()

# load the cell fraction estimates and use ICeDT results
prop = readRDS(prop_output_file)
names(prop)
prop = prop$ICeDT
dim(prop); head(prop)

table(dat$RNAseq_sample_id %in% rownames(prop))

prop = prop[match(dat$RNAseq_sample_id, rownames(prop)),]
dim(prop); head(prop)

# ----------------------------------------------------------------------
# first run PCA on residuals to help determin the number of PCs to use
# ----------------------------------------------------------------------

ctypes2use = which(colnames(prop) != "Exc")
plog = log(prop[,ctypes2use] + 1e-5)
plog = plog - log(prop[,which(colnames(prop) == "Exc")])
colnames(plog) = paste("log", colnames(plog), sep="_")
dat = cbind(dat, plog)
dim(dat)
dat[1:2,]

table(dat$gender, dat$Sex)

mod0.terms = c("gender","Institution","libclust","age_death","PMI",
               "RIN","RIN2", paste0("genoPC",1:5),"log_depth", 
               colnames(plog))
length(mod0.terms)

mod0.str = paste(mod0.terms, collapse = " + ")

# Looping over genes for lm() and residuals
GG = nrow(trec); NN = ncol(trec)
rr = matrix(NA,GG,NN) # residual matrix

for(gg in seq(GG)){
  # gg = 1
  if(gg %% 1e2 == 0) cat(".")
  if(gg %% 2e3 == 0 || gg == GG) cat(sprintf("%s out of %s\n",gg,GG))
  
  trec.gg = trec[gg,]
  
  lm_out = lm(formula(sprintf("trec.gg ~ %s", mod0.str)), data = dat)
  rr[gg,] = as.numeric(lm_out$residuals)
  aa = drop1(lm_out,.~.,test="F")
  
  if(gg == 1){
    pp = matrix(NA,GG,nrow(aa)-1)
    colnames(pp) = rownames(aa)[-1]
  }
  
  pp[gg,] = aa[-1,6]
  rm(lm_out)
}

pdf(file.path(cmc_dir,"figures/SCZ_expression_PCA_all_log_Prop.pdf"),
    height = 8,width = 8)

show_pvalue_hist(mat_pvalues = pp, test_type0 = 3)

rr2 = rr - rowMeans(rr,na.rm = TRUE)
cov_rr2 = t(rr2) %*% rr2 / nrow(rr2); pca_rr2 = eigen(cov_rr2)

show_screeplot(pca_rr2,main = "")
par(mfrow = c(2,1),mar = c(4,4,1,1),oma = c(0,0,2,0))
barplot(pca_rr2$values[3:21], names.arg =3:21, main = "", 
        xlab = "Index", ylab = "Eigen-value")
barplot(diff(-pca_rr2$values[3:22]), names.arg =3:21, main = "", 
        xlab = "Index", ylab = "Eigen-value[i] - Eigen-value[i+1]")

num_pcs = 5; pcs = smart_df(pca_rr2$vectors[,1:num_pcs])
names(pcs) = paste0("PC",seq(num_pcs))
show_pc_color(pcs,submain = "")

dev.off()

# ----------------------------------------------------------------------
# SVA
# ----------------------------------------------------------------------

mod.terms = c("Dx", mod0.terms)
length(mod.terms)

mod0 = model.matrix(as.formula(paste("~", paste(mod0.terms, collapse=" + "))), 
                    data=dat)
mod  = model.matrix(as.formula(paste("~", paste(mod.terms, collapse=" + "))), 
                    data=dat)

dim(mod0)
mod0[1:2,]

dim(mod)
mod[1:2,]

n.sv = num.sv(trec,mod,method="leek")
n.sv

n.sv = num.sv(trec,mod,method="be")
n.sv

sv5  = sva(trec, mod, mod0, n.sv=5)
sv12 = sva(trec, mod, mod0, n.sv=12)

# ----------------------------------------------------------------------
# check associatoin between PCs and svs
# ----------------------------------------------------------------------

pcs = pca_rr2$vectors[,1:20]
npc = 1:20
round(cor(sv5$sv, pcs[,1:5]),2)
round(cor(sv5$sv, sv12$sv),2)

r2 = matrix(NA, nrow=12, ncol=length(npc))
for(k in 1:12){
  svk = sv12$sv[,k]
  for(m in 1:length(npc)){
    lmk = summary(lm(svk ~ pcs[,1:npc[m]]))
    r2[k,m] = lmk$r.squared
  }
}

round(r2, 2)

pdf(file.path(cmc_dir,"figures/svs_vs_pcs.pdf"),
    height = 6,width = 6)
par(mar=c(5,4,1,1), bty="n")
for(i in 1:nrow(r2)){
  if(i==1){
    plot(npc, r2[1,], ylim=c(0,0.85), type="l", xlim=c(0,21), 
         xlab="# of PCs", ylab="R2 of surrogate variables explained by PCs")
  }else{
    lines(npc, r2[i,])
  }
}
abline(v=10, lty=2, col="darkred")
abline(v=12, lty=2, col="darkred")
abline(v=14, lty=2, col="darkred")
text(rep(21,12), r2[,20], labels=1:12)
dev.off()

# ----------------------------------------------------------------------
# R squared vs. number of surrogate variables to use
# ----------------------------------------------------------------------
ctypes2use = which(colnames(prop) != "Exc")
plog = log(prop[,ctypes2use] + 1e-5)
plog = plog - log(prop[,which(colnames(prop) == "Exc")])
colnames(plog) = paste("log", colnames(plog), sep="_")
colnames(sv12$sv) = paste0("sv", 1:12)
dat = cbind(cbind(dat, plog), sv12$sv)
dim(dat)
dat[1:2,]
log_trec = trec
final.mod0.terms = mod0.terms
final.mod0.str = paste(final.mod0.terms, collapse = " + ")
max_number_of_SV = 12
final.mod0.str.vector = rep(NA, max_number_of_SV + 1)
final.mod0.str.vector[1] = final.mod0.str
for (number_of_SV in (seq_len(max_number_of_SV))) {
  final.mod0.str.vector[1 + number_of_SV] = sprintf("%s + sv%d", final.mod0.str.vector[number_of_SV], number_of_SV)
}

GG = nrow(log_trec); NN = ncol(log_trec)
rr = matrix(NA,GG,NN) # residual matrix

rsquared = matrix(NA, nrow = GG, ncol = 1 + max_number_of_SV)
for(gg in seq_len(GG)){
  # gg = 1
  if(gg %% 1e2 == 0) cat(".")
  if(gg %% 2e3 == 0 || gg == GG) cat(sprintf("%s out of %s\n",gg,GG))
  
  log_trec.gg = log_trec[gg,]
  
  # starts from 0
  for (number_of_SV in (-1 + seq_len(max_number_of_SV + 1))) {
    lm_out = lm(formula(sprintf("log_trec.gg ~ %s", final.mod0.str.vector[1 + number_of_SV])), data = dat)
    rsquared[gg, 1 + number_of_SV] = summary(lm_out)$r.squared
    rm(lm_out)
  }
}
colnames(rsquared) = paste0("SV",  -1 + seq_len(max_number_of_SV + 1))

saveRDS(rsquared, file.path(cmc_dir, "results/SCZ_rsquared.rds"))

colMedians(rsquared, na.rm=TRUE)
colnames(rsquared) = sprintf("%s\n%.2f", paste0("SV",  -1 + seq_len(max_number_of_SV + 1)), colMedians(rsquared, na.rm=TRUE))

png(file.path(cmc_dir, "figures/SCZ_rsquared.png"),
    width=9, height=6, units="in", res=400)
par(mar = c(6.5, 3.0, 0.5, 1.0))
boxplot(rsquared, sub = paste(strwrap(paste("log_trec ~", final.mod0.str), 110), collapse = "\n"))
dev.off()


# ----------------------------------------------------------------------
# save results
# ----------------------------------------------------------------------

mod[1:2,]
colnames(mod)

modDat = cbind(mod[,2:22], sv12$sv)
colnames(modDat)[22:33] = paste0("sv", 1:12)
dim(modDat)
modDat[1:2,]

covariate.file.name = "data/dat_cavariates_scz_control_with_svs.rds"
saveRDS(modDat, file=file.path(cmc_dir, covariate.file.name))

trec.file.name = "data/trec_filtered_log_rd_corrected_scz_control.rds"
saveRDS(trec, file=file.path(cmc_dir, trec.file.name))


trec0 = readRDS(file=file.path(cmc_dir, "data/trec_filtered.rds"))
dim(trec0)
trec0[1:2,1:5]

table(dat$RNAseq_sample_id %in% colnames(trec0))
trec0 = trec0[,match(dat$RNAseq_sample_id, colnames(trec0))]
dim(trec0)
trec0[1:2,1:5]

trec0.file.name = "data/trec_filtered_scz_control.rds"
saveRDS(trec0, file=file.path(cmc_dir, trec0.file.name))

# mem_used()
gc()

sessionInfo()

