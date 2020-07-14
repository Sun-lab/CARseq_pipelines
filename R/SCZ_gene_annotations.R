library(SummarizedExperiment)

input_folder = file.path("..", "data")
results_folder = file.path("..", "results")
figures_folder = file.path("..", "figures")
if (!file.exists(results_folder)) {
  dir.create(results_folder)
}
if (!file.exists(figures_folder)) {
  dir.create(figures_folder, recursive = TRUE)
}

# input
rse_filtered_file = file.path(input_folder, "rse_filtered_SV.rds")
if (file.exists(rse_filtered_file)) {
  rse_filtered = readRDS(rse_filtered_file)
} else {
  # NOTE: see step3_test_methods.R for a fuller version of processing data.
  
  # read files provided by Wei Sun
  CMC_count = readRDS("../data/trec_filtered_scz_control.rds")  # 20788x527 matrix
  clinical_variables_raw = readRDS("../data/dat_cavariates_scz_control_with_svs.rds")    # 527x33 matrix

  # for code compatibility:
  CMC_clinical_merged2 = as.data.frame(clinical_variables_raw)
  CMC_clinical_merged2$Dx = ifelse(CMC_clinical_merged2$DxSCZ == 1, "SCZ", "Control")
  CMC_clinical_merged2$RNAseq.Sample_RNA_ID = colnames(CMC_count)
  # rescale age_death, PMI, RIN, RIN^2, log_depth
  CMC_clinical_merged2$scaled_age_death = scale(CMC_clinical_merged2$age_death)
  CMC_clinical_merged2$scaled_log_depth = scale(CMC_clinical_merged2$log_depth)
  CMC_clinical_merged2$scaled_PMI = scale(CMC_clinical_merged2$PMI)
  CMC_clinical_merged2$scaled_RIN = scale(CMC_clinical_merged2$RIN)
  CMC_clinical_merged2$scaled_RIN2 = scale(CMC_clinical_merged2$RIN2)
  # rescale PCs and SVs
  covariates_to_scale = grep("sv|PC", colnames(CMC_clinical_merged2))
  for (m in covariates_to_scale) {
    CMC_clinical_merged2[, m] = scale(CMC_clinical_merged2[, m])
  }
  rownames(CMC_clinical_merged2) = CMC_clinical_merged2$RNAseq.Sample_RNA_ID
  
  disease = CMC_clinical_merged2$Dx
  table(disease, useNA="ifany")
  
  # proportion estimates from either CIBERSORT or ICeDT
  cellular_frequencies = readRDS("../data/SCZ_prop.rds")
  
  H = ncol(cellular_frequencies[[1]])
  n_B = nrow(cellular_frequencies[[1]])
  
  library(SummarizedExperiment)
  rse_filtered = SummarizedExperiment(assays=SimpleList(couSCZ_gene_annotations.R:nts=CMC_count),
                                      colData=DataFrame(CMC_clinical_merged2))
  rse_filtered = as(rse_filtered, "RangedSummarizedExperiment")
  
  # only use about 15000 ~ 20000 highly expressed genes:
  rd75 = apply(CMC_count, 1, function(x) quantile(x, 0.75))
  
  # fetch gene names
  # obtain gene length from gtf (used in bulk expression quantification):
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
    gencode_gtf = gencode_gtf[match(rownames(rse_filtered), gencode_gtf$gene_id),
                              c("gene_id", "source", "gene_biotype", "gene_name")]
    gencode_gtf$gene_length = as.numeric(exonic_gene_sizes)
    
    save(file = genelength_file, exonic_gene_sizes, gencode_gtf)
  } else {
    load(genelength_file)
  }
  
  # Calculate TPM and save to file
  stopifnot(all(row.names(rse_filtered) %in% gencode_gtf$gene_id))
  gencode_gtf = gencode_gtf[match(row.names(rse_filtered), gencode_gtf$gene_id), ]
  gencode_gtf$rd75 = rd75
  rowData(rse_filtered) = 
    gencode_gtf[, c("gene_id", "source", "gene_biotype", "gene_name", "gene_length", "rd75")]
  # WARNING: the GRanges are only one exon of a gene, so it will not match gene_length (which is longer). 
  TPM = assays(rse_filtered)$counts / rowData(rse_filtered)$gene_length
  assays(rse_filtered)$TPM = t(t(TPM)*1e6/colSums(TPM))
  
  saveRDS(rse_filtered, rse_filtered_file)
}

