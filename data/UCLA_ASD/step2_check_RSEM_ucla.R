
library(data.table)

# ----------------------------------------------------------------------
# check RNAseq meta data
# ----------------------------------------------------------------------

# work.dir = "research/data/PsychENCODE/ASD"
# setwd(work.dir)

ucla.rna = read.csv("meta_data/UCLA_R01MH094714_ASD_Metadata_RNAseq.csv",
                    as.is=TRUE)
dim(ucla.rna)
ucla.rna[1:2,]

length(unique(ucla.rna$Individual_ID))
length(unique(ucla.rna$Sample_ID))

table(ucla.rna$BrodmannArea)
table(ucla.rna$BrainRegion)

table(ucla.rna$BrodmannArea, ucla.rna$BrainRegion)

w2kp = which(ucla.rna$BrodmannArea == "ba9")
ucla.rna = ucla.rna[w2kp,]
length(unique(ucla.rna$Individual_ID))
length(unique(ucla.rna$Sample_ID))

t2 = table(ucla.rna$Individual_ID)
t2 = t2[t2 > 1]
t2

ucla.rna[which(ucla.rna$Individual_ID %in% names(t2)),]

table(ucla.rna$Sequencing.Batch)

# ----------------------------------------------------------------------
# check clinical meta data
# ----------------------------------------------------------------------

fnm.ucla.cl = "UCLA_R01MH094714_ASD_Metadata_Clinical_August2016Release.csv"
ucla.clinic = read.csv(file.path("meta_data", fnm.ucla.cl), as.is=TRUE)
dim(ucla.clinic)
ucla.clinic[1:2,]

sapply(c("Grant","StudyName","BrainBank","Organism", "Sex", 
         "Diagnosis", "Agonal.State"), 
       function(xx) table(ucla.clinic[,xx], useNA="ifany"))

length(unique(ucla.clinic$Individual_ID))
table(unique(ucla.rna$Individual_ID) %in% ucla.clinic$Individual_ID)

w2kp = which(ucla.clinic$Individual_ID %in% ucla.rna$Individual_ID)
ucla.clinic = ucla.clinic[w2kp,]
dim(ucla.clinic)

sapply(c("Grant","StudyName","BrainBank","Organism", "Sex", 
         "Diagnosis", "Agonal.State"), 
       function(xx) table(ucla.clinic[,xx], useNA="ifany"))

# ----------------------------------------------------------------------
# read in gene expresion data
# ----------------------------------------------------------------------

fnms = list.files(path="UCLA_ASD", pattern="RSEM_Quant.genes.results", 
                  full.names=TRUE, recursive=TRUE)
length(fnms)
fnms[1:5]

dim(ucla.rna)
ucla.rna[1:2,1:10]

ucla.rna$Sample_ID = gsub("Sample_", "", ucla.rna$Sample_ID)
fnms2find = paste0(ucla.rna$Sample_ID, ".RSEM_Quant.genes.results")

# make sure the resutls are available for all the samples
table(grepl("_ba9", fnms))
table(grepl("_ba9", fnms2find))

table(sapply(fnms2find, function(xx){sum(grepl(xx, x=fnms))}))

# collect gene expression data
for(k in 1:length(fnms2find)){
  wi = grep(fnms2find[k], fnms); stopifnot(length(wi) == 1)
  di = fread(fnms[wi])
  
  if(k == 1){
    eDat = di$expected_count
    genes = di$gene_id
  }else{
    if(nrow(di) != length(genes)){
      stop("# of rows in di does not match with the # of rows of the 1st file\n")
    }
    
    if(any(di$gene_id != genes)){
      stop("gene ids of di do not match with the gene ids of the 1st file\n")
    }
    
    eDat = cbind(eDat, di$expected_count)
  }
}

dim(eDat)
eDat[1:2,1:5]
rownames(eDat) = genes
colnames(eDat) = ucla.rna$Sample_ID

# ----------------------------------------------------------------------
# check again the dupliated samples to decide which one to keep
# keep the one with larger RIN
# ----------------------------------------------------------------------

ucla.rna$Total.Expected.Count = colSums(eDat)

ucla.rna.dup = ucla.rna[which(ucla.rna$Individual_ID %in% names(t2)), ]
ucla.rna.dup = ucla.rna.dup[order(ucla.rna.dup$Individual_ID, ucla.rna.dup$RIN),]
dim(ucla.rna.dup)
ucla.rna.dup[,c(1:10,34)]

sam2rm = ucla.rna.dup$Sample_ID[seq(1,7,by=2)]
sam2rm

w2rm = which(ucla.rna$Sample_ID %in% sam2rm)
w2rm

ucla.rna = ucla.rna[-w2rm,]
table(table(ucla.rna$Individual_ID))

# ----------------------------------------------------------------------
# write out data 
# ----------------------------------------------------------------------

intersect(names(ucla.rna), names(ucla.clinic))
ucla.rna = merge(ucla.rna, ucla.clinic)
dim(ucla.rna)
ucla.rna[1:2,]

saveRDS(ucla.rna, file = "data/ucla_cDat.rds")
saveRDS(eDat, file = "data/ucla_eDat.rds")

# mem_used()
gc()

sessionInfo()
q(save="no")

