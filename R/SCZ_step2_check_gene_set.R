
library(ggplot2)
library(ggpubr)
library(tidyr)
theme_set(theme_bw())
library(tidyverse)
library(SummarizedExperiment)
library(qusage)
library(wesanderson)

# ------------------------------------------------------------------------
# read in SummarizedExperiment object that contains gene name mapping
# ------------------------------------------------------------------------

rse_filtered = readRDS("../data/rse_filtered_SV.rds")
dim(rse_filtered)
rse_filtered[1:2,1:5]

# ------------------------------------------------------------------------
# read in CARseq results
# ------------------------------------------------------------------------

res = readRDS("../results/SCZ_CARseq_ICeDT_SV2.rds")
dim(res)

names(res)

lfc = res$"shrunken_lfc"
dim(lfc)
lfc[1:2,]

rowData(rse_filtered)[1:2,]
table(rowData(rse_filtered)$gene_id == rownames(lfc))

# ------------------------------------------------------------------------
# read in genes within gene set
# ------------------------------------------------------------------------

react = read.gmt("../data/c2.cp.reactome.v7.1.symbols.gmt")
length(react)
names(react)[1:3]
summary(sapply(react, length))

# ------------------------------------------------------------------------
# illustrate pathways for neuron cells
# ------------------------------------------------------------------------

g1 = grep("NEURONAL_SYSTEM", names(react))
g1
names(react)[g1]

g2 = grep("SYNTHESIS_OF_PIPS_AT_THE_LATE_ENDOSOME_MEMBRANE", names(react))
g2
names(react)[g2]

g3 = grep("RAB_REGULATION_OF_TRAFFICKING", names(react))
g3
names(react)[g3]

g4 = grep("ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION", names(react))
g4
names(react)[g3]

g5 = grep("HSF1_DEPENDENT_TRANSACTIVATION", names(react))
g5
names(react)[g5]

g6 = grep("UNBLOCKING_OF_NMDA_RECEPTORS", names(react))
g6
names(react)[g6]

g7 = grep("NMDA_RECEPTOR_MEDIATED_NEURONAL_TRANSMISSION", names(react))
g7
names(react)[g7]

g1 = c(g1, g2, g3, g4, g5, g6, g7)

dfL = NULL

for(j in g1){
  genes = react[[j]]
  w2  = which(rowData(rse_filtered)$gene_name %in% genes)
  dfj = cbind(rowData(rse_filtered)[w2,], lfc[w2,], res$p[w2,])
  dfj$react_pathway = rep(names(react)[j], length(w2))
  dfL = rbind(dfL, dfj)
}

dim(dfL)
dfL[1:2,]

dfL$react_pathway = gsub("^REACTOME_", "", dfL$react_pathway)
dfL$react_pathway = gsub("_", " ", dfL$react_pathway)
table(dfL$react_pathway)

path1 = "NMDA RECEPTOR MEDIATED NEURONAL TRANSMISSION"
path2 = "NMDA RECEPTOR \nMEDIATED NEURONAL TRANSMISSION             "
dfL$react_pathway = gsub(path1, path2, dfL$react_pathway)

dfL$X = as.character(dfL$X)

dfL = as.data.frame(dfL)
names(dfL) = gsub(":", "_", names(dfL), fixed=TRUE)
dim(dfL)
dfL[1:2,]

tapply(dfL$Exc, dfL$react_pathway, summary)
tapply(dfL$Inh, dfL$react_pathway, summary)

p1 = ggplot(dfL[which(dfL$Exc < 0.05),], 
            aes(x=react_pathway, y=SCZ_vs_Control.Exc)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + 
  geom_jitter(shape=16, col="darkblue", alpha=0.7, 
              position=position_jitter(0.15)) + 
  ylab("log fold change, SCZ vs. control") + xlab("") + 
  labs(title="Excitatory neuron") + 
  geom_hline(yintercept = 0, color = "red", size=1.0)

p2 = ggplot(dfL[which(dfL$Inh < 0.05),], 
            aes(x=react_pathway, y=SCZ_vs_Control.Inh)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + 
  geom_jitter(shape=16, col="darkblue", alpha=0.7, 
              position=position_jitter(0.15)) + 
  ylab("log fold change, SCZ vs. control") + xlab("") + 
  labs(title="Inhibitory neuron") + 
  geom_hline(yintercept = 0, color = "red", size=1.0)

pdf("../figures/SCZ_log_fc_neurons.pdf", width=8, height=6)
ggarrange(p1, p2, ncol = 1, nrow = 2, labels="AUTO")
dev.off()

# grouped boxplot
dfL_neuron = rbind(dfL[which(dfL$Exc < 0.05), c("react_pathway","SCZ_vs_Control.Exc")] %>% 
                     rename(SCZ_vs_Control = SCZ_vs_Control.Exc) %>%
                     add_column(CellType = "Exc"),
                   dfL[which(dfL$Inh < 0.05), c("react_pathway","SCZ_vs_Control.Inh")] %>%
                     rename(SCZ_vs_Control = SCZ_vs_Control.Inh) %>%
                     add_column(CellType = "Inh"))
dfL_neuron$react_pathway = strwrap(dfL_neuron$react_pathway, 25, simplify=FALSE) %>%
sapply(function(x) paste(x, collapse="\n"))

pdf("../figures/SCZ_log_fc_neurons_grouped.pdf", width=5.5, height=5)
ggplot(dfL_neuron, 
       aes(x=react_pathway, y=SCZ_vs_Control, col=CellType)) +
  geom_boxplot(outlier.shape = NA, aes(col = CellType), 
               position = position_dodge(preserve="single")) +
  geom_point(aes(col = CellType), position=position_jitterdodge()) +
  coord_flip() + 
  scale_color_manual(name = "CellType", values = as.character(wes_palette(name="Darjeeling2", 6, type="continuous"))[c(2,3)])+
  ylab("log fold change, SCZ vs. control") + xlab("") +    
  labs(title="Neurons") + 
  geom_hline(yintercept = 0, color = "grey", size=1.0)
dev.off()

# ------------------------------------------------------------------------
# illustrate pathways for glial cells
# ------------------------------------------------------------------------

g1 = grep("ROBOS", names(react))
g1
names(react)[g1]

g2 = grep("INNATE_IMMUNE_SYSTEM", names(react))
g2
names(react)[g2]

g3 = grep("INFLUENZA", names(react))
g3
names(react)[g3]

g4 = which("REACTOME_CELL_CYCLE" == names(react))
g4
names(react)[g4]

g5 = grep("EUKARYOTIC_TRANSLATION_ELONGATION", names(react))
g5
names(react)[g5]

g6 = grep("METABOLISM_OF_RNA", names(react))
g6
names(react)[g6]

g7 = grep("RESPONSE_OF_EIF2AK4_GCN2", names(react))
g7
names(react)[g7]

g1 = c(g1, g2, g3, g4, g5, g6, g7)

dfL = NULL

for(j in g1){
  genes = react[[j]]
  w2  = which(rowData(rse_filtered)$gene_name %in% genes)
  dfj = cbind(rowData(rse_filtered)[w2,], lfc[w2,], res$p[w2,])
  dfj$react_pathway = rep(names(react)[j], length(w2))
  dfL = rbind(dfL, dfj)
}

dim(dfL)
dfL[1:2,]

dfL$react_pathway = gsub("^REACTOME_", "", dfL$react_pathway)
dfL$react_pathway = gsub("_", " ", dfL$react_pathway)
table(dfL$react_pathway)

dfL$X = as.character(dfL$X)

dfL = as.data.frame(dfL)
names(dfL) = gsub(":", "_", names(dfL), fixed=TRUE)
dim(dfL)
dfL[1:2,]

tapply(dfL$Micro, dfL$react_pathway, summary)
tapply(dfL$Oligo, dfL$react_pathway, summary)

p1 = ggplot(dfL[which(dfL$Micro < 0.05),], 
            aes(x=react_pathway, y=SCZ_vs_Control.Micro)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + 
  geom_jitter(shape=16, col="darkblue", alpha=0.7, 
              position=position_jitter(0.15)) + 
  ylab("log fold change, SCZ vs. control") + xlab("") + 
  labs(title="Microglia") + 
  geom_hline(yintercept = 0, color = "red", size=1.0)

p2 = ggplot(dfL[which(dfL$Oligo < 0.05),], 
            aes(x=react_pathway, y=SCZ_vs_Control.Oligo)) + 
  geom_boxplot(outlier.shape = NA) + coord_flip() + 
  geom_jitter(shape=16, col="darkblue", alpha=0.7, 
              position=position_jitter(0.15)) + 
  ylab("log fold change, SCZ vs. control") + xlab("") + 
  labs(title="Oligodendrocyte") + 
  geom_hline(yintercept = 0, color = "red", size=1.0)

pdf("../figures/SCZ_log_fc_glias.pdf", width=7, height=6)
ggarrange(p1, p2, ncol = 1, nrow = 2, labels="AUTO")
dev.off()


gc()

sessionInfo()
q(save="no")

