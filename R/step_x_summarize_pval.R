
library(data.table)
library(stringr)

library(ggplot2)
library(ggpubr)
library(ggpointdensity)
theme_set(theme_classic2())
library(ggcorrplot)

# ----------------------------------------------------------------------
# read in p-values for ASD
# ----------------------------------------------------------------------

pval.path = "../results/_pvalues"

pval.files = list.files(path=pval.path, pattern="ASD", full.names=TRUE)
pval.files

f1 = pval.files[1]
pval.ASD = fread(f1)
lb = str_extract(basename(f1), "(?<=ASD_)\\w+_\\w+(?=.txt)")
names(pval.ASD)[3] = lb

dim(pval.ASD)
pval.ASD[1:2,]

length(unique(pval.ASD$gene_id))
length(unique(pval.ASD$gene_name))

for(f1 in pval.files[-1]){
  d1 = fread(f1)
  if(any(d1$gene_id != pval.ASD$gene_id)){
    stop("gene_id do not match\n")
  }
  
  bn = basename(f1)
  lb = str_extract(bn, "(?<=ASD_)\\w+_\\w+(?=.txt)")
  
  pval.ASD[[lb]] = d1$pvalue
}

dim(pval.ASD)
pval.ASD[1:2,]

cor(pval.ASD[,3:8], pval.ASD[,9], use="pairwise.complete.obs")
cor(pval.ASD[,10:15], pval.ASD[,9], use="pairwise.complete.obs")

cr1 = cor(pval.ASD[,3:8], pval.ASD[,10:15], use="pairwise.complete.obs")
round(cr1,3)
summary(c(cr1))
rownames(cr1) = gsub("CARseq_", "", rownames(cr1))
colnames(cr1) = gsub("TOAST_", "", colnames(cr1))

gc1 = ggcorrplot(cr1, lab = TRUE) + 
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       limit=c(-0.2, 0.2))
gc1 = annotate_figure(gc1, bottom = text_grob("CARseq"), 
                      left=text_grob("TOAST", rot = 90))

pdf("../figures/step_x_ASD_CARseq_vs_TOAST_cor.pdf", width=4.5, height=4)
print(gc1)
dev.off()

# ----------------------------------------------------------------------
# read in p-values for SCZ
# ----------------------------------------------------------------------

pval.files = list.files(path=pval.path, pattern="SCZ", full.names=TRUE)
pval.files

f1 = pval.files[1]
pval.SCZ = fread(f1)
lb = str_extract(basename(f1), "(?<=SCZ_)\\w+_\\w+(?=.txt)")
names(pval.SCZ)[3] = lb

dim(pval.SCZ)
pval.SCZ[1:2,]

length(unique(pval.SCZ$gene_id))
length(unique(pval.SCZ$gene_name))

for(f1 in pval.files[-1]){
  d1 = fread(f1)
  if(any(d1$gene_id != pval.SCZ$gene_id)){
    stop("gene_id do not match\n")
  }
  
  bn = basename(f1)
  lb = str_extract(bn, "(?<=SCZ_)\\w+_\\w+(?=.txt)")
  
  pval.SCZ[[lb]] = d1$pvalue
}

dim(pval.SCZ)
pval.SCZ[1:2,]

cor(-log10(pval.SCZ[,3:8]), -log10(pval.SCZ[,9]), 
    use="pairwise.complete.obs")
cor(-log10(pval.SCZ[,10:15]), -log10(pval.SCZ[,9]), 
    use="pairwise.complete.obs")

cr1 = cor(pval.SCZ[,3:8], pval.SCZ[,10:15], use="pairwise.complete.obs")
round(cr1,3)
summary(c(cr1))

rownames(cr1) = gsub("CARseq_", "", rownames(cr1))
colnames(cr1) = gsub("TOAST_", "", colnames(cr1))

gc1 = ggcorrplot(cr1, lab = TRUE) + 
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       limit=c(-0.25, 0.25))
gc1 = annotate_figure(gc1, bottom = text_grob("CARseq"), 
                      left=text_grob("TOAST", rot = 90))

pdf("../figures/step_x_SCZ_CARseq_vs_TOAST_cor.pdf", width=4.5, height=4)
print(gc1)
dev.off()


gs1 = ggplot(pval.SCZ,aes(x=-log10(CARseq_Astro),y=-log10(TOAST_Astro))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval CARseq_Astro)", y= "-log10(pval TOAST_Astro)")

gs2 = ggplot(pval.SCZ,aes(x=-log10(CARseq_Exc),y=-log10(TOAST_Exc))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval CARseq_Exc)", y= "-log10(pval TOAST_Exc)")

gs3 = ggplot(pval.SCZ,aes(x=-log10(CARseq_Inh),y=-log10(TOAST_Inh))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval CARseq_Inh)", y= "-log10(pval TOAST_Inh)")

gs4 = ggplot(pval.SCZ,aes(x=-log10(CARseq_Micro),y=-log10(TOAST_Micro))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval CARseq_Micro)", y= "-log10(pval TOAST_Micro)")

gs5 = ggplot(pval.SCZ,aes(x=-log10(CARseq_Oligo),y=-log10(TOAST_Oligo))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval CARseq_Oligo)", y= "-log10(pval TOAST_Oligo)")

gs6 = ggplot(pval.SCZ,aes(x=-log10(CARseq_OPC),y=-log10(TOAST_OPC))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval CARseq_OPC)", y= "-log10(pval TOAST_OPC)")

pdf(file="../figures/step_x_SCZ_CARseq_vs_TOAST_scatter.pdf", width=13.5, height=7)
ggarrange(gs1, gs2, gs3, gs4, gs5, gs6, 
          ncol = 3, nrow = 2)
dev.off()


# ----------------------------------------------------------------------
# compare p-values between SCZ and ASD
# ----------------------------------------------------------------------

dim(pval.SCZ)
pval.SCZ[1:2,]

dim(pval.ASD)
pval.ASD[1:2,]

pval.ASD$gene_id = str_extract(pval.ASD$gene_id, "\\S+(?=\\.\\d+)")
dim(pval.ASD)
pval.ASD[1:2,]

table(pval.ASD$gene_id %in% pval.SCZ$gene_id)

pvals = merge(pval.SCZ, pval.ASD, by="gene_id", 
              suffixes = c(".SCZ", ".ASD"))
dim(pvals)
pvals[1:2,]

names(pvals)


cr1 = cor(-log10(pvals[,3:8]), -log10(pvals[,17:22]), 
          use="pairwise.complete.obs")

cr1.sp = cor(-log10(pvals[,3:8]), -log10(pvals[,17:22]), 
          use="pairwise.complete.obs", method="spearman")

cr2 = cor(-log10(pvals[,10:15]), -log10(pvals[,24:29]), 
          use="pairwise.complete.obs")

cr2.sp = cor(-log10(pvals[,10:15]), -log10(pvals[,24:29]), 
          use="pairwise.complete.obs", method="spearman")

rownames(cr1) = gsub("CARseq_", "", rownames(cr1))
colnames(cr1) = gsub("CARseq_", "", colnames(cr1))

rownames(cr1.sp) = gsub("CARseq_", "", rownames(cr1.sp))
colnames(cr1.sp) = gsub("CARseq_", "", colnames(cr1.sp))

rownames(cr2) = gsub("TOAST_", "", rownames(cr2))
colnames(cr2) = gsub("TOAST_", "", colnames(cr2))

rownames(cr2.sp) = gsub("TOAST_", "", rownames(cr2.sp))
colnames(cr2.sp) = gsub("TOAST_", "", colnames(cr2.sp))

cor.test(-log10(pvals[,CARseq_Micro.SCZ]), 
         -log10(pvals[,CARseq_Micro.ASD]))

cor.test(-log10(pvals[,TOAST_Micro.SCZ]), 
         -log10(pvals[,TOAST_Micro.ASD]))


cor.test(-log10(pvals[,CARseq_Micro.SCZ]), 
         -log10(pvals[,CARseq_Micro.ASD]), method="spearman")

cor.test(-log10(pvals[,TOAST_Micro.SCZ]), 
         -log10(pvals[,TOAST_Micro.ASD]), method="spearman")

summary(c(cr1))
summary(c(cr2))


round(cr1, 2)
round(cr1.sp, 2)

rownames(cr1) = gsub(".SCZ", "", rownames(cr1), fixed=TRUE)
colnames(cr1) = gsub(".ASD", "", colnames(cr1), fixed=TRUE)

gc1 = ggcorrplot(cr1, lab = TRUE) + 
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       limit=c(-0.2, 0.2))
gc1 = annotate_figure(gc1, bottom = text_grob("Schizophrenia"), 
                      left=text_grob("ASD", rot = 90))

pdf("../figures/step_x_SCZ_vs_ASD_CARseq_cor.pdf", width=4.5, height=4)
print(gc1)
dev.off()


rownames(cr1.sp) = gsub(".SCZ", "", rownames(cr1.sp), fixed=TRUE)
colnames(cr1.sp) = gsub(".ASD", "", colnames(cr1.sp), fixed=TRUE)

gc1 = ggcorrplot(cr1.sp, lab = TRUE) + 
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       limit=c(-0.25, 0.25))
gc1 = annotate_figure(gc1, bottom = text_grob("Schizophrenia"), 
                      left=text_grob("ASD", rot = 90))

pdf("../figures/step_x_SCZ_vs_ASD_CARseq_cor_spearman.pdf", width=4.5, height=4)
print(gc1)
dev.off()


round(cr2, 2)
round(cr2.sp, 2)

rownames(cr2) = gsub(".SCZ", "", rownames(cr2), fixed=TRUE)
colnames(cr2) = gsub(".ASD", "", colnames(cr2), fixed=TRUE)

gc2 = ggcorrplot(cr2, lab = TRUE) +   
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       limit=c(-0.2, 0.2))
gc2 = annotate_figure(gc2, bottom = text_grob("Schizophrenia"), 
                      left=text_grob("ASD", rot = 90))
pdf("../figures/step_x_SCZ_vs_ASD_TOAST_cor.pdf", width=4.5, height=4)
print(gc2)
dev.off()

gs1 = ggplot(pvals,aes(x=-log10(CARseq_Astro.SCZ),y=-log10(CARseq_Astro.ASD))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval Astro) for SCZ", y= "-log10(pval Astro) for ASD")

gs2 = ggplot(pvals,aes(x=-log10(CARseq_Exc.SCZ),y=-log10(CARseq_Exc.ASD))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval Exc) for SCZ", y= "-log10(pval Exc) for ASD")

gs3 = ggplot(pvals,aes(x=-log10(CARseq_Inh.SCZ),y=-log10(CARseq_Inh.ASD))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval Inh) for SCZ", y= "-log10(pval Inh) for ASD")

gs4 = ggplot(pvals,aes(x=-log10(CARseq_Micro.SCZ),y=-log10(CARseq_Micro.ASD))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval Micro) for SCZ", y= "-log10(pval Micro) for ASD")

gs5 = ggplot(pvals,aes(x=-log10(CARseq_Oligo.SCZ),y=-log10(CARseq_Oligo.ASD))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval Oligo) for SCZ", y= "-log10(pval Oligo) for ASD")

gs6 = ggplot(pvals,aes(x=-log10(CARseq_OPC.SCZ),y=-log10(CARseq_OPC.ASD))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "-log10(pval OPC) for SCZ", y= "-log10(pval OPC) for ASD")

pdf(file="../figures/step_x_SCZ_vs_ASD_CARseq_scatter.pdf", width=13.5, height=7)
ggarrange(gs1, gs2, gs3, gs4, gs5, gs6, 
          ncol = 3, nrow = 2)
dev.off()

# ----------------------------------------------------------------------
# compare DESeq2 results between ASD and SCZ
# ----------------------------------------------------------------------

cor.test(-log10(pvals[,DESeq2_bulk.SCZ]), -log10(pvals[,DESeq2_bulk.ASD]))

table(pvals[,DESeq2_bulk.SCZ] < 0.05, pvals[,DESeq2_bulk.ASD] < 0.05)

# ----------------------------------------------------------------------
# check effect sizes
# ----------------------------------------------------------------------

asd = readRDS("../results/ASD_CARseq_ICeDT_seqSV4.rds")
sapply(asd, dim)

scz = readRDS("../results/SCZ_CARseq_ICeDT_SV2.rds")
sapply(scz, dim)

asd$shrunken_lfc[1:2,]
scz$shrunken_lfc[1:2,]

asd$p[1:2,]
pval.ASD[1:2,]

summary(asd$p[,1] - pval.ASD$CARseq_Astro)

asd.lfc = asd$shrunken_lfc
rownames(asd.lfc) = pval.ASD$gene_id

scz.lfc = scz$shrunken_lfc

gene.ids = intersect(rownames(asd.lfc), rownames(scz.lfc))
length(gene.ids)

df1 = data.frame(cbind(asd.lfc[match(gene.ids, rownames(asd.lfc)),],
                       scz.lfc[match(gene.ids, rownames(scz.lfc)),]))
dim(df1)
df1[1:2,]

cr1 = cor(df1[,1:6], df1[7:12], use="pairwise.complete.obs")
round(cr1,3)

gs1 = ggplot(df1,aes(x=SCZ_vs_Control.Astro,y=ASD_vs_Control.Astro)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "lfc Astro for SCZ", y= "lfc Astro for ASD")

gs2 = ggplot(df1,aes(x=SCZ_vs_Control.Exc,y=ASD_vs_Control.Exc)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "lfc Exc for SCZ", y= "lfc Exc for ASD")

gs3 = ggplot(df1,aes(x=SCZ_vs_Control.Inh,y=ASD_vs_Control.Inh)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "lfc Inh for SCZ", y= "lfc Inh for ASD")

gs4 = ggplot(df1,aes(x=SCZ_vs_Control.Micro,y=ASD_vs_Control.Micro)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "lfc Micro for SCZ", y= "lfc Micro for ASD")

gs5 = ggplot(df1,aes(x=SCZ_vs_Control.Oligo,y=ASD_vs_Control.Oligo)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "lfc Oligo for SCZ", y= "lfc Oligo for ASD")

gs6 = ggplot(df1,aes(x=SCZ_vs_Control.OPC,y=ASD_vs_Control.OPC)) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x= "lfc OPC for SCZ", y= "lfc OPC for ASD")

pdf(file="../figures/step_x_SCZ_vs_ASD_CARseq_lfc.pdf", width=13.5, height=7)
ggarrange(gs1, gs2, gs3, gs4, gs5, gs6, 
          ncol = 3, nrow = 2)
dev.off()

# ----------------------------------------------------------------------
# check genes DE in both SCZ and ASD for Mcri
# ----------------------------------------------------------------------

table(pvals[,CARseq_Micro.SCZ] < 0.05, pvals[,CARseq_Micro.ASD] < 0.05)
c1 = chisq.test(pvals[,CARseq_Micro.SCZ] < 0.05, pvals[,CARseq_Micro.ASD] < 0.05)
c1
c1$expected

table(pvals[,CARseq_Micro.SCZ] < 0.01, pvals[,CARseq_Micro.ASD] < 0.01)
c1 = chisq.test(pvals[,CARseq_Micro.SCZ] < 0.01, pvals[,CARseq_Micro.ASD] < 0.01)
c1
c1$expected

df2 = cbind(pvals, df1)
select_cols = c("gene_id", "gene_name.SCZ", "CARseq_Micro.SCZ", "CARseq_Micro.ASD", 
                "SCZ_vs_Control.Micro", "ASD_vs_Control.Micro")

df3 = df2[CARseq_Micro.SCZ < 0.01 & CARseq_Micro.ASD < 0.01, ..select_cols]
df3[order(SCZ_vs_Control.Micro, ASD_vs_Control.Micro)]

df4 = df2[CARseq_Micro.SCZ < 0.05 & CARseq_Micro.ASD < 0.05, ..select_cols]
df4 = df4[order(SCZ_vs_Control.Micro, ASD_vs_Control.Micro)]

dim(df4)
df4[1:2,]
names(df4)[2] = "gene_name"

fwrite(df4, "../results/DEG_Micro_in_both_SCZ_ASD.csv")

# ----------------------------------------------------------------------
# check genes DE in both SCZ and ASD for Olig
# ----------------------------------------------------------------------

table(pvals[,CARseq_Oligo.SCZ] < 0.05, pvals[,CARseq_Oligo.ASD] < 0.05)
c1 = chisq.test(pvals[,CARseq_Oligo.SCZ] < 0.05, pvals[,CARseq_Oligo.ASD] < 0.05)
c1
c1$expected

table(pvals[,CARseq_Oligo.SCZ] < 0.02, pvals[,CARseq_Oligo.ASD] < 0.02)
c1 = chisq.test(pvals[,CARseq_Oligo.SCZ] < 0.02, pvals[,CARseq_Oligo.ASD] < 0.02)
c1
c1$expected

table(pvals[,CARseq_Oligo.SCZ] < 0.03, pvals[,CARseq_Oligo.ASD] < 0.03)
c1 = chisq.test(pvals[,CARseq_Oligo.SCZ] < 0.03, pvals[,CARseq_Oligo.ASD] < 0.03)
c1
c1$expected


df2 = cbind(pvals, df1)
select_cols = c("gene_id", "gene_name.SCZ", "CARseq_Oligo.SCZ", "CARseq_Oligo.ASD", 
                "SCZ_vs_Control.Oligo", "ASD_vs_Control.Oligo")

df3 = df2[CARseq_Oligo.SCZ < 0.02 & CARseq_Oligo.ASD < 0.02, ..select_cols]
df3[order(SCZ_vs_Control.Oligo, ASD_vs_Control.Oligo)]

df4 = df2[CARseq_Oligo.SCZ < 0.05 & CARseq_Oligo.ASD < 0.05, ..select_cols]
df4 = df4[order(SCZ_vs_Control.Oligo, ASD_vs_Control.Oligo)]

dim(df4)
df4[1:2,]
names(df4)[2] = "gene_name"

fwrite(df4, "../results/DEG_Oligo_in_both_SCZ_ASD.csv")

gc()

sessionInfo()
q(save="no")


