
library("readxl")
library(ggplot2)
theme_set(theme_minimal())

# ----------------------------------------------------------------------
# just to visualize the results of micro-specific DE genes SCZ vs. ASD 
# ----------------------------------------------------------------------

dat = read_excel("../results/GO_seq_SCZ_vs_ASD_micro.xlsx")
dim(dat)
dat[1:2,]
names(dat)[5]="qval"

dat = as.data.frame(dat)
dat$category = gsub()
dim(dat)
dat[1:2,]

w1 = which(dat$category=="SRP DEPENDENT COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE")
dat$category[w1] = "SRP DEPENDENT COTRANSLATIONAL PROTEIN TARGETING."
dat  = dat[which(dat$qval < 0.05),]
dat$category = as.factor(dat$category)
dat$category = reorder(dat$category, nrow(dat):1)

p = ggplot(data=dat, aes(x=category, y=-log10(qval), fill=category)) +
  geom_bar(stat="identity", color="black") + 
  coord_flip() + theme(legend.position='none') +
  scale_fill_manual(values=rep("#56B4E9",nrow(dat)))

pdf("../figures/GO_seq_SCZ_vs_ASD_micro.pdf", width=6, height=3)
p
dev.off()

gc()

sessionInfo()
q(save="no")


