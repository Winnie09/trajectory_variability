setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca')
order =readRDS('./result/ery/order.rds')

ct = readRDS('./data/HCA/proc/ct/sc.rds')

order = data.frame(order, Patient = gsub('_.*','', order$Cell), Celltype = ct[match(order$Cell, names(ct))])
v = sapply(order$Patient, function(i) ifelse(grepl('female', i)=='TRUE','Female','Male'))
order = data.frame(order, Gender = v)

library(ggplot2)
library(gridExtra)
pdf('./plot/time/time_patient_celltype_sampleplot_histogram.pdf',width=7,height=3)
ggplot() + geom_histogram(data=order,aes(x = Pseudotime, fill=Celltype), alpha=.3) + facet_wrap(~Patient, scales='free',nrow=2) + theme_classic()
dev.off()

pdf('./plot/time/time_patient_celltype_sampleplot_density.pdf',width=7,height=3)
ggplot() + geom_density(data=order,aes(x = Pseudotime, fill=Celltype), alpha=.3) + facet_wrap(~Patient, scales='free',nrow=2) + theme_classic()
dev.off()

pdf('./plot/time/time_patient_density.pdf',width=7,height=3)
ggplot() + geom_density(data=order,aes(x = Pseudotime, fill=Gender), alpha=.3) + facet_wrap(~Patient, nrow=2) + theme_classic()
dev.off()

pdf('./plot/time/time_patient_celltype_density_freescale.pdf',width=12,height=12)
ggplot() + geom_density(data=order,aes(x = Pseudotime, fill=Celltype), alpha=.3) + facet_grid(Celltype~Patient, scales='free') + theme_classic()
dev.off()

pdf('./plot/time/time_patient_celltype_density.pdf',width=12,height=12)
ggplot() + geom_density(data=order,aes(x = Pseudotime, fill=Celltype), alpha=.3) + facet_grid(Celltype~Patient) + theme_classic()
dev.off()

pdf('./plot/time/time_patient_celltype.pdf',width=12,height=12)
ggplot() + geom_histogram(data=order,aes(x = Pseudotime, fill=Celltype)) + facet_grid(Celltype~Patient) + theme_classic()
dev.off()


