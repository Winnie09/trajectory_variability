library(here)
setwd(here())
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/pt.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/case_design.rds')
identical(names(pt),cellanno[,1])

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res <- cellPropTest(cellanno=cellanno, pseudotime=pt, design=design)
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/caseprop/res.rds')

library(ggplot2)
df <- data.frame(pt=pt,sample=cellanno[,2],class=ifelse(names(pt) %in% cellanno[cellanno[,2]%in% rownames(design)[design[,2]==1],1],1,0),stringsAsFactors = F)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/plot/caseprop/prop.pdf',width=10,height=5)
ggplot(df,aes(pt,col=sample)) + geom_density() + theme_classic() + facet_wrap(~class) + theme(legend.position = 'none')
dev.off()

