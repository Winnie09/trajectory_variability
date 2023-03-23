pid <- 2
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res <- cellPropTest(cellanno=cellanno, pseudotime=pt, design=design)
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sexprop/pc2.rds')


library(ggplot2)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sexprop/pc2.pdf')
ggplot(data.frame(cellanno,pt=pt[cellanno[,1]],sex=ifelse(design[match(cellanno[,2],rownames(design)),2]==1,'female','male')),aes(pt,group=sample,col=sex)) + geom_density() + theme_classic() + scale_color_manual(values=c('orange','royalblue'))
dev.off()

