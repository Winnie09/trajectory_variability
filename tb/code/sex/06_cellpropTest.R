library(rhdf5)
pid <- 2
ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex_addBatch/'
pt <- readRDS(paste0(ddir,'ptpc',pid,'.rds'))
cellanno <- readRDS(paste0(ddir, 'cellanno.rds'))
design <- readRDS(paste0(ddir, 'design.rds'))
design = design[, c(1, 39, 2:38)]

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
af = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function', pattern = 'multi.R$')
for (f in af){
  source(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/', f))
}
res <- cellPropTest(cellanno=cellanno, pseudotime=pt, design=design)
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sexprop/pc2.rds')

pd = data.frame(cellanno,pt=pt[cellanno[,1]],sex=ifelse(design[match(cellanno[,2],rownames(design)),2]==1,'female','male'))
saveRDS(pd,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sexprop/pc2_pd.rds')



