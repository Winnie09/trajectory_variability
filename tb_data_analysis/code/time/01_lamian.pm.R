library(rhdf5)
pid <- 2
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')
## female = 1, male = 0

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/01_function.R')
res <- testpt(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/exprpc',pid,'.h5'), cellanno, pt, design=design,test.type='Time', ncores = 18)
saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/time/pc2/lamian_pm_test.rds'))


