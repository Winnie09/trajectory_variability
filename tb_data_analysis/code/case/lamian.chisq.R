library(rhdf5)

pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/pt.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/case_design.rds')
path <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/expr.h5'

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/01_function.R')
res <- testpt('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/expr.h5', cellanno, pt, design=design,test.method = 'chisq',test.type='Variable',EMitercutoff=10)
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/case/lamian_chisq.rds')
