mat <- readRDS('/home-4/zji4@jhu.edu/scratch/raisin/HCA/common/proc/matrix/rawcount.rds')
n <- readRDS('/home-4/zji4@jhu.edu/scratch/raisin/HCA/common/proc/qc/expressedgenenumber.rds')
c <- readRDS('/home-4/zji4@jhu.edu/scratch/raisin/HCA/common/proc/qc/totalreadcount.rds')
m <- readRDS('/home-4/zji4@jhu.edu/scratch/raisin/HCA/common/proc/qc/mitoproportion.rds')
mat <- mat[,n>500 & m < 0.1]
mat <- mat[rowSums(mat) > 0,]
saveRDS(mat,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcarev/data/matrix/count.rds')
