setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('./function/01_function.R')
order = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/result/HCA/order.rds')

mat = readRDS('./hca/data/HCA/proc/matrix/saver.rds')
mat = mat[,order$Cell]
ap = gsub('_.*','', colnames(mat))
p = ap[1]
paralist <- lapply(unique(ap), function(p){
  p = gsub(':','_',p)
  para = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/result/HCA//spline_coef_', p, '.rds'))
})
names(paralist) = unique(ap)
group = read.table('./simu/data/meta/group.txt', as.is = T)
g1 = group[,1]
g2 = group[,2]
paralist = paralist[c(g1,g2)]

f <- sapply(seq(1, nrow(paralist[[1]])), function(i){
  tmp <- t(sapply(names(paralist), function(p){
    paralist[[p]][i, ]  
  }))
  f <- ANOVA_model_f_stat(tmp, c(rep(1,length(g1)), rep(0, length(g2))))
})
names(f) = rownames(paralist[[1]])
saveRDS(f, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu2/result/f_statistics.rds')

