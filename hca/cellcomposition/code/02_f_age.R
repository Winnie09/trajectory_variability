setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
order = readRDS('./hca/result/ery/order.rds')
source('./function/01_function.R')
order = data.frame(order, Patient = gsub('_.*','', order$Cell))
ap = as.character(unique(order$Patient))
age = as.numeric(sapply(ap, function(i) strsplit(i,':')[[1]][2]))
g1 = ap[order(age)[1:4]]
g2 = ap[order(age)[5:8]]
f <- f_statistics_from_pseudotime(order, g1, g2)
saveRDS(f,'./hca/cellcomposition/result/f_statistics_from_pseudotime_age.rds')
f_pm = f_statistics_from_pseudotime_permute(order, g1, g2, num.permute=1e4)
saveRDS(f_pm,'./hca/cellcomposition/result/f_statistics_from_pseudotime_age_permute.rds')
