setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
order = readRDS('./hca/result/ery/order.rds')

source('./function/01_function.R')
order = data.frame(order, Patient = gsub('_.*','', order$Cell))
ap = as.character(unique(order$Patient))
g1 = ap[grepl('female', ap)]
g2 = ap[grepl(':male', ap)]
f <- f_statistics_from_pseudotime(order, g1, g2)
saveRDS(f,'./hca/cellcomposition/result/f_statistics_from_pseudotime_gender.rds')

f_pm = f_statistics_from_pseudotime_permute(order, g1, g2, num.permute=1e4)
saveRDS(f_pm,'./hca/cellcomposition/result/f_statistics_from_pseudotime_gender_permute.rds')
