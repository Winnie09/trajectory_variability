allf <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/permu/linear')
res <- t(sapply(allf, function(f){
  r <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/permu/linear/', f))
  r[['sensfdr']]
}))
colnames(res)[1] <- 'parameter'
res <- res[order(res[,1]),]

par(mfrow=c(1,2))
plot(res[,2]~res[,1], pch=20, xlab='Slope of Linear Signal', ylab='(area under Real_FDR ~ Reported_FDR)-0.25*0.25/2', main='Use Permutation P-value')

plot(res[,3]~res[,1], pch=20, xlab='Slope of Linear Signal', ylab='Area under sensitivity-fdr curve', main='Use Permutation P-value')


