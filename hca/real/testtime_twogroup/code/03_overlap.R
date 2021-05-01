library(here)
setwd(here())
ddir <- 'hca/real/testtime_twogroup/result/'
path = 'monocyte'

overlap <- list()
# method = 'EM_pm'
# Res1 <- readRDS(paste0(ddir, method, '/', path, '/1256/testtime_res.rds'))
# Res2 <- readRDS(paste0(ddir, method, '/', path, '/3478/testtime_res.rds'))
# s1 = Res1$statistics
# s2 = Res2$statistics
# int <- intersect(rownames(s1), rownames(s2))
# s1 = s1[int, ]
# s1 = s1[order(s1[,1], -s1[,3]), ]
# s2 = s2[int, ]
# s2 = s2[order(s2[,1], -s2[,3]), ]
# print(max(which(s1[,3]<0.05)))
# o <- sapply((1:100)*100, function(i){
#   length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
# })
# mean(o)


method = 'tradeSeq'
Res1 <- readRDS(paste0(ddir, method, '/', path, '/1256/testtime_res.rds'))
Res2 <- readRDS(paste0(ddir, method, '/', path, '/3478/testtime_res.rds'))
tmp <- Res1[['startVsEndTest']]
s1 <- tmp[complete.cases(tmp), ]
tmp <- Res2[['startVsEndTest']]
s2 <- tmp[complete.cases(tmp), ]
int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,1]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,1]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
overlap[[paste0(method, ':startVsEndTest')]] <- mean(o)

tmp <- Res1[['associationTest']]
s1 <- tmp[complete.cases(tmp), ]
tmp <- Res2[['associationTest']]
s2 <- tmp[complete.cases(tmp), ]
int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,1]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,1]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
overlap[[paste0(method, ':associationTest')]] <- mean(o)

##
method = 'tscan'
s1 <- readRDS(paste0(ddir, method, '/', path, '/1256/testtime_res.rds'))
s2 <- readRDS(paste0(ddir, method, '/', path, '/3478/testtime_res.rds'))
int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,2]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,2]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
overlap[[method]] <- mean(o)

## monocle2
method = 'monocle2'
s1 <- readRDS(paste0(ddir, method, '/', path, '/1256/testtime_res.rds'))
s2 <- readRDS(paste0(ddir, method, '/', path, '/3478/testtime_res.rds'))
int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,2]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,2]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
overlap[[method]] <- mean(o)

saveRDS(overlap, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_twogroup/result/perf/overlap.rds')
