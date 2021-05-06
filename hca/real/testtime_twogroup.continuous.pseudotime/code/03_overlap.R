library(here)
setwd(here())
ddir <- 'hca/real/testtime_twogroup/result/'
path = 'monocyte'

omat <- overlap <- list()
numdg2 <- numdg1 <- list()
method = 'EM_pm'
Res1 <- readRDS(paste0(ddir, method, '/', path, '/1256/testtime_res.rds'))
Res2 <- readRDS(paste0(ddir, method, '/', path, '/3478/testtime_res.rds'))
s1 = Res1$statistics
s2 = Res2$statistics
numdg1[[method]] <- sum(s1[,1] < 0.05)
numdg2[[method]] <- sum(s2[,1] < 0.05)
int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,1], -s1[,3]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,1], -s2[,3]), ]
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
mean(o)
omat[[method]] <- o
overlap[[method]] <- mean(o)

method = 'tradeSeq'
Res1 <- readRDS(paste0(ddir, method, '/', path, '/1256/testtime_res.rds'))
Res2 <- readRDS(paste0(ddir, method, '/', path, '/3478/testtime_res.rds'))
tmp <- Res1[['startVsEndTest']]
s1 <- tmp[complete.cases(tmp), ]
tmp <- Res2[['startVsEndTest']]
s2 <- tmp[complete.cases(tmp), ]
numdg1[[paste0(method, ':startVsEndTest')]] <- sum(s1[,3] < 0.05)
numdg2[[paste0(method, ':startVsEndTest')]] <- sum(s2[,3] < 0.05)

int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,1]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,1]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
omat[[paste0(method, ':startVsEndTest')]] <- o
overlap[[paste0(method, ':startVsEndTest')]] <- mean(o)

tmp <- Res1[['associationTest']]
s1 <- tmp[complete.cases(tmp), ]
tmp <- Res2[['associationTest']]
s2 <- tmp[complete.cases(tmp), ]
numdg1[[paste0(method, ':associationTest')]] <- sum(s1[,3] < 0.05)
numdg2[[paste0(method, ':associationTest')]] <- sum(s2[,3] < 0.05)

int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,1]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,1]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
omat[[paste0(method, ':associationTest')]] <- o
overlap[[paste0(method, ':associationTest')]] <- mean(o)

##
method = 'tscan'
s1 <- readRDS(paste0(ddir, method, '/', path, '/1256/testtime_res.rds'))
s2 <- readRDS(paste0(ddir, method, '/', path, '/3478/testtime_res.rds'))
numdg1[[method]] <- sum(s1[,3] < 0.05)
numdg2[[method]] <- sum(s2[,3] < 0.05)

int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,2]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,2]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
omat[[method]] <- o
overlap[[method]] <- mean(o)

## monocle2
method = 'monocle2'
s1 <- readRDS(paste0(ddir, method, '/', path, '/1256/testtime_res.rds'))
s2 <- readRDS(paste0(ddir, method, '/', path, '/3478/testtime_res.rds'))
numdg1[[method]] <- sum(s1[,3] < 0.05)
numdg2[[method]] <- sum(s2[,3] < 0.05)

int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,2]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,2]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply((1:100)*100, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
omat[[method]] <- o
overlap[[method]] <- mean(o)
saveRDS(overlap, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_twogroup/result/perf/overlap.rds')
saveRDS(numdg1, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_twogroup/result/perf/numdg1_1256.rds')
saveRDS(numdg2, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_twogroup/result/perf/numdg_3478.rds')

omat <- do.call(cbind, omat)
saveRDS(omat, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_twogroup/result/perf/overlap_all.rds')

# omat <-readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/real/testtime_twogroup/result/perf/overlap_all.rds')
library(reshape2)
pd <- melt(omat)
library(ggplot2)
library(RColorBrewer)
ggplot() + geom_freqpoly(data = pd, aes(value,color = Var2), fill = NA, binwidth = 0.01) + scale_color_brewer(palette = 'Set1') + theme_classic()

pdf('/Users/wenpinhou/Dropbox/trajectory_variability/hca/real/testtime_twogroup/plot/overlap_proportion_all.pdf', width = 6, height = 3.5)
ggplot(data = pd, aes(x=Var1, y = value, color = Var2)) + geom_line(size = 0.1) + 
  geom_point(size = 0.1) +
  scale_color_brewer(palette = 'Set1') + theme_classic() + 
  xlab('top n*10 genes') + ylab('overlap proportion')
dev.off()


