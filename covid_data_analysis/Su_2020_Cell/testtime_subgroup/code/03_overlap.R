ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testtime_subgroup/result/'
rdir <- paste0(ddir, 'perf/')
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testtime_subgroup/plot/'
dir.create(pdir)
dir.create(rdir)
compareid <- (1:100)*10
omat <- overlap <- list()
numdg2 <- numdg1 <- list()
method = 'EM_pm'
Res1 <- readRDS(paste0(ddir, method, '/testtime_res_subgroup1.rds'))
Res2 <- readRDS(paste0(ddir, method, '/testtime_res_subgroup2.rds'))

s1 = Res1$statistics
s2 = Res2$statistics
numdg1[[method]] <- sum(s1[,1] < 0.05)
numdg2[[method]] <- sum(s2[,1] < 0.05)
int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,1], -s1[,3]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,1], -s2[,3]), ]
o <- sapply(compareid, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
mean(o)
omat[[method]] <- o
overlap[[method]] <- mean(o)

method = 'tradeSeq'
Res1 <- readRDS(paste0(ddir, method, '/testtime_res_subgroup1.rds'))
Res2 <- readRDS(paste0(ddir, method, '/testtime_res_subgroup2.rds'))
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
o <- sapply(compareid, function(i){
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
o <- sapply(compareid, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
omat[[paste0(method, ':associationTest')]] <- o
overlap[[paste0(method, ':associationTest')]] <- mean(o)

##
method = 'tscan'
s1 <- readRDS(paste0(ddir, method, '/testtime_res_subgroup1.rds'))
s2 <- readRDS(paste0(ddir, method, '/testtime_res_subgroup2.rds'))
numdg1[[method]] <- sum(s1[,3] < 0.05)
numdg2[[method]] <- sum(s2[,3] < 0.05)

int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,2]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,2]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply(compareid, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
omat[[method]] <- o
overlap[[method]] <- mean(o)

## monocle2
method = 'monocle2'
s1 <- readRDS(paste0(ddir, method, '/testtime_res_subgroup1.rds'))
s2 <- readRDS(paste0(ddir, method, '/testtime_res_subgroup2.rds'))
numdg1[[method]] <- sum(s1[,3] < 0.05)
numdg2[[method]] <- sum(s2[,3] < 0.05)

int <- intersect(rownames(s1), rownames(s2))
s1 = s1[int, ]
s1 = s1[order(s1[,3], -s1[,2]), ]
s2 = s2[int, ]
s2 = s2[order(s2[,3], -s2[,2]), ]
print(max(which(s1[,3]<0.05)))
o <- sapply(compareid, function(i){
  length(intersect(rownames(s1)[1:i], rownames(s2)[1:i]))/i
})
omat[[method]] <- o
overlap[[method]] <- mean(o)

saveRDS(overlap, paste0(rdir, 'overlap.rds'))
saveRDS(numdg1, paste0(rdir,'numdg1.rds'))
saveRDS(numdg2, paste0(rdir, 'numdg2.rds'))
saveRDS(omat, paste0(rdir,'overlap_all.rds'))

df <- data.frame(overlap = unlist(overlap), numdg1 = unlist(numdg1), numdg2 = unlist(numdg2))
write.csv(df, paste0(pdir, 'overlap_numdg1_numdg2.csv'))

library(reshape2)
library(ggplot2)
library(RColorBrewer)
pd <- melt(do.call(cbind,omat))
pdf(paste0(pdir, 'overlap_proportion_all.pdf'), width = 6, height = 3.5)
ggplot(data = pd, aes(x=Var1, y = value, color = Var2)) + geom_line(size = 0.1) + 
  geom_point(size = 0.1) +
  scale_color_brewer(palette = 'Set1') + theme_classic() + 
  xlab('top n*10 genes') + ylab('overlap proportion')
dev.off()
