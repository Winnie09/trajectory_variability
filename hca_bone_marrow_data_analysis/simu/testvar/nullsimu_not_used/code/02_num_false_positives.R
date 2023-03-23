library(here)
setwd(here())
ddir <- 'hca/simu/testvar/nullsimu/result/'
rdir <- 'hca/simu/testvar/nullsimu/result/perf/'
dir.create(rdir)
res.tscan <- readRDS(paste0(ddir, 'tscan/res.rds'))
res.em <- readRDS(paste0(ddir, 'EM_pm/res.rds'))
res.limma <- readRDS(paste0(ddir, 'limma/res.rds'))
res.tradeseq <- readRDS(paste0(ddir, 'tradeSeq/res.rds'))

num <- c(sum(res.em[[1]][,1] < 0.05), 
         sum(res.limma[,5] < 0.05), 
         sum(res.tradeseq[['earlyDETest']][,3] < 0.05, na.rm = T),
         sum(res.tradeseq[['patternTest']][,3] < 0.05, na.rm = T),
         sum(res.tradeseq[['diffEndTest']][,3] < 0.05, na.rm = T))
names(num) <- c('Lamian', 'Limma', 'earlyDETest', 'patternTest', 'diffEndTest')
saveRDS(num, paste0(rdir, 'number_of_false_positives.rds'))


pdf('hca/simu/testvar/nullsimu/plot/pvalue_distribution.pdf', width = 7.5, height = 5)
par(mfrow=c(2,3))
hist(res.em[[1]][,2], main = 'Lamian', col = 'grey', xlab = 'p values')
hist(res.limma[,4], main = 'limma', col = 'grey', xlab = 'p values')
hist(res.tradeseq[['earlyDETest']][,2], main = 'tradeSeq_earlyDETest', col = 'grey', xlab = 'p values')
hist(res.tradeseq[['patternTest']][,2], main = 'tradeSeq_patternTest', col = 'grey', xlab = 'p values')
hist(res.tradeseq[['diffEndTest']][,2], main = 'tradeSeq_diffEndTest', col = 'grey', xlab = 'p values')
dev.off()


