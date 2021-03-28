library(here)
setwd(here('hca/simu/testvar/nullsimu/data/'))
source(here('function','01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')

### load saver, count matrix, and pseudotime
saverlog <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','hsc_mep_ery_saver.rds'))
cnt <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','hsc_mep_ery_count.rds'))
pt <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','pseudotime.rds'))
saverlog <- saverlog[, pt[,1]]
cnt <- cnt[, pt[,1]]

### permute the sample-cell relationsihp 
identical(colnames(saverlog), colnames(cnt))
sn <- sub(':.*', '', colnames(cnt))
cn <- sapply(1:ncol(cnt), function(i) sub(paste0(sn[i], ':'), '', colnames(cnt)[i]))
# plot(saverlog[grep('XIST', rownames(cnt)),], col = as.factor(sn %in% c(paste0('BM', c(1,2,5,6)))))
# plot(saverlog[grep('XIST', rownames(cnt)),], col = as.factor(sn))
set.seed(12345)
pt[,1] <- colnames(saverlog) <- colnames(cnt) <- paste0(sample(sn),':', cn)
pt.v <- pt[,2]
names(pt.v) <- pt[,1]
saveRDS(pt, 'pseudotime_pm_df.rds')
saveRDS(pt.v, 'pseudotime_pm.rds')
saveRDS(saverlog, 'saverlog_pm.rds')
saveRDS(cnt, 'cnt_pm.rds')

cellanno_pm <- data.frame(cell = colnames(cnt), sample = sub(':.*', '', colnames(cnt)), stringsAsFactors = F)
saveRDS(cellanno_pm, 'cellanno_pm.rds')

design <- matrix(c(rep(1,8), 1,1,0,0,1,1,0,0), nrow = 8)
dimnames(design) <- list(paste0('BM', 1:8), c('intercept', 'group'))
saveRDS(design, 'design.rds')
