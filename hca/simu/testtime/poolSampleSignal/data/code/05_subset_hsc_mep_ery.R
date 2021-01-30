setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testvar/data/plain/'
expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
rownames(expr) = sapply(rownames(expr), function(i) sub('_', '-', i))
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
saveRDS(expr, './testtime/data/data/null/hsc_mep_ery_saver.rds')

cnt <- readRDS('./hca/data/HCA/proc/matrix/count.rds')
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-', i))
a <- cnt[rownames(expr), colnames(expr)]
saveRDS(a, './testtime/data/data/null/hsc_mep_ery_count.rds')
