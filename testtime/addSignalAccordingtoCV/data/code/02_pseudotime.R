setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testtime/data/data/null/'
pc <- readRDS('./hca/data/HCA/proc/integrate/pca/pca.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
pc <- pc[names(ct[ct %in% c('HSC','MEP','Ery')]),]
saveRDS(pc, paste0(rdir, 'hsc_mep_ery_integrated_pca.rds'))

# pcvar <- apply(pc, 2, var)
# plot(pcvar, pch=19)

dr <- pc[,1:4]
## construct pseudotime 
library(TSCAN)
pseudotime <- rev(TSCANorder(exprmclust(t(dr),reduce = F,clusternum=3,clustermethod='kmeans'),orderonly = T))
psn <- 1:length(pseudotime)
names(psn) <- pseudotime
pd <- data.frame(x = dr[,1], y = dr[,2], ct = ct[rownames(dr)], pseudotime = psn[rownames(dr)], stringsAsFactors = F)
library(ggplot2)
ggplot() + geom_point(data=pd, aes(x=x, y=y, col=ct), size = 0.1, alpha = 0.1)
ggplot() + geom_point(data=pd, aes(x=x, y=y, col=pseudotime))
saveRDS(data.frame(cell = names(psn), pseudotime = psn, stringsAsFactors = F), paste0(rdir, 'pseudotime.rds'))
