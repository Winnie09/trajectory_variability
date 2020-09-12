setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testtime/data/data/null/'
ct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/HCA/proc/ct/sc.rds')
dr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/null/pr.rds')
## construct pseudotime 
pseudotime <- rev(TSCANorder(exprmclust(t(dr),reduce = F,clusternum=3,clustermethod='kmeans'),orderonly = T))
psn <- 1:length(pseudotime)
names(psn) <- pseudotime
pd <- data.frame(x = dr[,1], y = dr[,2], ct = ct[rownames(dr)], pseudotime = psn[rownames(dr)], stringsAsFactors = F)
ggplot() + geom_point(data=pd, aes(x=x, y=y, col=ct))
ggplot() + geom_point(data=pd, aes(x=x, y=y, col=pseudotime))
saveRDS(data.frame(cell = names(psn), pseudotime = psn, stringsAsFactors = F), paste0(rdir, 'pseudotime.rds'))
