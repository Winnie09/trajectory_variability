library(here)
setwd(here())
rdir <- './hca/real/build3traj/manual/result/erythroid/'
pdir <- './hca/real/build3traj/manual/plot/erythroid/'
dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
expr <- readRDS('./hca/data/proc/matrix/saver.rds') ## [1:22401, 1:32819] 
str(expr)
ct <- readRDS('./hca/data/proc/ct/sc.rds')
table(ct)
id <- which(ct %in% c('HSC','MPP','CMP', 'MEP', 'Ery')) ##
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9415, 1:8886] 
pr <- PCA(genebycellmat = expr, save.pca = T, plot.statistics=F, plot.dir = pdir, result.dir = rdir, PC_for_columns = TRUE, findVariableGenes = TRUE, maxVariableGenes = 3000, numPC = NULL, smoothFittingMethod = 'loess')
UMAP(samplebyfeature_mat = pr, save.umap = T, result.dir = rdir)

# pr <- readRDS('./hca/real/build3traj/manual/result/lymph/pr.rds')
# prvar <- apply(pr, 2, var)
# plot(prvar, prh=19)
set.seed(12345)
clu <- kmeans(pr, 3)$cluster
dr <- pr
## construct pseudotime 
library(TSCAN)
library(ggplot2)
pseudotime <- TSCANorder(exprmclust(t(dr), reduce = T, cluster = clu), MSTorder = c(3,2,1), orderonly = TRUE)
# qplot(ct[pseudotime],1:length(pseudotime))
psn <- 1:length(pseudotime)
names(psn) <- pseudotime
pd <- data.frame(x = dr[,1], y = dr[,2], ct = ct[rownames(dr)], pseudotime = psn[rownames(dr)],cluster = clu[rownames(dr)],  stringsAsFactors = F)

p1 <- ggplot() + geom_point(data=pd, aes(x=x, y=y, col=ct), size = 0.5, alpha = 0.5)+theme_classic() + xlab('PC1') + ylab('PC2')
p2 <- ggplot() + geom_point(data=pd, aes(x=x, y=y, col=pseudotime), size = 0.5)+theme_classic()+ xlab('PC1') + ylab('PC2')
p3 <- ggplot() + geom_point(data=pd, aes(x=x, y=y, col=cluster), size = 0.5)+theme_classic()+ xlab('PC1') + ylab('PC2')
pdf(paste0(pdir, 'pca.pdf'), width = 7, height = 2.5)
gridExtra::grid.arrange(p1,p2,p3, nrow = 1)
dev.off()
saveRDS(pseudotime, paste0(rdir, 'pseudotime.rds'))

