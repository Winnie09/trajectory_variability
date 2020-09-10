clusterType <- as.numeric(commandArgs(trailingOnly = T)[[1]])
pctGene <- as.numeric(commandArgs(trailingOnly = T)[[2]])

clusterType = 9
pctGene = 1
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu')
source('./function/01_function.R')
datadir <- './testtime/data/data/'
rdir <- './testtime/result/'
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/'
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)

## one group along pseudotime
# ------------------------------------------------
# see if  demean the data will change the results
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/'
pseudotime <- readRDS(paste0(ddir, 'null/pseudotime.rds'))
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/EM_SelectKnots/'
expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
design = matrix(rep(1,8), nrow=8)
dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
expr <- log2(expr + 1)
print(method)
design = cbind(1,design)
pt <- pseudotime[,2]
names(pt) <- pseudotime[,1]
pseudotime = pt

for (p in unique(cellanno[,2])){
  expr[, cellanno[,2] == p] = expr[, cellanno[,2] == p] - rowMeans(expr[, cellanno[,2] == p])
}

testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time')

# saveRDS(testres, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_testres.rds'))  
res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
rownames(res) <- names(testres$fdr)
res <- res[order(res[,1]),,drop=F]
final <- list()
final[['res']] <- res
final[['perll']] <- testres$perll
final[['knotnum']] <- testres$knotnum
# saveRDS(final, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
saveRDS(testres, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/EM/clusterType9_1.rds')

oldres <- readRDS(paste0(rdir, 'EM_SelectKnots/clusterType9_1_testres.rds'))
plot(oldres$fdr ~ testres$fdr, pch = 20)
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/selgene/selgene.rds')

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/plot/compare_demean_and_not_demean.pdf', width = 5, height = 3.5)
plot(oldres$fdr ~ testres$fdr, pch = 20, col = ifelse(names(oldres$fdr) %in% selgene, 'red', 'black'), xlab = 'fdr (centered data)', ylab = 'fdr (un-centered data)')
dev.off()

g <- names(which(testres$fdr < 0.05 & oldres$fdr > 0.5))
plotGene(testres, g[g %in% selgene])
plotGene(testres, g[!g %in% selgene])


g <- names(which(testres$fdr > 0.5 & oldres$fdr < 0.05))
expr.ori <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/null/hsc_mep_ery_saver.rds')

length(names(which(testres$fdr < 0.2 & oldres$fdr < 0.2)))
length(names(which(testres$fdr < 0.05 & oldres$fdr < 0.05)))
AreaUnderSensFdr(SensFdr(selgene, data.frame(fdr = testres$fdr, foldchange = testres$foldchange)))
