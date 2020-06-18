method <- as.character(commandArgs(trailingOnly = T)[[1]])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testtime/data/data/'
dir.create(paste0(rdir,method), recursive = T, showWarnings = F)
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
allp = sub(':.*','', colnames(expr))
sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
design <- cbind(rep(1,8))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'

### do not reorder cells for each patient
pt <- readRDS(paste0(rdir,'null/pseudotime.rds'))
pt <- data.frame(pt, sample = sub(':.*','', pt$cell), stringsAsFactors = F)
pseudotime = pt
design = cbind(1,design)
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)


if (method == 'EM_SelectKnots'){
  testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time')
  saveRDS(testres, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/result/NoReorderNoSignalNull/EM_SelectKnots_testres.rds')
}
  
if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
 
  counts <- round(exp(expr + 1))
  pdt <- data.frame(curve1 = pseudotime[,2], curve2 = pseudotime[,2])
  rownames(pdt) <- pseudotime[,1]
  pdt = pdt[colnames(counts), ]
  
  v <- (cellanno$sample %in% paste0('BM',c(1,2,5,6)) + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(counts)
  
  set.seed(12345)
  sce <- fitGAM(counts = counts, pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE)
  saveRDS(sce, paste0(rdir, method,'/sce.rds'))
    
  Final <- list()
  for (TestType in (c('diffEndTest', 'patternTest', 'earlyDETest'))){
    print(TestType)
    if (grepl('diffEndTest', TestType)){
      Res <- diffEndTest(sce)  
    } else if (grepl('patternTest', TestType)){
      Res <- patternTest(sce)  
    } else if (grepl('earlyDETest', TestType)){
      Res <- earlyDETest(sce, knots = c(1,2), global = TRUE, pairwise = TRUE)
    }
    res <- data.frame(waldStat = Res[,'waldStat'], P.Value = Res[,'pvalue'] ,adj.P.Val = p.adjust(Res$pvalue, method='fdr'))
    row.names(res) <- row.names(Res)
    res <- res[order(res[,3], -res[,1]), ]
    final <- list()
    final[['res']] <- res
    Final[[TestType]] <- final
  }
  saveRDS(Final, paste0(rdir, method,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/result/NoReorderNoSignalNull/tradeSeq_testres.rds'))  
}

###### check what statistics is most suitable for selecting the gene sets to add signals: the hightest signals should contribute to smallest fdr.  This analysis shows that, for genes with across each of the sample sd larger would have fdr smaller. And this relationship is stronger in sd than in cv or mean. Therefore, we need to use sd as the statistics where larger sd means larger signal. 
fdr <- testres$fdr
sd <- apply(expr, 1, sd)
mean <- rowMeans(expr)
cv <- sd/cv
cor(fdr, sd)
# [1] -0.2259065
cor(fdr, mean)
# [1] -0.1786241
cor(fdr, cv)
# [1] -0.1397419

sds <- sapply(unique(allp), function(p){
  tmp <- expr[, allp == p]
  sd <- apply(tmp, 1, sd)   
  rank(sd)
})
means <- sapply(unique(allp), function(p){
  tmp <- expr[, allp == p]
  mean <- rowMeans(tmp)
  rank(mean)
})
cvs <- sapply(unique(allp), function(p){
  tmp <- expr[, allp == p]
  cv <- apply(tmp, 1, sd)/rowMeans(tmp)
  rank(cv)
})
cor(rowMeans(sds), fdr)
# [1] -0.4480778
cor(rowMeans(means), fdr)
# [1] -0.2982547
cor(rowMeans(cvs), fdr)
# [1] -0.3191218


sds <- sapply(unique(allp), function(p){
  tmp <- expr[, allp == p]
  sd <- apply(tmp, 1, sd)   
})
means <- sapply(unique(allp), function(p){
  tmp <- expr[, allp == p]
  mean <- rowMeans(tmp)
})
cvs <- sapply(unique(allp), function(p){
  tmp <- expr[, allp == p]
  cv <- apply(tmp, 1, sd)/rowMeans(tmp)
})

> cor(rowMeans(sds), fdr)
# [1] -0.2445643
> cor(rowMeans(means), fdr)
# [1] -0.1768864
> cor(rowMeans(cvs), fdr)
# [1] NA
