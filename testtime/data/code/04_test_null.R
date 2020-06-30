# [whou10@jhu.edu@bigmem0032 code]$ cat 04_test_null.R 
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
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'

### reorder cells for each patient
pt <- readRDS(paste0(rdir,'null/pseudotime.rds'))
pt <- data.frame(pt, sample = sub(':.*','', pt$cell), stringsAsFactors = F)
ptlist <- lapply(unique(pt[,3]), function(p){
  tmp <- pt[as.character(pt[,3])==p, ]
  set.seed(12345)
  tmp[, 'pseudotime'] <- sample(tmp[, 'pseudotime'])
  tmp
})
pseudotime <- do.call(rbind, ptlist)
design = cbind(1,design)
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)

# > str(pseudotime)
# 'data.frame': 13269 obs. of  3 variables:
#  $ cell      : chr  "BM4:29:male_164021" "BM4:29:male_280662" "BM4:29:male_76213" "BM4:29:male_56540" ...
#  $ pseudotime: int  250 4362 4437 7898 4975 3665 11342 464 6586 10414 ...
#  $ sample    : chr  "BM4" "BM4" "BM4" "BM4" ...

if (method == 'EM_SelectKnots'){
  testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time')
  saveRDS(testres, paste0(rdir, method, '/testres.rds'))
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
  saveRDS(Final, paste0(rdir, method,'/testres.rds'))  
}

if (method == 'tscan'){
  psn <- as.numeric(pseudotime[,2])
  names(psn) <- pseudotime[,1]
  expr <- expr[, pseudotime[,2]]
  res <- TSCAN_time(expr=expr,pseudotime=psn)
  saveRDS(res, paste0(rdir, method,'/testres.rds'))  
}

if (method == 'monocle2'){
  psn <- as.numeric(pseudotime[,2])
  names(psn) <- pseudotime[,1]
  expr <- expr[, pseudotime[,2]]
  res <- monocle2_time(expr=expr,pseudotime=psn)
  saveRDS(res, paste0(rdir, method,'/testres.rds'))
}

if (method == 'monocle3'){
  psn <- as.numeric(pseudotime[,2])
  names(psn) <- pseudotime[,1]
  expr <- expr[, pseudotime[,2]]
  pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/null/hsc_mep_ery_integrated_pca.rds')
  res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
  saveRDS(res, paste0(rdir, method,'/testres.rds'))
}


