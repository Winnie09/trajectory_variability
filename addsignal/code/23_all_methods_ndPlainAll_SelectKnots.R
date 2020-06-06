rm(list=ls())
geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
method <- as.character(commandArgs(trailingOnly = T)[[4]])
## be reminded that the pca is done on cv > 0.5 (for all samples) genes. You should redo it !!!
# geneProp <- 0.05
# addSignalType <- 'constant'
# addSignalPara <-  2.1
# method <- 'EM_SelectKnots'

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('./function/01_function.R')
datadir <- './addsignal/data/ndPlainAll/1256/'
rdir <- './addsignal/result/ndPlainAll/'
if (!file.exists(paste0(rdir, method, '/', addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))){
  suppressMessages(library(parallel))
  suppressMessages(library(splines))
  suppressMessages(library(limma))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  dir.create(paste0(rdir, method, '/', addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
  res <- readRDS(paste0(datadir,addSignalType,'/', geneProp,'_', addSignalPara,'.rds')) #####
  expr <- res$expr
  selgene <- res$selgene
  pseudotime <- res$pseudotime
  design <- res$design
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  
  ## two group along pseudotime
  if (method == 'limma'){
    df <- t(sapply(0:20, function(i){
      print(i)
      dr <- readRDS('./addsignal/result/dr.rds')
      r <- exprdiff(expr=expr,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=FALSE)
      saveRDS(r, paste0(rdir, method, addSignalType,'/',geneProp,'_',i,'.rds'))
      r <- r[order(r[,3]),]
      sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
      res <- list()
      res[['res']] <- r
      res[['sensfdr']] <- c(i,AreaUnderSensFdr(sensfdr))
      print(c(i,AreaUnderSensFdr(sensfdr)))
      return(c(i,AreaUnderSensFdr(sensfdr)))
    }))
    colnames(df)[1] <- c('Parameter')
    saveRDS(df, paste0(rdir, method,addSignalType,'/',geneProp,'_Para_FdrDiff_Area_Group1256.rds'))
  } 
  
  if (method == 'permu'){
    if (!file.exists(paste0(rdir, method, addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))){
      dr <- readRDS('./addsignal/result/dr.rds')
      res <- exprdiff(expr=expr,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=TRUE, permutime=10000)
      r <- res[['res']]
      r <- r[order(r[,3]),]
      sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
      v <- c(addSignalPara,AreaUnderSensFdr(sensfdr))
      names(v)[1] = 'Parameter'
      res[['sensfdr']] <- v
      saveRDS(res, paste0(rdir, method, addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))
    }
  } 
  
  if (method == 'permu_IR'){ ## ignore repetitive
    if (!file.exists(paste0(rdir, method, addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))){
      dr <- readRDS('./addsignal/result/dr.rds')
      res <- exprdiff(expr=expr,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=TRUE, permutime=10000, IgnoreRepetitive = TRUE)
      r <- res[['res']]
      r <- r[order(r[,3]),]
      sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
      v <- c(addSignalPara,AreaUnderSensFdr(sensfdr))
      names(v)[1] = 'Parameter'
      res[['sensfdr']] <- v
      saveRDS(res, paste0(rdir, method, addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))
    }
  } 
  
  if (method == 'permu_IRGroup'){ ## ignore repetitive
    if (!file.exists(paste0(rdir, method, addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))){
      dr <- readRDS('./addsignal/result/dr.rds')
      res <- exprdiff(expr=expr,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=TRUE, permutime=10000, IgnoreRepetitive = TRUE)
      r <- res[['res']]
      r <- r[order(r[,3]),]
      sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
      v <- c(addSignalPara,AreaUnderSensFdr(sensfdr))
      names(v)[1] = 'Parameter'
      res[['sensfdr']] <- v
      saveRDS(res, paste0(rdir, method, addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))
    }
  } 
  
  if (grepl('tradeSeq', method)){
      suppressMessages(library(tradeSeq))
    # dir.create(paste0('./addsignal/result/', method, '/',addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
    if (!file.exists(paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))){
      counts <- round(exp(expr + 1))
      
      psn <- seq(1, length(pseudotime))
      names(psn) <- pseudotime
      pdt <- data.frame(curve1 = psn, curve2 = psn)
      rownames(pdt) <- names(psn)
      pdt = pdt[colnames(counts), ]
      
      v <- (cellanno$sample %in% paste0('BM',c(1,2,5,6)) + 0)
      v <- ifelse(v==1, 0.99, 0.01)
      cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
      rownames(cellWeights) <- colnames(counts)
      
      set.seed(12345)
      sce <- fitGAM(counts = counts, pseudotime = pdt, cellWeights = cellWeights,
                    nknots = 6, verbose = FALSE,parallel=TRUE)
      saveRDS(sce, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))
    } else {
      sce <- readRDS(paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))
    }
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
      sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
      final <- list()
      final[['res']] <- res
      final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
      Final[[TestType]] <- final
    }
    saveRDS(Final, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))  
  }
  
  if (method == 'EM'){
    print(method)
    design = cbind(1,design)
    
    library(gtools)
    per <- permutations(nrow(design),nrow(design))
    id <- which(!duplicated(apply(per,1,function(i) paste0(as.vector(design[i,]),collapse = '_'))))
    per <- per[id,]
    oriid <- which(apply(per,1,function(i) paste0(as.vector(design[i,]),collapse = '_'))==paste0(as.vector(design),collapse = '_'))
    library(parallel)
    # perll <- mclapply(1:2,function(i) {
    perll <- mclapply(1:nrow(per),function(i) {
      perdesign <- design[per[i,],,drop=F]
      row.names(perdesign) <- row.names(design)
      diffpt(expr=expr,design=perdesign,pseudotime=pseudotime,num.knot = 3,cellanno = cellanno, verbose = T)$logL
    }, mc.cores = 4)
    perll <- do.call(cbind,perll)
    pval <- sapply(1:nrow(perll), function(i) pnorm(perll[i,oriid],mean(perll[i,-oriid]),sd(perll[i,-oriid]),lower.tail = F))
    fdr <- p.adjust(pval,method='fdr')
    names(fdr) <- row.names(perll)
    saveRDS(fdr, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_fdr.rds'))  
    res <- data.frame(P.Value = pval, adj.P.Val = fdr, stringsAsFactors = F)
    rownames(res) <- names(fdr)
    res <- res[order(res[,2]),]
    sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
    final <- list()
    final[['res']] <- res
    final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
    final[['perll']] <- perll
    saveRDS(final, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))  
  }
  
  if (method == 'EM_SelectKnots'){
    print(method)
    design = cbind(1,design)
    psn <- 1:length(pseudotime)
    names(psn) <- pseudotime
    pseudotime <- psn
    testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100)
    saveRDS(testres$fdr, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_fdr.rds'))  
    saveRDS(testres, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_testres.rds'))  
    res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
    rownames(res) <- names(testres$fdr)
    res <- res[order(res[,1]),,drop=F]
    sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
    final <- list()
    final[['res']] <- res
    final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
    final[['perll']] <- testres$perll
    final[['knotnum']] <- testres$knotnum
    saveRDS(final, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))  
  }
  rm(list=ls())
}


