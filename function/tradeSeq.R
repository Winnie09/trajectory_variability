tradeSeq_testvar <- function(count, cellanno, design, pseudotime){
  ## count: gene by cell count matrix
  ## cellanno: 1st column cell name (chr), 2nc column sample name (chr)
  ## design: 1st column the group division (num), colname are the variable name (e.g. 'group'). rownames are sample names.
  ## pseudotime: 1st column cell name (chr), 2nd column pseudotime (num).
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))

  pdt <- data.frame(curve1 = pseudotime[,2], curve2 = pseudotime[,2])
  rownames(pdt) <- pseudotime[,1]
  pdt = pdt[colnames(count), ]
  v <- (cellanno[,2] %in% rownames(design)[which(design[,1]==1)] + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(counts)
  ### run test
  set.seed(12345)
  sce <- fitGAM(counts = count, pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE)
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
    Final[[TestType]] <- res[order(res[,3], -res[,1]), ]
  }
  return(Final)
}




