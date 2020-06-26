SensFdr <- function(TruePositive, Statistics){
  ## Statistics: a dataframe or matrix, should contain a column of fdr , for example c('FStat','P.Value','adj.P.Val')
  ## find the column of fdr
  fdrchar <- intersect(colnames(Statistics), c('adj.P.Val','adj.pvalue','fdr','FDR','Fdr','adj.p', 'adj.P','adj.Pval'))
  fdrcol <- which(colnames(Statistics) == fdrchar)
  ## if not ordered by significance, then rank by significance
  Statistics = Statistics[complete.cases(Statistics), , drop = FALSE]
  if (sum((diff(Statistics[, fdrcol])<0)+0) > 0){
    Statistics <- Statistics[order(Statistics[,fdrcol]), , drop = FALSE]
  }
  Order <- rownames(Statistics)
  ## calculate sensitivity, realfdr, reported fdr
  perf <- t(sapply(seq(1,length(Order)), function(i){
    num <- sum(Order[seq(1,i)] %in% TruePositive)
    c(num/length(TruePositive), (i - num)/i, Statistics[i, fdrcol])
  }))
  ## reorder real fdr so that it is monotonically increasing
  if (nrow(perf) > 1){
    for (i in (nrow(perf)):2){
      if (perf[i-1,2] > perf[i,2]) perf[i-1,2] <- perf[i,2]
    }
  }
  colnames(perf) <- c('Sensitivity','Real_FDR','Reported_FDR')
  rbind(c(0,0,0),perf)
}

