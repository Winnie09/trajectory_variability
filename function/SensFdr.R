SensFdr <- function(Order, TruePositive,statistics){
  ## statistics columns are: c('FStat','P.Value','adj.P.Val')
  perf <- t(sapply(seq(1,length(Order)), function(i){
    num <- sum(Order[seq(1,i)] %in% selgene)
    c(num/length(TruePositive), (i - num)/i, statistics$adj.P.Val[i])
  }))
  if (nrow(perf) > 1){
    for (i in seq(2, nrow(perf))){
      if (perf[i-1,2] > perf[i,2]) perf[i-1,2] <- perf[i,2]
    }
  }
  colnames(perf) <- c('Sensitivity','Real_FDR','Reported_FDR')
  rbind(c(0,0,0),perf)
}
