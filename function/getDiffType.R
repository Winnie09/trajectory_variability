getDiffType <- function(testobj, cutoff = 0.05){
  res <- testobj$statistics
  diffType <- sapply(rownames(res), function(g){
    if (res[g, 7] < cutoff){
      if (res[g, 1] < cutoff & res[g, 4] > cutoff){
        'meanOnly'
      } else if (res[g, 1] > cutoff & res[g, 4] < cutoff){
        'trendOnly'
      } else if (res[g, 1] < cutoff & res[g, 4] < cutoff){
        'both'
      } else {
        'unknown'
      }
    } else {
      'nonDiff'
    }
  })
  return(diffType)
}
