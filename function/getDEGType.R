getDEGType <- function(testobj, cutoff = 0.05){
  res <- testobj$statistics
  diffType <- sapply(rownames(res), function(g){
    if (res[g, 7] < cutoff){
      if (res[g, 1] < cutoff & res[g, 4] > cutoff){
        'meanSig'
      } else if (res[g, 1] > cutoff & res[g, 4] < cutoff){
        'trendSig'
      } else if (res[g, 1] < cutoff & res[g, 4] < cutoff){
        'bothSig'
      } else {
        'other'
      }
    } else {
      'nonDEG'
    }
  })
  return(diffType)
}
