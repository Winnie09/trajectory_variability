getXDEType <- function(testobj, cutoff = 0.05){
  res <- testobj$statistics
  if (toupper(testobj$test.type) == 'VARIABLE'){
    diffType <- sapply(rownames(res), function(g){
      if (res[g, grep('^fdr.*overall$', colnames(res))] < cutoff){
        if (res[g, grep('^fdr.*meanDiff$', colnames(res))] < cutoff & res[g, grep('^fdr.*trendDiff$', colnames(res))] > cutoff){
          'meanSig'
        } else if (res[g, grep('^fdr.*meanDiff$', colnames(res))] > cutoff & res[g, grep('^fdr.*trendDiff$', colnames(res))] < cutoff){
          'trendSig'
        } else if (res[g, grep('^fdr.*meanDiff$', colnames(res))] < cutoff & res[g, grep('^fdr.*trendDiff$', colnames(res))] < cutoff){
          'bothSig'
        } else {
          'other'
        }
      } else {
        'nonXDE'
      }
    })
  } else if (toupper(testobj$test.type) == 'TIME'){
    print('ConstantTimeTest does not lead to XDEType!')
  }
  return(diffType)
}



