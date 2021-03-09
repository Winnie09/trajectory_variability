getDEGType <- function(testobj, cutoff = 0.05, test.type = 'variable', test.method = 'permutation'){
  res <- testobj$statistics
  if (test.type == 'variable'){
    if (test.method == 'permutation'){
      diffType <- sapply(rownames(res), function(g){
        if (res[g, 1] < cutoff){
          if (res[g, 4] < cutoff & res[g, 7] > cutoff){
            'meanSig'
          } else if (res[g, 4] > cutoff & res[g, 7] < cutoff){
            'trendSig'
          } else if (res[g, 4] < cutoff & res[g, 7] < cutoff){
            'bothSig'
          } else {
            'other'
          }
        } else {
          'nonDEG'
        }
      })
    } else if (test.method == 'chisq'){
      diffType <- sapply(rownames(res), function(g){
        if (res[g, 1] < cutoff){
          if (res[g, 5] < cutoff & res[g, 3] > cutoff){
            'meanSig'
          } else if (res[g, 5] > cutoff & res[g, 3] < cutoff){
            'trendSig'
          } else if (res[g, 5] < cutoff & res[g, 3] < cutoff){
            'bothSig'
          } else {
            'other'
          }
        } else {
          'nonDEG'
        }
      })
    }
    
  } else if (test.type == 'time'){
    print('work on it now!')
  }
  
  return(diffType)
}


