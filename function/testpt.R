testpt <- function(expr, cellanno, pseudotime, design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type = 'Variable') {
  ## test.type = c('Variable', 'Time')
  if (ncol(pseudotime) > 1) {
    pseudotime <- pseudotime[order(pseudotime[,2]), ]
  } else {
    print('pseudotime should be a dataframe containing 1st column cell and 2nd column pseudotime...')
  }
  if (type == 'Variable'){
    print('Testing on variable ...')
    return(testpt_Variable(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, permuiter=permuiter, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores))
  } else if (type == 'Time'){
    print('Testing on time ... ')
    return(testpt_Time(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, permuiter=permuiter, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores))
  }
}


