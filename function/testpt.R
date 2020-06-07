testpt <- function(expr, cellanno, pseudotime, design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type = 'Variable') {
  ## test.type = c('Variable', 'Time')
  if (type == 'Variable'){
    print('Testing on variable ...')
    return(testpt_Variable(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, permuiter=permuiter, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores))
  } else if (type == 'Time'){
    print('Testing on time ... ')
    return(testpt_Time(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, permuiter=permuiter, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores))
  }
}

