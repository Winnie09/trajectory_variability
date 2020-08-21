ptest <- function(expr, cellanno, pseudotime, design=NULL, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type='Time', test.pattern = c('intercept', 'slope', 'overall'), test.position = 'all', fit.resolution = 1000){
  print('testing mean difference ...')
  res1 <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design=design, permuiter=permuiter, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores, type=type, test.pattern = 'intercept', test.position = test.position, fit.resolution = fit.resolution, return.all.data = TRUE)
  print('testing trend difference ...')
  res2 <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design=design, permuiter=permuiter, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores, type=type, test.pattern = 'slope', test.position = test.position, fit.resolution = fit.resolution, return.all.data = FALSE)
  g <- rownames(expr)
  res <- data.frame(meandiff.fdr = res1$fdr[g], 
                   meandiff.lfc = res1$foldchange[g], 
                   meandiff.diff = res1$meandiff[g],
                   trenddiff.fdr = res2$fdr[g],
                   trenddiff.lfc = res2$foldchange[g],
                   trenddiff.diff = res2$meandiff[g],
                   stringsAsFactors = FALSE)
  rownames(res) <- g
  return(list(res = res, 
              meandiff.parameter=res1$parameter, 
              meandiff.orill=res1$orill, 
              meandiff.perll = res1$perll, 
              trenddiff.parameter=res2$parameter, 
              trenddiff.orill=res2$orill, 
              trenddiff.perll = res2$perll, 
              knotnum = res1$knotnum,  
              predict.values = res1$pred[,colnames(expr)], 
              pseudotime = pseudotime[colnames(expr)], 
              design = design, 
              cellanno = cellanno, 
              expression = expr))
}


 
