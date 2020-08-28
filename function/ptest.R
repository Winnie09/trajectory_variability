ptest <- function(expr, cellanno, pseudotime, design=NULL, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type='Time', test.position = 'all', fit.resolution = 1000){
  ## expr: gene by cell matrix. values are log-transformed library-size-adjusted expression values. zero-expression genes should have been filtered out.
  ## cellanno: a dataframe whose first column is cells' names. second column is cells' samples (for example patients).
  ## pseudotime: a numeric vector of pseudotime for the cells, and the names of this vector entires are cells' names.
  ## design: a numeric dataframe or matrix whose first column is all 1, second column is covariate values. column names are intercept, covariate1, etc.. 
  print('testing mean difference ...')
  res1 <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design=design, permuiter=permuiter, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores, type=type, test.pattern = 'intercept', test.position = test.position, fit.resolution = fit.resolution, return.all.data = TRUE)
  print('testing trend difference ...')
  res2 <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design=design, permuiter=permuiter, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores, type=type, test.pattern = 'slope', test.position = test.position, fit.resolution = fit.resolution, return.all.data = FALSE)
  g <- rownames(expr)
  res <- data.frame(interceptdiff.fdr = res1$fdr[g], 
                   interceptdiff.lfc = res1$foldchange[g], 
                   interceptdiff.diff = res1$interceptdiff[g],
                   interceptdiff.pvalue = res1$pvalue[g],
                   trenddiff.fdr = res2$fdr[g],
                   trenddiff.lfc = res2$foldchange[g],
                   trenddiff.diff = res2$interceptdiff[g],
                   trenddiff.pvalue = res2$pvalue[g],
                   stringsAsFactors = FALSE)
  rownames(res) <- g
  return(list(res = res, 
              interceptdiff.parameter=res1$parameter, 
              interceptdiff.orill=res1$orill, 
              interceptdiff.perll = res1$perll, 
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


 
