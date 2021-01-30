meandiff <- function(expr, cellanno, design){
  ## expr: librarysize-adjusted log-normalized gene by cell expression matrix. Row names are genes, col names are cell names. 
  ## design: a numeric dataframe or matrix, each row is a sample, column 1 is all 1, column 2 is covariate that want to be tested on. rownames are sample names. 
  ## cellanno: a character dataframe. column 1 is the cell names. column 2 is the sample names.
  library(limma)
  if (sum(rownames(design) %in% unique(cellanno[,2])) != nrow(design)) print('Design matrix does not include all the samples in rows!')
  design = as.matrix(design)
  cellanno <- data.frame(cell = as.character(cellanno[,1]), sample = as.character(cellanno[,2]), stringsAsFactors = FALSE )
  agg <- vapply(rownames(design), function(i)
    rowMeans(expr[, cellanno[,2] == i, drop = FALSE], na.rm = TRUE), 
    numeric(nrow(expr)))
  agg <- agg[, rownames(design)]
  res <- topTable(eBayes(lmFit(agg, design)),n=nrow(agg),coef=2)
  return(res)
}
  
