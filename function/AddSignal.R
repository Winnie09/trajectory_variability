AddSignal <- function(expr, sample = NULL, SelectGene, SelectSample = NULL, pseudotime, method, parameter = 1, type = 'all'){
  if (ncol(pseudotime) > 1){
    pseudotime <- pseudotime[order(pseudotime[,2]), ]
    psn <- pseudotime[,2]
    names(psn) <- pseudotime[,1]
  } else {
    'pseudotime should be a dataframe containing two columns: cell, pseudotime.'
  }
  if (type == 'all'){
    AddSignal_AllSample(expr=expr, SelectGene = SelectGene, pseudotime =  pseudotime, method = method, parameter = parameter)
  } else {
    AddSignal_SelectSample(expr = expr, sample = sample, SelectGene = SelectGene, SelectSample, pseudotime = pseudotime, method = method, parameter = parameter)
  }
}

