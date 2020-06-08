AddSignal_AllSample <- function(expr, SelectGene, pseudotime, method, parameter=1){
  if (length(parameter) == 1) parameter[2] <- 0
  psn <- 1:length(pseudotime)
  names(psn) <- pseudotime
  expr <- expr[, names(psn)]
  if (method == 'constant'){
    expr[SelectGene,] <- t(t(expr[SelectGene,]) + parameter[1]  )
  } else if (method == 'linear'){
    expr[SelectGene,] <- t(t(expr[SelectGene,]) + psn/max(psn) * parameter[1] + parameter[2]) ###
  } else if (method == 'power'){
    expr[SelectGene,] <- t(t(expr[SelectGene,]) + (psn/max(psn))^(parameter[1]/10))
  } else if (method == 'log'){
    expr[SelectGene,] <- t(t(expr[SelectGene,]) + log(psn/max(psn) * parameter[1]))
  } else if (method == 'exponential'){
    expr[SelectGene,] <- t(t(expr[SelectGene,]) + exp((parameter[1]/10)* psn/max(psn)))
  }
  return(expr)
}

