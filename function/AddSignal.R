AddSignal <- function(expr, sample, SelectGene, SelectSample, pseudotime, method, parameter=1){
  if (length(parameter) == 1) parameter[2] <- 0
  psn <- 1:length(pseudotime)
  names(psn) <- pseudotime
  if (method == 'constant'){
    expr[SelectGene,sample %in% SelectSample] <- expr[SelectGene,sample %in% SelectSample] + parameter[1]  
  } else if (method == 'linear'){
    expr[SelectGene,sample %in% SelectSample] <- expr[SelectGene,sample %in% SelectSample] + psn[colnames(expr)[sample %in% SelectSample]]/max(psn) * parameter[1] + parameter[2]
  } else if (method == 'power'){
    expr[SelectGene,sample %in% SelectSample] <- expr[SelectGene,sample %in% SelectSample] + (psn[colnames(expr)[sample %in% SelectSample]]/max(psn))^(parameter[1]/10)
  } else if (method == 'log'){
    expr[SelectGene,sample %in% SelectSample] <- expr[SelectGene,sample %in% SelectSample] + log(psn[colnames(expr)[sample %in% SelectSample]]/max(psn) * parameter[1])
  } else if (method == 'exponential'){
    expr[SelectGene,sample %in% SelectSample] <- expr[SelectGene,sample %in% SelectSample] + exp((parameter[1]/10)* psn[colnames(expr)[sample %in% SelectSample]]/max(psn))
  }
  return(expr)
}

