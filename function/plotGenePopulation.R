plotGenePopulation <- function(testobj,
                               gene,
                               variable = NULL){
    ## testobj: object returned from testpt(). testvar is ready. testtime is not ready (should adjust the fit data structure a bit).
    ## gene: a character vector of gene names. It can be of length 1 or > 1.
    ## variable: a character (within the column names in design matrix) to get population pattern.
    if (is.null(variable)){
    print('variable = NULL, but testtime is not ready.')
  } else {
      fit <- sapply(gene, function(g){
        tmp <- get_population_fit(Res, variable = variable, gene = g)
      }, simplify = FALSE)
      
      pd <- sapply(1:length(fit), function(i){
            tmp <- reshape2::melt(fit[[i]])
            colnames(tmp) <- c('cell', 'type', 'expression')
            tmp <- data.frame(tmp, 
                             gene = names(fit)[i], 
                             pseudotime = 1:nrow(fit[[i]]), 
                             stringsAsFactors = FALSE)
      }, simplify = FALSE)
      pd <- do.call(rbind, pd)
      pd$gene <- as.factor(pd$gene)
      
      p <- ggplot2::ggplot(data= pd, aes(x = pseudotime, y = expression, group = type, color = type)) + 
              geom_point() +
              theme_classic() +
              facet_wrap(~gene, nrow = round(sqrt(length(gene))))+
              xlab('Pseudotime') +
              ylab('Expression') +
              labs(color = variable)
      if (length(unique(pd$type)) < 8)  {
        p <-  p + scale_color_brewer(palette = 'Dark2')
      } else {
        p <- p + RColorBrewer::colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(pd$type)))
      }
  }
  print(p)
}

