plotGenePopulation <- function(testobj,
                               gene,
                               type = 'time',
                               ylim = NA,
                               xlim = NA,
                               sep = NA,
                               free.scale = TRUE,
                               palette = 'Dark2'){
    ## testobj: object returned from testpt(). 
    ## gene: a character vector of gene names. It can be of length 1 or > 1.
    ## variable: a character (within the column names in design matrix) to get population pattern. If variable == NA, then return testtime fit. Other wise return testvar fit with the variable.
  library(ggplot2)
  library(RColorBrewer)
  a <- if (free.scale) 'free' else 'fixed'
    if (type == 'time'){
      fit <- sapply(gene, function(g){
        tmp <- getPopulationFit(testobj, gene = g, type = type)
      }, simplify = FALSE)
      pd <- sapply(1:length(fit), function(i){
            tmp <- reshape2::melt(fit[[i]])
            colnames(tmp) <- c('expression')
            tmp <- data.frame(tmp, 
                             gene = names(fit)[i], 
                             pseudotime = 1:length(fit[[i]]), 
                             stringsAsFactors = FALSE)
      }, simplify = FALSE)
      pd <- do.call(rbind, pd)
      if (!is.na(sep)) pd$gene <- sub(sep, '', pd$gene)
      pd$gene <- as.factor(pd$gene)
      p <- ggplot2::ggplot(data= pd, aes(x = pseudotime, y = expression, color = 'red')) + 
              geom_point() +
              theme_classic() +
              facet_wrap(~gene, nrow = round(sqrt(length(gene))), scales = a)+
              xlab('Pseudotime') +
              ylab('Expression') +
              labs(color = '')
  } else {
      fit <- getPopulationFit(testobj, gene = gene, type = type)
      pd <- sapply(1:length(fit), function(i){
            tmp <- reshape2::melt(fit[[i]])
            
            colnames(tmp) <- c('gene', 'pseudotime', 'expression')
            tmp <- data.frame(tmp, 
                             type = names(fit)[i],
                             stringsAsFactors = FALSE)
      }, simplify = FALSE)
      pd <- do.call(rbind, pd)
      if (!is.na(sep)) pd$gene <- sub(sep, '', pd$gene)
      pd$gene <- as.factor(pd$gene)
      
      p <- ggplot(data= pd, aes(x = pseudotime, y = expression, group = type, color = type)) + 
              geom_point() +
              theme_classic() +
              facet_wrap(~gene, nrow = round(sqrt(length(gene))), scales = a)+
              xlab('Pseudotime') +
              ylab('Expression') +
              labs(color = '')
      if (length(unique(pd$type)) < 8)  {
        p <-  p + scale_color_brewer(palette = palette)
      } else {
        p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(8, palette))(length(unique(pd$type))))
      }
  }
  
  if (!is.na(ylim)[1]) p <- p + ylim(ylim)
  if (!is.na(xlim)[1]) p <- p + xlim(xlim)
  print(p)
}

