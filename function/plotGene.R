 plotGene <- function(testptObj, gene, variable = NULL, variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = FALSE, line.alpha = 1, line.size = 1, point.alpha=1, point.size=0.5, continuous = TRUE, original.expr = FALSE){
  ## testptObj: the output of function testpt() which is a list containing fdr, etc..
  ## variable: character, the variable (covariate) to color the samples, should be null or one of the column names of design matrix. Default is NULL, meaning each sample is colored differently. Otherwise, samples are colored by the variable (covariate) values.
  ## variable.text: a character vector. The text for the legend of the plot, corresponding to each variable values.
  ## continuous: if TRUE, samples are colored using viridis continuous colors. If FALSE, RColorBrewer "Dark2" discrete palette.
  ## expression: a character ('demean',or 'original') to define the expression values shown on the plots. if "demean"(default), show demeaned expression. if 'original", show original gene expression.
  library(splines)
  library(ggplot2)
  library(gridExtra)
  library(viridis)   
  library(RColorBrewer) ##
  pseudotime <- testptObj[['pseudotime']]
  cellanno <- testptObj[['cellanno']]
  colnames(cellanno) <- c('Cell', 'Sample')
  if (original.expr){
    expression <- testptObj[['expr.ori']]
    predict.values <- predict_fitting(expr = expression, knotnum = testptObj$knotnum, design = testptObj$design, cellanno = cellanno, pseudotime = pseudotime[colnames(expression)])
  } else {
    expression <- testptObj[['expr.demean']]
    predict.values <- testptObj$predict.values
  }
  pseudotime = pseudotime[colnames(expression)]
  cellanno <- cellanno[match(colnames(expression), cellanno[,1]), ]
  predict.values <- predict.values[, colnames(expression)]
  knotnum <- testptObj$knotnum
  knotnum[knotnum==0] <- 1  ## in case the fitting of line would cause bugs
  design <- testptObj[['design']]
  
  
  cellanno <- data.frame(Cell = as.character(cellanno[,1]), 
                        Sample = as.character(cellanno[,2]), stringsAsFactors = FALSE)
  variable.d <- if(is.null(variable)) 1 else variable
  if (!is.null(variable.text) & exists('variable.d')) {
    design[,variable.d] <- ifelse(design[, variable.d] == 0, variable.text[1], variable.text[2])
  }
  a <- if (free.scale) 'free' else 'fixed'
  if (length(gene) == 1){
    print('plotting one gene ...')
    pd <- data.frame(expr = expression[gene, ], 
                     Sample = cellanno[,2], 
                     Variable = design[match(cellanno[,2], rownames(design)), variable.d],  ##
                     pseudotime = pseudotime[colnames(expression)])
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    linedlist <- lapply(unique(cellanno[,2]), function(p){
      tmpcell <- cellanno[cellanno[,2]==p,1]
      tmpdf <- data.frame(expr=predict.values[gene,tmpcell], 
                          Sample=p, 
                          Variable = design[rownames(design) == p, variable.d], ##
                          pseudotime=pseudotime[tmpcell])
    })
    ld <- do.call(rbind, linedlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    ld <- ld[order(ld$pseudotime), ] ## add 20200812
    if (is.null(variable)){
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Sample), alpha=point.alpha, size=point.size) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size) 
          
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)   
      }
    } else {
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Variable), alpha=point.alpha, size=point.size) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample),alpha=line.alpha, size=line.size) 
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample),alpha=line.alpha, size=line.size)   
      }   
    }
    p <- p + 
      theme_classic() +
      # ggtitle(paste0(sub(':.*','',gene),',adj.pvalue=', formatC(testptObj$fdr[gene], format = "e", digits = 2))) +
      ggtitle(sub(':.*','',gene)) +
      xlab('Pseudotime') + ylab('Expression') + 
      labs(color = variable)
    if (continuous){
      p <- p + scale_color_viridis(discrete = TRUE, direction = -1) 
    } else {
      p <- p + scale_color_brewer(palette = 'Dark2') 
    }
      if (facet.sample){
          print(p + facet_wrap(~Sample, scales=a))
      } else {
          print(p)
      }
  } else {
    print('plotting multiple genes ...')
    pdlist <- ldlist <- list()
    for (g in gene){
      pd <- data.frame(expr = expression[g, ], 
                       pseudotime = pseudotime[colnames(expression)],
                       Sample = cellanno[match(colnames(expression), cellanno[,1]),2], 
                       Variable = design[match(cellanno[,2], rownames(design)), variable.d],
                       g = g)
      pdlist[[g]] <- pd
      linedlist <- lapply(unique(testptObj$cellanno[,2]), function(p){
        tmpcell <- cellanno[cellanno[,2]==p, 1]
        tmpdf <- data.frame(expr=predict.values[g,tmpcell], 
                            pseudotime=pseudotime[tmpcell], 
                            Sample=p, 
                            Variable = design[rownames(design) == p, variable.d],
                            g=g)
      })
      ld = do.call(rbind, linedlist)
      ldlist[[g]] <- ld
    }
    pd <- do.call(rbind, pdlist)
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    ld <- do.call(rbind, ldlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    ld <- ld[order(ld$pseudotime), ] 
    if (is.null(variable)){
      if (plot.point){
         p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Sample), alpha=point.alpha, size=point.size) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)
          
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)   
      }
        
    } else {
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Variable), alpha=point.alpha, size=point.size) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample), alpha=line.alpha, size=line.size) 
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample), alpha=line.alpha, size=line.size) 
      }
    }
    p <- p + 
      theme_classic() + 
      xlab('Pseudotime') + ylab('Expression') + 
      labs(color = variable) +
      facet_wrap(~g, scales = a) 
    if (continuous){
      p <- p + scale_color_viridis(discrete = TRUE, direction = -1) 
    } else {
      p <- p + scale_color_brewer(palette = 'Dark2') 
    }
    print(p)
  }
}

