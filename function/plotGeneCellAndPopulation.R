plotGeneCellAndPopulation <- function(testobj,
                               gene = NA,
                               type = 'time',
                               ylim = NA,
                               xlim = NA,
                               sep = NA,
                               free.scale = TRUE,
                               palette = 'Dark2', 
                               ncol = NA,
                               subSampleNumber = NA,
                               line.size = 1, 
                               line.alpha = 1, 
                               point.size = 0.2,
                               point.alpha = 0.2,
                               axis.text.blank = F){
  ## testobj: object returned from testpt(). 
  ## gene: a character vector of gene names. It can be of length 1 or > 1.
  ## variable: a character (within the column names in design matrix) to get population pattern. If variable == NA, then return testtime fit. Other wise return testvar fit with the variable.
  library(ggplot2)
  library(RColorBrewer)
  if (is.na(ncol))  nrow = round(sqrt(length(gene))) else nrow = NA 
  a <- if (free.scale) 'free' else 'fixed'
  if ('populationFit' %in% names(testobj)) fit <- testobj$populationFit else fit = getPopulationFit(testobj, gene = gene, type = type)
  if ('expr.ori' %in% names(testobj)) expression <- testobj$expr.ori else expression <- testobj$expr
  pseudotime <- testobj$pseudotime
  cellanno <- testobj$cellanno
  design <- testobj$design
  pdlist <- list()
  if (toupper(type) == 'TIME') variable.d <- 1 else variable.d <- 2
  for (g in gene){
    pd <- data.frame(expression = expression[g, ], 
                     pseudotime = pseudotime[colnames(expression)],
                     Sample = cellanno[match(colnames(expression), cellanno[,1]),2], 
                     Variable = design[match(cellanno[,2], rownames(design)), variable.d],
                     gene = g, stringsAsFactors = F)
    pdlist[[g]] <- pd
  }
  pd <- do.call(rbind, pdlist)
  pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
  if (!is.na(sep)) {
    pd$gene <- gsub(sep, '', pd$gene)
    pd$gene <- factor(pd$gene, levels = gsub(sep, '', gene))
  } else {
    pd$gene <- factor(pd$gene, levels = gene)
  }
  
  if (toupper(type) == 'TIME'){
    if (is.na(gene)) gene <- rownames(fit)    
    ld <- sapply(gene, function(g){
      if (!is.na(sep)) g2 <- sub(sep, '', g) else g2 = g  
      tmp <- data.frame(gene = g2, expression = fit[g, ], pseudotime = testobj$pseudotime, stringsAsFactors = FALSE)
    }, simplify = F)
    ld <- do.call(rbind, ld)
    
    if (!is.na(sep)) {
      ld$gene <- factor(as.character(ld$gene), levels = sub(sep, '', gene))
    } else {
      ld$gene <- factor(as.character(ld$gene), levels = gene)
    }
    p <- ggplot2::ggplot() +
      geom_point(data= pd, aes(x = pseudotime, y = expression),color = 'black', size = point.size, alpha = point.alpha) + 
      geom_line(data = ld, aes(x = pseudotime, y = expression), color = 'red',size = line.size, alpha = line.alpha) +
      theme_classic() +
      xlab('Pseudotime') +
      ylab('Expression') +
      labs(color = '')
    if (is.na(ncol)){
      p <- p + facet_wrap(~gene, nrow = nrow, scales = a)
    } else {
      p <- p + facet_wrap(~gene, ncol = ncol, scales = a)
    }
  } else {
    if (is.na(gene)) gene <- rownames(fit[[1]])
    ld <- sapply(1:length(fit), function(i){
      
      tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
      if (!is.na(subSampleNumber)) {
        set.seed(12345)
        id <- sample(1:nrow(fit[[i]]), subSampleNumber)
        tmp <- reshape2::melt(fit[[i]][gene, id ,drop=F])
      } else {
        tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
      }
      
      colnames(tmp) <- c('gene', 'pseudotime', 'expression')
      if (!is.na(subSampleNumber)) {
        set.seed(12345)
        tmp <- tmp[sample(1:nrow(tmp), subSampleNumber), , drop=F]
      }
      tmp <- data.frame(tmp, 
                        type = names(fit)[i],
                        stringsAsFactors = FALSE)
    }, simplify = FALSE)
    ld <- do.call(rbind, ld)
    if (!is.na(sep)) {
      ld$gene <- sub(sep, '', ld$gene)
      ld$gene <- factor(as.character(ld$gene), levels = sub(sep, '', gene))
    } else{
      
      ld$gene <- factor(as.character(ld$gene), levels = gene)
      
    }
    p <- ggplot(data= ld, aes(x = pseudotime, y = expression, group = type, color = type)) + 
      geom_line(size = line.size) +
      theme_classic() +
      xlab('Pseudotime') +
      ylab('Expression') +
      labs(color = '') 
    if (is.na(ncol)){
      p <- p + facet_wrap(~gene, nrow = nrow, scales = a)
    } else {
      p <- p + facet_wrap(~gene, ncol = ncol, scales = a)
    }
    if (length(unique(ld$type)) < 8)  {
      p <-  p + scale_color_brewer(palette = palette)
    } else {
      p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(8, palette))(length(unique(ld$type))))
    }
  }
  
  if (!is.na(ylim)[1]) p <- p + ylim(ylim)
  if (!is.na(xlim)[1]) p <- p + xlim(xlim)
  if (axis.text.blank) {
    p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank())
  } else {
    p <- p + scale_x_continuous(breaks=c(min(ld$pseudotime),max(ld$pseudotime)))
  }
  print(p)
}


