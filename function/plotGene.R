plotGene <- function(testobj, gene, variable = NULL, variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = FALSE, line.alpha = 1, line.size = 1, point.alpha=1, point.size=0.5, continuous = TRUE, sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = F, cellProp = FALSE, x.lab = 'Pseudotime', y.lab = 'Expression', use.palette = FALSE, legend.position = 'right', axis.text.range.only = FALSE, ylim = NA, xlim = NA, axis.size = 8, title.size = 8, axis.text.size = 8){
  ## testobj: the output of function testpt() which is a list containing fdr, etc..
  ## variable: character, the variable (covariate) to color the samples, should be null or one of the column names of design matrix. Default is NULL, meaning each sample is colored differently. Otherwise, samples are colored by the variable (covariate) values.
  ## variable.text: a character vector. The text for the legend of the plot, corresponding to each variable values.
  ## continuous: if TRUE, samples are colored using viridis continuous colors. If FALSE, RColorBrewer "Dark2" discrete palette.
  ## expression: a character ('demean',or 'original') to define the expression values shown on the plots. if "demean"(default), show demeaned expression. if 'original", show original gene expression.
  ## ncol: only functional when plotting multiple genes. Used to define the number of columns.
  ## cellProp: if FALSE (default), plot gene expression. If TRUE, it is cell proportion.
  ## use.palette. If TRUE, forced to use color palette defined by the argument "palette". If FALSE (defult), only use "palette" when continuous = FALSE. 
  library(splines)
  library(ggplot2)
  library(gridExtra)
  library(viridis)   
  library(RColorBrewer) ##
  pseudotime <- testobj[['pseudotime']]
  cellanno <- testobj[['cellanno']]
  colnames(cellanno) <- c('Cell', 'Sample')
  if ('expr.ori' %in% names(testobj)) expression <- testobj[['expr.ori']] else expression <- testobj[['expr']]
  if (cellProp){
    ptw <- cut(pseudotime,seq(min(pseudotime),max(pseudotime),length.out = 100),include.lowest = T)
    ptdat <- table(ptw,cellanno[match(names(pseudotime),cellanno[,1]),2])
    ptdat <- t(t(ptdat)/colSums(ptdat)) ## divided by rowsum (rowsum = 1). interval * samples. 
    ptdat <- as.data.frame(ptdat)
    colnames(ptdat) <- c('pt','s','prop')
    ptdat[,1] <- match(ptdat[,1],levels(ptw))
    
    ptdat$cell <- paste0('cell',1:nrow(ptdat))
    ptexpr <- t(ptdat[,c('prop','prop'),drop=F])
    colnames(ptexpr) <- ptdat$cell
    
    ptpt <- ptdat$pt
    names(ptpt) <- ptdat$cell
    expr=ptexpr
    cellanno=data.frame(cell=ptdat$cell,sample=ptdat$s)
    pseudotime=ptpt
  }
  
  predict.values <- predict_fitting(testobj, gene= gene, test.type = testobj$test.type)
  pseudotime = pseudotime[colnames(expression)]
  cellanno <- cellanno[match(colnames(expression), cellanno[,1]), ]
  # predict.values <- predict.values[, colnames(expression),drop=F]
  knotnum <- testobj$knotnum
  knotnum[knotnum==0] <- 1  ## in case the fitting of line would cause bugs
  design <- testobj[['design']]
  
  
  cellanno <- data.frame(Cell = as.character(cellanno[,1]), 
                         Sample = as.character(cellanno[,2]), stringsAsFactors = FALSE)
  variable.d <- if(is.null(variable)) 1 else variable
  if (!is.null(variable.text) & variable.d != 1) {
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
      # tmpcell <- cellanno[cellanno[,2]==p,1]
      tmpcellid <- which(cellanno[,2]==p)
      if (toupper(testobj$test.type) == 'TIME'){  ##### add
        tmpdf <- data.frame(expr=predict.values[gene,tmpcellid], 
                            Sample=p, 
                            Variable = design[rownames(design) == p, variable.d], ##
                            pseudotime=pseudotime[tmpcellid])
      } else {
        tmpdf <- data.frame(
          expr = predict.values[[which(as.numeric(sub('.*_', '', names(predict.values))) == design[p, variable.d])]][gene, tmpcellid],
          Sample = p,
          Variable = design[rownames(design) == p, variable.d], ##
          pseudotime=pseudotime[tmpcellid])
      }
    })
    ld <- do.call(rbind, linedlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    ld <- ld[order(ld$pseudotime), ] ## add 20200812
    if (is.null(variable)){
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Sample), alpha=point.alpha, size=point.size, stroke = 0) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size) 
        
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)   
      }
    } else {
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Variable), alpha=point.alpha, size=point.size, stroke = 0) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample),alpha=line.alpha, size=line.size) 
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample),alpha=line.alpha, size=line.size)   
      }   
    }
    p <- p + 
      theme_classic() +
      # ggtitle(paste0(sub(':.*','',gene),',adj.pvalue=', formatC(testobj$fdr[gene], format = "e", digits = 2))) +
      xlab(x.lab) + ylab(y.lab) + 
      labs(color = variable) +
      theme(legend.spacing.y = unit(0.01, 'cm'), legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.1, "cm"), legend.position = legend.position) +
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))
    if (!is.na(sep)){
      p <- p + ggtitle(sub(sep,'',gene)) 
    } else {
      p <- p + ggtitle(gene)
    }
    if (continuous & !use.palette){
      p <- p + scale_color_viridis(discrete = TRUE, direction = -1)   
    } else {
      if (length(unique(ld[, 'Variable'])) > 8){
        p <- p + scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, palette)))(length(unique(ld[, 'Variable']))))  
      } else {
        p <- p + scale_color_manual(values = rev(brewer.pal(8, palette)[1:length(unique(ld[, 'Variable']))]))
      } 
    }
    
    if (!is.na(xlim)) p <- p + xlim(xlim)
    if (!is.na(ylim)) p <- p + ylim(ylim)
    
     if (axis.text.blank) {
      p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank())
    } else if (axis.text.range.only) {
      p <- p + scale_x_continuous(breaks=c(min(pd$pseudotime),max(pd$pseudotime)))
    }
    
    p <- p + theme(legend.spacing.y = unit(0.01, 'cm'), legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.1, "cm"), axis.title = element_text(size = axis.size), axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.size), axis.text.y = element_text(size = axis.text.size), plot.title = element_text(size = title.size), legend.position = legend.position) +
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) 
    if (facet.sample){
      if (is.null(ncol)){
        print(p + facet_wrap(~Sample, scales=a))  
      } else {
        print(p + facet_wrap(~Sample, scales=a, ncol = ncol))
      }
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
                       g = g, stringsAsFactors = F)
      pdlist[[g]] <- pd
      linedlist <- lapply(unique(testobj$cellanno[,2]), function(p){
        tmpcellid <- which(cellanno[,2]==p) 
        if (toupper(testobj$test.type) == 'TIME'){  ##### add
          tmpdf <- data.frame(expr=predict.values[g,tmpcellid], 
                              pseudotime=pseudotime[tmpcellid],
                              Sample=p, 
                              Variable = design[rownames(design) == p, variable.d], 
                              g = g, stringsAsFactors = F)
        } else {
          tmpdf <- data.frame(
            expr = predict.values[[which(as.numeric(sub('.*_', '', names(predict.values))) == design[p, variable.d])]][g, tmpcellid],
            pseudotime=pseudotime[tmpcellid],
            Sample = p,
            Variable = design[rownames(design) == p, variable.d], 
            g = g, stringsAsFactors = F)
        }
      })
      ld = do.call(rbind, linedlist)
      ldlist[[g]] <- ld
    }
    pd <- do.call(rbind, pdlist)
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    ld <- do.call(rbind, ldlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    ld <- ld[order(ld$pseudotime), ] 
    if (!is.na(sep)) {
      pd$g <- gsub(sep, '', pd$g)
      ld$g <- gsub(sep, '', ld$g)
      pd$g <- factor(pd$g, levels = gsub(sep, '', gene))
      ld$g <- factor(ld$g, levels = gsub(sep, '', gene))
    } else {
      pd$g <- factor(pd$g, levels = gene)
      ld$g <- factor(ld$g, levels = gene)
    }
    if (is.null(variable)){
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Sample), alpha=point.alpha, size=point.size, stroke = 0) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)
        
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)   
      }
      
    } else {
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Variable), alpha=point.alpha, size=point.size, stroke = 0) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample), alpha=line.alpha, size=line.size) 
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample), alpha=line.alpha, size=line.size) 
      }
    }
    
    p <- p + 
      theme_classic() + 
      xlab(x.lab) + ylab(y.lab) + 
      labs(color = variable) +
      facet_wrap(~g, scales = a, ncol = ncol) +
      theme(legend.spacing.y = unit(0.01, 'cm'), legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.1, "cm"), legend.position = legend.position, axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank())+
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))
    
    if (continuous & !use.palette){
      p <- p + scale_color_viridis(discrete = TRUE, direction = -1) 
    } else {
      if (length(unique(ld[, 'Variable'])) > 8){
        p <- p + scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, palette)))(length(unique(ld[, 'Variable']))))  
      } else {
        p <- p + scale_color_manual(values = rev(brewer.pal(8, palette)[1:length(unique(ld[, 'Variable']))]))
      } 
    }
    if (axis.text.blank) {
      p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank())
    } else if (axis.text.range.only) {
      p <- p + scale_x_continuous(breaks=c(min(pd$pseudotime),max(pd$pseudotime)))
    }
    if (!is.na(xlim)) p <- p + xlim(xlim)
    if (!is.na(ylim)) p <- p + ylim(ylim)
    print(p)
  }
  
}


