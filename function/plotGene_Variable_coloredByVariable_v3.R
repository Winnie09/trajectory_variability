plotGene_Variable_coloredByVariable_v3 <- function(testptObj, gene, expression, alpha=1, size=0.5, plot.points = FALSE, scale.free = FALSE, by.sample = FALSE, variable = 'Var2'){
  ## expression: gene by cell Matrix
  ## testptObj: the output of function testpt() which is a list containing fdr, etc..
  library(ggplot2)
  library(gridExtra)
  library(viridis)   
  ## extract the values we need for plotting
  order <- data.frame(cell = names(pseudotime), pseudotime = pseudotime, stringsAsFactors = FALSE)
  cellano = testptObj$cellanno
  colnames(cellano) <- c('cell', 'Sample')
  pseudotime = testptObj$pseudotime
  design = testptObj$design
  pred <- testptObj$predict.values
  knotnum <- testptObj$knotnum
  knotnum[knotnum==0] <- 1  ## in case the fitting of line would cause bugs
  ## make all cells the same order according to pseudotime
  expression <- expression[, names(pseudotime), drop=F]
  cellanno = cellanno[match(names(pseudotime), cellanno[,1]), ]
  pred <- pred[, names(pseudotime)]
  
  if (scale.free) {
    a  = 'free'
  } else {
    a  = 'fixed'
  }
  if (length(gene) == 1){
    print('plotting one gene ...')
    pd <- data.frame(Expr = expression[gene, ], cell = cellano[,1], Sample = cellano[,2], Variable = design[match(cellano[,2], rownames(design)), 1], stringsAsFactors = F)
    pd = cbind(pd, pseudotime = order[match(pd$cell, order$cell),'pseudotime'])
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    linedlist <- lapply(unique(pd$Sample), function(p){
    #   tMat = expression[gene, which(cellano[,2]==p), drop=F]
    #   trainX = order[match(colnames(tMat), order$cell),'pseudotime']      ### use time 
    #   pred <- get_spline_fit(tMat, trainX=trainX, fit.min=min(order$pseudotime), fit.max=max(order$pseudotime), num.base = knotnum[gene])
    #   tmpdf <- data.frame(Expr=pred[1,], pseudotime=seq(min(order$pseudotime),max(order$pseudotime),length.out = 1000), Sample=p, Variable = design[rownames(design) == p, 1])
        tmpdf = data.frame(Expr = pred[gene, gsub(':.*','', colnames(pred)) == p], pseudotime=seq(min(order$pseudotime),max(order$pseudotime), length.out = ncol(pred)/nrow(design)), Sample = p, Variable = design[rownames(design) == p, 1])
    })
    ld = do.call(rbind, linedlist)

    ##
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    if (plot.points){
      if (by.sample){
        ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=Expr, color=Variable, group = Sample), alpha=alpha, size=size)  +
          geom_line(data=ld, aes(x=pseudotime, y=Expr, color=Variable, group = Sample),alpha=1, size=.5) +
          theme_classic() +
          scale_color_viridis(discrete = TRUE, direction = -1) +
          ggtitle(paste0(sub(':.*','',gene),',adjP=', formatC(testptObj$fdr[gene], format = "e", digits = 2))) +
          xlab('Pseudotime') + ylab('Expression') + 
          labs(col = colnames(design)[1]) +
          facet_wrap(~Sample, scales=a)
      } else {
        ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=Expr, color=Variable, group = Sample), alpha=alpha, size=size)  +
          geom_line(data=ld, aes(x=pseudotime, y=Expr, color=Variable, group = Sample),alpha=1, size=.5) +
          theme_classic() +
          scale_color_viridis(discrete = TRUE, direction = -1)+
          ggtitle(paste0(sub(':.*','',gene),',adjP=', formatC(testptObj$fdr[gene], format = "e", digits = 2))) +
          xlab('Pseudotime') + ylab('Expression') + 
          labs(col = colnames(design)[1])  
      }
    } else {
      ggplot() + 
        geom_line(data=ld, aes(x=pseudotime, y=Expr, color=Variable,group=Sample),alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1)+
        ggtitle(paste0(sub(':.*','',gene),',adjP=', formatC(testptObj$fdr[gene], format = "e", digits = 2))) +
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(design)[1])
    }
  } else {
    print('plotting multiple genes ...')
    pdlist <- ldlist <- list()
    for (g in gene){
      pd <- data.frame(Expr = expression[g, ], cell = cellano[,1], Sample = cellano[,2], Variable = design[match(cellano[,2], rownames(design)), variable])
      pd = cbind(pd, pseudotime = pseudotime[pd$cell], g = g)
      pdlist[[g]] <- pd
      linedlist <- lapply(unique(pd$Sample), function(p){
      tmpdf = data.frame(Expr = pred[g, gsub(':.*','', colnames(pred)) == p], pseudotime=pseudotime[gsub(':.*','', names(pseudotime)) == p], Sample = p, Variable = design[rownames(design) == p, variable], g = g)     
      })
      ld = do.call(rbind, linedlist)
      
      ldlist[[g]] <- ld
    }
    ld <- do.call(rbind, ldlist)
    ld = ld[order(ld$pseudotime),]
    pd <- do.call(rbind, pdlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    if (plot.points){
      ggplot() + 
        geom_point(data=pd, aes(x=pseudotime, y=Expr, color=Variable, group = Sample), alpha=alpha, size=size)  +
        geom_line(data=ld, aes(x=pseudotime, y=Expr, color=Variable, group = Sample), alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1) +
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(design)[1]) + 
        facet_wrap(~g, scales = a)
    } else {
      ggplot() + 
        geom_line(data=ld, aes(x=pseudotime, y=Expr, color=Variable, group = Sample), alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1) +
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(design)[1]) + 
        facet_wrap(~g, scales = a)
    }
  }
}


