plotGene <- function(Gene, Mat, Order, Cellanno, Design, Stat, Alpha=1, Size=0.5, PlotPoints = FALSE, FreeScale = FALSE){
  ## Mat: gene by cell Matrix
  ## Order: dataframe, columns are: Cell, Pseudotime
  library(ggplot2)
  library(gridExtra)
  library(viridis)   
  colnames(Order) <- c('Cell', 'Pseudotime')
  colnames(Cellanno) <- c('Cell', 'Sample')
  Mat <- Mat[, Cellanno[,1]]
  if (FreeScale) {
    a  = 'free'
  } else {
    a  = 'fixed'
  }
  if (length(Gene) == 1){
    print('plotting one gene ...')
    pd <- data.frame(Expr = Mat[Gene, ], Cell = Cellanno[,1], Sample = Cellanno[,2], Variable = Design[match(Cellanno[,2], rownames(Design)), 1])
    pd = cbind(pd, Pseudotime = Order[match(pd$Cell, Order$Cell),'Pseudotime'])
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    linedlist <- lapply(unique(pd$Sample), function(p){
      tMat = Mat[g,grepl(p,colnames(Mat)),drop=F]
      trainX = Order[match(colnames(tMat), Order$Cell),'Pseudotime']      ### use time 
      pred <- get_spline_fit(tMat, trainX=trainX, fit.min=min(Order$Pseudotime), fit.max=max(Order$Pseudotime), num.base = 5)
      tmpdf <- data.frame(Expr=pred[1,], Pseudotime=seq(min(Order$Pseudotime),max(Order$Pseudotime),length.out = 1000), Sample=p, Variable = Design[rownames(Design) == p, 1])
    })
    ld = do.call(rbind, linedlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    if (PlotPoints){
      ggplot() + 
      geom_point(data=pd, aes(x=Pseudotime, y=Expr, color=Variable), alpha=Alpha, size=Size)  +
      geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Variable),alpha=1, size=.5) +
      theme_classic() +
      scale_color_viridis(discrete = TRUE, direction = -1)+
      ggtitle(paste0(sub(':.*','',Gene),',adjP=', formatC(Stat[g,'adj.P.Val'], format = "e", digits = 2))) +
      xlab('Pseudotime') + ylab('Expression') + 
      labs(col = colnames(Design)[1])
    } else {
      ggplot() + 
        geom_point(data=pd, aes(x=Pseudotime, y=Expr, color=Variable), alpha=Alpha, size=Size)  +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1)+
        ggtitle(paste0(sub(':.*','',Gene),',adjP=', formatC(Stat[g,'adj.P.Val'], format = "e", digits = 2))) +
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(Design)[1])
    }
  } else {
    print('plotting multiple genes ...')
    pdlist <- ldlist <- list()
    for (g in Gene){
      pd <- data.frame(Expr = Mat[g, ], Cell = Cellanno[,1], Sample = Cellanno[,2], Variable = Design[match(Cellanno[,2], rownames(Design)), 1])
      pd = cbind(pd, Pseudotime = Order[match(pd$Cell, Order$Cell),'Pseudotime'], g = g)
      pdlist[[g]] <- pd
      linedlist <- lapply(unique(pd$Sample), function(p){
        tMat = Mat[g,grepl(p,colnames(Mat)),drop=F]
        trainX = Order[match(colnames(tMat), Order$Cell),'Pseudotime']      ### use time 
        pred <- get_spline_fit(tMat, trainX=trainX, fit.min=min(Order$Pseudotime), fit.max=max(Order$Pseudotime), num.base = 3)
        tmpdf <- data.frame(Expr=pred[1,], Pseudotime=seq(min(Order$Pseudotime),max(Order$Pseudotime),length.out = 1000), Sample=p, Variable = Design[rownames(Design) == p, 1], g = g)
      })
      ld = do.call(rbind, linedlist)
      ldlist[[g]] <- ld
    }
    ld <- do.call(rbind, ldlist)
    pd <- do.call(rbind, pdlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    if (PlotPoints){
      ggplot() + 
        geom_point(data=pd, aes(x=Pseudotime, y=Expr, color=Variable), alpha=Alpha, size=Size)  +
        geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Variable), alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1)+
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(Design)[1]) + 
        facet_wrap(~g, scales = a)
    } else {
      ggplot() + 
        geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Variable), alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1)+
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(Design)[1]) + 
        facet_wrap(~g, scales = a)
    }
      
  }
}

