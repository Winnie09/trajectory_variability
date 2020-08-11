plotGene_Time <- function(testptObj, Gene, Mat, Pseudotime, Cellanno, Design,  Alpha=1, Size=0.5, PlotPoints = FALSE, FreeScale = FALSE, BySample = FALSE){
  ## Mat: gene by cell Matrix
  ## testptObj: the output of function testpt() which is a list containing fdr, etc..
  library(ggplot2)
  library(gridExtra)
  library(viridis)   
  Order <- data.frame(Cell = names(Pseudotime), Pseudotime = Pseudotime, stringsAsFactors = FALSE)
  colnames(Cellanno) <- c('Cell', 'Sample')
  Mat <- Mat[, Cellanno[,1]]
  knotnum <- testptObj$knotnum
  knotnum[knotnum==0] <- 1  ## in case the fitting of line would cause bugs
  
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
    linedlist <- lapply(unique(testptObj$cellanno[,2]), function(p){
      tmpcell <- testptObj$cellanno[testptObj$cellanno[,2]==p,1]
      tmpdf <- data.frame(Expr=testptObj$predict.values[Gene,tmpcell], Pseudotime=testptObj$pseudotime[tmpcell], Sample=p, Variable = Design[rownames(Design) == p, 1])
    })
    ld = do.call(rbind, linedlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    if (PlotPoints){
      if (BySample){
        ggplot() + 
        geom_point(data=pd, aes(x=Pseudotime, y=Expr, color=Sample), alpha=Alpha, size=Size)  +
        geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Sample),alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1) +
        ggtitle(paste0(sub(':.*','',Gene),',adjP=', formatC(testptObj$fdr[Gene], format = "e", digits = 2))) +
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(Design)[1]) +
        facet_wrap(~Sample, scales=a)
      } else {
        ggplot() + 
        geom_point(data=pd, aes(x=Pseudotime, y=Expr, color=Sample), alpha=Alpha, size=Size)  +
        geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Sample),alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1)+
        ggtitle(paste0(sub(':.*','',Gene),',adjP=', formatC(testptObj$fdr[Gene], format = "e", digits = 2))) +
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(Design)[1])  
      }
    } else {
      ggplot() + 
        geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Sample),alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1)+
        ggtitle(paste0(sub(':.*','',Gene),',adjP=', formatC(testptObj$fdr[Gene], format = "e", digits = 2))) +
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
    linedlist <- lapply(unique(testptObj$cellanno[,2]), function(p){
      tmpcell <- testptObj$cellanno[testptObj$cellanno[,2]==p,1]
      tmpdf <- data.frame(Expr=testptObj$predict.values[g,tmpcell], Pseudotime=testptObj$pseudotime[tmpcell], Sample=p, Variable = Design[rownames(Design) == p, 1],g=g)
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
        geom_point(data=pd, aes(x=Pseudotime, y=Expr, color=Sample), alpha=Alpha, size=Size)  +
        geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Sample), alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1) +
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(Design)[1]) + 
        facet_wrap(~g, scales = a)
    } else {
      ggplot() + 
        geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Sample), alpha=1, size=.5) +
        theme_classic() +
        scale_color_viridis(discrete = TRUE, direction = -1) +
        xlab('Pseudotime') + ylab('Expression') + 
        labs(col = colnames(Design)[1]) + 
        facet_wrap(~g, scales = a)
    }
  }
}
