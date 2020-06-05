plotGeneGenderDiff <- function(GeneSet, mat, order, sample, design, statistics, GroupNames=c('1','0')){
  ## mat: gene by cell matrix
  ## order: dataframe, columns are: Cell, Pseudotime
  library(ggplot2)
  library(gridExtra)
  plist <- list()
  for (g in GeneSet){
    pd1 = mat[g, sample %in% rownames(design)[design==1]]
    pd1 = data.frame(Expr=pd1, Cell=names(pd1), Patient = gsub(':.*','',names(pd1) ), Group=GroupNames[1])
    pd2 = mat[g, sample %in% rownames(design)[design==0]]
    pd2 = data.frame(Expr=pd2, Cell=names(pd2), Patient = gsub(':.*','',names(pd2) ), Group=GroupNames[2])
    pd = rbind(pd1, pd2)
    pd = cbind(pd, Pseudotime = order[match(pd$Cell, order$Cell),'Pseudotime'])
    linedlist <- lapply(unique(pd$Patient), function(p){
      tmat = mat[g,grepl(p,colnames(mat)),drop=F]
      trainX = order[match(colnames(tmat), order$Cell),'Pseudotime']      ### use time 
      pred <- get_spline_fit(tmat, trainX=trainX, fit.min=min(order$Pseudotime), fit.max=max(order$Pseudotime), num.base = 3)
      tmpdf <- data.frame(Expr=pred[1,], Pseudotime=seq(min(order$Pseudotime),max(order$Pseudotime),length.out = 1000), Patient=p, Group = ifelse(design[names(unlist(sapply(rownames(design), function(i) grep(i, p)))),1]==1,GroupNames[1],GroupNames[2]))
    })
    ld = do.call(rbind, linedlist)
    ld[,'Group'] = as.factor(ld[,'Group'])
    plist[[g]] <- ggplot() + geom_point(data=pd, aes(x=Pseudotime, y=Expr, color=Patient), alpha=.1, size=.2)  + 
      geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Patient),alpha=1, size=.5) +
      theme_classic() +
      ggtitle(paste0(sub(':.*','',g),',adjP=', round(statistics[g,'adj.P.Val'],4),',f=',round(statistics[g,'FStat'],2))) + theme(legend.position = 'none', plot.title = element_text(size=12)) + scale_color_manual(values=c(rep('darkblue',4),rep('orange',4)))
  }
  return(plist)
}

