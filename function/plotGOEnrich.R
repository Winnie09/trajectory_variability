plotGOEnrich <- function(goRes, n = 5, sortByFDR = TRUE, fdr.cutoff = 0.05, fc.cutoff = 2){
  ## goRes: output from GOEnrich. A list of GO enrichment result. Length same as number of clusters.
  ## n: number of top GO terms shown for each cluster
  ## sortByFDR: if TRUE (default), sort GO terms by FDR first, and then negative FC. If FALSE,  soft by negative FC and then FDR. 
  d <- lapply(names(goRes), function(i){  
    cbind(Cluster = i, goRes[[i]])
  })
  d <- do.call(rbind, d)    
  if (sortByFDR){
    d <- d[order(d$FDR,-d$FC),]  
  } else {
    d <- d[order(-d$FC, d$FDR),]
  }
  cd <- do.call(rbind,sapply(sort(unique(d$Cluster)),function(i) {
    tmp <- d[d$Cluster==i,]
    tmp <- tmp[tmp$FDR < fdr.cutoff & tmp$FC > fc.cutoff,,drop=F]
    if (nrow(tmp) > 0) {
      tmp[1:min(n,nrow(tmp)),]  
    } else {
      NULL
    }
  },simplify = F))
  ut <- unique(cd$Term)
  d <- d[d$Term %in% ut,c('Cluster','Term','FDR','FC')]
  d$Cluster <- as.numeric(as.character(d$Cluster))
  
  d <- d[d$FDR < fdr.cutoff & d$FC > fc.cutoff,,drop=F]
  d$enc <- TRUE
  fulld <- expand.grid(1:length(goRes),unique(d$Term))
  fulld <- fulld[!paste0(fulld[,1],fulld[,2]) %in% paste0(d[,1],d[,2]),]
  fulld <- data.frame(Cluster=fulld[,1],Term=as.character(fulld[,2]),FDR=1,FC=0,enc=FALSE)
  pd <- rbind(d,fulld)
  
  # dmat <- dcast(d,Term~Cluster)
  # rownames(dmat) <- dmat[,1]
  # dmat <- as.matrix(dmat[,-1,drop=F])
  # dmat <- is.na(dmat)
  # 
  # v <- sapply(1:nrow(dmat), function(i) which.min(dmat[i,]))
  # names(v) <- rownames(dmat)
  # 
  # pd <- melt(dmat)
  # colnames(pd) <- c('Term','Cluster','enc')
  # # pd$Term <- factor(as.character(pd$Term),levels=rownames(dmat)[hclust(dist(dmat))$order])
  pd$Term <- factor(as.character(pd$Term),levels=names(v[rev(order(v))]))
  pd$enc <- ifelse(pd$enc,'Significant','Non-significant')
  pd$Cluster <- as.character(pd$Cluster)
  library(RColorBrewer)
  p <- ggplot(pd,aes(x=Cluster,y=Term,fill=enc)) + geom_tile() + theme_classic() + 
    scale_fill_manual(values=c('darkblue', 'orange'))+
    theme(legend.position = 'right')+
    scale_y_discrete(position = "right")
  return(p)
}

