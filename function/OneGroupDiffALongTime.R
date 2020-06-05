OneGroupDiffALongTime <- function(expr,sample,dr=NULL,pseudotime=NULL,permutation=FALSE, fstatonly=FALSE,permutime=10000, num.base = 3) {
  n <- pseudotime
  pseudotime <- 1:length(pseudotime)
  names(pseudotime) <- n
  
  knots = seq(min(pseudotime),max(pseudotime),length.out=num.base+2)[2:(num.base+1)]
  base <- lapply(unique(sample),function(sg) {
    n <- names(sample)[sample==sg]
    cbind(1,bs(pseudotime[n],knots = knots))
  })
  names(base) <- unique(sample)
  singularid <- unique(unlist(lapply(base,function(i) which(apply(i,2,function(i) length(unique(i)))==1))))
  nonsigularid <- setdiff(1:ncol(base[[1]]),setdiff(singularid,1))
  coef <- lapply(unique(sample),function(sg) {
    n <- names(sample)[sample==sg]
    est <- expr[,n] %*% base[[sg]][,nonsigularid] %*% t(chol2inv(chol(crossprod(base[[sg]][,nonsigularid]))))
    est[,-1]
  })
  coef <- do.call(cbind,coef) 
  
  comdesign <- kronecker(rep(1,length(unique(sample))),diag(length(nonsigularid)-1))
  res <- topTable(eBayes(lmFit(coef,comdesign)),number=nrow(coef),sort.by='none')
  
  fstat <- res[row.names(expr),'F']
  adjPval <- res[row.names(expr),'adj.P.Val']
  if (fstatonly) {
    res <- data.frame(Fstat=fstat)
    row.names(res) <- row.names(expr)
    res[order(-res[,1]),,drop=F] ###
  } else {
    if (!permutation){
      res <- res[rownames(expr), c('F','P.Value','adj.P.Val')]
      colnames(res)[1] <- 'FStat'
      res[order(res[,3],-res[,1]),]
    } else {
      print('Permutation version will be developed.')
      # permuf <- mclapply(1:permutime,function(permu) {
      #   # limma may set seed by itself. Need to resetseed here
      #   set.seed(permu)
      #   #presetting
      #   simuid <- unlist(sapply(unique(sample),function(i) sample(names(sample)[sample==i],replace=T)))
      #   perdr <- dr[simuid,]
      #   perexpr <- expr[,simuid]
      #   persample <- sample[simuid]
      #   names(persample) <- colnames(perexpr) <- row.names(perdr) <- paste0('cell',1:nrow(perdr))
      #   # permudesign <- design[sample(1:nrow(design)),,drop=F]
      #   # row.names(permudesign) <- row.names(design)
      #   
      #   #redo pseudotime for bootstrapped samples
      #   ord <- lapply(2:10,function(i) {
      #     tmp <- TSCANorder(exprmclust(t(perdr),reduce = F,clusternum=i,clustermethod='kmeans'),orderonly = T)  
      #     n <- 1:length(tmp)
      #     names(n) <- tmp
      #     n
      #   })
      #   cv <- sapply(ord,function(i) {
      #     if (length(i) != length(pseudotime)) { ## there exists branches
      #       0
      #     } else {
      #       names(i) <- simuid[match(names(i),paste0('cell',1:nrow(perdr)))]
      #       int <- intersect(names(i),names(pseudotime))
      #       cor(i[int],pseudotime[int])
      #     }
      #   })
      #   pertime <- ord[[which.max(abs(cv))]]
      #   if (cv[which.max(abs(cv))] < 0) pertime <- rev(max(pertime) - pertime + 1)
      #   #fitting splines
      #   knots = seq(min(pertime),max(pertime),length.out=num.base+2)[2:(num.base+1)]
      #   base <- lapply(row.names(permudesign),function(sg) {
      #     n <- names(persample)[persample==sg]
      #     cbind(1,bs(pertime[n],knots = knots))
      #   })
      #   names(base) <- unique(sample) ##
      #   singularid <- unique(unlist(lapply(base,function(i) which(apply(i,2,function(i) length(unique(i)))==1))))
      #   nonsigularid <- setdiff(1:ncol(base[[1]]),setdiff(singularid,1))
      #   coef <- lapply(unique(sample),function(sg) {
      #     n <- names(persample)[persample==sg]
      #     perexpr[,n] %*% base[[sg]][,nonsigularid] %*% t(chol2inv(chol(crossprod(base[[sg]][,nonsigularid]))))
      #   })
      #   coef <- do.call(cbind,coef)
      #   coef <- coef[,as.vector(sapply(1:length(nonsigularid),function(i) i+c(0:(length(unique(sample))-1))*length(nonsigularid)))]
      #   compermudesign <- kronecker(diag(length(nonsigularid)),cbind(1,c(rep(1, length(unique(sample)) ))))
      #   coefid <- as.vector(t(sapply(which(colnames(design) %in% var) + 1, function(i) i + 0:(length(nonsigularid) - 1) * (ncol(design) + 1)))) ## to modify
      #   res <- topTable(eBayes(lmFit(coef,compermudesign)),coef=coefid,number=nrow(coef),sort.by='none') ## to modify
      #   res[row.names(expr),'F']
      # },mc.cores=detectCores())
      # permuf <- do.call(cbind,permuf)
      # pval <- sapply(1:length(fstat),function(i) mean(permuf[i,] > fstat[i]))
      # fdr <- p.adjust(pval,method='fdr')
      # res <- data.frame(FStat=fstat,P.Value=pval,adj.P.Val=fdr)
      # row.names(res) <- row.names(expr)
      # res <- res[order(res[,2],-res[,1]),]                   ######
      # list(res=res,permuf=permuf)
    }
  }
}


