predict_fitting <- function(testObj, gene = NULL, test.type = 'time'){
  ## make the cells order according to pseudotime order
  ## fitting
  if ('expr.ori' %in% names(testObj)) expr <- testObj$expr.ori else  expr <- testObj$expr
  knotnum = testObj$knotnum[gene]
  design = testObj$design
  cellanno = testObj$cellanno
  pseudotime = testObj$pseudotime[colnames(expr)]
  if (is.null(gene)) gene <- rownames(expr)
  philist <- lapply(sort(unique(knotnum)), function(num.knot) {
    if (num.knot==0) {
      # phi <- cbind(1,bs(pseudotime))
      phi <- bs(pseudotime, intercept = TRUE)
    } else {
      knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      # phi <- cbind(1,bs(pseudotime,knots = knots))  
      phi <- bs(pseudotime,knots = knots, intercept = TRUE)
    }
  })
  
  names(philist) <- as.character(sort(unique(knotnum)))
  as <- row.names(design)
  sname <- sapply(as,function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  
  pred <- lapply(unique(knotnum), function(num.knot){
    genesub <- names(knotnum)[knotnum == num.knot]
    B <- t(sapply(genesub, function(g){
      testObj$parameter[[g]]$beta
    }))
    
    omega <- t(sapply(genesub, function(g){
      testObj$parameter[[g]]$omega
    }))
     
    phi <- philist[[as.character(num.knot)]]
    phipred <- phi
    phi <- sapply(as,function(ss) phi[sname[[ss]],],simplify = F)
    
    if (toupper(test.type)  == 'TIME') {
      xs <- sapply(row.names(design), function(i) {
        kronecker(diag(num.knot + 4), design[i, 1, drop = F])
      }, simplify = F)
    } else if (toupper(test.type) == 'VARIABLE') {
      xs <- sapply(row.names(design), function(i) {
        kronecker(diag(num.knot + 4), design[i, ])
      }, simplify = F)
    }
    phiphi <- sapply(as,function(s) {
      t(phi[[s]]) %*% phi[[s]]
    },simplify = F)
    phiX <- sapply(as, function(s){
      phi[[s]] %*% t(xs[[s]])
    },simplify = F)
    
    predtmp <- sapply(as, function(s){
      sexpr <- expr[genesub,,drop=F] 
      sexpr_phibx <- sexpr[genesub, cellanno[,2]==s, drop=F]-B[genesub,] %*% t(phiX[[s]])
      
      nb <- num.knot + 4
      oinv <- sapply(genesub,function(g) {
        chol2inv(chol(matrix(omega[g,,drop=F],nrow=nb)))
      },simplify = F)
      
      Jchol <- sapply(genesub,function(g) {
        chol(phiphi[[s]] + oinv[[g]])
      },simplify = F)
      
      Jsolve <- sapply(genesub,function(g) {
        chol2inv(Jchol[[g]])
      })
      K <- tcrossprod(t(phi[[s]]),sexpr_phibx)
      JK <- rowsum((Jsolve*K[rep(1:nb,nb),,drop=FALSE]),rep(1:nb,each=nb)) ## u's poterior mean
      ## 20210515 <<<<<<<<<<<<<<<<<<<<
      if (toupper(test.type) == 'VARIABLE') {
      tmpfit <- t(phipred %*% JK) 
      tmpfit <- tmpfit[, colnames(tmpfit) %in% cellanno[cellanno[,2] == s,1], drop = FALSE] 
      } else {
        t(phipred %*% JK) 
      }
      ## 20210515 >>>>>>>>>>>>>>>>>>>
    }, simplify = F)
    if (toupper(test.type) == 'VARIABLE') {
      predtmp <- do.call(cbind, predtmp)
    } else {
      return(predtmp)
    }
  })
  
  if ('populationFit' %in% names(testObj))  populationFit = testObj$populationFit else 
    populationFit <- getPopulationFit(testObj,gene, type = testObj$test.type)
  
  if (toupper(test.type) == 'TIME'){
    predfit <- lapply(pred, function(i){
      for (j in 1:length(i)){
        if (ncol(populationFit) < ncol(i[[j]])){ ## cell proportion test case
          i[[j]] <- i[[j]][, cellanno[cellanno[,2] == names(i)[j], 1], drop = FALSE]
          i[[j]] <- i[[j]] + populationFit[rownames(i[[j]]), , drop = FALSE]
          
        } else { ## differential gene test case
          i[[j]] <- i[[j]] + populationFit[rownames(i[[j]]), , drop = FALSE]
          i[[j]] <- i[[j]][, cellanno[cellanno[,2] == names(i)[j], 1], drop = FALSE]
        }
          
      }
      tmp <- do.call(cbind, i)
      tmp[, names(pseudotime), drop = FALSE]
    })
    predfit <- do.call(rbind, predfit)
    return(predfit)
  } else {
    pred <- do.call(rbind, pred)
    pred <- pred[gene, colnames(expr), drop=FALSE]
    l <- lapply(populationFit, function(i){
      pred + i[gene, pseudotime , drop=F] 
    })
    return(l)
  }
}


