# ori.design = design
# test.pattern = 'overall'
# test.position = 'all'
# maxknotallowed=10; EMmaxiter=1000; EMitercutoff=0.01; verbose=F; ncores=1; model = 3
# test.pattern = 'overall'
fitpt <- function(expr, cellanno, pseudotime, design, maxknotallowed=10, EMmaxiter=1000, EMitercutoff=0.01, verbose=F, ncores=1, model = 3, knotnum = NULL) {
  # set.seed(12345)
  suppressMessages(library(Matrix))
  suppressMessages(library(parallel))
  suppressMessages(library(splines))
  suppressMessages(library(matrixcalc))
  ## expr: gene by cell matrix, entires and log-transformed expression
  ## pseudotime: a numeric vecotor of pseudotime, the names are ordered cell names
  ## design: design: sample by feature design matrix. rownames are sample names. first column is 1, second column is the group partition, currently only one variable.
  ## cellanno: dataframe, first column is cell names, second column is sample names.
  pseudotime <- pseudotime[colnames(expr)]
  cellanno <- cellanno[match(colnames(expr),cellanno[,1]),]
  sname <- sapply(row.names(design),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  design = as.matrix(design)
  
  philist <- lapply(0:maxknotallowed,function(num.knot) {
    if (num.knot==0) {
      phi <- cbind(1,bs(pseudotime))
    } else {
      knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      phi <- cbind(1,bs(pseudotime,knots = knots))  
    }
  })
  names(philist) <- as.character(0:maxknotallowed)
  
  maxknot <- 0
  testpos <- 1
  
  while (testpos & maxknot < maxknotallowed) {
    maxknot <- maxknot + 1
    phi <- philist[[as.character(maxknot)]]
    testpos <- mean(sapply(names(sname), function(ss) {
      is.positive.definite(crossprod(phi[sname[[ss]],,drop=F]))
    })) == 1
  }
  maxknot <- maxknot - 1
  sexpr <- sapply(names(sname),function(ss) expr[,sname[[ss]],drop=F],simplify = F)
  
  ## automatically select knotnum for each gene
  if (is.null(knotnum)){
    bicfunc <- function(num.knot) {
      phi <- philist[[as.character(num.knot)]]
      ll <- sapply(names(sname), function(ss) {
        phiss <- phi[sname[[ss]],,drop=F]
        dif2 <- sexpr[[ss]] - sexpr[[ss]] %*% (phiss %*% chol2inv(chol(crossprod(phiss)))) %*% t(phiss)
        dif2 <- rowSums(dif2 * dif2)
        s2 <- dif2/(length(sname[[ss]])-ncol(phi))
        log(2*pi*s2)*nrow(phiss) + dif2/s2
      })
      if (is.vector(ll)){
        sum(ll,na.rm=T) + log(nrow(phi))*((ncol(phi)+1)*sum(!is.na(ll)))
      } else {
        rowSums(ll,na.rm=T) + log(nrow(phi))*((ncol(phi)+1)*rowSums(!is.na(ll)))
      }    
    }
    
    if (ncores!=1) {
      bic <- mclapply(0:maxknot,bicfunc,mc.cores=ncores)
      bic <- do.call(cbind,bic)
    } else {
      bic <- sapply(0:maxknot,bicfunc)
    }
    
    rm('sexpr')
    if (is.vector(bic)){
      knotnum <- c(0:maxknot)[which.min(bic)]
    } else {
      knotnum <- c(0:maxknot)[apply(bic,1,which.min)]
    }
      
    names(knotnum) <- rownames(expr)
  } else {
    knotnum <- knotnum[rownames(expr)]
  }
  
  sfit <- function(num.knot) {
    gid <- names(which(knotnum==num.knot))
    sexpr <- expr[gid,,drop=F] ## !!! double check, should be list len =S
    phi <- philist[[as.character(num.knot)]]
    phicrossprod <- apply(phi,1,tcrossprod)
    phicrossprod <- sapply(names(sname),function(ss) phicrossprod[,sname[[ss]]],simplify = F)
    phi <- sapply(names(sname),function(ss) phi[sname[[ss]],],simplify = F) 
    
    if (model == 0){
      xs <- sapply(row.names(design), function(i){
        as.matrix(1, nrow = 1, ncol = 1)
      }, simplify = FALSE)
      phi <- sapply(phi, function(i){
        i[,1,drop=F]
      }, simplify = FALSE)
    } else if (model == 1) {
      xs <- sapply(row.names(design), function(i) {
        kronecker(diag(num.knot + 4), design[i, 1, drop = F]) 
      }, simplify = F)
    } else if (model == 2) {
      xs <- sapply(row.names(design), function(i) {
        tmp <- kronecker(diag(num.knot + 4), design[i, ]) 
        tmp <- tmp[-seq(4, nrow(tmp), 2), ] 
      }, simplify = F)
    } else if (model == 3) {
      xs <- sapply(row.names(design), function(i) {
        kronecker(diag(num.knot + 4), design[i, ]) 
      }, simplify = F)
    }
    
    as <- names(phi)
    nb <- ncol(phi[[1]])
    cn <- sapply(as,function(s) nrow(phi[[s]]))
    
    phiphi <- sapply(as,function(s) {
      t(phi[[s]]) %*% phi[[s]]
    },simplify = F)
    phiX <- sapply(as, function(s){
      phi[[s]] %*% t(xs[[s]])
    },simplify = F)
    
    phiXTphiX <- sapply(as, function(s){
      t(phiX[[s]]) %*% (phiX[[s]])
    },simplify = F)
    
    ## initialize
    B1 <- Reduce('+',phiXTphiX)
    B2 <- Reduce('+',sapply(as,function(s) t(phiX[[s]]) %*% t(sexpr[, cellanno[,2]==s, drop=F]),simplify = F))
    B <- t(solve(B1) %*% B2)
    
    indfit <- sapply(as,function(s) {
      sexpr[, cellanno[,2]==s, drop=F] %*% (phi[[s]] %*% chol2inv(chol(crossprod(phi[[s]]))))
    },simplify = F)
    
    s2 <- matrix(sapply(as,function(s) {
      tmp <- sexpr[, cellanno[,2]==s, drop=F]-indfit[[s]] %*% t(phi[[s]])
      n <- ncol(tmp)
      (rowMeans(tmp*tmp)-(rowMeans(tmp))^2)*(n)/(n-1)
    }),nrow=length(gid),dimnames=list(gid,as))
    alpha <- 1/apply(s2,1,var)*rowMeans(s2)^2+2
    eta <- (alpha-1)*rowMeans(s2)
    
    diffindfit <- sapply(as,function(s) {
      indfit[[s]] - B %*% xs[[s]]
    },simplify = F)
    
    if (model==0) {
      omega <- matrix(sapply(rownames(sexpr),function(g) {
        m <- sapply(diffindfit,function(i) i[g,,drop=F])
        var(m) + 0.01
      }),ncol=1,dimnames = list(rownames(sexpr),NULL))
    } else {
      omega <- t(sapply(rownames(sexpr),function(g) {
        m <- sapply(diffindfit,function(i) i[g,,drop=F])
        m <- m-rowMeans(m)
        tcrossprod(m)/(ncol(m)-1) + diag(nrow(m)) * 0.01
      }))
    }
    
    iter <- 0
    gidr <- rownames(sexpr)
    all <- matrix(-Inf, nrow=nrow(sexpr),ncol=1,dimnames = list(rownames=gidr))
    etalist <- alphalist <- omegalist <- Nlist <- Jslist <- list()
    while (iter < EMmaxiter && length(gidr) > 0) {
      sexpr_phibx <- sapply(as,function(s) {
        sexpr[, cellanno[,2]==s, drop=F][gidr,,drop=F]-B[gidr,] %*% t(phiX[[s]])
      },simplify = F)
      
      L <- sapply(as,function(s) {
        rowSums(sexpr_phibx[[s]] * sexpr_phibx[[s]])
      },simplify = F)
      
      Jsolve <- sapply(as,function(s) {
        sapply(gidr,function(g) {
          chol2inv(chol(phiphi[[s]] + chol2inv(chol(matrix(omega[g,],nrow=nb)))))
        },simplify = F)
      },simplify = F)
      
      K <- sapply(as,function(s) {
        tcrossprod(t(phi[[s]]),sexpr_phibx[[s]])
      },simplify = F)
      
      L2eKJK <- sapply(as,function(s) {
        2*eta[gidr]+L[[s]]-colSums(K[[s]][rep(1:nb,nb),,drop=FALSE] * K[[s]][rep(1:nb,each=nb),,drop=FALSE] * sapply(Jsolve[[s]],as.vector))
      },simplify = F)
      A <- matrix(sapply(as,function(s) {
        log(L2eKJK[[s]]/2)-digamma(alpha[gidr]+cn[s]/2)
      }),nrow=length(gidr),dimnames=list(gidr,as))
      N <- matrix(sapply(as,function(s) {
        (2*alpha[gidr]+cn[s])/L2eKJK[[s]]
      }),nrow=length(gidr),dimnames=list(gidr,as))
      
      ll <- rowSums(matrix(sapply(as,function(s) {
        dv <- sapply(gidr,function(g) {
          det(matrix(omega[g,],nrow=nb))/det(Jsolve[[s]][[g]])
        })
        alpha[gidr]*log(2*eta[gidr])+lgamma(cn[s]/2+alpha[gidr])-cn[s]*log(pi)/2-lgamma(alpha[gidr])-log(dv)/2-(cn[s]/2+alpha[gidr])*log(L2eKJK[[s]])
      }),nrow=length(gidr),dimnames = list(rownames=gidr,colnames=as)))
      
      JK <- sapply(as,function(s) {
        t(rowsum(t(t(sapply(Jsolve[[s]],as.vector))*t(K[[s]])[,rep(1:nb,nb),drop=FALSE]),rep(1:nb,each=nb)))
      },simplify = F)
      
      B1 <- Reduce('+', lapply(as, function(s){
        tcrossprod(N[,s],as.vector(phiXTphiX[[s]]))
      }))
      rownames(B1) <- gidr
      B2 <- Reduce('+', lapply(as, function(s){
        N[,s] * t(t(phiX[[s]]) %*%  t(sexpr[ ,cellanno[,2]==s, drop=F][gidr,] - t(phi[[s]] %*% t(JK[[s]]))))
      }))
      
      np <- nrow(xs[[1]])
      B[gidr,] <- t(sapply(gidr, function(g){  ## each column is a gene's all betas
        chol2inv(chol(matrix(B1[g,],nrow=np))) %*% B2[g,]
      }))
      
      ## M -step: 
      omega[gidr,] <- Reduce('+',sapply(as,function(s) {
        t(sapply(Jsolve[[s]],as.vector)) + N[,s]*JK[[s]][,rep(1:nb,nb)] * JK[[s]][,rep(1:nb,each=nb)]
      },simplify=F))/length(as)
      
      eta[gidr] <- sapply(gidr,function(g) {
        meanN <- mean(N[g,])
        meanA <- mean(A[g,])
        uniroot(function(eta) {digamma(eta * meanN)-log(eta)+meanA},c(1e-10,1e10))$root
      })
      alpha[gidr] <- eta[gidr] * rowMeans(N)
      
      iter <- iter + 1
      
      llv <- all[,ncol(all)]
      llv[gidr] <- ll
      all <- cbind(all,llv)
      gidr <- names(which(all[,ncol(all)] - all[,ncol(all)-1] > EMitercutoff))
      etalist[[iter]] <- eta
      alphalist[[iter]] <-alpha
      omegalist[[iter]] <-omega
      Nlist[[iter]] <- N
      Jslist[[iter]] <- Jsolve
      rm(list = c('L','Jsolve', 'K'))
    }
    # print(table(apply(all,1,function(i) mean(diff(i) >= 0))))  
    return(list(beta = B, alpha = alpha, eta = eta, omega = omega, logL = all))
  }
  
  if (ncores!=1) {
    allres <- mclapply(unique(knotnum),sfit,mc.cores=ncores)
  } else {
    allres <- sapply(unique(knotnum),sfit, simplify = FALSE)
  }
  
  para <- list()
  for (i in 1:length(allres)) {
    for (j in row.names(allres[[i]][[1]])) {
      para[[j]] <- list(beta=allres[[i]][[1]][j,],
                        alpha=allres[[i]][[2]][j],
                        eta=allres[[i]][[3]][j],
                        omega=allres[[i]][[4]][j,],
                        ll=allres[[i]][[5]][j,ncol(allres[[i]][[5]])])
    }
  }
  list(parameter=para[rownames(expr)],knotnum=knotnum[rownames(expr)])
}


