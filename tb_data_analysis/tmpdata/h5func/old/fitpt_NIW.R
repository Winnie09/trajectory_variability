# ori.design = design
# test.pattern = 'overall'
# test.position = 'all'
# maxknotallowed=10; EMmaxiter=1000; EMitercutoff=0.01; verbose=F; ncores=1; model = 3
# test.pattern = 'overall'
fitpt <- function(path, design, maxknotallowed=10, EMmaxiter=100, EMitercutoff=0.1, verbose=F, ncores=1, model = 3, knotnum = NULL) {
  suppressMessages(library(Matrix))
  suppressMessages(library(parallel))
  suppressMessages(library(splines))
  suppressMessages(library(matrixcalc))
  suppressMessages(library(rhdf5))
  ## expr: gene by cell matrix, entires and log-transformed expression
  ## pseudotime: a numeric vecotor of pseudotime, the names are ordered cell names
  ## design: design: sample by feature design matrix. rownames are sample names. first column is 1, second column is the group partition, currently only one variable.
  ## cellanno: dataframe, first column is cell names, second column is sample names.
  
  samp <- h5ls(path,recursive=F)$name
  pseudotime <- unlist(sapply(samp,function(s) as.vector(h5read(path,paste0(s,'/pseudotime')))))
  sname <- sapply(samp,function(s) as.vector(h5read(path,paste0(s,'/barcode'))))
  names(pseudotime) <- unlist(sname)
  pseudotime <- sort(pseudotime)
  gn <- h5read(path,paste0(samp[1],'/feature'))
  
  exprreadfunc <- function(s,geneid) {
    e <- h5read(path,paste0(s,'/expr'),index=list(geneid,NULL))
    g <- as.vector(h5read(path,paste0(s,'/feature'))[geneid])
    rownames(e) <- g
    e
  }
  
  philist <- lapply(0:maxknotallowed,function(num.knot) {
    if (num.knot==0) {
      # phi <- cbind(1,bs(pseudotime))
      phi <- bs(pseudotime, intercept = TRUE)
    } else {
      knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      # phi <- cbind(1,bs(pseudotime,knots = knots))  
      phi <- bs(pseudotime,knots = knots, intercept = TRUE)
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
  
  ## automatically select knotnum for each gene
  if (is.null(knotnum)){
    bicfunc <- function(num.knot) {
      phi <- philist[[as.character(num.knot)]]
      ll <- sapply(names(sname), function(ss) {
        sexpr <- exprreadfunc(ss,1:length(gn))
        phiss <- phi[sname[[ss]],,drop=F]
        dif2 <- sexpr - sexpr %*% (phiss %*% chol2inv(chol(crossprod(phiss)))) %*% t(phiss)
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
    
    if (is.vector(bic)){
      knotnum <- c(0:maxknot)[which.min(bic)]
    } else {
      knotnum <- c(0:maxknot)[apply(bic,1,which.min)]
    }
    
    names(knotnum) <- gn
  } else {
    knotnum <- knotnum[gn]
  }
  
  sfit <- function(num.knot) {
    #print(paste0('num.knot = ', num.knot))
    gid <- names(which(knotnum==num.knot))
    #print(gid)
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
      # xs <- sapply(row.names(design), function(i) {
      #   tmp <- kronecker(diag(num.knot + 4), design[i, ])
      #   tmp <- tmp[-seq(4, nrow(tmp), 2), ]
      # }, simplify = F)
      xs <- sapply(row.names(design), function(i) {  ## change X
        tmp <- kronecker(diag(num.knot + 4), c(1,1))
        tmp <- tmp[-seq(4, nrow(tmp), 2), ]
        tmp[1,] <- design[i,2]
        tmp
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
    B2 <- Reduce('+',sapply(as,function(s) exprreadfunc(s,match(gid,gn)) %*% phiX[[s]],simplify = F))
    B <- B2 %*% solve(B1)
    
    indfit <- sapply(as,function(s) {
      exprreadfunc(s,match(gid,gn)) %*% (phi[[s]] %*% chol2inv(chol(crossprod(phi[[s]]))))
    },simplify = F)
    
    s2 <- matrix(sapply(as,function(s) {
      tmp <- exprreadfunc(s,match(gid,gn))-indfit[[s]] %*% t(phi[[s]])
      n <- ncol(tmp)
      (rowMeans(tmp*tmp)-(rowMeans(tmp))^2)*(n)/(n-1)
    }),nrow=length(gid),dimnames=list(gid,as))
    alpha <- 1/apply(s2,1,var)*rowMeans(s2)^2+2
    eta <- (alpha-1)*rowMeans(s2)
    
    diffindfit <- sapply(as,function(s) {
      indfit[[s]] - B %*% xs[[s]]
    },simplify = F)
    
    omega <- t(sapply(gid,function(g) {
      m <- sapply(diffindfit,function(i) i[g,,drop=F])
      m <- m-rowMeans(m)
      tcrossprod(m)/(ncol(m)-1) + diag(nrow(m)) * 0.01
    }))
    omega <- omega/rowMeans(s2)
    
    iter <- 0
    gidr <- gid
    all <- matrix(-Inf, nrow=length(gidr),ncol=1,dimnames = list(rownames=gidr))
    #etalist <- alphalist <- omegalist <- Nlist <- Jslist <- list()
    while (iter < EMmaxiter && length(gidr) > 0) {
      
      oinv <- sapply(gidr,function(g) {
        chol2inv(chol(matrix(omega[g,],nrow=nb)))
      },simplify = F)
      
      omegadet <- sapply(gidr,function(g) {
        log(det(matrix(omega[g,],nrow=nb)))/2
      })
      
      ll <- NULL
      for (s in as) {
        sexpr_phibx <- exprreadfunc(s,match(gidr,gn))-B[gidr,] %*% t(phiX[[s]])
        
        L <- rowSums(sexpr_phibx * sexpr_phibx)
        
        Jchol <- sapply(gidr,function(g) {
          chol(phiphi[[s]] + oinv[[g]])
        },simplify = F)
        
        Jsolve <- sapply(gidr,function(g) {
          chol2inv(Jchol[[g]])
        })
        
        K <- tcrossprod(t(phi[[s]]),sexpr_phibx)
        
        L2eKJK <- 2*eta[gidr]+L-colSums(K[rep(1:nb,nb),,drop=FALSE] * K[rep(1:nb,each=nb),,drop=FALSE] * Jsolve)
        A <- log(L2eKJK/2)-digamma(alpha[gidr]+cn[s]/2)
        N <- (2*alpha[gidr]+cn[s])/L2eKJK
        
        JK <- t(rowsum((Jsolve*K[rep(1:nb,nb),,drop=FALSE]),rep(1:nb,each=nb)))
        
        logdv <- -2*colSums(log(sapply(Jchol,as.vector)[seq(1,nb*nb,nb+1),,drop=F]))
        if (is.null(ll)) {
          ll <- alpha[gidr]*log(2*eta[gidr])+lgamma(cn[s]/2+alpha[gidr])-cn[s]*log(pi)/2-lgamma(alpha[gidr])-omegadet+logdv/2-(cn[s]/2+alpha[gidr])*log(L2eKJK)
          B1 <- tcrossprod(N,as.vector(phiXTphiX[[s]]))
          rownames(B1) <- gidr
          B2 <- N * ((exprreadfunc(s,match(gidr,gn)) - t(phi[[s]] %*% t(JK))) %*% phiX[[s]])
          omegalist <- t(Jsolve) + N*JK[,rep(1:nb,nb)] * JK[,rep(1:nb,each=nb)]
          sumA <- A
          sumN <- N
        } else {
          ll <- ll + alpha[gidr]*log(2*eta[gidr])+lgamma(cn[s]/2+alpha[gidr])-cn[s]*log(pi)/2-lgamma(alpha[gidr])-omegadet+logdv/2-(cn[s]/2+alpha[gidr])*log(L2eKJK)
          B1 <- B1 + tcrossprod(N,as.vector(phiXTphiX[[s]]))
          B2 <- B2 + N * ((exprreadfunc(s,match(gidr,gn)) - t(phi[[s]] %*% t(JK))) %*% phiX[[s]])
          omegalist <- omegalist + t(Jsolve) + N*JK[,rep(1:nb,nb)] * JK[,rep(1:nb,each=nb)]
          sumA <- sumA + A
          sumN <- sumN + N
        }
      }
      
      np <- nrow(xs[[1]])
      B[gidr,] <- t(sapply(gidr, function(g){  ## each column is a gene's all betas
        chol2inv(chol(matrix(B1[g,],nrow=np))) %*% B2[g,]
      }))
      
      omega[gidr,] <- omegalist/length(as)
      
      rN <- sumN/length(as)
      rA <- sumA/length(as)
      eta[gidr] <- sapply(gidr,function(g) {
        meanN <- rN[g]
        meanA <- rA[g]
        # uniroot(function(eta) {digamma(eta * meanN)-log(eta)+meanA},c(1e-10,1e10))$root
        optim(eta[g],fn = function(eta) {(digamma(eta * meanN)-log(eta)+meanA)^2},gr = function(eta) {2*(digamma(eta * meanN)-log(eta)+meanA)*(trigamma(eta*meanN)*meanN-1/eta)},lower = 1e-10,method = 'L-BFGS-B')$par
      })
      alpha[gidr] <- eta[gidr] * rN
      
      iter <- iter + 1
      
      llv <- all[,ncol(all)]
      llv[gidr] <- ll
      all <- cbind(all,llv)
      gidr <- names(which(all[,ncol(all)] - all[,ncol(all)-1] > EMitercutoff))
      #etalist[[iter]] <- eta
      #alphalist[[iter]] <-alpha
      #omegalist[[iter]] <-omega
      #Nlist[[iter]] <- N
      #Jslist[[iter]] <- Jsolve
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

