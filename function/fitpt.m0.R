# ori.design = design
# test.pattern = 'overall'
# test.position = 'all'
# maxknotallowed=10; EMmaxiter=1000; EMitercutoff=0.01; verbose=F; ncores=1; model = 3
# test.pattern = 'overall'
fitpt.m0 <- function(expr, cellanno, pseudotime, design, EMmaxiter=1000, EMitercutoff=0.01, verbose=F, ncores=1) {
  # set.seed(12345)
  suppressMessages(library(Matrix))
  suppressMessages(library(parallel))
  suppressMessages(library(splines))
  suppressMessages(library(matrixcalc))
  ## This function is for model 0 (beta0 for intercept base only) fitting. 
  ## expr: gene by cell matrix, entires and log-transformed expression
  ## pseudotime: a numeric vecotor of pseudotime, the names are ordered cell names
  ## design: design: sample by feature design matrix. rownames are sample names. first column is 1, second column is the group partition, currently only one variable.
  ## cellanno: dataframe, first column is cell names, second column is sample names.
  ## test.pattern: c('slope', 'intercept', 'overall).
  pseudotime <- pseudotime[colnames(expr)]
  cellanno <- cellanno[match(colnames(expr),cellanno[,1]),]
  sname <- sapply(row.names(design),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  design = as.matrix(design)
  
  phi <- matrix(1, nrow = length(pseudotime), ncol = 1, dimnames = list(names(pseudotime), NULL))
  sexpr = expr
  phicrossprod <- apply(phi,1,tcrossprod)
  phicrossprod <- sapply(names(sname),function(ss) phicrossprod[sname[[ss]]],simplify = F)
  phi <- sapply(names(sname),function(ss) phi[sname[[ss]],,drop=F],simplify = F) 
  
  xs <- sapply(row.names(design), function(i){
    as.matrix(1, nrow = 1, ncol = 1)
  }, simplify = FALSE)
  phi <- sapply(phi, function(i){
    i[,1,drop=F]
  }, simplify = FALSE)

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
  }),nrow=nrow(expr),dimnames=list(rownames(expr),as))
  alpha <- 1/apply(s2,1,var)*rowMeans(s2)^2+2
  eta <- (alpha-1)*rowMeans(s2)
  
  diffindfit <- sapply(as,function(s) {
    indfit[[s]] - B %*% xs[[s]]
  },simplify = F)
  
  omega <- sapply(rownames(sexpr),function(g) {
    m <- sapply(diffindfit,function(i) i[g,,drop=F])
    var(m) + 0.01
  })
  
  iter <- 0
  EMitercutoff <- 0
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
        1/(phiphi[[s]] + 1/omega[g])
      })
    })
    
    K <- sapply(as,function(s) {
      tcrossprod(t(phi[[s]]),sexpr_phibx[[s]])
    },simplify = F)
    
    L2eKJK <- sapply(as,function(s) {
      2*eta[gidr] + L[[s]] - colSums(K[[s]] * K[[s]] * Jsolve[,s])
    },simplify = F)
    
    A <- matrix(sapply(as,function(s) {
      log(L2eKJK[[s]]/2)-digamma(alpha[gidr]+cn[s]/2)
    }),nrow=length(gidr),dimnames=list(gidr,as))
    
    N <- matrix(sapply(as,function(s) {
      (2*alpha[gidr]+cn[s])/L2eKJK[[s]]
    }),nrow=length(gidr),dimnames=list(gidr,as))
    
    ll <- rowSums(matrix(sapply(as,function(s) {
      dv <- sapply(gidr,function(g) {
        omega[g]/Jsolve[g,s]
      })
      alpha[gidr]*log(2*eta[gidr])+lgamma(cn[s]/2+alpha[gidr])-cn[s]*log(pi)/2-lgamma(alpha[gidr])-log(dv)/2-(cn[s]/2+alpha[gidr])*log(L2eKJK[[s]])
    }),nrow=length(gidr),dimnames = list(rownames=gidr,colnames=as)))
    
    JK <- sapply(as,function(s) {
      Jsolve[,s] * K[[s]]
    },simplify = F)
    
    ## -------------->
    B1 <- Reduce('+', lapply(as, function(s){
      tcrossprod(N[,s],as.vector(phiXTphiX[[s]]))
    }))
    rownames(B1) <- gidr
    B2 <- Reduce('+', lapply(as, function(s){
      N[,s] * t(t(phiX[[s]]) %*%  
                  (t(sexpr[ ,cellanno[,2]==s, drop=F][gidr,]) - as.vector(phi[[s]] %*% JK[[s]]))
      )
    }))
    
    np <- nrow(xs[[1]])
    B[gidr, ] <- t(sapply(gidr, function(g){  ## each column is a gene's all betas
      B2[g,]/B1[g,]
    }))
    
    ## M -step: 
    omega[gidr] <- rowMeans(sapply(as,function(s) {
      Jsolve[,s] + (N[,s]*JK[[s]] * JK[[s]])[1,]
    }))
    
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
  
  allres  <- list(beta = B, alpha = alpha, eta = eta, omega = omega, logL = all)
  
  para <- list()
  
  for (j in row.names(allres[[1]])) {
    para[[j]] <- list(beta=allres[[1]][j,],
                      alpha=allres[[2]][j],
                      eta=allres[[3]][j],
                      omega=allres[[4]][j],
                      ll=allres[[5]][j,ncol(allres[[5]])])
  }
  
  return(list(parameter=para))
}



