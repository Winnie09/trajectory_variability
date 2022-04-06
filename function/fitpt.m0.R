# ori.design = design
# test.pattern = 'overall'
# test.position = 'all'
# maxknotallowed=10; EMmaxiter=1000; EMitercutoff=0.01; verbose=F; ncores=1; model = 3
# test.pattern = 'overall'
fitpt.m0 <- function(expr, cellanno, pseudotime, design, EMmaxiter=100, EMitercutoff=0.05, verbose=F) {
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
  sname <- sapply(unique(cellanno[,2]),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  as <- names(sname)
  cn <- sapply(sname,length)
  
  ## initialize
  B <- rowMeans(expr)
  
  indfit <- sapply(as,function(s) {
    rowMeans(expr[, cellanno[,2]==s, drop=F])
  },simplify = F)
  
  s2 <- matrix(sapply(as,function(s) {
    tmp <- expr[, cellanno[,2]==s, drop=F]-indfit[[s]]
    n <- ncol(tmp)
    (rowMeans(tmp*tmp)-(rowMeans(tmp))^2)*(n)/(n-1)
  }),nrow=nrow(expr),dimnames=list(rownames(expr),as))
  alpha <- 1/apply(s2,1,var)*rowMeans(s2)^2+2
  eta <- (alpha-1)*rowMeans(s2)
  
  diffindfit <- sapply(as,function(s) {
    indfit[[s]] - B
  })
  
  s2[s2 < 0.001] <- 0.001
  omega <- rowMeans(diffindfit^2/s2)
  
  iter <- 0
  gidr <- rownames(expr)
  oldpara <- list(beta = B, alpha = alpha, eta = eta, omega = omega)
  
  while (iter < EMmaxiter && length(gidr) > 0) {
    expr_phibx <- sapply(as,function(s) {
      expr[, cellanno[,2]==s, drop=F][gidr,,drop=F]-B[gidr]
    },simplify = F)
    
    L <- sapply(as,function(s) {
      rowSums(expr_phibx[[s]] * expr_phibx[[s]])
    },simplify = F)
    
    Jsolve <- matrix(sapply(as,function(s) {
      1/(cn[s] + 1/omega[gidr])
    }), nrow = length(gidr), dimnames = list(gidr, as))
    
    K <- sapply(as,function(s) {
      rowSums(expr_phibx[[s]])
    },simplify = F)
    
    JK <- sapply(as,function(s) {
      Jsolve[,s,drop=F] * K[[s]]  ## debug here !!
    },simplify = F)
    
    L2eKJK <- sapply(as,function(s) {
      2*eta[gidr] + L[[s]] - K[[s]] * JK[[s]]
    },simplify = F)
    
    A <- matrix(sapply(as,function(s) {
      log(L2eKJK[[s]]/2)-digamma(alpha[gidr]+cn[s]/2)
    }),nrow=length(gidr),dimnames=list(gidr,as))
    
    N <- matrix(sapply(as,function(s) {
      (2*alpha[gidr]+cn[s])/L2eKJK[[s]]
    }),nrow=length(gidr),dimnames=list(gidr,as))
    
    B1 <- rowSums(matrix(sapply(as, function(s){
      N[,s] * cn[s]
    }), ncol = length(as), dimnames = list(gidr, as)))  ## length(gidr) * length(as) debug here !!!
    B2 <- rowSums(matrix(sapply(as, function(s){   ## debug here !!!
      N[,s,drop=F] * colSums(t(expr[ ,cellanno[,2]==s, drop=F][gidr,,drop=F]) - rep(JK[[s]],each=cn[s]))
    }), ncol = length(as), dimnames = list(gidr, as)))  ## vector of length(gidr) debug here !!!
    B[gidr] <- B2/B1
    
    ## M -step:
    omega[gidr] <- rowMeans(matrix(sapply(as,function(s) {
      Jsolve[,s,drop=F] + N[,s,drop=F]*JK[[s]] * JK[[s]]  ## debug here !!!
    }), ncol=length(as), dimnames = list(gidr, as)))
    
    eta[gidr] <- sapply(gidr,function(g) {
      meanN <- mean(N[g,])
      meanA <- mean(A[g,])
      uniroot(function(eta) {digamma(eta * meanN)-log(eta)+meanA},c(1e-10,1e10))$root
    })
    alpha[gidr] <- eta[gidr] * rowMeans(N)
    para <- list(beta = B, alpha = alpha, eta = eta, omega = omega)
    paradiff <- sapply(names(para),function(s) {
      if (is.vector(para[[s]])) {
        abs(para[[s]]-oldpara[[s]])/abs(oldpara[[s]])
      } else {
        rowSums(abs(para[[s]]-oldpara[[s]]))/rowSums(abs(oldpara[[s]]))
      }
    })
    if (is.vector(paradiff)) paradiff <- matrix(paradiff,nrow=1,dimnames=list(gid,NULL))    
    paradiff <- apply(paradiff,1,max)
    gidr <- names(which(paradiff > EMitercutoff))
    oldpara <- para
    iter <- iter + 1
    
    rm(list = c('L','Jsolve', 'K'))
  }
  
  ll <- rowSums(matrix(sapply(as,function(s) {
    expr_phibx <- sapply(as,function(s) {
      expr[, cellanno[,2]==s, drop=F]-B
    },simplify = F)
    
    L <- sapply(as,function(s) {
      rowSums(expr_phibx[[s]] * expr_phibx[[s]])
    },simplify = F)
    
    Jsolve <- matrix(sapply(as,function(s) {
      1/(cn[s] + 1/omega)
    }), nrow = nrow(expr), dimnames = list(rownames(expr), as))
    
    K <- sapply(as,function(s) {
      rowSums(expr_phibx[[s]])
    },simplify = F)
    
    JK <- sapply(as,function(s) {
      Jsolve[,s,drop=F] * K[[s]]  ## debug here !!
    },simplify = F)
    
    L2eKJK <- sapply(as,function(s) {
      2*eta + L[[s]] - K[[s]] * JK[[s]]
    },simplify = F)
    dv <- omega/Jsolve[,s,drop=F]
    alpha*log(2*eta)+lgamma(cn[s]/2+alpha)-cn[s]*log(pi)/2-lgamma(alpha)-log(dv)/2-(cn[s]/2+alpha)*log(L2eKJK[[s]])
  }),nrow=nrow(expr),dimnames=list(rownames(expr),NULL)))
  
  allres  <- list(beta = B, alpha = alpha, eta = eta, omega = omega, ll = ll)
  
  para <- list()
  
  for (j in rownames(expr)) {
    para[[j]] <- list(beta=allres[[1]][j],
                      alpha=allres[[2]][j],
                      eta=allres[[3]][j],
                      omega=allres[[4]][j],
                      ll=allres[[5]][j])
  }
  
  return(list(parameter=para))
}





