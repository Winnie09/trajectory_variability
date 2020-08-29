diffpt <- function(expr, cellanno, pseudotime, design, num.knot=3, maxiter=10, verbose=F) {
  library(Matrix)
  library(mvtnorm)
  library(parallel)
  library(splines)
  ## expr: gene by cell matrix, entires are log-transformed expression
  ## pseudotime: a  vecotor of ordered cell names
  ## design: sample by feature design matrix. rownames are sample names. first column is the group partition, currently only one column.
  ## cellanno: dataframe, first column is cell names, second column is sample names.
  psn <- 1:length(pseudotime)
  names(psn) <- pseudotime
  psn <- psn[colnames(expr)]
  cellanno <- cellanno[match(colnames(expr),cellanno[,1]),]
  knots = seq(min(psn),max(psn),length.out=num.knot+2)[2:(num.knot+1)]
  phi <- cbind(1,bs(psn,knots = knots))
  ## initialize data
  num.sample <- length(unique(cellanno$sample))  
  num.coef <- num.knot + 4
  
  phicrossprod <- apply(phi,1,tcrossprod)
  sname <- sapply(row.names(design),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  xs <- sapply(row.names(design),function(i) {
    kronecker(diag(num.knot + 4),design[i,])
  },simplify = F)
  
  ## initialize tau
  phi_xs <- sapply(names(xs),function(i) phi[sname[[i]],] %*% t(xs[[i]]))
  phi_xs_rbind <- do.call(rbind,phi_xs)
  beta <- expr[,row.names(phi_xs_rbind)] %*% (phi_xs_rbind %*% chol2inv(chol(t(phi_xs_rbind) %*% phi_xs_rbind)))
  
  ## initialize tau
  indfit <- sapply(names(xs), function(ss){
    expr[,sname[[ss]]] %*% (phi[sname[[ss]],] %*% chol2inv(chol(t(phi[sname[[ss]],]) %*% phi[sname[[ss]],])))
  },simplify = F)
  indfitdiff <- sapply(names(xs), function(ss){
    indfit[[ss]] - beta %*% xs[[ss]]
  },simplify = F)
  tau <- sapply(row.names(indfitdiff[[1]]),function(rid) {
    m <- sapply(names(xs),function(n) indfitdiff[[n]][rid,])
    m <- m-rowMeans(m)
    m %*% t(m) / ncol(m) + diag(nrow(m)) * 0.01
  })
  
  ## initialize gamma
  gamma <- sapply(names(xs), function(ss){
    diff <- expr[,sname[[ss]]] - indfit[[ss]] %*% t(phi[sname[[ss]],])
    rowMeans(diff * diff) - rowMeans(diff)^2
  })
  gamma[gamma < 0.01] <- 0.01
  
  for (dump in seq(1,maxiter)){
    E_phi_u_e <- M <- phi_xs_beta <- phi_tau_phi <- list()
    for (ss in names(sname)) {
      phi_tau_phi[[ss]] <- t(tau) %*% phicrossprod[,sname[[ss]]]
      phi_xs_beta[[ss]] <- beta %*% t(phi_xs[[ss]])
      sigma1s <- 1/(1/phi_tau_phi[[ss]] + 1/gamma[,ss])
      E_phi_u_e[[ss]] <- sigma1s * (phi_xs_beta[[ss]]/phi_tau_phi[[ss]] + expr[,sname[[ss]]]/gamma[,ss]) - phi_xs_beta[[ss]]
      M[[ss]] <- (expr[,sname[[ss]]] - phi_xs_beta[[ss]] - E_phi_u_e[[ss]])^2 + sigma1s
    }
    
    E_u_con_theta <- lapply(names(sname),function(ss) {
      m <- matrix(rowSums(phicrossprod[,sname[[ss]]]),nrow=ncol(phi))
      n <- ncol(phi)
      tmp <- lapply(colnames(tau),function(i) {
        chol2inv(chol(m/gamma[i,ss] + chol2inv(chol(matrix(tau[,i],nrow=n)))))
      })
      names(tmp) <- colnames(tau)
      tmp
    })
    names(E_u_con_theta) <- names(sname)
    
    k_s <- lapply(names(sname),function(ss){
      tmp <- (expr[,sname[[ss]]] - phi_xs_beta[[ss]]) %*% phi[sname[[ss]],]
      t(sapply(row.names(tmp),function(i) tmp[i,] %*% E_u_con_theta[[ss]][[i]])) / gamma[,ss]
    })
    names(k_s) <- names(sname)
    
    tmp2 <- (1/gamma) %*% t(sapply(names(xs), function(ss){
      as.vector(t(phi_xs[[ss]]) %*% phi_xs[[ss]])
    }))
    
    tmp3 <- lapply(names(xs), function(ss){
      (expr[,sname[[ss]]] - E_phi_u_e[[ss]]) %*% phi_xs[[ss]] / gamma[,ss]
    })
    tmp3 <- Reduce('+', tmp3)     
    newbeta <- t(sapply(row.names(tmp2),function(i) chol2inv(chol(matrix(tmp2[i,],nrow=ncol(beta)))) %*% tmp3[i,]))
    
    ### gamma
    newgamma <- sapply(M, rowMeans)
    
    ### tau
    tmp <- lapply(names(xs), function(ss) {
      sapply(names(E_u_con_theta[[ss]]),function(i) {
        E_u_con_theta[[ss]][[i]] + k_s[[ss]][i,] %*% t(k_s[[ss]][i,])  
      })
    })
    newtau <- Reduce('+', tmp)/length(xs)
    
    ### log L
    ll <- sapply(names(xs), function(ss) {
      tmp <- log(1/sqrt(2 * pi * (phi_tau_phi[[ss]] + gamma[,ss]))) - ((expr[,sname[[ss]]] - phi_xs_beta[[ss]])^2 / (2 * (phi_tau_phi[[ss]] + gamma[,ss])))
    })
    ll <- rowSums(sapply(ll,rowSums))
    if (verbose) print(mean(ll))
    beta <- newbeta
    gamma <- newgamma
    tau <- newtau
  }
  return(list(beta = beta, gamma = gamma, tau = tau, logL = ll))
}
