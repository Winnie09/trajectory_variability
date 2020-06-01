##### implement model
diffpt <- function(expr, cellanno, pseudotime, design, num.knot=3, maxiter=10, verbose=F) {
  library(Matrix)
  library(mvtnorm)
  library(parallel)
  library(splines)
  ## expr: gene by cell matrix, entires and log-transformed expression
  ## Pseudotime: a  vecotor of ordered cell names
  ## design: sample by feature design matrix. rownames are sample names. first column is the group partition, currently only one column.
  ## cellanno: dataframe, first column is cell names, second column is sample names.
  design <- as.matrix(design)
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
  phi_xs <- sapply(names(xs),function(i) phi[sname[[i]],,drop=F] %*% t(xs[[i]]))
  phi_xs_rbind <- do.call(rbind,phi_xs)
  phi_xs_rbind <- phi_xs_rbind[colnames(expr),]
  beta <- expr %*% (phi_xs_rbind %*% chol2inv(chol(crossprod(phi_xs_rbind))))
  
  ## initialize tau
  indfit <- sapply(names(xs), function(ss){
    print(ss)
    if (length(sname[[ss]]) > 1){
      expr[,sname[[ss]],drop=F] %*% (phi[sname[[ss]],] %*% chol2inv(chol(crossprod(phi[sname[[ss]],]))))
    } else {
      expr[,sname[[ss]],drop=F] %*% t(phi[sname[[ss]],] %*% chol2inv(chol(crossprod(phi[sname[[ss]],]))))
    }
  },simplify = F)
  indfitdiff <- sapply(names(xs), function(ss){
    indfit[[ss]] - beta %*% xs[[ss]]
  },simplify = F)
  tau <- sapply(row.names(indfitdiff[[1]]),function(rid) {
    m <- sapply(names(xs),function(n) indfitdiff[[n]][rid,])
    m <- m-rowMeans(m)
    tcrossprod(m) / ncol(m) + diag(nrow(m)) * 0.01
  })
  
  ## initialize gamma
  gamma <- sapply(names(xs), function(ss){
    diff <- expr[,sname[[ss]]] - indfit[[ss]] %*% t(phi[sname[[ss]],])
    rowMeans(diff * diff) - rowMeans(diff)^2
  })
  gamma[gamma < 0.01] <- 0.01
  
  for (dump in seq(1,maxiter)){
    E_phi_u_e <- phi_xs_beta <- phi_tau_phi <- list()
    M <- matrix(0,nrow=nrow(expr),ncol=nrow(design),dimnames=list(row.names(expr),names(sname)))
    for (ss in names(sname)) {
      phi_tau_phi[[ss]] <- t(tau) %*% phicrossprod[,sname[[ss]]]
      phi_xs_beta[[ss]] <- beta %*% t(phi_xs[[ss]])
      sigma1s <- 1/(1/phi_tau_phi[[ss]] + 1/gamma[,ss])
      E_phi_u_e[[ss]] <- sigma1s * (phi_xs_beta[[ss]]/phi_tau_phi[[ss]] + expr[,sname[[ss]]]/gamma[,ss]) - phi_xs_beta[[ss]]
      M[,ss] <- (rowSums((expr[,sname[[ss]]] - phi_xs_beta[[ss]] - E_phi_u_e[[ss]])^2) + rowSums(sigma1s))/ncol(sigma1s)
    }

    tauinv <- lapply(colnames(tau),function(i) {
      chol2inv(chol(matrix(tau[,i],nrow=ncol(phi))))
    })
    names(tauinv) <- colnames(tau)
    E_u_con_theta <- lapply(names(sname),function(ss) {
      m <- matrix(rowSums(phicrossprod[,sname[[ss]]]),nrow=ncol(phi))
      tmp <- lapply(colnames(tau),function(i) {
        chol2inv(chol(m/gamma[i,ss] + tauinv[[i]]))
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
    newgamma <- M
    
    ### tau
    tmp <- lapply(names(xs), function(ss) {
      sapply(E_u_con_theta[[ss]],as.vector) + apply(k_s[[ss]],1,tcrossprod)
    })
    newtau <- Reduce('+', tmp)/length(xs)
    
    ### log L
    ll <- sapply(names(xs), function(ss) {
      tmp <- 2 * (phi_tau_phi[[ss]] + gamma[,ss])
      rowSums(-log(pi * tmp)/2) - rowSums((expr[,sname[[ss]]] - phi_xs_beta[[ss]])^2 / tmp)
    })
    ll <- rowSums(ll)
    if (verbose) print(mean(ll))
    beta <- newbeta
    gamma <- newgamma
    tau <- newtau
  }
  return(list(beta = beta, gamma = gamma, tau = tau, logL = ll))
}

