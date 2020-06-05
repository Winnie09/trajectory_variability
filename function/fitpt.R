fitpt <- function(expr, cellanno, pseudotime, design, knotnum=NULL, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores()) {
  library(Matrix)
  library(parallel)
  library(splines)
  library(matrixcalc)
  ## expr: gene by cell matrix, entires and log-transformed expression
  ## Pseudotime: a  vecotor of ordered cell names
  ## design: design: sample by feature design matrix. rownames are sample names. first column is 1, second column is the group partition, currently only one variable.
  ## cellanno: dataframe, first column is cell names, second column is sample names.
  pseudotime <- pseudotime[colnames(expr)]
  cellanno <- cellanno[match(colnames(expr),cellanno[,1]),]
  sname <- sapply(row.names(design),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  design = as.matrix(design)
  if (is.null(knotnum)) {
    id <- split(sample(colnames(expr)),ceiling(seq_along(colnames(expr))/(ncol(expr)/5)))
    maxknot <- 0
    testpos <- 1
    while (testpos == 1 & maxknot < 50) {
      maxknot <- maxknot + 1
      knots = seq(min(pseudotime),max(pseudotime),length.out=maxknot+2)[2:(maxknot+1)]
      phi <- cbind(1,bs(pseudotime,knots = knots))
      testpos <- mean(sapply(1:5,function(cvid) {
        sapply(names(sname), function(ss){
          traincell <- intersect(sname[[ss]],unlist(id[-cvid]))
          is.positive.definite(crossprod(phi[traincell,]))
        })
      }))
    }
    print('selecting the optimal knots ...')  
    diff <- mclapply(0:maxknot,function(num.knot) {
      if (num.knot==0) {
        phi <- cbind(1,bs(pseudotime))
      } else {
        knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
        phi <- cbind(1,bs(pseudotime,knots = knots))  
      }
      rowSums(sapply(1:5,function(cvid) {
        rowSums(sapply(names(sname), function(ss){
          traincell <- intersect(sname[[ss]],unlist(id[-cvid]))
          testcell <- intersect(sname[[ss]],id[[cvid]])
          fit <- expr[,traincell] %*% (phi[traincell,] %*% chol2inv(chol(crossprod(phi[traincell,])))) %*% t(phi[testcell,])
          rowSums((fit-expr[,testcell])^2)
        }))
      }))
    },mc.cores=ncores)
    diff <- do.call(cbind,diff)
    knotnum <- c(0:50)[apply(diff,1,which.min)]
    names(knotnum) <- row.names(expr)  
  }
 
  sfit <- function(num.knot) {
    gid <- names(which(knotnum==num.knot))
    if (num.knot==0) {
      phi <- cbind(1,bs(pseudotime))
    } else {
      knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      phi <- cbind(1,bs(pseudotime,knots = knots))  
    }
    phicrossprod <- apply(phi,1,tcrossprod)
   
    xs <- sapply(row.names(design),function(i) {
      kronecker(diag(num.knot + 4),design[i,])
    },simplify = F)
   
    ## initialize tau
    phi_xs <- sapply(names(xs),function(i) phi[sname[[i]],] %*% t(xs[[i]]))
    phi_xs_rbind <- do.call(rbind,phi_xs)
    phi_xs_rbind <- phi_xs_rbind[colnames(expr),]
    beta <- expr[gid,,drop=F] %*% (phi_xs_rbind %*% chol2inv(chol(crossprod(phi_xs_rbind))))
   
    ## initialize tau
    indfit <- sapply(names(xs), function(ss){
      expr[gid,sname[[ss]],drop=F] %*% (phi[sname[[ss]],] %*% chol2inv(chol(crossprod(phi[sname[[ss]],]))))
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
      diff <- expr[gid,sname[[ss]],drop=F] - indfit[[ss]] %*% t(phi[sname[[ss]],])
      rowMeans(diff * diff) - rowMeans(diff)^2
    })
    if (is.vector(gamma)) {
      gamma <- t(gamma)
      colnames(gamma) <- names(xs)
      row.names(gamma) <- gid
    }
    gamma[gamma < 0.01] <- 0.01
   
    iter <- 0
    ll <- rep(-Inf,length(gid))
    names(ll) <- gid
    gidr <- gid
    while (iter < EMmaxiter && length(gidr) > 0) {
      llold <- ll
      E_phi_u_e <- phi_xs_beta <- phi_tau_phi <- list()
      M <- matrix(0,nrow=length(gidr),ncol=nrow(design),dimnames=list(gidr,names(sname)))
      for (ss in names(sname)) {
        phi_tau_phi[[ss]] <- t(tau[,gidr,drop=F]) %*% phicrossprod[,sname[[ss]]]
        phi_xs_beta[[ss]] <- beta[gidr,,drop=F] %*% t(phi_xs[[ss]])
        sigma1s <- 1/(1/phi_tau_phi[[ss]] + 1/gamma[gidr,ss])
        E_phi_u_e[[ss]] <- sigma1s * (phi_xs_beta[[ss]]/phi_tau_phi[[ss]] + expr[gidr,sname[[ss]]]/gamma[gidr,ss]) - phi_xs_beta[[ss]]
        M[,ss] <- (rowSums((expr[gidr,sname[[ss]],drop=F] - phi_xs_beta[[ss]] - E_phi_u_e[[ss]])^2) + rowSums(sigma1s))/ncol(sigma1s)
      }
     
      tauinv <- lapply(gidr,function(i) {
        chol2inv(chol(matrix(tau[,i],nrow=ncol(phi))))
      })
      names(tauinv) <- gidr
      E_u_con_theta <- lapply(names(sname),function(ss) {
        m <- matrix(rowSums(phicrossprod[,sname[[ss]]]),nrow=ncol(phi))
        tmp <- lapply(gidr,function(i) {
          chol2inv(chol(m/gamma[i,ss] + tauinv[[i]]))
        })
        names(tmp) <- gidr
        tmp
      })
      names(E_u_con_theta) <- names(sname)
     
      k_s <- lapply(names(sname),function(ss){
        tmp <- (expr[gidr,sname[[ss]],drop=F] - phi_xs_beta[[ss]]) %*% phi[sname[[ss]],]
        t(sapply(row.names(tmp),function(i) tmp[i,] %*% E_u_con_theta[[ss]][[i]])) / gamma[gidr,ss]
      })
      names(k_s) <- names(sname)
     
      tmp2 <- (1/gamma[gidr,,drop=F]) %*% t(sapply(names(xs), function(ss){
        as.vector(t(phi_xs[[ss]]) %*% phi_xs[[ss]])
      }))
     
      tmp3 <- lapply(names(xs), function(ss){
        (expr[gidr,sname[[ss]],drop=F] - E_phi_u_e[[ss]]) %*% phi_xs[[ss]] / gamma[gidr,ss]
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
      tmpll <- sapply(names(xs), function(ss) {
        tmp <- 2 * (phi_tau_phi[[ss]] + gamma[gidr,ss])
        rowSums(-log(pi * tmp)/2) - rowSums((expr[gidr,sname[[ss]],drop=F] - phi_xs_beta[[ss]])^2 / tmp)
      })
      if (is.vector(tmpll)) {
        tmpll <- t(tmpll)
        row.names(tmpll) <- gidr
      }
      ll[gidr] <- rowSums(tmpll)
      beta[gidr,] <- newbeta
      gamma[gidr,] <- newgamma
      tau[,gidr] <- newtau
      iter <- iter + 1
      gidr <- names(which(ll-llold > EMitercutoff))
    }
    return(list(beta = beta, gamma = gamma, tau = tau, logL = ll))
  }
  
  allres <- mclapply(unique(knotnum),sfit,mc.cores=ncores)
 
  para <- list()
  for (i in 1:length(allres)) {
    for (j in row.names(allres[[i]][[1]])) {
      para[[j]] <- list(beta=allres[[i]][[1]][j,],gamma=allres[[i]][[2]][j,],tau=matrix(allres[[i]][[3]][,j],nrow=sqrt(length(allres[[i]][[3]][,j]))),ll=allres[[i]][[4]][j])
    }
  }
  list(parameter=para,knotnum=knotnum)
}

