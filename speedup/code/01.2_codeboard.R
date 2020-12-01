expr=expr; cellanno=cellanno; pseudotime=pseudotime; design=design; EMmaxiter=EMmaxiter; EMitercutoff=EMitercutoff; verbose=verbose; ncores=1; test.pattern = test.pattern; test.position = test.position
ori.design = design; maxknotallowed=10; EMmaxiter=100; EMitercutoff=1; verbose=F; ncores=detectCores()

fitpt <- function(expr, cellanno, pseudotime, design, ori.design = design, test.pattern = 'overall', test.position = 'all',  maxknotallowed=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores()) {
  print('Running fitpt ...')
  suppressMessages(library(Matrix))
  suppressMessages(library(parallel))
  suppressMessages(library(splines))
  suppressMessages(library(matrixcalc))
  ## expr: gene by cell matrix, entires and log-transformed expression
  ## pseudotime: a numeric vecotor of pseudotime, the names are ordered cell names
  ## design: design: sample by feature design matrix. rownames are sample names. first column is 1, second column is the group partition, currently only one variable.
  ## cellanno: dataframe, first column is cell names, second column is sample names.
  ## test.pattern: c('slope', 'intercept', 'overall).
  pseudotime <- pseudotime[colnames(expr)]
  cellanno <- cellanno[match(colnames(expr),cellanno[,1]),]
  sname <- sapply(row.names(design),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  design = as.matrix(design)

  ### chunk 1 ----------------
   system.time({
  philist <- lapply(0:maxknotallowed,function(num.knot) {
    if (num.knot==0) {
      phi <- cbind(1,bs(pseudotime))
    } else {
      knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      phi <- cbind(1,bs(pseudotime,knots = knots))  
    }
  })
  names(philist) <- as.character(0:maxknotallowed)
   })
  #     user  system elapsed 
  # 0.310   0.129   0.439 
   ### chunk 2 ----------------
  system.time({
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
})
    
  #  user  system elapsed 
  # 1.155   0.666   1.823   
      
  ### chunk 3 ----------------
  system.time({
  bicfunc <- function(num.knot) {
    phi <- philist[[as.character(num.knot)]]
    ll <- sapply(names(sname), function(ss) {
      phiss <- phi[sname[[ss]],,drop=F]
      dif2 <- sexpr[[ss]] - sexpr[[ss]] %*% (phiss %*% chol2inv(chol(crossprod(phiss)))) %*% t(phiss)
      dif2 <- rowSums(dif2 * dif2)
      s2 <- dif2/(length(sname[[ss]])-ncol(phi))
      log(2*pi*s2)*nrow(phiss) + dif2/s2
    })
    rowSums(ll,na.rm=T) + log(nrow(phi))*((ncol(phi)+1)*rowSums(!is.na(ll)))
  }
  
  if (ncores!=1) {
    ## chunk 3.1 ------------
    system.time({
    bic <- mclapply(0:maxknot,bicfunc,mc.cores=ncores)
    bic <- do.call(cbind,bic)
    })
  } else {
    ## chunk 3.2 ------------
    system.time({
    bic <- sapply(0:maxknot,bicfunc)
    })
  }
  
  knotnum <- c(0:maxknot)[apply(bic,1,which.min)]
  names(knotnum) <- row.names(expr)
  })
  #  user  system elapsed 
  # 105.970  48.296  20.343  --> 13.602 
  
  ## chunk 4 ------------------
  system.time({
  sfit <- function(num.knot) {
    ## chunk 4.1 ---------------
    system.time({
    gid <- names(which(knotnum==num.knot))
    sexpr <- expr[gid,,drop=F]
    phi <- philist[[as.character(num.knot)]]
    phicrossprod <- apply(phi,1,tcrossprod)
    phicrossprod <- sapply(names(sname),function(ss) phicrossprod[,sname[[ss]]],simplify = F)
    phi <- sapply(names(sname),function(ss) phi[sname[[ss]],],simplify = F)
    })
    # user  system elapsed 
    # 0.853   0.609   1.462 
    
    # --------------
    # change here >>
    # --------------
    ## chunk 4.2 ------------------
    system.time({
    design = design[rownames(ori.design), ,drop=F]
    xs <- sapply(row.names(ori.design),function(i) {  
      kronecker(diag(num.knot + 4), ori.design[i,])
    },simplify = F)
    xs_test <- sapply(row.names(design),function(i) {  
      kronecker(diag(num.knot + 4),design[i,])
    },simplify = F)
    if (test.pattern == 'slope'){
      if (is.na(test.position) | test.position == 'all'){
        for (id in seq(1, length(xs))){
          xs[[id]][, seq(2, ncol(xs[[1]]))] <- xs_test[[id]][,seq(2, ncol(xs[[1]]))]
        }
      } else if (test.position == 'start'){  ## check
        for (id in seq(1, length(xs))){
          xs[[id]][,seq(2, ceiling(ncol(xs[[1]])/3))] <- xs_test[[id]][,seq(2, ceiling(ncol(xs[[1]])/3))]
        }
      } else if (test.position == 'middle'){## check
        for (id in seq(1, length(xs))){
          xs[[id]][,seq(ceiling(ncol(xs[[1]])/3+1), ceiling(ncol(xs[[1]])*2/3))] <- xs_test[[id]][,seq(ceiling(ncol(xs[[1]])/3+1), ceiling(ncol(xs[[1]])*2/3))]
        }
      } else if (test.position == 'end'){ ## check
        for (id in seq(1, length(xs))){
          xs[[id]][, seq(ceiling(ncol(xs[[1]])*2/3+1), ceiling(ncol(xs[[1]])))] <- xs_test[[id]][,seq(ceiling(ncol(xs[[1]])*2/3+1), ceiling(ncol(xs[[1]])))]
        }
      }
    } else if (test.pattern == 'intercept'){
        for (id in seq(1, length(xs))){
          xs[[id]][,1] <- xs_test[[id]][,1]
        }
    } else if (test.pattern == 'overall'){
        if (is.na(test.position) | test.position == 'all'){
          xs <- xs_test
        } else if (test.position == 'start'){
            for (id in seq(1, length(xs))){
              xs[[id]][, c(1, seq(2, ceiling(ncol(xs[[1]])/3)))] <- xs_test[[id]][ , c(1, seq(2, ceiling(ncol(xs[[1]])/3)))]
            }
        } else if (test.position == 'middle'){
            for (id in seq(1, length(xs))){
              xs[[id]][, c(1, seq(ceiling(ncol(xs[[1]])/3+1), ceiling(ncol(xs[[1]])*2/3)))] <- xs_test[[id]][, c(1, seq(ceiling(ncol(xs[[1]])/3+1), ceiling(ncol(xs[[1]])*2/3)))]
            }
        } else if (test.position == 'end'){
          for (id in seq(1, length(xs))){
            xs[[id]][, c(1, seq(ceiling(ncol(xs[[1]])*2/3+1), ceiling(ncol(xs[[1]]))))] <- xs_test[[id]][, c(1, seq(ceiling(ncol(xs[[1]])*2/3+1), ceiling(ncol(xs[[1]]))))]
          }
        }
    }
     })
    # user  system elapsed 
    # 0.008   0.000   0.008 
  
    ## chunk 4.3------------
    system.time({
    # --------------
    # change here <<
    # --------------
    ## chunk 4.3.1 ----------
    system.time({
    ## initialize tau
    phi_xs <- sapply(names(xs),function(ss) phi[[ss]] %*% t(xs[[ss]]))
    phi_xs_rbind <- do.call(rbind,phi_xs)
    phi_xs_rbind <- phi_xs_rbind[colnames(expr),]
    beta <- sexpr %*% (phi_xs_rbind %*% chol2inv(chol(crossprod(phi_xs_rbind))))
    })
    # user  system elapsed 
    # 1.098   0.001   1.098 

    ## chunk 4.3.2 ----------
    system.time({  
    sexpr <- sapply(names(xs),function(ss) sexpr[,sname[[ss]],drop=F],simplify = F)
    })
  #     user  system elapsed 
  # 0.515   0.498   1.015 
    ## chunk 4.3.3 ----------
    system.time({
    ## initialize tau
    indfit <- sapply(names(xs), function(ss){
      sexpr[[ss]] %*% (phi[[ss]] %*% chol2inv(chol(crossprod(phi[[ss]]))))
    },simplify = F)
    })
  #        user  system elapsed 
  # 1.077   0.000   1.077
    ## chunk 4.3.4 ----------
    system.time({  
    indfitdiff <- sapply(names(xs), function(ss){
      indfit[[ss]] - beta %*% xs[[ss]]
    },simplify = F)
    })
  #        user  system elapsed 
  # 0.008   0.000   0.007 
    ## chunk 4.3.5 ----------
    system.time({  
    tau <- sapply(row.names(indfitdiff[[1]]),function(rid) {
      system.time({  
      m <- sapply(names(xs),function(n) indfitdiff[[n]][rid,]) 
      }) ## 0.002 s
      system.time({  
      m <- m-rowMeans(m) 
      }) ## 0s
      system.time({  
      tcrossprod(m) / (ncol(m)-1) + diag(nrow(m)) * 0.01
      })  ## 0s
    })
    })
  # user  system elapsed 
  # 6.560   0.590   7.151 
    ## chunk 4.3.6 ----------
    system.time({  
    ## initialize gamma
    gamma <- sapply(names(xs), function(ss){
      diff <- sexpr[[ss]] - tcrossprod(indfit[[ss]],phi[[ss]])
      m <- rowMeans(diff)
      rowMeans(diff * diff) - m * m
    })
    if (is.vector(gamma)) {
      gamma <- t(gamma)
      colnames(gamma) <- names(xs)
      row.names(gamma) <- gid
    }
    gamma[gamma < 0.01] <- 0.01
    })
  #        user  system elapsed 
  # 2.595   0.521   3.115
    })
   # user  system elapsed 
   # 12.066   2.650  14.718
    
    ## chunk 4.4, *** very slow *** -------------
     system.time({
    iter <- 0
    ll <- rep(-Inf,length(gid))
    names(ll) <- gid
    gidr <- gid
    while (iter < EMmaxiter && length(gidr) > 0) {
      ## chunk 4.4.1 ------------------
      system.time({
      llold <- ll
      sexpr <- sapply(sexpr,function(i) i[gidr,,drop=F],simplify = F)
      E_phi_u_e <- phi_xs_beta <- phi_tau_phi <- list()
      M <- matrix(0,nrow=length(gidr),ncol=nrow(design),dimnames=list(gidr,names(sname)))
      for (ss in names(sname)) {
        phi_tau_phi[[ss]] <- crossprod(tau[,gidr,drop=F], phicrossprod[[ss]])
        phi_xs_beta[[ss]] <- tcrossprod(beta[gidr,,drop=F],phi_xs[[ss]])
        sigma1s <- 1/(1/phi_tau_phi[[ss]] + 1/gamma[gidr,ss])
        E_phi_u_e[[ss]] <- sigma1s * (phi_xs_beta[[ss]]/phi_tau_phi[[ss]] + sexpr[[ss]]/gamma[gidr,ss]) - phi_xs_beta[[ss]]
        phi_xs_beta[[ss]] <- sexpr[[ss]] - phi_xs_beta[[ss]]
        M[,ss] <- (rowSums((phi_xs_beta[[ss]] - E_phi_u_e[[ss]])^2) + rowSums(sigma1s))/ncol(sigma1s)
      }
      })
      # user  system elapsed 
      # 11.892   1.599  13.493 
      ## chunk 4.4.2 ------------------
      system.time({
      tauinv <- sapply(gidr,function(i) {
        chol2inv(chol(matrix(tau[,i],nrow=ncol(phi[[1]]))))
      },simplify = F)
      })
      # user  system elapsed 
      # 0.569   0.042   0.611 
      ## chunk 4.4.3 ------------------
      system.time({
      E_u_con_theta <- sapply(names(sname),function(ss) {
        m <- matrix(rowSums(phicrossprod[[ss]]),nrow=ncol(phi[[ss]]))
        tmp <- sapply(gidr,function(i) {
          chol2inv(chol(m/gamma[i,ss] + tauinv[[i]]))
        },simplify = F)
      },simplify = F)
      })
      # user  system elapsed 
      # 18.873   0.058  18.934 
      ## chunk 4.4.3.2 ------------------
      system.time({
      k_s <- sapply(names(sname),function(ss){
        tmp <- phi_xs_beta[[ss]] %*% phi[[ss]]
        rowsum(as.vector(t(tmp)) * do.call(rbind,E_u_con_theta[[ss]]),rep(1:nrow(tmp),each=ncol(tmp))) / gamma[gidr,ss]
      },simplify = F)
      })
      # user  system elapsed 
      # 1.386   0.100   1.487
      ## chunk 4.4.4 ------------------
      system.time({
      tmp2 <- (1/gamma[gidr,,drop=F]) %*% t(sapply(names(xs), function(ss){
        crossprod(phi_xs[[ss]])
      }))
      })
      # user  system elapsed 
      # 0.007   0.000   0.007
      ## chunk 4.4.5 ------------------
      system.time({
      tmp3 <- lapply(names(xs), function(ss){
        (sexpr[[ss]] - E_phi_u_e[[ss]]) %*% phi_xs[[ss]] / gamma[gidr,ss]
      })
      tmp3 <- Reduce('+', tmp3)
      newbeta <- t(sapply(row.names(tmp2),function(i) chol2inv(chol(matrix(tmp2[i,],nrow=ncol(beta)))) %*% tmp3[i,]))
      })
      # user  system elapsed 
      # 2.355   0.357   2.712 
      ### gamma
      newgamma <- M
      
      ## chunk 4.4.6 ------------------
      system.time({
      ### tau
      tmp <- lapply(names(xs), function(ss) {
        sapply(E_u_con_theta[[ss]],as.vector) + apply(k_s[[ss]],1,tcrossprod)
      })
      newtau <- Reduce('+', tmp)/length(xs)
      })
      #  user  system elapsed 
      # 2.248   0.001   2.250 
      ## chunk 4.4.7 ------------------
      system.time({
      ### log L
      tmpll <- sapply(names(xs), function(ss) {
        tmp <- 2 * (phi_tau_phi[[ss]] + gamma[gidr,ss])
        rowSums(-log(pi * tmp))/2 - rowSums(phi_xs_beta[[ss]]*phi_xs_beta[[ss]] / tmp)
      })
      })
     # user  system elapsed 
     # 12.262   1.108  13.372 
      ## chunk 4.4.8 ------------------
      system.time({
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
      })
      #   user  system elapsed 
      # 0.004   0.000   0.004
    }
    })
    
    return(list(beta = beta, gamma = gamma, tau = tau, logL = ll))
  }
  if (ncores!=1) {
    allres <- mclapply(unique(knotnum),sfit,mc.cores=ncores)
  } else {
    allres <- lapply(unique(knotnum),sfit)
  }
  })
  # user   system  elapsed 
  # 2957.508  311.263 3269.204 
  
  ## chunk 5 -------------
  system.time({
  para <- list()
  for (i in 1:length(allres)) {
    for (j in row.names(allres[[i]][[1]])) {
      para[[j]] <- list(beta=allres[[i]][[1]][j,],
                        gamma=allres[[i]][[2]][j,],
                        tau=matrix(allres[[i]][[3]][,j],nrow=sqrt(length(allres[[i]][[3]][,j]))),
                        ll=allres[[i]][[4]][j])
    }
  }
  })
  list(parameter=para,knotnum=knotnum)
}

