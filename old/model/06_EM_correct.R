##### E-M algorithm
rm(list=ls())
res <- readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/model/model_testing_data.rds')
expr = res$expr
selgene <- res$selgene
pseudotime  <- res$pseudotime
design = res$design

##### implement model
### full model

GeneByCellExpression=expr
Design=design
num.knot = 3
GeneName <- 'NCBP3'
DebugMode = FALSE
Cellanno <- data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)

# fullmodel <- function(GeneByCellExpression, Pseudotime, Design, num.knot = 3, GeneName, DebugMode = FALSE){
## GeneByCellExpression: cell names are pasted sample number in front with a ":" as seperation
## Pseudotime: a numeric vecotor, names the same as expression columns
## Design: sample by feature design matrix. rownames are sample names. first column is the group partition, currently only one column.

psn <- 1:length(pseudotime)
names(psn) <- pseudotime
knots = seq(min(psn),max(psn),length.out=num.knot+2)[2:(num.knot+1)]
base <- lapply(row.names(Design),function(sg) {  
  n <- Cellanno[Cellanno$sample==sg, 'cell']
  cbind(1,bs(psn[n],knots = knots))
})
base <- do.call(rbind, base)
base <- base[colnames(GeneByCellExpression), ]

#Y <- GeneByCellExpression[GeneName,]

## initialize data
num.sample <- length(unique(Cellanno$sample))  
num.coef <- num.knot + 4
set.seed(12345)

library(Matrix)
library(mvtnorm)
truetau <- tau <- diag(num.coef)/10
truegamma <- gamma <- abs(rnorm(nrow(design)))/10
names(gamma) <- row.names(Design)
truebeta <- beta <- rnorm(14)


id = grep(GeneName, rownames(GeneByCellExpression))
X_s <- Y_s <- Phi_s <- list() 
for (i in seq(1, nrow(Design))){
  if(Design[i,]==1){
    X_s[[i]] <- kronecker(diag(num.knot + 4),matrix(c(1,1),nrow=1)) 
  } else {
    X_s[[i]] <- kronecker(diag(num.knot + 4),matrix(c(1,0),nrow=1))
  }
  Phi_s[[i]] <- base[rownames(base)[grep(rownames(Design)[i], rownames(base))], ]
  a <- X_s[[i]] %*% beta + t(rmvnorm(1,mean=rep(0,nrow(X_s[[i]])),sigma = tau))
  epsilon <- rnorm(nrow(Phi_s[[i]]), mean = 0, sd = sqrt(gamma[[i]]))
  Y_s[[i]] <- Phi_s[[i]] %*% a + epsilon  
}
names(X_s) <- names(Phi_s) <- row.names(Design)

X <- do.call(rbind, X_s) 
Y <- do.call(rbind,Y_s)[,1]

phi <- as.matrix(.bdiag(Phi_s))
X_til <- phi %*% X
beta <- solve(t(X_til) %*% X_til) %*% t(X_til) %*% Y

## postier E_e_con_u_theta
get_mu1 <- function(mu0, sigma0s,x, sigmas){
  (1/(1/sigma0s + 1/sigmas)) * (mu0/sigma0s + x/sigmas)
}

get_sigma1s <- function(sigma0s,n, sigmas){
  as.numeric(1/(1/sigma0s + 1/sigmas))
}

get_E_phi_u_e_cs <- function(c, ss){
  phi_cs <- t(Phi_s[[ss]][c,,drop=F])
  x_s <- t(X_s[[ss]])
  mu1 <- get_mu1(mu0 = t(phi_cs) %*% t(x_s) %*% beta ,
                 sigma0s = t(phi_cs) %*% tau %*% phi_cs,
                 x = Y[c],
                 sigmas = gamma[ss])
  as.numeric(mu1 - t(phi_cs) %*% t(x_s) %*% beta)
}

get_M_cs <- function(c, ss){  
  phi_cs <- t(Phi_s[[ss]][c,,drop=F])
  E_phi_u_e_cs <- get_E_phi_u_e_cs(c,ss)
  sigma1s <- get_sigma1s(sigma0s = t(phi_cs) %*% tau %*% phi_cs,
                         n = 1,
                         sigmas = gamma[ss])
  E_phi_u_e_cs2 <- sigma1s + E_phi_u_e_cs^2  ###
  e_cs <- Y[c]
  phi_cs <- t(Phi_s[[ss]][c,,drop=F])
  x_s <- t(X_s[[ss]])
  phi_x_beta <- as.numeric(t(phi_cs) %*% t(x_s) %*% beta)
  tmp <- e_cs^2 + phi_x_beta^2 + E_phi_u_e_cs2 - 2* e_cs * phi_x_beta - 2 * (e_cs - phi_x_beta) * E_phi_u_e_cs
  as.numeric(tmp)
}

get_A_cs <- function(c, ss){
  M_cs <- get_M_cs(c, ss)
  log(2 * pi * gamma[ss]) + M_cs/(2 * gamma[ss])
}
### get posterior E_u_con_theta
get_Lambda_s <- function(ss){
  tmp <- lapply(Cellanno[Cellanno$sample==ss,'cell'], function(c){
    t(Phi_s[[ss]][c,,drop=F]) %*% Phi_s[[ss]][c,,drop=F]
  })
  tmp <- Reduce('+', tmp)
  return(solve(tmp/gamma[ss] + solve(tau)))
}

get_k_s <- function(ss){
  x_sbeta <- X_s[[ss]] %*% beta
  tmp <- lapply(Cellanno[Cellanno$sample==ss,'cell'], function(c){
    (Y[c] - sum(Phi_s[[ss]][c,,drop=F][1,] * x_sbeta[,1])) * Phi_s[[ss]][c,,drop=F]
  })
  return(t(Reduce('+', tmp) %*% get_Lambda_s(ss) / gamma[ss]))
}

get_E_posterior_u_tau_u <- function(ss){
  Lambda_s <- get_Lambda_s(ss)
  k_s <- get_k_s(ss)
  return(sum(diag(solve(tau) %*% Lambda_s)) + t(k_s) %*% solve(tau) %*% k_s)
}

betalist  <- taulist <- gammalist <- logLlist <- list()
# for (dump in seq(1,20)){
#   print('interation')
#   print(dump)
while(ifelse(length(logLlist) <= 1, TRUE, logLlist[[length(logLlist)]] - logLlist[[length(logLlist)-1]] >= 1)){
  ##### M step
  ### beta
  tmp2 <- lapply(names(X_s), function(ss){
    print(ss)
    tmp <- lapply(Cellanno[Cellanno$sample==ss,'cell'], function(c){
      phi_cs_tmp <- t(Phi_s[[ss]][c,,drop=F])
      x_s_tmp <- t(X_s[[ss]])
      x_s_tmp %*% phi_cs_tmp %*% t(phi_cs_tmp) %*% t(x_s_tmp) / gamma[ss] 
    })
    tmp <- Reduce('+', tmp)  
  })
  tmp2 <- Reduce('+', tmp2)     
  
  tmp3 <- lapply(names(X_s), function(ss){
    print(ss)
    tmp <- lapply(Cellanno[Cellanno$sample==ss,'cell'], function(c){
      phi_cs_x_s_tmp <- Phi_s[[ss]][c,,drop=F] %*% X_s[[ss]]
      e_cs_tmp <- Y[c]
      (e_cs_tmp * phi_cs_x_s_tmp - get_E_phi_u_e_cs(c,ss) * phi_cs_x_s_tmp) / gamma[ss]
    })
    tmp <- Reduce('+', tmp)  
  })
  tmp3 <- Reduce('+', tmp3)     
  newbeta <- solve(tmp2) %*% t(tmp3)
  
  ### gamma
  newgamma <- sapply(names(X_s), function(ss){   ### long
    tmp <- sapply(Cellanno[Cellanno$sample==ss,'cell'], function(c){
      get_M_cs(c, ss)
    })
    sum(unlist(tmp))/sum(Cellanno$sample==ss)
  })
  ### tau
  tmp <- lapply(names(X_s), function(ss){  ### long
    Lambda_s <- get_Lambda_s(ss)
    k_s <- get_k_s(ss)
    Lambda_s + k_s %*% t(k_s)
  })
  newtau <- Reduce('+', tmp)/length(X_s)
  ### log L
  logL <- sapply(names(X_s), function(ss){   ### long
    tmp <- sapply(Cellanno[Cellanno$sample==ss,'cell'], function(c){
      phi_cs_tmp <- t(Phi_s[[ss]][c,,drop=F])
      x_s_tmp <- t(X_s[[ss]])
      e_cs_tmp <- Y[c]
      phi_tau_phi <-  as.numeric(t(phi_cs_tmp) %*% tau %*% phi_cs_tmp)
      log(1/sqrt(2 * pi * (phi_tau_phi + gamma[ss]))) - ((e_cs_tmp - t(phi_cs_tmp) %*% t(x_s_tmp) %*% beta)^2 / (2 * (phi_tau_phi + gamma[ss])))
    })
    sum(tmp)
  })
  print(newlogL <- sum(logL))
  
  beta <- newbeta
  gamma <- newgamma
  tau <- newtau
  
  betalist[[length(betalist)+1]] <- beta
  taulist[[length(taulist)+1]] <- tau
  gammalist[[length(gammalist)+1]] <- gamma
  logLlist[[length(logLlist)+1]] <- newlogL
}
return(list(beta = betalist, gamma = gammalist, tau = taulist, logL = logLlist))
}


