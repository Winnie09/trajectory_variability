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

## initialize data
num.sample <- length(unique(Cellanno$sample))  
num.coef <- num.knot + 4
set.seed(12345)
tau <- diag(num.coef)
gamma <- rep(1, num.sample)
names(gamma) <- row.names(Design)
# beta <- rep(0, num.sample * 2)
id = grep(GeneName, rownames(GeneByCellExpression))
X_s <- Y_s <- Phi_s <- list() 
for (i in seq(1, nrow(Design))){
  if(Design[i,]==1){
    X_s[[i]] <- kronecker(diag(num.knot + 4),matrix(c(1,1),nrow=1)) 
  } else {
    X_s[[i]] <- kronecker(diag(num.knot + 4),matrix(c(1,0),nrow=1))
  }
  Y_s[[i]] <- GeneByCellExpression[grep(GeneName, rownames(GeneByCellExpression)), grep(rownames(Design)[i], colnames(GeneByCellExpression))]
  
  Phi_s[[i]] <- base[rownames(base)[grep(rownames(Design)[i], rownames(base))], ]
}
names(X_s) <- names(Y_s) <- names(Phi_s) <- row.names(Design)
X <- do.call(rbind, X_s)
Y <- matrix(unlist(Y_s), ncol=1)
library(Matrix)
phi <- as.matrix(.bdiag(Phi_s))
X_til <- phi %*% X
beta <- solve(t(X_til) %*% X_til) %*% t(X_til) %*% Y

betalist  <- taulist <- gammalist <- list()
for (dump in 1:5) {
  print('interative round: ')
  print(dump)
  print('\n')
  ## postier E_e_con_u_theta
  get_mu1 <- function(mu0, sigma0s,x, sigmas){
    1/(1/sigma0s + ncol(x)/sigmas) * (mu0/sigma0s + rowSums(x)/sigmas)
  }
  
  get_sigma1s <- function(sigma0s, n, sigmas){
    as.numeric(1/(1/sigma0s + n/sigmas))
  }
  
  get_E_phi_u_e_cs <- function(c, ss){
    phi_cs <- t(Phi_s[[ss]][c,,drop=F])
    x_s <- t(X_s[[ss]])
    mu1 <- get_mu1(mu0 = t(phi_cs) %*% t(x_s) %*% beta ,
                   sigma0s = t(phi_cs) %*% tau %*% phi_cs,
                   x = GeneByCellExpression[id, Cellanno$sample == ss, drop=F],
                   sigmas = gamma[ss])
    
    as.numeric(mu1 - t(phi_cs) %*% t(x_s) %*% beta)
  }
  
  get_E_phi_u_e_cs2 <- function(c, ss){
    sigma1s <- get_sigma1s(sigma0s = t(phi_cs) %*% tau %*% phi_cs,
                           n = ncol(GeneByCellExpression[id, Cellanno$sample == ss, drop=F]),
                           sigmas = gamma[ss])
    sigma1s - get_E_phi_u_e_cs(c,ss)^2
  }
  
  get_M_cs <- function(c, ss){  
    phi_cs <- t(Phi_s[[ss]][c,,drop=F])
    E_phi_u_e_cs <- get_E_phi_u_e_cs(c,ss)
    sigma1s <- get_sigma1s(sigma0s = t(phi_cs) %*% tau %*% phi_cs,
                           n = ncol(GeneByCellExpression[id, Cellanno$sample == ss, drop=F]),
                           sigmas = gamma[ss])
    E_phi_u_e_cs2 <- sigma1s + E_phi_u_e_cs^2
    e_cs <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
    phi_cs <- t(Phi_s[[ss]][c,,drop=F])
    x_s <- t(X_s[[ss]])
    tmp <- e_cs^2 + t(t(phi_cs) %*% t(x_s) %*% beta) %*% (t(phi_cs) %*% t(x_s) %*% beta) + E_phi_u_e_cs2 - 2* e_cs * (t(phi_cs) %*% t(x_s) %*% beta) - 2 * (e_cs - (t(phi_cs) %*% t(x_s) %*% beta)) * E_phi_u_e_cs
    as.numeric(tmp)
  }
  
  get_A_cs <- function(c, ss){
    M_cs <- get_M_cs(c, ss)
    log(2 * pi * gamma[ss]) + M_cs/(2 * gamma[ss])
  }
  
  # c <- colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F])[1]
  # phi_cs <- t(Phi_s[[ss]][c,,drop=F])
  # x_s <- t(X_s[[ss]])
  # e_cs <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
  # A <- sapply(names(X_s), function(ss){
  #   print(ss)          ##########
  #   sapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
  #     get_A_cs(c=c, ss=ss)    
  #   })
  # })
  # A <- -sum(unlist(A))/2
  # 
  ### get posterior E_u_con_theta
  get_Lambda_s <- function(ss){
    tmp <- lapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
      t(Phi_s[[ss]][c,,drop=F]) %*% Phi_s[[ss]][c,,drop=F]
    })
    tmp <- Reduce('+', tmp)
    return(solve(tmp/gamma[ss] + solve(tau)))
  }
  
  #original
  # get_k_s <- function(ss){
  #   x_s <- t(X_s[[ss]])
  #   tmp <- lapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
  #     phi_cs_tmp <- t(Phi_s[[ss]][c,,drop=F])
  #     e_cs_tmp <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
  #     (e_cs_tmp - t(phi_cs_tmp) %*% t(x_s) %*% beta) %*% t(phi_cs_tmp)
  #   })
  #   return(t(Reduce('+', tmp)))
  # }

  get_k_s <- function(ss){
    x_sbeta <- X_s[[ss]] %*% beta
    tmp <- lapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
      (GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c] - sum(Phi_s[[ss]][c,,drop=F][1,] * x_sbeta[,1])) * Phi_s[[ss]][c,,drop=F]
    })
    return(t(Reduce('+', tmp) %*% get_Lambda_s(ss) / gamma[ss]))
  }

  get_E_posterior_u_tau_u <- function(ss){
    Lambda_s <- get_Lambda_s(ss)
    k_s <- get_k_s(ss)
    return(sum(diag(solve(tau) %*% Lambda_s)) + t(k_s) %*% solve(tau) %*% k_s)
  }
  
  # B <- length(unique(Cellanno$sample)) * log( (2*pi)^(-num.coef/2) * (det(tau))^(-1/2)) -
  #   sum(sapply(unique(Cellanno$sample), function(ss){
  #     print(ss)          ##################
  #     get_E_posterior_u_tau_u(ss)
  #   }))/2  ###
  # Q = A + B  
  # 
  
  ##### M step
  ### beta
  
  tmp2 <- lapply(names(X_s), function(ss){
    print(ss)
    tmp <- lapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
      phi_cs_tmp <- t(Phi_s[[ss]][c,,drop=F])
      x_s_tmp <- t(X_s[[ss]])
      x_s_tmp %*% phi_cs_tmp %*% t(phi_cs_tmp) %*% t(x_s_tmp) / gamma[ss] 
    })
    tmp <- Reduce('+', tmp)  
  })
  tmp2 <- Reduce('+', tmp2)     
  
  tmp3 <- lapply(names(X_s), function(ss){
    print(ss)
    tmp <- lapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
      phi_cs_tmp <- t(Phi_s[[ss]][c,,drop=F])
      x_s_tmp <- t(X_s[[ss]])
      e_cs_tmp <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
      e_cs_tmp * t(phi_cs_tmp) %*% t(x_s_tmp) - get_E_phi_u_e_cs(c,ss) * (t(phi_cs_tmp) %*% t(x_s_tmp)) / gamma[ss]
    })
    tmp <- Reduce('+', tmp)  
  })
  tmp3 <- Reduce('+', tmp3)     
  newbeta <- solve(tmp2) %*% t(tmp3)
  
  ### gamma
  newgamma <- sapply(names(X_s), function(ss){
    tmp <- sapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
      get_M_cs(c, ss)
    })
    sum(unlist(tmp))/(2 * ncol(GeneByCellExpression[,Cellanno$sample==ss,drop=F]))
  })
  ### tau
  tmp <- lapply(names(X_s), function(ss){
    Lambda_s <- get_Lambda_s(ss)
    k_s <- get_k_s(ss)
    Lambda_s + k_s %*% t(k_s)
  })
  newtau <- Reduce('+', tmp)/length(X_s)
  
  beta <- newbeta
  gamma <- newgamma
  tau <- newtau
  betalist[[length(betalist)+1]] <- beta
  taulist[[length(taulist)+1]] <- tau
  gammalist[[length(gammalist)+1]] <- gamma
}  

##### plot
GeneName <- 'NCBP3'
pd = data.frame(sample=Cellanno$sample, e = expr[grepl(GeneName, rownames(expr)),], time = psn[Cellanno$cell], group=sapply(Cellanno$sample, function(i) ifelse(i %in% c('BM1','BM2','BM5','BM6'), 'group1','group2')))
ggplot() + geom_point(data=pd, aes(x=time, y = e, color=group), size=.3)


