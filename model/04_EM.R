##### E-M algorithm
rm(list=ls())
res <- readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/model/model_testing_data.rds')
expr = res$expr
selgene <- res$selgene
pseudotime  <- res$pseudotime
design = res$design
rm(res)
##### implement model
### full model
psn <- 1:length(pseudotime)
names(psn) <- pseudotime

GeneByCellExpression=expr
Pseudotime=psn
Design=design
num.knot = 3
GeneName = 'ACE2'
DebugMode = FALSE
Cellanno <- data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)))


fullmodel <- function(GeneByCellExpression, Pseudotime, Design, num.knot = 3, GeneName, DebugMode = FALSE){
  ## GeneByCellExpression: cell names are pasted sample number in front with a ":" as seperation
  ## Pseudotime: a numeric vecotor, names the same as expression columns
  ## Design: sample by feature design matrix. rownames are sample names. first column is the group partition, currently only one column.
  
  knots = seq(min(psn),max(psn),length.out=num.knot+2)[2:(num.knot+1)]
  base <- lapply(row.names(Design),function(sg) {  
    n <- Cellanno[Cellanno$sample==sg, 'cell']
    cbind(1,bs(psn[n],knots = knots))
  })
  base <- do.call(rbind, base)
  base <- base[colnames(GeneByCellExpression), ]
  # print('printing base:\n')
  # print(tail(base,))
  
## initialize data
num.sample <- length(unique(Cellanno$sample))  
num.coef <- num.knot + 4
set.seed(12345)
tau <- diag(num.coef)
gamma <- rep(1, num.sample)
names(gamma) <- row.names(Design)
beta <- rep(0,3)
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
  
  
  print('printing Y: \n')
  print(Y[1:20])  
  print(head(names(Y)))
  
    
  
  ## postier E_e_con_u_theta
  get_mu1 <- function(mu0, sigma0s,x, sigmas){
    1/(1/sigma0s + ncol(x)/sigmas) * (mu0/sigma0s + rowSums(x)/sigmas)
  }
  
  get_sigma1s <- function(sigma0s, n, sigmas){
    1/(1/sigma0s + n/sigmas)
  }
  
  get_A_cs <- function(c, ss){
    mu1 <- get_mu1(mu0 = t(phi_cs) %*% t(x_s) %*% beta ,
                   sigma0s = t(phi_cs) %*% tau %*% phi_cs,
                   x = GeneByCellExpression[id, alls == ss, drop=F],
                   sigmas = gamma[ss])
    sigma1s <- get_sigma1s(sigma0s = t(phi_cs) %*% tau %*% phi_cs,
                           n = ncol(GeneByCellExpression[id, alls == ss, drop=F]),
                           sigmas = gamma[ss])
    E_phi_u_e_cs <- mu1 - t(phi_cs) %*% t(x_s) %*% beta
    E_phi_u_e_cs2 <- sigma1s - E_phi_u_e_cs^2
    M_cs <- e_cs^2 + t(t(phi_cs) %*% t(x_s) %*% beta) %*% (t(phi_cs) %*% t(x_s) %*% beta) + E_phi_u_e_cs2 - 2* e_cs * (t(phi_cs) %*% t(x_s) %*% beta) - 2 * (e_cs - (t(phi_cs) %*% t(x_s) %*% beta)) * E_phi_u_e_cs
    log(2 * pi * gamma[ss]) + M_cs/(2 * gamma[ss])
  }
  c <- colnames(GeneByCellExpression[,alls==ss,drop=F])[1]
  phi_cs <- t(Phi_s[[ss]][c,,drop=F])
  x_s <- t(X_s[[ss]])
  e_cs <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
  A <- sapply(names(X_s), function(ss){
    sapply(colnames(GeneByCellExpression[,alls==ss,drop=F]), function(c){
      get_A_cs(c=c, ss=ss)    
    })
  })
  A <- -sum(unlist(A))/2
      
  ### get posterior E_u_con_theta
  get_E_posterior_u_tau_u <- function(ss){
    tmp <- lapply(colnames(GeneByCellExpression[,alls==ss,drop=F]), function(c){
      t(Phi_s[[ss]][c,,drop=F])
    })
    tmp <- Reduce('+', tmp)
    Lambda_s <- solve((tmp %*% t(tmp)) /gamma[ss] + solve(tau))
    
    c <- colnames(GeneByCellExpression[,alls==ss,drop=F])[1]
    phi_cs <- t(Phi_s[[ss]][c,,drop=F])
    e_cs <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
    
    tmp <- lapply(colnames(GeneByCellExpression[,alls==ss,drop=F]), function(c){
      phi_cs_tmp <- t(Phi_s[[ss]][c,,drop=F])
      e_cs_tmp <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
      (e_cs_tmp - t(phi_cs_tmp) %*% t(x_s) %*% beta) %*% t(phi_cs_tmp)
    })
    k_s <- t(Reduce('+', tmp))
    return(sum(diag(solve(tau) %*% Lambda_s)) + t(k_s) %*% solve(tau) %*% k_s)
  }
    
  B <- length(unique(alls)) * log( (2*pi)^(-num.coef/2) * (det(tau))^(-1/2)) -
    sum(sapply(unique(alls), function(ss){
    print(ss)
    get_E_posterior_u_tau_u(ss)
  }))/2  ###
  Q = A + B  
  
  
  ##### M step
  tmp <- lapply(colnames(GeneByCellExpression[,alls==ss,drop=F]), function(c){
      phi_cs_tmp <- t(Phi_s[[ss]][c,,drop=F])
      x_s_tmp <- t(X_s[[ss]])
      
      1/gamma[ss] * (t(phi_cs_tmp) %*% t(x_s_tmp) )%*% t((t(phi_cs_tmp) %*% t(x_s_tmp) ))
  })
    
    
     t(phi_cs_tmp) %*% t(x_s_tmp)
  
  
  
  
  
  
  
  
  
  #############
  ## get w_gs
  get_W_s <- function(beta){
    W_s <- list()
    for (i in seq(1, nrow(Design))){
      y_s <- matrix(Y_s[[i]], ncol=1)
      phi_s <- Phi_s[[i]]
      a <- t(y_s - phi_s %*% (X_s[[i]] %*% beta)) %*% (y_s - phi_s %*% (X_s[[i]] %*% beta))
      W_s[[i]] <- length(y_s)/a
    }
    return(W_s)
  }
  
  get_W_vec <- function(W_s){
    W_vec <- rep(as.vector(unlist(W_s)), sapply(Y_s,length))
  }
    
  ## get log-Likelihood
  get_logL <- function(W_s){
    log_L <- sum(sapply(seq(1,length(unique(sample))), function(i){
      sum(sapply(seq(1,length(W_s[[i]])), function(c){
        log(W_s[[i]][c] * 1/sqrt(2*pi)) - (1/W_s[[i]][c])/(2 * W_s[[i]][c])  
      }))
    }))
    return(log_L)
  }
    
  ## update beta
  update_beta <- function(W_vec){
    X_til_new <- t(sapply(seq(1,length(W_vec)), function(i){
      W_vec[i] %*% X_til[i,] 
    }))
    beta_new <- chol2inv(chol(t(X_til_new) %*% X_til)) %*% t(X_til_new) %*% Y
  }
  
  ## interative algo
  interative_beta <- function(beta_initial){
    beta <- beta_initial
    logl <- c(1, get_logL(get_W_s(beta)))
    while (abs(logl[length(logl)]-logl[length(logl)-1]) > 1e-10) {
      beta <- update_beta(get_W_vec(get_W_s(beta)))
      W_s <- get_W_s(beta)
      logl <- c(logl,get_logL(W_s))
      print(logl[length(logl)])
    }
    return(list(beta=beta, W=get_W_vec(W_s)))
  }
  
  if (DebugMode){
    interative_beta2 <- function(beta_initial){  
      beta_final <- beta <- beta_initial
      logl <- c(1, get_logL(get_W_s(beta)))
      while (abs(logl[length(logl)]-logl[length(logl)-1]) > 1e-10 & ncol(beta_final) < 20) {
        beta <- update_beta(get_W_vec(get_W_s(beta)))
        W_s <- get_W_s(beta)
        logl <- c(logl,get_logL(W_s))
        print(logl[length(logl)])
        beta_final <- cbind(beta_final, beta)
      } 
      return(beta_final)
    }
    return(interative_beta2(beta))
  } else {
    return(interative_beta(beta))
  }
}



b <- fullmodel(GeneByCellExpression=expr2, Pseudotime=pseudotime, Design=design, num.base=3, GeneName=sub(':.*','',selgene[5]), DebugMode=F)

##### simulate data for a double-check
set.seed(123)
beta <- rnorm(14)
num.base = 3
Pseudotime <- psn
Design <- design
knots = seq(min(Pseudotime),max(Pseudotime),length.out=num.base+2)[2:(num.base+1)]
base <- lapply(row.names(Design),function(sg) {
  n <- names(sample)[sample==sg]
  cbind(1,bs(Pseudotime[n],knots = knots))
})
base <- do.call(rbind, base)
GeneByCellExpression <- expr2
GeneName <- 'ACE2'
X_s <- Phi_s <- list() 
for (i in seq(1, nrow(design))){
  if(design[i,]==1){
    X_s[[i]] <- kronecker(diag(num.base + 4),matrix(c(1,1),nrow=1)) 
  } else {
    X_s[[i]] <- kronecker(diag(num.base + 4),matrix(c(1,0),nrow=1))
  }
  Phi_s[[i]] <- base[rownames(base)[grep(rownames(design)[i], rownames(base))], ]
}
X <- do.call(rbind, X_s)
library(Matrix)
phi <- as.matrix(.bdiag(Phi_s))

tau <- diag(nrow(X))
set.seed(1)
a <- X %*% beta + t(rmvnorm(1,mean=rep(0,nrow(X)),sigma = 0.2*tau))
epsilon <- as.vector(unlist(sapply(seq(1,nrow(design)), function(i) rnorm(table(sample)[rownames(design)[i]], mean = 0, sd = 0.5))))
Y <- t(phi %*% a + epsilon)
colnames(Y) <- colnames(expr2)
rownames(Y) <- 'gene1'
## get Y ordered by samples
Y_s <-  list() 
for (i in seq(1, nrow(Design))){
  Y_s[[i]] <- Y[, grep(rownames(Design)[i], colnames(Y))]
}
Y <- t(matrix(unlist(Y_s), ncol=1))
colnames(Y) <- unlist(sapply(Y_s, names))
rownames(Y) <- 'gene1'
res <- fullmodel(GeneByCellExpression=Y, Pseudotime=psn, Design=design, num.base=3, GeneName='gene1', DebugMode=F)
b = res$beta
W_out = res$W

a <- X %*% b + t(rmvnorm(1,mean=rep(0,nrow(X)),sigma=1e-6*tau))
Y_check <- t(phi %*% a + epsilon)
var(Y_check[1,])

truevar <- 0.2 * t(phi) %*% phi + (0.5)^2

1e-6 * tau



