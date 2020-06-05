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

fullmodel <- function(GeneByCellExpression, Cellanno=Cellanno, Pseudotime, Design, GeneName, num.knot=3){
  psn <- 1:length(Pseudotime)
  names(psn) <- Pseudotime
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
  id = grep(GeneName, rownames(GeneByCellExpression))
  X_s <- Y_s <- Phi_s <- list() 
  for (i in seq(1, nrow(Design))){
    if(Design[i,]==1){
      X_s[[i]] <- kronecker(diag(num.knot + 4),matrix(c(1,1),nrow=1)) 
    } else {
      X_s[[i]] <- kronecker(diag(num.knot + 4),matrix(c(1,0),nrow=1))
    }
    Y_s[[i]] <- GeneByCellExpression[grep(GeneName, rownames(GeneByCellExpression)), grep(rownames(Design)[i], colnames(GeneByCellExpression)), drop=F]
    
    Phi_s[[i]] <- base[rownames(base)[grep(rownames(Design)[i], rownames(base))], ]
  }
  names(X_s) <- names(Y_s) <- names(Phi_s) <- row.names(Design)
  
  X <- do.call(rbind, X_s) 
  Y <- matrix(unlist(Y_s), ncol=1)
  
  library(Matrix)
  phi <- as.matrix(.bdiag(Phi_s))
  X_til <- phi %*% X
  beta <- solve(t(X_til) %*% X_til) %*% t(X_til) %*% Y
  
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
  
  get_M_cs <- function(c, ss){  
    phi_cs <- t(Phi_s[[ss]][c,,drop=F])
    E_phi_u_e_cs <- get_E_phi_u_e_cs(c,ss)
    sigma1s <- get_sigma1s(sigma0s = t(phi_cs) %*% tau %*% phi_cs,
                           n = sum(Cellanno$sample == ss),
                           sigmas = gamma[ss])
    E_phi_u_e_cs2 <- sigma1s + E_phi_u_e_cs^2  ###
    e_cs <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
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
    tmp <- lapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
      t(Phi_s[[ss]][c,,drop=F]) %*% Phi_s[[ss]][c,,drop=F]
    })
    tmp <- Reduce('+', tmp)
    return(solve(tmp/gamma[ss] + solve(tau)))
  }
  
  #original
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

  betalist  <- taulist <- gammalist <- logLlist <- list()
  for (dump in seq(1,5)){
    print('interation')
    print(dump)
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
        phi_cs_x_s_tmp <- Phi_s[[ss]][c,,drop=F] %*% X_s[[ss]]
        e_cs_tmp <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
        (e_cs_tmp * phi_cs_x_s_tmp - get_E_phi_u_e_cs(c,ss) * phi_cs_x_s_tmp) / gamma[ss]
      })
      tmp <- Reduce('+', tmp)  
    })
    tmp3 <- Reduce('+', tmp3)     
    newbeta <- solve(tmp2) %*% t(tmp3)
    
    ### gamma
    newgamma <- sapply(names(X_s), function(ss){   ### long
      tmp <- sapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
        get_M_cs(c, ss)
      })
      sum(unlist(tmp))/(2 * ncol(GeneByCellExpression[,Cellanno$sample==ss,drop=F]))
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
      tmp <- sapply(colnames(GeneByCellExpression[,Cellanno$sample==ss,drop=F]), function(c){
          phi_cs_tmp <- t(Phi_s[[ss]][c,,drop=F])
          x_s_tmp <- t(X_s[[ss]])
          e_cs_tmp <- GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), c]
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



##### check u
e_s_tmp <- t(GeneByCellExpression[grep(GeneName,row.names(GeneByCellExpression)), Cellanno$sample==ss,drop=F])
solve(t(Phi_s[[ss]])%*% Phi_s[[ss]]) %*% t(Phi_s[[ss]]) %*% e_s_tmp - X_s[[ss]] %*% beta


##### plot
GeneName <- 'NCBP3'
pd = data.frame(sample=Cellanno$sample, e = expr[grepl(GeneName, rownames(expr)),], time = psn[Cellanno$cell], group=sapply(Cellanno$sample, function(i) ifelse(i %in% c('BM1','BM2','BM5','BM6'), 'group1','group2')))
ggplot() + geom_point(data=pd, aes(x=time, y = e, color=sample), size=.3) + facet_grid(~sample)


##### simulate data for a double-check
set.seed(123)
beta <- rnorm(14)
num.base = 3
Design <- design
psn = seq(1,length(pseudotime))
names(psn) <- pseudotime
Pseudotime <- psn
knots = seq(min(Pseudotime),max(Pseudotime),length.out=num.base+2)[2:(num.base+1)]
base <- lapply(row.names(Design),function(sg) {
  n <- Cellanno[Cellanno$sample==sg, 'cell']
  cbind(1,bs(Pseudotime[n],knots = knots))
})

base <- do.call(rbind, base)
GeneName <- 'NCBP3'
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
library(mvtnorm)
tau <- diag(nrow(X))
set.seed(1)
a <- X %*% beta + t(rmvnorm(1,mean=rep(0,nrow(X)),sigma = 0.001*tau))
epsilon <- as.vector(unlist(sapply(seq(1,nrow(design)), function(i) rnorm(table(Cellanno$sample)[rownames(design)[i]], mean = 0, sd = 0.001))))
Y <- t(phi %*% a + epsilon)
colnames(Y) <- colnames(expr)
rownames(Y) <- 'gene1'

## get Y ordered by samples
Y_s <-  list() 
for (i in seq(1, nrow(Design))){
  Y_s[[i]] <- Y[, grep(rownames(Design)[i], colnames(Y))]
}
Y <- t(matrix(unlist(Y_s), ncol=1))
colnames(Y) <- unlist(sapply(Y_s, names))
rownames(Y) <- 'gene1'

res <- fullmodel(GeneByCellExpression=Y, Pseudotime=pseudotime, Design=design, GeneName='gene1')


b = res$beta
W_out = res$W

a <- X %*% b + t(rmvnorm(1,mean=rep(0,nrow(X)),sigma=1e-6*tau))
Y_check <- t(phi %*% a + epsilon)
var(Y_check[1,])

truevar <- 0.2 * t(phi) %*% phi + (0.5)^2

1e-6 * tau


