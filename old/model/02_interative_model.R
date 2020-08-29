rm(list=ls())
# geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
# addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
# addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
# method <- as.character(commandArgs(trailingOnly = T)[[4]])
## be reminded that the pca is done on cv > 0.5 (for all samples) genes. You should redo it !!!
geneProp <- 0.05
addSignalType <- 'linear'
addSignalPara <-  2

# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')

## source functions
# source('./function/01_function.R')
allf <- list.files('./function')
allf = allf[!grepl('01_function.R',allf)]
for (f in allf){  
  source(paste0('./function/', f))
}

# source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
library(parallel)
library(splines)
library(limma)
library(RColorBrewer)
suppressMessages(library(SingleCellExperiment))

expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
allp = sub(':.*','', colnames(expr))

dr <- readRDS('./addsignal/result/dr.rds')
  
library(TSCAN)
pseudotime <- TSCANorder(exprmclust(t(dr),reduce = F,clusternum=4,clustermethod='kmeans'),orderonly = T)
psn <- 1:length(pseudotime)
names(psn) <- pseudotime
order = data.frame(Pseudotime = seq(1, length(pseudotime)), Cell=pseudotime)
set.seed(12345)
selgene <- sample(rownames(expr), round(geneProp*nrow(expr)))

sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'
  

expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara)

##### implement model
### full model

fullmodel <- function(GeneByCellExpression, Pseudotime, Design, num.base = 3, GeneName, DebugMode = FALSE){
  ## GeneByCellExpression: cell names are pasted sample number in front with a ":" as seperation
  ## Pseudotime: a numeric vecotor, names the same as expression columns
  ## Design: sample by feature design matrix. rownames are sample names. first column is the group partition, currently only one column.
  num.base = 3
  knots = seq(min(Pseudotime),max(Pseudotime),length.out=num.base+2)[2:(num.base+1)]
  base <- lapply(row.names(Design),function(sg) {
    n <- names(sample)[sample==sg]
    cbind(1,bs(Pseudotime[n],knots = knots))
  })
  base <- do.call(rbind, base)
  base <- base[colnames(GeneByCellExpression), ]
  # print('printing base:\n')
  # print(tail(base,))
  
  ## initialize beta 
  X_s <- Y_s <- Phi_s <- list() 
  for (i in seq(1, nrow(Design))){
    if(Design[i,]==1){
      X_s[[i]] <- kronecker(diag(num.base + 4),matrix(c(1,1),nrow=1)) 
    } else {
      X_s[[i]] <- kronecker(diag(num.base + 4),matrix(c(1,0),nrow=1))
    }
    Y_s[[i]] <- GeneByCellExpression[grep(GeneName, rownames(GeneByCellExpression)), grep(rownames(Design)[i], colnames(GeneByCellExpression))]
    # Phi_s[[i]] <- base[names(Y_s[[i]]), ]
    Phi_s[[i]] <- base[rownames(base)[grep(rownames(Design)[i], rownames(base))], ]
  }
  X <- do.call(rbind, X_s)
  Y <- matrix(unlist(Y_s), ncol=1)
  library(Matrix)
  phi <- as.matrix(.bdiag(Phi_s))
  X_til <- phi %*% X
  beta <- solve(t(X_til) %*% X_til) %*% t(X_til) %*% Y
  # print('printing X:\n')
  # print(X[1:10,1:10])
  # print('printing phi: \n')
  # print(phi[1:10, 1:10])
  print('printing Y: \n')
  print(Y[1:20])  
  print(head(names(Y)))
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
b <- fullmodel(GeneByCellExpression=expr2, Pseudotime=psn, Design=design, num.base=3, GeneName=sub(':.*','',selgene[5]), DebugMode=F)

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


