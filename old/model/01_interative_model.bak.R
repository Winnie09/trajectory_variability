rm(list=ls())
geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
method <- as.character(commandArgs(trailingOnly = T)[[4]])
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
suppressMessages(library(tradeSeq))
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)


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
num.base = 3
knots = seq(min(psn),max(psn),length.out=num.base+2)[2:(num.base+1)]
base <- lapply(row.names(design),function(sg) {
  n <- names(sample)[sample==sg]
  cbind(1,bs(psn[n],knots = knots))
})
base <- do.call(rbind, base)

## initialize beta 
X_s <- Y_s <- Phi_s <- list() 
for (i in seq(1, length(unique(sample)))){
  if(design[paste0('BM',i),]==1){
    X_s[[i]] <- kronecker(diag(7),matrix(c(1,1),nrow=1)) 
  } else {
    X_s[[i]] <- kronecker(diag(7),matrix(c(1,0),nrow=1))
  }
  Y_s[[i]] <- expr2[grep('ACE2', rownames(expr2)), grep(paste0('BM', i), colnames(expr2))]
  Phi_s[[i]] <- base[names(Y_s[[i]]), ]
}
X <- do.call(rbind, X_s)
Y <- matrix(unlist(Y_s), ncol=1)
library(Matrix)
phi <- as.matrix(.bdiag(Phi_s))
X_til <- phi %*% X
beta <- solve(t(X_til) %*% X_til) %*% t(X_til) %*% Y


## plot to check fitting
# v <- phi %*% X %*% beta 
# t = psn[names(unlist(Y_s))]
# plot(Y[,1]~t, pch=20, xlab='expression', ylab='t', .cex=.1, col=as.factor(sp))
# points(v~t, col='blue', lwd=.5)

## get w_gs
get_W_s <- function(beta){
  W_s <- list()
  for (i in seq(1, length(unique(sample)))){
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
  beta_new <- chol2inv(chol(t(X_til_new) %*% X_til_new)) %*% t(X_til_new) %*% Y
}

## interative algo
interative_beta <- function(beta_initial){
  beta <- beta_initial
  logl <- c(1, get_logL(get_W_s(beta)))
  while (abs(logl[length(logl)]-logl[length(logl)-1]) > 1e-10) {
    beta <- update_beta(get_W_vec(get_W_s(beta)))
    logl <- c(logl,get_logL(get_W_s(beta)))
    print(logl[length(logl)])
  }
  return(beta)
}
interative_beta(beta)


#### simulate data for a double-check
set.seed(123)
beta <- rnorm(14)
X_s <- list() 
for (i in seq(1, length(unique(sample)))){
  if(design[paste0('BM',i),]==1){
    X_s[[i]] <- kronecker(diag(7),matrix(c(1,1),nrow=1)) 
  } else {
    X_s[[i]] <- kronecker(diag(7),matrix(c(1,0),nrow=1))
  }
}
library(mvtnorm)
X <- do.call(rbind, X_s)
tau <- diag(nrow(X))
set.seed(1)
a <- X %*% beta + t(rmvnorm(1,mean=rep(0,nrow(X)),sigma=0.001*tau))
epsilon <- unlist(sapply(seq(1,length(Y_s)), function(i) rnorm(sapply(Y_s,length)[i], mean=0, sd=0.01)))
Y <- phi %*% a + epsilon

Y_s <- list() 
for (i in seq(1, length(unique(sample)))){
  if(design[paste0('BM',i),]==1){
    X_s[[i]] <- kronecker(diag(7),matrix(c(1,1),nrow=1)) 
  } else {
    X_s[[i]] <- kronecker(diag(7),matrix(c(1,0),nrow=1))
  }
  Y_s[[i]] <- expr2[grep('ACE2', rownames(expr2)), grep(paste0('BM', i), colnames(expr2))]
  Phi_s[[i]] <- base[names(Y_s[[i]]), ]
}
beta2 <- interative_beta(beta)
plot(beta~beta2)



