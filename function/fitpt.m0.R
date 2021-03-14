# ori.design = design
# test.pattern = 'overall'
# test.position = 'all'
# maxknotallowed=10; EMmaxiter=1000; EMitercutoff=0.01; verbose=F; ncores=1; model = 3
# test.pattern = 'overall'
fitpt.m0 <- function(expr, cellanno, pseudotime, design, EMmaxiter=100, EMitercutoff=0.1, verbose=F) {
  # set.seed(12345)
  suppressMessages(library(Matrix))
  suppressMessages(library(parallel))
  suppressMessages(library(splines))
  suppressMessages(library(matrixcalc))
## This function is for model 0 (beta0 for intercept base only) fitting.
## expr: gene by cell matrix, entires and log-transformed expression
## pseudotime: a numeric vecotor of pseudotime, the names are ordered cell names
## design: design: sample by feature design matrix. rownames are sample names. first column is 1, second column is the group partition, currently only one variable.
## cellanno: dataframe, first column is cell names, second column is sample names.
## test.pattern: c('slope', 'intercept', 'overall).
pseudotime <- pseudotime[colnames(expr)]
cellanno <- cellanno[match(colnames(expr),cellanno[,1]),]
sname <- sapply(unique(cellanno[,2]),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
as <- names(sname)
cn <- sapply(sname,length)

## initialize
B <- rowMeans(expr)

indfit <- sapply(as,function(s) {
  rowMeans(expr[, cellanno[,2]==s, drop=F])
},simplify = F)

s2 <- matrix(sapply(as,function(s) {
  tmp <- expr[, cellanno[,2]==s, drop=F]-indfit[[s]]
  n <- ncol(tmp)
  (rowMeans(tmp*tmp)-(rowMeans(tmp))^2)*(n)/(n-1)
}),nrow=nrow(expr),dimnames=list(rownames(expr),as))
alpha <- 1/apply(s2,1,var)*rowMeans(s2)^2+2
eta <- (alpha-1)*rowMeans(s2)

diffindfit <- sapply(as,function(s) {
  indfit[[s]] - B
})

omega <- apply(diffindfit,1,var)

print('start EM ......')
iter <- 0
EMitercutoff <- 0
gidr <- rownames(expr)
all <- matrix(-Inf, nrow=nrow(expr),ncol=1,dimnames = list(rownames=gidr))
etalist <- alphalist <- omegalist <- Nlist <- Jslist <- list()
while (iter < EMmaxiter && length(gidr) > 0) {
  print(paste0('iter ', iter))
  print(gidr)
  expr_phibx <- sapply(as,function(s) {
    expr[, cellanno[,2]==s, drop=F][gidr,,drop=F]-B[gidr]
  },simplify = F)
 
  L <- sapply(as,function(s) {
    rowSums(expr_phibx[[s]] * expr_phibx[[s]])
  },simplify = F)
 
  Jsolve <- matrix(sapply(as,function(s) {
    1/(cn[s] + 1/omega[gidr])
  }), nrow = length(gidr), dimnames = list(gidr, as))
 
  K <- sapply(as,function(s) {
    rowSums(expr_phibx[[s]])
  },simplify = F)
 
  JK <- sapply(as,function(s) {
    Jsolve[,s,drop=F] * K[[s]]  ## debug here !!
  },simplify = F)
 
  L2eKJK <- sapply(as,function(s) {
    2*eta[gidr] + L[[s]] - K[[s]] * JK[[s]]
  },simplify = F)
 
  A <- matrix(sapply(as,function(s) {
    log(L2eKJK[[s]]/2)-digamma(alpha[gidr]+cn[s]/2)
  }),nrow=length(gidr),dimnames=list(gidr,as))
 
  N <- matrix(sapply(as,function(s) {
    (2*alpha[gidr]+cn[s])/L2eKJK[[s]]
  }),nrow=length(gidr),dimnames=list(gidr,as))
 
  ll <- rowSums(matrix(sapply(as,function(s) {
    dv <- omega[gidr]/Jsolve[gidr,s,drop=F]
    alpha[gidr]*log(2*eta[gidr])+lgamma(cn[s]/2+alpha[gidr])-cn[s]*log(pi)/2-lgamma(alpha[gidr])-log(dv)/2-(cn[s]/2+alpha[gidr])*log(L2eKJK[[s]])
  }),nrow=length(gidr),dimnames = list(rownames=gidr,colnames=as)))
 
  ## -------------->
  B1 <- rowSums(sapply(as, function(s){
    N[,s] * cn[s]
  }))
  B2 <- rowSums(sapply(as, function(s){
    N[,s] * colSums(t(expr[ ,cellanno[,2]==s, drop=F][gidr,]) - rep(JK[[s]],each=cn[s]))
  }))
  B[gidr] <- B2/B1
 
  ## M -step:
  omega[gidr] <- rowMeans(sapply(as,function(s) {
    Jsolve[,s,drop=F] + N[,s]*JK[[s]] * JK[[s]]  ## debug here !!!
  }))
 
  eta[gidr] <- sapply(gidr,function(g) {
    meanN <- mean(N[g,])
    meanA <- mean(A[g,])
    uniroot(function(eta) {digamma(eta * meanN)-log(eta)+meanA},c(1e-10,1e10))$root
  })
  alpha[gidr] <- eta[gidr] * rowMeans(N)
 
  iter <- iter + 1
 
  llv <- all[,ncol(all)]
  llv[gidr] <- ll
  all <- cbind(all,llv)
  gidr <- names(which(all[,ncol(all)] - all[,ncol(all)-1] > EMitercutoff))
  etalist[[iter]] <- eta
  alphalist[[iter]] <-alpha
  omegalist[[iter]] <-omega
  Nlist[[iter]] <- N
  Jslist[[iter]] <- Jsolve
  rm(list = c('L','Jsolve', 'K'))
}
  
  print('list parameter ......')
  allres  <- list(beta = B, alpha = alpha, eta = eta, omega = omega, logL = all)
  
  para <- list()
  
  for (j in rownames(expr)) {
    para[[j]] <- list(beta=allres[[1]][j],
                      alpha=allres[[2]][j],
                      eta=allres[[3]][j],
                      omega=allres[[4]][j],
                      ll=allres[[5]][j,ncol(allres[[5]])])
  }
  
  return(list(parameter=para))
}



