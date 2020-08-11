predict_fitting_v2 <- function(Knotnum, Design, Cellanno, Pseudotime, Resolution = 1000){
  philist <- lapply(seq(0, max(Knotnum)), function(num.knot) {
    if (num.knot==0) {
      phi <- cbind(1,bs(Pseudotime))
    } else {
      knots = seq(min(Pseudotime),max(Pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      phi <- cbind(1,bs(Pseudotime,knots = knots))  
    }
  })
  names(philist) <- as.character(0:max(Knotnum))
  sname <- sapply(row.names(Design),function(i) Cellanno[Cellanno[,2]==i,1],simplify = F)
  sexpr <- sapply(names(sname),function(ss) expr[,sname[[ss]],drop=F],simplify = F)
  pred <- do.call(rbind,lapply(unique(Knotnum),function(knot.num) {
    phi <- philist[[as.character(knot.num)]]
    do.call(cbind,lapply(names(sname),function(ss) {
      phiss <- phi[sname[[ss]],,drop=F]
      sexpr[[ss]] %*% (phiss %*% chol2inv(chol(crossprod(phiss)))) %*% t(phiss)
    }))
  }))
  pred <- pred[rownames(expr), colnames(expr)]
}
