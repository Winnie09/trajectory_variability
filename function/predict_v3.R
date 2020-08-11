predict_fitting <- function(knotnum, design, cellanno, pseudotime){
  ## make the cells order according to pseudotime order
  pseudotime = pseudotime[order(pseudotime)]
  cellanno <- cellanno[match(names(pseudotime), cellanno[,1]), ]
  ## fitting
  philist <- lapply(seq(0, max(knotnum)), function(num.knot) {
    if (num.knot==0) {
      phi <- cbind(1,bs(pseudotime))
    } else {
      knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      phi <- cbind(1,bs(pseudotime,knots = knots))  
    }
  })
  names(philist) <- as.character(0:max(knotnum))
  sname <- sapply(row.names(design),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  sexpr <- sapply(names(sname),function(ss) expr[,sname[[ss]],drop=F],simplify = F)
  pred <- do.call(rbind,lapply(unique(knotnum),function(knot.num) {
    phi <- philist[[as.character(knot.num)]]
    do.call(cbind,lapply(names(sname),function(ss) {
      phiss <- phi[sname[[ss]],,drop=F]
      sexpr[[ss]] %*% (phiss %*% chol2inv(chol(crossprod(phiss)))) %*% t(phiss)
    }))
  }))
  pred <- pred[rownames(expr), colnames(expr)]
}
