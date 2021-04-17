mykmeans <- function(matrix, number.cluster = NA, maxclunum = 30, seed = 12345){
  ## cluster the rows
  library(parallel)
  if (is.na(number.cluster)){
    rss <- mclapply(1:maxclunum,function(clunum) {
      set.seed(12345)
      tmp <- kmeans(matrix,clunum,iter.max = 1000)
      tmp$betweenss/tmp$totss
    },mc.cores=30)
    rss <- unlist(rss)
    print(rss)
    # number.cluster <- which(diff(rss) < 1e-2)[1]
    x <- 2:maxclunum
    number.cluster <- x[which.min(sapply(1:length(x), function(i) {
      x2 <- pmax(0, x - i)
      sum(lm(rss[-1] ~ x + x2)$residuals^2)  ## check this
    }))]
  }
  print(number.cluster)
  set.seed(seed)
  clu <- kmeans(matrix, number.cluster)    
  return(clu)
}

