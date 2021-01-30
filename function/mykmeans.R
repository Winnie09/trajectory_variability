mykmeans <- function(matrix, number.cluster = NA, maxclunum = 50, seed = 12345){
  ## cluster the rows
  set.seed(seed)
  library(parallel)
  if (is.na(number.cluster)){
    rss <- mclapply(1:maxclunum,function(clunum) {
      set.seed(12345)
      tmp <- kmeans(matrix,clunum,iter.max = 1000)
      tmp$betweenss/tmp$totss
    },mc.cores=20)
    rss <- unlist(rss)
    x <- 1:maxclunum
    optclunum <- which.min(sapply(1:maxclunum, function(i) {
        x2 <- pmax(0, x - i)
        sum(lm(rss ~ x + x2)$residuals^2)  ## check this
    }))
    print(optclunum)
    clu <- kmeans(matrix,optclunum)
  } else {
    clu <- kmeans(matrix, number.cluster)    
  }
    return(clu)
}
