findbranch <- function(mst, order, origin){
  deg <- degree(mst)
  vertex <- names(deg[which(deg > 2 | deg == 1)])
  if (!origin %in% vertex) vertex <- c(origin, vertex)
  eg <- expand.grid(1:length(vertex), 1:length(vertex))
  eg <- eg[eg[,1]<eg[,2],]
  eg = data.frame(vertex[eg[,1]], vertex[eg[,2]], stringsAsFactors = FALSE)
  library(igraph)
  tmpbranch <- lapply(seq(1,nrow(eg)), function(i){
    sp <- shortest_paths(mst, from = eg[i,1], to = eg[i,2])$vpath[[1]]
    if (sum(vertex %in% sp) == 2) as.vector(sp)
  })
  tmpbranch <- tmpbranch[sapply(tmpbranch, length) >0]  
 
  allbranch <- gsub('backbone ', '', gsub('branch: ', '', names(order)))
  allbranch <- sapply(allbranch, function(i) strsplit(i, ',')[[1]])
  allbranch <- paste0(names(allbranch), collapse = ' ')
  newbranch <-sapply(tmpbranch, function(i) {
      tmp <- paste0(i, collapse = ',')
      if (!grepl(tmp, allbranch)){
        rev(i)
      } else {
        i
      }
  })
  return(newbranch)
}
