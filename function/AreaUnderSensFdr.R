AreaUnderSensFdr <- function(sensfdr){
  tmp <- sensfdr
  bound <- approx(x=tmp[,3],y=tmp[,2],xout=0.25)$y
  tmp <- rbind(tmp[tmp[,3] < 0.25,2:3],c(bound,0.25))
  tmp <- unique(tmp)
  diff <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm = T)-0.25*0.25/2   ## (area under Real_FDR ~ Reported_FDR)-0.25*0.25/2
  tmp <- sensfdr
  bound <- approx(x=tmp[,2],y=tmp[,1],xout=0.25)$y
  tmp <- rbind(tmp[tmp[,2] < 0.25,1:2],c(bound,0.25))
  area <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm=T)/0.25
  res <- c(diff, area)
  names(res) = c('Fdr.Diff','Area')
  res
}

