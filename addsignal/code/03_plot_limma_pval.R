allfd <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/limma/')
Res <- lapply(allfd, function(fd){
  allf <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/limma/', fd), pattern = 'Para_FdrDiff_Area')
  res <- lapply(allf, function(f){
    tmp <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/limma/', fd, '/', f))
    a <- data.frame(tmp, GeneProp=sub('_.*','', f), SignalType = fd) 
  })
  res <- do.call(rbind, res)  
})
res <- do.call(rbind, Res)  
res <- res[res[,1]!=0, ]


library(ggplot2)
library(gridExtra)
ggplot(res, aes(x = Parameter, y = Area, color = GeneProp)) + geom_point()  + geom_line() + theme_classic() + facet_wrap(~SignalType)
ggplot(res, aes(x = Parameter, y = Area, color = GeneProp)) + geom_point()  + geom_line() + theme_classic() + facet_wrap(~SignalType)

