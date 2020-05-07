allfd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/tradeSeq/')
final <- lapply(allfd, function(fd){
  allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/tradeSeq/', fd))
  allf <- allf[!grepl('sce', allf)]
  if (length(allf)>=1){
    Res <- lapply(allf, function(f){
      r <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/tradeSeq/', fd, '/', f))
      res <- t(sapply(seq(1,length(r)), function(i){
        v <- r[[i]][['sensfdr']]
        v[1] <- paste0('tradeSeq_',names(r)[i])
        v
      }))
      colnames(res)[1] <- 'Method'  
      res = data.frame(res, SignalType = fd, GeneProp = sub('_.*','', f), Parameter = gsub('.rds','',sub('.*_','', f)))
    })
    Res <- do.call(rbind, Res)
    Res
  } else return(NA)
})
final <- do.call(rbind, final)
final[!colnames(final)%in%c('Method','SignalType','GeneProp')] <- apply(final[!colnames(final)%in%c('Method','SignalType','GeneProp')], 2, as.numeric)
final = final[!is.na(final[,1]), ]    


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
res <- data.frame(res, Method = 'Ours')

pd = rbind(final, res[, colnames(final)])
library(ggplot2)
library(gridExtra)
p1 <- ggplot(pd, aes(x = Parameter, y = Fdr.Diff, shape = GeneProp, color=Method)) + geom_point()  + geom_line() + theme_classic() + facet_wrap(~SignalType)
p2 <- ggplot(pd, aes(x = Parameter, y = Area, shape = GeneProp, color=Method)) + geom_point()  + geom_line() + theme_classic() + facet_wrap(~SignalType)
grid.arrange(p1,p2, nrow=1)

