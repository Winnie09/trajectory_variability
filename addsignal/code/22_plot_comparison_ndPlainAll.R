allfd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/ndPlainAll/tradeSeq/')
final <- lapply(allfd, function(fd){
  allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/ndPlainAll/tradeSeq/', fd))
  allf <- allf[!grepl('sce', allf)]
  if (length(allf)>=1){
    Res <- lapply(allf, function(f){
      r <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/ndPlainAll/tradeSeq/', fd, '/', f))
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
Res1 = final[!is.na(final[,1]), ]    

allfd <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/ndPlainAll/EM/')
final <- lapply(allfd, function(fd){
  allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/ndPlainAll/EM/', fd))
  allf <- allf[!grepl('fdr', allf)]
  if (length(allf)>=1){
    Res <- t(sapply(allf, function(f){
      r <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/ndPlainAll/EM/', fd, '/', f))
      if (ncol(r$perll) == 70) {
        c(Method = 'EM', r[['sensfdr']][-1], SignalType = fd, GeneProp = sub('_.*','', f), Parameter = sub('.rds','',sub('.*_','',f)))  
      } else {
        c(Method = '', 'Fdr.Diff', 'Area', SignalType = '', GeneProp = '', Parameter = '')
      }
      
    }))
  } else return(NA)
})
final <- do.call(rbind, final)
final <- as.data.frame(final)
final[!colnames(final)%in%c('Method','SignalType','GeneProp')] <- apply(final[!colnames(final)%in%c('Method','SignalType','GeneProp')], 2, as.numeric)
Res2 = final[!is.na(final[,1]), ]

pd <- rbind(Res1, Res2)
pd <- pd[complete.cases(pd), ]
library(ggplot2)
library(gridExtra)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/plot/ndPlainAll/compare_fdr_diff.pdf',width=8,height=4)
ggplot(pd, aes(x = Parameter, y = Fdr.Diff, shape = GeneProp, color=Method)) + geom_point(size=1)  + geom_line(size=0.1) + theme_classic() + facet_wrap(~SignalType) + ylab('fdr.diff(real~reported - 0.25*0.25/2)') 
dev.off()
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/plot/ndPlainAllcompare_auc.pdf',width=8,height=4)
ggplot(pd, aes(x = Parameter, y = Area, shape = GeneProp, color=Method)) + geom_point(size=1)  + geom_line(size=0.1) + theme_classic() + facet_wrap(~SignalType)
dev.off()

### debug code
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
fd = 'linear'
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/ndPlainAll/EM/', fd))
  allf <- allf[!grepl('fdr', allf)]
f <- '0.05_2.4.rds'
r <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/result/ndPlainAll/EM/', fd, '/', f))
res <- r[['res']]
data <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/data/ndPlainAll/1256/linear/', f))
selgene <- data$selgene


sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
plot(sensfdr[,2]~sensfdr[,3], xlab='Reported Fdr', ylab='Real Fdr', main=paste0('linear signal, ',sub('.rds','',f), ', fdr.diff=',round(as.numeric(r[[2]][2]),3)), pch=20)
abline(0, 1, col='red')


plot(sensfdr[,1]~sensfdr[,3], xlab='Reported Fdr', ylab='Sensitivity', main=paste0('linear signal, ',sub('.rds','',f), ', AUC=',round(as.numeric(r[[2]][3]),3)), pch=20)
abline(0, 1, col='red')

ggplot(data.frame(x=data$pseudotime,y=data$expr[data$selgene[1],],p=sub(':.*','',colnames(data$expr))),aes(x=x,y=y)) + geom_point() + facet_wrap(~p)

