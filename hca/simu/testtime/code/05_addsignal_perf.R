library(here)
here()
source(here('function', '01_function.R'))
selgene <- readRDS(here('hca', 'data', 'simu', 'testtime', 'addMultiSignalUsingExpr','selgene','selgene.rds'))
rdir <- here('hca', 'simu','testtime','result','addsignal')

## EM_centered
perflist <- list()
af = list.files(paste0(rdir, '/EM_centered/'))
f = af[1]
if (length(af) > 0){
  df1<- t(sapply(af, function(f){
    print(f)
    Res = readRDS(paste0(rdir, '/EM_centered/', f))
    res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange)
    res = res[order(res[,1], -res[,2]), ]
    c(sub('.rds','',f), 'EM_centered', AreaUnderSensFdr(SensFdr(selgene, res)))
  }))
  perflist[['EM_centered']] <- df1
}
  
## EM_NOT_centered
rdir2 <- here('hca','simu','testtime','result','addsignal','EM_NOT_centered')
af = list.files(rdir2)
if (length(af) > 0){
  df6<- t(sapply(af, function(f){
    print(f)
    Res = readRDS(paste0(rdir2, '/',f))
    res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange)
    res = res[order(res[,1], -res[,2]), ]
    c(sub('.rds','',f), 'EM_NOT_centered', AreaUnderSensFdr(SensFdr(selgene, res)))
  }))
  perflist[['EM_NOT_centered']] <- df6
}
  
## tradeSeq
m = 'tradeSeq'
af = list.files(paste0(rdir, '/',m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('sce', af)]
if (length(af) > 0){
  df <- lapply(af, function(f){
    print(f)
    r = readRDS(paste0(rdir,'/', m, '/', f))
    tmp <- r[['startVsEndTest']]
    tmp <- tmp[complete.cases(tmp), ]
    a = c(sub('.rds','',f), 'tradeSeq_startVsEndTest', AreaUnderSensFdr(SensFdr(selgene, tmp)))
    tmp <- r[['associationTest']]
    tmp <- tmp[complete.cases(tmp), ]
    b = c(sub('.rds','',f), 'tradeSeq_associationTest', AreaUnderSensFdr(SensFdr(selgene, tmp)))
    rbind(a, b)
  })
  df2 <- do.call(rbind, df)
  perflist[['tradeSeq']] <- df2
}
  
## tscan
m = 'tscan'
af = list.files(paste0(rdir, '/',m, '/'))
if (length(af) > 0){
  df3 <- t(sapply(af, function(f){
    print(f)
    r = readRDS(paste0(rdir, '/', m, '/', f))
    c(sub('.rds','',f),  'tscan', AreaUnderSensFdr(SensFdr(selgene, r)))
  }))
  perflist[['tscan']] <- df3
}
  

## monocle2
m = 'monocle2'
af = list.files(paste0(rdir, '/', m, '/'))
if (length(af) > 0){
  df4 <- t(sapply(af, function(f){
    print(f)
    r = readRDS(paste0(rdir, '/', m, '/', f))
    c(sub('.rds','',f),  'monocle2', AreaUnderSensFdr(SensFdr(selgene, r)))
  }))
  perflist[['monocle2']] <- df4
}
  

## monocle3
m = 'monocle3'
af = list.files(paste0(rdir, '/', m, '/'))
if (length(af) > 0){
  df5 <- t(sapply(af, function(f){
    print(f)
    r = readRDS(paste0(rdir, '/', m, '/', f))
    c(sub('.rds','',f),  'monocle3', AreaUnderSensFdr(SensFdr(selgene, r)))
  }))
  perflist[['monocle3']] <- df5
}
    

## concatenate
perf <- do.call(rbind, perflist)
colnames(perf) <- c('Type', 'Method', 'Fdr.Diff', 'AUC')
saveRDS(perf, here('hca', 'simu','testtime','result','addsignal','perf','perf.rds'))
rm(list=ls())


