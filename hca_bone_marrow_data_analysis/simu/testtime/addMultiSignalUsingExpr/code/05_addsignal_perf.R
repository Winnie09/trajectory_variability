library(here)
here()
source(here('function', '01_function.R'))
selgene <- readRDS(here('hca','simu', 'testtime', 'addMultiSignalUsingExpr','data','selgene','selgene.rds'))
rdir <- here('hca', 'simu','testtime','addMultiSignalUsingExpr','result','addsignal')
perflist <- list()

## EM_NOT_centered
rdir2 <- here('hca','simu','testtime','addMultiSignalUsingExpr','result','addsignal','EM_pm')
af = list.files(rdir2)
if (length(af) > 0){
  df6<- sapply(af, function(f){
    print(f)
    Res = readRDS(paste0(rdir2, '/',f))
    stat = Res$statistics
    res = data.frame(fdr = stat$fdr.overall, zscore = stat$z.overall)
    rownames(res) = rownames(stat)
    res = res[order(res[,1], -res[,2]), ]
    c(sub('.rds','',f), 'Lamian.pm', AreaUnderSensFdr(SensFdr(selgene, res)))
  }, simplify = FALSE)
  df6 <- do.call(rbind, df6)
  perflist[['Lamian.pm']] <- df6
}

## EM chisq
m = 'EM_chisq'
af = list.files(paste0(rdir, '/',m, '/'))
if (length(af) > 0){
  df <- sapply(af, function(f){
    print(f)
    Res = readRDS(paste0(rdir, '/', m, '/', f))
    stat = Res$statistics
    res = data.frame(fdr = stat[,1], statistics = stat$llr/stat$df.diff)
    rownames(res) = rownames(stat)

    res = res[order(res[,1], -res[,2]), ]
    c(sub('.rds','',f),  'Lamian.chisq', AreaUnderSensFdr(SensFdr(selgene, res)))
  }, simplify = FALSE)
  df2 <- do.call(rbind, df)
  perflist[['Lamian.chisq']] <- df2
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
  df3 <- sapply(af, function(f){
    print(f)
    r = readRDS(paste0(rdir, '/', m, '/', f))
    c(sub('.rds','',f),  'tscan', AreaUnderSensFdr(SensFdr(selgene, r)))
  }, simplify = FALSE)
  df3 <- do.call(rbind, df3)
  perflist[['tscan']] <- df3
}

## monocle2
m = 'monocle2'
af = list.files(paste0(rdir, '/', m, '/'))
if (length(af) > 0){
  df4 <- sapply(af, function(f){
    print(f)
    r = readRDS(paste0(rdir, '/', m, '/', f))
    c(sub('.rds','',f),  'monocle2', AreaUnderSensFdr(SensFdr(selgene, r)))
  }, simplify = FALSE)
  df4 <- do.call(rbind, df4)
  perflist[['monocle2']] <- df4
}

## monocle3
m = 'monocle3'
af = list.files(paste0(rdir, '/', m, '/'))
if (length(af) > 0){
  df5 <- sapply(af, function(f){
    print(f)
    r = readRDS(paste0(rdir, '/', m, '/', f))
    c(sub('.rds','',f),  'monocle3', AreaUnderSensFdr(SensFdr(selgene, r)))
  }, simplify = FALSE)
  df5 <- do.call(rbind, df5)
  perflist[['monocle3']] <- df5
}

## concatenate
perf <- do.call(rbind, perflist)
colnames(perf) <- c('Type', 'Method', 'Fdr.Diff', 'AUC')
perf <- as.data.frame(perf)
perf$SignalStreghth <- as.numeric(perf$Type)
perf[,3] = as.numeric(as.character(perf[,3]))
perf[,4] = as.numeric(as.character(perf[,4]))
saveRDS(perf, here('hca', 'simu','testtime','addMultiSignalUsingExpr','result','addsignal','perf','perf.rds'))
rm(list=ls())





