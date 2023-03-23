library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')

### load saver, count matrix, and pseudotime
saverlog <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','hsc_mep_ery_saver.rds'))
cnt <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','hsc_mep_ery_count.rds'))
pt <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','pseudotime.rds'))
saverlog <- saverlog[, pt[,1]]
cnt <- cnt[, pt[,1]]

intgene <- intersect(rownames(cnt),rownames(saverlog))
cnt <- cnt[intgene,]
saverlog <- saverlog[intgene,]

# ### permute the sample-cell relationsihp 
# identical(colnames(saverlog), colnames(cnt))
# sn <- sub(':.*', '', colnames(cnt))
# cn <- sapply(1:ncol(cnt), function(i) sub(paste0(sn[i], ':'), '', colnames(cnt)[i]))
# plot(saverlog[grep('XIST', rownames(cnt)),], col = as.factor(sn %in% c(paste0('BM', c(1,2,5,6)))))
# plot(saverlog[grep('XIST', rownames(cnt)),], col = as.factor(sn))
seed <- as.numeric(commandArgs(trailingOnly = T))
set.seed(seed)
pt.v <- pt[,2]
names(pt.v) <- pt[,1]

cellanno_pm <- data.frame(cell = colnames(cnt), sample = sub(':.*', '', colnames(cnt)), stringsAsFactors = F)

design <- matrix(c(rep(1,8), 1,1,0,0,1,1,0,0), nrow = 8)
dimnames(design) <- list(paste0('BM', 1:8), c('intercept', 'group'))
g1g <- rownames(design)[design[,2]==1]
g0g <- rownames(design)[design[,2]==0]
window <- cut(pt.v,quantile(pt.v,seq(0,1,length.out=101)),include.lowest=T)
names(window) <- names(pt.v)


for (sw in unique(window)) {
  print(sw)
  sid <- names(window)[which(window==sw)]
  sample <- sub(':.*','',sid)
  ssaver <- saverlog[,sid]
  mm <- sapply(unique(sample),function(us) apply(ssaver[,sample==us,drop=F],1,median))
  g1m <- apply(mm[,intersect(colnames(mm),g1g)],1,median)
  g0m <- apply(mm[,intersect(colnames(mm),g0g)],1,median)
  g1lowgene <- names(which(g1m < g0m))
  g0lowgene <- names(which(g1m > g0m))
  
  g1lowgeneg1targetcell <- sapply(g1lowgene,function(g) names(which.min(abs(ssaver[g,sid[sample %in% g1g]]-g1m[g]))))
  g1lowgeneg0targetcell <- sapply(g1lowgene,function(g) names(which.min(abs(ssaver[g,sid[sample %in% g0g]]-g0m[g]))))
  g0lowgeneg0targetcell <- sapply(g0lowgene,function(g) names(which.min(abs(ssaver[g,sid[sample %in% g0g]]-g0m[g]))))
  g0lowgeneg1targetcell <- sapply(g0lowgene,function(g) names(which.min(abs(ssaver[g,sid[sample %in% g1g]]-g1m[g]))))
  
  # summary(abs(g1m[g1lowgene]-sapply(g1lowgene,function(g) ssaver[g,g1lowgeneg1targetcell[g]])))
  # summary(abs(g1m[g0lowgene]-sapply(g0lowgene,function(g) ssaver[g,g0lowgeneg1targetcell[g]])))
  # summary(abs(g0m[g1lowgene]-sapply(g1lowgene,function(g) ssaver[g,g1lowgeneg0targetcell[g]])))
  # summary(abs(g0m[g0lowgene]-sapply(g0lowgene,function(g) ssaver[g,g0lowgeneg0targetcell[g]])))
  
  saverlog[g1lowgene,sid[sample %in% g1g]] <- saverlog[g1lowgene,sid[sample %in% g1g]]-g1m[g1lowgene] + g0m[g1lowgene]
  saverlog[g0lowgene,sid[sample %in% g0g]] <- saverlog[g0lowgene,sid[sample %in% g0g]]-g0m[g0lowgene] + g1m[g0lowgene]
  
  cnt[g1lowgene,sid[sample %in% g1g]] <- cnt[g1lowgene,sid[sample %in% g1g]]-sapply(g1lowgene,function(g) cnt[g,g1lowgeneg1targetcell[g]]) + sapply(g1lowgene,function(g) cnt[g,g1lowgeneg0targetcell[g]])
  cnt[g0lowgene,sid[sample %in% g0g]] <- cnt[g0lowgene,sid[sample %in% g0g]]-sapply(g0lowgene,function(g) cnt[g,g0lowgeneg0targetcell[g]]) + sapply(g0lowgene,function(g) cnt[g,g0lowgeneg1targetcell[g]])
}

ptper <- pt.v
for (sw in unique(window)) {
  tmp <- pt.v[window==sw]
  samp <- sub(':.*','',names(tmp))
  for (u in unique(samp)) 
    if (sum(samp==u) > 1)
      tmp[samp==u] <- sample(tmp[samp==u])
  ptper[window==sw] <- tmp
}

ptper <- sort(ptper)
cellanno_pm <- cellanno_pm[match(names(ptper),cellanno_pm[,1]),]
saverlog <- saverlog[,names(ptper)]
cnt <- cnt[,names(ptper)]
savepath <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/repdata/data/',seed)
dir.create(savepath)
saveRDS(cellanno_pm, paste0(savepath,'cellanno_pm.rds'))
saveRDS(design, paste0(savepath,'design.rds'))
saveRDS(ptper, paste0(savepath,'pseudotime_pm.rds'))
saveRDS(saverlog, paste0(savepath,'saverlog_pm.rds'))
saveRDS(cnt, paste0(savepath,'cnt_pm.rds'))


# samp <- sub(':.*','',colnames(saverlog))
# mmm <- sapply(unique(samp),function(us) apply(saverlog[,samp==us],1,median))
# 
# ggplot(data.frame(e=saverlog['HBA2:ENSG00000188536',],pt=pt.v,samp=samp,group=samp %in% g1g),aes(x=pt,y=e,col=group,group=samp)) + geom_smooth()

