library(here)
setwd(here())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
get_ptprop <- function(d, cellanno){
  ## pseudotime density
  pseudotime = d[,4]
  names(pseudotime) = rownames(d)
  ptw <- cut(pseudotime,seq(min(pseudotime),max(pseudotime),length.out = 100),include.lowest = T)
  ptdat <- table(ptw,cellanno[match(names(pseudotime),cellanno[,1]),2])
  ptdat <- t(t(ptdat)/colSums(ptdat)) ## divided by rowsum (rowsum = 1). interval * samples. 
  ptdat <- as.data.frame(ptdat)
  colnames(ptdat) <- c('pt','s','prop')
  return(ptdat)
}

get_ccfit <- function(d, cellanno){
  ## cc score spline fitting
  cc.s = d[,1]
  cc.g2m = d[,2]
  cc.g1 = d[,3]
  pt = d[,4]
  names(pt) <- names(cc.s) <- names(cc.g2m) <- names(cc.g1) <- rownames(d)
  bs.all = splines::bs(pt)
  pdlist <- lapply(unique(cellanno[,2]), function(sample){
    ## subset the cells and pseudotime of this sample
    cell.tmp = cellanno[cellanno[,2] == sample, 1]
    pt.tmp = sort(pt[cell.tmp])
    ## cellcycle scores of s genes 
    cc.s.tmp = cc.s[names(pt.tmp)] 
    cc.g2m.tmp = cc.g2m[names(pt.tmp)]
    cc.g1.tmp = cc.g1[names(pt.tmp)]
    fitd = data.frame(cc.s.sample = cc.s.tmp, 
                      cc.g2m.sample = cc.g2m.tmp,
                      cc.g1.sample = cc.g1.tmp,
                      stringsAsFactors = FALSE)
    ## cellcycle scores of s genes
    fit = lm(cc.s.sample ~ bs.all[cell.tmp, ], data = fitd)
    summary(fit)
    bs.grid = cbind(1,bs.all[round(seq(1,nrow(bs.all),length.out = 100-1)),])
    cc.s.pred = (bs.grid %*% fit$coefficient)[,1]
    ## cellcycle scores of g2m genes
    fit = lm(cc.g2m.sample ~ bs.all[cell.tmp, ], data = fitd)
    cc.g2m.pred = (bs.grid %*% fit$coefficient)[,1]
    ## cellcycle scores of g1 genes
    fit = lm(cc.g1.sample ~ bs.all[cell.tmp, ], data = fitd)
    cc.g1.pred = (bs.grid %*% fit$coefficient)[,1]
    ## return the dataframe containing all the scores
    pd <- data.frame(cc.s.sample = cc.s.pred, 
                     cc.g2m.sample = cc.g2m.pred, 
                     cc.g1.sample = cc.g1.pred,
                     sample = sample)
  })
  pd <- do.call(rbind, pdlist)
  return(pd)
}

get_cor <- function(ptprop, ccfit){
  m <- t(sapply(unique(ptprop[,2]), function(s){
    ptprop.tmp = ptprop[ptprop[,2]==s, ]
    ccfit.tmp = ccfit[ccfit[,'sample']==s,]
    c(cor(ptprop.tmp[,3], ccfit.tmp[,1]), cor(ptprop.tmp[,3], ccfit.tmp[,2]), cor(ptprop.tmp[,3], ccfit.tmp[,3]))
  }))
  dimnames(m) = list(unique(ptprop[,2]), c('S','G2M','G1'))  
  return(m)
}

ddir <- 'hca/real/build_from_tree_variability/result/'
d <- readRDS(paste0('cellcycle/res/hca_lymph_cellcycle_score.rds'))
cellanno <- readRDS(paste0(ddir, 'lymph/input_cellanno.rds'))
ptdat = get_ptprop(d, cellanno)
ccfit = get_ccfit(d, cellanno)
mat.tb = get_cor(ptdat, ccfit)
pd = reshape2::melt(mat.tb)
colnames(pd) = c('sample', 'cellcycle', 'correlation')
pd$data = 'HCA_lymph'
pd11 <- pd

d <- readRDS(paste0('cellcycle/res/hca_erythroid_cellcycle_score.rds'))
cellanno <- readRDS(paste0(ddir, 'erythroid/input_cellanno.rds'))
ptdat = get_ptprop(d, cellanno)
ccfit = get_ccfit(d, cellanno)
mat.tb = get_cor(ptdat, ccfit)
pd = reshape2::melt(mat.tb)
colnames(pd) = c('sample','cellcycle', 'correlation')
pd$data = 'HCA_erythroid'
pd12 <- pd



d <- readRDS(paste0('cellcycle/res/hca_monocyte_cellcycle_score.rds'))
cellanno <- readRDS(paste0(ddir, 'monocyte/input_cellanno.rds'))
ptdat = get_ptprop(d, cellanno)
ccfit = get_ccfit(d, cellanno)
mat.tb = get_cor(ptdat, ccfit)
pd = reshape2::melt(mat.tb)
colnames(pd) = c('sample','cellcycle', 'correlation')
pd$data = 'HCA_monocyte'
pd13 <- pd


d  = readRDS('/home/whou10/scratch16/whou10/trajectory_variability/cellcycle/res/covid_Mod_Mi_cellcycle_score.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[rownames(d), ]
ptdat = get_ptprop(d, cellanno)
ccfit = get_ccfit(d, cellanno)
mat.tb = get_cor(ptdat, ccfit)
pd = reshape2::melt(mat.tb)
colnames(pd) = c('sample', 'cellcycle', 'correlation')
pd$data = 'COVID19'
pd2 = pd


d  = readRDS('/home/whou10/scratch16/whou10/trajectory_variability/cellcycle/res/tb_cellcycle_score.rds')
cellanno = readRDS('tb/data/sex/cellanno.rds')
ptdat = get_ptprop(d, cellanno)
ccfit = get_ccfit(d, cellanno)
mat.tb = get_cor(ptdat, ccfit)
pd = reshape2::melt(mat.tb)
colnames(pd) = c('sample', 'cellcycle', 'correlation')
pd$data = 'TB'
pd3 = pd
pd = rbind(pd11, pd12, pd13, pd2, pd3)

library(ggplot2)
pd$data = factor(pd$data, levels = c(unique(pd$data)[grepl('HCA', unique(pd$data))], setdiff(unique(pd$data), unique(pd$data)[grepl('HCA', unique(pd$data))])))

pd$grp = paste0(pd[,'data'], ';', pd[,'cellcycle'])
pd$grp = factor(pd$grp, levels = c(unique(pd$grp)[grepl('HCA', unique(pd$grp))], setdiff(unique(pd$grp), unique(pd$grp)[grepl('HCA', unique(pd$grp))])))
saveRDS(pd, '/home/whou10/scratch16/whou10/trajectory_variability/cellcycle/plot/cellcycle_pseudotime_correlation_pd.rds')

# pd <- readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/cellcycle/plot/cellcycle_pseudotime_correlation_pd.rds')

# source('/Users/wenpinhou/Dropbox/resource/ggplot_theme.R')
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('cellcycle/plot/cellcycle_pseudotime_correlation.pdf', width = 3.3, height = 2.6)
ggplot() + 
  geom_boxplot(data = pd, aes(x = data, y = correlation, fill = cellcycle), outlier.shape = NA, position=position_dodge(1), alpha = 0.3)+
  scale_fill_brewer(palette = 'Set1')+
  # geom_jitter(data = pd, aes(x = data, y = correlation, color = cellcycle),  width = 0.25, size = 0.01, stoke = 0) + 
  geom_point(data = pd, aes(x = data, y = correlation, color = cellcycle), position=position_jitterdodge(), size = 0.5, stroke = 0, alpha = 0.6) + 
  ylab('Correlation between cellcycle and density')+
  theme(axis.text.x = element_text(angle = 45))
dev.off()  
  



