library(Seurat)
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
#meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds') ## Mod = 0, Mi = 1

cc.genes[['g1.genes']] <- c('CDC23',	'CDC25C',	'CDC6',	'CDK10',	'CDK2',	'CDK6',	'CDKN1C',	'E2F1',	'FOXO4',	'GFI1B',	'MAP3K11',	'PRUNE2',	'RB1',	'TAF1',	'TBRG4')

cc.s = colMeans(expr[intersect(cc.genes[[1]], rownames(expr)),])
cc.g2m = colMeans(expr[intersect(cc.genes[[2]], rownames(expr)),])
cc.g1 = colMeans(expr[intersect(cc.genes[[3]], rownames(expr)),])

pdlist <- lapply(unique(rownames(design)), function(sample){
  ## subset the cells and pseudotime of this sample
  cell.tmp = cellanno[cellanno[,2] == sample, 1]
  pt.tmp = sort(pt[cell.tmp])
  
  ## cellcycle scores of s genes 
  cc.s.tmp = cc.s[names(pt.tmp)] 
  cc.g2m.tmp = cc.g2m[names(pt.tmp)]
  cc.g1.tmp = cc.g1[names(pt.tmp)]

  d = data.frame(cc.s.sample = cc.s.tmp, 
                 cc.g2m.sample = cc.g2m.tmp,
                 cc.g1.sample = cc.g1.tmp,
                 pt.sample = pt.tmp, 
                 sample = sample, 
                 severity = ifelse(design[sample, 2], 'Mi', 'Mod'),
                 stringsAsFactors = FALSE)
})
pd <- do.call(rbind, pdlist)
saveRDS(pd, 'cellcycle/res/covid_Mod_Mi_cellcycle_score.rds')


