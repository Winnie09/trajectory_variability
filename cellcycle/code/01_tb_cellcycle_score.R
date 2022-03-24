library(Seurat)
pid = 2
pt = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
expr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds')
cellanno = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds') ## ## female = 1, male = 0

cc.s = colSums(expr[intersect(cc.genes[[1]], rownames(expr)),])
cc.g2m = colSums(expr[intersect(cc.genes[[2]], rownames(expr)),])

bs.all = bs(pt)
pdlist <- lapply(unique(cellanno[,2]), function(sample){
  ## subset the cells and pseudotime of this sample
  cell.tmp = cellanno[cellanno[,2] == sample, 1]
  pt.tmp = sort(pt[cell.tmp])
  
  ## cellcycle scores of s genes 
  cc.s.tmp = cc.s[names(pt.tmp)] 
  cc.g2m.tmp = cc.g2m[names(pt.tmp)]
  d = data.frame(cc.s.sample = cc.s.tmp, 
                 cc.g2m.sample = cc.g2m.tmp,
                 pt.sample = pt.tmp, 
                 stringsAsFactors = FALSE)
  
  ## cellcycle scores of s genes
  fit = lm(cc.s.sample ~ bs.all[cell.tmp, ], data = d)
  summary(fit)
  bs.grid = cbind(1,bs.all[round(seq(1,nrow(bs.all),length.out = 1000)),])
  cc.s.pred = (bs.grid %*% fit$coefficient)[,1]
  pt.grid = pt[round(seq(1,nrow(bs.all),length.out = 1000))]
  
  ## cellcycle scores of g2m genes
  fit = lm(cc.g2m.sample ~ bs.all[cell.tmp, ], data = d)
  cc.g2m.pred = (bs.grid %*% fit$coefficient)[,1]
  
  ## return the dataframe containing all the scores
  data.frame(cc.s.sample = cc.s.pred, 
             cc.g2m.sample = cc.g2m.pred, 
             pt.sample = pt.grid, 
             sample = sample, 
             sex = ifelse(design[sample, 2], 'female', 'male'))
  #plot(pt.tmp, cc.tmp, col='grey', xlab = 'pt', ylab = 'cc')
  #points(pt.grid, pred, col = 'darkgreen')
})
pd <- do.call(rbind, pdlist)
saveRDS(pd, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/res/tb_cellcycle_score.rds')

