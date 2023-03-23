cc.genes <- readRDS('/home/whou10/scratch16/whou10/trajectory_variability/cellcycle/data/seurat_cc_genes.rds')
ddir <- '/home/whou10/scratch16/whou10/trajectory_variability/hca/real/build_from_tree_variability/result/'
library(splines)
cc.genes[['g1.genes']] <- c('CDC23',	'CDC25C',	'CDC6',	'CDK10',	'CDK2',	'CDK6',	'CDKN1C',	'E2F1',	'FOXO4',	'GFI1B',	'MAP3K11',	'PRUNE2',	'RB1',	'TAF1',	'TBRG4')

for (celltype in list.files(ddir)){
  print(celltype)
  cellanno <- readRDS(paste0(ddir, celltype, '/input_cellanno.rds'))
  design <- readRDS(paste0(ddir, celltype, '/input_design.rds'))
  expr <- readRDS(paste0(ddir, celltype, '/input_expr.rds'))
  rownames(expr) <- sub(':.*', '', rownames(expr))
  
  pt <- readRDS(paste0(ddir, celltype, '/input_pseudotime.rds'))
  cc.s = colMeans(expr[intersect(cc.genes[[1]], rownames(expr)),])
  cc.g2m = colMeans(expr[intersect(cc.genes[[2]], rownames(expr)),])
  cc.g1 = colMeans(expr[intersect(cc.genes[[3]], rownames(expr)),])
  bs.all = bs(pt)
  
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
                   sex = design[sample, 2], 
                   stringsAsFactors = FALSE)
    
  })
  pd <- do.call(rbind, pdlist)
  saveRDS(pd, paste0('/home/whou10/scratch16/whou10/trajectory_variability/cellcycle/res/hca_', celltype, '_cellcycle_score.rds'))
}


