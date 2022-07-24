cc.genes <- readRDS('/home/whou10/scratch16/whou10/trajectory_variability/cellcycle/data/seurat_cc_genes.rds')
ddir <- '/home/whou10/scratch16/whou10/trajectory_variability/hca/real/build_from_tree_variability/result/'
library(splines)

for (celltype in list.files(ddir)){
  cellanno <- readRDS(paste0(ddir, celltype, '/input_cellanno.rds'))
  design <- readRDS(paste0(ddir, celltype, '/input_design.rds'))
  expr <- readRDS(paste0(ddir, celltype, '/input_expr.rds'))
  rownames(expr) <- sub(':.*', '', rownames(expr))
  pt <- readRDS(paste0(ddir, celltype, '/input_pseudotime.rds'))
  cc.s = colSums(expr[intersect(cc.genes[[1]], rownames(expr)),])
  cc.g2m = colSums(expr[intersect(cc.genes[[2]], rownames(expr)),])
  bs.all = bs(pt)
  
  pdlist <- lapply(unique(rownames(design)), function(sample){
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
    a = data.frame(cc.s.sample = cc.s.pred, 
                   cc.g2m.sample = cc.g2m.pred, 
                   pt.sample = pt.grid, 
                   sample = sample, 
                   sex = design[sample, 2])
    #plot(pt.tmp, cc.tmp, col='grey', xlab = 'pt', ylab = 'cc')
    #points(pt.grid, pred, col = 'darkgreen')
  })
  pd <- do.call(rbind, pdlist)
  saveRDS(pd, paste0('/home/whou10/scratch16/whou10/trajectory_variability/cellcycle/res/hca_', celltype, '_cellcycle_score.rds'))
}
