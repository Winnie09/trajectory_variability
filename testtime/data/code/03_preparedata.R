geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
# geneProp <- 0.05
# addSignalType <- 'linear'
# addSignalPara <-  0.5

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testtime/data/data/'
if (!file.exists(paste0(rdir, addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))){
  source('./function/01_function.R')
  source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
  expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
  ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
  id <- which(ct %in% c('HSC','MEP','Ery'))
  expr <- expr[,id]
  expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
  allp = sub(':.*','', colnames(expr))
  sample <- sub(':.*','',colnames(expr))
  names(sample) <- colnames(expr)
  design <- cbind(c(1,1,0,0,1,1,0,0))
  row.names(design) <- paste0('BM',1:8)
  colnames(design) <- 'group'
  
  ### reorder cells for each patient
  pt <- readRDS(paste0(rdir,'null/pseudotime.rds'))
  pt <- data.frame(pt, sample = sub(':.*','', pt$cell), stringsAsFactors = F)
  ptlist <- lapply(unique(pt[,3]), function(p){
    tmp <- pt[as.character(pt[,3])==p, ]
    set.seed(12345)
    tmp[, 'pseudotime'] <- sample(tmp[, 'pseudotime'])
    tmp
  })
  pt <- do.call(rbind, ptlist)
  
  ### add signal
  set.seed(12345)
  selgene <- sample(row.names(expr), round(geneProp * nrow(expr)))
  celln <- sapply(colnames(expr), function(i) paste0(strsplit(i,':')[[1]][2],':',strsplit(i,':')[[1]][3]))
  expr2 <- AddSignal(expr = expr,SelectGene = selgene, pseudotime = pt, method = addSignalType, parameter=addSignalPara, type = 'all')
  # plot(expr[g,pt[,1]]~pt[,2], pch=20, cex=.1)
  # points(expr3[g,pt[,1]]~pt[,2], pch=20, col='red', cex=.1)

  dir.create(paste0(rdir, addSignalType), recursive = TRUE)
  saveRDS(list(expr = expr2, selgene = selgene, pseudotime = pt, design = design), paste0(rdir, addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))
}

