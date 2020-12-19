library(here)
here()
ddir <- here('hca','simu','testtime','result','addsignal')
pdir <- here('hca','simu','testtime','plot','addsignal', 'exampleGene')
fit <- readRDS(paste0(ddir,'/cluster','/4','/population_fit.rds'))
clu <- readRDS(paste0(ddir, '/cluster','/4','/cluster.rds'))
source(here('function/01_function.R'))

tradeseq <- readRDS(paste0(ddir, '/tradeSeq/2.rds'))
  tradeseq1 <- tradeseq[[1]]
  tradeseq1 <- rownames(tradeseq1[tradeseq1[,3] > 0.05,])
  tradeseq2 <- tradeseq[[2]]
  tradeseq2 <- rownames(tradeseq2[tradeseq2[,3] > 0.05,])
  
  tscan <- readRDS(paste0(ddir, '/tscan/2.rds'))
  tscan <- rownames(tscan)[tscan[,3] > 0.05]
  
  monocle2 <- readRDS(paste0(ddir, '/monocle2/2.rds'))
  monocle2 <- rownames(monocle2)[monocle2[,3] > 0.05]
  
  monocle3 <- readRDS(paste0(ddir, '/monocle3/2.rds'))
  monocle3 <- rownames(monocle3)[monocle3[,5] > 0.05]
  
  tab <- table(c(tradeseq1, tradeseq2, tscan, monocle2, monocle3))
  gl <- names(tab)[which(tab == 5)]
  
  selgene <- readRDS(here('hca/data/simu/testtime/addMultiSignalUsingExpr/selgene/selgene.rds'))
  gl <- intersect(gl, selgene) ## other methods cannot detect
  
for (clu.select in 1:max(clu)){
  print(clu.select)
  
  gl.select <- names(clu[clu == clu.select])
  
  res <- readRDS(paste0(ddir, '/EM_NOT_centered/2.rds'))
  fdr <- res$fdr
  gl.select <- names(sort(fdr[gl.select][fdr[gl.select] < 0.05]))
  print(gl.select[1])
  
  pdf(paste0(pdir, '/exampleGene_cluster', clu.select, '_sampleLevel.pdf' ), width = 3.8, height = 3)
  plotGene(res, gl.select[1], plot.point = TRUE, line.alpha = 1, line.size = 2, point.alpha=0.5, point.size=0.5)
  dev.off()
  pdf(paste0(pdir, '/exampleGene_cluster', clu.select, '_populationLevel.pdf' ), width = 3.8, height = 3)
  plotGenePopulation(res, gl.select[1])
  dev.off()
}

  
