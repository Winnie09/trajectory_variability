library(parallel)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
plotdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/test_type_position/plot/'
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/data/data/count/clusterType9_4.rds')
rownames(d) <- sub(':.*','',rownames(d))
m = log2(d + 1)
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/null/pseudotime.rds')
pseudotime = pt[,2]
names(pseudotime) = pt[,1]
ap <- sub(':.*', '', colnames(m))
design = cbind(1, c(1,1,0,0,1,1,0,0))
rownames(design) = paste0('BM', seq(1,8))
colnames(design) <- c('intersect', 'condition')
ca <- data.frame(Cell = colnames(m), Sample = ap)

##  test and plot
tmp <- mclapply(c('slope', 'all'), function(type){
  mclapply(c('all', 'start', 'middle', 'end'), function(pos){
      res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = type, test.position = pos)
      res1 <- res
      gene = names(sort(res$fdr))
      pdf(paste0(plotdir, type, '_', pos, '_genes.pdf'), widh = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      pdf(paste0(plotdir, type, '_', pos, 'l_fdr_foldchange.pdf'), widh = 4, height = 4)
      print(plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange'))
      dev.off()
      return(0)
  })
})
rm(list=ls())
# res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = 'slope', test.position = 'all')
# res1 <- res
# gene = names(sort(res$fdr))
# pdf(paste0(plotdir, 'slope_all_genes.pdf'), widh = 8, height = 8)
# plotGene(testptObj = res, gene = gene[1:16], variable = 'condition')
# dev.off()
# pdf(paste0(plotdir, 'slope_all_fdr_foldchange.pdf'), widh = 4, height = 4)
# plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange')
# dev.off()
# ##  2
# res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = 'slope', test.position = 'start')
# res2 <- res
# gene = names(sort(res$fdr))
# pdf(paste0(plotdir, 'slope_start_genes.pdf'), widh = 8, height = 8)
# plotGene(testptObj = res, gene = gene[1:16], variable = 'condition')
# dev.off()
# pdf(paste0(plotdir, 'slope_start_fdr_foldchange.pdf'), widh = 4, height = 4)
# plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange')
# dev.off()
# ## 3
# res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = 'slope', test.position = 'middle')
# res3 <- res
# gene = names(sort(res$fdr))
# pdf(paste0(plotdir, 'slope_middle_gene.pdf'), widh = 8, height = 8)
# plotGene(testptObj = res, gene = gene[1:16], variable = 'condition')
# dev.off()
# pdf(paste0(plotdir, 'slope_middle_fdr_foldchange.pdf'), widh = 4, height = 4)
# plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange')
# dev.off()
# 
# ## 4
# res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = 'slope', test.position = 'end')
# res4 <- res
# gene = names(sort(res$fdr))
# pdf(paste0(plotdir, 'slope_end_genes.pdf'), widh = 8, height = 8)
# plotGene(testptObj = res, gene = gene[1:16], variable = 'condition')
# dev.off()
# pdf(paste0(plotdir, 'slope_end_fdr_foldchange.pdf'), widh = 4, height = 4)
# plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange')
# dev.off()
# 
# ##  21
# res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = 'all', test.position = 'all')
# res5 <- res
# gene = names(sort(res$fdr))
# pdf(paste0(plotdir, 'all_all_genes.pdf'), widh = 8, height = 8)
# plotGene(testptObj = res, gene = gene[1:16], variable = 'condition')
# dev.off()
# pdf(paste0(plotdir, 'all_all_fdr_foldchange.pdf'), widh = 4, height = 4)
# plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange')
# dev.off()
# 
# ##  22
# res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = 'all', test.position = 'start')
# res6 <- res
# gene = names(sort(res$fdr))
# pdf(paste0(plotdir, 'all_start_genes.pdf'), widh = 8, height = 8)
# plotGene(testptObj = res, gene = gene[1:16], variable = 'condition')
# dev.off()
# pdf(paste0(plotdir, 'all_start_fdr_foldchange.pdf'), widh = 4, height = 4)
# plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange')
# dev.off()
# 
# ## 23
# res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = 'all', test.position = 'middle')
# res7 <- res
# gene = names(sort(res$fdr))
# pdf(paste0(plotdir, 'all_middle_genes.pdf'), widh = 8, height = 8)
# plotGene(testptObj = res, gene = gene[1:16], variable = 'condition')
# dev.off()
# pdf(paste0(plotdir, 'all_middle_fdr_foldchange.pdf'), widh = 4, height = 4)
# plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange')
# dev.off()
# 
# ## 24
# res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1,type='Variable', test.type = 'all', test.position = 'end')
# res8 <- res
# gene = names(sort(res$fdr))
# pdf(paste0(plotdir, 'all_end_genes.pdf'), widh = 8, height = 8)
# plotGene(testptObj = res, gene = gene[1:16], variable = 'condition')
# dev.off()
# pdf(paste0(plotdir, 'all_end_fdr_foldchange.pdf'), widh = 4, height = 4)
# plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange')
# dev.off()
# 
# # plotHeatmap(testptObj = res, Mat = m, Pseudotime = pt, Cellanno = ca, Design = design, Max.gene = 500, heatmap.x = 0.4, heatmap.y = 0.5, heatmap.width = 0.8, heatmap.height = 1.0, dend.x = 0.90, dend.y = 0.44, dend.width = 0.2, dend.height = 0.98)
# # pd <- get_heatmap_plots(testptObj = res, Mat = m, Pseudotime = pt, Cellanno = ca, Design = design[,2,drop=F], Max.gene = NULL)
# # plotQuickHeatmap(pd)

