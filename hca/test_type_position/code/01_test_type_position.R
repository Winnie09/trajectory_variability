library(parallel)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
plotdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/test_type_position/plot/testvar/clusterType9_1/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/test_type_position/result/testvar/clusterType9_1/'
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/data/data/count/clusterType9_1.rds')
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
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
##  test and plot
tmp <- mclapply(c('all', 'start', 'middle', 'end'), function(pos){
      ## slope only
      res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1, type='Variable', test.slope.only = TRUE, test.position = pos)
      saveRDS(res, paste0(rdir, 'slope_', pos, '.rds'))
      gene = names(rev(sort(abs(res$meandiff))))
      pdf(paste0(plotdir, 'slope_', pos, '_meandiff_genes.pdf'), width = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      
      gene = names(sort(res$fdr))
      pdf(paste0(plotdir, 'slope_', pos, '_fdr_genes.pdf'), width = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      pdf(paste0(plotdir, 'slope_', pos, '_fdr_foldchange.pdf'), width = 4, height = 4)
      print(plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange'))
      dev.off()
      pdf(paste0(plotdir, 'slope_', pos, '_fdr_meandiff.pdf'), width = 4, height = 4)
      print(plot(res$meandiff ~ res$fdr[names(res$meandiff)], pch = 20, xlab = 'fdr', ylab = 'group mean difference'))
      dev.off()
      
        
      ## all (intersept + slope)
      res <- testpt(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1, type='Variable', test.slope.only = FALSE, test.position = pos)
      saveRDS(res, paste0(rdir, 'all_', pos, '.rds'))
      gene = names(rev(sort(abs(res$meandiff))))
      pdf(paste0(plotdir, 'all_', pos, '_meandiff_genes.pdf'), width = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      
      gene = names(sort(res$fdr))
      pdf(paste0(plotdir, 'all_', pos, '_fdr_genes.pdf'), width = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      pdf(paste0(plotdir, 'all_', pos, '_fdr_foldchange.pdf'), width = 4, height = 4)
      print(plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange'))
      dev.off()
      pdf(paste0(plotdir, 'slope_', pos, '_fdr_meandiff.pdf'), width = 4, height = 4)
      print(plot(res$meandiff ~ res$fdr[names(res$meandiff)], pch = 20, xlab = 'fdr', ylab = 'group mean difference'))
      dev.off()
      
      return(0)
}, mc.cores = 8)
