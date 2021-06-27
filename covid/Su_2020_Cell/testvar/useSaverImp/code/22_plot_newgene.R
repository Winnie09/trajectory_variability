gl <- readLines('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/EM_pm/Mod_Mi/newgene/genelist.txt')
res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/Mod_Mi/numeric_res_with_clu.rds')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/EM_pm/Mod_Mi/'

png(paste0(pdir, 'newgene.png'), res = 200, width = 1600, height = 2500)
plotGene(res, gl, variable = 'type', ncol = 6)
dev.off()

selgene <- c('ATG14','PBRM1', 'VCP','SUCLG2','RNF11','RBPJ','PAK1','NCK2')
png(paste0(pdir, 'newgene_zeyu.png'), res = 200, width = 1000, height = 450)
plotGene(res, selgene, variable = 'type', ncol = 4)
dev.off()

selgene <- c('NCK2', 'PAK1','RBPJ', 'PBRM1')
png(paste0(pdir, 'newgene_zeyu4.png'), res = 300, width = 1000, height = 900)
plotGene(res, selgene, variable = 'type', ncol = 2, continuous = F, line.size = 0.1, line.alpha = 0.8)
dev.off()

pdf(paste0(pdir, 'newgene_zeyu4.pdf'), width = 4.1, height = 4)
plotGene(res, selgene, variable = 'type', ncol = 2, continuous = F, line.size = 0.1, line.alpha = 0.8)
dev.off()


