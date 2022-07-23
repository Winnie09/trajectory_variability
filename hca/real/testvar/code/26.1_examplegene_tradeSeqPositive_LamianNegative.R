lamian <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/monocyte/gender/gender_res.rds')
lamsig <- rownames(lamian[[1]])[lamian[[1]][,1] < 0.05]
tradeseq <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/tradeSeq/monocyte/gender/testvar_res.rds')

tradesig <- unique(unlist(sapply(tradeseq,function(i) rownames(i)[which(i$adj.P.Val < 0.05)])))
targ <- setdiff(tradesig,lamsig)

library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/monocyte/'

source('function/01_function.R')

lamian$populationFit <- getPopulationFit(lamian, rownames(lamian$expr), 'variable')
"BCLAF1:ENSG00000029363"
"CHPT1:ENSG00000111666"
set.seed(12345)
for (i in sample(targ,100)) {
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/examplegene_tradeSeqPositive_LamianNegative/',i,'.pdf'), width = 2.6, height = 2)
  # plotGeneSampleAndPopulation(lamian,i,variable='gender')
  plotGene(lamian, i, variable = 'gender', continuous = F, axis.text.blank = T, sep = ':.*')
  dev.off()
}

