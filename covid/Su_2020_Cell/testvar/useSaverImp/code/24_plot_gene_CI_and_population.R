lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/Mod_Mi/numeric_res.rds'))
lamiansig <- rownames(lamian[[1]])[lamian[[1]][,'fdr.overall'] < 0.05]
limma <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/limma/Mod_Mi.rds')
limma <- rownames(limma)[limma$adj.P.Val < 0.05]
tradeseq <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/tradeSeq/Mod_Mi/testvar_res.rds')
tradeseq <- unique(unlist(sapply(tradeseq,function(i) {rownames(i)[i$adj.P.Val < 0.05]})))
tradeseq <- tradeseq[!is.na(tradeseq)]

newgene <- setdiff(lamiansig,c(limma,tradeseq))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

for (g in newgene) {
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/addSE/lamianonly/',g,'.pdf'))
  print(plotGeneCIAndPopulation(lamian,g,variable='type',line.alpha=0.2))
  dev.off()
}

for (g in sample(setdiff(limma,lamiansig),50)) {
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/addSE/limmaonly/',g,'.pdf'))
  print(plotGeneCIAndPopulation(lamian,g,variable='type',line.alpha=0.2))
  dev.off()
}

