library(here)
setwd(here())
source('function/01_function.R')

lamian <- readRDS('hca/real/testvar/result/EM_pm/monocyte/gender/gender_res.rds')
lamian$populationFit <- getPopulationFit(lamian, gene = rownames(lamian$expr), type = 'variable', num.timepoint = length(lamian$pseudotime))
colnames(lamian$populationFit[[1]]) <- colnames(lamian$populationFit[[2]]) <- names(sort(lamian$pseudotime))

celltype = 'monocyte'
genes <- readRDS(paste0('hca/real/testvar/plot/perf/',celltype,'_venn_significant_genes.rds'))
length(genes)
genes[['tradeSeqDT']] <- genes[['tradeSeq']][[1]]
genes[['tradeSeqPT']] <- genes[['tradeSeq']][[2]]
genes[['tradeSeqET']] <- genes[['tradeSeq']][[3]]
genes[['tradeSeq']] = unique(unlist(genes[['tradeSeq']]))

siggene <- list()
i = 'Lamian.pm'
other = setdiff(names(genes), i)
for (j in other){
  siggene[[paste0(i, '_exc', j)]] <- setdiff(genes[[i]], genes[[j]])
  siggene[[paste0(j, '_exc', i)]] <- setdiff(genes[[j]], genes[[i]])
}

pdir <- 'hca/real/testvar/plot/examplegene/monocyte/'
dir.create(pdir, recursive = T, showWarnings = F)

source('/home/whou10/scratch16/whou10/trajectory_variability/package/Lamian/R/plotGene.R')
for (group in names(siggene)) {
  print(group)
  set.seed(12345)
  targ = sample(siggene[[group]], min(100, length(siggene[[group]])))
  dir.create(paste0(pdir, group), recursive = T, showWarnings = F)
  for (i in targ){
    pdf(paste0(pdir, group, '/',i,'.pdf'), width = 2.6, height = 2)
    # plotGeneSampleAndPopulation(lamian,i,variable='gender')
    plotGene(lamian, i, variable = 'gender', continuous = F, axis.text.blank = T, sep = ':.*')
    dev.off()
  }
}


