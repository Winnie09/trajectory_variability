library(here)
setwd(here())
ddir <- 'hca/real/testvar/result/EM_pm/'
path = 'erythroid'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source('function/01_function.R')

for (path in c('monocyte', 'lymph', 'erythroid')){
  print(path)
  Res <- readRDS(paste0(ddir, path, '/age/age_res.rds'))
  go <- GOEnrich(Res, use.clusters = F)
  go <- go[[1]]
  go <- go[go[,'FDR']<0.05, ,drop=F]
  if (nrow(go) > 0){
    write.csv(go, paste0(ddir, path, 'age/age_res_GO.csv'))
  }
}

  

