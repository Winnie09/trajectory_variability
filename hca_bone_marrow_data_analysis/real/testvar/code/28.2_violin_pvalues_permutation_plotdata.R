rm(list=ls())
library(here)
library(ggplot2)
library(RColorBrewer)
library(parallel)
# setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
setwd('/home/whou10/scratch16/whou10/')
# source('function/01_function.R')

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('resource/chrX_genename.rds')
u2 = readRDS('resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')
n.permute = 1e4

## "erythroid" "lymph"     "monocyte"
for (celltype in setdiff(list.files('trajectory_variability/hca/real/testvar/result/EM_pm'),'perf')) {
  genes <- readRDS(paste0('trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_violin_plotdata_genes.rds'))
  tradeseq <- readRDS(paste0('trajectory_variability/hca/real/testvar/result/tradeSeq/',celltype,'/gender/testvar_res.rds'))
  allg <- sub(':.*','',rownames(tradeseq[[1]]))
  permud <- reald <- NULL
  for (met in names(genes)) {
    print(met)
    gn <- genes[[met]]
    allg <- sub(':.*','',gn)
    v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
    v1 <- mean(v1[seq(1,length(v1)) %% 10 == 0])
    v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
    v2 <- mean(v2[seq(1,length(v2)) %% 10 == 0])
    # saveRDS(v1, paste0(rdir, path, '/gender/gender_chrX_overlap.rds'))
    # saveRDS(v2, paste0(rdir, path, '/gender/gender_chrY_overlap.rds'))
    ##### permute reported gene order
    v1_pm <- unlist(mclapply(seq(1,n.permute), function(myseed){
      set.seed(myseed+100)
      allg.pm = sample(allg)
      tmp = cumsum(allg.pm %in% u1)/seq(1,length(allg.pm))
      mean(tmp[seq(1,length(tmp)) %% 10 == 0])
    },mc.cores=detectCores()))
    
    v2_pm <- unlist(mclapply(seq(1,n.permute), function(myseed){
      set.seed(myseed+100)
      allg.pm = sample(allg)
      tmp = cumsum(allg.pm %in% u2)/seq(1,length(allg.pm))
      mean(tmp[seq(1,length(tmp)) %% 10 == 0])
    },mc.cores=detectCores()))
    
    permud <- rbind(permud,data.frame(per = c(v1_pm,v2_pm), type=rep(c('chrX','chrY'),each=n.permute),method=met, stringsAsFactors = FALSE))
    reald <- rbind(reald,data.frame(per=c(v1,v2),type=c('chrX','chrY'),pvalue=c(mean(v1_pm >= v1),mean(v2_pm >= v2)),method=met,stringsAsFactors = F))
  }
  saveRDS(list(permud = permud, reald = reald), paste0('trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_pvalue_violin_plotdata.rds'))
}

