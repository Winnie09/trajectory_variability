library(data.table)
af <- list.files('/home-4/zji4@jhu.edu/scratch/covid/data/EMTAB9357_pbmc/data/raw/pro/')
m <- sapply(af,function(f) {
  d <- fread(paste0('/home-4/zji4@jhu.edu/scratch/covid/data/EMTAB9357_pbmc/data/raw/pro/',f),data.table = F)  
  rownames(d) <- d[,1]
  d <- as.matrix(d[,-1])
  rownames(d) <- paste0(gsub('heathlab_dc_9_17_pbmc_pro_|.txt.gz','',f),':',sub(':.*','',rownames(d)))
  d
})
m <- do.call(rbind,m)
m <- t(m)
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/covid/data/EMTAB9357_pbmc/data/proc/prooriname.rds')
