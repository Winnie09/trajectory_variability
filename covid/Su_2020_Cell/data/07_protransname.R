d <- readRDS('/home-4/zji4@jhu.edu/scratch/covid/data/EMTAB9357_pbmc/data/proc/prooriname.rds')
cn <- colnames(d)
cell <- sub('.*:','',cn)
lib <- sub('library_','',sub(':.*','',cn))
cn <- paste0(cell,":",sub('Mix-donor1','Mix_donor1',sub('_','-',lib)))
n <- readRDS('/home-4/zji4@jhu.edu/scratch/covid/data/EMTAB9357_pbmc/data/proc/transname.rds')
tmp <- names(n)
tmp <- do.call(rbind,strsplit(tmp,':'))
tmp <- paste0(tmp[,1],':',tmp[,3])
names(n) <- tmp
colnames(d) <- n[cn]
saveRDS(d,'/home-4/zji4@jhu.edu/scratch/covid/data/EMTAB9357_pbmc/data/proc/pro.rds')
