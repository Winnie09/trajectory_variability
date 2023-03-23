## if we adjust for tb status (contrl/case) in XDE test, then we can only use the samples that have both sexes in each unique combination of batch and tb status, then we will have 216921 cells from 99 samples (47 case, 52 control. 48 female, 51 male)
ddir = '/home/whou10/scratch16/whou10/trajectory_variability/tb/data/sex/'
tb <- read.table('/home/whou10/scratch16/whou10/trajectory_variability/tb/data/raw/GSE158769_meta_data.txt', sep = '\t',header=T)
cc <- unique(tb[,c('batch','TB_status','donor')])
design <- readRDS(paste0(ddir, 'design.rds')) ## 1 female
rownames(design) <- sub('_.*', '', rownames(design))
cc <- cc[match(rownames(design),cc[,3]),]
cc$sex <- design[,'sex']
comb <- paste0(cc[,1],cc[,2])
unisex <- tapply(cc[,4],list(comb),function(i) length(unique(i)))

selcomb <- names(which(unisex==2))
selid <- cc$donor[comb%in%selcomb]  ## 99 samples 
saveRDS(selid,file='/home/whou10/scratch16/whou10/trajectory_variability/tb/res/sex_Batch/selid.rds')
