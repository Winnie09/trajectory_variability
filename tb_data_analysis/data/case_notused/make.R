library(data.table)
u <- fread('/hpc/group/jilab/zj/tb/data/raw/tbru_unsupervised_updated.csv',data.table=F)
expr <- readRDS(paste0('/hpc/group/jilab/zj/tb/combine/saverfilter1k5.rds'))
cn <- sub('.*:','',colnames(expr))
sn <- sub(':.*','',colnames(expr))
colnames(expr) <- cn
u <- u[u[,1] %in% cn,]
u <- u[order(u$innateness),]
ord <- u[,1]
pt <- 1:length(ord)
names(pt) <- ord
expr <- expr[rowMeans(expr>0.1)>0.1,]
expr <- expr[,names(pt)]
saveRDS(pt,file='/hpc/group/jilab/zj/tb/data/proc/pt.rds')
saveRDS(expr,file='/hpc/group/jilab/zj/tb/data/proc/expr.rds')

meta <- fread('/hpc/group/jilab/zj/tb/data/raw/GSE158769_meta_data.txt.gz',data.table=F)
meta <- meta[,c('TB_status','donor')]
meta <- unique(meta)
meta <- meta[meta[,2] %in% sn,]
meta[,1] <- as.numeric(meta[,1]=='CASE')

cellanno <- data.frame(cell=cn,sample=sn,stringsAsFactors = F,row.names = cn)
cellanno <- cellanno[match(names(pseudotime), cellanno[,1]), ]
usn <- unique(sn)
design <- data.frame(intercept=1,contrast=meta[,1],row.names = meta[,2])
identical(sort(rownames(design)), sort(unique(cellanno[,2])))
saveRDS(cellanno,file='/hpc/group/jilab/zj/tb/data/proc/cellanno.rds')
saveRDS(design,file='/hpc/group/jilab/zj/tb/data/proc/case_design.rds')
