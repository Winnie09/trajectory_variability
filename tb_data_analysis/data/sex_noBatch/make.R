library(data.table)
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/saver.rds')
u <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/raw/mRNA_NAM_PCs.csv',data.table=F)
s <- fread('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/tb/data/raw/GSE158769_meta_data.txt',data.table=F)
u[,3] <- paste0(s[match(u[,3],s[,1]),'donor'],':',u[,3])

int <- intersect(colnames(expr),u[,3])
expr <- expr[,int]
sn <- sub(':.*','',colnames(expr))
tab <- table(sn)
id <- which(sn%in%names(tab)[tab>=1000])
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.01,]
u <- u[match(colnames(expr),u[,3]),]

identical(u[,3],colnames(expr))

sex <- readRDS('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/tb/data/sex/sex_female.rds')

cellanno <- data.frame(cell=colnames(expr),sample=sub(':.*','',colnames(expr)),stringsAsFactors = F,row.names = colnames(expr))
print(mean(cellanno[,2] %in% names(sex)))
sex <- sex[unique(cellanno[,2])] ## female = 1, male = 0

design <- data.frame(intercept=1,contrast=as.numeric(sex),row.names = names(sex))

identical(sort(rownames(design)), sort(unique(cellanno[,2])))

saveRDS(cellanno,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
saveRDS(design,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')
saveRDS(expr,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds')


pt <- u[,4]
names(pt) <- u[,3]
pt <- sort(pt)
n <- names(pt)
pt <- 1:length(pt)
names(pt) <- n
saveRDS(pt,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc2.rds')

pt <- u[,5]
names(pt) <- u[,3]
pt <- sort(pt)
n <- names(pt)
pt <- 1:length(pt)
names(pt) <- n
saveRDS(pt,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc4.rds')
