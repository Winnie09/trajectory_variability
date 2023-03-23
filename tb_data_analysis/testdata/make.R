
ca <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/cellanno.rds')
de <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/case_design.rds')
co <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/count.rds')
ex <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/expr.rds')
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/pt.rds')
samp <- unique(ca[,2])
ss <- sample(samp,5)
ca <- ca[ca[,2]%in%ss,]
de <- de[ss,]
co <- co[,ca[,1]]
ex <- ex[,ca[,1]]
pt <- pt[ca[,1]]
ex <- ex[1:100,]
co <- co[1:100,]
saveRDS(ca,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/cellanno.rds')
saveRDS(de,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/case_design.rds')
saveRDS(co,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/count.rds')
saveRDS(ex,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/expr.rds')
saveRDS(pt,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/pt.rds')

