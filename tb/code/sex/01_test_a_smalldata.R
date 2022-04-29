library(rhdf5)
pid <- as.numeric(commandArgs(trailingOnly = T))
ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex_addBatch/'

pt <- readRDS(paste0(ddir,'ptpc',pid,'.rds'))
cellanno <- readRDS(paste0(ddir, 'cellanno.rds'))
design <- readRDS(paste0(ddir, 'design.rds'))
expr = readRDS(paste0(ddir, 'expr.rds'))


## gene by cell matrix
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function',pattern = 'multi')
for (f in af) source(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/',f))

t2c <- testpt(expr=expr[c(1:4,which(rownames(expr)=='XIST')),], cellanno=cellanno, pseudotime=pt, design=design,ncores=30, testvar=ncol(design),test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'chisq')
saveRDS(t2c, 't2c.rds')

