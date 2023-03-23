prop1 = as.numeric(commandArgs(trailingOnly = T)[[1]])
prop2 = as.numeric(commandArgs(trailingOnly = T)[[2]])
prop3 = as.numeric(commandArgs(trailingOnly = T)[[3]])
sampleID = as.numeric(commandArgs(trailingOnly = T)[[4]])
if ((prop1 + prop2 + prop3)!=1) {
  sum = prop1 + prop2 + prop3
  prop1 = prop1/sum
  prop2 = prop2/sum
  prop3 = prop3/sum
}
setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability')
mat = readRDS('./hca/data/HCA/proc/matrix/normcount.rds')
ct = readRDS('./hca/data/HCA/proc/ct/sc.rds')
ct <- ct[!is.na(ct)]
ct <- ct[grepl('HSC',ct) | grepl('MEP',ct) | grepl('Ery',ct)]
mat = mat[,names(ct)]
ap = gsub(':.*','',colnames(mat))
p = ap[sampleID]
tab = table(ct[colnames(mat)[ap == p]])
cttmp = ct[colnames(mat)[ap == p]]
numCell = floor(min(tab['HSC']/prop1,tab['MEP']/prop2,tab['Ery']/prop3) * 0.8)
set.seed(12345)
ct.alloc = sample(c('a','b','c'), numCell, prob = c(prop1,prop2,prop3), replace=T)
ct1 = names(sample(cttmp[cttmp=='HSC'], length(ct.alloc[ct.alloc=='a'])))
ct2 = names(sample(cttmp[cttmp=='MEP'], length(ct.alloc[ct.alloc=='b'])))
ct3 = names(sample(cttmp[cttmp=='Ery'], length(ct.alloc[ct.alloc=='c'])))
smat = mat[, c(ct1,ct2,ct3)]
dir.create(paste0('./simu/data/hscMepEry_',round(prop1,2),'_',round(prop2,2),'_',round(prop3,2),'/',p,'/matrix/'),recursive = T,showWarnings = F)
saveRDS(smat, paste0('./simu/data/hscMepEry_',round(prop1,2),'_',round(prop2,2),'_',round(prop3,2),'/',p,'/matrix/normcount.rds'))
ct = ct[colnames(smat)]
dir.create(paste0('./simu/data/hscMepEry_',round(prop1,2),'_',round(prop2,2),'_',round(prop3,2),'/',p,'/ct/'), recursive = T, showWarnings = F)
saveRDS(smat, paste0('./simu/data/hscMepEry_',round(prop1,2),'_',round(prop2,2),'_',round(prop3,2),'/',p,'/ct/ct.rds'))
rm(list=ls())

