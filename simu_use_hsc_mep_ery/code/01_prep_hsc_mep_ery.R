mat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/HCA/proc/matrix/normcount.rds')
ct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/HCA/proc/ct/sc.rds')
ct = ct[colnames(mat)]
smat = mat[,ct == 'HSC' | ct == 'MEP' | ct == 'Ery']
saveRDS(smat,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/hscMepEry/matrix/normcount.rds')
ct  = ct[colnames(smat)]
saveRDS(ct,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/hscMepEry/ct/ct.rds')

