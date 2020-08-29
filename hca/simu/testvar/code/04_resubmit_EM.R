setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/code')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/testvar/result/'
method = 'EM_SelectKnots'

for (clusterType in 1:10){
  print(clusterType)
  for (pctGene in 1:4){
    fn <- paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds')
    if (file.exists(fn)) {
      break
    } else {
      system(paste0('qsub run02.EM.sh ', clusterType, ' ', pctGene))
    }
  }
}
  
