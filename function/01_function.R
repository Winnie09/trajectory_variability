allf <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function')
allf = allf[!grepl('01_function.R',allf)]
res <- sapply(allf, function(f){
  source(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/', f))
})
rm(allf)
rm(res)

