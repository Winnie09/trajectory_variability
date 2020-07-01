allfiles <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function')
allfiles = allfiles[!grepl('01_function.R',allfiles)]
allfiles <- allfiles[grepl('.R$', allfiles)]
for (f in allfiles){
  source(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/', f))
}
rm(allfiles)
