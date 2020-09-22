allfiles <- list.files('/dcl02/hongkai/data/whou/trajectory_variability/function')
allfiles = allfiles[!grepl('01_function.R',allfiles)]
allfiles <- allfiles[grepl('.R$', allfiles)]
for (f in allfiles){
  source(paste0('/dcl02/hongkai/data/whou/trajectory_variability/function/', f))
}
rm(allfiles)


