# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_chisq')
# numDDG <- sapply(c('lymph', 'erythroid', 'monocyte'), function(path){
#   print(path)
#   sapply(c('gender_res.rds', 'age_res.rds'), function(fn){
#     print(fn)
#     res = readRDS(paste0(path, '/', fn))
#     s = res$statistics
#     sum(s[, grepl('^fdr.*overall$', colnames(s))] < 0.05)
#   })
#   
# })
    


setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm')
numDDG <- sapply(c('lymph', 'erythroid', 'monocyte'), function(path){
  print(path)
  sapply(c('gender_res.rds', 'age_res.rds'), function(fn){
    print(fn)
    if (file.exists(paste0(path, '/', sub('_.*','',fn), '/', fn))){
      res = readRDS(paste0(path, '/', sub('_.*','',fn),'/',fn))
      s = res$statistics
      sum(s[, grepl('^fdr.*overall$', colnames(s))] < 0.05)
    } else {
      NA
    }
  })
})
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/perf/', recursive = T, showWarnings = F)
write.csv(numDDG, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/perf/numDDG.csv')
  

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar_not_NIW/result/')
numDDG <- sapply(c('monocyte', 'lymph', 'erythroid'), function(path){
  print(path)
  sapply(list.files(path, pattern = '.rds'), function(fn){
    print(fn)
    if (file.exists(paste0(path, '/', fn))){
      res = readRDS(paste0(path, '/', fn))
      if ('adj.P.Val' %in% colnames(res)){
        sum(res[, grepl('adj.P.Val', colnames(res))] < 0.05)
      } else {
        sum(res[[1]][, grepl('both.fdr', colnames(res[[1]]))] < 0.05)
      }
        
    } else {
      NA
    }  
  })
})
write.csv(numDDG, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar_not_NIW/result/perf/numDDG.csv')





