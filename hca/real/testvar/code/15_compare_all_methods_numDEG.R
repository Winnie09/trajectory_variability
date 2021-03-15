setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_chisq')
numDEG <- sapply(c('lymph', 'erythroid', 'monocyte'), function(path){
  print(path)
  sapply(c('gender_res.rds', 'age_res.rds'), function(fn){
    print(fn)
    res = readRDS(paste0(path, '/', fn))
    s = res$statistics
    sum(s[, grepl('^fdr.*overall$', colnames(s))] < 0.05)
  })
  
})
    


setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm')
numDEG <- sapply(c('lymph', 'erythroid', 'monocyte'), function(path){
  print(path)
  sapply(c('gender_res.rds', 'age_res.rds'), function(fn){
    print(fn)
    if (file.exists(paste0(path, '/', fn))){
      res = readRDS(paste0(path, '/', fn))
      s = res$statistics
      sum(s[, grepl('^fdr.*overall$', colnames(s))] < 0.05)
    } else {
      NA
    }
      
  })
  
})
    


setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar_not_NIW/result/')
numDEG <- sapply(c('lymph', 'erythroid', 'monocyte'), function(path){
  print(path)
  sapply(c('meandiff_gender_res.rds', 'meandiff_age_res.rds'), function(fn){
    print(fn)
    if (file.exists(paste0(path, '/', fn))){
      res = readRDS(paste0(path, '/', fn))
      print(head(res))
      print(summary(res[, grepl('adj.P.Val', colnames(res))]))
      sum(res[, grepl('adj.P.Val', colnames(res))] < 0.05)
    } else {
      NA
    }  
  })
})


