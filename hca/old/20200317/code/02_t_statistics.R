datapath = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/ery'
af = sub('.rds','',list.files(datapath))
af = gsub('spline_coef_','',af)
sex = sub('.*_','',af)
get_coef_list <- function(samples, datapath){
  m = list()
  for (j in 1:length(samples)){
    fn = paste0(datapath,'/spline_coef_',samples[j],'.rds')
    if (file.exists(fn)){
        coef = readRDS(fn)
        for (i in 1:ncol(coef)){
          if (j == 1) m[[i]] = coef[,i]
          else m[[i]] = cbind(m[[i]], coef[,i])  
        }  
    }
  }
  m  
}
g1 = as.character(unique(sex)[1])
ap1 <- af[sex==g1]
m1 = get_coef_list(ap1,datapath)

g2 = as.character(unique(sex)[2])
ap2 <- af[sex==g2]
m2 = get_coef_list(ap2,datapath)

t <- sapply(1:length(m1), function(i){
  sapply(rownames(m1[[1]]), function(gene){
    v1 = m1[[i]][gene,]
    v2 = m2[[i]][gene,]
    (mean(v1) - mean(v2))/(sqrt(var(v1)/length(v1) + var(v2)/length(v2)))
  })
})

saveRDS(t, paste0(datapath,'/t_statistics.rds'))
