input.seed = as.numeric(commandArgs(trailingOnly = T)[1])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/')
mat = readRDS('./proc/matrix/saver.rds')
vargene = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/variable_genes_cv0.5.rds')
mat = mat[vargene,]
order = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/order.rds')
mat = mat[, order$Cell]
datapath = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/ery'
t_true = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/t_statistics.rds')
t_true = t_true[vargene,]
pval = 0
dir.create(paste0(datapath,'/',input.seed))
mat.bak = mat
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
for (myseed in (input.seed*100) : (input.seed*100 + 99)){    ################ CHANGE TO 99
  print(myseed)
  ## <----- randomization 
  set.seed(myseed)
  id = sample(1:nrow(order),nrow(order), replace=T) #######
  mat = mat.bak[, order$Cell[id]]
  ## --------------------->
  res <- sapply(unique(sub('_.*','',colnames(mat))), function(s){
    trainData = mat[, grepl(s, colnames(mat))]
    trainX = order[match(colnames(trainData),order$Cell), 'Pseudotime']
    para = get_spline_coefficient(trainData, trainX, fit.min = min(order$Pseudotime), fit.max = max(order$Pseudotime))
    saveRDS(para, paste0(datapath,'/', input.seed,'/spline_coef_', gsub(':','_',s), '.rds'))
    return(0)
  })
        
  ########## a new file
  af = list.files(paste0(datapath,'/',input.seed))
  af = af[grepl('spline_coef_',af)]
  af = sub('.rds','',gsub('spline_coef_','',af))
  ## <------   permute sample gender lables 
  set.seed(myseed)
  sex = sample(sub('.*_','',af))
  ## ------------------------------------>
  get_coef_list <- function(samples, datapath){
      m = list()
      for (j in 1:length(samples)){
        fn = paste0(datapath,'/',input.seed,'/spline_coef_',samples[j],'.rds')
        if (file.exists(fn)){
            coef = readRDS(fn)
            for (i in 1:ncol(coef)){
              if (j == 1) m[[i]] = coef[,i,drop=F]
              else m[[i]] = cbind(m[[i]], coef[,i,drop=F])  
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
  pval <- ((abs(t)>t_true)-0) + pval
  str(pval)
}
saveRDS(pval, paste0(datapath,'/',input.seed,'/pvalue.rds'))
