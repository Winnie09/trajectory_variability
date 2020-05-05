seed = as.numeric(commandArgs(trailingOnly = T)[[1]])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA')
mat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/proc/matrix/saver.rds')
order =readRDS('./order.rds')
mat = mat[,order$Cell]
order = cbind(order, Sample = sub('_.*','',order$Cell), Gender = gsub('.*:','',sub('_.*','',order$Cell)))
ap = sub('_.*','',colnames(mat))
ap1 = as.character(unique(order[which(order[,'Gender']=='female'),'Sample']))
ap2 = as.character(unique(order[which(order[,'Gender']=='male'),'Sample']))

df = readRDS('./significant_genes.rds')

mat.bak = mat
get_plotdata <- function(samples, g, seed){
  myseed = seed * 100 + 99          ###### adapt to permutation times
  ## <--------randomization 
  set.seed(myseed)
  id = sample(1:nrow(order),nrow(order), replace=T) 
  mat = mat.bak[, order$Cell[id]]
  ## --------------------->
  pd = NULL
  pd <- lapply(samples, function(i){
    if (file.exists(paste0('./ery/',seed,'/spline_coef_',gsub(':','_',i),'.rds'))){
      print(i)
      tmp = mat[g,grepl(i,colnames(mat))]
      time = order[match(names(tmp),order$Cell),'Pseudotime']
      coef = readRDS(paste0('./ery/',seed,'/spline_coef_',gsub(':','_',i),'.rds')) ## comment
      trainX = time
      num.base = 10
      knots = seq(min(order$Pseudotime),max(order$Pseudotime),length.out=num.base+2)[2:(num.base+1)]
      library(splines)
      base = cbind(1,bs(trainX,knots = knots))
      colidx = NULL
      for (ii in seq(2,ncol(base))){
        if (length(unique(base[,ii])) == 1) colidx = c(colidx, ii)
      }
      if (length(colidx)) base = base[,-colidx]
      # para = chol2inv(chol(crossprod(base))) %*% t(base) %*% tmp ### fit here instead of using coef
      # pred = t(base %*% para)
      pred = base %*% coef[g,]
      data.frame(time = time, expr = tmp, pred = pred, sample = i, stringsAsFactors = F)  
    }
  })
  pd = do.call(rbind,pd)
}

for (i in 1:20){
  g = rownames(df)[df[,2]>0][i]
  print(g)
  pd1 = get_plotdata(ap1,g,seed)
  pd2 = get_plotdata(ap2,g,seed)
  library(ggplot2)
  dir.create(paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/plot/HCA/ery_pm/',seed,'/'), showWarnings = F, recursive = T)
  pdf(paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/plot/HCA/ery_pm/',seed,'/',g,'_1.pdf'),width=12,height=3)
  print(ggplot(data=pd1) + geom_point(aes(x=time, y = expr),size=0.5,alpha=0.5) + geom_line(aes(x=time, y=pred),size=2,color='red',alpha=0.5) + facet_grid(~sample)) + ggtitle(g)
  dev.off()
  
  pdf(paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/plot/HCA/ery_pm/',seed,'/',g,'_2.pdf'),width=12,height=3)
  print(ggplot(data=pd2) + geom_point(aes(x=time, y = expr),size=0.5,alpha=0.5) + geom_line(aes(x=time, y=pred),size=2,color='red',alpha=0.5) + facet_grid(~sample))+ ggtitle(g)
  dev.off()
}

