setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA')
mat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/proc/matrix/saver.rds')
order =readRDS('./order.rds')
mat = mat[,order$Cell]
order = cbind(order, Sample = sub('_.*','',order$Cell), Gender = gsub('.*:','',sub('_.*','',order$Cell)))
ap = sub('_.*','',colnames(mat))
ap1 = as.character(unique(order[which(order[,'Gender']=='female'),'Sample']))
m1 = mat[, order[order[,'Sample']%in%ap1,'Cell']]
ap2 = as.character(unique(order[which(order[,'Gender']=='male'),'Sample']))
m2 = mat[, order[order[,'Sample']%in%ap2,'Cell']]

df = readRDS('./ery/significant_genes.rds')

get_plotdata <- function(samples, g){
  pd = NULL
  pd <- lapply(samples, function(i){
    if (file.exists(paste0('./spline_coef_',gsub(':','_',i),'.rds'))){
      print(i)
      tmp = mat[g, grepl(i,ap)]
      time = order[match(names(tmp),order$Cell),'Pseudotime']
      coef = readRDS(paste0('./spline_coef_',gsub(':','_',i),'.rds'))
      
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
      pred = base %*% coef[g,]
      data.frame(time = time, expr = tmp, pred = pred, sample = i, stringsAsFactors = F)  
    }
  })
  # pd = pd[!sapply(pd, function(i) is.null(nrow(i)))]
  pd = do.call(rbind,pd)
}

allg = rownames(df)[df[,2]>0]
for (i in 1:length(allg)){
  g = rownames(df)[df[,2]>0][i]
  print(g)
  pd1 = get_plotdata(ap1,g)
  pd2 = get_plotdata(ap2,g)
  library(ggplot2)
  pdf(paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/plot/HCA/ery/',g,'_1.pdf'),width=12,height=3)
  print(ggplot(data=pd1) + geom_point(aes(x=time, y = expr),size=0.5,alpha=0.5) + geom_line(aes(x=time, y=pred),size=2,color='red',alpha=0.5) + facet_grid(~sample)) + ggtitle(g)
  dev.off()
  
  pdf(paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/plot/HCA/ery/',g,'_2.pdf'),width=12,height=3)
  print(ggplot(data=pd2) + geom_point(aes(x=time, y = expr),size=0.5,alpha=0.5) + geom_line(aes(x=time, y=pred),size=2,color='red',alpha=0.5) + facet_grid(~sample))+ ggtitle(g)
  dev.off()
}




for (i in 1:50){
  g = rownames(df)[df[,2]<0][i]
  print(g)
  pd1 = get_plotdata(ap1,g)
  pd2 = get_plotdata(ap2,g)
  library(ggplot2)
  pdf(paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/plot/HCA/ery_insig/',g,'_1.pdf'),width=12,height=3)
  print(ggplot(data=pd1) + geom_point(aes(x=time, y = expr),size=0.5,alpha=0.5) + geom_line(aes(x=time, y=pred),size=2,color='red',alpha=0.5) + facet_grid(~sample)) + ggtitle(g)
  dev.off()
  
  pdf(paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/plot/HCA/ery_insig/',g,'_2.pdf'),width=12,height=3)
  print(ggplot(data=pd2) + geom_point(aes(x=time, y = expr),size=0.5,alpha=0.5) + geom_line(aes(x=time, y=pred),size=2,color='red',alpha=0.5) + facet_grid(~sample))+ ggtitle(g)
  dev.off()
}
