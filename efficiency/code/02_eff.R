library(chron)
library(reshape2)
timefunc <- function(i) {
  if (!grepl('-',i)) {
    t = as.numeric(times(i))
  } else {
    t = as.numeric(sub('-.*','',i))+as.numeric(times(sub('.*-','',i)))
  }
  round(24 * t,3)
}

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/code/')
af = list.files(getwd(), pattern = '.txt$')

eff <- lapply(af,function(f){
  print(f)
  d=readLines(f)
  
  a = sapply(d[grep('whou10',d)+1],function(i) {
    tmp <- strsplit(i,' ')[[1]]
    tmp <- tmp[nchar(tmp) > 0]
    if ('COMPLETED' %in% tmp) {
      tmp[(length(tmp)-1):length(tmp)]  
    } else if ('FAILED' %in% tmp | 'RUNNING' %in% tmp | 'CANCELLED' %in% tmp) {
      c(NA,NA)
    }
  },USE.NAMES = F, simplify = F)
  if (is.list(a)){
    a = a[!is.na(a)]
    a = do.call(cbind,a)
  }
  a = a[,order(as.numeric(sub('K','',a[2,]))), drop = F]
  a[1,] <- sapply(a[1,],timefunc)
  rownames(a) = c('time(h)', 'memory(GB)')
  a.mt = melt(a)
  df = data.frame(evaluation=as.character(a.mt[,1]), id=as.character(a.mt[,2]), value=as.character(a.mt[,3]), method = sub('_.*', '', f), data = strsplit(f, '_')[[1]][2], test = strsplit(f, '_')[[1]][3], stringsAsFactors = F)
  df[,3] = sapply(seq(1,nrow(df)), function(ii) {
    if(df[ii, 1] == 'memory(GB)'){
      if (grepl('K',df[ii, 3])){
      as.numeric(sub('K','', df[ii, 3]))/1e6
      } else {
      as.numeric(df[ii, 3])/1e9  
      }
    } else {
        df[ii, 3]
      }  
  })
  df
})
pd = do.call(rbind, eff)
pd[,1] = as.character(pd[,1])
pd[,2] = as.character(pd[,2])
pd[,3] = as.numeric(pd[,3])


time = pd[pd[,1] == 'time(h)', ]
mem = pd[pd[,1] == 'memory(GB)', ]
time = time[time[,4] != 'phenopath3', ]
time = time[time[,4] != 'phenopath500', ]
mem = mem[mem[,4] != 'phenopath500', ]
mem = mem[mem[,4] != 'phenopath3', ]


saveRDS(time, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/time.rds')
saveRDS(mem, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/memory.rds')

write.csv(time, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/time.csv')
write.csv(mem, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/memory.csv')


