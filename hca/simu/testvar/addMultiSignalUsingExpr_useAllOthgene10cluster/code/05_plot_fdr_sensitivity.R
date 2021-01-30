rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/perf/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'

af = paste0(1:4, '.rds')
allm = c('EM', 'EM_meandiff', 'tradeSeq', 'monocle2', 'monocle3', 'tscan')
pd <- lapply(allm, function(m){
  print(m)
  tmp <- lapply(af, function(f){
    print(f)
    r = readRDS(paste0(ddir, m, '/', f))
    if (m == 'EM') {
      res = data.frame(fdr = r$fdr, foldchange = r$foldchange, stringsAsFactors = F)
      res = res[order(res[,1], -res[,2]), , drop=F]
      data.frame(SensFdr(selgene, res), Method = m, SignalStrength = sub('.rds', '', f), stringsAsFactors = FALSE)
    } else if (m == 'EM_meandiff'){
      res = r[order(r[,5], -abs(r[,1])), ]
      data.frame(SensFdr(selgene, res), Method = m, SignalStrength = sub('.rds', '', f), stringsAsFactors = FALSE)
    } else if (m == 'tradeSeq'){
      res = r[['earlyDETest']]
      a = c(sub('.rds','',f), 'tradeSeq_earlyDETest', AreaUnderSensFdr(SensFdr(selgene, res[['res']])))
      a = data.frame(SensFdr(selgene, res[['res']]), Method = 'tradeSeq_earlyDETest', SignalStrength = sub('.rds', '', f), stringsAsFactors = FALSE)
      res = r[['patternTest']]
      b = data.frame(SensFdr(selgene, res[['res']]), Method = 'tradeSeq_patternTest', SignalStrength = sub('.rds', '', f), stringsAsFactors = FALSE)
      res = r[['diffEndTest']]
      c = data.frame(SensFdr(selgene, res[['res']]), Method = 'tradeSeq_diffEndTest', SignalStrength = sub('.rds', '', f), stringsAsFactors = FALSE)
      rbind(a, b, c)
    } else if (m == 'monocle2' | m == 'tscan'){
      res = r[order(r[,3]), ]
      data.frame(SensFdr(selgene, res), Method = m, SignalStrength = sub('.rds', '', f), stringsAsFactors = FALSE)
    } else if (m == 'monocle3'){
      res = r[order(r[,3], -abs(r[,1])), ]
      data.frame(SensFdr(selgene, res), Method = m, SignalStrength = sub('.rds', '', f), stringsAsFactors = FALSE)
    } 
  })
  tmp <- do.call(rbind, tmp)  
})
pd2 <- do.call(rbind, pd)  
saveRDS(pd2, paste0(rdir, 'real_reported_fdr.rds'))

pd2[pd2$Method == 'EM', 'Method'] <- 'EM_trenddiff'

pdf(paste0(pdir, 'real_fdr_vs_reported_fdr.pdf'), width = 9,height = 7)
ggplot(data = pd2, aes(x = Reported_FDR, y = Real_FDR, color = Method))+
  geom_point(size = 0.2) +
  facet_wrap(~SignalStrength) +
  scale_color_brewer(palette = 'Dark2') +
  theme_classic()
dev.off()

pdf(paste0(pdir, 'sensitivity_vs_real_fdr.pdf'), width = 9,height = 7)
ggplot(data = pd2, aes(x = Real_FDR, y = Sensitivity, color = Method))+
  geom_point(size = 0.2) +
  facet_wrap(~SignalStrength) +
  scale_color_brewer(palette = 'Dark2') +
  theme_classic()
dev.off()


