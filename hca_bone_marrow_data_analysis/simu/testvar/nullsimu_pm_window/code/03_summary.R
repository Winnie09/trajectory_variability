rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/'

n = NULL

# phe.fit = readRDS(paste0(rdir, 'phenopath/fit.rds'))
# phe.res = readRDS(paste0(rdir, 'phenopath/sig.rds'))
# n['phenopath'] = sum(phe.res == 'TRUE')

mon.res = readRDS(paste0(rdir, 'monocle2_trajtest/res.rds'))
n['monocle2_trajtest'] = sum(mon.res$fdr < 0.05)

mon.res = readRDS(paste0(rdir, 'condiments/cond_genes.rds'))
mon.res$pvalue[is.na(mon.res$pvalue)] <- 1
n['condiments'] = sum(p.adjust(mon.res$pvalue,method='fdr') < 0.05)

lam.res = readRDS(paste0(rdir, 'lamian/testvar_res.rds'))
n['Lamian'] = sum(lam.res[[1]][,1] < 0.05)

saveRDS(n, paste0(rdir, 'summary/num_false_positives.rds'))
write.csv(n, paste0(rdir, 'summary/num_false_positives.csv'))
names(phe.res) = names(mon.res)
int = intersect(names(phe.res)[phe.res], names(mon.res)[mon.res<0.05])
str(int)

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

sg = names(sort(mon.res[int])[1:32])
fit = getPopulationFit(lam.res, sg, type = 'variable')
plotGeneSampleAndPopulation(lam.res, gene = sg[2], variable = 'group', line.size= 0.5, line.alpha = 0.5, sep = ':.*', continuous = F, plot.point = T, point.alpha = 0.2, point.size = 0.1)

g = sg[2]
plot(fit[[1]][g,]~lam.res$pseudotime)
expr = lam.res$expr
plot(expr[g, ]~lam.res$pseudotime, col = as.factor(sub(':.*', '', names(lam.res$pseudotime))))


g.sd = sapply(int, function(g){
  mean(tapply(expr[g,], sub(':.*', '', names(lam.res$pseudotime)), sd))
})

pdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/plot/'  
png(paste0(pdir, 'false_positive1.png'), res = 300, height = 1000, width = 1200)
plotGeneSampleAndPopulation(lam.res, gene = names(sort(g.sd, decreasing = T)[1:4]), variable = 'group', line.size= 0.5, line.alpha = 0.5, sep = ':.*', continuous = F, plot.point = T, point.alpha = 0.1, point.size = 0.01)
dev.off()


png(paste0(pdir, 'false_positive1.png'), res = 300, height = 4500, width = 5400)
plotGeneSampleAndPopulation(lam.res, gene = names(sort(g.sd, decreasing = T)[1:36]), variable = 'group', line.size= 0.3, line.alpha = 0.5, sep = ':.*', continuous = F, plot.point = T, point.alpha = 0.2, point.size = 0.5)
dev.off()


int = names(sort(g.sd, decreasing = T))
df = data.frame(phenopath.fdr = phe.res[int], monocle2.fdr = mon.res[int], lam.res[[1]][int, 1], stringsAsFactors = F)

write.csv(df, paste0(rdir, 'summary/false_positive_fdr.csv'))


