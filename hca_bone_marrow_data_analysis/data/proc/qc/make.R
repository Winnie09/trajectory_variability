ddir <- '/home/whou10/scratch16/whou10/trajectory_variability/hca/data/proc/qc/'
n <- readRDS(paste0(ddir, 'expressedgenenumber.rds'))
m <- readRDS(paste0(ddir, 'mitoproportion.rds'))
c <- readRDS(paste0(ddir, 'totalreadcount.rds'))

pdir <- '/home/whou10/scratch16/whou10/trajectory_variability/hca/qcplot/'

pdf(paste0(pdir, 'genenumber.pdf'))
plot(density(n),main='Number of expressed genes')
abline(v=1000,col='red')
dev.off()

pdf(paste0(pdir, 'mitoprop.pdf'))
plot(density(m),main='Mito gene proportion')
abline(v=0.1,col='red')
dev.off()

pdf(paste0('readnumber.pdf'))
plot(density(c),main='Number of reads')
abline(v=5000,col='red')
dev.off()

library(scattermore)
pdf(paste0(pdir, 'genenumber_readnumber.pdf'), width = 4, height = 4.5)
scattermoreplot(n, c, col=heat.colors(length(n), alpha=.1), main='', xlab = 'Number of expressed genes', ylab = 'Number of reads')
abline(h = 5000, col='red')
abline(v = 1000, col='blue')
dev.off()

pdf(paste0(pdir, 'mito_readnumber.pdf'), width = 4, height = 4.5)
scattermoreplot(m, c, col=heat.colors(length(n), alpha=.1), main='', xlab = 'Mito gene proportion', ylab = 'Number of reads')
abline(h = 5000, col='red')
abline(v = 0.1, col='blue')
dev.off()

pdf(paste0(pdir, 'genenumber_mito.pdf'), width = 4, height = 4.5)
scattermoreplot(n, m, col=heat.colors(length(n), alpha=.1), main='', xlab = 'Number of expressed genes', ylab = 'Mito gene proportion')
abline(h = 0.1, col='red')
abline(v = 1000, col='blue')
dev.off()


