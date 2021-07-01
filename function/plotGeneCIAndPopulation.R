plotGeneCIAndPopulation <- function(testobj, gene, variable = NULL, variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = FALSE,  point.alpha=1, point.size=0.5, ribbon.alpha = 0.3, ribbon.size = 1, line.size = 1, line.alpha = 1, continuous = TRUE, sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = F, x.lab='Pseudotime', y.lab = 'Expression', legend.position = 'right', axis.size = 8, title.size = 8, axis.text.size = 8){
  ## testobj: the output of function testpt() which is a list containing fdr, etc..
  ## variable: character, the variable (covariate) to color the samples, should be null or one of the column names of design matrix. Default is NULL, meaning each sample is colored differently. Otherwise, samples are colored by the variable (covariate) values.
  ## variable.text: a character vector. The text for the legend of the plot, corresponding to each variable values.
  ## continuous: if TRUE, samples are colored using viridis continuous colors. If FALSE, RColorBrewer "Dark2" discrete palette.
  ## expression: a character ('demean',or 'original') to define the expression values shown on the plots. if "demean"(default), show demeaned expression. if 'original", show original gene expression.
  ## ncol: only functional when plotting multiple genes. Used to define the number of columns. 
  library(splines)
  library(ggplot2)
  library(gridExtra)
  library(viridis)   
  library(reshape2)
  library(RColorBrewer) ##
  pseudotime <- testobj[['pseudotime']]
  cellanno <- testobj[['cellanno']]
  colnames(cellanno) <- c('Cell', 'Sample')
  if ('expr.ori' %in% names(testobj)) expression <- testobj[['expr.ori']] else expression <- testobj[['expr']]
  if ('populationFit' %in% names(testobj)) fit <- testobj$populationFit else fit = getPopulationFit(testobj, gene = gene, type = testobj$test.type)
  predict.values <- predict_fitting(testobj, gene= gene, test.type = testobj$test.type)
  
  pseudotime = pseudotime[colnames(expression)]
  cellanno <- cellanno[match(colnames(expression), cellanno[,1]), ]
  # predict.values <- predict.values[, colnames(expression),drop=F]
  knotnum <- testobj$knotnum
  knotnum[knotnum==0] <- 1  ## in case the fitting of line would cause bugs
  design <- testobj[['design']]
  
  
  cellanno <- data.frame(Cell = as.character(cellanno[,1]), 
                         Sample = as.character(cellanno[,2]), stringsAsFactors = FALSE)
  variable.d <- if(is.null(variable)) 1 else variable
  if (!is.null(variable.text) & variable.d != 1) {
    design[,variable.d] <- ifelse(design[, variable.d] == 0, variable.text[1], variable.text[2])
  }
  a <- if (free.scale) 'free' else 'fixed'
  
  pd <- data.frame(expr = expression[gene, ], 
                   Sample = cellanno[,2], 
                   Variable = design[match(cellanno[,2], rownames(design)), variable.d],  ##
                   pseudotime = pseudotime[colnames(expression)])
  pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
  linedlist <- lapply(unique(cellanno[,2]), function(p){
    # tmpcell <- cellanno[cellanno[,2]==p,1]
    tmpcellid <- which(cellanno[,2]==p)
    if (toupper(testobj$test.type) == 'TIME'){  ##### add
      tmpdf <- data.frame(expr=predict.values[gene,tmpcellid], 
                          Sample=p, 
                          Variable = design[rownames(design) == p, variable.d], ##
                          pseudotime=pseudotime[tmpcellid])
    } else {
      tmpdf <- data.frame(
        expr = predict.values[[which(as.numeric(sub('.*_', '', names(predict.values))) == design[p, variable.d])]][gene, tmpcellid],
        Sample = p,
        Variable = design[rownames(design) == p, variable.d], ##
        pseudotime=pseudotime[tmpcellid])
    }
  })
  ld <- do.call(rbind, linedlist)
  ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
  ld <- ld[order(ld$pseudotime), ] ## add 20200812
  
  ### add population Line >>>>>>>>>>>>>>
  pd.pop <- sapply(1:length(fit), function(i){
    tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
    tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
    colnames(tmp) <- c('g', 'pseudotime', 'expr')
    tmp <- data.frame(tmp, 
                      type = names(fit)[i],
                      stringsAsFactors = FALSE)
  }, simplify = FALSE)
  pd.pop <- do.call(rbind, pd.pop)
  if (!is.na(sep)) {
    pd.pop$g <- sub(sep, '', pd.pop$g)
    pd.pop$g <- factor(as.character(pd.pop$g), levels = sub(sep, '', gene))
  } else{
    pd.pop$g <- factor(as.character(pd.pop$g), levels = gene)
  }
  pd.pop$type <- sub('.*_', '', pd.pop$type)
  
  
  ldim <- t(sapply(unique(ld$Sample),function(i) {
    tmp <- ld[ld$Sample==i,]
    approx(tmp$pseudotime,tmp$expr,xout=1:max(ld$pseudotime),rule=2)$y
  }))
  rownames(ldim) <- unique(ld$Sample)
  
  ldg <- do.call(rbind,lapply(unique(design[,variable]),function(s) {
    tmp <- pd.pop[pd.pop$type==s,]
    mv <-approx(tmp$pseudotime,tmp$expr,xout=1:max(ld$pseudotime),rule=2)$y
    tmp <- ldim[names(which(design[,variable]==s)),]
    sdv <- apply(tmp,2,sd)/sqrt(nrow(tmp))
    # data.frame(ymin=mv+qnorm(0.025)*sdv,ymax=mv+qnorm(0.975)*sdv,Variable=s,pseudotime=1:max(ld$pseudotime))
    data.frame(ymin=mv-1.96*sdv,ymax=mv+1.96*sdv,Variable=s,pseudotime=1:max(ld$pseudotime))
  }))
  ldg$type <- as.factor(ldg$Variable)
  
  apd <- pd.pop
  apd$ymin <- ldg$ymin
  apd$ymax <- ldg$ymax
  
  ### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  p <- ggplot(data=apd) + 
    geom_line(aes(x = pseudotime, y = expr, group = type, color = type), size = line.size, alpha = line.alpha) + 
    geom_ribbon(aes(x=pseudotime, ymin=ymin,ymax=ymax, fill=type), alpha=ribbon.alpha, size=ribbon.size, show.legend = F) +
    scale_color_brewer(palette = 'Dark2', direction = -1) +
    scale_fill_manual(values = brewer.pal(3,'Dark2')[2:1]) + ggtitle(gene)
  p <- p + theme_classic() +
    # ggtitle(paste0(sub(':.*','',gene),',adj.pvalue=', formatC(testobj$fdr[gene], format = "e", digits = 2))) +
    xlab(x.lab) + ylab(y.lab) + 
    theme(legend.spacing.y = unit(0.01, 'cm'), legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.1, "cm"), axis.title = element_text(size = axis.size), axis.text.x = element_text(angle = 45, hjust = 1, size = axis.text.size), axis.text.y = element_text(size = axis.text.size), plot.title = element_text(size = title.size), legend.position = legend.position) +
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) 
  if (axis.text.blank) {
    p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank())
  } else {
    p 
  }
}



