plotQuickHeatmap <- function(heatmap_plots, Max.gene = NULL, heatmap.x = 0.4, heatmap.y = 0.5, heatmap.width = 0.8, heatmap.height = 1.0, dend.x = 0.90, dend.y = 0.44, dend.width = 0.2, dend.height = 0.98){
  library(grid)
  heatmap.plot = heatmap_plots[[1]]
  dendro.plot = heatmap_plots[[2]]
  # All together
  grid.newpage()
  # print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
  # print(dendro.plot, vp = viewport(x = 0.90, y = 0.44, width = 0.2, height = 0.98))
  print(heatmap.plot, vp = viewport(x = heatmap.x, y = heatmap.y, width = heatmap.width, height = heatmap.height))
  print(dendro.plot, vp = viewport(x = dend.x, y = dend.y, width = dend.width, height = dend.height))
}
  

