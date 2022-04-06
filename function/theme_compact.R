theme_compact <- function(){ 
  font <- "sans"   #assign font family up front
  
  theme_classic() %+replace%    #replace elements we want to change
    
    theme(
      strip.background = element_rect(color='white',fill='white'),
      legend.spacing.x =unit(0,'cm'),
      
      strip.text = element_text(
        family = font,
        size = 10,
        color = 'black'),
      
      strip.text.x = element_text(
        family = font,
        size = 10,
        color = 'black'),
      
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      # axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 10,                #set font size
        face = 'bold',            #bold typeface
        color = 'black',
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        color = 'black',
        size = 10),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        color = 'black',
        size = 10,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        color = 'black',
        size = 10),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        color = 'black',
        size = 8),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10))
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    ) 
  # guides(color = guide_legend(byrow = TRUE,override.aes = list(size=2,alpha=1)))
}

