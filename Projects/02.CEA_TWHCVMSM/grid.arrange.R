# loading library 
# if not exisiting, install.packages("PACKAGE NAME")
library(cowplot)
library(ggplot2)
library(gridExtra)
library(ggpubr)

heatmap_A <- ggplot(XXX) 
heatmap_B <- ggplot(XXX) 


legend <- cowplot::get_legend(heatmap_A) # get common lenged 

heatmap_A= heatmap_A + ggplot2::theme(legend.position = 'none') # remove lenged from the plot 
heatmap_B= heatmap_B + ggplot2::theme(legend.position = 'none') # remove lenged from the plot 

# arrange plot common axis 
gridExtra::grid.arrange(grobs = list(heatmap_A, heatmap_B, , , ,legend),
                        ncol = 3,  # number of plot per column 
                        left = ggpubr::text_grob("Your y-Axis", 
                                                 rot = 90, 
                                                 vjust = 1),
                        bottom = ggpubr::text_grob("your X-Axis", 
                                                   rot = 0, 
                                                   vjust = 0)
