plotOpts <- theme_bw() + theme(plot.title = element_text(face = "bold",
                                                         size = rel(1.2), 
                                                         hjust = 0.5),
                               text = element_text(face = "bold",size=14,
                                                   colour="black"),
                               axis.text.x = element_text(face = "bold", 
                                                          angle = 90,
                                                          size=14,
                                                          colour="black", 
                                                          hjust = 0, vjust = 0.5),
                               axis.text.y = element_text(face = "bold",
                                                          size=14,
                                                          colour="black"),
                               axis.line.x = element_line(colour="black"),
                               axis.line.y = element_line(colour="black"),
                               axis.ticks = element_line(colour="black"),
                               legend.background = element_rect(),
                               legend.key = element_blank(),
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(), 
                               panel.background = element_blank(), 
                               panel.border = element_blank(),
                               axis.line = element_line(colour = "black"),
                               strip.background = element_blank(),
                               plot.margin =unit(c(2.5,2.5,2.5,2.5),"mm"),
                               strip.text = element_text(face="bold")) +
  theme(plot.title = element_text(hjust = 0.5))





theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size) +
    theme(plot.title = element_text(face = "bold",
                                      size = rel(0.8), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle=90, vjust = 0.5,hjust = -1,
                                       face = "bold",size = rel(1)), 
            axis.text.y = element_text(
                                     face = "bold",size = rel(1)),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = rel(1)),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="bold", size= rel(1)),
            #plot.margin =unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


theme_Publication_facet <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size) + 
      theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border     = element_rect(fill = NA, colour = "black", 
                                            size = NA),
            panel.grid       = element_line(colour = "black"),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle=90, vjust = 0.5,hjust = -1,
                                       face = "bold",size = rel(1)), 
            axis.text.y = element_text(
              face = "bold",size = rel(1)),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = rel(1)),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="bold", size= rel(1)),
            plot.margin = unit(c(10,5,5,5),"mm"),
            strip.background = element_blank(),
            strip.text = element_text(face="bold")
            
      ))
  
}
 # 
#strip.background = element_blank(),
#strip.text = element_text(face="bold")