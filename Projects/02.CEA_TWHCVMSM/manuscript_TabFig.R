# this script generated the figures for the manuscript of CEAMSM 
rm(list = ls()) 

library(here)
here()

# Load useful libraries
library("dplyr")
library("tidyr")
library("ggthemes")
library("stringr")
library("viridis")
library("data.table")
library("formattable")
# Setup directories after setting working directory to source file 
# directory 

codep <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/TWHCV-model"

com_codeP <- "/Users/jjwu/Documents/Simplified-HCV-testing-model/03. Code"

epidatapath <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Taiwan-MSM-HCV-model"

project_codep <- "/Users/jjwu/Documents/Simplified-HCV-testing-model/Projects/02.CEA_TWHCVMSM"

#file path of "TWHCV-model" project
Rcode <- file.path(codep, '03. Code')

scenariopath <- file.path(epidatapath, '01. DATA/model input/Scenarios')

# save and cost data folder 
dtp <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/02.CEA_TWHCVMSM"

DataFolder <- file.path(dtp, "01. DATA/model input")

outputdt <- file.path(dtp, "02. Output/RDA")

outputfig <- file.path(dtp, "02. Output/Figs")

# source 

source(file.path(com_codeP, "Functions/plotFunctions.R"))

source(file.path(com_codeP, "Functions/plotOptions.R"))

source(file.path(com_codeP, "Functions/Scenarios_wrangling.R"))

source(file.path(com_codeP, "Functions/plotManuscript.R"))

source(file.path(project_codep, "AggregateRes.R"))

load(file.path(paste0(outputdt, "/" ,"measureOutcome.rda", sep ="")))

projectFile <- file.path(epidatapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile)

Currency_converter <- 14.07

TW_GDP2022 <- 32811

# color plate 
colpal <- scales::viridis_pal()(3)
Inci_pop$PrEPsame$population
####HCV incidence from 2022-2030 ####
# entire MSM population 
Inci_plot <- list()
tags <- c("A", "B", "C", "D")
for(i in seq_along(names(Inci))){ 
  
  Inci_plot[[i]] <- Inci[[i]]%>%
    mutate(year = (HCV$cabY + year -1),
           `Testing strategies` = factor(testing, labels = c("Current practice", 
                                                "Point-of-care antibody testing",
                                                "Reflex RNA testing",
                                                "Point-of-care RNA testing")))%>%
    ggplot(., aes(x = year, y = best)) + 
    geom_line(aes(colour = `Testing strategies`), linewidth = 1.5 ) + 
    labs(x = "Year", y = "HCV incidence (/1000)", tag = tags[i]) + 
    scale_color_viridis(discrete = TRUE, option =  "D") + 
    scale_x_continuous(limits = c(2022, 2030)) + 
    scale_y_continuous(limits = c(0, 15), breaks = seq(0,15,1))+
    theme_Publication(base_size = 20) +
    theme(legend.position = "none")
  
}
# number of infection and infection averted 
num_plot <- list()
for(i in seq_along(names(Inci))){ 
  num_plot[[i]] <- NumInf[[i]]%>%
    mutate(year = (HCV$cabY + year -1),
           `Testing strategies` = factor(testing, 
                                         labels = c("Current practice", 
                                                    "Point-of-care antibody testing",
                                                    "Reflex RNA testing",
                                                    "Point-of-care RNA testing")))%>%
    ggplot(., aes(x = year, y = best)) + 
    geom_line(aes(colour = `Testing strategies`), linewidth = 1.5 ) + 
    labs(x = "Year", y = "Number of new infections", tag = tags[i]) + 
    scale_color_viridis(discrete = TRUE, option =  "D") + 
    scale_x_continuous(limits = c(2022, 2030)) + 
    scale_y_continuous(limits = c(0, 5000), breaks = seq(0,5000,500))+
    theme_Publication(base_size = 20) +
    theme(legend.position = "none")
  
}


fig1x <- ggarrange(plotlist = Inci_plot, 
                                  common.legend = TRUE, legend="top",
                  nrow = 1, 
                  font.label = list(size = 20))+ 
  theme(legend.title = element_text(size = 20, face ="bold"),
        legend.text = element_text(size = 20, face ="bold")) 

fig1xx <- ggarrange(plotlist = num_plot,
                   nrow = 1, 
                   font.label = list(size = 20)) +
  theme(legend.title = element_text(size = 20, face ="bold"),
        legend.text = element_text(size = 20, face ="bold")) 


fig1 <- ggarrange(plotlist = list(fig1x, fig1xx), nrow = 2, 
                  common.legend = TRUE,  legend="top") + 
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
ggsave(file=file.path(outputfig, paste0("fig_inci_all.jpeg")), 
                      fig1, width = 16, height = 10, dpi = 300)



## sub pop 
Incisubpop_plot <- list()
tags <- c("A", "B", "C", "D")
for(i in seq_along(names(Inci_pop))){ 
  
  Incisubpop_plot[[i]] <- Inci_pop[[i]]%>%
    mutate(year = (HCV$cabY + year -1),
           `Testing strategies` = factor(testing, labels = c("Current practice", 
                                                             "Point-of-care antibody testing",
                                                             "Reflex RNA testing",
                                                             "Point-of-care RNA testing")))%>%
    filter(population =="HIV-PrEP")%>%
    ggplot(., aes(x = year, y = best)) + 
    geom_line(aes(colour = `Testing strategies`), linewidth = 1.5 ) + 
    labs(x = "Year", y = "HCV incidence (/1000)", tag = tags[i]) + 
    scale_color_viridis(discrete = TRUE, option =  "D") + 
    scale_x_continuous(limits = c(2022, 2030)) + 
    scale_y_continuous(limits = c(0, 6), breaks = seq(0,6,1))+
    theme_Publication(base_size = 20) +
    theme(legend.position = "none")
  
}

# number of infection and infection averted 
numsubpop_plot <- list()
for(i in seq_along(names(Inci))){ 
  numsubpop_plot[[i]] <- NumInf_pop[[i]]%>%
    mutate(year = (HCV$cabY + year -1),
           `Testing strategies` = factor(testing, 
                                         labels = c("Current practice", 
                                                    "Point-of-care antibody testing",
                                                    "Reflex RNA testing",
                                                    "Point-of-care RNA testing")))%>%
    filter(population == "HIV-PrEP")%>%
    ggplot(., aes(x = year, y = best)) + 
    geom_line(aes(colour = `Testing strategies`), linewidth = 1.5 ) + 
    labs(x = "Year", y = "Number of new infections", tag = tags[i]) + 
    scale_color_viridis(discrete = TRUE, option =  "D") + 
    scale_x_continuous(limits = c(2022, 2030)) + 
    scale_y_continuous(limits = c(0, 150), breaks = seq(0,150,10))+
    theme_Publication(base_size = 20) +
    theme(legend.position = "none")
  
}

fig1x_prep <- ggarrange(plotlist = Incisubpop_plot, 
                   common.legend = TRUE, legend="top",
                   nrow = 1, 
                   font.label = list(size = 20))+ 
  theme(legend.title = element_text(size = 20, face ="bold"),
        legend.text = element_text(size = 20, face ="bold")) 

fig1xx_prep <- ggarrange(plotlist = numsubpop_plot,
                    nrow = 1, 
                    font.label = list(size = 20)) +
  theme(legend.title = element_text(size = 20, face ="bold"),
        legend.text = element_text(size = 20, face ="bold")) 


fig1_prep <- ggarrange(plotlist = list(fig1x_prep, fig1xx_prep), nrow = 2, 
                  common.legend = TRUE,  legend="top") + 
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
ggsave(file=file.path(outputfig, paste0("fig_inci_PrEP.jpeg")), 
       fig1_prep, width = 16, height = 10, dpi = 300)


# MSM living with HIV 

Incisubpop_plot <- list()
tags <- c("A", "B", "C", "D")
for(i in seq_along(names(Inci_pop))){ 
  
  Incisubpop_plot[[i]] <- Inci_pop[[i]]%>%
    mutate(year = (HCV$cabY + year -1),
           `Testing strategies` = factor(testing, labels = c("Current practice", 
                                                             "Point-of-care antibody testing",
                                                             "Reflex RNA testing",
                                                             "Point-of-care RNA testing")))%>%
    filter(population =="HIV+d")%>%
    ggplot(., aes(x = year, y = best)) + 
    geom_line(aes(colour = `Testing strategies`), linewidth = 1.5 ) + 
    labs(x = "Year", y = "HCV incidence (/1000)", tag = tags[i]) + 
    scale_color_viridis(discrete = TRUE, option =  "D") + 
    scale_x_continuous(limits = c(2022, 2030)) + 
    scale_y_continuous(limits = c(0, 50), breaks = seq(0,50,10))+
    theme_Publication(base_size = 20) +
    theme(legend.position = "none")
  
}

# number of infection and infection averted 
numsubpop_plot <- list()
for(i in seq_along(names(Inci))){ 
  numsubpop_plot[[i]] <- NumInf_pop[[i]]%>%
    mutate(year = (HCV$cabY + year -1),
           `Testing strategies` = factor(testing, 
                                         labels = c("Current practice", 
                                                    "Point-of-care antibody testing",
                                                    "Reflex RNA testing",
                                                    "Point-of-care RNA testing")))%>%
    filter(population == "HIV+d")%>%
    ggplot(., aes(x = year, y = best)) + 
    geom_line(aes(colour = `Testing strategies`), linewidth = 1.5 ) + 
    labs(x = "Year", y = "Number of new infections", tag = tags[i]) + 
    scale_color_viridis(discrete = TRUE, option =  "D") + 
    scale_x_continuous(limits = c(2022, 2030)) + 
    scale_y_continuous(limits = c(0, 3000), breaks = seq(0,3000,500))+
    theme_Publication(base_size = 20) +
    theme(legend.position = "none")
  
}


fig1x_hivd <- ggarrange(plotlist = Incisubpop_plot, 
                        common.legend = TRUE, legend="top",
                        nrow = 1, 
                        font.label = list(size = 20))+ 
  theme(legend.title = element_text(size = 20, face ="bold"),
        legend.text = element_text(size = 20, face ="bold")) 

fig1xx_hivd <- ggarrange(plotlist = numsubpop_plot,
                         nrow = 1, 
                         font.label = list(size = 20)) +
  theme(legend.title = element_text(size = 20, face ="bold"),
        legend.text = element_text(size = 20, face ="bold")) 


fig1_hivd <- ggarrange(plotlist = list(fig1x_hivd , fig1xx_hivd ), nrow = 2, 
                       common.legend = TRUE,  legend="top") + 
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
ggsave(file=file.path(outputfig, paste0("fig_inci_HIVD.jpeg")), 
       fig1_hivd, width = 16, height = 10, dpi = 300)



####========================================================================####
test_order <- c("Status Quo", "POC_antibody", "Reflex_RNA", "POC_RNA")

Cost_plot <- list()

for(i in seq_along(names(HCVcost))){ 
  Cost_plot[["Status Quo"]][[i]] <- 
    HCVcost[[i]]$`Status Quo`%>%as.data.frame()%>%
    group_by(year, testing, Casca)%>%
    summarize(best = sum(best)/Currency_converter)%>%
    filter(year<= 2090)%>%
    mutate(`Cost categories` = factor(Casca, labels = c("Diagnosis", "Management", "Treatment")))%>%
    ggplot(., aes(x=year, y=best, fill=`Cost categories`)) + 
    geom_area() + 
    scale_fill_viridis(discrete = TRUE, option = "mako") + 
    scale_x_continuous(expand = c(0,0), limits = c(2022, 2090), 
                       breaks = c(2022, 2030, seq(2030, 2090,10))) + 
    scale_y_continuous(expand = c(0,0), limits = c(0, 10000000), 
                       breaks = seq(0, 10000000, 1000000),
                       labels = seq(0, 10000000, 1000000)/1000000) + 
    labs(x = "", y = "") +

    ggtitle("Current practice of HCV testing") +
    theme_Publication(base_size = 14) + 
    theme(legend.position = "right") + 
    guides(fill=guide_legend(direction='vertical', override.aes = list(size=5)))
  
  Cost_plot[["POC_antibody"]][[i]] <- 
    HCVcost[[i]]$POC_antibody%>%as.data.frame()%>%
    group_by(year, testing, Casca)%>%
    summarize(best = sum(best)/Currency_converter)%>%
    filter(year<= 2090)%>%
    mutate(`Cost categories` = factor(Casca, labels = c("Diagnosis", "Management", "Treatment")))%>%
    ggplot(., aes(x=year, y=best, fill=`Cost categories`)) + 
    geom_area() + 
    scale_fill_viridis(discrete = TRUE, option = "mako") + 
    scale_x_continuous(expand = c(0,0), limits = c(2022, 2090), 
                       breaks = c(2022, 2030, seq(2030, 2090,10))) + 
    scale_y_continuous(expand = c(0,0), limits = c(0, 10000000), 
                       breaks = seq(0, 10000000, 1000000),
                       labels = seq(0, 10000000, 1000000)/1000000) + 
    labs(x = "", y = "") +
    ggtitle("Point-of-care antibody testing") +
    theme_Publication(base_size = 14) + 
    theme(legend.position = "none")
  
  
  
  Cost_plot[["Reflex_RNA"]][[i]] <- 
    HCVcost[[i]]$Reflex_RNA%>%as.data.frame()%>%
    group_by(year, testing, Casca)%>%
    summarize(best = sum(best)/Currency_converter)%>%
    filter(year<= 2090)%>%
    mutate(`Cost categories` = factor(Casca, labels = c("Diagnosis", "Management", "Treatment")))%>%
    ggplot(., aes(x=year, y=best, fill=`Cost categories`)) + 
    geom_area() + 
    scale_fill_viridis(discrete = TRUE, option = "mako") + 
    scale_x_continuous(expand = c(0,0),limits = c(2022, 2090), 
                       breaks = c(2022, 2030, seq(2030, 2090,10))) + 
    scale_y_continuous(expand = c(0,0), limits = c(0, 10000000), 
                       breaks = seq(0, 10000000, 1000000),
                       labels = seq(0, 10000000, 1000000)/1000000) + 
    labs(x = "", y = "") +
    ggtitle("Reflex RNA testing") +
    theme_Publication(base_size = 14) + 
    theme(legend.position = "none")
  
  Cost_plot[["POC_RNA"]][[i]] <- 
    HCVcost[[i]]$POC_RNA%>%as.data.frame()%>%
    group_by(year, testing, Casca)%>%
    summarize(best = sum(best)/Currency_converter)%>%
    filter(year<= 2090)%>%
    mutate(`Cost categories` = factor(Casca, labels = c("Diagnosis", "Management", "Treatment")))%>%
    ggplot(., aes(x=year, y=best, fill=`Cost categories`)) + 
    geom_area() + 
    scale_fill_viridis(discrete = TRUE, option = "mako") + 
    scale_x_continuous(expand = c(0,0),limits = c(2022, 2090), 
                       breaks = c(2022, 2030, seq(2030, 2090,10))) + 
    scale_y_continuous(expand = c(0,0), limits = c(0, 10000000), 
                       breaks = seq(0, 10000000, 1000000),
                       labels = seq(0, 10000000, 1000000)/1000000) + 
    labs(x = "", y = "") +
    ggtitle("Point-of-care RNA testing") +
    theme_Publication(base_size = 14) + 
    theme(legend.position = "")
  
  
}


#### get legend function 
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

x_l <- g_legend(Cost_plot[["Status Quo"]][[1]])
Cost_plot[["Status Quo"]][[1]] <- Cost_plot[["Status Quo"]][[1]] + 
  theme(legend.position = "")
Cost_plot[["Status Quo"]][[2]] <- Cost_plot[["Status Quo"]][[2]] + 
  theme(legend.position = "")
Cost_plot[["Status Quo"]][[3]] <- Cost_plot[["Status Quo"]][[3]] + 
  theme(legend.position = "")
Cost_plot[["Status Quo"]][[4]] <- Cost_plot[["Status Quo"]][[4]] + 
  theme(legend.position = "")

# arrange beside ggarrange 
lay <- rbind(c(1,1,1,1,1,1, NA),
             c(2,2,2,2,2,2,5),
             c(3,3,3,3,3,3,5),
             c(4,4,4,4,4,4,NA)) 



xtt <- ggarrange(plotlist = list(Cost_plot[["Status Quo"]][[1]], 
                                 Cost_plot$POC_antibody[[1]],
                                 Cost_plot$Reflex_RNA[[1]],
                                 Cost_plot$POC_RNA[[1]]), labels = c(tags[1]),
                 nrow = 1, font.label = list(size = 14), 
                 common.legend = FALSE)

xtt2 <- ggarrange(plotlist = list(Cost_plot[["Status Quo"]][[2]], 
                                  Cost_plot$POC_antibody[[2]],
                                  Cost_plot$Reflex_RNA[[2]],
                                  Cost_plot$POC_RNA[[2]]), labels = c(tags[2]), 
                  font.label = list(size = 14), nrow = 1) 

xtt3 <- ggarrange(plotlist = list(Cost_plot[["Status Quo"]][[3]], 
                                  Cost_plot$POC_antibody[[3]],
                                  Cost_plot$Reflex_RNA[[3]],
                                  Cost_plot$POC_RNA[[3]]), labels = c(tags[3]), 
                  font.label = list(size = 14), nrow = 1)  

xtt4 <- ggarrange(plotlist = list(Cost_plot[["Status Quo"]][[4]], 
                                  Cost_plot$POC_antibody[[4]],
                                  Cost_plot$Reflex_RNA[[4]],
                                  Cost_plot$POC_RNA[[4]]), labels = c(tags[4]),
                  font.label = list(size = 14),
                  nrow = 1)  


xxx <- grid.arrange(grobs = list(xtt, xtt2, xtt3, xtt4, x_l), layout_matrix = lay)
# Annotate the figure by adding a common labels

fig3 <- annotate_figure(xxx,
                       bottom = text_grob("Year",,
                                          hjust = 1, vjust = -1, x = 0.45,  
                                          size = 14, face ="bold"),
                       left = text_grob("Cost (millions, US$)",  rot = 90, 
                                        size = 14, face = "bold")
)


ggsave(file=file.path(outputfig, paste0("fig_cost_all.jpeg")), 
       fig3, width = 14, height = 10, 
       dpi = 300)



# number output for results 
Inci_num <- list()
for(i in names(Inci)){ 
  Inci_num[[i]] <- Inci[[i]]%>%
    mutate(year = (HCV$cabY + year -1))%>%
    filter(year == 2030)%>%select(year, population, testing, best, q5, q95)
  }

Inci_num <- bind_rows(Inci_num, .id = "Sec")
  do.call("rbind", Inci_num, by = "Sce")
