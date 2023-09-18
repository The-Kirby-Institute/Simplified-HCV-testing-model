# cost effectiveness frontier plot 

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
library("dampack")
library("kableExtra")
library("dplyr")
library("writexl")
library("ggrepel")
library("magick")
library("ggpmisc")
# Setup directories after setting working directory to source file 
# directory 

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

load(file.path(paste0(outputdt, "/" ,"CostBenfit.rda", sep ="")))
load(file.path(paste0(outputdt, "/" ,"CostBenfit_incre.rda", sep ="")))

# table numbers 
tab_costben <- dt_tab_char%>%filter(indicator %in% c("Lifetime cost (discounted, millions US$)",
                                                       "Lifetime QALY (discounted, thousands)"))%>%
  arrange(scenario, testing, population)%>%
  mutate(M_PI = paste0(round(Med, digits = 0), "(", round(q5,digits = 0), "-", round(q95,digits = 0), ")"))%>%
  select(population, scenario, testing, indicator, M_PI)


tab_incre <- dt_incre%>%filter(indicator %in% c("Incremental lifetime cost (discounted, millions US$)",
                                                "QALY gained (discounted)"))%>%
  arrange(scenario, testing, population)%>%
  mutate(M_PI = paste0(round(Med, digits = 0), "(", round(q5,digits = 0), "-", round(q95,digits = 0), ")"))%>%
  select(population, scenario, testing, indicator, M_PI)

write_xlsx(tab_costben,paste0(outputfig,"\table.xlsx"))
write_xlsx(tab_incre,paste0(outputfig,"\table_incre.xlsx"))


dt_tab_char <- dt_tab_char%>%mutate(front_lab = case_when(scenario == "PrEP coverage and HIV treatment coverage remains unchanged" ~ 
                                                            paste0("S1", "_", testing),
                                                          scenario == "HIV treatment coverages increase to 95% by 2030" ~
                                                            paste0("S3", "_", testing),
                                                          scenario == "PrEP coverage increases to 20% by 2030" ~
                                                            paste0("S2", "_", testing),
                                                          scenario == "PrEP coverage increases to 20% and HIV treatment increased to 95% by 2030" ~
                                                            paste0("S4", "_", testing)
                                                          ))
dt_tab_char$testing
dt_tab_char$front_lab <- gsub("Status Quo", "lab-based testing", dt_tab_char$front_lab)

dt_tab_all <- list()
for(i in unique(dt_tab_char$population)){ 
  dt_tab_all[[i]] <- dt_tab_char%>%filter(population == i & 
                                       indicator %in% c("Lifetime cost (discounted, millions US$)",
                                                        "Lifetime QALY (discounted)"))%>%
    select(front_lab, indicator, Med)%>%
    spread(indicator, Med)
  
  dt_tab_all[[i]] <- dt_tab_all[[i]]%>%mutate(Cost = `Lifetime cost (discounted, millions US$)`,
                                    effect = `Lifetime QALY (discounted)`)
  
  }

icer_hcv <- list()
for(i in unique(dt_tab_char$population)){
  icer_hcv[[i]] <- calculate_icers(cost=dt_tab_all[[i]]$Cost, 
                              effect=dt_tab_all[[i]]$effect, 
                              strategies=dt_tab_all[[i]]$front_lab)
  
  
}

icer_hcv$`Entire population of MSM`%>%
  kable() %>%
  kable_styling()
icer_hcv$`MSM who are HIV diagnosed (and on treatment)` %>%
  kable() %>%
  kable_styling()

icer_hcv$`HIV-negative MSM on PrEP` %>%
  kable() %>%
  kable_styling()
dt_tab_char$population
dt_tab_HIVD <- dt_tab_char%>%filter(population %in% c("MSM who are HIV diagnosed (and on treatment)") & 
                                     indicator %in% c("Lifetime cost (discounted, millions US$)",
                                                      "Lifetime QALY (discounted)",
                                                      ""))%>%
  select(front_lab, indicator, Med)%>%
  spread(indicator, Med)

dt_tab_HIVD <- dt_tab_HIVD%>%mutate(Cost = `Lifetime cost (discounted, millions US$)`,
                                  effect = `Lifetime QALY (discounted, thousands)`)


icer_hcv <- calculate_icers(cost=dt_tab_HIVD$Cost, 
                            effect=dt_tab_HIVD$effect, 
                            strategies=dt_tab_HIVD$front_lab)
icer_hcv %>%
  kable() %>%
  kable_styling()
unique(icer_hcv$`Entire population of MSM`$`Testing strategy`)
icer_hcv$`Entire population of MSM` <- icer_hcv$`Entire population of MSM`%>%
  mutate(`Testing strategy` = substr(Strategy, 4, nchar(Strategy)),
         `HIV prevention and care scenarios` = substr(Strategy, 1, 2))%>%
  mutate(`Testing strategy` = factor(`Testing strategy`, 
                                     levels = c("lab-based testing",
                                                "Point-of-care antibody testing",
                                                "Reflex RNA testing",
                                                "Point-of-care RNA testing"
                                                ),
                                     
                                     labels = c("Current practice of HCV testing",
                                                "Point-of-care antibody testing",
                                                "Reflex RNA testing",
                                                "Point-of-care RNA testing"     
                                                )))

# replacing inc_cost and inc_effect compared to S1_current practice 



base <- icer_hcv$`Entire population of MSM`%>%filter(Strategy == "S1_lab-based testing")


icer_hcv$`Entire population of MSM` <- icer_hcv$`Entire population of MSM`%>%
  mutate(Inc_Cost = Cost - base$Cost,
         Inc_Effect = Effect - base$Effect)



colpal <- scales::viridis_pal()(4)
show_col(colpal)
x <- icer_hcv$`Entire population of MSM`%>%
  ggplot(data = ., 
         aes(x = Inc_Effect, 
             y = Inc_Cost, shape = Strategy, colour = Strategy, group = Strategy)) + 
  geom_point(aes(color = Strategy, shape = Strategy), size = 2) + 
  geom_line(data = icer_hcv$`Entire population of MSM`%>%filter(Status == "ND" ),
            aes(x=Inc_Effect, y = Inc_Cost),color = "black")  +
  theme_Publication(base_size=14)  +   
  scale_shape_manual(name = "Scenarios",
                     labels = c("S1: Current practice of HCV testing", 
                                "S1: Point-of-care antibody testing", 
                                "S1: Reflex RNA testing",
                                "S1: Point-of-care RNA testing",
                                "S2: Current practice of HCV testing", 
                                "S2: Point-of-care antibody testing", 
                                "S2: Reflex RNA testing",
                                "S2: Point-of-care RNA testing",
                                "S3: Current practice of HCV testing", 
                                "S3: Point-of-care antibody testing", 
                                "S3: Reflex RNA testing",
                                "S3: Point-of-care RNA testing",
                                "S4: Current practice of HCV testing", 
                                "S4: Point-of-care antibody testing", 
                                "S4: Reflex RNA testing",
                                "S4: Point-of-care RNA testing"),
                     values=c(15,19,17,18,
                              15,19,17,18,
                              15,19,17,18,
                              15,19,17,18)) + 
  scale_colour_manual(name = "Scenarios",
                      labels = c("S1: Current practice of HCV testing", 
                                 "S1: Point-of-care antibody testing", 
                                 "S1: Reflex RNA testing",
                                 "S1: Point-of-care RNA testing",
                                 "S2: Current practice of HCV testing", 
                                 "S2: Point-of-care antibody testing", 
                                 "S2: Reflex RNA testing",
                                 "S2: Point-of-care RNA testing",
                                 "S3: Current practice of HCV testing", 
                                 "S3: Point-of-care antibody testing", 
                                 "S3: Reflex RNA testing",
                                 "S3: Point-of-care RNA testing",
                                 "S4: Current practice of HCV testing", 
                                 "S4: Point-of-care antibody testing", 
                                 "S4: Reflex RNA testing",
                                 "S4: Point-of-care RNA testing"),
                      values = c("#440154FF","#440154FF", "#440154FF","#440154FF",
                                 "#31688EFF","#31688EFF","#31688EFF","#31688EFF",
                                 "#35B779FF","#35B779FF","#35B779FF","#35B779FF",
                                 "#FDE725FF","#FDE725FF","#FDE725FF","#FDE725FF"
                                 )) +
  guides(fill=guide_legend(nrow=4)) + 
  theme(legend.text = element_text(size=10))

x 

x
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
xl <- g_legend(x)


# main plot: excluding S4 and explorating the frontier line. 
x_main <- icer_hcv$`Entire population of MSM`%>%
  filter(`Testing strategy` != "Point-of-care RNA testing")
# get the base case 
x_main[x_main$Inc_Cost==0, "Status"] <- "ND"  
x_main[x_main$Strategy=="S1_Point-of-care antibody testing", "ICER"] <-  
  x_main[x_main$Strategy=="S1_Point-of-care antibody testing", "Inc_Cost"]/
  x_main[x_main$Strategy=="S1_Point-of-care antibody testing", "Inc_Effect"]

exp_intercept <- c(5864953.8 - 4467034.4)/c(7615.584-7355.095)
# adding the explorated dot (invisible) 
x_main_exp <- x_main[x_main$Inc_Cost==0,]
# manuualy calculation exploration 

endp <- exp_intercept* (8000 - 7615.584) + 5864953.8 
x_main_exp$Inc_Cost <- endp
x_main_exp$Inc_Effect <- 8000


x_mainp <- ggplot(data = x_main, 
         aes(x = Inc_Effect, 
             y = Inc_Cost )) + 
  geom_point(aes(color = `HIV prevention and care scenarios`, shape = `Testing strategy`), size = 2) + 
  geom_line(data = icer_hcv$`Entire population of MSM`%>%filter(Status == "ND" ),
            aes(x=Inc_Effect, y = Inc_Cost),color = "black")  +
  geom_segment(x = 7615.584, y = 5864953.8, xend = 8000, yend = endp, 
               linetype = 2) + 
  geom_label_repel(data = x_main%>%filter(Status == "ND" & Inc_Cost!= 0),
                   aes(label=paste0 ("$",round(ICER, digits = 1),"/QALYs")), 
                   angle = 45, face = "bold", segment.size=0.5, vjust = 0.5, 
                   position = position_dodge(0.5))+
  theme_Publication(base_size=14) + 
  scale_x_continuous(limits = c(-1000,  9000), 
                     breaks = seq(-1000, 8000, 1000)) + 
  scale_y_continuous(limits = c(-2000000, 10000000), 
                     breaks = seq(-2000000, 8000000, 1000000),
                     labels = seq(-2000000, 8000000, 1000000)/1000000) + 
  scale_shape_manual(values=c(15,19,17,18)) + 
  scale_colour_manual(values=c("#440154FF","#31688EFF", 
                               "#35B779FF","#FDE725FF")) + 
  theme(legend.position = "") + 
  theme(axis.text.x = element_text(angle = 360, hjust = 0.5)) + 
  labs(x = "Incremental QALYs", y = "Incremental Cost (millions)")
  

# combine main plot and legend
lay <- lay <- rbind(c(1,1,1,1,1,1,1,1),
                    c(1,1,1,1,1,1,1,1),
                    c(1,1,1,1,1,1,1,1),
                    c(1,1,1,1,1,1,1,1),
                    c(2,2,2,2,2,2,2,2)) 
frontierplot_main <- grid.arrange(grobs = list(x_mainp, xl), layout_matrix = lay) 


frontierplot_main <- frontierplot_main + annotate(geom = "figure",
                                                x = 9,
                                                y = 3)

ggsave(path = outputfig, file="frontierplot_main.png", frontierplot_main, 
       height = 5, width = 11, dpi = 800)


x_POCRNA <- ggplot(data = icer_hcv$`Entire population of MSM`%>%
                     filter(`Testing strategy` == "Point-of-care RNA testing"), 
                  aes(x = Inc_Effect, 
                      y = Inc_Cost )) + 
  geom_point(aes(color = `HIV prevention and care scenarios`, shape = `Testing strategy`), size = 5) + 
  geom_line(data = icer_hcv$`Entire population of MSM`%>%filter(Status == "ND" ),
            aes(x=Inc_Effect, y = Inc_Cost),color = "black") +   
  theme_Publication(base_size=14) + 
  scale_x_continuous(limits = c(3000,  10000), 
                     breaks = seq(3000, 10000, 1000)) + 
  scale_y_continuous(limits = c(80000000, 140000000), 
                     breaks = seq(80000000, 140000000, 10000000),
                     labels = seq(80000000, 140000000, 10000000)/1000000) + 
  scale_shape_manual(values=c(18)) + 
  scale_colour_manual(values=c("#440154FF","#31688EFF", 
                               "#35B779FF","#FDE725FF")) + 
  theme(legend.position = "") + 
  theme(axis.text.x = element_text(angle = 360, hjust = 0.5)) + 
  labs(x = "Incremental QALYs", y = "")

x_POCRNA 
ggsave(path = outputfig, file="frontierplot_POCRNA.png", x_POCRNA , 
       height = 3
         , width = 6, dpi = 300)


ggsave(path = outputfig, file="frontierplot.png", x, height = 6, width = 15, 
       dpi = 800)

a <- plot(icer_hcv$`Entire population of MSM`, currency = "USD", 
          effect_units = "quality-adjusted life-years"#, label="all"
          ) + 
  theme_Publication(base_size=20) + 
  
  theme(axis.text.x = element_text(angle = 360, hjust = 0.5),
        text = element_text(size = 28))

ggsave(path = outputfig, file="frontierplot.png", a, height = 6, width = 9, dpi = 800)
