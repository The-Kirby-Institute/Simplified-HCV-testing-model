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
# Setup directories after setting working directory to source file 
# directory 

#file path of "TWHCV-model" project
codepath <- file.path(here() %>% dirname()%>%dirname(), '03. Code/Functions')

DataFolder <- file.path(here(), "01. DATA/model input")

# Rda file path 
# load the .rda file of base estimate 
rdapath <- file.path(here()%>%dirname()%>%dirname(), "Taiwan-MSM-HCV-model")

# output file path 
# dir.create("02. Output") # create subdircetory 
# dir.create("02. Output/RDA") 
# dir.create("02. Output/Results")
# dir.create("02. Output/Figs")
outputdt <- here("02. Output/RDA")

outputfig <- here("02. Output/Figs")

# source 
source(file.path(codepath, "plotOptions.R"))

source(file.path(codepath, "plotFunctions.R"))

source(file.path(codepath, "plotManuscript.R"))

source(file.path(here(),"AggregateRes.R"))

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
                                                        "Lifetime QALY (discounted, thousands)"))%>%
    select(front_lab, indicator, Med)%>%
    spread(indicator, Med)
  
  dt_tab_all[[i]] <- dt_tab_all[[i]]%>%mutate(Cost = `Lifetime cost (discounted, millions US$)`,
                                    effect = `Lifetime QALY (discounted, thousands)`)
  
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

dt_tab_HIVD <- list()
icer_hcv <- list()
for(i in unique(dt_tab_char$population)){ 
  dt_tab_HIVD[[i]] <- dt_tab_char%>%filter(population == i & 
                                        indicator %in% c("Lifetime cost (discounted, millions US$)",
                                                         "Lifetime QALY (discounted, thousands)",
                                                         ""))%>%
    select(front_lab, indicator, Med)%>%
    spread(indicator, Med)%>%mutate(Cost = `Lifetime cost (discounted, millions US$)`,
                                    effect = `Lifetime QALY (discounted, thousands)`)
  
  icer_hcv[[i]] <- calculate_icers(cost=dt_tab_HIVD[[i]]$Cost, 
                              effect=dt_tab_HIVD[[i]]$effect, 
                              strategies=dt_tab_HIVD[[i]]$front_lab)
  
  
  }

# function for ceiling and floor 
library(plyr)

dtlim_cost <- list()
dtlim_effect <- list()
for(i in names(icer_hcv)){ 
  dtlim_cost[[i]] <- c(round_any(max(icer_hcv[[i]]$Cost, na.rm = T), 1000000, f = ceiling), 
                       round_any(min(icer_hcv[[i]]$Cost, na.rm = T), 1000000))
  dtlim_effect[[i]] <- c(floor_dec(min(icer_hcv[[i]]$Effect, na.rm = T), 100000), 
                         ceiling_dec(max(icer_hcv[[i]]$Effect, na.rm = T), 100000))
  }
icer_hcv$`HIV-negative MSM on PrEP`
x <- ggplot(data = icer_hcv$`Entire population of MSM`, 
            aes(y = Effect, x = Cost, shape = Status)) + 
  geom_point(aes(shape = Status)) + 
  geom_line(data = icer_hcv$`Entire population of MSM`%>%filter(Status == "ND"),
            aes(y=Effect,x=Cost),color = "black") + 
  geom_label_repel(data = icer_hcv$`Entire population of MSM`%>%filter(Status == "ND"),
    aes(label=Strategy), 
    nudge_x = 0.25, nudge_y = 1) +
  theme_Publication(base_size=20) + 
    scale_x_continuous(limits = round(c(dtlim_cost[[1]]), 1), 
                     breaks = seq(dtlim_cost[[1]][1], dtlim_cost[[1]][2], 5000000),
                     labels = round(seq(dtlim_cost[[1]][1]/1000000, 
                                        dtlim_cost[[1]][2]/1000000, 5000000/1000000),1)) +
  labs(x = "Cost, millions (US$)", y = "QALYs gained") 
  
ggsave(path = outputfig, file="frontierplot.png", x, height = 6, width = 9, dpi = 800)
a <- plot(icer_hcv$`Entire population of MSM`, currency = "USD", 
          effect_units = "quality-adjusted life-years"#, label="all"
          ) + 
  theme_Publication(base_size=20) + 
  
  theme(axis.text.x = element_text(angle = 360, hjust = 0.5),
        text = element_text(size = 28))

ggsave(path = outputfig, file="frontierplot.png", a, height = 6, width = 9, dpi = 800)
