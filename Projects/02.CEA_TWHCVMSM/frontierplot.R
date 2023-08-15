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
dt_tab_char$population
dt_tab_HIVD <- dt_tab_char%>%filter(population %in% c("MSM who are HIV diagnosed (and on treatment)") & 
                                     indicator %in% c("Lifetime cost (discounted, millions US$)",
                                                      "Lifetime QALY (discounted, thousands)",
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

icer_hcv$`Entire population of MSM`$Status

x <- ggplot(data = icer_hcv$`Entire population of MSM`, aes(x = Effect, y = Cost, shape = Status)) + 
  geom_point(aes(shape = Status)) + 
  geom_line(data = icer_hcv$`Entire population of MSM`%>%filter(Status == "ND"),
            aes(x=Effect,y=Cost),color = "black") + 
  geom_label_repel(data = icer_hcv$`Entire population of MSM`%>%filter(Status == "ND"),
    aes(label=Strategy), 
    nudge_x = 0.25, nudge_y = 1) +
  theme_Publication(base_size=20) + 
  scale_y_continuous(limits = c(50000000, 150000000), 
                     breaks = seq(50000000, 150000000, 10000000),
                     labels = seq(50, 150, 10)) +
  labs(y = "Cost, millions (US$)", x = "QALYs gained")
  
ggsave(path = outputfig, file="frontierplot.png", x, height = 6, width = 9, dpi = 800)
a <- plot(icer_hcv$`Entire population of MSM`, currency = "USD", 
          effect_units = "quality-adjusted life-years"#, label="all"
          ) + 
  theme_Publication(base_size=20) + 
  
  theme(axis.text.x = element_text(angle = 360, hjust = 0.5),
        text = element_text(size = 28))

ggsave(path = outputfig, file="frontierplot.png", a, height = 6, width = 9, dpi = 800)
