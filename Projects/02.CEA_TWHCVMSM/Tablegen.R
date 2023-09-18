#### This script is for  tables and figures generartion 

## import data 

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
library("grid")
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

# name required packages
list.of.packages <- c("gapminder", "gt", "tidyverse")

# install required packages, if necessary, and load them ----
{
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages, require, character.only = TRUE)
}

# table default settings 
n = 0
c_col = c("#1e3048", "#274060", "#2f5375", "#4073a0", "#5088b9")
c_col_light_blue = c("#edf2fb", "#e2eafc", "#d7e3fc", "#ccdbfd", "#c1d3fe")
c_container_width = px(800)
c_table_width = px(650)
c_rn = 30
c_save = TRUE
c_format = "html"


colset <- c(paste0("set", seq(1, HCV$numberSamples,1)))
####preparing the dataset for table generation ####
##### each scenario #####
dt_name <- c("PrEPsame","PrEP", "HIVD", "PrEPnHIVD")



num_table <- list()



# cumulative HCV new infections (2022-2030) in each scenario 
for(i in dt_name){ 
  # cumulative new infections for subpops
  num_table[["cumu_newinf"]][[i]] <- NumInf_pop[[i]]%>%                  
    mutate(year = year + HCV$cabY - 1)%>%
    filter(year>=2022 & year<=2030)%>%
    group_by(population, testing)%>%                                            # grouping by population and testing for cumsum 
    mutate(across(c(best, colset), ~cumsum(.x)))%>%
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca = NULL)%>%   # get the aggegrated value: min, max..., etc.
    filter(population %in% c("HIV-PrEP", "HIV+d"))%>%ungroup()                                        # get rid of testing scneairo we don't need for this project
  
  # cumulative reinfection for subpops 
  num_table[["cumu_newreinf"]][[i]] <- NumreInf_pop[[i]]%>%                  
    mutate(year = year + HCV$cabY - 1)%>%
    filter(year>=2022 & year<=2030)%>%
    group_by(population, testing)%>%                                            # grouping by population and testing for cumsum 
    mutate(across(c(best, colset), ~cumsum(.x)))%>%
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca = NULL)%>%   # get the aggegrated value: min, max..., etc.
    filter(population %in% c("HIV-PrEP", "HIV+d"))%>%ungroup() 
  
  # cumulative reinfection for subpops 
  num_table[["cumu_newpriinf"]][[i]] <- NumpriInf_pop[[i]]%>%                  
    mutate(year = year + HCV$cabY - 1)%>%
    filter(year>=2022 & year<=2030)%>%
    group_by(population, testing)%>%                                            # grouping by population and testing for cumsum 
    mutate(across(c(best,colset), ~cumsum(.x)))%>%
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca = NULL)%>%   # get the aggegrated value: min, max..., etc.
    filter(population %in% c("HIV-PrEP", "HIV+d"))%>%ungroup() 
  
  # cumulative life time total cost for subpops 
  num_table[["cumu_lifecost"]][[i]] <- HCVcost_ScenCum_dis[[i]]%>%                  
    filter(year<=2090)%>%                                                       # assuming lifetime: 2090
    group_by(year, population, testing)%>%                                      # grouping by population and testing for summarize the cost of different stages 
    summarize(across(c(best, colset), ~sum(.x))/Currency_converter)%>%
    ungroup()%>%group_by(population, testing)%>%
    filter(population %in% c("HIV-PrEP", "HIV+d"))%>%ungroup()%>%                                     # get rid of testing scneairo we don't need for this project
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca = NULL)
  
  # cumulative life time cost of stages for subpops 
  num_table[["cumu_lifecost_stage"]][[i]] <- HCVcost_ScenCum_dis[[i]]%>%                  
    filter(year<=2090)%>%                                                       
    group_by(population, testing, Casca)%>%                                    
    mutate(across(c(best, colset), ~.x/Currency_converter))%>%
    filter(population %in% c("HIV-PrEP", "HIV+d"))%>%ungroup()%>%                                     
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca = "y")
  
  # cumulative QALYs for subpops 
  num_table[["cumu_QALY"]][[i]] <- HCVQALYCum[[i]]%>%                  
    filter(year<=2090)%>%                                                       
    filter(population %in% c("HIV-PrEP", "HIV+d"))%>%ungroup()%>%                                     
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca =  NULL)
  
}



# cumulative HCV new infections (2022-2030) for all pops 
num_table_all <- list()

for(i in dt_name){ 
  # cumulative new infections for overall pop
  num_table_all[["cumu_newinf"]][[i]] <- NumInf[[i]]%>%                  
    mutate(year = year + HCV$cabY - 1)%>%
    filter(year>=2022 & year<=2030)%>%
    group_by(population, testing)%>%                                            # grouping by population and testing for cumsum 
    mutate(across(c(best, colset), ~cumsum(.x)))%>%
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca = NULL)%>%
    ungroup()  
  
  
  num_table_all[["cumu_newreinf"]][[i]] <- NumreInf[[i]]%>%                  
    mutate(year = year + HCV$cabY - 1)%>%
    filter(year>=2022 & year<=2030)%>%
    group_by(population, testing)%>%                                            # grouping by population and testing for cumsum 
    mutate(across(c(best, colset), ~cumsum(.x)))%>%
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca = NULL)%>%
    ungroup() 
  
  
  num_table_all[["cumu_newpriinf"]][[i]] <- NumpriInf[[i]]%>%                  
    mutate(year = year + HCV$cabY - 1)%>%
    filter(year>=2022 & year<=2030)%>%
    group_by(population, testing)%>%                                            # grouping by population and testing for cumsum 
    mutate(across(c(best, colset), ~cumsum(.x)))%>%
    popRange(., target_year = NULL, pop = NULL, test = NULL, Casca = NULL)%>%
    ungroup()

  num_table_all[["cumu_lifecost"]][[i]] <- HCVcost_ScenCum_dis[[i]]%>%                  
    filter(year<=2090)%>%                                                       # assuming lifetime: 2090
    group_by(year, testing)%>%                                                  # grouping by population and testing for summarize the cost of different stages 
    summarize(across(c(best, colset), ~sum(.x)/Currency_converter))%>%
    ungroup()%>%group_by(testing)%>%ungroup()%>%                                     # get rid of testing scneairo we don't need for this project
    popRange(., target_year = NULL, pop = "all", test = NULL, Casca = NULL)%>%
    mutate(population = "All")
  
  
  
  num_table_all[["cumu_lifecost_stage"]][[i]] <- HCVcost_ScenCum_dis[[i]]%>%                  
    filter(year<=2090)%>%                                                       
    group_by(year, testing, Casca)%>% 
    summarize(across(c(best, colset), ~sum(.x)/Currency_converter))%>%
    ungroup()%>%group_by(testing, Casca)%>%ungroup()%>%                                     
    popRange(., target_year = NULL, pop = "all", test = NULL, Casca = "y")%>%
    mutate(population = "All")
  
  
  num_table_all[["cumu_QALY"]][[i]] <- HCVQALYCum[[i]]%>%                  
    filter(year<=2090)%>%
    group_by(year, testing)%>%
    summarize(across(c(best, colset), ~sum(.x)))%>%
    ungroup()%>%                                     
    popRange(., target_year = NULL, pop = "all", test = NULL, Casca =  NULL)%>%
    mutate(population = "All")
  }


# combine all pops and subpops 
tt <- list()
for(nam in names(num_table)){
  for(m in names(num_table[[1]])){ 
    
    tt[[nam]][[m]] <- rbind(num_table_all[[nam]][[m]], num_table[[nam]][[m]])%>%
      mutate(scenario = m,
             indicator = nam)
    }
}


# making format consistent 
for(m in names(num_table[[1]])){ 
  
  tt$cumu_newinf[[m]] <- tt$cumu_newinf[[m]]%>%filter(year ==2030)%>%
    mutate(Casca = "total")%>%
    select(year, population, testing, scenario, indicator, Casca, best, colset, 
           min, max, Med, Mu, q5, q25, q75, q95)
  
  tt$cumu_newreinf[[m]] <- tt$cumu_newreinf[[m]]%>%filter(year ==2030)%>%
    mutate(Casca = "total")%>%
    select(year, population, testing, scenario, indicator, Casca, best, colset, 
           min, max, Med, Mu, q5, q25, q75, q95)
  
  tt$cumu_newpriinf[[m]] <- tt$cumu_newpriinf[[m]]%>%filter(year == 2030)%>%
    mutate(Casca = "total")%>%
    select(year, population, testing, scenario, indicator, Casca, best, colset, 
           min, max, Med, Mu, q5, q25, q75, q95)
    
  
  tt$cumu_lifecost[[m]] <- tt$cumu_lifecost[[m]]%>%filter(year == 2090)%>%
    mutate(Casca = "total")%>%
    select(year, population, testing, scenario, indicator, Casca, best, colset, 
           min, max, Med, Mu, q5, q25, q75, q95)
  
  tt$cumu_lifecost_stage[[m]] <- tt$cumu_lifecost_stage[[m]]%>%
    filter(year == 2090)%>%
    select(year, population, testing, scenario, indicator, Casca, best, colset, 
           min, max, Med, Mu, q5, q25, q75, q95)
  
  tt$cumu_QALY[[m]] <- tt$cumu_QALY[[m]]%>%filter(year == 2090)%>%
    mutate(Casca = "total")%>%
    select(year, population, testing, scenario, indicator, Casca, best, colset, 
           min, max, Med, Mu, q5, q25, q75, q95)
  }


# subset dtframe 
dt_tab <- list()

# combind as a dataframe for generating table 
dt_tab <- lapply(tt, function(x) rbindlist(x))%>%rbindlist(.)

dt_tab <- dt_tab%>%arrange(indicator, scenario, testing, population)

dt_tab <- dt_tab%>%mutate(testing = 
                            case_when(testing == "Status Quo" ~ "Status Quo",
                                      testing == "StatusQuo" ~ "Status Quo",
                                      testing == "POC_antibody" ~ "POC_antibody",
                                      testing == "Reflex_RNA" ~ "Reflex_RNA",
                                      testing == "POC_RNA" ~ "POC_RNA"))

# cost stage # 

cost_S <- dt_tab%>%filter(population == "All" & indicator %in% c("cumu_lifecost_stage", 
                                                                   "cumu_lifecost"))%>%
  select(population, testing, scenario, indicator, Casca, Med)%>%
  arrange(population, testing, scenario, Casca)%>%
  mutate(SS = paste0(indicator, "_", Casca))%>%select(!c(Casca, indicator))%>%
  spread(SS, Med)%>%
  mutate(PER = round(cumu_lifecost_stage_diag/cumu_lifecost_total*100, digits = 1))%>%
  arrange(cumu_lifecost_total)

cost_S <- dt_tab%>%filter(population == "All" & indicator %in% c("cumu_lifecost_stage", 
                                                                 "cumu_lifecost"))%>%
  select(population, testing, scenario, indicator, Casca, Med)%>%
  arrange(population, testing, scenario, Casca)

cost_S$cumu_lifecost_total

# tidy format and labels 
dt_tab_char <- dt_tab%>%
  mutate(testing = factor(testing, 
                          levels = c(unique(dt_tab$testing)), 
                          label = c("Status Quo", 
                                    "Point-of-care antibody testing",
                                    "Reflex RNA testing",
                                    "Point-of-care RNA testing"
                                    )),
         scenario = factor(scenario, 
                           levels = c("PrEPsame", "PrEP", "HIVD", "PrEPnHIVD"),
                           label = c("PrEP coverage and HIV treatment coverage remains unchanged",
                                     "PrEP coverage increases to 20% by 2030",
                                     "HIV treatment coverages increase to 95% by 2030",
                                     "PrEP coverage increases to 20% and HIV treatment increased to 95% by 2030")),
         
         population = factor(population, 
                             levels = c("All", "HIV-PrEP", "HIV+d"),
                             label = c("Entire population of MSM", 
                                       "HIV-negative MSM on PrEP", 
                                       "MSM who are HIV diagnosed (and on treatment)")),
         
         indicator = factor(indicator, 
                            levels = c( "cumu_newinf", "cumu_newpriinf", "cumu_newreinf" ,"cumu_lifecost", "cumu_QALY" ),
                            label = c("New HCV infections (2022-2030)",
                                      "New HCV primary infections (2022-2030)",
                                      "New HCV reinfections (2022-2030)",
                                      "Lifetime cost (discounted, millions US$)",
                                      "Lifetime QALY (discounted)"
                                      #,
                                      #"Incremental cost effectiveness (Cost per QALY gained (US$))"
                                      ))
         )
save(dt_tab_char,
     file = file.path(outputdt, "CostBenfit.rda"))
# output table for lifetime cost and QALY
x <- dt_tab_char %>% 
  select(scenario, indicator, population, testing, 
         Med, q5, q95)%>%
  filter(indicator!="cumu_lifecost_stage")%>%
  mutate(Med = case_when(indicator %in% c("Lifetime cost (discounted, millions US$)")~ round(Med/1000000, digits = 1),
                          indicator %in% c("Lifetime QALY (discounted)")~ round(Med, digits = 0),
                          indicator %in% c("New HCV infections (2022-2030)", 
                                           "New HCV primary infections (2022-2030)",
                                           "New HCV reinfections (2022-2030)")~ round(Med, digits =  -1),
                          .default = Med
                          ),
         q5= case_when(indicator %in% c("Lifetime cost (discounted, millions US$)")~ round(q5/1000000, digits = 1),
                          indicator %in% c("Lifetime QALY (discounted)") ~ round(q5, digits = 0),
                       indicator %in% c("New HCV infections (2022-2030)",
                                        "New HCV primary infections (2022-2030)",
                                        "New HCV reinfections (2022-2030)")~ round(q5, digits =  -1),
                       .default = q5
                       ),
         q95 = case_when(indicator %in% c("Lifetime cost (discounted, millions US$)")~ round(q95/1000000, digits = 1),
                          indicator %in% c("Lifetime QALY (discounted)")~ round(q95,digits = 0),
                         indicator %in% c("New HCV infections (2022-2030)",
                                          "New HCV primary infections (2022-2030)",
                                          "New HCV reinfections (2022-2030)")~ round(q95, digits =  -1),
                         .default = q95))%>%
  mutate(Med = formatC(Med,  format = "fg", big.mark = ","),
                      q5 = formatC(q5,  format = "fg", big.mark = ","),
                      q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
  
  mutate(Med = case_when(indicator %in% c("Lifetime cost (discounted, millions US$)", 
                                           "Incremental cost effectiveness (Cost per QALY gained (US$))")~ paste0("$", as.character(Med)),
                          .default = as.character(Med)),
         q5 = case_when(indicator %in% c("Lifetime cost (discounted, millions US$)", 
                                           "Incremental cost effectiveness (Cost per QALY gained (US$))")~ paste0("$", as.character(q5)),
                          .default = as.character(q5)),
         q95 = case_when(indicator %in% c("Lifetime cost (discounted, millions US$)", 
                                           "Incremental cost effectiveness (Cost per QALY gained (US$))")~ paste0("$", as.character(q95)),
                          .default = as.character(q95)))%>%
  mutate(vv = paste0(Med, "\n", "(", q5, "-", q95, ")"))%>%
  select(-c(Med, q5, q95))%>%
  # reorder columns to present in the right order
  arrange(scenario, indicator)%>%
  tidyr::pivot_wider(names_from = c("indicator","population"), 
                     values_from = c("vv"), 
                     names_glue = "{indicator}X{population}")
  # grouping rows
x%>%gt(groupname_col = "scenario",
     rowname_col = "testing")%>%
  tab_spanner_delim(delim = "X")%>%
  gtsave(., file = file.path(outputfig, paste0("mainRes_0913.docx")))


#### Incremental table  #### 
 
Iavert_bench <- list()
 
Iavert_scen <- list()

prIavert_bench <- list()

prIavert_scen <- list()

reIavert_bench <- list()

reIavert_scen <- list()
 
cost_incre_bench <- list()
 
cost_incre_scen <- list()
 
cost_stage_incre_bench <- list()
 
cost_stage_incre_scen <- list()
 
qaly_incre_bench <- list()
 
qaly_incre_scen <- list()
 
# rename "Status quo" 

for(ind in names(tt)){ 
  for(i in dt_name){
    tt[[ind]][[i]] <- tt[[ind]][[i]]%>%
      mutate(testing = 
               case_when(testing == "Status Quo" ~ "Status Quo",
                         testing == "StatusQuo" ~ "Status Quo",
                         testing == "POC_antibody" ~ "POC_antibody",
                         testing == "Reflex_RNA" ~ "Reflex_RNA",
                         testing == "POC_RNA" ~ "POC_RNA"))
    
  }
}


t_lab <- unique(tt$cumu_newinf$PrEP$testing) # remove status quo
 
t_lab <- t_lab[t_lab!= "Status Quo"] 
for(i in dt_name){ 
   # bench: status quo
   Iavert_bench[[i]] <- tt$cumu_newinf[[i]]%>%
     filter(testing == "Status Quo")%>%
     select(year, testing, population, scenario, indicator, Casca, best,
            paste0("set", seq(1,HCV$numberSamples,1)))%>%
     arrange(year, population)
   
   prIavert_bench[[i]] <- tt$cumu_newpriinf[[i]]%>%
     filter(testing == "Status Quo")%>%
     select(year, testing, population, scenario, indicator, Casca, best,
            paste0("set", seq(1,HCV$numberSamples,1)))%>%
     arrange(year, population)
   
   reIavert_bench[[i]] <- tt$cumu_newreinf[[i]]%>%
     filter(testing == "Status Quo")%>%
     select(year, testing, population, scenario, indicator, Casca, best,
            paste0("set", seq(1,HCV$numberSamples,1)))%>%
     arrange(year, population)
   
   cost_incre_bench[[i]] <- tt$cumu_lifecost[[i]]%>%
     filter(testing == "Status Quo")%>%
     select(year, testing, population, scenario, indicator, Casca, best,
            paste0("set", seq(1,HCV$numberSamples,1)))%>%
     arrange(year, population)
   
   cost_stage_incre_bench[[i]] <- tt$cumu_lifecost_stage[[i]]%>%
     filter(testing == "Status Quo")%>%
     select(year, testing, population, scenario, indicator, Casca, best,
            paste0("set", seq(1,HCV$numberSamples,1)))%>%
     arrange(year, population)
   
   qaly_incre_bench[[i]] <- tt$cumu_QALY[[i]]%>%
     filter(testing == "Status Quo")%>%
     select(year, testing, population, scenario, indicator, Casca, best,
            paste0("set", seq(1,HCV$numberSamples,1)))%>%
     arrange(year, population)
   
   for(t in t_lab){ 
     
     Iavert_scen[[i]][[t]] <-  tt$cumu_newinf[[i]]%>%
       filter(testing == t)%>%
       select(year, testing, population, scenario, indicator, Casca, best,
              paste0("set", seq(1,HCV$numberSamples,1)))%>%
       arrange(year, population)
     
     prIavert_scen[[i]][[t]] <-  tt$cumu_newpriinf[[i]]%>%
       filter(testing == t)%>%
       select(year, testing, population, scenario, indicator, Casca, best,
              paste0("set", seq(1,HCV$numberSamples,1)))%>%
       arrange(year, population)
     
     reIavert_scen[[i]][[t]] <-  tt$cumu_newreinf[[i]]%>%
       filter(testing == t)%>%
       select(year, testing, population, scenario, indicator, Casca, best,
              paste0("set", seq(1,HCV$numberSamples,1)))%>%
       arrange(year, population)
     
     cost_incre_scen[[i]][[t]] <- tt$cumu_lifecost[[i]]%>%
       filter(testing == t)%>%
       select(year, testing, population, scenario, indicator, Casca, best,
              paste0("set", seq(1,HCV$numberSamples,1)))%>%
       arrange(year, population)
     
     cost_stage_incre_scen[[i]][[t]] <- tt$cumu_lifecost_stage[[i]]%>%
       filter(testing == t)%>%
       select(year, testing, population, scenario, indicator, Casca, best,
              paste0("set", seq(1,HCV$numberSamples,1)))%>%
       arrange(year, population)
     
     qaly_incre_scen[[i]][[t]] <- tt$cumu_QALY[[i]]%>%
       filter(testing == t)%>%
       select(year, testing, population, scenario, indicator, Casca, best,
              paste0("set", seq(1,HCV$numberSamples,1)))%>%
       arrange(year, population)
     
   }
 }

 
Iavert <- list()

prIavert <- list()

reIavert <- list()
 
cost_incre <- list()
 
cost_stage_incre <- list()
 
qaly_incre <- list()
 
 for(i in dt_name){ 
   for(t in t_lab){ 
     
     Iavert[[i]][[t]] <- cbind(year = Iavert_bench[[i]]$year,
                               population = Iavert_bench[[i]]$population,
                               scenario = Iavert_bench[[i]]$scenario,
                               indicator = Iavert_bench[[i]]$indicator,
                               Casca = Iavert_bench[[i]]$Casca,
                               testing = t,
                               as.data.frame(
                                 (Iavert_bench[[i]][, -c(1:6)] - 
                                    Iavert_scen[[i]][[t]][,-c(1:6)])))%>%
       popRange(., Casca = "y")                                                 # Casca != NULL >> arrange all columns in character format  
     
     prIavert[[i]][[t]] <- cbind(year = prIavert_bench[[i]]$year,
                               population = prIavert_bench[[i]]$population,
                               scenario = prIavert_bench[[i]]$scenario,
                               indicator = prIavert_bench[[i]]$indicator,
                               Casca = prIavert_bench[[i]]$Casca,
                               testing = t,
                               as.data.frame(
                                 (prIavert_bench[[i]][, -c(1:6)] - 
                                    prIavert_scen[[i]][[t]][,-c(1:6)])))%>%
       popRange(., Casca = "y")  
     
     reIavert[[i]][[t]] <- cbind(year = reIavert_bench[[i]]$year,
                                 population = reIavert_bench[[i]]$population,
                                 scenario = reIavert_bench[[i]]$scenario,
                                 indicator = reIavert_bench[[i]]$indicator,
                                 Casca = reIavert_bench[[i]]$Casca,
                                 testing = t,
                                 as.data.frame(
                                   (reIavert_bench[[i]][, -c(1:6)] - 
                                      reIavert_scen[[i]][[t]][,-c(1:6)])))%>%
       popRange(., Casca = "y")
     
     
     
     cost_incre[[i]][[t]] <- cbind(year = cost_incre_bench[[i]]$year,
                                   population = cost_incre_bench[[i]]$population,
                                   scenario = cost_incre_bench[[i]]$scenario,
                                   indicator = cost_incre_bench[[i]]$indicator,
                                   Casca = cost_incre_bench[[i]]$Casca,
                                   testing = t,
                                   as.data.frame(
                                     (cost_incre_scen[[i]][[t]][,-c(1:6)] - 
                                        cost_incre_bench[[i]][, -c(1:6)])))%>%
       popRange(., Casca = "y")
     
     
     cost_stage_incre[[i]][[t]] <- cbind(year = cost_stage_incre_bench[[i]]$year,
                                         population = cost_stage_incre_bench[[i]]$population,
                                         scenario = cost_stage_incre_bench[[i]]$scenario,
                                         indicator = cost_stage_incre_bench[[i]]$indicator,
                                         Casca = cost_stage_incre_bench[[i]]$Casca,
                                         testing = t,
                                         as.data.frame(
                                           (cost_stage_incre_scen[[i]][[t]][,-c(1:6)] - 
                                              cost_stage_incre_bench[[i]][, -c(1:6)])))%>%
       popRange(., Casca = "y")
     
     
     qaly_incre[[i]][[t]] <- cbind(year = qaly_incre_bench[[i]]$year,
                                   population = qaly_incre_bench[[i]]$population,
                                   scenario = qaly_incre_bench[[i]]$scenario,
                                   indicator = qaly_incre_bench[[i]]$indicator,
                                   Casca = qaly_incre_bench[[i]]$Casca,
                                   testing = t,
                                   as.data.frame(
                                     (qaly_incre_scen[[i]][[t]][,-c(1:6)] - 
                                        qaly_incre_bench[[i]][, -c(1:6)])))%>%
       popRange(., Casca = "y")
     
     
     
   }
   
 }

CEA <- list()

for(i in dt_name){ 
  for(t in t_lab){ 
    
    cost_incre[[i]][[t]] <- cost_incre[[i]][[t]]%>%arrange(year, population)
    
    qaly_incre[[i]][[t]] <- qaly_incre[[i]][[t]]%>%arrange(year, population)
    
    CEA[[i]][[t]] <- cbind(year = qaly_incre[[i]][[t]]$year,
                           population = qaly_incre[[i]][[t]]$population,
                           scenario = qaly_incre[[i]][[t]]$scenario,
                           indicator = "ICER",
                           Casca = qaly_incre[[i]][[t]]$Casca,
                           testing = t,
                           as.data.frame(
                             (cost_incre[[i]][[t]][,-c(1:6)]/ 
                                qaly_incre[[i]][[t]][, -c(1:6)])))%>%
      select(year, population, scenario, indicator, Casca, testing, best, 
             paste0("set", seq(1, HCV$numberSamples,1)))%>%
      popRange(., Casca = "y")
    

  }
} 


for(i in dt_name){ 
  for(t in t_lab){ 
    CEA[[i]][[t]] <- CEA[[i]][[t]]%>%
      mutate(cost_p_med  = ifelse(cost_incre[[i]][[t]]$Med>0,1,0), 
             qaly_p_med  = ifelse(qaly_incre[[i]][[t]]$Med>0,1,0),
             cost_p_q5  = ifelse(cost_incre[[i]][[t]]$q5>0,1,0), 
             qaly_p_q5  = ifelse(qaly_incre[[i]][[t]]$q5>0,1,0),
             cost_p_q95  = ifelse(cost_incre[[i]][[t]]$q95>0,1,0), 
             qaly_p_q95  = ifelse(qaly_incre[[i]][[t]]$q95>0,1,0)
             )%>%
      mutate(dominant_Med = case_when((Med <0 & cost_p_med == 0)~"dominant",
                                       (Med <0 & cost_p_med == 1)~"dominanted",
                                       .default = scales::dollar(Med, prefix = " $")),
             
             
             dominant_q5 = case_when((q5 <0 & cost_p_q5 == 0 )~ "dominant",
                                      (q5 <0 & cost_p_q5 == 1) ~ "dominanted",
                                      .default = scales::dollar(q5, prefix = " $")),
             
             dominant_q95 = case_when((q95 <0 & cost_p_q95 == 0) ~ "dominant",
                                     (q95 <0 & cost_p_q95 == 1) ~ "dominanted",
                                     .default = scales::dollar(q95, prefix = " $"))
             
             )
  }
}

################ CEA table color ###############################################
CEA_table <- list()

for(i in dt_name){ 
  
  CEA_table[[i]] <-do.call("rbind", CEA[[i]])%>%
    select(year, population, scenario,indicator, testing, dominant_Med, 
           dominant_q5, dominant_q95, Med)
  
  }

CEA_table <- do.call("rbind", CEA_table)
CEA_table <- CEA_table%>%mutate(testing = factor(testing, 
                          levels = c(unique(CEA_table$testing)), 
                          label = c("Point-of-care antibody testing",
                                    "Clinic-based reflex RNA testing",
                                    "Single visit point-of-care RNA testing")),
         scenario = factor(scenario, 
                           levels = c("PrEPsame", "PrEP", "HIVD", "PrEPnHIVD"),
                           label = c("PrEP coverage and HIV treatment coverage remains unchanged",
                                     "PrEP coverage increases to 20% by 2030",
                                     "HIV treatment coverages increase to 95% by 2030",
                                     "PrEP coverage increases to 20% and HIV treatment increased to 95% by 2030")),
         
         population = factor(population, 
                             levels = c("All", "HIV-PrEP", "HIV+d"),
                             label = c("Entire population of MSM", 
                                       "HIV-negative MSM on PrEP", 
                                       "MSM who are HIV diagnosed (and on treatment)"))
         )
row.names(CEA_table) <- NULL

CEA_table <- CEA_table%>%as.data.frame()%>%
  mutate(ICER = paste0(dominant_Med, "\n", "(", dominant_q5, "-", dominant_q95, ")"))%>%
  gather(indicator, vv, -c(population, scenario, testing,year, dominant_Med, 
                        dominant_q5, dominant_q95))%>%
  select(-c(year,dominant_Med, 
            dominant_q5, dominant_q95))%>%
  filter(indicator %in% c("ICER", "Med"))

# row bind to a dataframe 
# combine stage cost as separate indicators 
 dt_incre <- list()
 for(i in dt_name){
 
   dt_incre[[i]] <- rbind(do.call("rbind", Iavert[[i]]),
                          do.call("rbind", prIavert[[i]]),
                          do.call("rbind", reIavert[[i]]),
                          do.call("rbind", cost_incre[[i]]),
                          do.call("rbind", cost_stage_incre[[i]]),
                          do.call("rbind", qaly_incre[[i]])
                          )%>%as_tibble()%>%
     mutate(Casca = ifelse(Casca == "total", "", Casca),
            indicator = ifelse(Casca!= "", 
                               paste0(indicator, "_", Casca, sep = ""),
                               indicator))%>%
     select(year, population, scenario, indicator, testing, 
            Med, q5, q95)
   
 }
 

 dt_incre <- do.call("rbind", dt_incre)
 save( dt_incre,
      file = file.path(outputdt, "CostBenfit_incre.rda"))

 dt_incre <-  dt_incre%>%mutate(testing = factor(testing, 
                         levels = c(unique(dt_incre$testing)), 
                         label = c("Point-of-care antibody testing",
                                   "Clinic-based reflex RNA testing",
                                   "Single visit point-of-care RNA testing")),
        scenario = factor(scenario, 
                          levels = c("PrEPsame", "PrEP", "HIVD", "PrEPnHIVD"),
                          label = c("PrEP coverage and HIV treatment coverage remains unchanged",
                                    "PrEP coverage increases to 20% by 2030",
                                    "HIV treatment coverages increase to 95% by 2030",
                                    "PrEP coverage increases to 20% and HIV treatment increased to 95% by 2030")),
        
        population = factor(population, 
                            levels = c("All", "HIV-PrEP", "HIV+d"),
                            label = c("Entire population of MSM", 
                                      "HIV-negative MSM on PrEP", 
                                      "MSM who are HIV diagnosed (and on treatment)")),
     
        indicator = factor(indicator, 
                           levels = c( "cumu_newinf", "cumu_newpriinf", "cumu_newreinf","cumu_lifecost",
                                       "cumu_lifecost_stage_diag",
                                       "cumu_lifecost_stage_treat",
                                       "cumu_lifecost_stage_management",
                                       "cumu_QALY" ),
                           label = c("HCV infections averted (2022-2030)",
                                     "HCV primary infections averted (2022-2030)",
                                     "HCV reinfections averted (2022-2030)",
                                     "Incremental lifetime cost (discounted, thousands US$)",
                                     "Incremental lifetime cost of HCV diagnosis (discounted, thousands US$)",
                                     "Incremental lifetime cost of HCV treatment (discounted, thousands US$)",
                                     "Incremental lifetime cost of HCV management (discounted, thousands US$)",
                                     "QALY gained (discounted)"
                                     
                           ))
 ) 
 
 dt_incre <- dt_incre %>% 
   select(scenario, indicator, population, testing, 
          Med, q5, q95)%>%
  filter(!indicator %in% c("Incremental lifetime cost of HCV diagnosis (discounted, thousands US$)",
                           "Incremental lifetime cost of HCV treatment (discounted, thousands US$)",
                           "Incremental lifetime cost of HCV management (discounted, thousands US$)"))%>%
   mutate(Med = case_when(indicator %in% 
                            c("Incremental lifetime cost (discounted, thousands US$)")~ round(Med/1000, digits = 1),
                           indicator %in% c("QALY gained (lifetime)")~ round(Med, digits = 0),
                           indicator %in% c("HCV infections averted (2022-2030)",
                                            "HCV primary infections averted (2022-2030)",
                                            "HCV reinfections averted (2022-2030)")~ round(Med, digits =  1),
                           .default = Med
   ),
   q5= case_when(indicator %in% c("Incremental lifetime cost (discounted, thousands US$)")~ round(q5/1000, digits = 1),
                 indicator %in% c("QALY gained (lifetime)") ~ round(q5, digits = 0),
                 indicator %in% c("HCV infections averted (2022-2030)")~ round(q5, digits =  1),
                 .default = q5
   ),
   q95 = case_when(indicator %in% c("Incremental lifetime cost (discounted, thousands US$)")~ round(q95/1000, digits = 1),
                   indicator %in% c("QALY gained (lifetime)")~ round(q95,digits = 0),
                   indicator %in% c("HCV infections averted (2022-2030)",
                                    "HCV primary infections averted (2022-2030)",
                                    "HCV reinfections averted (2022-2030)")~ round(q95, digits =  1),
                   .default = q95))%>%
   mutate(Med = formatC(Med, format = "fg", big.mark = ","),
          q5 = formatC(q5,format = "fg", big.mark = ","),
          q95 = formatC(q95,format = "fg", big.mark = ","))%>%
   mutate(vv = paste0(Med, "\n", "(", q5, "-", q95, ")"))%>%
   select(-c(Med, q5, q95))%>%
   # reorder columns to present in the right order
   arrange(indicator)

 # bind CEA table 
dt_incre <- rbind(dt_incre, CEA_table)

# long form 
dt_incre <- dt_incre%>%tidyr::pivot_wider(names_from = c("indicator","population"), 
                      values_from = c("vv"), 
                      names_glue = "{indicator}X{population}")
   # grouping rows
dt_incre%>%gt(groupname_col = "scenario",
      rowname_col = "testing")%>%
   tab_spanner_delim(delim = "X")%>%
   tab_style(
     style = list(
       cell_fill(color = "#DEEBF7"),
       cell_text(weight = "bold")),
     locations = cells_body(
       columns = c(`ICERXHIV-negative MSM on PrEP`,
       ),
       rows = c(as.numeric(`MedXHIV-negative MSM on PrEP`) <= 1*TW_GDP2022))
   )%>%
   tab_style(
     style = list(
       cell_fill(color = "#4292C6"),
       cell_text(weight = "bold")
     ),
     locations = cells_body(
       columns = c(`ICERXHIV-negative MSM on PrEP`),
       rows = c(as.numeric(`MedXHIV-negative MSM on PrEP`) > 1*TW_GDP2022)))%>%
   tab_style(
     style = list(
       cell_fill(color = "#DEEBF7"),
       cell_text(weight = "bold")),
     locations = cells_body(
       columns = c(`ICERXEntire population of MSM`,
       ),
       rows = c(as.numeric(`MedXEntire population of MSM`) <= 1*TW_GDP2022))
   )%>%
   tab_style(
     style = list(
       cell_fill(color = "#4292C6"),
       cell_text(weight = "bold")
     ),
     locations = cells_body(
       columns = c(`ICERXEntire population of MSM`),
       rows = c(as.numeric(`MedXEntire population of MSM`) > 1*TW_GDP2022
       )))%>%
   tab_style(
     style = list(
       cell_fill(color = "#DEEBF7"),
       cell_text(weight = "bold")),
     locations = cells_body(
       columns = c(`ICERXMSM who are HIV diagnosed (and on treatment)`),
       rows = c(`MedXMSM who are HIV diagnosed (and on treatment)`) <= 1*TW_GDP2022))%>%
   tab_style(
     style = list(
       cell_fill(color = "#4292C6"),
       cell_text(weight = "bold")
     ),
     locations = cells_body(
       columns = c(`ICERXMSM who are HIV diagnosed (and on treatment)`),
       rows = c(as.numeric(`MedXMSM who are HIV diagnosed (and on treatment)`) > 1*TW_GDP2022)
       )
     )%>%
   cols_hide(columns = c(`MedXMSM who are HIV diagnosed (and on treatment)`, 
                         `MedXEntire population of MSM`, 
                         `MedXHIV-negative MSM on PrEP`))%>%
  gtsave(., file = file.path(outputfig, paste0("mainRes_incre_0913.docx"))) 


#### 2023/09/13 #### 
# reorganize table 
##### ICER table #####
# main table in main text 
x_total <- x%>%select(scenario, testing, 
           `Lifetime cost (discounted, millions US$)XEntire population of MSM`,
           `Lifetime QALY (discounted)XEntire population of MSM`)%>%
  mutate(testing = factor(testing, level = c(unique(testing)),
                          labels = c("Status Quo", "Point-of-care antibody testing",
                                     "Clinic-based reflex RNA testing",
                                     "Single visit point-of-care RNA testing")))


x_incre <- dt_incre%>%
  select(scenario, testing,
         `Incremental lifetime cost (discounted, thousands US$)XEntire population of MSM`,
         `QALY gained (discounted)XEntire population of MSM`,
         `MedXEntire population of MSM`,
         `ICERXEntire population of MSM`)

incre_sq <- x_incre%>%filter(testing =="Point-of-care antibody testing")%>%
  mutate(testing = "Status Quo")

incre_sq[, 3:ncol(incre_sq)] <- "NA"

incre_wsq <- rbind(x_incre, incre_sq)
  
xxxt <- merge(x_total,incre_wsq, by = c("testing", "scenario"), all = TRUE)%>%
  mutate(testing = factor(testing, levels = c("Status Quo", "Point-of-care antibody testing",
                                              "Clinic-based reflex RNA testing",
                                              "Single visit point-of-care RNA testing")))%>%
  arrange(scenario, testing)

xxxt%>%gt(groupname_col = "scenario",
          rowname_col = "testing")%>%
  tab_spanner_delim(delim = "X")%>%
  tab_style(
    style = list(
      cell_fill(color = "#DEEBF7"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = c(`ICERXEntire population of MSM`,
      ),
      rows = c(as.numeric(`MedXEntire population of MSM`) <= 1*TW_GDP2022))
  )%>%
  tab_style(
    style = list(
      cell_fill(color = "#4292C6"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = c(`ICERXEntire population of MSM`),
      rows = c(as.numeric(`MedXEntire population of MSM`) > 1*TW_GDP2022
      )))%>%
  cols_hide(columns = c(`MedXEntire population of MSM`))%>%
  gtsave(., file = file.path(outputfig, paste0("mainRes_Tab_0913.docx")))

# subpops in supplementary 
# main table in main text 
# MSM on HIV treatment 

x_total <- x%>%select(scenario, testing, 
                      `Lifetime cost (discounted, millions US$)XMSM who are HIV diagnosed (and on treatment)`,
                      `Lifetime QALY (discounted)XMSM who are HIV diagnosed (and on treatment)`)%>%
  mutate(testing = factor(testing, level = c(unique(testing)),
                          labels = c("Status Quo", "Point-of-care antibody testing",
                                     "Clinic-based reflex RNA testing",
                                     "Single visit point-of-care RNA testing")))


x_incre <- dt_incre%>%
  select(scenario, testing,
         `Incremental lifetime cost (discounted, thousands US$)XMSM who are HIV diagnosed (and on treatment)`,
         `QALY gained (discounted)XMSM who are HIV diagnosed (and on treatment)`,
         `MedXMSM who are HIV diagnosed (and on treatment)`,
         `ICERXMSM who are HIV diagnosed (and on treatment)`)

incre_sq <- x_incre%>%filter(testing =="Point-of-care antibody testing")%>%
  mutate(testing = "Status Quo")

incre_sq[, 3:ncol(incre_sq)] <- "NA"

incre_wsq <- rbind(x_incre, incre_sq)

xxxt <- merge(x_total,incre_wsq, by = c("testing", "scenario"), all = TRUE)%>%
  mutate(testing = factor(testing, levels = c("Status Quo", "Point-of-care antibody testing",
                                              "Clinic-based reflex RNA testing",
                                              "Single visit point-of-care RNA testing")))%>%
  arrange(scenario, testing)

xxxt%>%gt(groupname_col = "scenario",
          rowname_col = "testing")%>%
  tab_spanner_delim(delim = "X")%>%
  tab_style(
    style = list(
      cell_fill(color = "#DEEBF7"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = c(`ICERXMSM who are HIV diagnosed (and on treatment)`,
      ),
      rows = c(as.numeric(`MedXMSM who are HIV diagnosed (and on treatment)`) <= 1*TW_GDP2022))
  )%>%
  tab_style(
    style = list(
      cell_fill(color = "#4292C6"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = c(`ICERXMSM who are HIV diagnosed (and on treatment)`),
      rows = c(as.numeric(`MedXMSM who are HIV diagnosed (and on treatment)`) > 1*TW_GDP2022
      )))%>%
  cols_hide(columns = c(`MedXMSM who are HIV diagnosed (and on treatment)`))%>%
  gtsave(., file = file.path(outputfig, paste0("suppleRes_Tab_HIVMSM_0913.docx")))

# MSM on PrEP
x_total <- x%>%select(scenario, testing, 
                      `Lifetime cost (discounted, millions US$)XHIV-negative MSM on PrEP`,
                      `Lifetime QALY (discounted)XHIV-negative MSM on PrEP`)%>%
  mutate(testing = factor(testing, level = c(unique(testing)),
                          labels = c("Status Quo", "Point-of-care antibody testing",
                                     "Clinic-based reflex RNA testing",
                                     "Single visit point-of-care RNA testing")))


x_incre <- dt_incre%>%
  select(scenario, testing,
         `Incremental lifetime cost (discounted, thousands US$)XHIV-negative MSM on PrEP`,
         `QALY gained (discounted)XHIV-negative MSM on PrEP`,
         `MedXHIV-negative MSM on PrEP`,
         `ICERXHIV-negative MSM on PrEP`)

incre_sq <- x_incre%>%filter(testing =="Point-of-care antibody testing")%>%
  mutate(testing = "Status Quo")

incre_sq[, 3:ncol(incre_sq)] <- "NA"

incre_wsq <- rbind(x_incre, incre_sq)

xxxt <- merge(x_total,incre_wsq, by = c("testing", "scenario"), all = TRUE)%>%
  mutate(testing = factor(testing, levels = c("Status Quo", "Point-of-care antibody testing",
                                              "Clinic-based reflex RNA testing",
                                              "Single visit point-of-care RNA testing")))%>%
  arrange(scenario, testing)

xxxt%>%gt(groupname_col = "scenario",
          rowname_col = "testing")%>%
  tab_spanner_delim(delim = "X")%>%
  tab_style(
    style = list(
      cell_fill(color = "#DEEBF7"),
      cell_text(weight = "bold")),
    locations = cells_body(
      columns = c(`ICERXHIV-negative MSM on PrEP`,
      ),
      rows = c(as.numeric(`MedXHIV-negative MSM on PrEP`) <= 1*TW_GDP2022))
  )%>%
  tab_style(
    style = list(
      cell_fill(color = "#4292C6"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = c(`ICERXHIV-negative MSM on PrEP`),
      rows = c(as.numeric(`MedXHIV-negative MSM on PrEP`) > 1*TW_GDP2022
      )))%>%
  cols_hide(columns = c(`MedXHIV-negative MSM on PrEP`))%>%
  gtsave(., file = file.path(outputfig, paste0("suppleRes_Tab_PrEPMSM_0913.docx")))




#### infection averted 

##### Infection and infection averted table 
# total infection 
x_total <- x%>%select(scenario, testing, 
                      `New HCV infections (2022-2030)XEntire population of MSM`,
                      `New HCV primary infections (2022-2030)XEntire population of MSM`,
                      `New HCV reinfections (2022-2030)XEntire population of MSM`
                      )%>%
  mutate(testing = factor(testing, level = c(unique(testing)),
                          labels = c("Status Quo", "Point-of-care antibody testing",
                                     "Clinic-based reflex RNA testing",
                                     "Single visit point-of-care RNA testing")))

Inf_averted <- dt_incre%>%
  select(scenario, testing,
         `HCV infections averted (2022-2030)XEntire population of MSM`,
         `HCV primary infections averted (2022-2030)XEntire population of MSM`,
         `HCV reinfections averted (2022-2030)XEntire population of MSM`
         )

incre_sq <- Inf_averted %>%filter(testing =="Point-of-care antibody testing")%>%
  mutate(testing = "Status Quo")

incre_sq[, 3:ncol(incre_sq)] <- "NA"

incre_wsq <- rbind(Inf_averted , incre_sq)

xxxt <- merge(x_total,incre_wsq, by = c("testing", "scenario"), all = TRUE)%>%
  mutate(testing = factor(testing, levels = c("Status Quo", "Point-of-care antibody testing",
                                              "Clinic-based reflex RNA testing",
                                              "Single visit point-of-care RNA testing")))%>%
  arrange(scenario, testing)

xxxt%>%gt(groupname_col = "scenario",
          rowname_col = "testing")%>%
  tab_spanner_delim(delim = "X")%>%
  gtsave(., file = file.path(outputfig, paste0("suppleRes_Tab_infAvrALL_0913.docx")))


# MSM on PrEP 
x_total <- x%>%select(scenario, testing, 
                      `New HCV infections (2022-2030)XHIV-negative MSM on PrEP`,
                      `New HCV primary infections (2022-2030)XHIV-negative MSM on PrEP`,
                      `New HCV reinfections (2022-2030)XHIV-negative MSM on PrEP`,
)%>%
  mutate(testing = factor(testing, level = c(unique(testing)),
                          labels = c("Status Quo", "Point-of-care antibody testing",
                                     "Clinic-based reflex RNA testing",
                                     "Single visit point-of-care RNA testing")))

Inf_averted <- dt_incre%>%
  select(scenario, testing,
         `HCV infections averted (2022-2030)XHIV-negative MSM on PrEP`,
         `HCV primary infections averted (2022-2030)XHIV-negative MSM on PrEP`,
         `HCV reinfections averted (2022-2030)XHIV-negative MSM on PrEP`
  )

incre_sq <- Inf_averted %>%filter(testing =="Point-of-care antibody testing")%>%
  mutate(testing = "Status Quo")

incre_sq[, 3:ncol(incre_sq)] <- "NA"

incre_wsq <- rbind(Inf_averted , incre_sq)

xxxt <- merge(x_total,incre_wsq, by = c("testing", "scenario"), all = TRUE)%>%
  mutate(testing = factor(testing, levels = c("Status Quo", "Point-of-care antibody testing",
                                              "Clinic-based reflex RNA testing",
                                              "Single visit point-of-care RNA testing")))%>%
  arrange(scenario, testing)

xxxt%>%gt(groupname_col = "scenario",
          rowname_col = "testing")%>%
  tab_spanner_delim(delim = "X")%>%
  gtsave(., file = file.path(outputfig, paste0("suppleRes_Tab_infAvrMSMPrEP_0913.docx")))

# MSM on HIV treatment
x_total <- x%>%select(scenario, testing, 
                      `New HCV infections (2022-2030)XMSM who are HIV diagnosed (and on treatment)`,
                      `New HCV primary infections (2022-2030)XMSM who are HIV diagnosed (and on treatment)`,
                      `New HCV reinfections (2022-2030)XMSM who are HIV diagnosed (and on treatment)`
)%>%
  mutate(testing = factor(testing, level = c(unique(testing)),
                          labels = c("Status Quo", "Point-of-care antibody testing",
                                     "Clinic-based reflex RNA testing",
                                     "Single visit point-of-care RNA testing")))

Inf_averted <- dt_incre%>%
  select(scenario, testing,
         `HCV infections averted (2022-2030)XMSM who are HIV diagnosed (and on treatment)`,
         `HCV primary infections averted (2022-2030)XMSM who are HIV diagnosed (and on treatment)`,
         `HCV reinfections averted (2022-2030)XMSM who are HIV diagnosed (and on treatment)`
  )

incre_sq <- Inf_averted %>%filter(testing =="Point-of-care antibody testing")%>%
  mutate(testing = "Status Quo")

incre_sq[, 3:ncol(incre_sq)] <- "NA"

incre_wsq <- rbind(Inf_averted , incre_sq)

xxxt <- merge(x_total,incre_wsq, by = c("testing", "scenario"), all = TRUE)%>%
  mutate(testing = factor(testing, levels = c("Status Quo", "Point-of-care antibody testing",
                                              "Clinic-based reflex RNA testing",
                                              "Single visit point-of-care RNA testing")))%>%
  arrange(scenario, testing)

xxxt%>%gt(groupname_col = "scenario",
          rowname_col = "testing")%>%
  tab_spanner_delim(delim = "X")%>%
  gtsave(., file = file.path(outputfig, paste0("suppleRes_Tab_infAvrMSMHIV_0913.docx")))


##### Cost categories ##### 
cost_tab <- list()
for(i in seq_along(tt$cumu_lifecost_stage)){ 
  cost_tab[[i]] <- tt$cumu_lifecost_stage[[i]]%>%
    mutate(testing = factor(testing, 
                          levels = c("Status Quo", "POC_antibody", "Reflex_RNA",
                                     "POC_RNA"), 
                          label = c("Status Quo",
                                    "Point-of-care anibody testing",
                                    "Clinic-based reflex RNA testing",
                                    "Single visit point-of-care RNA testing"
                          )),
         population = factor(population, 
                             levels = c("All", "HIV-PrEP", "HIV+d"),
                             label = c("Entire population of MSM", 
                                       "HIV-negative MSM on PrEP", 
                                       "MSM who are HIV diagnosed (and on treatment)")),
         Casca = factor(Casca, 
                        levels = c("diag", "management", "treat"),
                        label = c("HCV diagnosis", 
                                  "HCV management and follow up", 
                                  "HCV treatment"))
  )
}

x_cost <- cost_tab[[1]]%>%select(year, population, testing, scenario, 
                       indicator, Casca, Med, q5, q95)%>%
  mutate()


#### haven't test 2023/06/07================================================#### 


####cost of cascade stage plot  ####

plot_costcasca <- list()

scenario_lab <- c("PrEP coverage and HIV treatment coverage remains unchanged",
                  "PrEP coverage increases to 20% by 2030",
                  "HIV treatment coverages increase to 95% by 2030",
                  "PrEP coverage increases to 20% and HIV treatment increased to 95% by 2030")
plot_tag <- c("A", "B", "C", "D")

lifecost_stage <- list()
maxvalue <- list()

for(i in seq_along(tt$cumu_lifecost_stage)){ 
  
  lifecost_stage[[i]] <-tt$cumu_lifecost_stage[[i]]%>%
    mutate(testing = factor(testing, 
                            levels = c("Status Quo", "POC_antibody", "Reflex_RNA",
                                       "POC_RNA"), 
                            label = c("Current practice of HCV testing",
                                      "Point-of-care anibody testing",
                                      "Reflex RNA testing",
                                      "Point-of-care RNA testing"
                                      )),
           population = factor(population, 
                               levels = c("All", "HIV-PrEP", "HIV+d"),
                               label = c("Entire population of MSM", 
                                         "HIV-negative MSM on PrEP", 
                                         "MSM who are HIV diagnosed (and on treatment)")),
           Casca = factor(Casca, 
                            levels = c("diag", "management", "treat"),
                            label = c("HCV diagnosis", 
                                      "HCV management", 
                                      "HCV treatment"))
           ) 
  
  
  maxvalue[[i]] <- lifecost_stage[[i]]%>%group_by(population)%>%
    summarize(value = case_when(max(best/10000000) < 1 ~ round(max(best/10000000), digits = 1),
                                max(best/10000000) >=1  ~ ceiling(max(best/10000000)
                                )))


  plot_costcasca[[i]] <- lifecost_stage[[i]]%>%
    filter(year == 2090 & population == "Entire population of MSM")%>%
    ggplot(., aes(x= testing , y = best/(Currency_converter*1000000), 
                  fill = Casca)) + 
    geom_bar(stat = "identity", width = 0.6) + 
    scale_fill_viridis(discrete = TRUE, option = "mako") +
     
    labs(fill = "HCV Cascade", tag = plot_tag[i], x = "", y = "") + 
    ggtitle(scenario_lab[i])  +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(panel.background = element_rect(colour = "white"),
          plot.background = element_rect(colour = "white"),
          panel.border     = element_rect(fill = NA, colour = "black", 
                                          size = NA),
          plot.title = element_text(face = "bold",
                                    size = 14, hjust = 0.5),
          text = element_text(),
          axis.title = element_text(face = "bold",size = 14),
          axis.title.y = element_text(angle=90,vjust =1),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.x = element_text(angle=45, vjust = 1,hjust = 1,
                                     face = "bold",size = 12, colour = "black"), 
          axis.text.y = element_text(
            face = "bold",size = 14,colour = "black"),
          strip.text.x = element_text(size=14, color="black",
                                      face="bold"),
          legend.text = element_text(size = 14, face = "bold"),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(face="bold", size= 14),
          plot.margin = unit(c(10,5,5,5),"mm")
          
    ) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0.001, .1)), 
                       limits = 
                         c(0,as.numeric(maxvalue[[i]][1, 2])))
  
}


gplotlist_costcasca <- ggarrange(plotlist = plot_costcasca, 
                                 common.legend = TRUE, ncol  = 2, nrow = 2) + 
  theme(legend.position = "bottom")

figs11 <- annotate_figure(gplotlist_costcasca,
                          bottom = text_grob("Testing strategies",,
                                             hjust = 1, vjust = -1, x = 0.55,  
                                             size = 14, face ="bold"),
                          left = text_grob("Cost (Millions in USD$)",  rot = 90, 
                                           size = 14, face = "bold")
)



ggsave(path = outputfig, file="cost_cas_all.jpeg", 
       figs11, width = 16, 
       height = 13, dpi = 300) 

# subgroup 
plot_costcasca_prep <- list()
plot_costcasca_hivd <- list()
for(i in seq_along(lifecost_stage)){ 
  plot_costcasca_prep[[i]] <- lifecost_stage[[i]]%>%
    filter(year == 2090 & population == "HIV-negative MSM on PrEP")%>%
    ggplot(., aes(x= testing , y = best/(Currency_converter*1000000), 
                  fill = Casca)) + 
    geom_bar(stat = "identity", width = 0.6) + 
    scale_fill_viridis(discrete = TRUE, option = "mako") +
    
    labs(fill = "HCV Cascade", tag = plot_tag[i], x = "", y = "") + 
    ggtitle(scenario_lab[i])  +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(panel.background = element_rect(colour = "white"),
          plot.background = element_rect(colour = "white"),
          panel.border     = element_rect(fill = NA, colour = "black", 
                                          size = NA),
          plot.title = element_text(face = "bold",
                                    size = 12, hjust = 0.5),
          text = element_text(),
          axis.title = element_text(face = "bold",size = 14),
          axis.title.y = element_text(angle=90,vjust =1),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.x = element_text(angle=45, vjust = 1,hjust = 1,
                                     face = "bold",size = 12, colour = "black"), 
          axis.text.y = element_text(
            face = "bold",size = 14,colour = "black"),
          strip.text.x = element_text(size=14, color="black",
                                      face="bold"),
          legend.text = element_text(size = 14, face = "bold"),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(face="bold", size= 14),
          plot.margin = unit(c(10,5,5,5),"mm")
          
    ) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0.001, .1)), 
                       limits = 
                         c(0,as.numeric(maxvalue[[i]][2, 2])))
  
  plot_costcasca_hivd[[i]] <- lifecost_stage[[i]]%>%
    filter(year == 2090 & population == "MSM who are HIV diagnosed (and on treatment)")%>%
    ggplot(., aes(x= testing , y = best/(Currency_converter*1000000), 
                  fill = Casca)) + 
    geom_bar(stat = "identity", width = 0.6) + 
    scale_fill_viridis(discrete = TRUE, option = "mako") +
    
    labs(fill = "HCV Cascade", tag = plot_tag[i], x = "", y = "") + 
    ggtitle(scenario_lab[i])  +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(panel.background = element_rect(colour = "white"),
          plot.background = element_rect(colour = "white"),
          panel.border     = element_rect(fill = NA, colour = "black", 
                                          size = NA),
          plot.title = element_text(face = "bold",
                                    size = 12, hjust = 0.5),
          text = element_text(),
          axis.title = element_text(face = "bold",size = 14),
          axis.title.y = element_text(angle=90,vjust =1),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.x = element_text(angle=45, vjust = 1,hjust = 1,
                                     face = "bold",size = 12, colour = "black"), 
          axis.text.y = element_text(
            face = "bold",size = 14,colour = "black"),
          strip.text.x = element_text(size=14, color="black",
                                      face="bold"),
          legend.text = element_text(size = 14, face = "bold"),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(face="bold", size= 14),
          plot.margin = unit(c(10,5,5,5),"mm")
          
    ) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0.001, .1)), 
                       limits = 
                         c(0,as.numeric(maxvalue[[i]][3, 2])))
  
  
  }

gplotlist_costcasca_prep <- ggarrange(
  plotlist = list(plot_costcasca_prep[[1]], 
                  plot_costcasca_prep[[2]],
                  plot_costcasca_prep[[3]],
                  plot_costcasca_prep[[4]]), nrow = 2, ncol = 2,
  common.legend = TRUE) 

gplotlist_costcasca_hivd <- ggarrange(plotlist = list(
  plot_costcasca_hivd[[1]], 
  plot_costcasca_hivd[[2]],
  plot_costcasca_hivd[[3]],
  plot_costcasca_hivd[[4]]), nrow = 2, ncol = 2,
  common.legend = TRUE) + theme(legend.position = "")



figs12 <- annotate_figure(gplotlist_costcasca_prep,
                          bottom = text_grob("Testing strategies",,
                                             hjust = 1, vjust = -1, x = 0.55,  
                                             size = 14, face ="bold"),
                          left = text_grob("Cost (Millions in USD$)",  rot = 90, 
                                           size = 14, face = "bold")
)

figs13 <- annotate_figure(gplotlist_costcasca_hivd,
                          bottom = text_grob("Testing strategies",,
                                             hjust = 1, vjust = -1, x = 0.55,  
                                             size = 14, face ="bold"),
                          left = text_grob("Cost (Millions in USD$)",  rot = 90, 
                                           size = 14, face = "bold")
)

ggsave(path = outputfig, file="cost_cas_prep.jpeg", 
       figs12, width = 14, 
       height = 12, dpi = 300) 

ggsave(path = outputfig, file="cost_cas_hivd.jpeg", 
       figs13, width = 14, 
       height = 12, dpi = 300) 

#### numbers for results #### 
xt_num <- do.call("rbind",lifecost_stage)%>%
  select(year, population, testing, scenario, Casca, Med, q5, q95)
xt_numSum <- xt_num%>%group_by(year,population, testing,scenario)%>%
  summarize(Med = sum(Med),
            q5 = sum(q5),
            q95 = sum(q95))%>%mutate(Casca = "Total")%>%
  select(year,population, testing,scenario, Casca, Med, q5, q95)%>%
  arrange(population,scenario, testing)

ttt <- xt_num%>%
  filter(Casca == "HCV diagnosis" & population == "Entire population of MSM")%>%
  arrange(Med)
tttt$`HCV treatment`
tttt <- rbind(xt_num, xt_numSum)%>%
  arrange(year, population, testing, scenario)%>%
  pivot_wider(names_from = Casca, values_from = Med)%>%
  select(!c(`HCV treatment`, `HCV management and follow up`))

tttt_diag <- tttt%>%select(!Total)
tttt_diag <- tttt_diag%>%na.omit()
tttt_total <- tttt%>%select(Total)%>%na.omit()

tttt_p <- cbind(tttt_diag, tttt_total)%>%mutate(per = `HCV diagnosis`/Total)%>%
  arrange(per)%>%filter(population == "Entire population of MSM")

#### plot_costcasa for entire MSM only #### 
test3[1, 7]
test3 <- do.call("rbind",lifecost_stage)%>%
  mutate(scenario = factor(scenario, 
                           levels = c("PrEPsame", "PrEP", "HIVD", "PrEPnHIVD"),
                           labels = c("S1: 1.6% PrEP;84% HIV treatment", 
                                      "S2: 20% PrEP;84% HIV treatment",
                                      "S3: 1.6% PrEP;95% HIV treatment",
                                      "S4: 20% PrEP;95% HIV treatment")))
plot_tag
plot_ihea1 <- test3%>%
  filter(year == 2090 & population == "Entire population of MSM")%>%
  ggplot(., aes(x= Casca, y = best/(Currency_converter*1000000), 
                fill = Casca)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  xlab("Testing") + ylab("Cost (Millions in USD$)") + 
  labs(fill = "HCV Cascade")  +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.background = element_rect(colour = "white"),
        plot.background = element_rect(colour = "white"),
        panel.border     = element_rect(fill = NA, colour = "black", 
                                        size = NA),
        plot.title = element_text(face = "bold",
                                  size = 14, hjust = 0.5),
        text = element_text(),
        axis.title = element_text(face = "bold",size = 14),
        axis.title.y = element_text(angle=90,vjust =1),
        axis.title.x = element_text(vjust = -0.2),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(
          face = "bold",size = 14,colour = "black"),
        strip.text.x = element_text(size=14, color="black",
                                    face="bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(face="bold", size= 14),
        plot.margin = unit(c(10,5,5,5),"mm")
        
  ) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  facet_grid(scenario ~ testing, scale = "free",
             labeller = label_wrap_gen(width=15)) + 
  scale_y_continuous(limits = c(0,20))
  

plot_ihea1

posterdt <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/11. Documents_Joyce's PhD/General/04. Conferences/iHEA2023"
ggsave(path = posterdt, file=paste0("totalcost.png"), 
       plot_ihea1, width = 9, 
       height = 6)  

#### heatmap for ICER 

CEA_table <- CEA_table%>%mutate(scenario = factor(scenario, 
                             levels = c("PrEP coverage and HIV treatment coverage remains unchanged",
                                        "PrEP coverage increases to 20% by 2030",
                                        "HIV treatment coverages increase to 95% by 2030",
                                        "PrEP coverage increases to 20% and HIV treatment increased to 95% by 2030")))
labdata <- CEA_table%>%filter(indicator != "Med" & population == "Entire population of MSM")

heatmap(A, Rowv = NA, Colv = NA, 
        col = colorRampPalette(brewer.pal(8,"Blues"))(4))
show_col(c("#F7FBFF" "#B8D5EA" "#4F9BCB" "#084594"))

heat_ICER <- ggplot(CEA_table%>%filter(indicator == "Med" & 
                            population == "Entire population of MSM"), 
       aes(x = testing, y = rev(scenario), fill = as.numeric(vv))) +
  geom_tile(colour = "white", linewidth = 1) +  
  scale_fill_gradient(low = "#4F9BCB",
                      
                      high = "#084594",
                      
                      guide = "legend") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_text(aes(label = labdata$vv), color = "white", size = 5) +
  theme(text = element_text(size = 18)
        ,axis.title = element_blank()
        ,axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) 

ggsave(path = posterdt, file=paste0("heat_ICER.png"), 
       heat_ICER, width = 9, 
       height = 6)  
######
test_dt <- list()

test_dtQALY <- list()

test_all <- list()

for(i in names(cost_incre)){ 
  test_dt[[i]] <- do.call("rbind", cost_incre[[i]])%>%
    mutate(indicator = "Incremental Cost")
  
  test_dtQALY[[i]] <- do.call("rbind", qaly_incre[[i]])%>%
    mutate(indicator = "QALY gained")

  test_all[[i]] <- rbind(test_dt[[i]], test_dtQALY[[i]])%>%
    arrange(year, population, testing)%>%select(year, population, testing, 
                                                indicator,best, colset)%>%
    gather(., "simulation", "estimate", c(best, colset))%>%
    spread(indicator, estimate)
  
}


# finding limits for axis

# roundup function 
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

maxValue_incre <- list()

test_wide <- list()
for(i in names(test_all)){ 
  
  maxValue_incre[[i]] <- test_all[[i]]%>%filter(year == 2090)%>%
    dplyr::group_by(population)%>%
    mutate(`Incremental Cost` = 
                 roundUpNice(max(abs(`Incremental Cost`), na.rm = TRUE)),
           `QALY gained` = roundUpNice(max(abs(`QALY gained`), na.rm = TRUE)))%>%
    filter(simulation == "best" & testing == "POC_RNA")%>%
    mutate(cost_lab = `Incremental Cost`/1000000, 
           population = factor(population,
                               levels = c("All", "HIV-PrEP", "HIV+d"),
                               label = c("Entire population of MSM", 
                                         "HIV-negative MSM on PrEP", 
                                         "MSM who are HIV diagnosed (and on treatment)")))%>%
    mutate(cost_labbrk = case_when(cost_lab <=10 ~cost_lab/10,
                                  cost_lab <=100 ~ 5,
                                  cost_lab <=1000 ~ 50,
                                  cost_lab <= 10000 ~ 500),
           QALY_labbrk = case_when(`QALY gained` <=100 ~`QALY gained`/10,
                                   `QALY gained` <=1000 ~ 100,
                                   `QALY gained` <=10000 ~ 1000,
                                   `QALY gained` <=100000 ~ 10000,
                                   `QALY gained` >100000 ~ 50000))%>%
    arrange(population)
  
  
  
  test_wide[[i]] <- test_all[[i]]%>%filter(year == 2090)%>%
    filter(!population %in% c("HIV-", "HIV+"))%>%
    mutate(`Incremental Cost` = `Incremental Cost`)%>%
    mutate(population = factor(population, 
                               levels = c("All", "HIV-PrEP", "HIV+d"),
                               label = c("Entire population of MSM", 
                                         "HIV-negative MSM on PrEP", 
                                         "MSM who are HIV diagnosed (and on treatment)")),
           testing = factor(testing, 
                            levels = c( "POC_antibody", "Reflex_RNA",
                                        "POC_RNA"), 
                            label = c("Point-of-care anibody testing",
                                      "Clinic-based reflex RNA testing",
                                      "Single visit point-of-care RNA testing"
                            )
           ))
  
  

}

cost_lab <- list()
cost_limit <- list()
cost_break <- list()
qaly_limit <- list()
qaly_break <- list()
for(i in names(test_all)){ 
  for(m in seq_along(unique(maxValue_incre$PrEP$population))){
    cost_limit[[i]][[m]] <- c(-1*maxValue_incre[[i]]$`Incremental Cost`[m], 
                             maxValue_incre[[i]]$`Incremental Cost`[m])
    
    cost_break[[i]][[m]] <- seq(-1*maxValue_incre[[i]]$`Incremental Cost`[m], 
                              maxValue_incre[[i]]$`Incremental Cost`[m], 
                              maxValue_incre[[i]]$cost_labbrk[m]*1000000)
    
    cost_lab[[i]][[m]] <- seq(-1*maxValue_incre[[i]]$cost_lab[m], 
                              maxValue_incre[[i]]$cost_lab[m], 
                              maxValue_incre[[i]]$cost_labbrk[m]*1000000)
    
    qaly_limit[[i]][[m]] <- c(-1*maxValue_incre[[i]]$`QALY gained`[m], 
                              maxValue_incre[[i]]$`QALY gained`[m])
    
    
    qaly_break[[i]][[m]] <- seq(-1*maxValue_incre[[i]]$`QALY gained`[m], 
                                maxValue_incre[[i]]$`QALY gained`[m], 
                                maxValue_incre[[i]]$QALY_labbrk[m])
    
  }
}

# function to show lenged of geom_abline as "---" rather than "\" 
GeomAbline$draw_key <- function(data, params, size) 
{
  segmentsGrob(0, 0.5, 1, 0.5,
               gp = gpar(col = alpha(data$colour, 
                                     data$alpha), 
                         lwd = 0.5* .pt, lty = data$linetype, 
                         lineend = "butt"))
}


#### save data for plotting PSA ####
dtforPSA <- test_wide
save(dtforPSA,
     file = file.path(outputdt, "dtforPSA.rda"))



#### plot PSA: total #### 
plot_PSA <- list()
dt_psa <- list()
colpal <- scales::viridis_pal()(3)
thres <- data.frame(thres = c("1GDP"), value = c(32811))
popSelect <- unique(test_wide$PrEPsame$population)
for(i in seq_along(test_all)){ 
  
    dt_psa[[i]] <- test_wide[[i]]%>%
      mutate(testing = factor(testing, labels =  c(
                                                   "Point-of-care antibody testing",
                                                   "Reflex RNA testing",
                                                   "Point-of-care RNA testing")))%>%
      group_split(population)
      
    names(dt_psa[[i]]) <- unique(test_wide[[i]]$population)
} 
xt <- list()

for(i in seq_along(dt_psa)){ 
  xt[[i]] <- ggplot(dt_psa[[i]][[3]], 
               aes(y = `Incremental Cost`, x = `QALY gained`)) +
    geom_point(aes(colour = testing), size = 0.8) + 
    scale_color_viridis(discrete = TRUE, option =  "D") +
    geom_hline(yintercept=0, color = "black", linewidth = 1) +
    geom_vline(xintercept=0, color = "black", linewidth = 1) +
    stat_ellipse(aes(colour = testing), 
                 alpha = 0.8,
                 show.legend = FALSE, 
                 level = 0.95) +
    labs(colour = "Testing", tag = plot_tag[i], x = "", y = "") + 
    ggtitle(scenario_lab[i])  +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(panel.background = element_rect(colour = "white"),
          plot.background = element_rect(colour = "white"),
          panel.border     = element_rect(fill = NA, colour = "black", 
                                          size = NA),
          plot.title = element_text(face = "bold",
                                    size = 12, hjust = 0.5),
          text = element_text(),
          axis.title = element_text(face = "bold",size = 14),
          axis.title.y = element_text(angle=90,vjust =1),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,
                                     face = "bold",size = 12, colour = "black"), 
          axis.text.y = element_text(
            face = "bold",size = 14,colour = "black"),
          strip.text.x = element_text(size=14, color="black",
                                      face="bold"),
          legend.text = element_text(size = 14, face = "bold"),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(face="bold", size= 14),
          plot.margin = unit(c(10,5,5,5),"mm")) + 
    guides(color =guide_legend(direction='horizontal', override.aes = list(size=2))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
    scale_y_continuous(limits = cost_limit[[i]][[3]],
                       breaks = cost_break[[i]][[3]],
                       labels = cost_break[[i]][[3]]/1000000) + 
    scale_x_continuous(limits = qaly_limit[[i]][[3]],
                       breaks = qaly_break[[i]][[3]]) +
    geom_abline(aes(slope = 32811, intercept = 0,linetype = "1GDP"), colour = "black") + 
    scale_linetype_manual(name = "Threshold", values = c(2), 
                          guide = guide_legend(override.aes = list(color = c("black")))) 
    
  
  }



gplotlist_psa <- ggarrange(plotlist = xt, common.legend = TRUE, legend = "top")

gplotlist_psa <- annotate_figure(gplotlist_psa ,
                bottom = text_grob("Incremental QALY gained",,
                                   hjust = 1, vjust = -1, x = 0.55,  
                                   size = 14, face ="bold"),
                left = text_grob("Incremental cost (Millions in US$)",  rot = 90, 
                                 size = 14, face = "bold")
)
ggsave(path = outputfig, file=paste0("PSA_HIVD.pdf"), 
       gplotlist_psa, width = 16, 
       height = 12, dpi = 300)

#### bugs here, couldn't loop ####
for(i in seq_along(test_all)){ 
  plot_PSA[[i]] <- list()
    for(n in 1: length(popSelect)){
    plot_PSA[[i]][[n]] <-  ggplot(dt_psa[[i]][[n]], 
                                  aes(y = `Incremental Cost`, x = `QALY gained`)) +
      geom_point(aes(colour = testing), size = 0.8) + 
      scale_color_viridis(discrete = TRUE, option =  "D") +
      geom_hline(yintercept=0, color = "black", linewidth = 1) +
      geom_vline(xintercept=0, color = "black", linewidth = 1) +
      stat_ellipse(aes(colour = testing), 
                   alpha = 0.8,
                   show.legend = FALSE, 
                   level = 0.95) +
      ylab("Incremental cost (Millions in US$)") + xlab("Incremental QALY gained") + 
      labs(colour = "Testing", tag = plot_tag[i]) + 
      ggtitle(scenario_lab[i])  +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(panel.background = element_rect(colour = "white"),
            plot.background = element_rect(colour = "white"),
            panel.border     = element_rect(fill = NA, colour = "black", 
                                            size = NA),
            plot.title = element_text(face = "bold",
                                      size = 14, hjust = 0.5),
            text = element_text(),
            axis.title = element_text(face = "bold",size = 14),
            axis.title.y = element_text(angle=90,vjust =1),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,
                                       face = "bold",size = 12, colour = "black"), 
            axis.text.y = element_text(
              face = "bold",size = 14,colour = "black"),
            strip.text.x = element_text(size=14, color="black",
                                        face="bold"),
            legend.text = element_text(size = 14, face = "bold"),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title = element_text(face="bold", size= 14),
            plot.margin = unit(c(10,5,5,5),"mm")) + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
      scale_y_continuous(limits = cost_limit[[i]][[n]],
                         breaks = cost_break[[i]][[n]],
                         labels = cost_break[[i]][[n]]/1000000) + 
      scale_x_continuous(limits = qaly_limit[[i]][[n]],
                         breaks = qaly_break[[i]][[n]])
    
    
    
    
    plot_PSA[[i]][[n]] <- plot_PSA[[i]][[n]] + 
      geom_abline(aes(slope = 32811, intercept = 0,linetype = "1GDP"), colour = "black") + 
      geom_abline(aes(slope = 32811*3, intercept = 0,linetype = "3GDP"), colour = "black" ) + 
      scale_linetype_manual(name = "Threshold", values = c(2, 3), 
                            guide = guide_legend(override.aes = list(color = c("black", "black"))))
  }
  
  names(plot_PSA[[i]]) <- popSelect
}

names(plot_PSA) <- names(test_all)

gplotlist_PSA <- ggarrange(plotlist = plot_PSA, common.legend = TRUE, 
                           ncol  = 2, nrow = 2) 

ggsave(path = outputfig, file=paste0("PSA_allMSM.pdf"), 
       gplotlist_PSA, width = 25, dpi = 800,
       height = 32)


#### plot PSA: total #### 
plot_PSA <- list()
colpal <- scales::viridis_pal()(3)
thres <- data.frame(thres = c("1GDP"), value = c(32811))
popSelect <- unique(test_wide$PrEPsame$population)
for(i in seq_along(test_all)){ 
  plot_PSA[[i]] <-  ggplot(test_wide[[i]], 
                           aes(y = `Incremental Cost`, x = `QALY gained`)) +
    geom_point(aes(colour = testing), size = 0.8) + 
    scale_color_viridis(discrete = TRUE, option =  "D") +
    geom_hline(yintercept=0, color = "black", linewidth = 1) +
    geom_vline(xintercept=0, color = "black", linewidth = 1) +
    stat_ellipse(aes(colour = testing), 
                 alpha = 0.8,
                 show.legend = FALSE, 
                 level = 0.95) +
    ylab("Incremental cost (Millions in US$)") + xlab("Lifetime QALY gained") + 
    labs(colour = "Testing", tag = plot_tag[i]) + 
    ggtitle(scenario_lab[i])  +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(panel.background = element_rect(colour = "white"),
          plot.background = element_rect(colour = "white"),
          panel.border     = element_rect(fill = NA, colour = "black", 
                                          size = NA),
          plot.title = element_text(face = "bold",
                                    size = 14, hjust = 0.5),
          text = element_text(),
          axis.title = element_text(face = "bold",size = 14),
          axis.title.y = element_text(angle=90,vjust =1),
          axis.title.x = element_text(vjust = -0.2),
          axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,
                                     face = "bold",size = 12, colour = "black"), 
          axis.text.y = element_text(
            face = "bold",size = 14,colour = "black"),
          strip.text.x = element_text(size=14, color="black",
                                      face="bold"),
          legend.text = element_text(size = 14, face = "bold"),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(face="bold", size= 14),
          plot.margin = unit(c(10,5,5,5),"mm")) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
    facet_custom(~population, ncol = 3, scales = "free", 
                 scale_overrides = 
                   list(scale_new(1, 
                                  scale_y_continuous(limits = cost_limit[[i]][[1]],
                                                     breaks = cost_break[[i]][[1]],
                                                     labels = cost_break[[i]][[1]]/1000000)),
                        scale_new(2,
                                  scale_y_continuous(limits = cost_limit[[i]][[2]],
                                                     breaks = cost_break[[i]][[2]],
                                                     labels = format(cost_break[[i]][[2]]/1000000, 
                                                                     digits =1))),
                        scale_new(3,
                                  scale_y_continuous(limits = cost_limit[[i]][[3]],
                                                     breaks = cost_break[[i]][[3]],
                                                     labels = cost_break[[i]][[3]]/1000000)),
                        scale_new(1,
                                  scale_x_continuous(limits = qaly_limit[[i]][[1]],
                                                     breaks = qaly_break[[i]][[1]])),
                        scale_new(2,
                                  scale_x_continuous(limits = qaly_limit[[i]][[2]],
                                                     breaks = qaly_break[[i]][[2]])),
                        scale_new(3,
                                  scale_x_continuous(limits = qaly_limit[[i]][[3]],
                                                     breaks = qaly_break[[i]][[3]])))
    )
    
  plot_PSA[[i]] <- plot_PSA[[i]] + 
    geom_abline(aes(slope = 32811, intercept = 0,linetype = "1GDP"), colour = "black") + 
    geom_abline(aes(slope = 32811*3, intercept = 0,linetype = "3GDP"), colour = "black" ) + 
    scale_linetype_manual(name = "Threshold", values = c(2, 3), 
                          guide = guide_legend(override.aes = list(color = c("black", "black"))))
  
  
  }


  gplotlist_PSA <- ggarrange(plotlist = plot_PSA, common.legend = TRUE, nrow  = 2) 
  
  ggsave(path = outputfig, file=paste0("PSA.pdf"), 
         gplotlist_PSA, width = 25, 
         height = 32)

  
#### PSA with density plot: test ####
  
  
  
#### for IAS poster #### 
 
  plot_PSA_test <- list()
  colpal <- scales::viridis_pal()(3)
  thres <- data.frame(thres = c("1GDP", "3GDP"), value = c(32811, 32811*3))
  
  for(i in seq_along(test_all)){ 
    plot_PSA_test[[i]] <-  ggplot(test_wide[[i]]%>%filter(population == "Entire population of MSM"), 
                             aes(y = `Incremental Cost`, x = `QALY gained`)) +
      geom_point(aes(colour = testing), size = 0.8) + 
      scale_color_viridis(discrete = TRUE, option =  "D") +
      geom_hline(yintercept=0, color = "black", linewidth = 1) +
      geom_vline(xintercept=0, color = "black", linewidth = 1) +
      stat_ellipse(color = "gray",
                   linetype = "dashed",
                   linewidth = 0.8,
                   alpha = 0.8,
                   show.legend = FALSE, 
                   level = 0.95) +
      ylab("Incremental cost (Millions in US$)") + xlab("Lifetime QALY gained") + 
      labs(colour = "Testing", tag = plot_tag[i]) + 
      ggtitle(scenario_lab[i])  +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(panel.background = element_rect(colour = "white"),
            plot.background = element_rect(colour = "white"),
            panel.border     = element_rect(fill = NA, colour = "black", 
                                            size = NA),
            plot.title = element_text(face = "bold",
                                      size = 14, hjust = 0.5),
            text = element_text(),
            axis.title = element_text(face = "bold",size = 14),
            axis.title.y = element_text(angle=90,vjust =1),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,
                                       face = "bold",size = 12, colour = "black"), 
            axis.text.y = element_text(
              face = "bold",size = 14,colour = "black"),
            strip.text.x = element_text(size=14, color="black",
                                        face="bold"),
            legend.text = element_text(size = 14, face = "bold"),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title = element_text(face="bold", size= 14),
            plot.margin = unit(c(10,5,5,5),"mm")) + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
      facet_custom(~testing, ncol = 3, scales = "free", 
                   scale_overrides = 
                     list(scale_new(1, 
                                    scale_y_continuous(limits = cost_limit[[i]][[1]]/5,
                                                       breaks = seq(-20000000, 20000000, 5000000),
                                                       labels = seq(-20, 20, 5))),
                          scale_new(2,
                                    scale_y_continuous(limits = cost_limit[[i]][[1]],
                                                       breaks = seq(-100000000, 100000000, 50000000),
                                                       labels = seq(-100, 100, 50))),
                          scale_new(3,
                                    scale_y_continuous(limits = cost_limit[[i]][[1]],
                                                       breaks = seq(-100000000, 100000000, 50000000),
                                                       labels = seq(-100, 100, 50))),
                          scale_new(1,
                                    scale_x_continuous(limits = qaly_limit[[i]][[1]]/4,
                                                       breaks = qaly_break[[i]][[1]]/4)),
                          scale_new(2,
                                    scale_x_continuous(limits = qaly_limit[[i]][[1]],
                                                       breaks = qaly_break[[i]][[1]])),
                          scale_new(3,
                                    scale_x_continuous(limits = qaly_limit[[i]][[1]],
                                                       breaks = qaly_break[[i]][[1]]))
                          
                          
                     )
      )
    
    
    
    plot_PSA_test [[i]] <-  plot_PSA_test [[i]] + 
      geom_abline(aes(slope = 32811, intercept = 0,linetype = "1GDP"), colour = "black") + 
      geom_abline(aes(slope = 32811*3, intercept = 0,linetype = "3GDP"), colour = "black" ) + 
      scale_linetype_manual(name = "Threshold", values = c(2, 3), 
                            guide = guide_legend(override.aes = list(color = c("black", "black"))))
  }
  
  plot_PSA_test[[1]]
  gplotlist_PSA <- ggarrange(plotlist = plot_PSA, common.legend = TRUE, ncol  = 1) 
  posterdt <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/11. Documents_Joyce's PhD/General/04. Conferences/IAS 2023"
  ggsave(path = posterdt, file=paste0("PSA_ellispe.pdf"), 
         plot_PSA_test[[1]], width = 16, 
         height = 9)  


  
  