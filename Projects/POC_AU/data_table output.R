# tidy up outcomes 
rm(list = ls())

project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Documents/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
# Load useful libraries

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")
library("ggplot2")
library("viridis") 
library("openxlsx")
library(gt)
library(dplyr)
Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")
OutputFig <- file.path(OutputFolder, "Figs")
load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "Res_flowcost.rda")))
load(file.path(OutputFolder, paste0(project_name,"epiRes_timestep.rda")))

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
endY <- 100
year_obs <- c(POC_AU$simY  +5 - 1  , POC_AU$simY + 10 - 1, POC_AU$simY + 20 - 1)

#### cumulative: number of infections: overall in 5, 10, 20y ####
# starting at year of 2022
# tidy up data in dt to yearly and cumulative 
Resflow_cum_all <- list()

col_lst <- colnames(Resflow_year_all$`Status quo`)[-c(1,2)]
for(i in names(Resflow_year_all)){ 
  Resflow_cum_all[[i]] <- Resflow_year_all[[i]]%>%
    mutate(year = year + POC_AU$cabY - 1)%>%
    filter(year>= POC_AU$simY)%>%
    group_by(year)%>%
    summarise_at(col_lst, "sum")%>%ungroup()%>%
    mutate(across(col_lst, list(cum=cumsum), .names = "{fn}_{col}"))
  
}

# avert number of transitions 
Resflow_cum_all_avert <- list()
col_lst <- colnames(Resflow_cum_all$`Status quo`)[-c(1)]
for(i in names(Resflow_cum_all)){ 
  Resflow_cum_all_avert[[i]] <- cbind(
    year = Resflow_cum_all[[i]]$year, 
    data.frame(Resflow_cum_all[[1]][, col_lst] - Resflow_cum_all[[i]][, col_lst]))%>%
    rename_at(vars(col_lst), ~ paste0("avert_", col_lst))
  
}

Resflow_sc_cum_all <- list()

col_lst <- colnames(Resflow_sc_year_all$`Status quo`)[-c(1,2)]
for(i in names(Resflow_year_all)){ 
  Resflow_sc_cum_all[[i]] <- Resflow_sc_year_all[[i]]%>%
    mutate(year = year + POC_AU$cabY - 1)%>%
    filter(year>= POC_AU$simY)%>%
    group_by(year)%>%
    summarise_at(col_lst, "sum")%>%ungroup()%>%
    mutate(across(col_lst, list(cum=cumsum), .names = "{fn}_{col}"))
  
}
# avert number of transitions 
Resflow_sc_cum_all_avert <- list()
col_lst <- colnames(Resflow_sc_cum_all$`Status quo`)[-c(1)]
for(i in names(Resflow_sc_cum_all)){ 
  Resflow_sc_cum_all_avert[[i]] <- cbind(
    year = Resflow_sc_cum_all[[i]]$year, 
    data.frame(Resflow_sc_cum_all[[1]][, col_lst] - Resflow_sc_cum_all[[i]][, col_lst]))%>%
    rename_at(vars(col_lst), ~ paste0("avert_", col_lst))
  
}


# function for generating plots 
plot_pocau <- function(pj,dt,indicator, type, year_obs){ 
  # type: new, cum, and avert  
  
  if(isTRUE(type == "new")){ 
    
    a <- dt%>%filter(index == indicator)%>%ggplot(data = ., 
                                                  aes(x = year, y = best, colour = scenario)) + 
      geom_line() + 
      scale_x_continuous(limits = c(pj$simY,2050), breaks = c(seq(pj$simY,2050, 5), 2050))+
      scale_color_viridis(discrete = TRUE,option = "D") + 
      theme_linedraw() + 
      ggtitle("Yearly number")
    
  }
  else if(isTRUE(type =="cum")){ 
    a <- dt%>%filter(index == indicator & year %in% year_obs)%>%
      ggplot(data =. , aes(x = as.character(year), y = best, fill = scenario)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      labs(x = "", y = indicator) + 
      scale_x_discrete(labels = c(paste0((year_obs - pj$simY + 1), "-Year", sep = ""))) + 
      scale_fill_viridis(discrete = TRUE,option = "D") + 
      theme_linedraw() +
      ggtitle("Cumulative number")
    
    
  }
  
  else if(isTRUE(type =="avert")){ 
    a <- dt%>%filter(index == indicator & year %in% year_obs)%>%
      ggplot(data =. , aes(x = as.character(year), y = best, fill = scenario)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      labs(x = "", y = indicator) + 
      scale_x_discrete(labels = c(paste0((year_obs - pj$simY + 1), "-Year", sep = ""))) + 
      scale_fill_viridis(discrete = TRUE,option = "D") + 
      theme_linedraw() +
      ggtitle("Number averted (cumulative)")
    
    
  }
  
  return(a)
}

# combind list to one dataframe 
# convert to long form 
convert_long <- function(dt){ 
  x <- dplyr::bind_rows(dt, .id = 'scenario')%>%
    mutate(scenario = factor(scenario, levels = c(names(dt)))) 
  x.long <- x%>%gather(index, best, -c(year, scenario))
  
  return(x.long)
}


floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)




# number of dc, hcc, lt 
# sum to yearly and cumulative & output long form 
num_sum <- function(dt, endY){ 
  dt_year <- list()
  
  for(i in names(dt)){ 
    
    dt_year[[i]] <- dt[[i]]%>%ungroup()%>%
      mutate(year = rep(seq(1, endY -1, 1), each = POC_AU$npops*(1/POC_AU$timestep)))%>%
      mutate(year = year + POC_AU$cabY - 1)%>%
      filter(year>= POC_AU$simY)%>%
      group_by(year)%>%
      summarise(best = mean(best))%>%ungroup()%>%
      mutate(cum_best = cumsum(best))
    
    
    dt_year[[i]] <- cbind(dt_year[[i]], 
                          avertcum_best = dt_year[[1]]$cum_best - dt_year[[i]]$cum_best)
    
  }
  
  dt_year_bind <- dplyr::bind_rows(dt_year, .id= 'scenario')%>%
    mutate(scenario = factor(scenario, levels = c(names(dt)))) 
  
  long.dt_year_bind <- dt_year_bind%>%gather(index, best, -c(scenario, year))
  
  return(long.dt_year_bind)
}

#####################################
# cost data plot 
# extract yearly value 
cap <- 200000000
cost_year_all <- list()


for(i in names(Rescost_yearcum_all)){ 
  
  cost_year_all[[i]] <- Rescost_yearcum_all[[i]]%>%group_by(year)%>%slice(n())%>%
    mutate(cost_totalDAAcap_ycum = ifelse(cost_totalDAA_ycum>=cap, cap, cost_totalDAA_ycum))%>%
    mutate(cost_totalDAAcap_disycum = ifelse(cost_totalDAAcap_ycum<cap, cost_totalDAAcap_ycum/discount, cap))%>%
    mutate(cost_totalcap_ycum = cost_compartment_ycum + cost_ab_ycum + 
             cost_RNA_ycum +  cost_POCT_ycum + cost_totalDAA_ycum + 
             cost_TreatOther_ycum + cost_RetreatOther_ycum + cost_Cured_ycum)%>%
    mutate(cost_totalcap_disycum = cost_compartment_disycum + cost_ab_disycum + 
             cost_RNA_disycum +  cost_POCT_disycum + cost_totalDAA_disycum + 
             cost_TreatOther_disycum + cost_RetreatOther_disycum + cost_Cured_disycum)%>%
    select(year, discount, matches('_ycum'), cost_totalcap_ycum, 
           cost_totalcap_disycum,cost_totalDAAcap_ycum, cost_totalDAAcap_disycum, 
           matches('_disycum'))%>%
    filter(year>= POC_AU$simY)
  
}

col_lst <- colnames(cost_year_all[[1]])[-c(1,2)]
cost_cum_all <- list()
cost_cumavert_all <- list()
for(i in names(cost_year_all)){ 
  cost_cum_all[[i]] <- cost_year_all[[i]]%>%ungroup()%>%
    mutate_at(col_lst, "cumsum")%>%
    select(year, discount ,col_lst)
  
  cost_cumavert_all[[i]] <- cbind(year =  cost_cum_all[[i]]$year, 
                                  discount = cost_cum_all[[i]]$year, 
                                  data.frame(cost_cum_all[[1]][, col_lst] - 
                                               cost_cum_all[[i]][, col_lst]))%>%
    rename_at(vars(col_lst), ~ paste0("avert_", col_lst))
}


cost_stage <- list()
cost_stage_avert <- list()


for(i in names(cost_cum_all)){ 
  
  cost_stage[[i]] <- cost_cum_all[[i]]%>%
    mutate(diagnosis_ycum = cost_ab_ycum + cost_RNA_ycum + cost_POCT_ycum,
           diagnosis_disycum = cost_ab_disycum + cost_RNA_disycum + cost_POCT_disycum,
           management_ycum = cost_TreatOther_ycum + cost_RetreatOther_ycum + 
             cost_compartment_ycum + cost_Cured_ycum,
           management_disycum = cost_TreatOther_disycum + cost_RetreatOther_disycum + 
             cost_compartment_disycum + cost_Cured_disycum,,
           treatment_ycum  = cost_totalDAA_ycum, 
           treatment_disycum = cost_totalDAA_disycum)%>%
    select(year, diagnosis_ycum, management_ycum,treatment_ycum, 
           diagnosis_disycum,management_disycum,treatment_disycum )
}

cost_stagesave <- list()

cost_stagesave_year_obs <- list()
col_lst <- c("diagnosis_ycum", "management_ycum", 
             "treatment_ycum","diagnosis_disycum", "management_disycum", 
             "treatment_disycum" )
for(i in names(cost_stage)){ 
  cost_stagesave_year_obs[[i]] <- cost_stage[[i]]%>%filter(year%in% year_obs)
  cost_stagesave[[i]] <- cbind(year = cost_stagesave_year_obs[[i]]$year, 
                               data.frame(cost_stagesave_year_obs[[1]][, col_lst] - 
                                            cost_stagesave_year_obs[[i]][, col_lst]))
  
}



##### Table generation #####
# number of indicators and cost 

# number of infections & death 
# select infections and death columns and combine each scenarios and year <= 2050 
library(stringr)

inf_death <- list()
Test_treat <- list()
for(i in names(Resflow_cum_all)){ 
  inf_death[[i]] <- Resflow_cum_all[[i]]%>%filter(year<= 2050)%>%
    select(year, newInfections, HCVdeath)%>%ungroup()%>%
    mutate(newInfections = round(newInfections, digits = 1),
           HCVdeath = round(HCVdeath, digits = 1))%>%
    gather(index, best, - c("year"))
  
  Test_treat[[i]] <- Resflow_cum_all[[i]]%>%filter(year<= 2050)%>%
    select(year, Treatment, Retreat, Treatment_sc, 
           Testing_ab, Testing_ab_neg, Testing_ab_sc, Testing_ab_sc_neg,
           Testing_RNA, Testing_RNA_neg, Testing_RNA_sc, Testing_RNA_sc_neg,
           Testing_POCT, Testing_POCT_neg, Testing_POCT_sc, Testing_POCT_sc_neg)%>%ungroup()%>%
    mutate(Treatment_all = (Treatment + Retreat + Treatment_sc), 
           Testing_ab_all = (Testing_ab + Testing_ab_neg + Testing_ab_sc + Testing_ab_sc_neg),
           Testing_RNA_all = (Testing_RNA + Testing_RNA_neg + Testing_RNA_sc + Testing_RNA_sc_neg),
           Testing_POCT_all = (Testing_POCT + Testing_POCT_neg + Testing_POCT_sc + Testing_POCT_sc_neg),
           Treatment_sc = (Treatment_sc), 
           Testing_ab_sc = (Testing_ab_sc + Testing_ab_sc_neg),
           Testing_RNA_sc = (Testing_RNA_sc + Testing_RNA_sc_neg),
           Testing_POCT_sc = (Testing_POCT_sc + Testing_POCT_sc_neg))%>%
    select(year, Treatment_all, Testing_ab_all, Testing_RNA_all, Testing_POCT_all,
           Treatment_sc, Testing_ab_sc, Testing_RNA_sc, Testing_POCT_sc)%>%
    gather(index, best, - c("year"))
  
  }

Num_dc_year_bind <- num_sum(Num_dc, endY =100)%>%mutate(best = round(best, digits = 1))%>%
  spread(scenario, best)%>%filter(year <= 2050)%>%filter(index == "best")%>%
  mutate( index = "DC")

Num_hcc_year_bind <- num_sum(Num_hcc, endY =100)%>%mutate(best = round(best, digits = 1))%>%
  spread(scenario, best)%>%filter(year <= 2050)%>%filter(index == "best")%>%
  mutate( index = "HCC")


Num_lt_year_bind <- num_sum(Num_lt, endY =100)%>%mutate(best = round(best, digits = 1))%>%
  spread(scenario, best)%>%filter(year <= 2050)%>%filter(index == "best")%>%
  mutate( index = "Liver transplant")




inf_death_bind <- dplyr::bind_rows(inf_death, .id = "scenario")%>%
  spread(scenario, best)

num_inf_adv <- rbind(Num_dc_year_bind, Num_hcc_year_bind,Num_lt_year_bind, inf_death_bind)

Test_treat_bind <- dplyr::bind_rows(Test_treat, .id = "scenario")%>%
  spread(scenario, best)%>%arrange(index)%>%
  select(year, index, names(Resflow_cum_all))%>%gt()

# for spreadsheet output 
Res_bind <- list()
Res_total <- list()
Res_sc <- list()
for(n in names(Resflow_cum_all)){ 
  Res_bind[[n]] <- Resflow_cum_all[[n]]%>%mutate(
    Treatment_all = (Treatment + Retreat + Treatment_sc), 
    Testing_ab_all = (Testing_ab + Testing_ab_neg + Testing_ab_sc + Testing_ab_sc_neg),
    Testing_RNA_all = (Testing_RNA + Testing_RNA_neg + Testing_RNA_sc + Testing_RNA_sc_neg),
    Testing_POCT_all = (Testing_POCT + Testing_POCT_neg + Testing_POCT_sc + Testing_POCT_sc_neg),
    Testing_ab_sc = (Testing_ab_sc + Testing_ab_sc_neg),
    Testing_RNA_sc = (Testing_RNA_sc + Testing_RNA_sc_neg),
    Testing_POCT_sc = (Testing_POCT_sc + Testing_POCT_sc_neg),
  )
 
  Res_total[[n]] <- Res_bind[[n]]%>%select(year, Testing_ab_all, Testing_RNA_all, 
                                           Testing_POCT_all, Treatment_all, Cured)
  
  Res_sc[[n]] <- Res_bind[[n]]%>%select(year, Testing_ab_sc, Testing_RNA_sc, 
                                           Testing_POCT_sc, Treatment_sc)
  
  
  }
Res_bind$`Status quo`$Cured
View(Res_bind$`Status quo`%>%select(year,Treatment_all, Testing_ab_all, Testing_RNA_all,
                                    Testing_POCT_all))

Res_total$

str(Res_bind)
Resflow_cum_all[[1]]
Resflow_cum_all[[1]]
Res_bind[[1]]

Resflow_cum_all$dfList_NP_2023%>%mutate(
  Treatment_all = (Treatment + Retreat + Treatment_sc), 
  Testing_ab_all = (Testing_ab + Testing_ab_neg + Testing_ab_sc + Testing_ab_sc_neg),
  Testing_RNA_all = (Testing_RNA + Testing_RNA_neg + Testing_RNA_sc + Testing_RNA_sc_neg),
  Testing_POCT_all = (Testing_POCT + Testing_POCT_neg + Testing_POCT_sc + Testing_POCT_sc_neg),
  Testing_ab_sc = (Testing_ab_sc + Testing_ab_sc_neg),
  Testing_RNA_sc = (Testing_RNA_sc + Testing_RNA_sc_neg),
  Testing_POCT_sc = (Testing_POCT_sc + Testing_POCT_sc_neg),
)



for(i in names(Resflow_cum_all)){ 
  Res_sheet_sc[[i]] <- Res_sheet[[i]]%>% 
    select(year, Treatment_sc, Testing_ab_sc, Testing_RNA_sc, Testing_POCT_sc)
    
} 
cost_year_all
View($dfList_NP_2023)
Res_cost_all <- list()
Res_cost_sc <- list()
for(i in names(Resflow_cum_all)){
  Res_cost_all[[i]] <- cbind(
    year = Res_bind[[i]]$year,
    cost_total = cost_year_all[[i]]$cost_total_ycum,
    num_ab = Res_bind[[i]]$Testing_ab_all, 
    cost_ab = cost_year_all[[i]]$cost_ab_ycum + cost_year_all[[i]]$cost_ab_sc_ycum, 
    num_RNA = Res_bind[[i]]$Testing_RNA_all, 
    cost_RNA = cost_year_all[[i]]$cost_RNA_ycum + cost_year_all[[i]]$cost_RNA_sc_ycum,
    num_POCT = Res_bind[[i]]$Testing_POCT_all, 
    cost_POCT = cost_year_all[[i]]$cost_POCT_ycum + cost_year_all[[i]]$cost_POCT_sc_ycum,
    num_Treat = Res_bind[[i]]$Treatment_all, 
    cost_DAA = cost_year_all[[i]]$cost_totalDAA_ycum, 
    cost_DAAcap = cost_year_all[[i]]$cost_totalDAAcap_ycum,
    cost_TreatOther = cost_year_all[[i]]$cost_TreatOther_ycum + cost_year_all[[i]]$cost_RetreatOther_ycum,
    num_Cured = Res_bind[[i]]$Cured,
    cost_cured = cost_year_all[[i]]$cost_Cured_ycum,
    cost_compartment = cost_year_all[[i]]$cost_compartment_ycum)%>%
    as_tibble()%>%
    mutate(cost_total_cap = cost_ab + cost_RNA + cost_POCT + cost_DAAcap + 
             cost_TreatOther + cost_cured + cost_compartment)
  
  
  
  Res_cost_sc[[i]] <- cbind(
    year = Res_bind[[i]]$year,
    num_ab_sc = Res_bind[[i]]$Testing_ab_sc,
    cost_ab_sc = cost_year_all[[i]]$cost_ab_sc_ycum,
    num_RNA_sc = Res_bind[[i]]$Testing_RNA_sc, 
    cost_RNA_sc =  cost_year_all[[i]]$cost_RNA_sc_ycum,
    num_POCT_sc = Res_bind[[i]]$Testing_POCT_sc, 
    cost_POCT_sc = cost_year_all[[i]]$cost_POCT_sc_ycum,
    num_Treat_sc = Res_bind[[i]]$Treatment_sc, 
    cost_TreatOther_sc = cost_year_all[[i]]$cost_TreatOther_sc_ycum)%>%
    as_tibble()%>%
    mutate(NP_cost = cost_ab_sc + cost_RNA_sc + cost_POCT_sc + cost_TreatOther_sc)%>%
    select(year, NP_cost, num_ab_sc, cost_ab_sc, num_RNA_sc, cost_RNA_sc, num_POCT_sc, 
           cost_POCT_sc, num_Treat_sc, cost_TreatOther_sc)
  
}

write.xlsx(Res_cost_sc, file = file.path(OutputFolder, paste0("POC_AU_Res_cost_NP.xlsx")), 
           append=TRUE) 
write.xlsx(Res_cost_all, file = file.path(OutputFolder, paste0("POC_AU_Res_cost_total.xlsx")), 
           append=TRUE)

