# this script imports the datasets tidy up in Resultssummary.R 
# {Res_dt.rda} and {epiRes_timestep.rda} 
# we aggregate the number to annual numbers and apply cost discount in this script. 
# three datasets we work on in this script: {Num_box}, {Resflow_dt}, {Rescost_dt}
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
library("openxlsx")

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")
# project specific code path 
Proj_code <- file.path(codefun_path, paste0("projects/", project_name))

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "Res_dt.rda")))
load(file.path(OutputFolder, paste0(project_name, "epiRes_timestep.rda")))
source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Proj_code, "/model_timestep.R")) 

AUdiscount <- 0.05
unitC_DAA <- 21085
cap <- 200000000
endY <- 100
# summarize the number of flow for each pop in each year 
Resflow_year_pop <- list()
for(i in names(Resflow_dt)){ 
  Resflow_year_pop[[i]] <- Resflow_dt[[i]]%>%as_tibble()%>%arrange(population, year)%>%
    group_by(year, population)%>%
    mutate(newInf_year = cumsum(newInfections), 
           HCVdeath_year = cumsum(HCVdeath), 
           Treatment_year = cumsum(Treatment), 
           Retreat_year = cumsum(Retreat), 
           Testing_ab_year = cumsum(Testing_ab), 
           Testing_RNA_year = cumsum(Testing_RNA), 
           Testing_POCT_year = cumsum(Testing_POCT), 
           Cured_year = cumsum(Cured))%>%
    # finding the last row in each year and pop
    mutate(index = c(round((timestep%%1)/POC_AU$timestep, digits = 0)))%>%
    mutate(newInf_ycum = ifelse(index == 11, newInf_year, NA), 
           HCVdeath_ycum = ifelse(index == 11, HCVdeath_year, NA),
           Treatment_ycum = ifelse(index == 11, Treatment_year, NA),
           Retreat_ycum = ifelse(index == 11, Retreat_year, NA),
           Testing_ab_ycum = ifelse(index == 11, Testing_ab_year, NA),
           Testing_RNA_ycum = ifelse(index == 11, Testing_RNA_year, NA),
           Testing_POCT_ycum = ifelse(index == 11,  Testing_POCT_year, NA),
           Cured_ycum = ifelse(index == 11,  Cured_year, NA),
           )%>%
    select(year, timestep, population, newInfections, newInf_ycum, 
           HCVdeath, HCVdeath_ycum, Treatment, Treatment_ycum, 
           Retreat, Retreat_ycum, Cured, Cured_ycum, 
           Testing_ab, Testing_ab_ycum, Testing_RNA, Testing_RNA_ycum, 
           Testing_POCT, Testing_POCT_ycum)%>%
    arrange(timestep)

}


Resflow_year_all <- list()
colnames(Resflow_year_pop[[1]])
for(i in names(Resflow_dt)){ 
  Resflow_year_all[[i]] <- Resflow_year_pop[[i]]%>%ungroup()%>%
    arrange(timestep)%>%group_by(timestep)%>%
    summarise(
      newInfections = sum(newInfections, na.rm = TRUE), 
      HCVdeath = sum(HCVdeath, na.rm = TRUE), 
      Treatment = sum(Treatment, na.rm = TRUE),
      Retreat = sum(Retreat, na.rm = TRUE), 
      Cured = sum(Cured, na.rm = TRUE),
      Testing_ab = sum(Testing_ab, na.rm = TRUE),
      Testing_RNA = sum(Testing_RNA, na.rm = TRUE), 
      Testing_POCT = sum(Testing_POCT, na.rm = TRUE))%>%
    arrange(timestep)%>%
    mutate(year = rep(seq(1, endY - 1, 1),each =(1/POC_AU$timestep)))%>%
    select(year, timestep, newInfections, HCVdeath, Treatment, Retreat, Cured, Testing_ab, 
           Testing_RNA, Testing_POCT)
}

# separate cost of DAA and other costs related to treatment initiation  


# treatment and retreat 
RescostDAA <- list()

for(i in names(Rescost_dt)){ 
  Rescost_dt[[i]] <- Rescost_dt[[i]]%>%arrange(timestep)
  RescostDAA[[i]] <- cbind(year = Resflow_year_pop[[i]]$year, 
                                   timestep = Rescost_dt[[i]]$timestep,
                                   population = Rescost_dt[[i]]$population,
                                   newTreat = Resflow_year_pop[[i]]$Treatment, 
                                   costTreat = Rescost_dt[[i]]$cost_Treatment,
                                   newRetreat = Resflow_year_pop[[i]]$Retreat,
                                   costRetreat = Rescost_dt[[i]]$cost_Retreat)%>%
    as_tibble()%>%
    mutate(TreatDAA_cost = newTreat*unitC_DAA,
           RetreatDAA_cost = newRetreat*unitC_DAA)%>%
    mutate(TreatOther_cost = costTreat - TreatDAA_cost, 
           RetreatOther_cost = costRetreat - RetreatDAA_cost)%>%
    mutate(totalDAA_cost = TreatDAA_cost + RetreatDAA_cost)%>%
    ungroup()
  
}

RescostDAA_totalpop <- list()


for(i in names(RescostDAA)){
  
  RescostDAA_totalpop[[i]] <- RescostDAA[[i]]%>%ungroup()%>%
    arrange(timestep)%>%group_by(timestep)%>%
    summarise(
      newTreat = sum(newTreat, na.rm = TRUE), 
      TreatDAA_cost = sum(TreatDAA_cost, na.rm = TRUE),
      TreatOther_cost = sum(TreatOther_cost, na.rm = TRUE), 
      newRetreat = sum(newTreat, na.rm = TRUE),
      RetreatDAA_cost = sum(RetreatDAA_cost, na.rm = TRUE), 
      RetreatOther_cost = sum(RetreatOther_cost, na.rm = TRUE), 
      totalDAA_cost = sum(totalDAA_cost, na.rm = TRUE))%>%ungroup()
}

# cost_total DAA, exculdeDAA replaced cost_Treatment & cost_Retreat 
# valide with cost_total, 
# calculating total_cost_cap  
# extract the columns to replace in cost_treatment& cost_retreat 
# "cost_totalDAA_y": "cost_totalDAAcap_discum"
rep_col <- c("TreatDAA_cost", "TreatOther_cost", "RetreatDAA_cost", 
             "RetreatOther_cost", "totalDAA_cost")

Rescost_DAA_dt <- list()

for(i in names(Rescost_dt)){ 
  Rescost_DAA_dt[[i]] <- Rescost_dt[[i]]%>%
    select(-c("cost_Treatment", 
              "cost_Retreat"))
  Rescost_DAA_dt[[i]] <- 
    cbind(Rescost_DAA_dt[[i]], 
          RescostDAA[[i]][, (names(RescostDAA[[i]]) %in% rep_col)])
}

# summarize the costs in each years by pops
Rescost_year <- list()

for(i in names(Rescost_dt)){ 
  Rescost_year[[i]] <- Rescost_DAA_dt[[i]]%>%as_tibble()%>%
    mutate(year = c(rep(seq(1, endY - 2, 1), each = POC_AU$npops*1/POC_AU$timestep),
                    c(rep(endY - 1, POC_AU$npops*(1/POC_AU$timestep)))))%>%
    group_by(year, population)%>%
    mutate(cost_total_year = cumsum(cost_total), 
           cost_compartment_year = cumsum(cost_compartment), 
           cost_ab_year = cumsum(cost_ab), 
           cost_RNA_year = cumsum(cost_RNA), 
           cost_POCT_year = cumsum(cost_POCT), 
           cost_TreatDAA_year = cumsum(TreatDAA_cost), 
           cost_RetreatDAA_year = cumsum(RetreatDAA_cost), 
           cost_TreatOther_year = cumsum(TreatOther_cost),
           cost_RetreatOther_year = cumsum(RetreatOther_cost),
           cost_totalDAA_year = cumsum(totalDAA_cost),
           cost_Cured_year = cumsum(cost_Cured))%>%
    # finding the last row in each year and pop
    mutate(index = c(round((timestep%%1)/POC_AU$timestep, digits = 0)))%>%
    mutate(cost_total_ycum = ifelse(index == 11, cost_total_year, NA), 
           cost_compartment_ycum = ifelse(index == 11, cost_compartment_year, NA),
           cost_ab_ycum = ifelse(index == 11, cost_ab_year, NA),
           cost_RNA_ycum = ifelse(index == 11, cost_RNA_year, NA),
           cost_POCT_ycum = ifelse(index == 11, cost_POCT_year, NA),
           cost_TreatDAA_ycum = ifelse(index == 11, cost_TreatDAA_year, NA), 
           cost_RetreatDAA_ycum = ifelse(index == 11, cost_RetreatDAA_year, NA),
           cost_TreatOther_ycum = ifelse(index == 11, cost_TreatOther_year, NA), 
           cost_RetreatOther_ycum = ifelse(index == 11, cost_RetreatOther_year, NA),
           cost_totalDAA_ycum = ifelse(index == 11, cost_totalDAA_year, NA),
           cost_Cured_ycum = ifelse(index == 11,  cost_Cured_year, NA)
    )%>%
    mutate(cost_TreatDAA = TreatDAA_cost, 
           cost_RetreatDAA = RetreatDAA_cost, 
           cost_TreatOther = TreatOther_cost, 
           cost_RetreatOther = RetreatOther_cost, 
           cost_totalDAA = totalDAA_cost)%>%
    select(year, timestep, population, cost_total, cost_total_ycum, 
           cost_compartment, cost_compartment_ycum, cost_ab, cost_ab_ycum, 
           cost_RNA, cost_RNA_ycum, cost_POCT, cost_POCT_ycum, 
           cost_TreatDAA, cost_TreatDAA_ycum, 
           cost_RetreatDAA, cost_RetreatDAA_ycum, 
           cost_TreatOther, cost_TreatOther_ycum, 
           cost_RetreatOther, cost_RetreatOther_ycum, 
           cost_totalDAA, cost_totalDAA_ycum,
           cost_Cured, cost_Cured_ycum)
}


# cumulative yearly cost by population  & discount value 
Rescost_yearcum_pop <- list()

for(i in names(Rescost_dt)){ 
  Rescost_yearcum_pop[[i]] <- Rescost_year[[i]]%>%as_tibble()%>%arrange(population)%>%
   group_by(population)%>%
    mutate(cost_total_cum = cumsum(cost_total), 
           cost_compartment_cum = cumsum(cost_compartment), 
           cost_ab_cum = cumsum(cost_ab), 
           cost_RNA_cum = cumsum(cost_RNA),
           cost_POCT_cum = cumsum(cost_POCT),
           cost_TreatDAA_cum = cumsum(cost_TreatDAA), 
           cost_RetreatDAA_cum = cumsum(cost_RetreatDAA),
           cost_TreatOther_cum = cumsum(cost_TreatOther), 
           cost_RetreatOther_cum = cumsum(cost_RetreatOther), 
           cost_totalDAA_cum = cumsum(cost_totalDAA), 
           cost_Cured_cum = cumsum(cost_Cured))%>%
    mutate(cost_total_cum = ifelse(is.na(cost_total_ycum), NA, cost_total_cum),
           cost_compartment_cum = ifelse(is.na(cost_compartment_ycum), NA, cost_compartment_cum),
           cost_ab_cum = ifelse(is.na(cost_ab_ycum), NA, cost_ab_cum),
           cost_RNA_cum = ifelse(is.na(cost_RNA_ycum), NA, cost_RNA_cum),
           cost_POCT_cum = ifelse(is.na(cost_POCT_ycum), NA, cost_POCT_cum),
           cost_TreatDAA_cum = ifelse(is.na(cost_TreatDAA_ycum), NA, cost_TreatDAA_cum),
           cost_RetreatDAA_cum = ifelse(is.na(cost_RetreatDAA_ycum), NA, cost_RetreatDAA_cum),
           cost_TreatOther_cum = ifelse(is.na(cost_TreatOther_ycum), NA, cost_TreatOther_cum),
           cost_RetreatOther_cum = ifelse(is.na(cost_RetreatOther_ycum), NA, cost_RetreatOther_cum),
           cost_totalDAA_cum = ifelse(is.na(cost_totalDAA_ycum), NA, cost_totalDAA_cum),
           cost_Cured_cum = ifelse(is.na(cost_Cured_ycum), NA, cost_Cured_cum)
           )%>%
    # discount index 
    mutate(year = year + POC_AU$cabY - 1,
           id = year - POC_AU$simY, 
           discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))%>%
    # discount yearly value 
    mutate(cost_total_disycum = ifelse(is.na(cost_total_ycum)|is.na(discount), 0, cost_total_ycum/discount), 
           cost_compartment_disycum = ifelse(is.na(cost_compartment_ycum)|is.na(discount), 0, cost_compartment_ycum/discount), 
           cost_ab_disycum  = ifelse(is.na(cost_ab_ycum)|is.na(discount), 0, cost_ab_ycum/discount), 
           cost_RNA_disycum  = ifelse(is.na(cost_RNA_ycum)|is.na(discount), 0, cost_RNA_ycum/discount), 
           cost_POCT_disycum  = ifelse(is.na(cost_POCT_ycum)|is.na(discount), 0, cost_POCT_ycum/discount), 
           cost_TreatDAA_disycum  = ifelse(is.na(cost_TreatDAA_ycum)|is.na(discount), 0, cost_TreatDAA_ycum/discount), 
           cost_RetreatDAA_disycum  = ifelse(is.na(cost_RetreatDAA_ycum)|is.na(discount), 0, cost_RetreatDAA_ycum/discount), 
           cost_TreatOther_disycum  = ifelse(is.na(cost_TreatOther_ycum)|is.na(discount), 0, cost_TreatOther_ycum/discount), 
           cost_RetreatOther_disycum  = ifelse(is.na(cost_RetreatOther_ycum)|is.na(discount), 0, cost_RetreatOther_ycum/discount),
           cost_totalDAA_disycum  = ifelse(is.na(cost_totalDAA_ycum)|is.na(discount), 0, cost_totalDAA_ycum/discount), 
           cost_Cured_disycum = ifelse(is.na(cost_Cured_ycum)|is.na(discount), 0, cost_Cured_ycum/discount)
           )%>%
    # discount cumulative value 
    mutate(cost_total_discum = cumsum(cost_total_disycum), 
           cost_compartment_discum = cumsum(cost_compartment_disycum), 
           cost_ab_discum  = cumsum(cost_ab_disycum), 
           cost_RNA_discum  = cumsum(cost_RNA_disycum), 
           cost_POCT_discum  = cumsum(cost_POCT_disycum), 
           cost_TreatDAA_discum  = cumsum(cost_TreatDAA_disycum), 
           cost_RetreatDAA_discum  = cumsum(cost_RetreatDAA_disycum), 
           cost_TreatOther_discum  = cumsum(cost_TreatOther_disycum), 
           cost_RetreatOther_discum  = cumsum(cost_RetreatOther_disycum), 
           cost_totalDAA_discum  = cumsum(cost_totalDAA_disycum),
           cost_Cured_discum = cumsum(cost_Cured_disycum)
    )%>%
    select(year, timestep, population,discount, 
           
           cost_total, cost_total_ycum, cost_total_cum, cost_total_disycum, 
           cost_total_discum, 
           
           cost_compartment, cost_compartment_ycum, cost_compartment_cum, 
           cost_compartment_disycum, cost_compartment_discum,
           
           cost_ab, cost_ab_ycum, cost_ab_cum, cost_ab_disycum, cost_ab_discum, 
           
           cost_RNA, cost_RNA_ycum, cost_RNA_cum, cost_RNA_disycum, 
           cost_RNA_discum, 
           
           cost_POCT, cost_POCT_ycum, cost_POCT_cum, cost_POCT_disycum, 
           cost_POCT_discum, 
           
           cost_TreatDAA, cost_TreatDAA_ycum, cost_TreatDAA_cum, 
           cost_TreatDAA_disycum, cost_TreatDAA_discum,
           
           cost_RetreatDAA, cost_RetreatDAA_ycum, cost_RetreatDAA_cum, 
           cost_RetreatDAA_disycum, cost_RetreatDAA_discum,
           
           cost_TreatOther, cost_TreatOther_ycum, cost_TreatOther_cum, 
           cost_TreatOther_disycum, cost_TreatOther_discum,
           
           cost_RetreatOther, cost_RetreatOther_ycum, cost_RetreatOther_cum, 
           cost_RetreatOther_disycum, cost_RetreatOther_discum,
           
           cost_totalDAA, cost_totalDAA_ycum, cost_totalDAA_cum, 
           cost_totalDAA_disycum, cost_totalDAA_discum,
           
           cost_Cured, cost_Cured_ycum, cost_Cured_cum, cost_Cured_disycum, 
           cost_Cured_discum
           )
  }

# summarize the costs in each years overall 
Rescost_year_all <- list()

for(i in names(Rescost_dt)){ 
  Rescost_year_all[[i]] <- Rescost_year[[i]]%>%ungroup()%>%group_by(timestep)%>%
    summarise(cost_total = sum(cost_total),
              cost_compartment = sum(cost_compartment), 
              cost_ab = sum(cost_ab), 
              cost_RNA = sum(cost_RNA), 
              cost_POCT = sum(cost_POCT), 
              cost_TreatDAA = sum(cost_TreatDAA),
              cost_RetreatDAA = sum(cost_RetreatDAA),
              cost_TreatOther = sum(cost_TreatOther),
              cost_RetreatOther = sum(cost_RetreatOther),
              cost_totalDAA = sum(cost_totalDAA),
              cost_Cured = sum(cost_Cured))%>%
    ungroup()%>%
    mutate(year = c(rep(seq(1, endY - 2, 1), each = 1/POC_AU$timestep),
                    c(rep(endY - 1, (1/POC_AU$timestep)))))%>%
    group_by(year)%>%
    mutate(cost_total_year = cumsum(cost_total),
           cost_compartment_year = cumsum(cost_compartment), 
           cost_ab_year = cumsum(cost_ab), 
           cost_RNA_year = cumsum(cost_RNA), 
           cost_POCT_year = cumsum(cost_POCT), 
           cost_TreatDAA_year = cumsum(cost_TreatDAA), 
           cost_RetreatDAA_year = cumsum(cost_RetreatDAA), 
           cost_TreatOther_year = cumsum(cost_TreatOther),
           cost_RetreatOther_year = cumsum(cost_RetreatOther),
           cost_totalDAA_year = cumsum(cost_totalDAA),
           cost_Cured_year = cumsum(cost_Cured))%>%
    # finding the last row in each year and pop
    mutate(index = c(round((timestep%%1)/POC_AU$timestep, digits = 0)))%>%
    mutate(cost_total_ycum = ifelse(index == 11, cost_total_year, NA), 
           cost_compartment_ycum = ifelse(index == 11, cost_compartment_year, NA),
           cost_ab_ycum = ifelse(index == 11, cost_ab_year, NA),
           cost_RNA_ycum = ifelse(index == 11, cost_RNA_year, NA),
           cost_POCT_ycum = ifelse(index == 11, cost_POCT_year, NA),
           cost_TreatDAA_ycum = ifelse(index == 11, cost_TreatDAA_year, NA), 
           cost_RetreatDAA_ycum = ifelse(index == 11, cost_RetreatDAA_year, NA),
           cost_TreatOther_ycum = ifelse(index == 11, cost_TreatOther_year, NA), 
           cost_RetreatOther_ycum = ifelse(index == 11, cost_RetreatOther_year, NA),
           cost_totalDAA_ycum = ifelse(index == 11, cost_totalDAA_year, NA),
           cost_Cured_ycum = ifelse(index == 11,  cost_Cured_year, NA)
    )%>%
    mutate(cost_totalDAA_ycumCap = ifelse(cost_totalDAA_ycum>= cap, cap,cost_totalDAA_ycum))%>%
    mutate(cost_total_ycumCap = cost_compartment_ycum + cost_ab_ycum + cost_RNA_ycum + cost_POCT_ycum + 
             cost_totalDAA_ycumCap + cost_TreatOther_ycum + cost_RetreatOther_ycum + cost_Cured_ycum)%>%
    
    select(year, timestep, cost_total, cost_total_ycum, cost_total_ycumCap,
           cost_compartment, cost_compartment_ycum, cost_ab, cost_ab_ycum, 
           cost_RNA, cost_RNA_ycum, cost_POCT, cost_POCT_ycum, 
           cost_TreatDAA, cost_TreatDAA_ycum, 
           cost_RetreatDAA, cost_RetreatDAA_ycum, 
           cost_TreatOther, cost_TreatOther_ycum, 
           cost_RetreatOther, cost_RetreatOther_ycum, 
           cost_totalDAA, cost_totalDAA_ycum, cost_totalDAA_ycumCap,
           cost_Cured, cost_Cured_ycum)
}


# cumulative cost each year for overall pops
Rescost_yearcum_all <- list()

for(i in names(Rescost_dt)){ 
  Rescost_yearcum_all[[i]] <- Rescost_year_all[[i]]%>%ungroup()%>%
    mutate(cost_total_cum = cumsum(cost_total), 
           cost_compartment_cum = cumsum(cost_compartment), 
           cost_ab_cum = cumsum(cost_ab), 
           cost_RNA_cum = cumsum(cost_RNA),
           cost_POCT_cum = cumsum(cost_POCT), 
           cost_TreatDAA_cum = cumsum(cost_TreatDAA), 
           cost_RetreatDAA_cum = cumsum(cost_RetreatDAA), 
           cost_TreatOther_cum = cumsum(cost_TreatOther),
           cost_RetreatOther_cum = cumsum(cost_RetreatOther),
           cost_totalDAA_cum = cumsum(cost_totalDAA), 
           cost_Cured_cum = cumsum(cost_Cured))%>%
    mutate(cost_total_cum = ifelse(is.na(cost_total_ycum), NA, cost_total_cum),
           cost_compartment_cum = ifelse(is.na(cost_compartment_ycum), NA, cost_compartment_cum),
           cost_ab_cum = ifelse(is.na(cost_ab_ycum), NA, cost_ab_cum),
           cost_RNA_cum = ifelse(is.na(cost_RNA_ycum), NA, cost_RNA_cum),
           cost_POCT_cum = ifelse(is.na(cost_POCT_ycum), NA, cost_POCT_cum),
           cost_TreatDAA_cum = ifelse(is.na(cost_TreatDAA_ycum), NA, cost_TreatDAA_cum), 
           cost_RetreatDAA_cum = ifelse(is.na(cost_RetreatDAA_ycum), NA, cost_RetreatDAA_cum), 
           cost_TreatOther_cum = ifelse(is.na(cost_TreatOther_ycum), NA, cost_TreatOther_cum), 
           cost_RetreatOther_cum = ifelse(is.na(cost_RetreatOther_ycum), NA, cost_RetreatOther_cum), 
           cost_totalDAA_ycum = ifelse(is.na(cost_totalDAA_ycum), NA, cost_totalDAA_ycum), 
           cost_Cured_cum = ifelse(is.na(cost_Cured_ycum), NA, cost_Cured_cum)
    )%>%
    # discount index 
    mutate(year = year + POC_AU$cabY - 1,
           id = year - POC_AU$simY, 
           discount = ifelse(id>=0, (1 + AUdiscount)^id, NA))%>%
    
    # discount yearly value 
    mutate(cost_total_disycum = ifelse(is.na(cost_total_ycum)|is.na(discount), 0, cost_total_ycum/discount), 
           cost_compartment_disycum = ifelse(is.na(cost_compartment_ycum)|is.na(discount), 0, cost_compartment_ycum/discount), 
           cost_ab_disycum  = ifelse(is.na(cost_ab_ycum)|is.na(discount), 0, cost_ab_ycum/discount), 
           cost_RNA_disycum  = ifelse(is.na(cost_RNA_ycum)|is.na(discount), 0, cost_RNA_ycum/discount), 
           cost_POCT_disycum  = ifelse(is.na(cost_POCT_ycum)|is.na(discount), 0, cost_POCT_ycum/discount), 
           cost_TreatDAA_disycum = ifelse(is.na(cost_TreatDAA_ycum)|is.na(discount), 0, cost_TreatDAA_ycum/discount),
           cost_RetreatDAA_disycum = ifelse(is.na(cost_RetreatDAA_ycum)|is.na(discount), 0, cost_RetreatDAA_ycum/discount),
           cost_TreatOther_disycum = ifelse(is.na(cost_TreatOther_ycum)|is.na(discount), 0, cost_TreatOther_ycum/discount),
           cost_RetreatOther_disycum = ifelse(is.na(cost_RetreatOther_ycum)|is.na(discount), 0, cost_RetreatOther_ycum/discount),
           cost_totalDAA_disycum = ifelse(is.na(cost_totalDAA_ycum)|is.na(discount), 0,  cost_totalDAA_ycum/discount),
           cost_Cured_disycum = ifelse(is.na(cost_Cured_ycum)|is.na(discount), 0, cost_Cured_ycum/discount)
    )%>%
    # discount cumulative value 
    mutate(cost_total_discum = cumsum(cost_total_disycum), 
           cost_compartment_discum = cumsum(cost_compartment_disycum), 
           cost_ab_discum  = cumsum(cost_ab_disycum), 
           cost_RNA_discum  = cumsum(cost_RNA_disycum), 
           cost_POCT_discum  = cumsum(cost_POCT_disycum), 
           cost_TreatDAA_discum = cumsum(cost_TreatDAA_disycum),
           cost_RetreatDAA_discum = cumsum(cost_RetreatDAA_disycum),
           cost_TreatOther_discum = cumsum(cost_TreatOther_disycum),
           cost_RetreatOther_discum = cumsum(cost_RetreatOther_disycum),
           cost_totalDAA_discum = cumsum(cost_totalDAA_disycum),
           cost_Cured_discum = cumsum(cost_RetreatDAA_disycum)
    )%>%
    select(year, timestep, discount, 
           
           cost_total, cost_total_ycum, cost_total_cum, cost_total_disycum, 
           cost_total_discum, 
           
           cost_compartment, cost_compartment_ycum, cost_compartment_cum, 
           cost_compartment_disycum, cost_compartment_discum,
           
           cost_ab, cost_ab_ycum, cost_ab_cum, cost_ab_disycum, cost_ab_discum, 
           
           cost_RNA, cost_RNA_ycum, cost_RNA_cum, cost_RNA_disycum, 
           cost_RNA_discum, 
           
           cost_POCT, cost_POCT_ycum, cost_POCT_cum, cost_POCT_disycum, 
           cost_POCT_discum, 
           
           cost_TreatDAA, cost_TreatDAA_ycum, cost_TreatDAA_cum, 
           cost_TreatDAA_disycum, cost_TreatDAA_discum, 
           
           cost_RetreatDAA, cost_RetreatDAA_ycum, cost_RetreatDAA_cum, 
           cost_RetreatDAA_disycum, cost_RetreatDAA_discum, 
           
           cost_TreatOther, cost_TreatOther_ycum, cost_TreatOther_cum, 
           cost_TreatOther_disycum, cost_TreatOther_discum, 
           
           cost_RetreatOther, cost_RetreatOther_ycum, cost_RetreatOther_cum,
           cost_RetreatOther_disycum, cost_RetreatOther_discum, 
           
           cost_totalDAA, cost_totalDAA_ycum, cost_totalDAA_cum, 
           cost_totalDAA_disycum, cost_totalDAA_discum, 
           
           cost_Cured, cost_Cured_ycum, cost_Cured_cum, cost_Cured_disycum, 
           cost_Cured_discum
    )

}

View(Rescost_yearcum_all$dfList_NP_2024)

# save rda files 
save(Resflow_year_pop, Resflow_year_all, 
     Rescost_yearcum_pop, Rescost_yearcum_all,
     Rescost_DAA_dt,
     AUdiscount,
     unitC_DAA,
     cap,
     file = file.path(OutputFolder,
                      paste0(project_name,"Res_flowcost" ,".rda"))) 

# further generate the the csv/xlsx files to check numbers in each time step 
# (1) epi indicators by pops or all 
# (2) epi_cost combination 

for(i in names(Num_box)){ 
  Num_box[[i]] <- Num_box[[i]]%>%select(year, timestep, population,state, 
                                        cascade, disease_prog, best)%>%
    mutate(year = year + POC_AU$cabY - 1)
  
  Resflow_year_pop[[i]] <- Resflow_year_pop[[i]]%>%
    mutate(year = year + POC_AU$cabY)
  
  Resflow_year_all[[i]] <- Resflow_year_all[[i]]%>%
    mutate(year = year + POC_AU$cabY - 1 )
}

write.xlsx(Num_box, file = file.path(OutputFolder, paste0("POC_AU_epi_compartment.xlsx")), 
           append=TRUE) 
write.xlsx(Resflow_year_pop, file = file.path(OutputFolder, paste0("POC_AU_epi_flow_pop.xlsx")), 
           append=TRUE) 
write.xlsx(Resflow_year_all, file = file.path(OutputFolder, paste0("POC_AU_epi_flow_all.xlsx")), 
           append=TRUE) 

epicost_pop <- list()


epi_index <- c("Treatment", "Treatment_ycum", "Retreat", "Retreat_ycum", 
               "Cured", "Cured_ycum", "Testing_ab", "Testing_ab_ycum", 
               "Testing_RNA", "Testing_RNA_ycum", "Testing_POCT", "Testing_POCT_ycum")

for(i in names(Rescost_yearcum_pop)){ 
  Rescost_yearcum_pop[[i]] <- Rescost_yearcum_pop[[i]]%>%arrange(timestep)
  Resflow_year_pop[[i]] <- Resflow_year_pop[[i]]%>%arrange(timestep)
  epicost_pop[[i]] <- cbind(Rescost_yearcum_pop[[i]],
                            Resflow_year_pop[[i]][, (names(Resflow_year_pop[[i]]) %in% epi_index)])%>%
    select(year, timestep, population, discount, 
           cost_total, cost_total_ycum, cost_total_cum, cost_total_disycum, cost_total_discum,
           cost_compartment, cost_compartment_ycum, cost_compartment_cum,   
           cost_compartment_disycum, cost_compartment_discum, 
           Testing_ab, Testing_ab_ycum,
           cost_ab, cost_ab_ycum, cost_ab_cum, cost_ab_disycum, cost_ab_discum, 
           Testing_RNA, Testing_RNA_ycum,
           cost_RNA, cost_RNA_ycum, cost_RNA_cum, cost_RNA_disycum, cost_RNA_discum, 
           Testing_POCT, Testing_POCT_ycum,
           cost_POCT, cost_POCT_ycum, cost_POCT_cum, cost_POCT_disycum, cost_POCT_discum,
           Treatment, Treatment_ycum,
           cost_TreatDAA, cost_TreatDAA_ycum,  cost_TreatDAA_cum,
           cost_TreatDAA_disycum,cost_TreatDAA_discum, 
           Retreat, Retreat_ycum,
           cost_RetreatDAA, cost_RetreatDAA_ycum,  cost_RetreatDAA_cum,
           cost_RetreatDAA_disycum,cost_RetreatDAA_discum,
           cost_TreatOther, cost_TreatOther_ycum,  cost_TreatOther_cum,
           cost_TreatOther_disycum,cost_TreatOther_discum, 
           cost_RetreatOther, cost_RetreatOther_ycum,  cost_RetreatOther_cum,
           cost_RetreatOther_disycum,cost_RetreatOther_discum,
           cost_totalDAA, cost_totalDAA_ycum,  cost_totalDAA_cum,
           cost_totalDAA_disycum,cost_totalDAA_discum, 
           Cured, Cured_ycum,
           cost_Cured, cost_Cured_ycum, cost_Cured_cum,
           cost_Cured_disycum,cost_Cured_discum
           )
  
}
write.xlsx(epicost_pop, file = file.path(OutputFolder, paste0("POC_AU_epicost_pop.xlsx")), 
           append=TRUE) 

write.xlsx(Rescost_yearcum_pop, file = file.path(OutputFolder, paste0("POC_AU_cost_pop.xlsx")), 
           append=TRUE)

write.xlsx(Rescost_yearcum_all, file = file.path(OutputFolder, paste0("POC_AU_cost_all.xlsx")), 
           append=TRUE)

