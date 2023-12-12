# figures & table generation 
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

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "coverage_sq.rda")))
load(file.path(OutputFolder, paste0(project_name, "coverage_np.rda")))
load(file.path(OutputFolder, paste0(project_name, "coverage_npscale.rda")))
load(file.path(OutputFolder, paste0(project_name, "ResSum.rda")))
load(file.path(OutputFolder, paste0(project_name, "cost_np.rda")))
load(file.path(OutputFolder, paste0(project_name, "S_NPscale_test.rda")))
load(file.path(OutputFolder, paste0(project_name, "S_NP_test.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 

# return of investment 
ROI_np <- cost_np$`National program`%>%
  mutate(ROI = (cost_np$`status quo`$discountValue_cum - cost_np$`National program`$discountValue_cum)/cost_np$`National program`$discountValue_cum)%>%
  mutate(scenario = paste0("National program: community ", 
                           round(Ccal$C*100, digits = 1),"%" ,", prisons", round(Ccal$P*100, digits = 1), "%"))%>%
  select(year, scenario, discountValue, discountValue_cum, ROI)%>%
  filter(year <= 2050)

ROI_npscale <- cost_np$`National program scale up`%>%
  mutate(ROI = (cost_np$`status quo`$discountValue_cum - cost_np$`National program scale up`$discountValue_cum)/cost_np$`National program scale up`$discountValue_cum)%>%
  mutate(scenario = paste0("National program scale up: community ", 
                           round(Ccal_scale$C*100, digits = 1),"%" ,", prisons", round(Ccal_scale$P*100, digits = 1), "%"))%>%
  select(year, scenario, discountValue, discountValue_cum, ROI)%>%
  filter(year <= 2050)



## this is using the fixed budget to explore the return of investment 
# ROI_np <- cost_np$`National program`%>%
#  mutate(ROI = (cost_np$`status quo`$discountValue_cum - cost_np$`National program`$discountValue_cum)/6500000)%>%
#  mutate(scenario = paste0("National program: community ", 
#                           round(Ccal$C*100, digits = 1),"%" ,", prisons", round(Ccal$P*100, digits = 1), "%"))%>%
#  select(year, scenario, discountValue, discountValue_cum, ROI)%>%
#  filter(year <= 2050)

# ROI_npscale <- cost_np$`National program scale up`%>%
#  mutate(ROI = (cost_np$`status quo`$discountValue_cum - cost_np$`National program scale up`$discountValue_cum)/15000000)%>%
#  mutate(scenario = paste0("National program scale up: community ", 
#                           round(Ccal_scale$C*100, digits = 1),"%" ,", prisons", round(Ccal_scale$P*100, digits = 1), "%"))%>%
#  select(year, scenario, discountValue, discountValue_cum, ROI)%>%
#  filter(year <= 2050)
 ROI <- rbind(ROI_np, ROI_npscale)

ggplot(data = ROI_npscale) + 
  geom_line(aes(x = year, y = ROI))

# save csv

# number of infection & number of infection averted 
Num_avert_sub <- list()
Num_avert_all <- list()
Num_avert_setting <- list()
for(n in names(Sce_flow_sub$`status quo`)){ 
  Num_avert_sub[[n]] <- rbind(Sce_flow_sub$`status quo`[[n]],
                     Sce_flow_sub$`National program`[[n]],
                     Sce_flow_sub$`National program scale up`[[n]])%>%
    mutate(num_avert = Sce_flow_sub$`status quo`[[n]]$best - best)
  
  Num_avert_all[[n]] <- rbind(Sce_flow_all$`status quo`[[n]],
                              Sce_flow_all$`National program`[[n]],
                              Sce_flow_all$`National program scale up`[[n]])%>%
    mutate(num_avert = Sce_flow_all$`status quo`[[n]]$best - best)
  for(i in names(Sce_flow_setting)){ 
    Num_avert_setting[[i]][[n]] <- 
      rbind(Sce_flow_setting[[i]]$`status quo`[[n]],
            Sce_flow_setting[[i]]$`National program`[[n]],
            Sce_flow_setting[[i]]$`National program scale up`[[n]])%>%
      mutate(num_avert = Sce_flow_setting[[i]]$`status quo`[[n]]$best - best)
    
    }
  }

Num_avert_sub_bind <- dplyr::bind_rows(Num_avert_sub, .id = 'indicator')%>%
  mutate(year = year + POC_AU$cabY - 1)
Num_avert_all_bind <- dplyr::bind_rows(Num_avert_all, .id = 'indicator')%>%
  mutate(year = year + POC_AU$cabY - 1)

Num_avert_setting_bind <- list()
for(i in names(Sce_flow_setting)){ 
  Num_avert_setting_bind[[i]] <- 
    dplyr::bind_rows(Num_avert_setting[[i]], .id = 'indicator')%>%
    mutate(year = year + POC_AU$cabY - 1)
  
}
Num_state_sub <- list("Num_diag" = Num_diag,
                      "Num_diag" = Num_diag_ab, 
                      "Num_Treated" = Num_Treated, 
                      "Num_chronic_cured" = Num_chronic_cured, 
                      "Num_curInf" = Num_curInf, 
                      "Num_dc" = Num_dc, 
                      "Num_hcc" = Num_hcc, 
                      "Num_lt" = Num_lt)
  
Num_state_all <- list("Num_diag" = Num_diag_all,
                   "Num_diag" = Num_diag_ab_all, 
                   "Num_Treated" = Num_Treated_all, 
                   "Num_chronic_cured" = Num_chronic_cured_all, 
                   "Num_curInf" = Num_curInf_all, 
                   "Num_dc" = Num_dc_all, 
                   "Num_hcc" = Num_hcc_all, 
                   "Num_lt" = Num_lt_all)

Num_state_setting <- list("Num_diag" = Num_diag_setting,
                      "Num_diag" = Num_diag_ab_setting, 
                      "Num_Treated" = Num_Treated_setting, 
                      "Num_chronic_cured" = Num_chronic_cured_setting, 
                      "Num_curInf" = Num_curInf_setting, 
                      "Num_dc" = Num_dc_setting, 
                      "Num_hcc" = Num_hcc_setting, 
                      "Num_lt" = Num_lt_setting)





Num_state_avert_sub <- list()
Num_state_avert_all <- list()
Num_state_avert_setting <- list()
Num_state_avert_sub$Num_diag$`status quo`[97, ]
Num_state_avert_sub$Num_diag$`National program`[97, ]

for(n in names(Num_state_sub)){ 
  for(i in names(Num_state_sub[[1]])){ 
    Num_state_avert_sub[[n]][[i]] <- Num_state_sub[[n]][[i]]%>%as_tibble()%>%
      mutate(num_avert = Num_state_sub[[n]][["status quo"]]$best - Num_state_sub[[n]][[i]]$best)%>%
      mutate(scenario = i)
    
    
    Num_state_avert_all[[n]][[i]] <- Num_state_all[[n]][[i]]%>%as_tibble()%>%
      mutate(num_avert = Num_state_all[[n]][["status quo"]]$best - Num_state_all[[n]][[i]]$best )%>%
      mutate(scenario =  i)
  }
  
    for(s in names(Num_state_setting[[1]])){
      Num_state_avert_setting[[s]][[n]][[i]] <- Num_state_setting[[n]][[s]][[i]]%>%as_tibble()%>%
        mutate(num_avert = Num_state_setting[[n]][[s]][["status quo"]]$best - 
                 Num_state_setting[[n]][[s]][[i]]$best)%>%
        mutate(scenario = i)
    }
    
  
} 

previnc_subpop <- list("tempPrev" = tempPrev_subpop,
                       "tempPrevRNA" = tempPrevRNA_subpop,
                       "HCVInc" = HCVInc_subpop)

previnc_setting <- list("tempPrev" = tempPrev_setting,
                        "tempPrevRNA" = tempPrevRNA_setting,
                        "HCVInc" = HCVInc_setting)
previnc_subpop_redp <- list()
previnc_setting_redp <- list()
for(n in names(previnc_subpop)){ 
  for(i in names(previnc_subpop[[1]])){ 
    for(s in names(tempPrev_setting$`status quo`)){ 
      
      previnc_setting_redp[[n]][[i]][[s]] <- 
        previnc_setting[[n]][[i]][[s]]%>%as_tibble()%>%
        mutate(propo_reduce = 
                 100*(previnc_setting[[n]][["status quo"]][[s]]$best - 
                        previnc_setting[[n]][[i]][[s]]$best)/previnc_setting[[n]][["status quo"]][[s]]$best
        )%>%
        mutate( year = POC_AU$cabY + year - 1)%>%as_tibble()
      }
    
    previnc_subpop_redp[[n]][[i]] <- previnc_subpop[[n]][[i]]%>%as_tibble()%>%
      mutate(propo_reduce = 
               100*(previnc_subpop[[n]][["status quo"]]$best - 
                  previnc_subpop[[n]][[i]]$best)/previnc_subpop[[n]][["status quo"]]$best
      )%>%
      mutate( year = POC_AU$cabY + year - 1)%>%as_tibble()
    }
} 



previnc_subpop_bind <- list()
previnc_setting_bind <- list()
for(i in names(previnc_subpop_redp)){ 
  previnc_subpop_bind[[i]] <- dplyr::bind_rows(previnc_subpop_redp[[i]], .id = 'scenario')
  for(i in names(previnc_setting_redp)){ 
    previnc_subpop_bind[[i]] <- dplyr::bind_rows(previnc_subpop_redp[[i]], .id = 'scenario')
    }
  }



save(ROI,cost_stage,
  Num_avert_sub,
     Num_avert_all,
     Num_avert_setting, 
     Num_state_avert_sub, Num_state_avert_all, Num_state_avert_setting,
     previnc_setting_redp,previnc_setting_redp,
     file = file.path(OutputFolder ,
                      paste0(project_name,"TabGen" ,".rda")))


Num_avert_sub_bind%>%filter(year == 2023, indicator == "newTestingAg")


