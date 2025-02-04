#This script using the Rda file {Res_flowcost.rda} from {Res_aggregate.R} to 
# generate the Tables and Figures 
# Tables and figures are generated based on overall pops, all subpops, and settings(community and prisons)
# the indicators include epi indicators and cost indicators 
# epi indicators include 
# 1. number of infections, number of infections averted 
# 2. number of advanced liver disease and averted 
#    (a) number of dc, (b) number of hcc, (c) number of lt and its averted cases 
# 3. number of HCV deaths and its averted cases 
# 4. Number of people living with HCV: excluding "s" and "cured" 
# 5. Number of people living with HCV get tested
# 6. Number of people living with diagnosed HCV
# 7. Number of reinfection 
# cost 
gc()
rm(list = ls())
gc()
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
Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")
OutputFig <- file.path(OutputFolder, "Figs")

OutputFig_y_cum_avert <- file.path(OutputFig , "y_cum_avert")
load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "Res_flowcost.rda")))

load(file.path(OutputFolder, paste0(project_name, "Res_numbox.rda")))

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
endY <- 100
year_obs <- c(POC_AU$simY  +5 - 1  , POC_AU$simY + 10 - 1, POC_AU$simY + 20 - 1)
par_col <- c("best", paste0("set", seq(1,1000,1)))
# sum up the total number of cascade 
for(i in names(Resflow_year_all)){ 
  for(n in names(Resflow_year_all[[1]])){ 
    Resflow_year_all[[i]][[n]] <- Resflow_year_all[[i]][[n]]%>%ungroup()
    Resflow_year_all[[i]][[n]][is.na(Resflow_year_all[[i]][[n]])] <- 0
    Resflow_year_pop[[i]][[n]] <-  Resflow_year_pop[[i]][[n]]%>%ungroup()
    Resflow_year_pop[[i]][[n]][is.na(Resflow_year_pop[[i]][[n]])]  <- 0 
      
    }
  }

for(i in names(Resflow_year_all)){ 
  
  Resflow_year_all[[i]][["Tot_Treatment"]] <- 
    cbind(year = Resflow_year_all[[i]][["Treatment"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Treatment"]][, par_col] + 
                          Resflow_year_all[[i]][["Retreat"]][, par_col] + 
                          Resflow_year_all[[i]][["Treatment_sc"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_ab"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_ab"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_ab"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_ab_neg"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_ab_sc"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_ab_sc_neg"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_RNA"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_RNA"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_RNA"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_RNA_neg"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_RNA_sc"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_POCT"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          dplyr::bind_cols(Resflow_year_all[[i]][["Testing_POCT"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_POCT_neg"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_POCT_sc"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )
}

for(i in names(Resflow_year_all)){
  Resflow_year_pop[[i]][["Tot_Treatment"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Treatment"]]$year, 
          population = Resflow_year_pop[[i]][["Treatment"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Treatment"]][, par_col] + 
                          Resflow_year_pop[[i]][["Retreat"]][, par_col] + 
                          Resflow_year_pop[[i]][["Treatment_sc"]][, par_col])
    )
  
  Resflow_year_pop[[i]][["Tot_Testing_ab"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_ab"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_ab"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_ab"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_ab_neg"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_ab_sc"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_ab_sc_neg"]][, par_col])
    )

    
  Resflow_year_pop[[i]][["Tot_Testing_RNA"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_RNA"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_RNA"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_RNA"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_RNA_neg"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_RNA_sc"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )
  
  Resflow_year_pop[[i]][["Tot_Testing_POCT"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          dplyr::bind_cols(Resflow_year_pop[[i]][["Testing_POCT"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_POCT_neg"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_POCT_sc"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )
  
}
Resflow_year_setting <- list()

for(i in names(Resflow_year_pop)){ 
  for(indic in names(Resflow_year_pop[[1]])){ 
    Resflow_year_setting[[i]][[indic]] <- 
      Resflow_year_pop[[i]][[indic]]%>%
      mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), 
                              "commu", "prisons"))%>%
      group_by(year, setting)%>%
      summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
      ungroup()%>%
      mutate(population = setting)%>%
      select(-setting)%>%
      mutate(population = factor(population, 
                    levels = c("commu", "prisons"), 
                    labels = c("Community", "Prisons")))%>%
      select(year, population, par_col)
      
     
    
    }
}


Resflow_year_setting_range <- list()
Resflow_year_pop_range <- list()
Resflow_year_all_range <- list()
for(i in names(Resflow_year_pop)){ 
  for(indic in names(Resflow_year_pop[[1]])){ 
    Resflow_year_setting_range[[i]][[indic]] <- 
      Resflow_year_setting[[i]][[indic]]%>%ungroup()%>%
      popResults_range(POC_AU, ., Population = c("Community", "Prisons"),
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)

    Resflow_year_pop_range[[i]][[indic]] <- 
      Resflow_year_pop[[i]][[indic]]%>%ungroup()%>%
      popResults_range(POC_AU, ., Population = POC_AU$popNames,
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
    
    Resflow_year_all_range[[i]][[indic]] <- 
      Resflow_year_all[[i]][[indic]]%>%ungroup()%>%
      popResults_range(POC_AU, ., Population = NULL,
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
    }
  }

#### cumulative: number of infections: overall in 5, 10, 20y ####
# starting at year of 2022
# tidy up data in dt to yearly and cumulative 
Resflow_cum_all <- list()
par_col <- c("best", paste0("set", seq(1, POC_AU$numberSamples,1)))

for(i in names(Resflow_year_all)){ 
  for(indic in names(Resflow_year_all[[1]])){ 
    Resflow_cum_all[[i]][[indic]] <- Resflow_year_all[[i]][[indic]]%>%
      ungroup()%>%
      mutate(year = year + POC_AU$cabY)%>%
      filter(year>= POC_AU$simY )%>%
      ungroup()%>%
      mutate(across(par_col, list(cum=cumsum), .names = "{col}"))
    }
}


# avert number of transitions 
Resflow_cum_all_avert <- list()

for(i in names(Resflow_cum_all)){ 
  for(indic in names(Resflow_cum_all[[1]])){
    Resflow_cum_all_avert[[i]][[indic]] <- cbind(
      year = Resflow_cum_all[[i]][[indic]]$year, 
      data.frame(Resflow_cum_all[["Status quo"]][[indic]][, c(par_col)] - 
                   Resflow_cum_all[[i]][[indic]][, c(par_col)]))
  }
}

Resflow_sc_cum_all <- list()

for(i in names(Resflow_sc_year_all)){ 
  for(indic in names(Resflow_sc_year_all[[1]])){
    Resflow_sc_cum_all[[i]][[indic]] <- Resflow_sc_year_all[[i]][[indic]]%>%
      ungroup()%>%
      mutate(year = year + POC_AU$cabY)%>%
      filter(year>= POC_AU$simY - 1 )%>%
      mutate(across(par_col, list(cum=cumsum), .names = "{col}"))
  }
}
# avert number of transitions 
Resflow_sc_cum_all_avert <- list()
for(i in names(Resflow_sc_cum_all)){ 
  for(indic in names(Resflow_sc_cum_all[[1]])){
    Resflow_sc_cum_all_avert[[i]][[indic]] <- cbind(
      year = Resflow_sc_cum_all[[i]][[indic]]$year, 
      dplyr::bind_cols(Resflow_sc_cum_all[["Status quo"]][[indic]][, c(par_col)] - 
                   Resflow_sc_cum_all[[i]][[indic]][, c(par_col)]))
  }
}



# aggregate q2.5 - q97.5
# Resflow_year_all, Resflow_sc_year_all
# Resflow_cum_all, Resflow_cum_all_avert, 
# Resflow_sc_cum_all, Resflow_sc_cum_all_avert 



Resflow_year_all_range <- list() 

for(i in names(Resflow_year_all)){ 
  for(indic in names(Resflow_year_all[[1]])){
    Resflow_year_all_range[[i]][[indic]] <- Resflow_year_all[[i]][[indic]]%>%
      mutate(year = year + POC_AU$cabY)%>%
      filter(year>= POC_AU$simY - 1)%>%
      popResults_range(POC_AU, . , Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}

Resflow_sc_year_all_range <- list() 

for(i in names(Resflow_sc_year_all)){ 
  for(indic in names(Resflow_sc_year_all[[1]])){
    Resflow_sc_year_all_range[[i]][[indic]] <- Resflow_sc_year_all[[i]][[indic]]%>%
      mutate(year = year + POC_AU$cabY)%>%
      filter(year>= POC_AU$simY - 1 )%>%
      popResults_range(POC_AU, . , Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
      
  }
}

Resflow_cum_all_range <- list() 

for(i in names(Resflow_cum_all)){ 
  for(indic in names(Resflow_cum_all[[1]])){
    Resflow_cum_all_range[[i]][[indic]] <- 
      popResults_range(POC_AU, Resflow_cum_all[[i]][[indic]], Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
    }
  }

Resflow_cum_all_avert_range <- list() 

for(i in names(Resflow_cum_all_avert)){ 
  for(indic in names(Resflow_cum_all_avert[[1]])){
    Resflow_cum_all_avert_range[[i]][[indic]] <- 
      popResults_range(POC_AU, Resflow_cum_all_avert[[i]][[indic]], Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}  

Resflow_sc_cum_all_range <- list()
for(i in names(Resflow_sc_cum_all)){ 
  for(indic in names(Resflow_sc_cum_all[[1]])){
    Resflow_sc_cum_all_range[[i]][[indic]] <- 
      popResults_range(POC_AU, Resflow_sc_cum_all[[i]][[indic]], Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}  

Resflow_sc_cum_all_avert_range <- list()
for(i in names(Resflow_sc_cum_all)){ 
  for(indic in names(Resflow_sc_cum_all[[1]])){
    Resflow_sc_cum_all_avert_range[[i]][[indic]] <- 
      popResults_range(POC_AU, Resflow_sc_cum_all_avert[[i]][[indic]], Population = NULL, 
                       Disease_prog = NULL, Cascade = NULL, end_Y = endY - 1)
  }
}  

# turn list inside out and tidy dt for plots and table generation 
# combine the cascade numbers 

Resflow_all <- list("Resflow_year" = Resflow_year_all_range, 
                    "Resflow_cum" = Resflow_cum_all_range, 
                    "Resflow_cum_avert" = Resflow_cum_all_avert_range)

Resflow_all_lst <- list()
for(i in names(Resflow_all)){ 
  Resflow_all_lst[[i]] <- Resflow_all[[i]]%>%purrr::transpose()%>%
    lapply(., function(x) dplyr::bind_rows(x, .id = 'scenario'))
  
  for(n in names(Resflow_all_lst[[i]])){ 
    Resflow_all_lst[[i]][[n]] <- Resflow_all_lst[[i]][[n]]%>%
      mutate(scenario = factor(scenario, 
                               levels = c("Status quo", "dfList_NP_2023", 
                                          "dfList_NP_2024", "dfList_NPexp_A", 
                                          "dfList_NPexp_B", "dfList_NPexp_C",
                                          "dfList_NPexp_D"), 
                               labels = c("No National Program", "Past program", 
                                          "Current program", "Continued program", 
                                          "NP expand 2025", "NP expand 2026", 
                                          "Expanded program")))
    }
  }

# extracting the scenario wanna present 
Resflow_all_lst_subsce <- list()
for(i in names(Resflow_all_lst)){
  for(n in names(Resflow_all_lst[[1]])){ 
    Resflow_all_lst_subsce[[i]][[n]] <- Resflow_all_lst[[i]][[n]]%>%
      filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
    
    }
}



#### plot for number of averted/cumulative infections by 5, 10, 15 years ####
# ceiling and floor function 

# function for generating plots 
plot_pocau <- function(pj,dt,indicator, type, year_obs, UI = NULL){ 
  # type: new, cum, and avert  
  if(length(unique(dt[[indicator]]$scenario)) == 2){ 
    col_pal <- c("#000000", "#E69F00")
    
  } 
  else{col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  }
  
  if(isTRUE(type == "new")){ 
    
    a <- dt[[indicator]]%>%ggplot(data = ., 
           aes(x = year, colour = scenario)) + 
      geom_line(aes(y = best, colour = scenario, linetype = scenario), size = 1) + 
      scale_x_continuous(expand = c(0,0), limits = c(pj$simY - 1,2050), breaks = c(seq(pj$simY - 1,2050, 5), 2050))+
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt[[indicator]]$scenario)) - 1))) +
      theme_Publication()
    
  }
  else if(isTRUE(type =="cum")){ 
   a <- dt[[indicator]]%>%filter(year %in% year_obs)%>%
     ggplot(data =. , aes(x = as.character(year),  colour = scenario)) + 
     geom_errorbar(aes(ymin = q5, ymax = q95, colour = scenario), width = 0.5,
                   size = 1, 
                   position = "dodge") + 
     geom_point(aes(y = best, colour = scenario),position = position_dodge(width = 0.5),size = 1.5) +
     labs(x = "", y = indicator) + 
     scale_x_discrete(labels = c(paste0((year_obs - pj$simY + 1), "-Year", sep = ""))) + 
     scale_color_manual(name = "Scenarios", values = col_pal ) + 
     scale_fill_manual(name = "Scenarios", values = col_pal ) + 
     scale_linetype_manual(name = "Scenarios", 
                           values = c("dashed", rep("solid", length(unique(dt[[indicator]]$scenario)) - 1))) +
     theme_Publication()
    
    
  }
  
  else if(isTRUE(type =="avert")){ 
    a <- dt[[indicator]]%>%filter( year %in% year_obs)%>%
      ggplot(data =. , aes(x = as.character(year), colour = scenario)) + 
      geom_errorbar(aes(ymin = q5, ymax = q95, colour = scenario),
                    size = 1, position = "dodge" , width = 0.5) + 
      geom_point(aes(y = best, colour = scenario), position = position_dodge(width = 0.5) , size = 1.5) + 
      labs(x = "", y = indicator) + 
      scale_x_discrete(labels = c(paste0((year_obs - pj$simY + 1), "-Year", sep = ""))) + 
      scale_color_manual(name = "Scenarios", values = col_pal[-1] ) + 
      scale_fill_manual(name = "Scenarios", values = col_pal[-1] ) + 

      theme_Publication()
    
  } 
  a <- a + theme(legend.key.size = unit(1,"line"))
  
  return(a)
}
lim_ident <- function(dt, year_range){ 
  
  lim <- dt%>%as.data.frame()%>%ungroup()%>%
    filter(year %in% year_range)%>%
    summarise(x = max(q95, na.rm = TRUE))%>%
    mutate(lim = case_when( 
      x <1 ~5, 
      x>=1 & x <10 ~ 10, 
      x>=10 & x< 30 ~ 30, 
      x>=30 & x<60 ~ 60, 
      x>=60 & x<80 ~ 80, 
      x>80 & x<=100 ~100,
      x>100 & x<=1000 ~ (x%/%100 + 1)*100,
      x >1000 & x<=10000 ~ (x%/%1000 + 1)*1000,
      x >10000 & x <= 100000  ~ (x%/%10000 + 1)*10000,
      x >100000 & x <= 1000000  ~ (x%/%100000 + 1)*100000,
      x >1000000 & x <= 10000000 ~ (x%/%1000000 + 1)*1000000
    ))
  
  return(lim)  
}


for(i in names(Resflow_all_lst_subsce$Resflow_cum_avert)){ 
  Resflow_all_lst_subsce$Resflow_cum_avert[[i]] <- 
    Resflow_all_lst_subsce$Resflow_cum_avert[[i]]%>%
    filter(scenario != "No National Program")
  }


p_pocau_y <- list()
p_pocau_cum <- list()
p_pocau_avert <- list()
lim_dt_y <- list()
lim_dt_cum <- list()
lim_dt_avert <- list()
ylab_name <- list()
ylab_name <- list("New HCV Infections", 
                  "Number of HCV releated deaths",
                  "Number of Treatment",
                  "Number of Retreat",
                  "Number of antibody testing\n(HCV positive, out of National Program)",
                  "Number of RNA testing\n(two-step, HCV positive, out of National Program)",
                  "Number of point-of-care RNA testing\n(HCV positive, out of National Program)",
                  "Number of antibody testing\n(HCV negative, out of National Program)",
                  "Number of RNA testing\n(two-step, HCV negative, out of National Program)",
                  "Number of point-of-care RNA testing\n(HCV negative, out of National Program)",
                  "Number of SVR achieved",
                  "Number of treatment initiated\nvia National Program",
                  "Number of antibody testing via National Program\n(HCV positive)",
                  "Number of RNA testing via\nNational Program(HCV positive)",
                  "Number of point-of-care RNA testing\nvia National Program(HCV positive)",
                  "Number of antibody testing via\nNational Program(HCV negative)",
                  "Number of RNA testing via\nNational Program(HCV negative)",
                  "Number of point-of-care RNA testing via\nNational Program(HCV negative)",
                  "Number of total treatment initiated",
                  "Number of total antibody testing",
                  "Number of total two-step RNA testing",
                  "Number of total point-of-care RNA testing")

yavert_lab_name <- paste0("Averted ", tolower(ylab_name))
names(ylab_name) <- names(Resflow_all_lst_subsce[[1]])
names(yavert_lab_name) <- names(Resflow_all_lst_subsce[[1]])


for(indic in names(Resflow_all_lst_subsce[[1]])){
  
  lim_dt_y[[indic]] <- lim_ident(Resflow_all_lst_subsce$Resflow_year[[indic]], seq(2021, 2050, 1))
  lim_dt_cum[[indic]] <- lim_ident(Resflow_all_lst_subsce$Resflow_cum[[indic]], seq(2022, year_obs[3], 1))
  lim_dt_avert[[indic]] <- lim_ident(Resflow_all_lst_subsce$Resflow_cum_avert[[indic]], seq(2022, year_obs[3], 1))
  
  p_pocau_y[[indic]] <- 
    plot_pocau(POC_AU, Resflow_all_lst_subsce$Resflow_year, indicator = indic,
               type = "new", 5) + 
    labs( y = ylab_name[[indic]], x = "Year") + 
    theme(axis.title = element_text()) + 
    scale_y_continuous(limits = c(0, as.numeric(lim_dt_y[[indic]][, "lim"])),
                       breaks = seq(0, as.numeric(lim_dt_y[[indic]][, "lim"]),
                                    (as.numeric(lim_dt_y[[indic]][, "lim"] - 0))/10)) + 
    theme(legend.position = "right", legend.direction="vertical")
  p_pocau_cum[[indic]] <- 
    plot_pocau(POC_AU, Resflow_all_lst_subsce$Resflow_cum, type = "cum", 
               indicator = indic, year_obs = year_obs) + 
    ylab(ylab_name[[indic]]) + 
    theme(axis.title = element_text()) + 
    scale_y_continuous(limits = c(0, as.numeric(lim_dt_cum[[indic]][, "lim"])),
                       breaks = seq(0, as.numeric(lim_dt_cum[[indic]][, "lim"]),
                                    (as.numeric(lim_dt_cum[[indic]][, "lim"] - 0))/10)) + 
    theme(legend.position = "right", legend.direction="vertical")
  p_pocau_avert[[indic]] <- 
    plot_pocau(POC_AU, 
               Resflow_all_lst_subsce$Resflow_cum_avert, type = "avert", 
               indicator = indic, year_obs = year_obs) + 
    ylab(yavert_lab_name[[indic]]) + 
    scale_y_continuous(limits = c(0, as.numeric(lim_dt_avert[[indic]][, "lim"])),
                       breaks = seq(0, as.numeric(lim_dt_avert[[indic]][, "lim"]),
                                    (as.numeric(lim_dt_avert[[indic]][, "lim"] - 0))/10)) + 
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 

    theme(axis.title = element_text()) + 
    theme(legend.position = "right", legend.direction="vertical")
    
}

View(Resflow_all_lst_subsce$Resflow_year$Tot_Testing_ab)
p_pocau_y$HCVdeath <- p_pocau_y$HCVdeath + 
  scale_x_continuous(expand = c(0,0),limits = c(2021, 2041),breaks = seq(2021, 2041, 5))
p_pocau_y$Tot_Treatment +scale_y_continuous(limits = c(0, 10000), 
                                            breaks = seq(0, 10000, 1000)) 
p_pocau_avert$Treatment <- p_pocau_avert$Treatment + 
  scale_y_continuous(limits = c(-10000, 30000), 
                     breaks = seq(-10000, 30000, 5000)) 

p_pocau_avert$Retreat <- p_pocau_avert$Retreat + 
  scale_y_continuous(limits = c(-3000, 2000), 
                     breaks = seq(-3000, 2000, 500))

p_pocau_avert$Testing_ab <- p_pocau_avert$Testing_ab + 
  scale_y_continuous(limits = c(0, 50000), 
                     breaks = seq(0, 50000, 5000))

p_pocau_avert$Testing_ab_neg <- p_pocau_avert$Testing_ab_neg + 
  scale_y_continuous(limits = c(-20000, 10000), 
                     breaks = seq(-20000, 10000, 1000))

p_pocau_avert$Testing_RNA <- p_pocau_avert$Testing_RNA + 
  scale_y_continuous(limits = c(-1000, 30000), 
                     breaks = seq(-1000, 30000, 1000))



p_pocau_avert$Testing_RNA_neg <- p_pocau_avert$Testing_RNA_neg + 
  scale_y_continuous(limits = c(-8000, 8000), 
                     breaks = seq(-8000, 8000, 1000))

p_pocau_avert$Testing_POCT <- p_pocau_avert$Testing_POCT + 
  scale_y_continuous(limits = c(-20000, 20000), 
                     breaks = seq(-20000, 20000, 1000))

p_pocau_avert$Testing_POCT_neg <- p_pocau_avert$Testing_POCT_neg + 
  scale_y_continuous(limits = c(-8000, 0), 
                     breaks = seq(-8000, 0, 1000))



p_pocau_avert$Cured <- p_pocau_avert$Cured + 
  scale_y_continuous(limits = c(-20000, 20000), 
                     breaks = seq(-20000, 20000, 5000))

p_pocau_avert$Treatment_sc <- p_pocau_avert$Treatment_sc + 
  scale_y_continuous(limits = c(-15000, 0), 
                     breaks = seq(-15000, 0, 3000))

p_pocau_avert$Testing_ab_sc <- p_pocau_avert$Testing_ab_sc + 
  scale_y_continuous(limits = c(-15000, 0), 
                     breaks = seq(-15000, 0, 1000))

p_pocau_avert$Testing_RNA_sc <- p_pocau_avert$Testing_RNA_sc + 
  scale_y_continuous(limits = c(-10000, 0), 
                     breaks = seq(-10000, 0, 1000))

p_pocau_avert$Testing_POCT_sc <- p_pocau_avert$Testing_POCT_sc + 
  scale_y_continuous(limits = c(-20000, 0), 
                     breaks = seq(-2000, 0, 1000))
p_pocau_avert$Testing_ab_sc_neg <- p_pocau_avert$Testing_ab_sc_neg + 
  scale_y_continuous(limits = c(-80000, 0), 
                     breaks = seq(-80000, 0, 5000))

p_pocau_avert$Testing_RNA_sc_neg <- p_pocau_avert$Testing_RNA_sc_neg + 
  scale_y_continuous(limits = c(-40000,100), 
                     breaks = seq(-40000, 100, 5000))

p_pocau_avert$Testing_POCT_sc_neg <- p_pocau_avert$Testing_POCT_sc_neg + 
  scale_y_continuous(limits = c(-300000,0), 
                     breaks = seq(-300000, 0, 10000))

p_pocau_avert$Tot_Treatment <- p_pocau_avert$Tot_Treatment + 
  scale_y_continuous(limits = c(-10000,20000), 
                     breaks = seq(-10000,20000, 2500))

p_pocau_avert$Tot_Testing_ab <- p_pocau_avert$Tot_Testing_ab + 
  scale_y_continuous(limits = c(-100000,10000), 
                     breaks = c(seq(-100000,10000, 5000), c(2000)))
p_pocau_avert$Tot_Testing_RNA <- p_pocau_avert$Tot_Testing_RNA + 
  scale_y_continuous(limits = c(-50000,20000), 
                     breaks = seq(-50000,20000, 10000))

p_pocau_avert$Tot_Testing_POCT <- p_pocau_avert$Tot_Testing_POCT + 
  scale_y_continuous(limits = c(-300000,0), 
                     breaks = seq(-300000,0, 10000))


for(i in names(p_pocau_avert)){ 
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0(i,"_avert" ,".png")), 
         p_pocau_avert[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0(i,"_cum" ,".png")), 
         p_pocau_cum[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0(i,"_y" ,".png")), 
         p_pocau_y[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
}

# number of treatment, testing plot with calibration points 
Resflow_year_setting_range <- Resflow_year_setting_range%>%
  purrr::transpose()%>%
  lapply(., function(x) dplyr::bind_rows(x, .id = 'scenario'))


Resflow_year_setting_range_trajectory <- lapply(Resflow_year_setting_range, function(x)
  x%>%mutate(scenario = factor(scenario, levels = c("Status quo", "dfList_NP_2023", "dfList_NP_2024", 
                                                             "dfList_NPexp_A", "dfList_NPexp_B", "dfList_NPexp_C",
                                                             "dfList_NPexp_D"), 
                               labels = c("No National Program", "Past program", 
                                          "Current program", "Continued program", 
                                          "NP expand 2025", "NP expand 2026", 
                                          "Expanded program")))%>%

    filter(! scenario %in% c("NP expand 2025", "NP expand 2026"))
)
# plot function for generating cascade numbers 
Cas_num_plot <- function(pj, dt, obdt =NULL, xlimits, UI = NULL, population = NULL){ 
  # @ UI the name of scenario for showing UI range
  if(length(unique(dt$scenario)) == 2){ 
    col_pal <- c("#000000", "#E69F00")
    
  } 
  else{col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  }
  
  if(is.null(obdt) & is.null(UI) & !is.null(population)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario), size = 1) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      theme(panel.spacing = unit(2, "lines")) + theme_Publication() + 
      theme(legend.key.size = unit(1,"line"))
  }
  
  else if(is.null(obdt) & is.null(UI)){ 
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
  }
  else if(is.null(obdt) & !is.null(UI)){ 
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      geom_ribbon(data =dt%>%filter(scenario == UI) ,aes(ymin = q5, ymax = q95), fill = "#000000", alpha = 0.2) +
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = "#000000" ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
    
  }
  
  else if(!is.null(obdt) & is.null(UI)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      geom_point(data=obdt, aes(y=realPop, x = time), 
                 colour = "black", size = 1) +
      geom_segment(data = obdt, 
                   aes ( y = low, yend = up, x = time, xend = time)) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
    
    
  }
  else if(!is.null(obdt) & !is.null(UI) & is.null(population)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      geom_ribbon(data =dt%>%filter(scenario == UI), aes(ymin = q5, ymax = q95), fill = "#000000", alpha = 0.2) +
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = "#000000" ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      geom_point(data=obdt, aes(y=realPop, x = time), 
                 colour = "black", size = 1) +
      geom_segment(data = obdt, 
                   aes ( y = low, yend = up, x = time, xend = time)) +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication_facet() + 
      theme(legend.key.size = unit(1,"line"))
  }else if(!is.null(obdt) & !is.null(UI) & !is.null(population)){
    traj_plot <- ggplot(dt, aes(x = year, y = best)) + 
      geom_line(aes(colour = scenario, linetype = scenario)) + 
      geom_ribbon(data =dt%>%filter(scenario == UI), aes(ymin = q5, ymax = q95), fill = "#000000", alpha = 0.2) +
      
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = "#000000" ) + 
      scale_linetype_manual(name = "Scenarios", 
                            values = c("dashed", rep("solid", length(unique(dt$scenario)) - 1))) + 
      coord_cartesian(xlim = xlimits[1:2]) +
      scale_x_continuous(expand = c(0.01, 0), limits =c(xlimits[1], xlimits[2]) ,
                         breaks = seq(xlimits[1], xlimits[2], 
                                      by = xlimits[3]),
                         labels = seq(pj$cabY + xlimits[1] - 1, 
                                      (pj$cabY + xlimits[2] - 1),
                                      xlimits[3])) + 
      geom_point(data=obdt, aes(y=realPop, x = time), 
                 colour = "black", size = 1)  +
      theme(panel.spacing = unit(2, "lines")) + theme_Publication() + 
      theme(legend.key.size = unit(1,"line"))
  }
  
  return(traj_plot)
}

HCVtreatinitN_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))
HCVtreatinitN_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower)

HCVtreatinitN_NP_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_NP_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower, 
                           population = factor(population, 
                                               levels = c("Total", "Natioanl Program"), 
                                               labels = c("Total", "National Program")))


Resflow_year_setting_range_trajectory$Tot_Treatment  <- 
  Resflow_year_setting_range_trajectory$Tot_Treatment%>%mutate(year = year + 1)

total_treatm <- Resflow_year_setting_range_trajectory$Tot_Treatment%>%
  group_by(year, scenario)%>%summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))
total_treatm_lst <- list()
for(i in unique(total_treatm$scenario)){ 
  
  total_treatm_lst[[i]] <- total_treatm%>%filter(scenario == i)%>%
    popResults_range(POC_AU, .)
  }

total_treatm_lst <- total_treatm_lst%>%bind_rows(., .id = 'scenario')%>%
  mutate(scenario = factor(scenario, levels = c("No National Program", "Past program", 
                                                "Current program", "Continued program", 
                                                 
                                                "Expanded program")))



p_T_num <- Cas_num_plot(POC_AU,Resflow_year_setting_range_trajectory$Tot_Treatment%>%
                          filter(scenario %in% c("No National Program", "Past program")) , 
             obdt =HCVtreatinitN_setting_fit, 
             xlimits = c(1, 16, 5), UI = "No National Program") + 
  labs(x = "Year", y = "Number of total treatment initiated") 
p_T_num <- p_T_num + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 40000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 8000)))
                  )) +
  theme(legend.position = "right", legend.direction="vertical")
Resflow_sc_year_all_range$dfList_NP_2023$Treatment_sc%>%select(year, best, q5, q95)
# total treatment number 
p_T_num_total <- Cas_num_plot(POC_AU,total_treatm_lst , 
                        obdt =HCVtreatinitN_fit, 
                        xlimits = c(1, 31, 5), UI = "No National Program", 
                        population = "n") + 
  labs(x = "Year", y = "Number of total treatment initiated") +
  theme(legend.position = "",
        legend.direction = "vertical")

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numTreat_all", ".png")), 
       p_T_num_total , 
       width = 8, height = 6, bg = "white", dpi = 300)



# new inf and treatment number for main text 

p_T_num_total_maintext <- Cas_num_plot(POC_AU,total_treatm_lst , 
                              obdt = NULL, 
                              xlimits = c(7, 16, 1), UI = NULL, 
                              population = "n") + 
  labs(x = "Year", y = "Treatment initiations", tag ="B") +
  theme(legend.position = "",
        legend.direction = "vertical") + 
  scale_y_continuous(limits = c(0, 9000), breaks = seq(0, 9000,500))

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("T_num_total_maintext", ".png")), 
       p_T_num_total_maintext, 
       width = 8, height = 6, bg = "white", dpi = 300)

p_newinf_num_total_maintext <- p_pocau_y$newInfections + 
  scale_x_continuous(expand = c(0,0), limits = c(2021, 2030), breaks = seq(2021, 2030, 1)) + 
  scale_y_continuous(limits = c(0, 9000), breaks = seq(0, 9000, 500)) + 
  labs(x = "Year", y = "New HCV infections", tag ="A") + 
  theme(legend.position = "",
        legend.direction = "vertical") 


ggsave(file=file.path(OutputFig_y_cum_avert, paste0("newinf_num_total_maintext", ".png")), 
       p_newinf_num_total_maintext , 
       width = 8, height = 6, bg = "white", dpi = 300)

p_newinfnT_maintext  <- ggarrange(plotlist = list(p_newinf_num_total_maintext, 
                                                  p_T_num_total_maintext), ncol = 2, nrow = 1, 
                       common.legend = FALSE)

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("newinfnT_maintext ", ".png")), 
       p_newinfnT_maintext , 
       width = 10, height = 6, bg = "white", dpi = 300)


p_HCVdeath_supp <- p_pocau_y$HCVdeath + 
  scale_x_continuous(expand = c(0,0), limits = c(2021, 2041), breaks = seq(2021, 2041, 1)) + 
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 50)) + 
  labs(x = "Year", y = "Number of HCV-related deaths") + 
  theme(legend.position = "bottom",
        legend.direction = "vertical") 


ggsave(file=file.path(OutputFig_y_cum_avert, paste0("HCVdeath_supp", ".png")), 
       p_HCVdeath_supp , 
       width = 8, height = 6, bg = "white", dpi = 300)
# testing 
xt_np_ab <- cbind(year = Resflow_year_setting_range_trajectory$Testing_ab_sc$year, 
                  population = Resflow_year_setting_range_trajectory$Testing_ab_sc$population, 
            scenario = Resflow_year_setting_range_trajectory$Testing_ab_sc$scenario, 
            as.data.frame(Resflow_year_setting_range_trajectory$Testing_ab_sc[, par_col] + 
                            Resflow_year_setting_range_trajectory$Testing_ab_sc_neg[, par_col]))%>%
  as.data.frame()%>%
  popResults_range(POC_AU, . , Population = c("Community", "Prisons"),
                   Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
  mutate(NP = "NP")

xt_pnp_ab <- cbind(year = Resflow_year_setting_range_trajectory$Testing_ab$year, 
                   population = Resflow_year_setting_range_trajectory$Testing_ab_sc$population, 
                  scenario = Resflow_year_setting_range_trajectory$Testing_ab$scenario, 
                  as.data.frame(Resflow_year_setting_range_trajectory$Testing_ab[, par_col] + 
                                  Resflow_year_setting_range_trajectory$Testing_ab_neg[, par_col]))%>%
  as.data.frame()%>%
  popResults_range(POC_AU, . , Population = c("Community", "Prisons"),
                   Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
  mutate(NP = "out of NP")

xt_np_rna <- cbind(year = Resflow_year_setting_range_trajectory$Testing_RNA_sc$year, 
                  population = Resflow_year_setting_range_trajectory$Testing_RNA_sc$population, 
                  scenario = Resflow_year_setting_range_trajectory$Testing_RNA_sc$scenario, 
                  as.data.frame(Resflow_year_setting_range_trajectory$Testing_RNA_sc[, par_col] + 
                                  Resflow_year_setting_range_trajectory$Testing_RNA_sc_neg[, par_col] + 
                                  Resflow_year_setting_range_trajectory$Testing_POCT_sc[, par_col] + 
                                  Resflow_year_setting_range_trajectory$Testing_POCT_sc_neg[, par_col]))%>%
  as.data.frame()%>%
  popResults_range(POC_AU, . , Population = c("Community", "Prisons"),
                   Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
  mutate(NP = "NP")

xt_pnp_rna <- cbind(year = Resflow_year_setting_range_trajectory$Testing_ab$year, 
                   population = Resflow_year_setting_range_trajectory$Testing_ab_sc$population, 
                   scenario = Resflow_year_setting_range_trajectory$Testing_ab$scenario, 
                   as.data.frame(Resflow_year_setting_range_trajectory$Testing_RNA[, par_col] + 
                                   Resflow_year_setting_range_trajectory$Testing_RNA_neg[, par_col] + 
                                   Resflow_year_setting_range_trajectory$Testing_POCT[, par_col] + 
                                   Resflow_year_setting_range_trajectory$Testing_POCT_neg[, par_col]))%>%
  as.data.frame()%>%
  popResults_range(POC_AU, . , Population = c("Community", "Prisons"),
                   Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)%>%
  mutate(NP = "out of NP")



xt_ab <- rbind(xt_np_ab, xt_pnp_ab)
xt_rna <- rbind(xt_np_rna, xt_pnp_rna)
xt_ab <- xt_ab%>%mutate(NP = factor(NP, levels = c("out of NP", "NP"),
                                    labels = c("Out of the National Program", 
                                               "National Program")))

xt_rna <- xt_rna%>%mutate(NP = factor(NP, levels = c("out of NP", "NP"),
                                    labels = c("Out of the National Program", 
                                               "National Program")))
xt_ab[, c(1:10)]
HCVNP_ab_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVNPab_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY , 
                           best = ab,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))%>%
  mutate(scenario = ifelse(year + POC_AU$cabY == 2024, "dfList_NP_2024" ,"dfList_NP_2023"),
         NP = "NP")%>%
  
  mutate( NP = factor(NP, 
                     labels = c(
                       "National program")))

HCVNP_rna_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVNPab_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year - POC_AU$cabY , 
                           best = RNA,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))%>%
  mutate(scenario = ifelse(year + POC_AU$cabY == 2024, "dfList_NP_2024" ,"dfList_NP_2023"),
         NP = "NP")%>%
  
  mutate(NP = factor(NP, 
                     labels = c(
                       "National Program")))
HCVNP_ab_setting_fit_all <- HCVNP_ab_setting_fit%>%ungroup()%>%
  group_by(year, scenario)%>%summarise(best= sum(best))%>%ungroup()
View(HCVNP_ab_setting_fit_all)
HCVNP_ab_setting_fit_lst <- list()
HCVNP_ab_setting_fit_lst_all <- list()
HCVNP_rna_setting_fit_lst <- list()
HCVNP_rna_setting_fit_lst_all <- list()
for(i in c("Past program", "Current program", 
           "Continued program","Expanded program")){ 
  HCVNP_ab_setting_fit_lst[[i]] <- as.data.frame(HCVNP_ab_setting_fit)%>%
    mutate(scenario = i)
 
  HCVNP_rna_setting_fit_lst[[i]] <- as.data.frame(HCVNP_rna_setting_fit)%>%
    mutate(scenario = i)

  }
HCVNP_ab_setting_fit_all <- list()
HCVNP_rna_setting_fit_all <- list()
for(i in names(HCVNP_ab_setting_fit_lst)){ 
  HCVNP_ab_setting_fit_all[[i]] <- HCVNP_ab_setting_fit_lst[[i]]%>%
    group_by(year, scenario)%>%summarise(best= sum(best)) 
  
  HCVNP_rna_setting_fit_all[[i]] <- HCVNP_rna_setting_fit_lst[[i]]%>%
  group_by(year, scenario)%>%summarise(best= sum(best))
  
  }


xt_ab <- xt_ab%>%
  mutate(scenario = factor(scenario, 
                           levels = c("No National Program", 
                                      "Past program", 
                                      "Current program", 
                                      "Continued program",
                                      "Expanded program")))%>%
  arrange(scenario)

xt_ab[is.na(xt_ab)] <- 0
xt_rna[is.na(xt_rna)] <- 0
parea_ab <- list()
parea_ab[["No National Program"]] <- 
  ggplot(data = xt_ab%>%filter(scenario == "No National Program"), aes(x = year, y = best)) +
  geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
            alpha=0.6 , size=0.2, colour="black") +
  scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
  facet_wrap(~population, scales='free_y') + 
  theme_Publication_facet() + 
  scale_x_continuous(limits = c(1,16), breaks = seq(1,16,5), 
                     labels = seq(1,16,5) + POC_AU$cabY - 1)  + 
  labs(x = "Year", y = "Number of antibody testing", title = "No National Program") + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 20000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 20000)))))

parea_rna <- list()
parea_rna[["Out of the National Program"]] <- 
  ggplot(data = xt_rna%>%filter(scenario == "No National Program"), aes(x = year, y = best)) +
  geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
            alpha=0.6 , size=0.2, colour="black") +
  scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
  facet_wrap(~population, scales='free_y') + 
  theme_Publication_facet() + 
  scale_x_continuous(limits = c(1,16), breaks = seq(1,16,5), 
                     labels = seq(1,16,5) + POC_AU$cabY - 1)  + 
  labs(x = "Year", y = "Number of RNA testing", title = "No National Program") + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 30000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 25000)))))
for(i in names(HCVNP_ab_setting_fit_lst)){ 
  parea_ab[[i]] <- ggplot(data = xt_ab%>%filter(scenario == i), 
                          aes(x = year, y = best)) +
    geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
              alpha=0.6 , size=0.2, colour="black") +
    scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
    geom_point(data=HCVNP_ab_setting_fit_lst[[i]], aes(y=best, x = time), 
               colour = "black", size = 1) + 
    facet_wrap(~population, scales='free_y') + 
    theme_Publication_facet() + 
    scale_x_continuous(limits = c(1,16), breaks = seq(1,16,5), 
                       labels = seq(1,16,5) + POC_AU$cabY - 1)  + 
    labs(x = "Year", y = "Number of antibody testing", title = i) + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 20000))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 20000)))))
  
  parea_rna[[i]] <- ggplot(data = xt_rna%>%filter(scenario == i), 
                          aes(x = year, y = best)) +
    geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
              alpha=0.6 , size=0.2, colour="black") +
    scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
    geom_point(data=HCVNP_rna_setting_fit_lst[[i]], aes(y=best, x = time), 
               colour = "black", size = 1) + 
    facet_wrap(~population, scales='free_y') + 
    theme_Publication_facet() + 
    scale_x_continuous(limits = c(1,16), breaks = seq(1,16,5), 
                       labels = seq(1,16,5) + POC_AU$cabY - 1)  + 
    labs(x = "Year", y = "Number of RNA testing", title = i) + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 30000))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 25000)))))
    
}
for(i in names(parea_ab)){ 
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numAb_area_", i,".png")), 
         parea_ab[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numRNA_area_", i ,".png")), 
         parea_rna[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
}
# remove y label 
parea_ab <- lapply(parea_ab, function(x) x + rremove("ylab") + rremove("xlab"))
parea_ab <- list(parea_ab[[1]],parea_ab[[2]], parea_ab[[3]])
parea_ab_ggarrange <- ggarrange(plotlist = parea_ab, ncol = 3, nrow = 1, 
                                common.legend = TRUE)
# adding common x, y label 
parea_ab_ggarrange <- annotate_figure(parea_ab_ggarrange, 
                left = textGrob("Number of antibody testing", 
                                rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Year", gp = gpar(cex = 1.3))
                )
ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numAb_arrange", ".png")), 
       parea_ab_ggarrange , 
       width = 12, height = 8, bg = "white", dpi = 300)
# remove y label 
parea_rna <- lapply(parea_rna, function(x) x + rremove("ylab") + rremove("xlab"))

parea_rna <- list(parea_rna[[1]],parea_rna[[2]], parea_rna[[3]])
parea_rna_ggarrange <- ggarrange(plotlist = parea_rna, ncol = 3, nrow = 1, 
                                common.legend = TRUE)
# adding common x, y label 
parea_rna_ggarrange <- annotate_figure(parea_rna_ggarrange, 
                                      left = textGrob("Number of RNA testing", 
                                                      rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                      bottom = textGrob("Year", gp = gpar(cex = 1.3)))


ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numRNA_arrange", ".png")), 
       parea_rna_ggarrange, 
       width = 12, height = 8, bg = "white", dpi = 300)


# total tests done in np 








ggsave(file=file.path(OutputFig_y_cum_avert, paste0("T_num" ,".png")), 
       p_T_num, 
       width = 9, height = 6, bg = "white", dpi = 300)

# total test 
xt_ab <- xt_ab%>%arrange(year, population, scenario, NP)%>%ungroup()
xt_rna <- xt_rna%>%arrange(year, population, scenario, NP)%>%ungroup()
xt_toltest <- cbind(year = xt_ab$year, 
                    population = xt_ab$population,
                    scenario = xt_ab$scenario, 
                    NP = xt_ab$NP,
                    as.data.frame(xt_ab[, c(par_col)] + 
                                    xt_rna[, c(par_col)]))%>%as.data.frame()%>%
  group_by(year, scenario, NP)%>%
  summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))

xt_toltest
xt_toltest_lst <- list()
for(i in unique(xt_toltest$scenario)){ 
  xt_toltest_lst[[i]] <- xt_toltest%>%filter(scenario == i)
  
}
xt_toltest_lst
xt_toltest_lst$`No National Program` <- 
  xt_toltest_lst$`No National Program`%>%ungroup()%>%
  mutate(dt = NA, dt_exp = NA)
xt_toltest_lst$`Past program` <- 
  xt_toltest_lst$`Past program`%>%ungroup()%>%
  mutate(dt = ifelse(year == 7 & NP != "National Program", 6615, 
                     ifelse(year == 8 & NP != "National Program", 11276, NA)), 
         dt_exp = NA)

xt_toltest_lst$`Current program` <- 
  xt_toltest_lst$`Current program`%>%ungroup()%>%
  mutate(dt = ifelse(year == 7 & NP != "National Program", 6615, 
                     ifelse(year == 8 & NP != "National Program", 11276, NA)), 
         dt_exp = ifelse(year %in% c(9) & NP != "National Program", 19480, NA))

xt_toltest_lst$`Continued program` <- 
  xt_toltest_lst$`Continued program`%>%ungroup()%>%
  mutate(dt = ifelse(year == 7 & NP != "National Program", 6615, 
                     ifelse(year == 8 & NP != "National Program", 11276, NA)), 
         dt_exp = ifelse(year %in% c(9) & NP != "National Program", 19480, NA))
                     
xt_toltest_lst$`Expanded program` <- 
  xt_toltest_lst$`Expanded program`%>%ungroup()%>%
  mutate(dt = ifelse(year == 7 & NP != "National Program", 6615, 
                     ifelse(year == 8 & NP != "National Program", 11276, NA)), 
         dt_exp = ifelse(year %in% c(9) & NP != "National Program", 19480,
                         ifelse(year %in% c(10) & NP != "National Program", 30000, 
                                ifelse(year %in% c(11) & NP != "National Program", 40000,
                                       ifelse(year %in% c(12) & NP != "National Program", 50000, NA)))))


parea_tol <- list()
for(i in names(HCVNP_ab_setting_fit_lst)){ 
  if(i == "Past program"){ 
    parea_tol[[i]] <- ggplot(xt_toltest_lst[[i]]) +
      geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
                alpha=0.6 , size=0.2, colour="black") +
      scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
      geom_point(aes(y=dt, x = year), colour = "black", size = 2
      )  +
      theme_Publication() + 
      scale_x_continuous(limits = c(1,16), breaks = seq(1,16,5), 
                         labels = seq(1,16,5) + POC_AU$cabY - 1)  + 
      labs(x = "Year", y = "Number of total tests (thousands)", title = i) + 
      scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 10000),
                         labels = seq(0, 100000, 10000)/1000)
  } 
  else{ 
    parea_tol[[i]] <- ggplot(xt_toltest_lst[[i]]) +
      geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
                alpha=0.6 , size=0.2, colour="black") +
      scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
      geom_point(aes(y=dt, x = year), colour = "black", size = 2
      ) + 
      geom_point(aes(y=dt_exp, x = year) , shape = 1, colour = "black", size = 2
      ) +
      theme_Publication() + 
      scale_x_continuous(limits = c(1,16), breaks = seq(1,16,5), 
                         labels = seq(1,16,5) + POC_AU$cabY - 1)  + 
      labs(x = "Year", y = "Number of total tests (thousands)", title = i) + 
      scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 10000),
                         labels = seq(0, 100000, 10000)/1000)
    }
    
  }

parea_tol <- lapply(parea_tol, function(x) x + rremove("ylab") + rremove("xlab"))


parea_tol_ggarrange <- ggarrange(plotlist = parea_tol, ncol = 4, nrow = 1, 
                                 common.legend = TRUE)
# adding common x, y label 
parea_tol_ggarrange <- annotate_figure(parea_tol_ggarrange, 
                                       left = textGrob("Number of total tests (thousands)", 
                                                       rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                       bottom = textGrob("Year", gp = gpar(cex = 1.3)))


ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numtotal_arrange", ".png")), 
       parea_tol_ggarrange, 
       width = 12, height = 8, bg = "white", dpi = 300)


################################################################################
# number of advanced liver diseases 
Res_Numbox_y <- list()
Res_Numbox_cum <- list()
Res_Numbox_avert <- list()
# output in timestep finding mid-year 
for(i in names(Res_numbox)){ 
  Res_Numbox_y[[i]][["DC"]] <- Res_numbox[[i]]%>%filter(disease_prog == "dc")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
    filter(timestep%%1 == 0.5)%>%mutate(year = timestep%/%1)%>%
    select(year, par_col)
  
  Res_Numbox_y[[i]][["HCC"]] <- Res_numbox[[i]]%>%filter(disease_prog == "hcc")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
    filter(timestep%%1 == 0.5)%>%mutate(year = timestep%/%1)%>%
    select(year, par_col)
  
  Res_Numbox_y[[i]][["LT"]] <- Res_numbox[[i]]%>%filter(disease_prog == "lt")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%
    filter(timestep%%1 == 0.5)%>%mutate(year = timestep%/%1)%>%
    select(year, par_col)
  
  Res_Numbox_y[[i]][["PLT"]] <- Res_numbox[[i]]%>%filter(disease_prog == "plt")%>%
    group_by(timestep)%>%
    summarise(across(c(par_col),~ mean(.x, na.rm = FALSE)))%>%
    filter(timestep%%1 == 0.5)%>%mutate(year = timestep%/%1)%>%
    select(year, par_col)
  
  }

for(i in names(Res_Numbox_y)){
  for(indic in names(Res_Numbox_y[[1]])){ 
    Res_Numbox_cum[[i]][[indic]] <- Res_Numbox_y[[i]][[indic]]%>%
      mutate(year = year + POC_AU$cabY - 1)%>%
      filter(year >= POC_AU$simY )%>%ungroup()%>%
      mutate(across(c(par_col), cumsum, .names = "{col}"))
      
  }
}

for(i in names(Res_Numbox_y)){ 
  for(indic in names(Res_Numbox_y[[1]])){
    Res_Numbox_avert[[i]][[indic]] <- cbind(
      year = Res_Numbox_cum[[i]][[indic]]$year, 
      data.frame(Res_Numbox_cum[["Status quo"]][[indic]][, c(par_col)] - 
                   Res_Numbox_cum[[i]][[indic]][, c(par_col)]))
  }
}
Res_Numbox_y_range <- list()
Res_Numbox_cum_range <- list()
Res_Numbox_avert_range <- list()
for(i in names(Res_Numbox_y)){ 
  for(indic in names(Res_Numbox_y[[1]])){
    Res_Numbox_y_range[[i]][[indic]] <- Res_Numbox_y[[i]][[indic]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    Res_Numbox_cum_range[[i]][[indic]] <- Res_Numbox_cum[[i]][[indic]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    Res_Numbox_avert_range[[i]][[indic]] <- Res_Numbox_avert[[i]][[indic]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  }

}

Res_Numbox_y_range <- Res_Numbox_y_range%>%purrr::transpose()
Res_Numbox_cum_range <- Res_Numbox_cum_range%>%purrr::transpose()
Res_Numbox_avert_range <- Res_Numbox_avert_range%>%purrr::transpose()

for(i in names(Res_Numbox_y_range)){ 
  Res_Numbox_y_range[[i]] <- dplyr::bind_rows(Res_Numbox_y_range[[i]], .id = "scenario")%>%
    mutate(year = year + POC_AU$cabY - 1)%>%
    mutate(scenario = factor(scenario, 
                             levels = c("Status quo", "dfList_NP_2023", 
                                        "dfList_NP_2024", "dfList_NPexp_A", 
                                        "dfList_NPexp_B", "dfList_NPexp_C",
                                        "dfList_NPexp_D"), 
                             labels = c("No National Program", "Past program", 
                                        "Current program", "Continued program", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "Expanded program")))%>%
    filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
  
  Res_Numbox_cum_range[[i]] <- dplyr::bind_rows(Res_Numbox_cum_range[[i]], .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = c("Status quo", "dfList_NP_2023", 
                                        "dfList_NP_2024", "dfList_NPexp_A", 
                                        "dfList_NPexp_B", "dfList_NPexp_C",
                                        "dfList_NPexp_D"), 
                             labels = c("No National Program", "Past program achievement", 
                                        "Current program achievement", "Continued program achievement", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "Expanded program achievement")))%>%
    filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
  Res_Numbox_avert_range[[i]] <- dplyr::bind_rows(Res_Numbox_avert_range[[i]], .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = c("Status quo", "dfList_NP_2023", 
                                        "dfList_NP_2024", "dfList_NPexp_A", 
                                        "dfList_NPexp_B", "dfList_NPexp_C",
                                        "dfList_NPexp_D"), 
                             labels = c("No National Program", "Past program", 
                                        "Current program", "Continued program", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "Expanded program")))%>%
    filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
    
  }

for(i in names(Res_Numbox_avert_range)){ 
  Res_Numbox_avert_range[[i]] <- Res_Numbox_avert_range[[i]]%>%
    filter(scenario != "No National Program")
  }

lim_adliver_y <- list()
lim_adliver_cum <- list()
lim_adliver_avert <- list()
for(i in names(Res_Numbox_y_range)){ 
  
  lim_adliver_y[[i]] <- lim_ident(Res_Numbox_y_range[[i]], seq(2021, 2050, 1))
  lim_adliver_cum[[i]] <- lim_ident(Res_Numbox_cum_range[[i]], seq(2021, 2050, 1))
  lim_adliver_avert[[i]] <- lim_ident(Res_Numbox_avert_range[[i]], seq(2021, 2050, 1))
  }


p_num_adliver_y <- list()
p_num_adliver_cum <- list()
p_num_adliver_avert <- list()
y_adlab_name <- list("Decompensated cirrhosis", 
                     "Hepatocellular carcinoma", 
                     "Liver transplant", 
                     "Post-liver transplant")
names(y_adlab_name) <- names(Res_Numbox_y_range)

for(i in names(Res_Numbox_y_range)){ 
  
  p_num_adliver_y[[i]] <- plot_pocau(POC_AU, 
                                     Res_Numbox_y_range, type = "new", 
                                     indicator = i) +
    scale_x_continuous(limits = c(2022, 2041), breaks = c(seq(2022, 2037, 5), 2041)) + 
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_y[[i]][, "lim"]))) + 
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]], x = "Year")  + 
    theme(theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  p_num_adliver_cum[[i]] <- plot_pocau(POC_AU, Res_Numbox_cum_range, type = "cum", 
             indicator = i, year_obs = year_obs) + 
    
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_cum[[i]][, "lim"]))) +
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]], x = "Year") 
  
  p_num_adliver_avert[[i]] <- plot_pocau(POC_AU, Res_Numbox_avert_range, type = "avert", 
                                       indicator = i, year_obs = year_obs) + 
    
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_avert[[i]][, "lim"]))) +
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]] , x = "Year") 
  }

p_num_adliver_y$DC <- p_num_adliver_y$DC + scale_y_continuous(limits = c(0, 2000))
p_num_adliver_y$HCC <- p_num_adliver_y$HCC + scale_y_continuous(limits = c(0, 500))
p_num_adliver_y$LT <- p_num_adliver_y$LT + scale_y_continuous(limits = c(0, 100))
p_num_adliver_y$PLT <- p_num_adliver_y$PLT + scale_y_continuous(limits = c(0, 100))

p_num_adliver_avert$DC <- p_num_adliver_avert$DC + 
  scale_y_continuous(limits = c(0, 2000)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") + 
  labs(y = "Number averted of Decompensated cirrhosis")

p_num_adliver_avert$HCC <- p_num_adliver_avert$HCC + 
  scale_y_continuous(limits = c(-5, 400)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  labs(y = "Number averted of Hepatocellular carcinoma")

p_num_adliver_avert$LT <- p_num_adliver_avert$LT + 
  scale_y_continuous(limits = c(-5, 100)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  labs(y = "Number averted of liver transplant")

p_num_adliver_avert$PLT <- p_num_adliver_avert$PLT + 
  scale_y_continuous(limits = c(-5, 50)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  labs(y = "Number averted of post-liver transplant") 
  
p_y_adliver_arrange <- ggarrange(plotlist = p_num_adliver_y, ncol = 2, nrow = 2, 
                                 common.legend = TRUE, legend = "bottom")
p_cum_adliver_arrange <- ggarrange(plotlist = p_num_adliver_cum, ncol = 2, nrow = 2, 
                                 common.legend = TRUE, legend = "bottom")
p_avert_adliver_arrange <- ggarrange(plotlist = p_num_adliver_avert, ncol = 2, nrow = 2, 
                                   common.legend = TRUE, legend = "bottom")

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("num_adliver_y"  ,".png")), 
       p_y_adliver_arrange, 
       width = 8, height = 10, bg = "white", dpi = 300)

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("cumnum_adliver",".png")), 
       p_cum_adliver_arrange, 
       width = 8, height = 10, bg = "white", dpi = 300)  

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("avert_adliver_", ".png")), 
       p_avert_adliver_arrange , 
       width = 8, height = 11, bg = "white", dpi = 300)


# save individual figures 
for(i in names(p_num_adliver_y)){ 
  ggsave(file=file.path(OutputFig_y_cum_avert, 
                        paste0("num_adliver_",i  ,".png")), 
         p_num_adliver_y[[i]], 
         width = 10, height = 6, bg = "white", dpi = 300)
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("cumnum_adliver_",i  ,".png")), 
         p_num_adliver_cum[[i]], 
         width = 10, height = 6, bg = "white", dpi = 300)
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("avert_adliver_",i  ,".png")), 
         p_num_adliver_avert[[i]], 
         width = 10, height = 6, bg = "white", dpi = 300)
  
  }

################################################################################

####                Simple and easy to understand plots                     #### 


################################################################################
# plots: x: 5y, 10y, 20y; y: the averted number of {a,b,c,d} in each scenario
# (a) newinfections; (b) HCV-related Deathes; (a) DC; (b) HCC; (c) LT; (d) PLT 
# (a) & (b) flows 
# (c) to (f) Num_box

# (a) to (b) 
# first select the list of dtset we would like to generate the plots, which is aligned with the sce_lab below 
# chose the indicator want to generate plots: new infections and HCV releated death 
Resflow_cumavert <- Resflow_cum_all_avert_range[names(Resflow_cum_all_avert_range)%in% 
                                   c("dfList_NP_2023", "dfList_NP_2024", "dfList_NPexp_A",
                                    "dfList_NPexp_D")]
nindic <- list("newInfections" = "New HCV Infections", "HCVdeath" = "HCV-related Deathes")
plot_num_avert <- list()
nlim_avert <- list("newInfections" = 20000, "HCVdeath" = 500)
nnumtag <- list("a", "b")
names(nnumtag) <- names(nindic)
for(indic in names(nindic)){ 
  for(s in names(Resflow_cumavert)){ 
    plot_num_avert[[s]][[indic]] <- 
      ggplot(Resflow_cumavert[[s]][[indic]]%>%
               filter(year %in% (year_obs)), 
             aes(x = as.character(year), y = best)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=round(best, digits = 0)), position = position_dodge(width = 0.55),
                hjust = 0.5, vjust = -0.5) + 
      scale_x_discrete(labels = c(paste0((year_obs - POC_AU$simY + 1), "-Year", sep = ""))) + 
      labs(y = "Number averted", x = "Time frame", tag = nnumtag[[indic]]) + 
      ggtitle(nindic[[indic]]) + 
      scale_y_continuous(limits = c(0,nlim_avert [[indic]])) +
      theme_Publication() 
    }
  
  }

sce_lab <- c(unique(Res_Numbox_avert_range$DC$scenario))

names(plot_num_avert) <- sce_lab
ntitle_avert <- list("Decompensated cirrhosis", "Hepatocellular carcinoma", 
                     "Liver transplant", "Post-liver transplant") 
ntag <- list("a", "b", "c", "d")
names(ntag) <- names(Res_Numbox_avert_range)
names(ntitle_avert) <- names(Res_Numbox_avert_range)
lim_avert <- list(30, 10, 5, 5 )
names(lim_avert) <- names(Res_Numbox_avert_range)
# (c) to (f)

for(indic in names(Res_Numbox_avert_range)){
  for(s in sce_lab){ 
    plot_num_avert[[s]][[indic]] <- 
      ggplot(Res_Numbox_avert_range[[indic]]%>%
               filter(scenario == s & year %in% (year_obs)), 
             aes(x = as.character(year), y = best)) + 
      geom_bar(stat="identity") + 
      geom_text(aes(label=round(best, digits = 0)),position = position_dodge(width = 0.55),
                hjust = 0.5, vjust = -0.5) + 
      scale_x_discrete(labels = c(paste0((year_obs - POC_AU$simY + 1), "-Year", sep = ""))) + 
      labs(y = "Number averted", x = "Time frame", tags = ntag[[indic]]) + 
      ggtitle(ntitle_avert[[indic]]) +
      theme_Publication() + 
      scale_y_continuous(limits = c(0,lim_avert[[indic]]))
    
  }
} 
plot_num_avert$`Expanded program`$DC + scale_y_continuous(limits = c(0,1000))
# adjust achievement 2023 
plot_num_avert$`Past program`$newInfections <- 
  plot_num_avert$`Past program`$newInfections + scale_y_continuous(limits = c(0,5000))
plot_num_avert$`Past program`$HCVdeath <- 
  plot_num_avert$`Past program`$HCVdeath + scale_y_continuous(limits = c(0,100))
for(i in names(plot_num_avert)){ 
  plot_num_avert[[i]]$DC <- plot_num_avert[[i]]$DC + scale_y_continuous(limits = c(0,1000))
  
  plot_num_avert[[i]]$HCC <- plot_num_avert[[i]]$HCC + scale_y_continuous(limits = c(0,100))
  
  plot_num_avert[[i]]$LT <- plot_num_avert[[i]]$LT + scale_y_continuous(limits = c(0,50))
  
  plot_num_avert[[i]]$PLT <- plot_num_avert[[i]]$PLT + scale_y_continuous(limits = c(0,5))
  }

## arrange the plot (a) to (f) adding common title and tag 

plot_numavert <- list() 
plot_flowavert <- list()
for(i in names(plot_num_avert)){ 
  plot_flowavert[[i]] <- ggarrange(plotlist = 
                                     list(plot_num_avert[[i]]$newInfections, 
                                       plot_num_avert[[i]]$HCVdeath),
                                  ncol = 2, nrow = 1, common.legend = TRUE,
                                  legend="bottom")
  
  plot_flowavert[[i]] <- annotate_figure(plot_flowavert[[i]], top = text_grob(i, 
                                                              color = "black", 
                                                              face = "bold", size = 14))
  
  plot_numavert[[i]] <- ggarrange(plotlist = 
                                     list(plot_num_avert[[i]]$DC, 
                                          plot_num_avert[[i]]$HCC, 
                                          plot_num_avert[[i]]$LT, 
                                          plot_num_avert[[i]]$PLT),
                                   ncol = 2, nrow = 2, common.legend = TRUE,
                                   legend="bottom")
  
  plot_numavert[[i]] <- annotate_figure(plot_numavert[[i]], 
                                        top = text_grob(i, 
                                                        color = "black", 
                                                        face = "bold", size = 14))
  
}

plot_numavert$`Past program achievement`
#dir.create(file.path(paste0(OutputFig, "/Reports")))
for(i in names(plot_flowavert)){ 
  
  ggsave(file=file.path(OutputFig, paste0("Reports/flow_cumavert_",i,".png")), 
         plot_flowavert[[i]], 
         width = 8 , height = 6, bg = "white")
  ggsave(file=file.path(OutputFig, paste0("Reports/box_cumavert_",i,".png")), 
         plot_numavert[[i]], 
         width =12, height = 8, bg = "white")
  
  }


################################################################################

# ggarrange plot 
#plot_flow_arrange <- list()
#plot_flow_arrange$Treatment_sc
#for(i in 1: length(plot_flow_year)){ 
#  plot_flow_arrange[[i]] <- ggpubr::ggarrange(
#    plotlist = list(plot_flow_year[[i]], 
#                    plot_flow_cum[[i]],
#                    plot_flow_cumavert[[i]]),
#    common.legend = TRUE, 
#    nrow=1)  + 
#    theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm")) 
#  
#  plot_flow_arrange[[i]] <- annotate_figure(plot_flow_arrange[[i]], 
#                                            top = text_grob(names(plot_flow_year)[i], color = "black", 
#                                                            face = "bold", size = 14))
#  
#  }

#####################################
# cost data plot 
# extract yearly value 
cap <- 200000000
cost_year_all <- list()


names(Rescost_year_all$dfList_NP_2023)

cost_y_categories <- list()
cost_disyear_categories <- list()

  
for(i in names(Rescost_year_all)){ 
  
  Rescost_year_all[[i]]$cost_POCT[is.na(Rescost_year_all[[i]]$cost_POCT)]  <- 0

    cost_y_categories[[i]][["Diagnosis"]] <- cbind(
      year = Rescost_year_all[[i]]$cost_ab$year, 
      as.data.frame(Rescost_year_all[[i]]$cost_ab[, par_col],
                         Rescost_year_all[[i]]$cost_RNA[, par_col],
                         Rescost_year_all[[i]]$cost_POCT[, par_col]))
 
    cost_y_categories[[i]][["Diagnosis"]][cost_y_categories[[i]][["Diagnosis"]] == 0] <- NA 
    cost_y_categories[[i]][["Diagnosis"]] <- cost_y_categories[[i]][["Diagnosis"]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  
    cost_y_categories[[i]][["Treatment_cap"]] <- Rescost_year_all[[i]]$cost_totalDAA_Cap%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
    cost_y_categories[[i]][["Treatment"]] <- Rescost_year_all[[i]]$cost_totalDAA%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
    
    Rescost_year_all[[i]]$cost_compartment[is.na(Rescost_year_all[[i]]$cost_compartment)]  <- 0
    Rescost_year_all[[i]]$cost_Cured[Rescost_year_all[[i]]$cost_Cured <0 ] <- 0
    Rescost_year_all[[i]]$cost_TreatOther[Rescost_year_all[[i]]$cost_TreatOther <0 ] <- 0
    Rescost_year_all[[i]]$cost_RetreatOther[Rescost_year_all[[i]]$cost_RetreatOther <0 ] <- 0
    cost_y_categories[[i]][["Management"]] <- cbind(
      year = Rescost_year_all[[i]]$cost_compartment$year, 
      as.data.frame(Rescost_year_all[[i]]$cost_compartment[, par_col] + 
                         Rescost_year_all[[i]]$cost_Cured[, par_col] +
                         Rescost_year_all[[i]]$cost_TreatOther[, par_col] + 
                         Rescost_year_all[[i]]$cost_RetreatOther[, par_col])) 
    
    
    cost_y_categories[[i]][["Management"]][cost_y_categories[[i]][["Management"]] == 0] <- NA  
    cost_y_categories[[i]][["Management"]] <- cost_y_categories[[i]][["Management"]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
    cost_y_categories[[i]] <- cost_y_categories[[i]]%>%
      dplyr::bind_rows(., .id = "Categories")
    
}
for(i in names(Rescost_year_all)){     
    Rescost_disyear_all[[i]]$cost_ab[is.na(Rescost_disyear_all[[i]]$cost_ab)]  <- 0
    Rescost_disyear_all[[i]]$cost_RNA[is.na(Rescost_disyear_all[[i]]$cost_RNA)]  <- 0
    Rescost_disyear_all[[i]]$cost_POCT[is.na(Rescost_disyear_all[[i]]$cost_POCT)]  <- 0
    
    cost_disyear_categories[[i]][["Diagnosis"]] <- cbind(
      year = Rescost_disyear_all[[i]]$cost_ab$year, 
      as.data.frame(Rescost_disyear_all[[i]]$cost_ab[, par_col] + 
                         Rescost_disyear_all[[i]]$cost_RNA[, par_col] +
                         Rescost_disyear_all[[i]]$cost_POCT[, par_col]))
    
    cost_disyear_categories[[i]][["Diagnosis"]][cost_disyear_categories[[i]][["Diagnosis"]] == 0] <- NA  
    
    cost_disyear_categories[[i]][["Diagnosis"]] <-  cost_disyear_categories[[i]][["Diagnosis"]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
    
    cost_disyear_categories[[i]][["Treatment_cap"]] <- Rescost_disyear_all[[i]]$cost_totalDAA_Cap%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
    cost_disyear_categories[[i]][["Treatment"]] <- Rescost_disyear_all[[i]]$cost_totalDAA%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
    Rescost_disyear_all[[i]]$cost_compartment[is.na(Rescost_disyear_all[[i]]$cost_compartment)]  <- 0
    Rescost_disyear_all[[i]]$cost_Cured[Rescost_disyear_all[[i]]$cost_Cured <0 ] <- 0
    Rescost_disyear_all[[i]]$cost_TreatOther[Rescost_disyear_all[[i]]$cost_TreatOther <0 ] <- 0
    Rescost_disyear_all[[i]]$cost_RetreatOther[Rescost_disyear_all[[i]]$cost_RetreatOther <0 ] <- 0
    
    cost_disyear_categories[[i]][["Management"]] <- cbind(
      year = Rescost_disyear_all[[i]]$cost_compartment$year, 
      as.data.frame(Rescost_disyear_all[[i]]$cost_compartment[, par_col] + 
                         Rescost_disyear_all[[i]]$cost_Cured[, par_col] +
                         Rescost_disyear_all[[i]]$cost_TreatOther[, par_col] + 
                         Rescost_disyear_all[[i]]$cost_RetreatOther[, par_col]))


    cost_disyear_categories[[i]][["Management"]][cost_disyear_categories[[i]][["Management"]] == 0] <- NA
    
    cost_disyear_categories[[i]][["Management"]] <- cost_disyear_categories[[i]][["Management"]]%>%
      popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
    cost_disyear_categories[[i]] <- cost_disyear_categories[[i]]%>%
      dplyr::bind_rows(., .id = "Categories")
}


y_cost_disyear_categories <- cost_disyear_categories%>%
  dplyr::bind_rows(., .id = "scenario")%>%ungroup()%>%
  group_by(scenario, year)%>%
  dplyr::summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))

y_cost_disyear_categories_range <- y_cost_disyear_categories%>%
  gather("simulation", "estimate", -c(year, scenario))%>%
  group_by(year, scenario)%>%
  summarise(min = min(estimate, na.rm = TRUE),
            max = max(estimate, na.rm = TRUE),
            Med = median(estimate, na.rm = TRUE),
            Mu = mean(estimate, na.rm = TRUE),
            q5 = quantile(estimate, prob = 0.025, na.rm = TRUE),
            q25 = quantile(estimate, prob = 0.25, na.rm = TRUE),
            q75 = quantile(estimate, prob = 0.75, na.rm = TRUE),
            q95 = quantile(estimate, prob = 0.975, na.rm = TRUE)
  )%>%ungroup

ggplot(y_cost_disyear_categories, aes(x = year, colour = scenario)) + 
  geom_line(aes(y = best)) + 
  scale_x_continuous(expand = c(0,0),limits = c(2022, 2045), breaks = seq(2022, 2045, 1)) + 
  geom_ribbon(data = y_cost_disyear_categories_range, aes(x = year, ymin = q25, ymax = q75, 
                                                          fill = scenario), 
              alpha =0.2, colour = NA) + 
  scale_y_continuous(limits = c(0, 1000000000), breaks = seq(0, 1000000000,100000000), 
                     labels = seq(0, 1000000000,100000000)/1000000000) + 
  theme_Publication()
ref_sce <- y_cost_disyear_categories%>%filter(scenario == "Status quo") 

y_cost_disyear_categories <- y_cost_disyear_categories%>%mutate(best_turning = best - ref_sce$best)
y_cost_disyear_categories <- y_cost_disyear_categories%>%select(scenario, year, best_turning,
                                                             par_col)

xt <- y_cost_disyear_categories%>%group_by(scenario)%>%
  filter(best_turning <0)%>%slice(1)

unique(y_cost_disyear_categories$scenario)
y_cost_disyear_categories <- y_cost_disyear_categories%>%
  mutate(scenario = factor(scenario, 
                           levels = c("Status quo", "dfList_NP_2023", 
                                      "dfList_NP_2024", "dfList_NPexp_A", 
                                      "dfList_NPexp_B", "dfList_NPexp_C",
                                      "dfList_NPexp_D"), 
                           labels = c("No National Program", "Past program", 
                                      "Current program", "Continued program", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "Expanded program")))%>%
  filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
# adding dashed line to indicate the turning point of year 
# change colours for scenarios 
col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
p_cost_y_turning <- ggplot(y_cost_disyear_categories, aes(x = year, colour = scenario)) + 
  geom_line(aes(y = best)) + 
  scale_x_continuous(expand = c(0,0),limits = c(2022, 2045), breaks = seq(2022, 2045, 1)) + 
  
  scale_y_continuous(expand = c(0,0), limits = c(0, 700000000), breaks = seq(0, 700000000,100000000), 
                     labels = seq(0, 700000000,100000000)/1000000000) + 
  theme_Publication() + 
  
  scale_colour_manual(name = "Scenarios", values = col_pal) + 
  labs( x = "Year", y = "Annual discounted costs, billions")

ggsave(file=file.path(OutputFig_y_cum_avert, paste0("p_cost_y_turning",".png")), 
       p_cost_y_turning, 
       width = 12, height = 8, bg = "white", dpi = 300)
# plot of cost categories to 2080 
# undiscount no cap 
cost_y_categories <- dplyr::bind_rows(cost_y_categories, .id = "scenario")

cost_y_categories <- cost_y_categories%>%
  mutate(scenario = factor(scenario, 
                           levels = c("Status quo", "dfList_NP_2023", 
                                      "dfList_NP_2024", "dfList_NPexp_A", 
                                      "dfList_NPexp_B", "dfList_NPexp_C",
                                      "dfList_NPexp_D"), 
                           labels = c("No National Program", "Past program", 
                                      "Current program", "Continued program", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "Expanded program")))%>%
  filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))

cost_ydaanocap_categories <- cost_y_categories%>%filter(Categories != "Treatment_cap")


# undiscount cap 

cost_ydaacap_categories <- cost_y_categories%>%filter(Categories != "Treatment")

# discount no cap 
cost_disyear_categories <- dplyr::bind_rows(cost_disyear_categories, .id = "scenario")



cost_disydaanocap_categories <- cost_disyear_categories%>%filter(Categories != "Treatment_cap")

# discount cap 
cost_disydaacap_categories <- cost_disyear_categories%>%filter(Categories != "Treatment")

View(cost_disyear_categories)


write.xlsx(cost_ydaacap_categories%>%
             select(scenario, Categories, year, best, min, max, 
                    Med, Mu, q5, q25, q75, q95), file = file.path(OutputFig, paste0("Reports/cost_y_daacap.xlsx")), 
           append=TRUE) 
write.xlsx(cost_disydaanocap_categories%>%
             select(scenario, Categories, year, best, min, max, 
                    Med, Mu, q5, q25, q75, q95), file = file.path(OutputFig, paste0("Reports/cost_disy_daacap.xlsx")), 
           append=TRUE) 

write.xlsx(cost_ydaanocap_categories%>%
             select(scenario, Categories, year, best, min, max, 
                    Med, Mu, q5, q25, q75, q95), file = file.path(OutputFig, paste0("Reports/cost_y_nocap.xlsx")), 
           append=TRUE) 
# gt_table: 4 tables by categories
# categories yearly cost and discount yearly cost to 2022- 2080 
# columns: scenarios 
4948698700	- 5019337375
# output excel files 

# benefit: Lifetime cost averted = total lifetime cost_ref -  total lifetime cost_program 
# cost: program cost: 5y

# diagnosis cost 
x_catcost <- cost_disydaacap_categories%>%
  filter(year>= 2022)%>%
  group_by(scenario, Categories)%>%
  mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                  "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
  arrange(scenario)%>%
  select(scenario, Categories, year, best, q5, q95)%>%
  mutate(scenario = factor(scenario, levels = c("Status quo", "dfList_NP_2023", 
                                                "dfList_NP_2024", "dfList_NPexp_A", 
                                                "dfList_NPexp_B", "dfList_NPexp_C",
                                                "dfList_NPexp_D"), 
                           labels = c("No National Program", "Past program", 
                                      "Current program", "Continued program", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "Expanded program")))%>%
  filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))

# cumulative_cost plots 
col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
x_catcost%>%
  ggplot(. , aes(x = year, y = best)) +
  geom_line(aes(colour = scenario)) + 
  scale_colour_manual( values = col_pal) + 
  scale_x_continuous(limits = c(2022, 2052), breaks = seq(2022,2052,5)) + 
  theme_Publication_facet() + 
  facet_custom (~Categories,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 70000000),
                                                 breaks = seq(0, 70000000, 7000000),
                                                 labels = seq(0, 70000000, 7000000)/1000000)),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 1800000000),
                                                 breaks = seq(0, 1800000000, 100000000),
                                                 labels = seq(0, 1800000000, 100000000)/1000000)),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 3000000000),
                                                 breaks = seq(0, 3000000000, 100000000),
                                                 labels = seq(0, 3000000000, 100000000)/1000000))
                  ))  + 
  theme(axis.line = element_line(),
        panel.margin = unit(2, "lines")) + 
  labs(y = "Costs (discounted, millions)") 
x_totalcost <- x_catcost%>%group_by(scenario, year)%>%
  summarise(across(c("best", "q5", "q95"), ~ sum(.x, na.rm = FALSE)))%>%
  mutate(Categories ="Total")

# barchart 
pcatcost <- list()
title_name <- c("5-Year: 2022-2026", 
                "10-Year: 2022-2031",
                "20-Year: 2022-2041")

for(i in seq_along(year_obs)){ 
  pcatcost[[i]]<- x_catcost%>%filter(year == year_obs[i] )%>%arrange(Categories)%>%
    ggplot(., aes(fill = Categories, y = best, x = scenario, 
                  label = round(best/1000000, digits = 1))) + 
    geom_bar(position="stack", stat="identity") + 
    theme(panel.spacing = unit(0, 'lines')) +
    scale_fill_manual(values = c( "grey10", "grey40","grey80")) + 
    theme_Publication(base_size = 16) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    scale_y_continuous(limit = c(0, 5000000000), 
                       breaks = seq(0, 5000000000, 500000000),
                       labels = seq(0, 5000000000, 500000000)/1000000) + 
    labs(y = "Cost (discounted, millions)") + 
    geom_text(aes(x = scenario, y = best + 50000000, 
                  label = paste0(format(round(best/1000000, digits = 1), nsmall = 1), "m"), 
                  group = Categories),
              position = position_stack(vjust = 0.5), size = 6) + 
    ggtitle(title_name[i])
  }

ggsave(file=file.path(OutputFig, paste0("Reports/cost_catego_20y",".png")), 
       pcatcost[[3]], 
       width = 10, height = 8, bg = "white", dpi = 300)  
pcatcost <- lapply(pcatcost, function(x) x + rremove("ylab") + rremove("xlab"))

pcatcost  <- ggarrange(plotlist = pcatcost , ncol = 3, nrow = 1, 
                                 common.legend = TRUE)
# adding common x, y label 
pcatcost <- annotate_figure(pcatcost, 
                                       left = textGrob("Cost (discounted, millions)", 
                                                       rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                       bottom = textGrob("Scenarios", gp = gpar(cex = 1.3)))

ggsave(file=file.path(OutputFig, paste0("Reports/cost_catego",".png")), 
       pcatcost, 
       width = 18, height = 10, bg = "white", dpi = 300)  


################################################################################ 

####                                  animation                             ####
# total cost turning point 
# x_catcost_total <- x_catcost%>%group_by(year, scenario)%>%
#  summarise(across(c("best", "q5", "q95"),~ sum(.x, na.rm = FALSE)))

# x_catcost_total%>%filter(year <2030 )%>%
#  ggplot(. ,aes(x = scenario, y = best/1000000000, colour = scenario)) + 
#  geom_point() + 
#  facet_wrap(~year, scale = "free")
# View(x_catcost_total)

# x_catcost_total <- x_catcost_total%>%mutate(year = as.integer(year))%>%
#  filter(year <=2050)
# p_plot <- ggplot(x_catcost_total,aes(x = year, y = best, group = scenario, color = scenario)) +
#  geom_line() + 
#  geom_point( position = position_dodge(width = 0.5)) +
#  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
#  labs(x = "Year", y = "Cumulative discounted cost") +
#  theme(legend.position = "top") + theme_Publication()
# p_plot <- p_plot  + transition_reveal(year)

  
# animate(p_plot, end_pause = 10, width=1000, height=600,fps = 5) 
  
#  ggplot(x_catcost_total,
#       aes(x=scenario, y=best, label =scenario, color = scenario)) + 
#  geom_point(stat = "identity", size = 2) 
# p_plot + 
#  transition_time(year) +
#  labs(title = "Year: {frame_time}")
# p <- p_plot + theme(legend.position = "none") + #removing legend
#  labs (title = "Year: {frame_time}", x= "Year", y = "Discounted cumulative cost") + # The animation codes starts here
#  transition_time(year) +
# ease_aes('linear') + theme_bw()
# animate(p, duration = 40, fps = 20, width = 500, height = 500, renderer = gifski_renderer())

################################################################################ 
######################## incremental line ######################################
x_total_ref <- cost_disydaacap_categories%>%
  filter(year>= 2022)%>%
  group_by(scenario, Categories)%>%
  mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                  "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
  arrange(scenario)%>%group_by(year, scenario)%>%
  summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))%>%filter(scenario == "Status quo") 

x_catcost_total <- cost_disydaacap_categories%>%
  filter(year>= 2022)%>%
  group_by(scenario, Categories)%>%
  mutate(across(c(par_col, "min", "max", "Med", "Mu", "q5", 
                  "q25", "q75", "q95"), cumsum, .names = "{col}"))%>%ungroup()%>%
  arrange(scenario)%>%group_by(year, scenario)%>%
  summarise(across(c(par_col),~ sum(.x, na.rm = FALSE)))
x_catcost_total_incre <- list()
x <- list()
for(i in unique(x_catcost_total$scenario)){ 
  x[[i]] <- x_catcost_total%>%filter(scenario == i)
  x_catcost_total_incre[[i]] <- 
    cbind(year = x[[i]]$year, scenario = x[[i]]$scenario, 
          as.data.frame(x[[i]][, c(par_col)] - x_total_ref[, c(par_col)]))%>%
    as.data.frame()
  
  x_catcost_total_incre[[i]][x_catcost_total_incre[[i]] == 0] <- NA
  x_catcost_total_incre[[i]] <- x_catcost_total_incre[[i]]%>%
    popResults_range(POC_AU, ., end_Y = 100)
  

}

x_total_incre <- x_catcost_total_incre%>%dplyr::bind_rows(., .id = "scenario")%>%
  mutate(scenario = factor(scenario, levels = c("Status quo", "dfList_NP_2023", 
                                                "dfList_NP_2024", "dfList_NPexp_A", 
                                                "dfList_NPexp_B", "dfList_NPexp_C",
                                                "dfList_NPexp_D"), 
                           labels = c("No National Program", "Past program", 
                                      "Current program", "Continued program", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "Expanded program")))%>%
  filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
  
incremental_cost <- ggplot(x_total_incre, aes(x = year, colour = scenario) ) + 
  geom_line(aes(x = year, y = best, colour = scenario,
                linetype = scenario), size = 1
            ) + 
  theme_Publication() + 
  scale_x_continuous(expand = c(0,0), limits = c(2021,2041), breaks = seq(2021, 2041, 1)) + 
  scale_y_continuous(limits = c(-300000000, 50000000), 
                     breaks = seq(-300000000, 50000000, 10000000), 
                     labels = seq(-300000000, 50000000, 10000000)/1000000) + 
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) + 
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid", "solid")) + 
  labs( y = "Incremental cost (in millions)", x = "Year")

incremental_cost <- incremental_cost + 
  geom_hline(linetype = "dashed", yintercept = 0, size = 1)
ggsave(file=file.path(OutputFig, paste0("Reports/incremental_cost ",".png")), 
       incremental_cost , 
       width = 16, height = 8, bg = "white", dpi = 300)   




#### Cost saving 5, 10, 20 y
# total 
p_cost_saving <- list()
col_pal <- list("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
names(col_pal) <- c(unique(x_total_incre$scenario))

for(i in unique(x_total_incre$scenario)){ 
  
  p_cost_saving[[i]] <- 
    ggplot(x_total_incre%>%filter(scenario == i & year %in% (year_obs)), 
           aes(x = as.character(year), y = best)) + 
    geom_bar(stat="identity") + 
    geom_text(aes(label=paste0(round(best/1000000, digits = 1), "m")),
              position = position_dodge(width = 0.55),
              hjust = 0.5, vjust = -0.5) + 
    scale_x_discrete(labels = c(paste0((year_obs - POC_AU$simY + 1), "-Year", sep = ""))) + 
    labs(y = "Incremental costs (millions)", x = "Time frame")  + 
    scale_y_continuous(limits = c(-300000000, 100000000), 
                       breaks = seq(-300000000, 100000000, 10000000), 
                       labels = seq(-300000000, 100000000, 10000000)/1000000) + 
    theme_Publication() + 
    theme(panel.grid.major = element_line(color = "gray80",
                                    size = 0.1,
                                    linetype = 1)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle(i)
}


for(i in names(p_cost_saving)){ 
  ggsave(file=file.path(OutputFig, paste0("Reports/Tot_cost_saving", i,".png")), 
         p_cost_saving[[i]] , 
         width = 10, height = 8, bg = "white", dpi = 300)   
  }

# by categories 
catcost_cum <- cost_disydaacap_categories%>%
  filter(year>= 2022)%>%
  group_by(scenario, Categories)%>%
  mutate(across(c(par_col), cumsum, .names = "{col}"))%>%ungroup()%>%
  arrange(scenario)

catcost_cum_ref <- catcost_cum%>%filter(scenario == "Status quo")
x <- list()
catcost_incre <- list()
for(i in unique(catcost_cum$scenario)){ 
  x[[i]] <- catcost_cum%>%filter(scenario == i)
  catcost_incre[[i]] <- cbind(scenario = x[[i]]$scenario, 
                              Categories = x[[i]]$Categories, 
                              year = x[[i]]$year, 
                              dplyr::bind_cols(x[[i]][, par_col] - catcost_cum_ref[, par_col]))%>%
    as.data.frame()%>%
    popResults_range(POC_AU, .)%>%
    filter(Categories != "Treatment") 
  }

catcost_incre <- dplyr::bind_rows(catcost_incre, .id = "Scenario")
catcost_incre <- catcost_incre%>%
  mutate(Scenario = factor(Scenario, levels = c("Status quo", "dfList_NP_2023", 
                                                "dfList_NP_2024", "dfList_NPexp_A", 
                                                "dfList_NPexp_B", "dfList_NPexp_C",
                                                "dfList_NPexp_D"), 
                           labels = c("No National Program", "Past program", 
                                      "Current program", "Continued program", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "Expanded program")))%>%
  filter(!Scenario %in% c("NP expand 2025", "NP expand 2026"))
catcost_incre <- catcost_incre%>%
  mutate(Categories = factor(Categories, 
                             levels = c("Diagnosis", "Management", "Treatment_cap"),
                             labels = c("Diagnosis", "Management", "Treatment")))

p_catcost_saving <- list()
lim_catcost_saving <- list(c(-20000000, 50000000 ), c(-80000000, 30000000), c(-250000000, 20000000))
bek_catcost_saving <- list(seq(-20000000, 50000000 , 5000000), 
                           seq(-80000000, 30000000, 5000000), 
                           seq(-250000000, 20000000, 10000000))
for(i in 1: length(year_obs)){ 
  
  p_catcost_saving[[i]] <- 
    ggplot(catcost_incre%>%
             filter(Scenario != "No National Program" & year %in% year_obs[i])%>%
             arrange(Categories), 
           aes(x = Scenario, y = best, fill = Scenario)) + 
    geom_bar(stat = "identity", width = 0.8) + 
    geom_text(aes(y = best + 2 * sign(best), 
                  label=paste0(round(best/1000000, digits = 1), "m")),
              position = position_stack(vjust = 0.5), size = 6) + 
    scale_fill_manual(values = c(col_pal[2:5])) + 
    facet_wrap(~Categories, nrow = 3) + 
    theme_Publication_facet() + 
    
    labs(y = "Incremental costs (millions)", x = "Scenarios")  + 
    scale_y_continuous(limits = lim_catcost_saving[[i]], 
                       breaks = bek_catcost_saving[[i]], 
                       labels = bek_catcost_saving[[i]]/1000000) + 
    theme_Publication_facet() + 
    theme(panel.grid.major = element_line(color = "gray80",
                                          size = 0.1,
                                          linetype = 1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    ggtitle(paste0(POC_AU$simY, "-", year_obs[i], 
                   " (",year_obs[i] - POC_AU$simY + 1 ,"-Year)"))
}


for(i in 1: length(year_obs)){ 
  ggsave(file=file.path(OutputFig, paste0("Reports/cat_cost_saving", i,".png")), 
         p_catcost_saving[[i]] , 
         width = 8, height = 6, bg = "white", dpi = 300)   
}
p_cat_cost_supple <- p_catcost_saving[[3]] + facet_custom (~Categories,
                                      scales = "free", ncol = 1,
                                      scale_overrides = 
                                        list(
                                          scale_new(1,
                                                    scale_y_continuous(limits = 
                                                                         c(-5000000, 15000000),
                                                                       breaks = seq(-5000000, 15000000, 5000000),
                                                                       labels = seq(-5000000, 15000000, 5000000)/1000000)),
                                          scale_new(2,
                                                    scale_y_continuous(limits = 
                                                                         c(-20000000, 0),
                                                                       breaks = seq(-20000000, 0, 5000000),
                                                                       labels = seq(-20000000, 0, 5000000)/1000000)),
                                          
                                          scale_new(3,
                                                    scale_y_continuous(limits = 
                                                                         c(-300000000, 0 ),
                                                                       breaks = seq(-300000000, 0, 50000000),
                                                                       labels = seq(-300000000, 0, 50000000)/1000000))
                                        )) + 
  labs(caption = "*HCV management cost in the expanded National Program scenario: A$35,818")+
  theme_Publication_facet(base_size = 22) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) + 
  theme(legend.direction = "vertical") + 
  theme(plot.caption.position = "plot",
        plot.caption = element_text(hjust = 1, size = 20))
  
  

ggsave(file=file.path(OutputFig, paste0("Reports/cat_cost_saving","cat_cost_supple",".png")), 
       p_cat_cost_supple , 
       width = 12, height = 14, bg = "white", dpi = 300)    

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

legend <- g_legend(p_cat_cost_supple) 

grid.newpage()
leg_x <- grid.draw(legend) 
ggsave(file=file.path(OutputFig, paste0("Reports/cat_cost_saving","legend",".png")), 
       legend  , 
       width = 6, height = 4, bg = "white", dpi = 300) 
################################################################################

####                            Tables generation                           ####

################################################################################

# generating tables 
# Resflow all indicators 
# DC, HCC, LT, PLT 
# best (95% PI) 
# substract years 
# making functions 
# for senstivity analysis: 
# load(file.path(OutputFolder, paste0(project_name, "Res_flowcost_discount_0.03.rda")))
# load(file.path(OutputFolder, paste0(project_name, "Res_flowcost_discount_0.07.rda")))

load(file.path(OutputFolder, paste0(project_name, "Res_flowcost_discount_0.07.rda")))
# test  
tab_epi <- Resflow_all_lst
cost_qaly_range <- list() 
cost_qaly_range_disy <- list()
cost_qaly_range_ycum <- list()
cost_qaly_range_disycum <- list()
for(i in names(Rescost_year_all)){ 
  for(n in names(Rescost_year_all[[1]])){ 
    Rescost_year_all[[i]][[n]][Rescost_year_all[[i]][[n]] ==0 ]  <- NA
    
    Rescost_disyear_all[[i]][[n]][Rescost_disyear_all[[i]][[n]] ==0 ]  <- NA
   
    cost_qaly_range[[i]][[n]] <- 
      popResults_range(POC_AU, Rescost_year_all[[i]][[n]], end_Y = 100-1)%>%
      as_tibble()
    
    cost_qaly_range_disy[[i]][[n]] <- 
      popResults_range(POC_AU, Rescost_disyear_all[[i]][[n]], end_Y = 100-1)%>%
      as_tibble()
    
    cost_qaly_range_ycum[[i]][[n]] <- Rescost_year_all[[i]][[n]]%>%
      filter(year>=2022)%>%
      mutate(across(par_col, list(cum=cumsum), .names = "{col}"))%>%
      popResults_range(POC_AU, ., end_Y = 100-1)%>%as_tibble()
    
    cost_qaly_range_disycum[[i]][[n]] <- Rescost_disyear_all[[i]][[n]]%>%
      filter(year>=2022)%>%
      mutate(across(par_col, list(cum=cumsum), .names = "{col}"))%>%
               popResults_range(POC_AU, ., end_Y = 100-1)%>%
               as_tibble()
            
    
  }
}

tab_costqaly <- cost_qaly_range%>%purrr::transpose()


tab_costqaly_disy <- cost_qaly_range_disy%>%purrr::transpose()
tab_costqaly_ycum <- cost_qaly_range_ycum%>%purrr::transpose()
tab_costqaly_disycum <- cost_qaly_range_disycum%>%purrr::transpose() 

# function to factors scenarios and bind list 
bindfac <- function(dt){ 
  x <- lapply(dt, function(x){ 
    x%>%dplyr::bind_rows(., .id = "scenario")%>%
      mutate(scenario = factor(scenario, 
                               levels = c("Status quo", "dfList_NP_2023", 
                                          "dfList_NP_2024", "dfList_NPexp_A", 
                                          "dfList_NPexp_B", "dfList_NPexp_C",
                                          "dfList_NPexp_D"), 
                               labels = c("No National Program", "Past program", 
                                          "Current program", "Continued program", 
                                          "NP expand 2025", "NP expand 2026", 
                                          "Expanded program")))
  } )
  
  return(x)
  }

tab_costqaly <- bindfac(tab_costqaly)
tab_costqaly_disy <- bindfac(tab_costqaly_disy)
tab_costqaly_ycum <- bindfac(tab_costqaly_ycum)
tab_costqaly_disycum <- bindfac(tab_costqaly_disycum)

tab_costqaly_lst <- list("year" = tab_costqaly, 
                         "disy" = tab_costqaly_disy, 
                         "ycum" = tab_costqaly_ycum, 
                         "disycum" = tab_costqaly_disycum)
#### table generation ####  

library("data.table")
library("formattable")
library("gt")
library("writexl")
library("gtsummary")
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



# grouping rows
tab_x <- list()
for(i in names(tab_epi)){ 
  tab_x[[i]] <- tab_epi[[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
    filter(year %in% seq(2022,2041,1))%>%
    select(Indicators, scenario, year, best, q5, q95)%>%
    mutate(best = formatC(best,  format = "fg", big.mark = ","),
           q5 = formatC(q5,  format = "fg", big.mark = ","),
           q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
    mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
    select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
  
  tab_x[[i]]%>%group_by(Indicators)%>%
    gt(groupname_col = "Indicators",
       rowname_col = "year")%>%
    gtsave(., file = file.path(OutputFig, paste0("Reports/", i, ".docx")))
}

tab_cost <- list()

for(i in names(tab_costqaly_lst)){ 
  
  tab_cost[[i]] <- tab_costqaly_lst[[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
    filter(year %in% seq(2022,2041,1))%>%
    select(Indicators, scenario, year, best, q5, q95)%>%
    mutate(best = formatC(best,  format = "fg", big.mark = ","),
           q5 = formatC(q5,  format = "fg", big.mark = ","),
           q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
    mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
    select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
  
  tab_cost[[i]]%>%group_by(Indicators)%>%
    gt(groupname_col = "Indicators",
       rowname_col = "year")%>%
    gtsave(., file = file.path(OutputFig, paste0("Reports/costqaly_", i, ".docx")))
  
  }

# numbox
tab_numbox = list("year" = Res_Numbox_y_range, 
                  "cumyear" = Res_Numbox_cum_range, 
                  "cumavert" = Res_Numbox_avert_range)

for(i in names(tab_numbox)){ 
  tab_numbox[[i]] <- tab_numbox[[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
    filter(year %in% seq(2022,2041,1))%>%
    select(Indicators, scenario, year, best, q5, q95)%>%
    mutate(best = formatC(best,  format = "fg", big.mark = ","),
           q5 = formatC(q5,  format = "fg", big.mark = ","),
           q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
    mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
    select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
  
  tab_numbox[[i]]%>%group_by(Indicators)%>%
    gt(groupname_col = "Indicators",
       rowname_col = "year")%>%
    gtsave(., file = file.path(OutputFig, paste0("Reports/numbox_", i, ".docx")))
}


################################################################################

####                                CEA                                     ####
# c("QALY_compartment", "cost_total_Cap", "cost_total")




CEAanalysis <- list()
timeframe <- c(5, 10, 20, 30, 40, 50, 60)
timeframe_name <- c("5y", "10y", "20y", "30y", "40y", "50y", "60y")
for(i in 1: length(timeframe)){ 
  CEAanalysis[[timeframe_name[i]]][["QALY"]] <- tab_costqaly_lst$disycum$QALY_compartment%>%
    filter(year == POC_AU$simY + timeframe[i] - 1)%>%
    split(.$scenario)%>%map(~.x %>% select(-scenario))%>%
    replace(is.na(.), 0)
  
  CEAanalysis[[timeframe_name[i]]][["Cost"]] <- tab_costqaly_lst$disycum$cost_total%>%
    filter(year == POC_AU$simY + timeframe[i] - 1)%>%
    split(.$scenario)%>%map(~.x %>% select(-scenario))%>%
    replace(is.na(.), 0)
  
  CEAanalysis[[timeframe_name[i]]][["Cost_cap"]] <- tab_costqaly_lst$disycum$cost_total_Cap%>%
    filter(year == POC_AU$simY + timeframe[i] - 1)%>%
    split(.$scenario)%>%map(~.x %>% select(-scenario))%>%
    replace(is.na(.), 0)
  
  }

tab_costqaly_lst$disycum$cost_total[, c(1:6)]
Incre <- list()
for(m in names(CEAanalysis)){
  for(i in names(CEAanalysis[[1]])){ 
    for(n in names(CEAanalysis[[1]][[1]])){ 
      Incre[[m]][[i]][[n]] <- cbind(year = CEAanalysis[[m]][[i]][[n]]$year, 
                                    dplyr::bind_cols(CEAanalysis[[m]][[i]][[n]][, par_col] - 
                                                    CEAanalysis[[m]][[i]][["No National Program"]][, par_col]))%>%
        as_tibble()%>%
        popResults_range(POC_AU, .)
    }
  }
  
}


Incre<- lapply(Incre, function(x) x%>%purrr::transpose())


CEA <- list()
CEA_cap <- list()
for(m in names(Incre)){
  for(i in names(Incre[[1]])){ 
    
    CEA[[m]][[i]] <- 
      cbind(year = Incre[[m]][[i]][["Cost"]]$year,
            dplyr::bind_cols(Incre[[m]][[i]][["Cost"]][, par_col]/Incre[[m]][[i]][["QALY"]][, par_col]))
    
    
    CEA_cap[[m]][[i]] <- cbind(year = Incre[[m]][[i]][["Cost_cap"]]$year,
                               dplyr::bind_cols(Incre[[m]][[i]][["Cost_cap"]][, par_col]/Incre[[m]][[i]][["QALY"]][, par_col]))
  }
  
  }

for(i in names(CEA)){ 
  for(n in names(CEA[[1]])){ 
    CEA[[i]][[n]] <- CEA[[i]][[n]]%>%as_tibble()%>%popResults_range(POC_AU, .)
    CEA_cap[[i]][[n]] <- CEA_cap[[i]][[n]]%>%as_tibble()%>%popResults_range(POC_AU, .)
    }
  }

te <- lapply(Incre, function(x) lapply(x, function(y) dplyr::bind_rows(y, .id = "indicator")))


te <- lapply(te, function(x) x%>%dplyr::bind_rows(., .id = "scenario"))


x_cap <- lapply(te, function(x) x%>%filter(indicator != "Cost"))
x_nocap <- lapply(te, function(x) x%>%filter(indicator != "Cost_cap"))
x_cap <- lapply(x_cap, function(x) x[, c(1:1004)])

x <- CEAanalysis$`20y`$QALY$`No National Program`
x_sce <- list()
for(i in names(CEAanalysis$`20y`$QALY)){ 
  x_sce[[i]] <- cbind(year = x$year, 
                      as.data.frame(CEAanalysis$`20y`$QALY[[i]][, par_col] - x[, par_col]))
  
  x_sce[[i]] <- x_sce[[i]]%>%popResults_range(POC_AU, .)
  
  x_sce[[i]] <- x_sce[[i]]%>%select(year, best, q5, q95)
  }
x_sce$`Expanded program achievement`

x_nocap[["20y"]]
for(i in names(x_cap)){ 
  x_cap[[i]] <- x_cap[[i]]%>%gather(sim, val, -c(indicator, scenario, year))%>%
    spread(indicator, val)
}
  
View(x_cap[["20y"]])
PSA_dt <- x_cap[["20y"]]%>%
  filter(!scenario %in% c("No National Program", "NP expand 2025", "NP expand 2026") & Cost_cap != 0)%>%
  mutate(scenario = factor(scenario, level = c("Past program", 
                                               "Current program", "Continued program", 
                                               "Expanded program")))


PSA_dt<- PSA_dt%>%mutate(outline = ifelse(sim == "best", 1,0))%>%
  mutate(category = paste0(scenario, outline))%>%
  mutate(category = factor(category, level = c("Past program0", 
                           "Current program0", "Continued program0", 
                           "Expanded program0", 
                           "Past program1", 
                           "Current program1", "Continued program1", 
                           "Expanded program1"))) 

PSA <- ggplot(PSA_dt, 
              aes(y = `Cost_cap`, x = QALY)) + 
  geom_point(aes(colour = scenario))  +
  facet_wrap(~scenario) +
  scale_color_manual(name = "Scenarios", 
                     values = c("#E69F00", "#56B4E9","#009E73", "#F0E442")) +
  geom_point(data = PSA_dt%>%filter(outline == 1), color = "gray50")  +
  facet_wrap(~scenario) +
  
  geom_hline(yintercept=0, color = "black", linewidth = 1) +
  geom_vline(xintercept=0, color = "black", linewidth = 1)   +
  labs(colour = "Scenarios",  x = "QALY", 
       y = "Costs (DAA cap, millions)") + 
  theme_bw() + 
    scale_y_continuous(limits = c(-500000000, 200000000),
                       breaks = seq(-500000000, 200000000, 100000000), 
                       labels = seq(-500000000, 200000000, 100000000)/1000000) + 
  scale_x_continuous(limits = c(0, 3000)) + 
  theme(plot.title = element_text(hjust = 0.5)) +  
  theme(panel.background = element_rect(colour = "white"),
        plot.background = element_rect(colour = "white"),
        panel.border     = element_rect(fill = NA, colour = "black", 
                                        size = NA),
        plot.title = element_text(face = "bold",
                                  size = 12, hjust = 0.5),
        text = element_text(),
        axis.title = element_text(face = "bold",size = 12),
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
        legend.direction = "vertical",
        legend.title = element_text(face="bold", size= 14),
        plot.margin = unit(c(10,5,5,5),"mm")) + 
  guides(color =guide_legend(direction='vertical',
                             override.aes = list(size=2)),
         fill = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_abline(aes(slope = 50000, intercept = 0,linetype = "WTP: A50,000"), colour = "black") + 
  scale_linetype_manual(name = "", values = c(2), 
                        guide = guide_legend(override.aes = list(color = c("black")))) +
  stat_ellipse(color = "gray50", 
               alpha = 0.7,
               linewidth = 0.6,
               show.legend = FALSE, 
               level = 0.95) 
PSA


ggsave(file=file.path(OutputFig, paste0("Reports/PSA.png")), 
       PSA , 
       width = 12, height = 8, bg = "white", dpi = 300)  


# CEA table 

CEA_cap <- lapply(CEA_cap, function(x) x%>%dplyr::bind_rows(., .id = "Scenario"))
CEA_cap <- CEA_cap%>%dplyr::bind_rows(., .id = "Timeframe")
CEA_cap <- CEA_cap%>%mutate(Timframe = factor(Timeframe, 
                                              levels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y"), 
                                              labels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y")))
CEA_cap <- CEA_cap%>%select(Timeframe, Scenario,best, q5, q95)%>%
  filter(!Scenario %in% c("No National Program", "NP expand 2025", "NP expand 2026"))%>%
  mutate(best = formatC(best,  format = "fg", big.mark = ","),
         q5 = formatC(q5,  format = "fg", big.mark = ","),
         q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
  mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
  select(-c(best, q5, q95))%>%ungroup()%>%spread(Scenario, vv)

CEA_cap%>%
  gt(
     rowname_col = "Timeframe")%>%
  gtsave(., file = file.path(OutputFig, paste0("Reports/costqaly_CEA.docx")))


# one-way sensitivity analysis 
CEA_nocap <- lapply(CEA, function(x) x%>%dplyr::bind_rows(., .id = "Scenario"))
CEA_nocap <- CEA_nocap%>%dplyr::bind_rows(., .id = "Timeframe")
CEA_nocap <- CEA_nocap%>%mutate(Timframe = factor(Timeframe, 
                                              levels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y"), 
                                              labels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y")))
CEA_nocap <- CEA_nocap%>%select(Timeframe, Scenario,best, q5, q95)%>%
  filter(!Scenario %in% c("No National Program", "NP expand 2025", "NP expand 2026"))%>%
  mutate(best = formatC(best,  format = "fg", big.mark = ","),
         q5 = formatC(q5,  format = "fg", big.mark = ","),
         q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
  mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
  select(-c(best, q5, q95))%>%ungroup()%>%spread(Scenario, vv)

CEA_nocap%>%
  gt(
    rowname_col = "Timeframe")%>%
  gtsave(., file = file.path(OutputFig, paste0("Reports/costqaly_CEA_nocap.docx")))

# # discount rate: 3%, 7% 
# # CEA_cap_discout_lower <- lapply(CEA_cap, function(x) x%>%dplyr::bind_rows(., .id = "Scenario"))
# # CEA_cap_discout_lower<- CEA_cap_discout_lower%>%dplyr::bind_rows(., .id = "Timeframe")
# # CEA_cap_discout_lower <- CEA_cap_discout_lower%>%mutate(Timframe = factor(Timeframe, 
# #                                               levels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y"), 
# #                                               labels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y")))
# # CEA_cap_discout_lower <- CEA_cap_discout_lower%>%select(Timeframe, Scenario,best, q5, q95)%>%
# #   filter(!Scenario %in% c("No National Program", "NP expand 2025", "NP expand 2026"))%>%
# #   mutate(best = formatC(best,  format = "fg", big.mark = ","),
# #          q5 = formatC(q5,  format = "fg", big.mark = ","),
# #          q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
# #   mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
# #   select(-c(best, q5, q95))%>%ungroup()%>%spread(Scenario, vv)
# # 
# # CEA_cap_discout_lower%>%
# #   gt(
# #     rowname_col = "Timeframe")%>%
# #   gtsave(., file = file.path(OutputFig, paste0("Reports/costqaly_CEA_discount_.docx")))
# 
# # 7% 
# # discount rate: 3%, 7% 
# CEA_cap_discout_upper <- lapply(CEA_cap, function(x) x%>%dplyr::bind_rows(., .id = "Scenario"))
# CEA_cap_discout_upper <- CEA_cap_discout_upper%>%dplyr::bind_rows(., .id = "Timeframe")
# CEA_cap_discout_upper <- CEA_cap_discout_upper%>%mutate(Timframe = factor(Timeframe, 
#                                                                           levels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y"), 
#                                                                           labels = c("5y", "10y", "20y", "30y", "40y", "50y", "60y")))
# CEA_cap_discout_upper <- CEA_cap_discout_upper%>%select(Timeframe, Scenario,best, q5, q95)%>%
#   filter(!Scenario %in% c("No National Program", "NP expand 2025", "NP expand 2026"))%>%
#   mutate(best = formatC(best,  format = "fg", big.mark = ","),
#          q5 = formatC(q5,  format = "fg", big.mark = ","),
#          q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
#   mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
#   select(-c(best, q5, q95))%>%ungroup()%>%spread(Scenario, vv)
# 
# CEA_cap_discout_upper%>%
#   gt(
#     rowname_col = "Timeframe")%>%
#   gtsave(., file = file.path(OutputFig, paste0("Reports/costqaly_CEA_discountupper_.docx")))
# 
# 
# 
# # tornado plot for one-way sensitivity analysis 
# 
# # import the data 
# 
# 
# library(readxl)
# library(ggplot2)
# library(plyr)
# library(dplyr)
# library(tidyverse)
# library(here)
# 
# # this is throwing some warnings in my computer, but it is reading the data frame correctly
# df <-read_xlsx(path = file.path(OutputFolder, paste0("tornado_dt.xlsx")))%>%
#   as.data.frame()%>%mutate(Diff = abs(Upper - Lower))%>%
#   mutate(scenario = factor(scenario, levels = c("Past program achievement", 
#                                                 "Current program achievement",
#                                                 "Continued program achievement",
#                                                 "Expanded program achievement"),
#                            labels = c("Past program", 
#                                       "Cureent program",
#                                       "Continued program",
#                                       "Expanded program")))
# 
# tornado_dt <- list()
# base_value <- list()
# order.parameters <- list()
# tornado_dt_order <- list()
# # width of columns in plot (value between 0 and 1)
# width <- 0.45
# for(i in unique(df$scenario)){ 
#   tornado_dt[[i]] <- df%>%filter(scenario == i)
#   base_value[[i]] <- tornado_dt[[i]]%>%filter(Paramters == "DAA cap"  )%>%
#     select(Lower)%>%as.numeric() 
#   
#   order.parameters[[i]] <- tornado_dt[[i]] %>% arrange(Diff) %>%
#     mutate(Paramters=factor(x=Paramters, levels=Paramters)) %>%
#     select(Paramters) %>% unlist() %>% levels()
#   
#   # get data frame in shape for ggplot and geom_rect
#   tornado_dt_order[[i]] <- tornado_dt[[i]] %>% 
#     # gather columns Lower_Bound and Upper_Bound into a single column using gather
#     gather(key='type', value='output.value', Lower:Upper) %>%
#     # just reordering columns
#     select(Paramters, type, output.value, Diff) %>%
#     # create the columns for geom_rect
#     mutate(Paramters=factor(Paramters, levels=order.parameters[[i]]),
#            ymin=pmin(output.value, base_value[[i]]),
#            ymax=pmax(output.value, base_value[[i]]),
#            xmin=as.numeric(Paramters)-width/2,
#            xmax=as.numeric(Paramters)+width/2)
#   tornado_dt_order[[i]] <- tornado_dt_order[[i]][order(tornado_dt_order[[i]]$Diff), ]
#   
# }
# mylabel <- list()
# mylabel[[1]] <- c("Discounted rate\n(3% vs. 7%)",
#              "Timeframe\n(10-year vs. 30-year)",
#              "DAA cap\n(Capped vs. no cap)")
# mylabel[[2]] <- c("Discounted rate\n(3% vs. 7%)",
#                   "Timeframe\n(10-year vs. 30-year)",
#                   "DAA cap\n(Capped vs. no cap)")
# 
# mylabel[[3]] <- c("Discounted rate\n(3% vs. 7%)",
#                   "Timeframe\n(10-year vs. 30-year)",
#                   "DAA cap\n(Capped vs. no cap)")
# 
# mylabel[[4]] <- c("Discounted rate\n(3% vs. 7%)",
#                   "DAA cap\n(Capped vs. no cap)", 
#                   "Timeframe\n(10-year vs. 30-year)")
# names(mylabel) <- names(tornado_dt_order)
# # get order of parameters according to size of intervals
# # (I use this to define the ordering of the factors which I then use to define the positions in the plot)
# p_tornado <- list()
# for(i in names(tornado_dt_order)){ 
#   p_tornado[[i]] <-  
#   ggplot() + 
#     geom_rect(data = tornado_dt_order[[i]], 
#               aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill=type)) +
#     theme_classic() + 
#     theme(axis.title.y = element_blank(), legend.position = 'bottom',
#           legend.title = element_blank(),
#           axis.line.y =element_blank(),
#           axis.ticks.y = element_blank()
#     ) + 
#     
#     theme(axis.text.x = element_text(colour="Black",
#                                      face="bold", 
#                                      size=14,
#                                      angle = 90,
#                                      hjust = 0.5,
#                                      vjust = 0.5))+
#     theme(axis.text.y = element_text(colour="Black",
#                                      face="bold", 
#                                      size=14, 
#                                      angle = 0,
#                                      hjust = 1,
#                                      vjust = 0.5))+
#     theme(axis.title.x = element_text(colour="Black",
#                                      face="bold", 
#                                      size=14, 
#                                      angle = 0,
#                                      hjust = 0.5)) +
#     theme(legend.text =  element_text(colour="Black",
#                                       face="bold", 
#                                       size=14, 
#                                       angle = 360,
#                                       hjust = 1,
#                                       vjust = 0.5))+
#     theme(plot.title = element_text(hjust = 0.5, colour="Black",
#                                     face="bold", 
#                                     size=14)) + 
#     geom_hline(yintercept=base_value[[i]],colour="black") + 
#     geom_hline(yintercept=0,colour="gray", linetype= "dashed") +
#     scale_x_continuous(breaks = c(1:length(order.parameters[[i]])), 
#                        labels = mylabel[[i]])+ 
#     
#     scale_y_continuous(limits = c(-140000,0), 
#                        breaks=c(-140000, base_value[[i]], seq(-100000 ,0, 20000))) +
#     scale_fill_manual(labels = c("Lower estimate", "Upper estimate"), 
#                       values = c("#779b9b", "#2f4f4e"))+
#     ylab("ICERs") + 
#    
#     theme(plot.background = element_rect(fill = "transparent",colour = NA))+
#     coord_flip() + 
#     ggtitle(i) 
#     
#   }
# 
# 
# p_tornado[[4]] <- p_tornado[[4]] + 
#   scale_y_continuous(limits = c(-140000,0), 
#                      breaks=c(-140000,-120000, base_value[[4]], seq(-80000 ,0, 20000)))
# 
# p_tornado_supple  <- ggarrange(plotlist = p_tornado, ncol = 2, nrow = 2, 
#                                   common.legend = TRUE)
# 
# ggsave(file=file.path(OutputFig_y_cum_avert, paste0("tornado_supple", ".png")), 
#        p_tornado_supple , 
#        width = 12, height = 10, bg = "white", dpi = 300)
# # create plot
# # (use scale_x_continuous to change labels in y axis to name of parameters)
# ## label
# 
# 
# 
# 
# ggplot() + 
#   geom_rect(data = df.2, 
#             aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill=type)) +
#   theme_classic() + 
#   theme(axis.title.y = element_blank(), legend.position = 'bottom',
#         legend.title = element_blank(),
#         axis.line.y =element_blank(),
#         axis.ticks.y = element_blank()
#   ) + 
#   
#   theme(axis.text.x = element_text(colour="Black",
#                                    face="bold", 
#                                    size=10,
#                                    angle = 0,
#                                    hjust = 0.5,
#                                    vjust = 1))+
#   theme(axis.text.y = element_text(colour="Black",
#                                    face="bold", 
#                                    size=10, 
#                                    angle = 360,
#                                    hjust = 1,
#                                    vjust = 0.5))+
#   theme(legend.text =  element_text(colour="Black",
#                                     face="bold", 
#                                     size=10, 
#                                     angle = 360,
#                                     hjust = 1,
#                                     vjust = 0.5))+
#   geom_hline(yintercept=base_value,colour="black") +
#   scale_x_continuous(breaks = c(1:length(order.parameters)), 
#                      labels = mylabel)+
#   scale_y_continuous(expand = c(0, 0),limits = c(0,21), 
#                      breaks=c(0,5,base_value,10,15,20)) +
#   scale_fill_manual(labels = c("Lower estimate", "Upper estimate"), 
#                     values = c("#779b9b", "#2f4f4e"))+
#   ylab("Ratio") + 
#   
#   
#   
#   annotate("text", x=9, y=18, label= "17.4", 
#            colour="Black",fontface="bold",  size=labsize) + 
#   annotate("text", x=8, y=5.9, label= "6.4",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=8, y=10.4, label= "9.9",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=7, y=3.3, label= "3.8",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=6, y=9, label= "8.5",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=6, y=5, label= "5.5",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=5, y=8.4, label= "7.9",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=5, y=4.7, label= "5.2",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=4, y=5.4, label= "5.9",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=4, y=8.2, label= "7.7",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=3, y=8, label= "7.5",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=3, y=6.2, label= "6.7",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=2, y=6.5, label= "6.9",
#            colour="Black", fontface="bold",  size=labsize) +
#   annotate("text", x=1, y=6.5, label= "7.0",
#            colour="Black", fontface="bold",  size=labsize) +
#   theme(plot.background = element_rect(fill = "transparent",colour = NA))+
#   coord_flip()
# 
# 
# ggsave(path = '/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/13. Other research project/PrEP modelling', 
#        filename="CEA_Sensitivity.png", 
#        width = 12, height = 8, device='png', dpi=900)
# 
# 
# 
# 
# ################################################################################
# #                                 USED                              
# ########################## plot animation ######################################
# install.packages("gganimate")
# library("gganimate")
# x_catcost_total <- x_catcost%>%group_by(year, scenario)%>%
#   summarise(across(c("best", "q5", "q95"),~ sum(.x, na.rm = FALSE)))%>%
#   mutate(year = as.integer(year))%>%filter(year <= 2045)
# View(x_catcost_total)
# library(forcats)
# # creating pause timeline 
# df2 <- x_catcost_total %>%
#   group_by(year) %>%
#   # The * 1 makes it possible to have non-integer ranks while sliding
#   mutate(ordering = min_rank(-best) * 1) %>%
#   ungroup()%>%
#   mutate(text_lab = paste0(format(round(best/1000000000, digits = 3), nsmall = 3), "b"),
#          text_y = best)%>%
#   mutate(text_lab =as.factor(text_lab))%>% 
#   mutate(show_time = case_when(year %in% c(2026,2028,2030, 2031, 2032) ~ 10,
#                                TRUE           ~ 1)) %>%
#   # uncount is a tidyr function which copies each line 'n' times
#   uncount(show_time) %>%
#   group_by(year, scenario) %>%
#   mutate(reveal_time = row_number()) %>%
#   ungroup()
# 
# p <- ggplot(df2, aes(ordering, group = scenario, 
#                           fill = as.factor(scenario), color = as.factor(scenario))) +
#   geom_tile(aes(y = best/2,
#                 height = best,
#                 width = 1), alpha = 1, color = NA) +
#   # text in x-axis (requires clip = "off" in coord_*)
#   # paste(country, " ")  is a hack to make pretty spacing, since hjust > 1 
#   #   leads to weird artifacts in text spacing.
#   geom_text(aes(y = 0, label = paste(scenario, " ")), vjust = 0.2, hjust = 1, 
#             size = 4, fontface = "bold") +
#   geom_text(aes(y = text_y, label = text_lab, 
#                  group=text_y), color = "white", vjust = 0.2, hjust = 1, 
#             size = 8, fontface='bold') + 
# 
#   coord_flip(clip = "off", expand = FALSE) +
#   scale_y_continuous(labels = scales::comma) +
#   scale_x_reverse() + 
# 
# 
#   scale_fill_manual(labels = c("", "", "", "", ""),
#                     values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
#   scale_color_manual(labels = c("", "", "", "", ""),
#                     values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
#   theme_Publication() + 
#   theme(legend.position = "") + 
#   guides(color = FALSE, fill = FALSE) +
# 
#   labs(title='{closest_state}', x = "", y = "Cumulative discounted cost") +
#   theme(plot.title = element_text(hjust = 0, size = 22),
#         axis.ticks.y = element_blank(),  # These relate to the axes post-flip
#         axis.text.y  = element_blank(),  # These relate to the axes post-flip
#         plot.margin = margin(1,1,1,5, "cm")) +
#   
#   transition_states(year, 
#                     transition_length = c(rep(1,4),3,1,3,3,3,3,3, 
#                                           rep(1, 12), 12), state_length = 1) +
#   ease_aes('cubic-in-out') 
# 
# anim <- animate(p, fps = 25, duration = 50, width = 800, height = 600)
# magick::image_write(anim, path= paste0(OutputFig, "totalcost_turningyear.gif"))
# # make different timeframe animation then compress together 
# anim
# 
# 
# 
# p <- ggplot(x_catcost_total,
#             aes(x=scenario, y=best, label = scenario, color = scenario)) + 
#   geom_point(stat = "identity", size = 3) + 
#   geom_segment(aes(
#     y = 5,
#     x = scenario,
#     yend = best,
#     xend = scenario)
#   ) + 
#   coord_flip() + 
#   theme_Publication() + 
#   scale_y_continuous(limits = c(0, 4500000000), breaks = seq(0,4500000000, 500000000), 
#                      labels = seq(0,4500000000, 500000000)/1000000000) + 
#   labs(y = "Cumulative discounted cost", x = "Scenarios")
#   transition_states(year, wrap = FALSE) 
# 
# 
# animate(p)
# 
# 
# 
# ################################################################################
# 
#                          # Used #
# 
# ################################################################################
# 
# 
# 
# 
# for(i in names(Rescost_yearcum_all)){ 
#   
#   cost_year_all[[i]] <- Rescost_yearcum_all[[i]]%>%group_by(year)%>%slice(n())%>%
#     mutate(cost_totalDAAcap_ycum = ifelse(cost_totalDAA_ycum>=cap, cap, cost_totalDAA_ycum))%>%
#     mutate(cost_totalDAAcap_disycum = ifelse(cost_totalDAAcap_ycum<cap, cost_totalDAAcap_ycum/discount, cap))%>%
#     mutate(cost_totalcap_ycum = cost_compartment_ycum + cost_ab_ycum + 
#              cost_RNA_ycum +  cost_POCT_ycum + cost_totalDAA_ycum + 
#              cost_TreatOther_ycum + cost_RetreatOther_ycum + cost_Cured_ycum)%>%
#     mutate(cost_totalcap_disycum = cost_compartment_disycum + cost_ab_disycum + 
#              cost_RNA_disycum +  cost_POCT_disycum + cost_totalDAA_disycum + 
#              cost_TreatOther_disycum + cost_RetreatOther_disycum + cost_Cured_disycum)%>%
#     select(year, discount, matches('_ycum'), cost_totalcap_ycum, 
#            cost_totalcap_disycum,cost_totalDAAcap_ycum, cost_totalDAAcap_disycum, 
#            matches('_disycum'))%>%
#     filter(year>= POC_AU$simY)
# 
# }
# 
# col_lst <- colnames(cost_year_all[[1]])[-c(1,2)]
# cost_cum_all <- list()
# cost_cumavert_all <- list()
# for(i in names(cost_year_all)){ 
#   cost_cum_all[[i]] <- cost_year_all[[i]]%>%ungroup()%>%
#     mutate_at(col_lst, "cumsum")%>%
#     select(year, discount ,col_lst)
#   
#   cost_cumavert_all[[i]] <- cbind(year =  cost_cum_all[[i]]$year, 
#                                   discount = cost_cum_all[[i]]$year, 
#                                   data.frame(cost_cum_all[[1]][, col_lst] - 
#                                             cost_cum_all[[i]][, col_lst]))%>%
#   rename_at(vars(col_lst), ~ paste0("avert_", col_lst))
# }
# 
# cost_cum_all_bind <- dplyr::bind_rows(cost_cum_all, .id ='scenario')%>%
#   mutate(scenario = factor(scenario, levels = c(names(cost_cum_all))))
# 
# long.cost_cum_all_bind <- gather(cost_cum_all_bind, index, best, 
#                                  -c(scenario, year, discount))
# 
# cost_cumavert_all_bind <- dplyr::bind_rows(cost_cumavert_all, .id ='scenario')%>%
#   mutate(scenario = factor(scenario, levels = c(names(cost_cum_all))))
# 
# long.cost_cumavert_all_bind <- gather(cost_cumavert_all_bind, index, best, 
#                                  -c(scenario, year, discount))
# 
# 
# plot_cost_cum <- list()
# plot_cost_cumavert <- list()
# # extracting cumsum releated index
# col_cum <- unique(long.cost_cum_all_bind $index)
# 
# ceil_avert <- list()
# ceil_cum <- list()
# floor <- list()
# long.cost_cum_all_bind <- long.cost_cum_all_bind%>%mutate(best = best/1000000000)
# for(i in col_cum){ 
#   
#   ceil_cum[[i]] <- long.cost_cum_all_bind%>%
#     filter(index == i & year%in% year_obs)%>%select(best)%>%
#     c()%>%unlist()%>%max()%>%ceiling_dec(.,-1)
#   floor[[i]] <- long.cost_cum_all_bind%>%filter(index == i)%>%select(best)%>%
#     c()%>%unlist()%>%min()%>%floor_dec(.,-1)
#   
#   if(floor[[i]]>=0){ 
#     plot_cost_cum[[i]] <- plot_pocau(POC_AU, long.cost_cum_all_bind, type = "cum", 
#                                      indicator = i, year_obs = year_obs) + 
#       labs(y = paste0(i, "(billions)")) + 
#       ggtitle("Cumulative cost") + 
#       scale_y_continuous(limits = c(0, ceil_cum[[i]]), breaks = seq(0, ceil_cum[[i]], 
#                                                                     ceil_cum[[i]]/10))
#   }
#   else{ 
#     
#     plot_cost_cum[[i]] <- plot_pocau(POC_AU, long.cost_cum_all_bind, type = "cum", 
#                                      indicator = i, year_obs = year_obs) + 
#       labs(y = paste0(i, "(billions)")) + 
#       ggtitle("Cumulative cost") + 
#       scale_y_continuous(limits = c(floor[[i]], ceil_cum[[i]]), breaks = seq(floor[[i]], ceil_cum[[i]], 
#                                                                              (ceil_cum[[i]]- floor[[i]])/10))
#   }
#   
# }
# 
# plot_cost_cum$cost_total_ycum
# 
# col_cum_avert <- c(paste0("avert_", col_cum))
# floor_avert <- list()
# plot_cost_cumavert <- list() 
# 
# col_lst <- colnames(cost_year_all[[1]])[-c(1,2)]
# cost_cum_all_year_obs <- list()
# cost_cumavert_all <- list()
# for(i in names(cost_cum_all)){
#   
#   cost_cum_all_year_obs[[i]] <- cost_cum_all[[i]]%>%filter(year %in% year_obs)
#   cost_cumavert_all[[i]] <- cbind(year = cost_cum_all_year_obs[[i]]$year, 
#                                   discount = cost_cum_all_year_obs[[i]]$discount, 
#                                   data.frame(cost_cum_all_year_obs[[1]][, col_lst] - 
#                                                cost_cum_all_year_obs[[i]][, col_lst]))%>%
#     rename_at(vars(col_lst), ~ paste0("avert_", col_lst))
#   
#   
# }
# 
# cost_cumavert_all_bind <- dplyr::bind_rows(cost_cumavert_all, .id ='scenario')%>%
#   mutate(scenario = factor(scenario, levels = c(names(cost_cum_all))))
# 
# long.cost_cumavert_all_bind <- gather(cost_cumavert_all_bind, index, best, 
#                                       -c(scenario, year, discount))
# for(i in col_cum_avert){ 
#   
#   ceil_avert[[i]] <- long.cost_cumavert_all_bind%>%filter(index == i & year%in% year_obs)%>%
#     select(best)%>%c()%>%unlist()%>%max()%>%ceiling_dec(.,-2)
#   
#   floor_avert[[i]] <- long.cost_cumavert_all_bind%>%filter(index == i & year%in% year_obs)%>%
#     select(best)%>%c()%>%unlist()%>%min()%>%ceiling_dec(.,-2)
#   
#   if(floor_avert[[i]]>0){
#     plot_cost_cumavert[[i]] <- plot_pocau(POC_AU, long.cost_cumavert_all_bind, 
#                                           type = "avert", 
#                                           indicator = i, year_obs = year_obs) + 
#       labs(y = i) + 
#       ggtitle("Cost saving") +  
#       scale_y_continuous(limits = c(0, ceil_avert[[i]]), breaks = seq(0, ceil_avert[[i]], 
#                                                                       ceil_avert[[i]]/10))
#     
#   }
#   else{ 
#     plot_cost_cumavert[[i]] <- plot_pocau(POC_AU, long.cost_cumavert_all_bind, 
#                                           type = "avert", 
#                                           indicator = i, year_obs = year_obs) + 
#       labs(y = paste0(i)) + 
#       ggtitle("Cost saving") + 
#       scale_y_continuous(limits = c(floor_avert[[i]], ceil_avert[[i]]), breaks = seq(floor_avert[[i]], ceil_avert[[i]], 
#                                                                                      (ceil_avert[[i]]-  floor_avert[[i]])/10))
#     
#   }
#   
# }
# 
# cost_stage <- list()
# cost_stage_avert <- list()
# 
# 
# for(i in names(cost_cum_all)){ 
#   
#   cost_stage[[i]] <- cost_cum_all[[i]]%>%
#     mutate(diagnosis_ycum = cost_ab_ycum + cost_RNA_ycum + cost_POCT_ycum,
#            diagnosis_disycum = cost_ab_disycum + cost_RNA_disycum + cost_POCT_disycum,
#            management_ycum = cost_TreatOther_ycum + cost_RetreatOther_ycum + 
#              cost_compartment_ycum + cost_Cured_ycum,
#            management_disycum = cost_TreatOther_disycum + cost_RetreatOther_disycum + 
#              cost_compartment_disycum + cost_Cured_disycum,,
#            treatment_ycum  = cost_totalDAA_ycum, 
#            treatment_disycum = cost_totalDAA_disycum)%>%
#     select(year, diagnosis_ycum, management_ycum,treatment_ycum, 
#            diagnosis_disycum,management_disycum,treatment_disycum )
# }
# 
# cost_stagesave <- list()
# 
# cost_stagesave_year_obs <- list()
# col_lst <- c("diagnosis_ycum", "management_ycum", 
#              "treatment_ycum","diagnosis_disycum", "management_disycum", 
#              "treatment_disycum" )
# for(i in names(cost_stage)){ 
#   cost_stagesave_year_obs[[i]] <- cost_stage[[i]]%>%filter(year%in% year_obs)
#   cost_stagesave[[i]] <- cbind(year = cost_stagesave_year_obs[[i]]$year, 
#                                data.frame(cost_stagesave_year_obs[[1]][, col_lst] - 
#                                             cost_stagesave_year_obs[[i]][, col_lst]))
#   
# }
# 
# cost_stage_bind <- dplyr::bind_rows(cost_stage, .id = 'scenario')%>%
#   mutate(scenario = factor(scenario, levels = c(names(cost_stage))))%>%
#   gather(index, best, -c(year, scenario))
# 
# cost_stagesave_bind <- dplyr::bind_rows(cost_stagesave, .id = 'scenario')%>%
#   mutate(scenario = factor(scenario, levels = c(names(cost_stage))))%>%
#   gather(index, best, -c(year, scenario))
# 
# p_cost_stage <- list()
# p_cost_stage[["cumulative"]] <-  cost_stage_bind%>%filter(year%in% year_obs & index %in% 
#                            c("diagnosis_ycum", "management_ycum", 
#                              "treatment_ycum"))%>%
#   ggplot(., aes(x = as.character(year), y = best/1000000000, fill = scenario)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   
#   scale_fill_viridis(discrete = TRUE, option = "D") + 
#   theme_linedraw() +
#   labs(x = "Year", y = "Cost, billions")
# 
# p_cost_stage[["discount cumulative"]] <-  cost_stage_bind%>%
#   filter(year%in% year_obs & index %in% c("diagnosis_disycum", "management_disycum", 
#                                                               "treatment_disycum"))%>%
#   ggplot(., aes(x = as.character(year), y = best/1000000000, fill = scenario)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   
#   scale_fill_viridis(discrete = TRUE, option = "D") + 
#   theme_linedraw() +
#   labs(x = "Year", y = "Cost, billions")
# p_cost_stage[["cumulative"]]  <- p_cost_stage[["cumulative"]] + facet_custom (~index,
#                 scales = "free", ncol = 3,
#                 scale_overrides = 
#                   list(
#                     scale_new(1,
#                               scale_y_continuous(limits = c(0,5))),
#                     scale_new(2,
#                               scale_y_continuous(limits = c(0, 5))),
#                     scale_new(3,
#                           scale_y_continuous(limits = c(0, 5))))) + 
#   scale_x_discrete(labels = c(paste0((year_obs - 
#                                         POC_AU$simY + 1), "-Year", sep = "")))
# 
# 
# p_cost_stage[["discount cumulative"]]  <- p_cost_stage[["discount cumulative"]] + 
#   facet_custom (~index,
#                 scales = "free", ncol = 3,
#                 scale_overrides = 
#                   list(
#                     scale_new(1,
#                               scale_y_continuous(limits = c(0, 3))),
#                     scale_new(2,
#                               scale_y_continuous(limits = c(0, 3))),
#                     scale_new(3,
#                               scale_y_continuous(limits = c(0, 3))))) +
#   scale_x_discrete(labels = c(paste0((year_obs - 
#                                         POC_AU$simY + 1), "-Year", sep = "")))
# 
# for(i in names(p_cost_stage)){ 
#   ggsave(file=file.path(OutputFig, paste0("plot_cost_stage_",i,".png")), p_cost_stage[[i]], 
#          width = 10 , height = 6, bg = "white")
#   }
# 
# # saving 
# p_cost_stagesave <- list()
# p_cost_stagesave[["cumulative"]] <-  
#   cost_stagesave_bind%>%
#   filter(year%in% year_obs & index %in% c("diagnosis_ycum", "management_ycum", 
#                                                               "treatment_ycum"))%>%
#   ggplot(., aes(x = as.character(year), y = best/1000000, fill = scenario)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   
#   scale_fill_viridis(discrete = TRUE, option = "D") + 
#   theme_linedraw() +
#   labs(x = "Year", y = "Cost, billions") + 
#   scale_x_discrete(labels = c(paste0((year_obs - 
#                                         POC_AU$simY + 1), "-Year", sep = "")))
# 
# p_cost_stagesave[["discount cumulative"]] + facet_wrap(~index, scale = "free")
# 
# p_cost_stagesave[["discount cumulative"]] <-  cost_stagesave_bind%>%
#   filter(year%in% year_obs & index %in% c("diagnosis_disycum", "management_disycum", 
#                                           "treatment_disycum"))%>%
#   ggplot(., aes(x = as.character(year), y = best/1000000, fill = scenario)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   
#   scale_fill_viridis(discrete = TRUE, option = "D") + 
#   theme_linedraw() +
#   labs(x = "Year", y = "Cost, millions") + 
#   scale_x_discrete(labels = c(paste0((year_obs - 
#                                         POC_AU$simY + 1), "-Year", sep = "")))
# 
# p_cost_stagesave[["cumulative"]]  <- p_cost_stagesave[["cumulative"]] + 
#   facet_custom (~index,
#                 scales = "free", ncol = 3,
#                 scale_overrides = 
#                   list(
#                     scale_new(1,
#                               scale_y_continuous(limits = c(-25, 0), 
#                               breaks = seq(-25, 0, 5))),
#                     scale_new(2,
#                               scale_y_continuous(limits = c(0, 100), 
#                               breaks = seq(0, 100, 10))),
#                     scale_new(3,
#                               scale_y_continuous(limits = c(-200, 200), 
#                                                  breaks = seq(-200, 200, 10))))) + 
#   labs(x = "Year", y = "Cost, millions")
# 
# 
# p_cost_stagesave[["discount cumulative"]]  <- p_cost_stagesave[["discount cumulative"]] + 
#   facet_custom (~index,
#                 scales = "free", ncol = 3,
#                 scale_overrides = 
#                   list(
#                     scale_new(1,
#                               scale_y_continuous(limits = c(-20, 0), 
#                                                  breaks = seq(-20, 0, 5))),
#                     scale_new(2,
#                               scale_y_continuous(limits = c(0, 60), 
#                                                  breaks = seq(0, 60, 10))),
#                     scale_new(3,
#                               scale_y_continuous(limits = c(-200, 10), 
#                                                  breaks = seq(-200, 10, 10)))))
# 
# for(i in names(p_cost_stage)){ 
#   ggsave(file=file.path(OutputFig, paste0("plot_cost_stage_",i,".png")), 
#          p_cost_stage[[i]], 
#          width = 10 , height = 6, bg = "white")
#   
#   ggsave(file=file.path(OutputFig, paste0("plot_cost_stagesave_",i,".png")), 
#          p_cost_stagesave[[i]], 
#          width = 10 , height = 6, bg = "white")
# }
# 
# plot_cost_arrange <- list()
# 
# ylab <- c("Total", "Management at eahc state", "Antibody testing", "RNA testing", 
#           "Immediate RNA testing", "Treat_DAA", "Retreat_DAA", 
#           "Other costs related to treatment", "Other costs related to Retreat", 
#           "Total DAA costs", "Cured", "Total DAA with cap", "Total cost with cap", 
#           "Discounted total cost with cap",  "Discounted total DAA with cap", 
#           "Discounted total", "Discounted Management at eahc state", 
#           "Discounted antibody testing", "Discounted RNA testing", 
#           "Discounted immediate RNA testing", "Discounted treat_DAA", 
#           "Discounted retreat_DAA", 
#           "Discounted other costs related to treatment", 
#           "Discounted other costs related to Retreat",
#           "Discounted total DAA costs", 
#           "Discounted cured"
#           )
# ylab_billion <- paste(ylab, "(billions)")
# for(i in 1: length(names(plot_cost_cum))){ 
#   
#   plot_cost_cum[[i]] <- plot_cost_cum[[i]] + 
#     labs(y = ylab_billion[i]) 
#   
#   plot_cost_cumavert[[i]] <- plot_cost_cumavert[[i]] + 
#     labs(y = ylab[i])
#   
# }
# 
# 
# 
# plot_cost_cum[[1]] <- plot_cost_cum[[1]] + 
#   scale_y_continuous(limits = c(0, 6))
# plot_cost_cum[[2]] <- plot_cost_cum[[2]] + 
#   scale_y_continuous(limits = c(0, 3))
# plot_cost_cum[[3]] <- plot_cost_cum[[3]] + 
#   scale_y_continuous(limits = c(0, 0.15))
# 
# plot_cost_cum[[4]] <- plot_cost_cum[[4]] + 
#   scale_y_continuous(limits = c(0, 0.5))
# 
# plot_cost_cum[[5]] <- plot_cost_cum[[5]] + 
#   scale_y_continuous(limits = c(0, 0.01), breaks = seq(0, 0.01, 0.001), 
#                      labels = seq(0, 0.01, 0.001)) 
# 
# plot_cost_cum[[6]] <- plot_cost_cum[[6]] + 
#   scale_y_continuous(limits = c(0, 3))
# plot_cost_cum[[7]] <- plot_cost_cum[[7]] + 
#   scale_y_continuous(limits = c(0, 0.03))
# 
# plot_cost_cum[[8]] <- plot_cost_cum[[8]] + 
#   scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.01), 
#                      labels = seq(0, 0.1, 0.01)) 
# 
# plot_cost_cum[[9]] <- plot_cost_cum[[9]] +  
#   scale_y_continuous(limits = c(0, 0.0015), breaks = seq(0, 0.0015, 0.0005), 
#                      labels = seq(0, 0.0015, 0.0005)) 
# 
# plot_cost_cum[[10]] <- plot_cost_cum[[10]] + 
#   scale_y_continuous(limits = c(0, 3))
# 
# plot_cost_cum[[11]] <- plot_cost_cum[[11]] + 
#   scale_y_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01), 
#                      labels = seq(0, 0.04, 0.01))
# 
# plot_cost_cum[[12]] <- plot_cost_cum[[12]] + 
#   scale_y_continuous(limits = c(0, 3))
# plot_cost_cum[[13]] <- plot_cost_cum[[13]] + 
#   scale_y_continuous(limits = c(0, 6))
# plot_cost_cum[[14]] <- plot_cost_cum[[14]] + 
#   scale_y_continuous(limits = c(0, 6))
# plot_cost_cum[[15]] <- plot_cost_cum[[15]] + 
#   scale_y_continuous(limits = c(0, 2))
# 
# plot_cost_cum[[16]] <- plot_cost_cum[[16]] + 
#   scale_y_continuous(limits = c(0, 4))
# 
# 
# plot_cost_cum[[17]] <- plot_cost_cum[[17]] + 
#   scale_y_continuous(limits = c(0, 2))
# 
# plot_cost_cum[[18]] <- plot_cost_cum[[18]] + 
#   scale_y_continuous(limits = c(0, 0.25)) 
# plot_cost_cum[[19]] <- plot_cost_cum[[19]] + 
#   scale_y_continuous(limits = c(0, 0.5))
# 
# plot_cost_cum[[20]] <- plot_cost_cum[[20]] + 
#   scale_y_continuous(limits = c(0, 0.007), breaks = seq(0, 0.007, 0.001), 
#                      labels = seq(0, 0.007, 0.001)) 
# plot_cost_cum[[21]] <- plot_cost_cum[[21]] + 
#   scale_y_continuous(limits = c(0, 3))
# 
# plot_cost_cum[[22]] <- plot_cost_cum[[22]] + 
#   scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.01), 
#                      labels = seq(0, 0.02, 0.01)) 
# 
# plot_cost_cum[[23]] <- plot_cost_cum[[23]] + 
#   scale_y_continuous(limits = c(0, 0.1))
# 
# plot_cost_cum[[24]] <- plot_cost_cum[[24]] + 
#   scale_y_continuous(limits = c(0, 0.001))
# 
# plot_cost_cum[[25]] <- plot_cost_cum[[25]] + 
#   scale_y_continuous(limits = c(0, 2))
# 
# plot_cost_cum[[26]] <- plot_cost_cum[[26]] + 
#   scale_y_continuous(limits = c(0, 0.03))
# 
# plot_cost_cumavert[[1]] <- plot_cost_cumavert[[1]] + 
#   scale_y_continuous(limits = c(-90000000,100000000), 
#                      breaks = seq(-90000000,100000000, 10000000),
#                      labels = seq(-90000000,100000000, 10000000)/1000000) + 
#   labs(y = "Total cost (million)")
# 
# plot_cost_cumavert[[2]] <- plot_cost_cumavert[[2]] + 
#   scale_y_continuous(limits = c(0,50000000), 
#                      breaks = seq(0,50000000, 5000000),
#                      labels = seq(0,50000000, 5000000)/1000000) + 
#   labs(y = paste0(ylab [2], "(million)"))
# 
# plot_cost_cumavert[[3]] <- plot_cost_cumavert[[3]] + 
#   scale_y_continuous(limits = c(-4000000,0), 
#                      breaks = seq(-4000000,0, 500000),
#                      labels = seq(-4000000,0, 500000)/1000000) +
#   labs(y = paste0(ylab [3], "(million)"))
# 
# plot_cost_cumavert[[4]] <- plot_cost_cumavert[[4]] + 
#   scale_y_continuous(limits = c(-20000000,0), 
#                      breaks = seq(-20000000,0, 1000000),
#                      labels = seq(-20000000,0, 1000000)/1000000) +
#   labs(y = paste0(ylab [4], "(million)"))
# 
# plot_cost_cumavert[[5]] <- plot_cost_cumavert[[5]] + 
#   scale_y_continuous(limits = c(-9000000,0), 
#                      breaks = seq(-9000000,0, 1000000),
#                      labels = seq(-9000000,0, 1000000)/1000000) +
#   labs(y = paste0(ylab [5], "(million)"))
# 
# plot_cost_cumavert[[6]] <- plot_cost_cumavert[[6]] + 
#   scale_y_continuous(limits = c(-80000000,80000000), 
#                      breaks = seq(-80000000,80000000, 5000000),
#                      labels = seq(-80000000,80000000, 5000000)/1000000) +
#   labs(y = paste0(ylab [6], "(million)"))
# 
# 
# plot_cost_cumavert[[7]] <- plot_cost_cumavert[[7]] + 
#   scale_y_continuous(limits = c(-1000000,500000), 
#                      breaks = seq(-1000000,500000, 100000),
#                      labels = seq(-1000000,500000, 100000)/1000000) +
#   labs(y = paste0(ylab [7], "(million)"))
# 
# plot_cost_cumavert[[8]] <- plot_cost_cumavert[[8]] + 
#   scale_y_continuous(limits = c(-4000000,3000000), 
#                      breaks = seq(-4000000,3000000, 500000),
#                      labels = seq(-4000000,3000000, 500000)/1000000) +
#   labs(y = paste0(ylab [8], "(million)"))
# 
# plot_cost_cumavert[[9]] <- plot_cost_cumavert[[9]] + 
#   scale_y_continuous(limits = c(-40000,20000), 
#                      breaks = seq(-40000,20000, 5000)) +
#   labs(y = paste0(ylab [9]))
# 
# plot_cost_cumavert[[10]] <- plot_cost_cumavert[[10]] + 
#   scale_y_continuous(limits = c(-80000000,60000000), 
#                      breaks = seq(-80000000,60000000, 5000000),
#                      labels = seq(-80000000,60000000, 5000000)/1000000) +
#   labs(y = paste0(ylab [10], "(million)"))
# 
# plot_cost_cumavert[[11]] <- plot_cost_cumavert[[11]] + 
#   scale_y_continuous(limits = c(-1000000,800000), 
#                      breaks = seq(-1000000,800000, 100000),
#                      labels = seq(-1000000,800000, 100000)/1000000) +
#   labs(y = paste0(ylab [11], "(million)"))
# 
# plot_cost_cumavert[[12]] <- plot_cost_cumavert[[12]] + 
#   scale_y_continuous(limits = c(-70000000,60000000), 
#                      breaks = seq(-70000000,60000000, 5000000),
#                      labels = seq(-70000000,60000000, 5000000)/1000000) +
#   labs(y = paste0(ylab [12], "(million)"))
# 
# plot_cost_cumavert[[13]] <- plot_cost_cumavert[[13]] + 
#   scale_y_continuous(limits = c(-90000000,90000000), 
#                      breaks = seq(-90000000,90000000, 10000000),
#                      labels = seq(-90000000,90000000, 10000000)/1000000) +
#   labs(y = paste0(ylab [13], "(million)"))
# 
# plot_cost_cumavert[[14]] <- plot_cost_cumavert[[14]] + 
#   scale_y_continuous(limits = c(-80000000,10000000), 
#                      breaks = seq(-80000000,10000000, 10000000),
#                      labels = seq(-80000000,10000000, 10000000)/1000000) +
#   labs(y = paste0(ylab [14], "(million)"))
# 
# plot_cost_cumavert[[15]] <- plot_cost_cumavert[[15]] + 
#   scale_y_continuous(limits = c(-60000000,10000000), 
#                      breaks = seq(-60000000,10000000, 5000000),
#                      labels = seq(-60000000,10000000, 5000000)/1000000) +
#   labs(y = paste0(ylab [15], "(million)"))
# 
# plot_cost_cumavert[[16]] <- plot_cost_cumavert[[16]] + 
#   scale_y_continuous(limits = c(-80000000,10000000), 
#                      breaks = seq(-80000000,10000000, 10000000),
#                      labels = seq(-80000000,10000000, 10000000)/1000000) +
#   labs(y = paste0(ylab [16], "(million)"))
# 
# plot_cost_cumavert[[17]] <- plot_cost_cumavert[[17]] + 
#   scale_y_continuous(limits = c(0,30000000), 
#                      breaks = seq(0,30000000, 5000000),
#                      labels = seq(0,30000000, 5000000)/1000000) +
#   labs(y = paste0(ylab [17], "(million)"))
# 
# plot_cost_cumavert[[18]] <- plot_cost_cumavert[[18]] + 
#   scale_y_continuous(limits = c(-3000000,0), 
#                      breaks = seq(-3000000,0, 500000),
#                      labels = seq(-3000000,0, 500000)/1000000) +
#   labs(y = paste0(ylab [18], "(million)"))
# 
# plot_cost_cumavert[[19]] <- plot_cost_cumavert[[19]] + 
#   scale_y_continuous(limits = c(-12000000,0), 
#                      breaks = seq(-12000000,0, 1000000),
#                      labels = seq(-12000000,0, 1000000)/1000000) +
#   labs(y = paste0(ylab [19], "(million)"))
# 
# plot_cost_cumavert[[20]] <- plot_cost_cumavert[[20]] + 
#   scale_y_continuous(limits = c(-8000000,0), 
#                      breaks = seq(-8000000,0, 1000000),
#                      labels = seq(-8000000,0, 1000000)/1000000) +
#   labs(y = paste0(ylab [20], "(million)"))
# 
# plot_cost_cumavert[[21]] <- plot_cost_cumavert[[21]] + 
#   scale_y_continuous(limits = c(-65000000,6000000), 
#                      breaks = seq(-65000000,6000000, 5000000),
#                      labels = seq(-65000000,6000000, 5000000)/1000000) +
#   labs(y = paste0(ylab [21], "(million)"))
# 
# plot_cost_cumavert[[22]] <- plot_cost_cumavert[[22]] + 
#   scale_y_continuous(limits = c(-600000,50000), 
#                      breaks = c(seq(-600000,0, 50000), 50000),
#                      labels = c(seq(-600000,0, 50000), 50000)/1000000) +
#   labs(y = paste0(ylab [22], "(million)"))
# 
# plot_cost_cumavert[[23]] <- plot_cost_cumavert[[23]] + 
#   scale_y_continuous(limits = c(-3000000,300000), 
#                      breaks = seq(-3000000,300000, 200000),
#                      labels = seq(-3000000,300000, 200000)/1000000) +
#   labs(y = paste0(ylab [23], "(million)"))
# 
# plot_cost_cumavert[[24]] <- plot_cost_cumavert[[24]] + 
#   scale_y_continuous(limits = c(-30000,5000), 
#                      breaks = seq(-30000,5000, 5000)) +
#   labs(y = paste0(ylab [24]))
# 
# plot_cost_cumavert[[25]] <- plot_cost_cumavert[[25]] + 
#   scale_y_continuous(limits = c(-65000000,5500000), 
#                      breaks = seq(-65000000,5500000, 5000000),
#                      labels = seq(-65000000,5500000, 5000000)/1000000) +
#   labs(y = paste0(ylab [25], "(million)"))
# 
# plot_cost_cumavert[[26]] <- plot_cost_cumavert[[26]] + 
#   scale_y_continuous(limits = c(-900000,100000), 
#                      breaks = seq(-900000,100000, 50000),
#                      labels = seq(-900000,100000, 50000)/1000000) +
#   labs(y = paste0(ylab [26], "(million)"))
# 
# 
# 
# for(i in 1:length(names(plot_cost_cum))){ 
#   
#   plot_cost_arrange[[i]] <- ggpubr::ggarrange(
#     plotlist = list(plot_cost_cum[[i]], plot_cost_cumavert[[i]]),
#     common.legend = TRUE, 
#     nrow=1)  + 
#     theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))  
#   
#   ggsave(file=file.path(OutputFig, paste0("plot_cost_each",names(plot_cost_cum)[i],".png")),  
#          plot_cost_arrange[[i]], 
#          width = 10 , height = 6, bg = "white")
#   
#   }
# 
# # combine plot together 
# # plot_cost_cum, plot_cost_cumavert
# # p_cost_stage
# t_string <- names(plot_cost_cum)
# 
# t <- str_match(t_string, "cost_\\s*(.*?)\\s*_ycum")
# t <- t[c(1:13) ,]
# t
# t_discum <- str_match(t_string, "cost_\\s*(.*?)\\s*_disycum")
# 
# 
# 
# 
# 
# 
# 
# 
# # save cost data 
# write.xlsx(cost_year_all, file = file.path(OutputFolder, paste0("POC_AU_cost_y.xlsx")), 
# append=TRUE) 
# write.xlsx(cost_cum_all, file = file.path(OutputFolder, paste0("POC_AU_costcum_y.xlsx")), 
#            append=TRUE) 
# 
# 


