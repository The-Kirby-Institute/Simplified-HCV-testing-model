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
  Resflow_year_all[[i]][["Tot_Treatment"]] <- 
    cbind(year = Resflow_year_all[[i]][["Treatment"]]$year, 
          as.data.frame(Resflow_year_all[[i]][["Treatment"]][, par_col] + 
                          Resflow_year_all[[i]][["Retreat"]][, par_col] + 
                          Resflow_year_all[[i]][["Treatment_sc"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_ab"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_ab"]]$year, 
          as.data.frame(Resflow_year_all[[i]][["Testing_ab"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_ab_neg"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_ab_sc"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_ab_sc_neg"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_RNA"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_RNA"]]$year, 
          as.data.frame(Resflow_year_all[[i]][["Testing_RNA"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_RNA_neg"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_RNA_sc"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )
  
  Resflow_year_all[[i]][["Tot_Testing_POCT"]] <- 
    cbind(year = Resflow_year_all[[i]][["Testing_POCT"]]$year, 
          as.data.frame(Resflow_year_all[[i]][["Testing_POCT"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_POCT_neg"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_POCT_sc"]][, par_col] + 
                          Resflow_year_all[[i]][["Testing_POCT_sc_neg"]][, par_col])
    )
}

for(i in names(Resflow_year_all)){
  Resflow_year_pop[[i]][["Tot_Treatment"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Treatment"]]$year, 
          population = Resflow_year_pop[[i]][["Treatment"]]$population,
          as.data.frame(Resflow_year_pop[[i]][["Treatment"]][, par_col] + 
                          Resflow_year_pop[[i]][["Retreat"]][, par_col] + 
                          Resflow_year_pop[[i]][["Treatment_sc"]][, par_col])
    )
  
  Resflow_year_pop[[i]][["Tot_Testing_ab"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_ab"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_ab"]]$population,
          as.data.frame(Resflow_year_pop[[i]][["Testing_ab"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_ab_neg"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_ab_sc"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_ab_sc_neg"]][, par_col])
    )
  
  Resflow_year_pop[[i]][["Tot_Testing_RNA"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_RNA"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_RNA"]]$population,
          as.data.frame(Resflow_year_pop[[i]][["Testing_RNA"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_RNA_neg"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_RNA_sc"]][, par_col] + 
                          Resflow_year_pop[[i]][["Testing_RNA_sc_neg"]][, par_col])
    )
  
  Resflow_year_pop[[i]][["Tot_Testing_POCT"]] <- 
    cbind(year = Resflow_year_pop[[i]][["Testing_POCT"]]$year, 
          population = Resflow_year_pop[[i]][["Testing_POCT"]]$population,
          as.data.frame(Resflow_year_pop[[i]][["Testing_POCT"]][, par_col] + 
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
      data.frame(Resflow_sc_cum_all[["Status quo"]][[indic]][, c(par_col)] - 
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
  Resflow_all_lst[[i]] <- Resflow_all[[i]]%>%transpose()%>%
    lapply(., function(x) dplyr::bind_rows(x, .id = 'scenario'))
  
  for(n in names(Resflow_all_lst[[i]])){ 
    Resflow_all_lst[[i]][[n]] <- Resflow_all_lst[[i]][[n]]%>%
      mutate(scenario = factor(scenario, 
                               levels = c("Status quo", "dfList_NP_2023", 
                                          "dfList_NP_2024", "dfList_NPexp_A", 
                                          "dfList_NPexp_B", "dfList_NPexp_C",
                                          "dfList_NPexp_D"), 
                               labels = c("Pre national program", "Achievement 2023", 
                                          "Achievement 2024", "NP expand 2024", 
                                          "NP expand 2025", "NP expand 2026", 
                                          "NP expand 2027")))
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
    summarise(x = max(q95))%>%
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
      x >100000 & x <= 1000000  ~ (x%/%100000 + 1)*100000
    ))
  
  return(lim)  
}


for(i in names(Resflow_all_lst_subsce$Resflow_cum_avert)){ 
  Resflow_all_lst_subsce$Resflow_cum_avert[[i]] <- 
    Resflow_all_lst_subsce$Resflow_cum_avert[[i]]%>%
    filter(scenario != "Pre national program")
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
    ylab(ylab_name[[indic]]) + 
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

p_pocau_avert$Treatment <- p_pocau_avert$Treatment + 
  scale_y_continuous(limits = c(-10000, 20000), 
                     breaks = seq(-10000, 20000, 5000)) 

p_pocau_avert$Retreat <- p_pocau_avert$Retreat + 
  scale_y_continuous(limits = c(-100, 100), 
                     breaks = seq(-100, 100, 20))

p_pocau_avert$Testing_ab <- p_pocau_avert$Testing_ab + 
  scale_y_continuous(limits = c(0, 30000), 
                     breaks = seq(0, 30000, 5000))

p_pocau_avert$Testing_ab_neg <- p_pocau_avert$Testing_ab_neg + 
  scale_y_continuous(limits = c(-6000, 0), 
                     breaks = seq(-6000, 0, 1000))

p_pocau_avert$Testing_RNA_neg <- p_pocau_avert$Testing_RNA_neg + 
  scale_y_continuous(limits = c(-8000, 0), 
                     breaks = seq(-8000, 0, 1000))

p_pocau_avert$Testing_POCT_neg <- p_pocau_avert$Testing_POCT_neg + 
  scale_y_continuous(limits = c(-8000, 0), 
                     breaks = seq(-8000, 0, 1000))

p_pocau_avert$Cured <- p_pocau_avert$Cured + 
  scale_y_continuous(limits = c(-9000, 12000), 
                     breaks = seq(-9000, 12000, 3000))

p_pocau_avert$Treatment_sc <- p_pocau_avert$Treatment_sc + 
  scale_y_continuous(limits = c(-15000, 0), 
                     breaks = seq(-15000, 0, 3000))

p_pocau_avert$Testing_ab_sc <- p_pocau_avert$Testing_ab_sc + 
  scale_y_continuous(limits = c(-10000, 0), 
                     breaks = seq(-10000, 0, 1000))

p_pocau_avert$Testing_RNA_sc <- p_pocau_avert$Testing_RNA_sc + 
  scale_y_continuous(limits = c(-6000, 0), 
                     breaks = seq(-6000, 0, 1000))

p_pocau_avert$Testing_POCT_sc <- p_pocau_avert$Testing_POCT_sc + 
  scale_y_continuous(limits = c(-15000, 0), 
                     breaks = seq(-15000, 0, 3000))
p_pocau_avert$Testing_ab_sc_neg <- p_pocau_avert$Testing_ab_sc_neg + 
  scale_y_continuous(limits = c(-80000, 0), 
                     breaks = seq(-80000, 0, 5000))

p_pocau_avert$Testing_RNA_sc_neg <- p_pocau_avert$Testing_RNA_sc_neg + 
  scale_y_continuous(limits = c(-40000,100), 
                     breaks = seq(-40000, 100, 5000))

p_pocau_avert$Testing_POCT_sc_neg <- p_pocau_avert$Testing_POCT_sc_neg + 
  scale_y_continuous(limits = c(-250000,0), 
                     breaks = seq(-250000, 0, 50000))

p_pocau_avert$Tot_Treatment <- p_pocau_avert$Tot_Treatment + 
  scale_y_continuous(limits = c(-10000,10000), 
                     breaks = seq(-10000,10000, 2500))

p_pocau_avert$Tot_Testing_ab <- p_pocau_avert$Tot_Testing_ab + 
  scale_y_continuous(limits = c(-70000,2000), 
                     breaks = c(seq(-70000,0, 5000), c(2000)))
p_pocau_avert$Tot_Testing_RNA <- p_pocau_avert$Tot_Testing_RNA + 
  scale_y_continuous(limits = c(-40000,10000), 
                     breaks = seq(-40000,10000, 5000))

p_pocau_avert$Tot_Testing_POCT <- p_pocau_avert$Tot_Testing_POCT + 
  scale_y_continuous(limits = c(-250000,0), 
                     breaks = seq(-250000,0, 50000))


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
  transpose()%>%
  lapply(., function(x) dplyr::bind_rows(x, .id = 'scenario'))


Resflow_year_setting_range_trajectory <- lapply(Resflow_year_setting_range, function(x)
  x%>%mutate(scenario = factor(scenario, levels = c("Status quo", "dfList_NP_2023", "dfList_NP_2024", 
                                                             "dfList_NPexp_A", "dfList_NPexp_B", "dfList_NPexp_C",
                                                             "dfList_NPexp_D"), 
                               labels = c("Pre national program", "Achievement 2023", 
                                          "Achievement 2024", "NP expand 2024", 
                                          "NP expand 2025", "NP expand 2026", 
                                          "NP expand 2027")))%>%

    filter(! scenario %in% c("NP expand 2025", "NP expand 2026"))
)
# plot function for generating cascade numbers 
Cas_num_plot <- function(pj, dt, obdt =NULL, xlimits, UI = NULL){ 
  # @ UI the name of scenario for showing UI range
  if(length(unique(dt$scenario)) == 2){ 
    col_pal <- c("#000000", "#E69F00")
    
  } 
  else{col_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  }
  
  
  
  if(is.null(obdt) & is.null(UI)){ 
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
  else if(!is.null(obdt) & !is.null(UI)){
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
Resflow_year_setting_range_trajectory$Tot_Treatment  <- 
  Resflow_year_setting_range_trajectory$Tot_Treatment%>%mutate(year = year + 1)
p_T_num <- Cas_num_plot(POC_AU,Resflow_year_setting_range_trajectory$Tot_Treatment , 
             obdt =HCVtreatinitN_setting_fit, 
             xlimits = c(1, 31, 5), UI = "Pre national program") + 
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
                                    labels = c("Out of national program", 
                                               "Naitonal program")))

xt_rna <- xt_rna%>%mutate(NP = factor(NP, levels = c("out of NP", "NP"),
                                    labels = c("Out of national program", 
                                               "Naitonal program")))

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
                       "Naitonal program")))

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
                       "Naitonal program")))

HCVNP_ab_setting_fit_lst <- list()
HCVNP_rna_setting_fit_lst <- list()
for(i in c("Achievement 2023", 
           "Achievement 2024", "NP expand 2024", 
           "NP expand 2027")){ 
  HCVNP_ab_setting_fit_lst[[i]] <- as.data.frame(HCVNP_ab_setting_fit)%>%
    mutate(scenario = i)
  
  HCVNP_rna_setting_fit_lst[[i]] <- as.data.frame(HCVNP_rna_setting_fit)%>%
    mutate(scenario = i)

  }


xt_ab <- xt_ab%>%
  mutate(scenario = factor(scenario, 
                           levels = c("Pre national program", "Achievement 2023", 
                                      "Achievement 2024", "NP expand 2024",  "NP expand 2027")))%>%
  arrange(scenario)

parea_ab <- list()
parea_ab[["Pre national program"]] <- 
  ggplot(data = xt_ab%>%filter(scenario == "Pre national program"), aes(x = year, y = best)) +
  geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
            alpha=0.6 , size=0.2, colour="black") +
  scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
  facet_wrap(~population, scales='free_y') + 
  theme_Publication_facet() + 
  scale_x_continuous(limits = c(1,31), breaks = seq(1,31,5), 
                     labels = seq(1,31,5) + POC_AU$cabY - 1)  + 
  labs(x = "Year", y = "Number of antibody testing", title = "Pre national program") + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 20000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 10000)))))

parea_rna <- list()
parea_rna[["Pre national program"]] <- 
  ggplot(data = xt_rna%>%filter(scenario == "Pre national program"), aes(x = year, y = best)) +
  geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
            alpha=0.6 , size=0.2, colour="black") +
  scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
  facet_wrap(~population, scales='free_y') + 
  theme_Publication_facet() + 
  scale_x_continuous(limits = c(1,31), breaks = seq(1,31,5), 
                     labels = seq(1,31,5) + POC_AU$cabY - 1)  + 
  labs(x = "Year", y = "Number of RNA testing", title = "Pre national program") + 
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
    scale_x_continuous(limits = c(1,31), breaks = seq(1,31,5), 
                       labels = seq(1,31,5) + POC_AU$cabY - 1)  + 
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
                                                     c(0, 10000)))))
  
  parea_rna[[i]] <- ggplot(data = xt_rna%>%filter(scenario == i), 
                          aes(x = year, y = best)) +
    geom_area(aes(x = year, y = best, fill = NP,colour = scenario), 
              alpha=0.6 , size=0.2, colour="black") +
    scale_fill_manual(" ", values=c ("#FFC57D","#1a979d"))   +
    geom_point(data=HCVNP_rna_setting_fit_lst[[i]], aes(y=best, x = time), 
               colour = "black", size = 1) + 
    facet_wrap(~population, scales='free_y') + 
    theme_Publication_facet() + 
    scale_x_continuous(limits = c(1,31), breaks = seq(1,31,5), 
                       labels = seq(1,31,5) + POC_AU$cabY - 1)  + 
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

parea_ab_ggarrange <- ggarrange(plotlist = parea_ab, ncol = 5, nrow = 1, 
                                common.legend = TRUE)
# adding common x, y label 
parea_ab_ggarrange <- annotate_figure(parea_ab_ggarrange, 
                left = textGrob("Number of antibody testing", 
                                rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Year", gp = gpar(cex = 1.3))
                )
ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numAb_arrange", ".png")), 
       parea_ab_ggarrange , 
       width = 14, height = 8, bg = "white", dpi = 300)
# remove y label 
parea_rna <- lapply(parea_rna, function(x) x + rremove("ylab") + rremove("xlab"))

parea_rna_ggarrange <- ggarrange(plotlist = parea_rna, ncol = 5, nrow = 1, 
                                common.legend = TRUE)
# adding common x, y label 
parea_rna_ggarrange <- annotate_figure(parea_rna_ggarrange, 
                                      left = textGrob("Number of RNA testing", 
                                                      rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                      bottom = textGrob("Year", gp = gpar(cex = 1.3)))


ggsave(file=file.path(OutputFig_y_cum_avert, paste0("NP_numRNA_arrange", ".png")), 
       parea_rna_ggarrange, 
       width = 14, height = 8, bg = "white", dpi = 300)

p_num_ab_line <-ggplot(data = xt_ab, aes(x = year, y = best)) +
  geom_line(aes(x = year, y = best, colour = scenario)) + 
  geom_ribbon(data = xt_ab%>%filter((NP == "Out of national program" & scenario == "Pre national program")),
              aes(x = year, ymin = q5, ymax = q95), alpha = 0.2) + 
   
  geom_point(data=HCVNP_ab_setting_fit, aes(y=best, x = time), 
             colour = "black", size = 1) + 
  facet_grid( population~ NP, scales='free_y') + 
  theme_Publication_facet() + 
  scale_x_continuous(limits = c(1,31), breaks = seq(1,31,5), 
                     labels = seq(1,31,5) + POC_AU$cabY - 1)  + 
  labs(x = "Year", y = "Number of antibody testing") 
  
p_num_ab_line <- p_num_ab_line + 
  facet_custom (population~ NP,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 12000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 15000))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 18000))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 6000)))
                    
                  )) +
  theme(legend.position = "right", legend.direction="vertical")

p_num_rna_line <-ggplot(data = xt_rna, aes(x = year, y = best)) +
  geom_line(aes(x = year, y = best, colour = scenario)) + 
  geom_ribbon(data = xt_rna%>%filter((NP == "Out of national program" & scenario == "Pre national program")),
              aes(x = year, ymin = q5, ymax = q95), alpha = 0.2) + 
  
  geom_point(data=HCVNP_rna_setting_fit, aes(y=best, x = time), 
             colour = "black", size = 1) + 
  facet_grid( population~ NP, scales='free_y') + 
  theme_Publication_facet() + 
  scale_x_continuous(limits = c(1,31), breaks = seq(1,31,5), 
                     labels = seq(1,31,5) + POC_AU$cabY - 1)  + 
  labs(x = "Year", y = "Number of RNA testing") 

p_num_rna_line <- p_num_rna_line + 
  facet_custom (population~ NP,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 30000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 15000))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 20000))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 30000)))
                    
                  )) +
  theme(legend.position = "right", legend.direction="vertical")




ggsave(file=file.path(OutputFig_y_cum_avert, paste0("T_num" ,".png")), 
       p_T_num, 
       width = 9, height = 6, bg = "white", dpi = 300)



ggsave(file=file.path(OutputFig_y_cum_avert, paste0("Ab_num_line" ,".png")), 
       p_num_ab_line, 
       width = 12, height = 6, bg = "white", dpi = 300)
ggsave(file=file.path(OutputFig_y_cum_avert, paste0("RNA_num_line" ,".png")), 
       p_num_rna_line, 
       width = , height = 6, bg = "white", dpi = 300)

################################################################################
# number of advanced liver diseases 
Res_Numbox_y <- list()
Res_Numbox_cum <- list()
Res_Numbox_avert <- list()
for(i in names(Res_numbox)){ 
  Res_Numbox_y[[i]][["DC"]] <- Res_numbox[[i]]%>%filter(disease_prog == "dc")%>%
    group_by(year)%>%
    summarise(across(c(par_col),~ mean(.x, na.rm = FALSE)))
  Res_Numbox_y[[i]][["HCC"]] <- Res_numbox[[i]]%>%filter(disease_prog == "hcc")%>%
    group_by(year)%>%
    summarise(across(c(par_col),~ mean(.x, na.rm = FALSE)))
  Res_Numbox_y[[i]][["LT"]] <- Res_numbox[[i]]%>%filter(disease_prog == "lt")%>%
    group_by(year)%>%
    summarise(across(c(par_col),~ mean(.x, na.rm = FALSE)))
  Res_Numbox_y[[i]][["PLT"]] <- Res_numbox[[i]]%>%filter(disease_prog == "plt")%>%
    group_by(year)%>%
    summarise(across(c(par_col),~ mean(.x, na.rm = FALSE)))
  
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

Res_Numbox_y_range <- Res_Numbox_y_range%>%transpose()
Res_Numbox_cum_range <- Res_Numbox_cum_range%>%transpose()
Res_Numbox_avert_range <- Res_Numbox_avert_range%>%transpose()

for(i in names(Res_Numbox_y_range)){ 
  Res_Numbox_y_range[[i]] <- dplyr::bind_rows(Res_Numbox_y_range[[i]], .id = "scenario")%>%
    mutate(year = year + POC_AU$cabY - 1)%>%
    mutate(scenario = factor(scenario, 
                             levels = c("Status quo", "dfList_NP_2023", 
                                        "dfList_NP_2024", "dfList_NPexp_A", 
                                        "dfList_NPexp_B", "dfList_NPexp_C",
                                        "dfList_NPexp_D"), 
                             labels = c("Pre national program", "Achievement 2023", 
                                        "Achievement 2024", "NP expand 2024", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "NP expand 2027")))%>%
    filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
  
  Res_Numbox_cum_range[[i]] <- dplyr::bind_rows(Res_Numbox_cum_range[[i]], .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = c("Status quo", "dfList_NP_2023", 
                                        "dfList_NP_2024", "dfList_NPexp_A", 
                                        "dfList_NPexp_B", "dfList_NPexp_C",
                                        "dfList_NPexp_D"), 
                             labels = c("Pre national program", "Achievement 2023", 
                                        "Achievement 2024", "NP expand 2024", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "NP expand 2027")))%>%
    filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
  Res_Numbox_avert_range[[i]] <- dplyr::bind_rows(Res_Numbox_avert_range[[i]], .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = c("Status quo", "dfList_NP_2023", 
                                        "dfList_NP_2024", "dfList_NPexp_A", 
                                        "dfList_NPexp_B", "dfList_NPexp_C",
                                        "dfList_NPexp_D"), 
                             labels = c("Pre national program", "Achievement 2023", 
                                        "Achievement 2024", "NP expand 2024", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "NP expand 2027")))%>%
    filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
    
  }

for(i in names(Res_Numbox_avert_range)){ 
  Res_Numbox_avert_range[[i]] <- Res_Numbox_avert_range[[i]]%>%
    filter(scenario != "Pre national program")
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
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_y[[i]][, "lim"]))) + 
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]])
  
  p_num_adliver_cum[[i]] <- plot_pocau(POC_AU, Res_Numbox_cum_range, type = "cum", 
             indicator = i, year_obs = year_obs) + 
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_cum[[i]][, "lim"]))) +
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]])
  
  p_num_adliver_avert[[i]] <- plot_pocau(POC_AU, Res_Numbox_avert_range, type = "avert", 
                                       indicator = i, year_obs = year_obs) + 
    scale_y_continuous(limits = c(0, as.numeric(lim_adliver_avert[[i]][, "lim"]))) +
    theme(legend.position = "right", legend.direction="vertical") + 
    labs(y = y_adlab_name[[i]])
  
  }

p_num_adliver_y$DC <- p_num_adliver_y$DC + scale_y_continuous(limits = c(0, 80))
p_num_adliver_y$HCC <- p_num_adliver_y$HCC + scale_y_continuous(limits = c(0, 20))
p_num_adliver_y$LT <- p_num_adliver_y$LT + scale_y_continuous(limits = c(0, 4))
p_num_adliver_y$PLT <- p_num_adliver_y$PLT + scale_y_continuous(limits = c(0, 80))

p_num_adliver_avert$DC <- p_num_adliver_avert$DC + 
  scale_y_continuous(limits = c(0, 80)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") + 
  labs(y = "Number averted of Decompensated cirrhosis")

p_num_adliver_avert$HCC <- p_num_adliver_avert$HCC + 
  scale_y_continuous(limits = c(-5, 20)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  labs(y = "Number averted of Hepatocellular carcinoma")

p_num_adliver_avert$LT <- p_num_adliver_avert$LT + 
  scale_y_continuous(limits = c(-5, 5)) + 
  geom_hline(yintercept = 0, size = 0.8, linetype = "dashed") +
  labs(y = "Number averted of liver transplant")

p_num_adliver_avert$PLT <- p_num_adliver_avert$PLT + 
  scale_y_continuous(limits = c(-5, 25)) + 
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
         width = , height = 6, bg = "white", dpi = 300)
  
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("cumnum_adliver_",i  ,".png")), 
         p_num_adliver_cum[[i]], 
         width = , height = 6, bg = "white", dpi = 300)
  ggsave(file=file.path(OutputFig_y_cum_avert, paste0("avert_adliver_",i  ,".png")), 
         p_num_adliver_avert[[i]], 
         width = , height = 6, bg = "white", dpi = 300)
  
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
nlim_avert <- list("newInfections" = 15000, "HCVdeath" = 200)
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

# adjust achievement 2023 
plot_num_avert$`Achievement 2023`$newInfections <- 
  plot_num_avert$`Achievement 2023`$newInfections + scale_y_continuous(limits = c(0,5000))
plot_num_avert$`Achievement 2023`$HCVdeath <- 
  plot_num_avert$`Achievement 2023`$HCVdeath + scale_y_continuous(limits = c(0,100))
plot_num_avert$`Achievement 2023`$DC <- 
  plot_num_avert$`Achievement 2023`$DC + scale_y_continuous(limits = c(0,10))
plot_num_avert$`Achievement 2023`$HCC <- 
  plot_num_avert$`Achievement 2023`$HCC + scale_y_continuous(limits = c(0,5))
plot_num_avert$`Achievement 2023`$LT <- 
  plot_num_avert$`Achievement 2023`$LT + scale_y_continuous(limits = c(0,2))
plot_num_avert$`Achievement 2023`$PLT <- 
  plot_num_avert$`Achievement 2023`$PLT + scale_y_continuous(limits = c(0,2))
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
Rescost_year_all$dfList_NP_2023$cost_total_Cap[, c(1:6)]

names(Rescost_year_all$dfList_NP_2023)

cost_y_categories <- list()
cost_disyear_categories <- list()
for(i in names(Rescost_year_all)){ 
  cost_y_categories[[i]][["Diagnosis"]] <- cbind(
    year = Rescost_year_all[[i]]$cost_ab$year, 
    as.data.frame(Rescost_year_all[[i]]$cost_ab[, par_col] + 
                    Rescost_year_all[[i]]$cost_RNA[, par_col] +
                    Rescost_year_all[[i]]$cost_POCT[, par_col]))%>%
    popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  
  cost_y_categories[[i]][["Treatment_cap"]] <- Rescost_year_all[[i]]$cost_totalDAA_Cap%>%
    popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  
  cost_y_categories[[i]][["Treatment"]] <- Rescost_year_all[[i]]$cost_totalDAA%>%
    popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
  cost_y_categories[[i]][["Management"]] <- cbind(
    year = Rescost_year_all[[i]]$cost_compartment$year, 
    as.data.frame(Rescost_year_all[[i]]$cost_compartment[, par_col] + 
                    Rescost_year_all[[i]]$cost_Cured[, par_col] +
                    Rescost_year_all[[i]]$cost_TreatOther[, par_col] + 
                    Rescost_year_all[[i]]$cost_RetreatOther[, par_col]))%>%
    popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  
  cost_y_categories[[i]] <- cost_y_categories[[i]]%>%
    dplyr::bind_rows(., .id = "Categories")
  
  
  cost_disyear_categories[[i]][["Diagnosis"]] <- cbind(
    year = Rescost_disyear_all[[i]]$cost_ab$year, 
    as.data.frame(Rescost_disyear_all[[i]]$cost_ab[, par_col] + 
                    Rescost_disyear_all[[i]]$cost_RNA[, par_col] +
                    Rescost_disyear_all[[i]]$cost_POCT[, par_col]))%>%
    popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
    
  
  cost_disyear_categories[[i]][["Treatment_cap"]] <- Rescost_disyear_all[[i]]$cost_totalDAA_Cap%>%
    popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  
  cost_disyear_categories[[i]][["Treatment"]] <- Rescost_disyear_all[[i]]$cost_totalDAA%>%
    popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  
  cost_disyear_categories[[i]][["Management"]] <- cbind(
    year = Rescost_disyear_all[[i]]$cost_compartment$year, 
    as.data.frame(Rescost_disyear_all[[i]]$cost_compartment[, par_col] + 
                    Rescost_disyear_all[[i]]$cost_Cured[, par_col] +
                    Rescost_disyear_all[[i]]$cost_TreatOther[, par_col] + 
                    Rescost_disyear_all[[i]]$cost_RetreatOther[, par_col]))%>%
    popResults_range(POC_AU, ., Population = NULL, end_Y = endY-1)
  
  cost_disyear_categories[[i]] <- cost_disyear_categories[[i]]%>%
    dplyr::bind_rows(., .id = "Categories")
    
}

y_cost_disyear_categories <- cost_disyear_categories%>% dplyr::bind_rows(., .id = "scenario")%>%
  group_by(scenario, year)%>%
  summarise(across(c(par_col),~ sum(.x, na.rm = TRUE))) 


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
  )

ggplot(y_cost_disyear_categories, aes(x = year, colour = scenario)) + 
  geom_line(aes(y = best)) + 
  scale_x_continuous(expand = c(0,0),limits = c(2022, 2045), breaks = seq(2022, 2045, 1)) + 
  geom_ribbon(data = y_cost_disyear_categories_range, aes(x = year, ymin = q25, ymax = q75, 
                                                          fill = scenario), 
              alpha =0.2, colour = NA) + 
  scale_y_continuous(limits = c(0, 800000000), breaks = seq(0, 800000000,100000000), 
                     labels = seq(0, 800000000,100000000)/1000000000) + 
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
                           labels = c("Pre national program", "Achievement 2023", 
                                      "Achievement 2024", "NP expand 2024", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "NP expand 2027")))%>%
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
                           labels = c("Pre national program", "Achievement 2023", 
                                      "Achievement 2024", "NP expand 2024", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "NP expand 2027")))%>%
  filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))

cost_ydaanocap_categories <- cost_y_categories%>%filter(Categories != "Treatment_cap")


# undiscount cap 

cost_ydaacap_categories <- cost_y_categories%>%filter(Categories != "Treatment")

# discount no cap 
cost_disyear_categories <- dplyr::bind_rows(cost_disyear_categories, .id = "scenario")
str(cost_disyear_categories)


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
                           labels = c("Pre national program", "Achievement 2023", 
                                      "Achievement 2024", "NP expand 2024", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "NP expand 2027")))%>%
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
                  label = round(best/1000000000, digits = 3))) + 
    geom_bar(position="stack", stat="identity") + 
    theme(panel.spacing = unit(0, 'lines')) +
    scale_fill_manual(values = c( "grey10", "grey40","grey80")) + 
    theme_Publication() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    scale_y_continuous(limit = c(0, 4000000000), 
                       breaks = seq(0, 4000000000, 500000000),
                       labels = seq(0, 4000000000, 500000000)/1000000000) + 
    labs(y = "Cost (discounted, billions)") + 
    geom_text(aes(x = scenario, y = best + 50000000, 
                  label = paste0(format(round(best/1000000000, digits = 3), nsmall = 3), "b"), 
                  group = Categories),
              position = position_stack(vjust = 0.5)) + 
    ggtitle(title_name[i])
  }



pcatcost <- lapply(pcatcost, function(x) x + rremove("ylab") + rremove("xlab"))

pcatcost  <- ggarrange(plotlist = pcatcost , ncol = 3, nrow = 1, 
                                 common.legend = TRUE)
# adding common x, y label 
pcatcost <- annotate_figure(pcatcost, 
                                       left = textGrob("Cost (discounted, billions)", 
                                                       rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                       bottom = textGrob("Scenarios", gp = gpar(cex = 1.3)))

ggsave(file=file.path(OutputFig, paste0("Reports/cost_catego",".png")), 
       pcatcost, 
       width = 18, height = 10, bg = "white", dpi = 300)  

x_catcost_total%>%filter(year == 2027)

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
          as.data.frame(x_total_ref[, c(par_col)] - x[[i]][, c(par_col)]))%>%
    as.data.frame()%>%
    popResults_range(POC_AU, ., end_Y = 100)
  

}

x_total_incre <- x_catcost_total_incre%>%dplyr::bind_rows(., .id = "scenario")%>%
  mutate(scenario = factor(scenario, levels = c("Status quo", "dfList_NP_2023", 
                                                "dfList_NP_2024", "dfList_NPexp_A", 
                                                "dfList_NPexp_B", "dfList_NPexp_C",
                                                "dfList_NPexp_D"), 
                           labels = c("Pre national program", "Achievement 2023", 
                                      "Achievement 2024", "NP expand 2024", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "NP expand 2027")))%>%
  filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
  
incremental_cost <- ggplot(x_total_incre, aes(x = year, colour = scenario) ) + 
  geom_line(aes(x = year, y = best, colour = scenario,
                linetype = scenario) 
            ) + 
  theme_Publication() + 
  scale_x_continuous(expand = c(0,0), limits = c(2022,2045), breaks = seq(2022, 2045, 1)) + 
  scale_y_continuous(limits = c(-100000000, 300000000), 
                     breaks = seq(-100000000, 300000000, 50000000), 
                     labels = seq(-100000000, 300000000, 50000000)/1000000) + 
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) + 
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid", "solid")) + 
  labs( y = "Cost savings, in millions")
  
ggsave(file=file.path(OutputFig, paste0("Reports/incremental_cost ",".png")), 
       incremental_cost , 
       width = 12, height = 8, bg = "white", dpi = 300)   

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
    labs(y = "Cost savings (millions)", x = "Time frame")  + 
    scale_y_continuous(limits = c(-30000000, 150000000), 
                       breaks = seq(-30000000, 150000000, 10000000), 
                       labels = seq(-30000000, 150000000, 10000000)/1000000) + 
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
         width = 6, height = 8, bg = "white", dpi = 300)   
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
                              as.data.frame(catcost_cum_ref[, par_col] - x[[i]][, par_col]))%>%
    as.data.frame()%>%
    popResults_range(POC_AU, .)%>%
    filter(Categories != "Treatment") 
  }

catcost_incre <- dplyr::bind_rows(catcost_incre, .id = "scenario")
catcost_incre <- catcost_incre%>%
  mutate(scenario = factor(scenario, levels = c("Status quo", "dfList_NP_2023", 
                                                "dfList_NP_2024", "dfList_NPexp_A", 
                                                "dfList_NPexp_B", "dfList_NPexp_C",
                                                "dfList_NPexp_D"), 
                           labels = c("Pre national program", "Achievement 2023", 
                                      "Achievement 2024", "NP expand 2024", 
                                      "NP expand 2025", "NP expand 2026", 
                                      "NP expand 2027")))%>%
  filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
catcost_incre <- catcost_incre%>%
  mutate(Categories = factor(Categories, 
                             levels = c("Diagnosis", "Management", "Treatment_cap"),
                             labels = c("Diagnosis", "Management", "Treatment")))
p_catcost_saving <- list()
lim_catcost_saving <- list(c(-50000000, 20000000), c(-50000000, 50000000), c(-20000000, 150000000))
bek_catcost_saving <- list(seq(-50000000, 20000000, 5000000), 
                           seq(-50000000, 50000000, 5000000), 
                           seq(-50000000, 150000000, 10000000))
for(i in 1: length(year_obs)){ 
  
  p_catcost_saving[[i]] <- 
    ggplot(catcost_incre%>%
             filter(scenario != "Pre national program" & year %in% year_obs[i])%>%
             arrange(Categories), 
           aes(x = scenario, y = best, fill = Categories)) + 
    geom_bar(stat = "identity", width = 0.8) + 
    geom_text(aes(y = best + 2 * sign(best), 
                  label=paste0(round(best/1000000, digits = 1), "m")),
              position = position_stack(vjust = 0.5)) + 
    scale_fill_manual(values = c( "grey30", "grey50","grey80")) + 
    theme_Publication() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(y = "Cost savings (millions)", x = "Time frame")  + 
    scale_y_continuous(limits = lim_catcost_saving[[i]], 
                       breaks = bek_catcost_saving[[i]], 
                       labels = bek_catcost_saving[[i]]/1000000) + 
    theme_Publication() + 
    theme(panel.grid.major = element_line(color = "gray80",
                                          size = 0.1,
                                          linetype = 1)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    ggtitle(paste0(POC_AU$simY, "-", year_obs[i], 
                   " (",year_obs[i] - POC_AU$simY + 1 ,"-Year)"))
}
for(i in 1: length(year_obs)){ 
  ggsave(file=file.path(OutputFig, paste0("Reports/cat_cost_saving", i,".png")), 
         p_catcost_saving[[i]] , 
         width = 8, height = 6, bg = "white", dpi = 300)   
}

################################################################################

####                            Tables generation                           ####

################################################################################

# generating tables 
# Resflow all indicators 
# DC, HCC, LT, PLT 
# best (95% PI) 
# substract years 
# making functions 

# test  
tab_epi <- Resflow_all_lst
cost_qaly_range <- list() 
cost_qaly_range_disy <- list()
for(i in names(Rescost_year_all)){ 
  for(n in names(Rescost_year_all[[1]])){ 
    cost_qaly_range[[i]][[n]] <- 
      popResults_range(POC_AU, Rescost_year_all[[i]][[n]], end_Y = 100-1)%>%as_tibble()
    
    cost_qaly_range_disy[[i]][[n]] <- 
      popResults_range(POC_AU, Rescost_disyear_all[[i]][[n]], end_Y = 100-1)%>%as_tibble()
  }
}

tab_costqaly <- cost_qaly_range%>%purrr::transpose()
tab_costqaly_disy <- cost_qaly_range_disy%>%purrr::transpose()

tab_costqaly <- lapply(tab_costqaly, function(x){ 
  x%>%dplyr::bind_rows(., .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = c("Status quo", "dfList_NP_2023", 
                                        "dfList_NP_2024", "dfList_NPexp_A", 
                                        "dfList_NPexp_B", "dfList_NPexp_C",
                                        "dfList_NPexp_D"), 
                             labels = c("Pre national program", "Achievement 2023", 
                                        "Achievement 2024", "NP expand 2024", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "NP expand 2027")))
} )

tab_costqaly_disy <- lapply(tab_costqaly_disy, function(x){ 
  x%>%dplyr::bind_rows(., .id = "scenario")%>%
    mutate(scenario = factor(scenario, 
                             levels = c("Status quo", "dfList_NP_2023", 
                                        "dfList_NP_2024", "dfList_NPexp_A", 
                                        "dfList_NPexp_B", "dfList_NPexp_C",
                                        "dfList_NPexp_D"), 
                             labels = c("Pre national program", "Achievement 2023", 
                                        "Achievement 2024", "NP expand 2024", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "NP expand 2027")))
} )




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
    filter(year %in% seq(2022,2030,1))%>%
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


tab_cost <- tab_costqaly%>%dplyr::bind_rows(., .id ="Indicators")%>%
  filter(year %in% seq(2022,2030,1))%>%
  select(Indicators, scenario, year, best, q5, q95)%>%
  mutate(best = formatC(best,  format = "fg", big.mark = ","),
         q5 = formatC(q5,  format = "fg", big.mark = ","),
         q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
  mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
  select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)
tab_cost%>%group_by(Indicators)%>%
  gt(groupname_col = "Indicators",
     rowname_col = "year")%>%
  gtsave(., file = file.path(OutputFig, paste0("Reports/year_cost.docx")))

tab_cost_disy <- tab_costqaly_disy%>%dplyr::bind_rows(., .id ="Indicators")%>%
  filter(year %in% seq(2022,2030,1))%>%
  select(Indicators, scenario, year, best, q5, q95)%>%
  mutate(best = formatC(best,  format = "fg", big.mark = ","),
         q5 = formatC(q5,  format = "fg", big.mark = ","),
         q95 = formatC(q95,  format = "fg", big.mark = ","))%>%
  mutate(vv = paste0(best, "\n", "(", q5, "-", q95, ")"))%>%
  select(-c(best, q5, q95))%>%ungroup()%>%spread(scenario, vv)

tab_cost_disy%>%group_by(Indicators)%>%
  gt(groupname_col = "Indicators",
     rowname_col = "year")%>%
  gtsave(., file = file.path(OutputFig, paste0("Reports/disyear_cost.docx")))

# numbox
tab_numbox = list("year" = Res_Numbox_y_range, 
                  "cumyear" = Res_Numbox_cum_range, 
                  "cumavert" = Res_Numbox_avert_range)

for(i in names(tab_numbox)){ 
  tab_numbox[[i]] <- tab_numbox[[i]]%>%dplyr::bind_rows(., .id ="Indicators")%>%
    filter(year %in% seq(2022,2030,1))%>%
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


########################## plot animation ######################################
install.packages("gganimate")
library("gganimate")
x_catcost_total <- x_catcost%>%group_by(year, scenario)%>%
  summarise(across(c("best", "q5", "q95"),~ sum(.x, na.rm = FALSE)))%>%
  mutate(year = as.integer(year))%>%filter(year <= 2045)
View(x_catcost_total)
library(forcats)
# creating pause timeline 
df2 <- x_catcost_total %>%
  group_by(year) %>%
  # The * 1 makes it possible to have non-integer ranks while sliding
  mutate(ordering = min_rank(-best) * 1) %>%
  ungroup()%>%
  mutate(text_lab = paste0(format(round(best/1000000000, digits = 3), nsmall = 3), "b"),
         text_y = best)%>%
  mutate(text_lab =as.factor(text_lab))%>% 
  mutate(show_time = case_when(year %in% c(2026,2028,2030, 2031, 2032) ~ 10,
                               TRUE           ~ 1)) %>%
  # uncount is a tidyr function which copies each line 'n' times
  uncount(show_time) %>%
  group_by(year, scenario) %>%
  mutate(reveal_time = row_number()) %>%
  ungroup()

p <- ggplot(df2, aes(ordering, group = scenario, 
                          fill = as.factor(scenario), color = as.factor(scenario))) +
  geom_tile(aes(y = best/2,
                height = best,
                width = 1), alpha = 1, color = NA) +
  # text in x-axis (requires clip = "off" in coord_*)
  # paste(country, " ")  is a hack to make pretty spacing, since hjust > 1 
  #   leads to weird artifacts in text spacing.
  geom_text(aes(y = 0, label = paste(scenario, " ")), vjust = 0.2, hjust = 1, 
            size = 4, fontface = "bold") +
  geom_text(aes(y = text_y, label = text_lab, 
                 group=text_y), color = "white", vjust = 0.2, hjust = 1, 
            size = 8, fontface='bold') + 

  coord_flip(clip = "off", expand = FALSE) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_reverse() + 


  scale_fill_manual(labels = c("", "", "", "", ""),
                    values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
  scale_color_manual(labels = c("", "", "", "", ""),
                    values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
  theme_Publication() + 
  theme(legend.position = "") + 
  guides(color = FALSE, fill = FALSE) +

  labs(title='{closest_state}', x = "", y = "Cumulative discounted cost") +
  theme(plot.title = element_text(hjust = 0, size = 22),
        axis.ticks.y = element_blank(),  # These relate to the axes post-flip
        axis.text.y  = element_blank(),  # These relate to the axes post-flip
        plot.margin = margin(1,1,1,5, "cm")) +
  
  transition_states(year, 
                    transition_length = c(rep(1,4),3,1,3,3,3,3,3, 
                                          rep(1, 12), 12), state_length = 1) +
  ease_aes('cubic-in-out') 

anim <- animate(p, fps = 25, duration = 50, width = 800, height = 600)
magick::image_write(anim, path= paste0(OutputFig, "totalcost_turningyear.gif"))
# make different timeframe animation then compress together 
anim



p <- ggplot(x_catcost_total,
            aes(x=scenario, y=best, label = scenario, color = scenario)) + 
  geom_point(stat = "identity", size = 3) + 
  geom_segment(aes(
    y = 5,
    x = scenario,
    yend = best,
    xend = scenario)
  ) + 
  coord_flip() + 
  theme_Publication() + 
  scale_y_continuous(limits = c(0, 4500000000), breaks = seq(0,4500000000, 500000000), 
                     labels = seq(0,4500000000, 500000000)/1000000000) + 
  labs(y = "Cumulative discounted cost", x = "Scenarios")
  transition_states(year, wrap = FALSE) 


animate(p)



################################################################################

                         # Used #

################################################################################




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

cost_cum_all_bind <- dplyr::bind_rows(cost_cum_all, .id ='scenario')%>%
  mutate(scenario = factor(scenario, levels = c(names(cost_cum_all))))

long.cost_cum_all_bind <- gather(cost_cum_all_bind, index, best, 
                                 -c(scenario, year, discount))

cost_cumavert_all_bind <- dplyr::bind_rows(cost_cumavert_all, .id ='scenario')%>%
  mutate(scenario = factor(scenario, levels = c(names(cost_cum_all))))

long.cost_cumavert_all_bind <- gather(cost_cumavert_all_bind, index, best, 
                                 -c(scenario, year, discount))


plot_cost_cum <- list()
plot_cost_cumavert <- list()
# extracting cumsum releated index
col_cum <- unique(long.cost_cum_all_bind $index)

ceil_avert <- list()
ceil_cum <- list()
floor <- list()
long.cost_cum_all_bind <- long.cost_cum_all_bind%>%mutate(best = best/1000000000)
for(i in col_cum){ 
  
  ceil_cum[[i]] <- long.cost_cum_all_bind%>%
    filter(index == i & year%in% year_obs)%>%select(best)%>%
    c()%>%unlist()%>%max()%>%ceiling_dec(.,-1)
  floor[[i]] <- long.cost_cum_all_bind%>%filter(index == i)%>%select(best)%>%
    c()%>%unlist()%>%min()%>%floor_dec(.,-1)
  
  if(floor[[i]]>=0){ 
    plot_cost_cum[[i]] <- plot_pocau(POC_AU, long.cost_cum_all_bind, type = "cum", 
                                     indicator = i, year_obs = year_obs) + 
      labs(y = paste0(i, "(billions)")) + 
      ggtitle("Cumulative cost") + 
      scale_y_continuous(limits = c(0, ceil_cum[[i]]), breaks = seq(0, ceil_cum[[i]], 
                                                                    ceil_cum[[i]]/10))
  }
  else{ 
    
    plot_cost_cum[[i]] <- plot_pocau(POC_AU, long.cost_cum_all_bind, type = "cum", 
                                     indicator = i, year_obs = year_obs) + 
      labs(y = paste0(i, "(billions)")) + 
      ggtitle("Cumulative cost") + 
      scale_y_continuous(limits = c(floor[[i]], ceil_cum[[i]]), breaks = seq(floor[[i]], ceil_cum[[i]], 
                                                                             (ceil_cum[[i]]- floor[[i]])/10))
  }
  
}

plot_cost_cum$cost_total_ycum

col_cum_avert <- c(paste0("avert_", col_cum))
floor_avert <- list()
plot_cost_cumavert <- list() 

col_lst <- colnames(cost_year_all[[1]])[-c(1,2)]
cost_cum_all_year_obs <- list()
cost_cumavert_all <- list()
for(i in names(cost_cum_all)){
  
  cost_cum_all_year_obs[[i]] <- cost_cum_all[[i]]%>%filter(year %in% year_obs)
  cost_cumavert_all[[i]] <- cbind(year = cost_cum_all_year_obs[[i]]$year, 
                                  discount = cost_cum_all_year_obs[[i]]$discount, 
                                  data.frame(cost_cum_all_year_obs[[1]][, col_lst] - 
                                               cost_cum_all_year_obs[[i]][, col_lst]))%>%
    rename_at(vars(col_lst), ~ paste0("avert_", col_lst))
  
  
}

cost_cumavert_all_bind <- dplyr::bind_rows(cost_cumavert_all, .id ='scenario')%>%
  mutate(scenario = factor(scenario, levels = c(names(cost_cum_all))))

long.cost_cumavert_all_bind <- gather(cost_cumavert_all_bind, index, best, 
                                      -c(scenario, year, discount))
for(i in col_cum_avert){ 
  
  ceil_avert[[i]] <- long.cost_cumavert_all_bind%>%filter(index == i & year%in% year_obs)%>%
    select(best)%>%c()%>%unlist()%>%max()%>%ceiling_dec(.,-2)
  
  floor_avert[[i]] <- long.cost_cumavert_all_bind%>%filter(index == i & year%in% year_obs)%>%
    select(best)%>%c()%>%unlist()%>%min()%>%ceiling_dec(.,-2)
  
  if(floor_avert[[i]]>0){
    plot_cost_cumavert[[i]] <- plot_pocau(POC_AU, long.cost_cumavert_all_bind, 
                                          type = "avert", 
                                          indicator = i, year_obs = year_obs) + 
      labs(y = i) + 
      ggtitle("Cost saving") +  
      scale_y_continuous(limits = c(0, ceil_avert[[i]]), breaks = seq(0, ceil_avert[[i]], 
                                                                      ceil_avert[[i]]/10))
    
  }
  else{ 
    plot_cost_cumavert[[i]] <- plot_pocau(POC_AU, long.cost_cumavert_all_bind, 
                                          type = "avert", 
                                          indicator = i, year_obs = year_obs) + 
      labs(y = paste0(i)) + 
      ggtitle("Cost saving") + 
      scale_y_continuous(limits = c(floor_avert[[i]], ceil_avert[[i]]), breaks = seq(floor_avert[[i]], ceil_avert[[i]], 
                                                                                     (ceil_avert[[i]]-  floor_avert[[i]])/10))
    
  }
  
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

cost_stage_bind <- dplyr::bind_rows(cost_stage, .id = 'scenario')%>%
  mutate(scenario = factor(scenario, levels = c(names(cost_stage))))%>%
  gather(index, best, -c(year, scenario))

cost_stagesave_bind <- dplyr::bind_rows(cost_stagesave, .id = 'scenario')%>%
  mutate(scenario = factor(scenario, levels = c(names(cost_stage))))%>%
  gather(index, best, -c(year, scenario))

p_cost_stage <- list()
p_cost_stage[["cumulative"]] <-  cost_stage_bind%>%filter(year%in% year_obs & index %in% 
                           c("diagnosis_ycum", "management_ycum", 
                             "treatment_ycum"))%>%
  ggplot(., aes(x = as.character(year), y = best/1000000000, fill = scenario)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  
  scale_fill_viridis(discrete = TRUE, option = "D") + 
  theme_linedraw() +
  labs(x = "Year", y = "Cost, billions")

p_cost_stage[["discount cumulative"]] <-  cost_stage_bind%>%
  filter(year%in% year_obs & index %in% c("diagnosis_disycum", "management_disycum", 
                                                              "treatment_disycum"))%>%
  ggplot(., aes(x = as.character(year), y = best/1000000000, fill = scenario)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  
  scale_fill_viridis(discrete = TRUE, option = "D") + 
  theme_linedraw() +
  labs(x = "Year", y = "Cost, billions")
p_cost_stage[["cumulative"]]  <- p_cost_stage[["cumulative"]] + facet_custom (~index,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = c(0,5))),
                    scale_new(2,
                              scale_y_continuous(limits = c(0, 5))),
                    scale_new(3,
                          scale_y_continuous(limits = c(0, 5))))) + 
  scale_x_discrete(labels = c(paste0((year_obs - 
                                        POC_AU$simY + 1), "-Year", sep = "")))


p_cost_stage[["discount cumulative"]]  <- p_cost_stage[["discount cumulative"]] + 
  facet_custom (~index,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = c(0, 3))),
                    scale_new(2,
                              scale_y_continuous(limits = c(0, 3))),
                    scale_new(3,
                              scale_y_continuous(limits = c(0, 3))))) +
  scale_x_discrete(labels = c(paste0((year_obs - 
                                        POC_AU$simY + 1), "-Year", sep = "")))

for(i in names(p_cost_stage)){ 
  ggsave(file=file.path(OutputFig, paste0("plot_cost_stage_",i,".png")), p_cost_stage[[i]], 
         width = 10 , height = 6, bg = "white")
  }

# saving 
p_cost_stagesave <- list()
p_cost_stagesave[["cumulative"]] <-  
  cost_stagesave_bind%>%
  filter(year%in% year_obs & index %in% c("diagnosis_ycum", "management_ycum", 
                                                              "treatment_ycum"))%>%
  ggplot(., aes(x = as.character(year), y = best/1000000, fill = scenario)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  
  scale_fill_viridis(discrete = TRUE, option = "D") + 
  theme_linedraw() +
  labs(x = "Year", y = "Cost, billions") + 
  scale_x_discrete(labels = c(paste0((year_obs - 
                                        POC_AU$simY + 1), "-Year", sep = "")))

p_cost_stagesave[["discount cumulative"]] + facet_wrap(~index, scale = "free")

p_cost_stagesave[["discount cumulative"]] <-  cost_stagesave_bind%>%
  filter(year%in% year_obs & index %in% c("diagnosis_disycum", "management_disycum", 
                                          "treatment_disycum"))%>%
  ggplot(., aes(x = as.character(year), y = best/1000000, fill = scenario)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  
  scale_fill_viridis(discrete = TRUE, option = "D") + 
  theme_linedraw() +
  labs(x = "Year", y = "Cost, millions") + 
  scale_x_discrete(labels = c(paste0((year_obs - 
                                        POC_AU$simY + 1), "-Year", sep = "")))

p_cost_stagesave[["cumulative"]]  <- p_cost_stagesave[["cumulative"]] + 
  facet_custom (~index,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = c(-25, 0), 
                              breaks = seq(-25, 0, 5))),
                    scale_new(2,
                              scale_y_continuous(limits = c(0, 100), 
                              breaks = seq(0, 100, 10))),
                    scale_new(3,
                              scale_y_continuous(limits = c(-200, 200), 
                                                 breaks = seq(-200, 200, 10))))) + 
  labs(x = "Year", y = "Cost, millions")


p_cost_stagesave[["discount cumulative"]]  <- p_cost_stagesave[["discount cumulative"]] + 
  facet_custom (~index,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = c(-20, 0), 
                                                 breaks = seq(-20, 0, 5))),
                    scale_new(2,
                              scale_y_continuous(limits = c(0, 60), 
                                                 breaks = seq(0, 60, 10))),
                    scale_new(3,
                              scale_y_continuous(limits = c(-200, 10), 
                                                 breaks = seq(-200, 10, 10)))))

for(i in names(p_cost_stage)){ 
  ggsave(file=file.path(OutputFig, paste0("plot_cost_stage_",i,".png")), 
         p_cost_stage[[i]], 
         width = 10 , height = 6, bg = "white")
  
  ggsave(file=file.path(OutputFig, paste0("plot_cost_stagesave_",i,".png")), 
         p_cost_stagesave[[i]], 
         width = 10 , height = 6, bg = "white")
}

plot_cost_arrange <- list()

ylab <- c("Total", "Management at eahc state", "Antibody testing", "RNA testing", 
          "Immediate RNA testing", "Treat_DAA", "Retreat_DAA", 
          "Other costs related to treatment", "Other costs related to Retreat", 
          "Total DAA costs", "Cured", "Total DAA with cap", "Total cost with cap", 
          "Discounted total cost with cap",  "Discounted total DAA with cap", 
          "Discounted total", "Discounted Management at eahc state", 
          "Discounted antibody testing", "Discounted RNA testing", 
          "Discounted immediate RNA testing", "Discounted treat_DAA", 
          "Discounted retreat_DAA", 
          "Discounted other costs related to treatment", 
          "Discounted other costs related to Retreat",
          "Discounted total DAA costs", 
          "Discounted cured"
          )
ylab_billion <- paste(ylab, "(billions)")
for(i in 1: length(names(plot_cost_cum))){ 
  
  plot_cost_cum[[i]] <- plot_cost_cum[[i]] + 
    labs(y = ylab_billion[i]) 
  
  plot_cost_cumavert[[i]] <- plot_cost_cumavert[[i]] + 
    labs(y = ylab[i])
  
}



plot_cost_cum[[1]] <- plot_cost_cum[[1]] + 
  scale_y_continuous(limits = c(0, 6))
plot_cost_cum[[2]] <- plot_cost_cum[[2]] + 
  scale_y_continuous(limits = c(0, 3))
plot_cost_cum[[3]] <- plot_cost_cum[[3]] + 
  scale_y_continuous(limits = c(0, 0.15))

plot_cost_cum[[4]] <- plot_cost_cum[[4]] + 
  scale_y_continuous(limits = c(0, 0.5))

plot_cost_cum[[5]] <- plot_cost_cum[[5]] + 
  scale_y_continuous(limits = c(0, 0.01), breaks = seq(0, 0.01, 0.001), 
                     labels = seq(0, 0.01, 0.001)) 

plot_cost_cum[[6]] <- plot_cost_cum[[6]] + 
  scale_y_continuous(limits = c(0, 3))
plot_cost_cum[[7]] <- plot_cost_cum[[7]] + 
  scale_y_continuous(limits = c(0, 0.03))

plot_cost_cum[[8]] <- plot_cost_cum[[8]] + 
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.01), 
                     labels = seq(0, 0.1, 0.01)) 

plot_cost_cum[[9]] <- plot_cost_cum[[9]] +  
  scale_y_continuous(limits = c(0, 0.0015), breaks = seq(0, 0.0015, 0.0005), 
                     labels = seq(0, 0.0015, 0.0005)) 

plot_cost_cum[[10]] <- plot_cost_cum[[10]] + 
  scale_y_continuous(limits = c(0, 3))

plot_cost_cum[[11]] <- plot_cost_cum[[11]] + 
  scale_y_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01), 
                     labels = seq(0, 0.04, 0.01))

plot_cost_cum[[12]] <- plot_cost_cum[[12]] + 
  scale_y_continuous(limits = c(0, 3))
plot_cost_cum[[13]] <- plot_cost_cum[[13]] + 
  scale_y_continuous(limits = c(0, 6))
plot_cost_cum[[14]] <- plot_cost_cum[[14]] + 
  scale_y_continuous(limits = c(0, 6))
plot_cost_cum[[15]] <- plot_cost_cum[[15]] + 
  scale_y_continuous(limits = c(0, 2))

plot_cost_cum[[16]] <- plot_cost_cum[[16]] + 
  scale_y_continuous(limits = c(0, 4))


plot_cost_cum[[17]] <- plot_cost_cum[[17]] + 
  scale_y_continuous(limits = c(0, 2))

plot_cost_cum[[18]] <- plot_cost_cum[[18]] + 
  scale_y_continuous(limits = c(0, 0.25)) 
plot_cost_cum[[19]] <- plot_cost_cum[[19]] + 
  scale_y_continuous(limits = c(0, 0.5))

plot_cost_cum[[20]] <- plot_cost_cum[[20]] + 
  scale_y_continuous(limits = c(0, 0.007), breaks = seq(0, 0.007, 0.001), 
                     labels = seq(0, 0.007, 0.001)) 
plot_cost_cum[[21]] <- plot_cost_cum[[21]] + 
  scale_y_continuous(limits = c(0, 3))

plot_cost_cum[[22]] <- plot_cost_cum[[22]] + 
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.01), 
                     labels = seq(0, 0.02, 0.01)) 

plot_cost_cum[[23]] <- plot_cost_cum[[23]] + 
  scale_y_continuous(limits = c(0, 0.1))

plot_cost_cum[[24]] <- plot_cost_cum[[24]] + 
  scale_y_continuous(limits = c(0, 0.001))

plot_cost_cum[[25]] <- plot_cost_cum[[25]] + 
  scale_y_continuous(limits = c(0, 2))

plot_cost_cum[[26]] <- plot_cost_cum[[26]] + 
  scale_y_continuous(limits = c(0, 0.03))

plot_cost_cumavert[[1]] <- plot_cost_cumavert[[1]] + 
  scale_y_continuous(limits = c(-90000000,100000000), 
                     breaks = seq(-90000000,100000000, 10000000),
                     labels = seq(-90000000,100000000, 10000000)/1000000) + 
  labs(y = "Total cost (million)")

plot_cost_cumavert[[2]] <- plot_cost_cumavert[[2]] + 
  scale_y_continuous(limits = c(0,50000000), 
                     breaks = seq(0,50000000, 5000000),
                     labels = seq(0,50000000, 5000000)/1000000) + 
  labs(y = paste0(ylab [2], "(million)"))

plot_cost_cumavert[[3]] <- plot_cost_cumavert[[3]] + 
  scale_y_continuous(limits = c(-4000000,0), 
                     breaks = seq(-4000000,0, 500000),
                     labels = seq(-4000000,0, 500000)/1000000) +
  labs(y = paste0(ylab [3], "(million)"))

plot_cost_cumavert[[4]] <- plot_cost_cumavert[[4]] + 
  scale_y_continuous(limits = c(-20000000,0), 
                     breaks = seq(-20000000,0, 1000000),
                     labels = seq(-20000000,0, 1000000)/1000000) +
  labs(y = paste0(ylab [4], "(million)"))

plot_cost_cumavert[[5]] <- plot_cost_cumavert[[5]] + 
  scale_y_continuous(limits = c(-9000000,0), 
                     breaks = seq(-9000000,0, 1000000),
                     labels = seq(-9000000,0, 1000000)/1000000) +
  labs(y = paste0(ylab [5], "(million)"))

plot_cost_cumavert[[6]] <- plot_cost_cumavert[[6]] + 
  scale_y_continuous(limits = c(-80000000,80000000), 
                     breaks = seq(-80000000,80000000, 5000000),
                     labels = seq(-80000000,80000000, 5000000)/1000000) +
  labs(y = paste0(ylab [6], "(million)"))


plot_cost_cumavert[[7]] <- plot_cost_cumavert[[7]] + 
  scale_y_continuous(limits = c(-1000000,500000), 
                     breaks = seq(-1000000,500000, 100000),
                     labels = seq(-1000000,500000, 100000)/1000000) +
  labs(y = paste0(ylab [7], "(million)"))

plot_cost_cumavert[[8]] <- plot_cost_cumavert[[8]] + 
  scale_y_continuous(limits = c(-4000000,3000000), 
                     breaks = seq(-4000000,3000000, 500000),
                     labels = seq(-4000000,3000000, 500000)/1000000) +
  labs(y = paste0(ylab [8], "(million)"))

plot_cost_cumavert[[9]] <- plot_cost_cumavert[[9]] + 
  scale_y_continuous(limits = c(-40000,20000), 
                     breaks = seq(-40000,20000, 5000)) +
  labs(y = paste0(ylab [9]))

plot_cost_cumavert[[10]] <- plot_cost_cumavert[[10]] + 
  scale_y_continuous(limits = c(-80000000,60000000), 
                     breaks = seq(-80000000,60000000, 5000000),
                     labels = seq(-80000000,60000000, 5000000)/1000000) +
  labs(y = paste0(ylab [10], "(million)"))

plot_cost_cumavert[[11]] <- plot_cost_cumavert[[11]] + 
  scale_y_continuous(limits = c(-1000000,800000), 
                     breaks = seq(-1000000,800000, 100000),
                     labels = seq(-1000000,800000, 100000)/1000000) +
  labs(y = paste0(ylab [11], "(million)"))

plot_cost_cumavert[[12]] <- plot_cost_cumavert[[12]] + 
  scale_y_continuous(limits = c(-70000000,60000000), 
                     breaks = seq(-70000000,60000000, 5000000),
                     labels = seq(-70000000,60000000, 5000000)/1000000) +
  labs(y = paste0(ylab [12], "(million)"))

plot_cost_cumavert[[13]] <- plot_cost_cumavert[[13]] + 
  scale_y_continuous(limits = c(-90000000,90000000), 
                     breaks = seq(-90000000,90000000, 10000000),
                     labels = seq(-90000000,90000000, 10000000)/1000000) +
  labs(y = paste0(ylab [13], "(million)"))

plot_cost_cumavert[[14]] <- plot_cost_cumavert[[14]] + 
  scale_y_continuous(limits = c(-80000000,10000000), 
                     breaks = seq(-80000000,10000000, 10000000),
                     labels = seq(-80000000,10000000, 10000000)/1000000) +
  labs(y = paste0(ylab [14], "(million)"))

plot_cost_cumavert[[15]] <- plot_cost_cumavert[[15]] + 
  scale_y_continuous(limits = c(-60000000,10000000), 
                     breaks = seq(-60000000,10000000, 5000000),
                     labels = seq(-60000000,10000000, 5000000)/1000000) +
  labs(y = paste0(ylab [15], "(million)"))

plot_cost_cumavert[[16]] <- plot_cost_cumavert[[16]] + 
  scale_y_continuous(limits = c(-80000000,10000000), 
                     breaks = seq(-80000000,10000000, 10000000),
                     labels = seq(-80000000,10000000, 10000000)/1000000) +
  labs(y = paste0(ylab [16], "(million)"))

plot_cost_cumavert[[17]] <- plot_cost_cumavert[[17]] + 
  scale_y_continuous(limits = c(0,30000000), 
                     breaks = seq(0,30000000, 5000000),
                     labels = seq(0,30000000, 5000000)/1000000) +
  labs(y = paste0(ylab [17], "(million)"))

plot_cost_cumavert[[18]] <- plot_cost_cumavert[[18]] + 
  scale_y_continuous(limits = c(-3000000,0), 
                     breaks = seq(-3000000,0, 500000),
                     labels = seq(-3000000,0, 500000)/1000000) +
  labs(y = paste0(ylab [18], "(million)"))

plot_cost_cumavert[[19]] <- plot_cost_cumavert[[19]] + 
  scale_y_continuous(limits = c(-12000000,0), 
                     breaks = seq(-12000000,0, 1000000),
                     labels = seq(-12000000,0, 1000000)/1000000) +
  labs(y = paste0(ylab [19], "(million)"))

plot_cost_cumavert[[20]] <- plot_cost_cumavert[[20]] + 
  scale_y_continuous(limits = c(-8000000,0), 
                     breaks = seq(-8000000,0, 1000000),
                     labels = seq(-8000000,0, 1000000)/1000000) +
  labs(y = paste0(ylab [20], "(million)"))

plot_cost_cumavert[[21]] <- plot_cost_cumavert[[21]] + 
  scale_y_continuous(limits = c(-65000000,6000000), 
                     breaks = seq(-65000000,6000000, 5000000),
                     labels = seq(-65000000,6000000, 5000000)/1000000) +
  labs(y = paste0(ylab [21], "(million)"))

plot_cost_cumavert[[22]] <- plot_cost_cumavert[[22]] + 
  scale_y_continuous(limits = c(-600000,50000), 
                     breaks = c(seq(-600000,0, 50000), 50000),
                     labels = c(seq(-600000,0, 50000), 50000)/1000000) +
  labs(y = paste0(ylab [22], "(million)"))

plot_cost_cumavert[[23]] <- plot_cost_cumavert[[23]] + 
  scale_y_continuous(limits = c(-3000000,300000), 
                     breaks = seq(-3000000,300000, 200000),
                     labels = seq(-3000000,300000, 200000)/1000000) +
  labs(y = paste0(ylab [23], "(million)"))

plot_cost_cumavert[[24]] <- plot_cost_cumavert[[24]] + 
  scale_y_continuous(limits = c(-30000,5000), 
                     breaks = seq(-30000,5000, 5000)) +
  labs(y = paste0(ylab [24]))

plot_cost_cumavert[[25]] <- plot_cost_cumavert[[25]] + 
  scale_y_continuous(limits = c(-65000000,5500000), 
                     breaks = seq(-65000000,5500000, 5000000),
                     labels = seq(-65000000,5500000, 5000000)/1000000) +
  labs(y = paste0(ylab [25], "(million)"))

plot_cost_cumavert[[26]] <- plot_cost_cumavert[[26]] + 
  scale_y_continuous(limits = c(-900000,100000), 
                     breaks = seq(-900000,100000, 50000),
                     labels = seq(-900000,100000, 50000)/1000000) +
  labs(y = paste0(ylab [26], "(million)"))



for(i in 1:length(names(plot_cost_cum))){ 
  
  plot_cost_arrange[[i]] <- ggpubr::ggarrange(
    plotlist = list(plot_cost_cum[[i]], plot_cost_cumavert[[i]]),
    common.legend = TRUE, 
    nrow=1)  + 
    theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))  
  
  ggsave(file=file.path(OutputFig, paste0("plot_cost_each",names(plot_cost_cum)[i],".png")),  
         plot_cost_arrange[[i]], 
         width = 10 , height = 6, bg = "white")
  
  }

# combine plot together 
# plot_cost_cum, plot_cost_cumavert
# p_cost_stage
t_string <- names(plot_cost_cum)

t <- str_match(t_string, "cost_\\s*(.*?)\\s*_ycum")
t <- t[c(1:13) ,]
t
t_discum <- str_match(t_string, "cost_\\s*(.*?)\\s*_disycum")








# save cost data 
write.xlsx(cost_year_all, file = file.path(OutputFolder, paste0("POC_AU_cost_y.xlsx")), 
append=TRUE) 
write.xlsx(cost_cum_all, file = file.path(OutputFolder, paste0("POC_AU_costcum_y.xlsx")), 
           append=TRUE) 




