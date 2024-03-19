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

# output to spreadsheets 
write.xlsx(Resflow_cum_all, file = file.path(OutputFolder, paste0("POC_AU_epiy.xlsx")), 
           append=TRUE) 
write.xlsx(Resflow_cum_all_avert, file = file.path(OutputFolder, paste0("POC_AU_epiy_avert.xlsx")), 
           append=TRUE)
write.xlsx(Resflow_sc_cum_all, file = file.path(OutputFolder, paste0("POC_AU_epiy_sc.xlsx")), 
           append=TRUE) 
write.xlsx(Resflow_sc_cum_all_avert, file = file.path(OutputFolder, paste0("POC_AU_epiy_sc_avert.xlsx")), 
           append=TRUE)


#### plot for number of averted/cumulative infections by 5, 10, 15 years ####
# ceiling and floor function 



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

long.Resflow_cum_all <- convert_long(Resflow_cum_all)
long.Resflow_cum_all_avert <- convert_long(Resflow_cum_all_avert)

floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

# specific indexs 
col_lst <- colnames(Resflow_year_all$`Status quo`)[-c(1,2)]

plot_flow_year <- list()
ceil <- list()


for(i in col_lst){ 
  ceil[[i]] <- long.Resflow_cum_all%>%filter(index == i)%>%select(best)%>%
    c()%>%unlist()%>%max()%>%ceiling_dec(.,-2)
  plot_flow_year[[i]] <- plot_pocau(POC_AU, long.Resflow_cum_all, type = "new", 
                               indicator = i) + labs(y = i) + 
    scale_y_continuous(limits = c(0, ceil[[i]]), breaks = seq(0, ceil[[i]], 
                                                              ceil[[i]]/10))
}

plot_flow_cum <- list()
plot_flow_cumavert <- list()
# extracting cumsum releated index
col_cum <- unique(long.Resflow_cum_all$index)
col_cum <- col_cum[! col_cum %in% col_lst]
ceil_avert <- list()
ceil_cum <- list()
floor <- list()
for(i in col_cum){ 
  ceil_cum[[i]] <- long.Resflow_cum_all%>%filter(index == i & year%in% year_obs)%>%select(best)%>%
    c()%>%unlist()%>%max()%>%ceiling_dec(.,-2)
  floor[[i]] <- long.Resflow_cum_all%>%filter(index == i)%>%select(best)%>%
    c()%>%unlist()%>%min()%>%floor_dec(.,-2)
  
  if(floor[[i]]>=0){ 
    plot_flow_cum[[i]] <- plot_pocau(POC_AU, long.Resflow_cum_all, type = "cum", 
                                     indicator = i, year_obs = year_obs) + 
      labs(y = i) + 
      scale_y_continuous(limits = c(0, ceil_cum[[i]]), breaks = seq(0, ceil_cum[[i]], 
                                                                ceil_cum[[i]]/10))
  }
  else{ 
    
    plot_flow_cum[[i]] <- plot_pocau(POC_AU, long.Resflow_cum_all, type = "cum", 
                                     indicator = i, year_obs = year_obs) + 
      labs(y = i) + 
      scale_y_continuous(limits = c(floor[[i]], ceil_cum[[i]]), breaks = seq(floor[[i]], ceil_cum[[i]], 
                                                                (ceil_cum[[i]]- floor[[i]])/10))
    }
  
}
col_cum_avert <- c(paste0("avert_", col_cum))
floor_avert <- list()
plot_flow_cumavert <- list()
for(i in col_cum_avert){ 
  
  ceil_avert[[i]] <- long.Resflow_cum_all_avert%>%filter(index == i & year%in% year_obs)%>%
    select(best)%>%c()%>%unlist()%>%max()%>%ceiling_dec(.,-2)
  
  floor_avert[[i]] <- long.Resflow_cum_all_avert%>%filter(index == i & year%in% year_obs)%>%
    select(best)%>%c()%>%unlist()%>%min()%>%floor_dec(.,-2)
  
  if(floor_avert[[i]]>0){
    plot_flow_cumavert[[i]] <- plot_pocau(POC_AU, long.Resflow_cum_all_avert, 
                                          type = "avert", 
                                          indicator = i, year_obs = year_obs) + 
      labs(y = i) + 
      scale_y_continuous(limits = c(0, ceil_avert[[i]]), breaks = seq(0, ceil_avert[[i]], 
                                                                      ceil_avert[[i]]/10))
    
  }
  else{ 
    plot_flow_cumavert[[i]] <- plot_pocau(POC_AU, long.Resflow_cum_all_avert, 
                                          type = "avert", 
                                          indicator = i, year_obs = year_obs) + 
      labs(y = i) + 
      scale_y_continuous(limits = c(floor_avert[[i]], ceil_avert[[i]]), breaks = seq(floor_avert[[i]], ceil_avert[[i]], 
                                                                      (ceil_avert[[i]]-  floor_avert[[i]])/10))
    
    }
  
  }

# ggarrange plot 
plot_flow_arrange <- list()
plot_flow_arrange$Treatment_sc
for(i in 1: length(plot_flow_year)){ 
  plot_flow_arrange[[i]] <- ggpubr::ggarrange(
    plotlist = list(plot_flow_year[[i]], 
                    plot_flow_cum[[i]],
                    plot_flow_cumavert[[i]]),
    common.legend = TRUE, 
    nrow=1)  + 
    theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm")) 
  
  plot_flow_arrange[[i]] <- annotate_figure(plot_flow_arrange[[i]], 
                                            top = text_grob(names(plot_flow_year)[i], color = "black", 
                                                            face = "bold", size = 14))
  
  }

names(plot_flow_arrange) <- names(plot_flow_year)

for(i in names(plot_flow_arrange)){ 
  ggsave(file=file.path(OutputFig, paste0("plot_flow",i,".png")), plot_flow_arrange[[i]], 
         width = 10 , height = 6, bg = "white")
  
  }

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

long.Num_dc_year_bind <- num_sum(Num_dc, endY =100)
long.Num_hcc_year_bind <- num_sum(Num_hcc, endY =100)
long.Num_lt_year_bind <- num_sum(Num_lt, endY =100)
long.num.ad <- list("Decompensated cirrhosis" = long.Num_dc_year_bind, 
                    "Hepatocellular carcinoma" = long.Num_hcc_year_bind, 
                    "Liver transplant" = long.Num_lt_year_bind)
ceil_new <- list()
ceil_cum <- list()
ceil_avert <- list()
floor_avert <- list()
year_obs
for(i in names(long.num.ad)){ 
  ceil_new[[i]] <- long.num.ad[[i]]%>%filter(index == "best")%>%select(best)%>%
    c()%>%unlist()%>%max()%>%ceiling_dec(.,-2)
  
  ceil_cum[[i]] <- long.num.ad[[i]]%>%
    filter(index == "cum_best" & year %in% year_obs)%>%select(best)%>%
    c()%>%unlist()%>%max()%>%ceiling_dec(.,-2)
  
  ceil_avert[[i]] <- long.num.ad[[i]]%>%
    filter(index == "avertcum_best" & year %in% year_obs)%>%select(best)%>%
    c()%>%unlist()%>%max()%>%ceiling_dec(.,-1)
  
  floor_avert[[i]] <- long.num.ad[[i]]%>%
    filter(index == "avertcum_best" & year %in% year_obs)%>%select(best)%>%
    c()%>%unlist()%>%min()%>%floor_dec(.,-1)
  
}


p_num_adliver <- list()
for(i in names(long.num.ad)){ 
  p_num_adliver[[i]][["yearly"]] <-  plot_pocau(POC_AU, long.num.ad[[i]], type = "new", 
                                              indicator = "best") + 
    labs(y = i) + 
    scale_y_continuous(limits = c(0, ceil_new[[i]]), breaks = seq(0, ceil_new[[i]], 
                                                                  ceil_new[[i]]/10))
  
  p_num_adliver[[i]][["cumulative"]] <-  plot_pocau(POC_AU, long.num.ad[[i]], type = "cum", 
                                              indicator = "cum_best", year_obs = year_obs) + 
    labs(y = i) + 
    scale_y_continuous(limits = c(0, ceil_cum[[i]]), breaks = seq(0, ceil_cum[[i]], 
                                                                  ceil_cum[[i]]/10))
  
  if(isTRUE(floor_avert[[i]]>=0)){ 
    p_num_adliver[[i]][["avert_cumulative"]] <- plot_pocau(POC_AU, long.num.ad[[i]], type = "avert", 
                                                           indicator = "avertcum_best", 
                                                           year_obs = year_obs) + 
      labs(y = i) + 
      scale_y_continuous(limits = c(0, ceil_avert[[i]]), breaks = seq(0, ceil_avert[[i]], 
                                                                    1))
  }
  else{ 
    p_num_adliver[[i]][["avert_cumulative"]] <- plot_pocau(POC_AU, long.num.ad[[i]], type = "avert", 
                                                           indicator = "avertcum_best", 
                                                           year_obs = year_obs) + 
      labs(y = i) + 
      scale_y_continuous(limits = c(floor_avert[[i]], ceil_avert[[i]]), 
                         breaks = seq(floor_avert[[i]], ceil_avert[[i]],
                                      1))
    
    
    }
  

  
} 

p_num_adliver$`Decompensated cirrhosis`[["avert_cumulative"]] + 
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 10))

plot_num_arrange <- list()

for(i in names(p_num_adliver)){ 
  plot_num_arrange[[i]] <- ggarrange(plotlist = p_num_adliver[[i]], nrow= 1, common.legend = TRUE)
  ggsave(file=file.path(OutputFig, paste0("plot_Num",i,".png")), plot_num_arrange[[i]], 
         width = 10 , height = 6, bg = "white")
  
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




