# this script is for visulizing the calibraiton outcomes 
# once check the code is running then turn it to markdown file to ouput documents 

rm(list = ls())

library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)

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

load(file.path(OutputFolder, paste0(project_name, "cali_init", ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_param", ".rda")))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 

# setting end Year 

endY <- 40
#### tidy up datapoint #### 

datapoint <- list()
datapoint[["N"]] <- 
  cbind.data.frame(year = rep(POC_AU$simY, length(POC_AU$popNames) +3), 
                   indicator = c(POC_AU$popNames, "Community", "Prison", "Total"), 
                   realpop = c(80000, 400000, 6993, 12696, 20311, 480000, 40000, 520000),
                   low= c(60000, 30000, 6163, 11819, 19332, 360000, 37314, 397314),
                   up = c(100000, 600000, 7869, 13525, 21292, 700000, 42686, 742686))

datapoint[["frac"]] <- 
  cbind.data.frame(year = rep(POC_AU$simY, 3), 
                   indicator = c("commu_proP_fit", "prison_proP_fit", 
                                 "prison_profP_fit"), 
                   realpop = c(16.67, 17.48, 31.74),
                   low= c(16.56, 17.11, 31.28),
                   up = c(16.77, 17.85, 32.7))

datapoint[["flow"]] <- 
  cbind.data.frame(year = rep(POC_AU$simY, 3), 
                   indicator = c("incar", "release", 
                                 "commu_stopinj"), 
                   realpop = c(63753, 63113, 5.4),
                   low= c(NA, NA, 4.3),
                   up = c(NA, NA, 6.7))



#### Result:  population ####
#### N: each subpop/ total  #### 
# subpop as list

for(i in 1: 10){ 
  calibrateInit <- list()
  calibrateInit <- cab[[i]]
  subpop_N <- lapply(POC_AU$popNames, function(x){ 
    
    a <- popResults_MidYear(POC_AU, calibrateInit,
                            Population = x,
                            Disease_prog = NULL, 
                            Cascade = NULL, param = NULL, 
                            endYear = endY)%>%ungroup()
  })
  
  names(subpop_N) <- POC_AU$popNames
  
  # all subpop in one list 
  pop_N <- dplyr::bind_rows(subpop_N, .id = 'population')
  
  # total number 
  
  total_N <- pop_N%>%group_by(year)%>%summarise(best = sum(best))
  
  #### N: community ####
  commu_N <- pop_N%>%filter(population %in% c("C_PWID", "C_fPWID"))%>%
    dplyr::group_by(year)%>%summarise_at("best",sum)
  
  #### N: prison ####
  prison_N <- pop_N%>%filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%
    dplyr::group_by(year)%>%summarise_at("best",sum)
  
  #### N: prisonPWID ####
  prisonPWID_N <- pop_N%>%filter(population %in% c("P_PWID", "P_fPWID"))%>%
    dplyr::group_by(year)%>%summarise_at("best",sum)
  
  #### %: PWID in community/prison ####
  commu_proP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                      as.data.frame(pop_N[pop_N$population =="C_PWID",3] / commu_N[ ,-1])*100)%>%
    tibble::as_tibble()
  
  prison_proP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                       as.data.frame(pop_N[pop_N$population =="P_PWID",3] / prison_N[ ,-1])*100)%>%
    tibble::as_tibble()
  
  #### % fPWID in prison #### 
  prison_profP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                        as.data.frame(pop_N[pop_N$population =="P_fPWID",3] / prison_N[ ,-1])*100)%>%
    tibble::as_tibble()
  
  
  
  #### number of annual leaving in each subpop ####
  # this is for getting annual number of released in non-PWID in prison 
  rel <- indicatorResult_uno(POC_AU, calibrateInit, "newLeave",
                             populations = POC_AU$popNames, endYear= endY)%>%
    mutate(year = year + POC_AU$cabY - 1) 
  
  #### number of annual entry model in each subpop ####
  # to get the annual number of incarceration in non-PWID in prison 
  entry <- indicatorResult_uno(POC_AU, calibrateInit, "newEntry",
                               populations = POC_AU$popNames, endYear= endY) %>%
    mutate(year = year + POC_AU$cabY - 1)
  
  
  
  #### annual number of incarceration/release/injection relapse/ stopping injection #### 
  
  PopTransition <- as.data.frame.table(calibrateInit$newpop_tran)%>%
    mutate(timestep = c(rep(seq(POC_AU$startYear, endY-POC_AU$timestep,
                                POC_AU$timestep),each = POC_AU$npops*POC_AU$npops)),
           from = Var1,
           To = Var2)%>%dplyr::select(-c(Var1, Var2, Var3))%>%
    filter(timestep != 1)
  # giving the average number to first time step 
  #impute <- PopTransition%>%filter(timestep >0 &timestep <2)%>%
  #  group_by(from, To)%>%
  #  mutate(total_pop = mean(Freq))%>%
  #  dplyr::select(total_pop)%>%ungroup()
  
  #impute <- impute[c(1:as.numeric(POC_AU$npops*POC_AU$npops)),]
  
  
  #PopTransition[c(1:as.numeric(POC_AU$npops*POC_AU$npops)), "Freq"] <- impute$total_pop 
  
  PopTransition_all <- cbind.data.frame(timestep = PopTransition$timestep, 
                                        from = PopTransition$from,  
                                        to = PopTransition$To, 
                                        best = PopTransition$Freq)%>%
    as_tibble()%>%
    mutate(year = c(rep(seq(POC_AU$startYear, endY - 2, 1), 
                        each = (1/POC_AU$timestep)*POC_AU$npops*POC_AU$npops),
                    rep(endY - 1 , each = (1/POC_AU$timestep -1)*POC_AU$npops*POC_AU$npops)))
  
  
  PPTranTo <- PopTransition_all%>%
    group_by(year, from, to)%>%summarize(best = sum(best))
  
  # incarceration 
  incarce <- list()
  
  incarce[["PWID"]] <- PPTranTo%>%filter(from == "C_PWID" & to == "P_PWID")
  
  incarce[["fPWID"]] <- PPTranTo%>%filter(from == "C_fPWID" & to == "P_fPWID")
  
  incarce_bind <- dplyr::bind_rows(incarce, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()%>%
    dplyr::select(year, population, best)
  
  entry_nonPWID <- entry%>%filter(population =="P_nPWID") 
  
  # total incarceration = incarceration in C_PWID + incarceration in C_fPWID + incarceration in non-PWID (entry:non-PWID)
  incar_total <- rbind(incarce_bind, entry_nonPWID)%>%group_by(year)%>%
    summarize(best = sum(best))
  
  # release 
  release <- list()
  
  release[["PWID"]] <- PPTranTo%>%filter(from == "P_PWID" & to == "C_PWID")
  
  release[["fPWID"]] <- PPTranTo%>%filter(from == "P_fPWID" & to == "C_fPWID")
  
  release_bind <- dplyr::bind_rows(release, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()%>%
    dplyr::select(year, population, best)
  
  release_nonPWID <- rel%>%filter(population =="P_nPWID") 
  
  release_total <- rbind(release_bind, release_nonPWID)%>%group_by(year)%>%
    summarize(best = sum(best))
  
  # injection relapse 
  inj_relap <- list()
  
  inj_relap[["community"]] <- PPTranTo%>%filter(from == "C_fPWID" & to == "C_PWID")
  
  inj_relap[["prison"]] <- PPTranTo%>%filter(from == "P_fPWID" & to == "P_PWID")
  
  inj_relap_bind <- dplyr::bind_rows(inj_relap, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()
  
  
  # stopping injection 
  inj_stop <- list()
  
  inj_stop[["community"]] <- PPTranTo%>%filter(from == "C_PWID" & to == "C_fPWID")
  
  inj_stop[["prison"]] <- PPTranTo%>%filter(from == "P_PWID" & to == "P_fPWID")
  
  inj_stop_bind <- dplyr::bind_rows(inj_stop, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()
  
  
  
  ####  incidence of stopping injecting #### 
  inj_stop_inc <- list()
  
  inj_stop_inc[["community"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(inj_stop[["community"]][ , "best"] / 
                          subpop_N[["C_PWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  
  inj_stop_inc[["prison"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(inj_stop[["prison"]][ , "best"] / 
                          subpop_N[["P_PWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  inj_stop_inc_bind <- dplyr::bind_rows(inj_stop_inc, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()
  ####  incidence of injection relapse #### 
  inj_relap_inc <- list()
  
  inj_relap_inc[["community"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(inj_relap[["community"]][ , "best"] / 
                          subpop_N[["C_fPWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  inj_relap_inc[["prison"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(inj_relap[["prison"]][ , "best"] / 
                          subpop_N[["P_fPWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  inj_relap_inc_bind <- dplyr::bind_rows(inj_relap_inc, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()
  
  #### incidence of incarceration  ####
  incar_inc <- list()
  
  incar_inc[["PWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                               as.data.frame(incarce[["PWID"]][ , "best"] / 
                                               subpop_N[["C_PWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  incar_inc[["fPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                as.data.frame(incarce[["fPWID"]][ , "best"] / 
                                                subpop_N[["C_fPWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  incar_inc_bind <- dplyr::bind_rows(incar_inc, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()
  
  #### incidence of release  ####
  release_inc <- list()
  
  release_inc[["PWID"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(release[["PWID"]][ , "best"] / 
                          subpop_N[["P_PWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  release_inc[["fPWID"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(release[["fPWID"]][ , "best"] /
                          subpop_N[["P_fPWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  release_inc_bind <- dplyr::bind_rows(release_inc, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()
  
  #==============================================================================#
  
  #==============================================================================#
  
  
  ####Plot: population ####
  #####plot: number of total pop  #####
  
  totalPop_plot <- indicatorPlot(POC_AU, total_N, 
                                 ylabel = "Number (thousands)",
                                 xlimits = c(POC_AU$startYear, 
                                             POC_AU$startYear+30, 5),
                                 calibration_Y = POC_AU$cabY,
                                 rangeun = NULL, 
                                 groupPlot = NULL, 
                                 facetPlot = NULL,
                                 observationData = NULL, 
                                 simulateYear = NULL) + theme_bw() +
    scale_y_continuous(limits = c(300000,825000), 
                       breaks = seq(300000,825000, 100000), 
                       labels = seq(300000,825000, 100000)/1000) + 
    ggtitle("Number of total population") + 
    geom_point(aes(
      y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Total", "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), 
      colour = "black", shape = 5, size = 2)  + 
    geom_segment(
      aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Total", "low"], 
          yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Total", "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  #####plot: number of pop in community/prison   #####
  commuPop_plot <- indicatorPlot(POC_AU, commu_N, 
                                 ylabel = "Number",
                                 xlimits = c(POC_AU$startYear, 
                                             POC_AU$startYear+30, 5),
                                 calibration_Y = POC_AU$cabY,
                                 rangeun = NULL, 
                                 groupPlot = NULL, 
                                 facetPlot = NULL,
                                 observationData = NULL, 
                                 simulateYear = NULL) + theme_bw() +
    ggtitle("Number of total population") + 
    geom_point(aes(
      y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Community", "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), 
      colour = "black", shape = 5, size = 2) + 
    geom_segment(
      aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Community", "low"], 
          yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Community", "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1))) + 
    scale_y_continuous(limits = c(0,750000), breaks = seq(0, 750000, 50000)) 
  
  prisonPop_plot <- indicatorPlot(POC_AU, prison_N, 
                                  ylabel = "Number",
                                  xlimits = c(POC_AU$startYear, 
                                              POC_AU$startYear+30, 5),
                                  calibration_Y = POC_AU$cabY,
                                  rangeun = NULL, 
                                  groupPlot = NULL, 
                                  facetPlot = NULL,
                                  observationData = NULL, 
                                  simulateYear = NULL) + theme_bw() +
    ggtitle("Number of total population") + 
    geom_point(aes(
      y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Prison", "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), 
      colour = "black", shape = 5, size = 2) + 
    geom_segment(
      aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Prison", "low"], 
          yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Prison", "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1))) + 
    scale_y_continuous(limits = c(0,60000), breaks = seq(0, 60000, 10000))
  
  #####plot: number of subpop   ####
  # sub pop full names 
  Namelab <- c("PWID in community", "former PWID in community",
               "PWID in prison", "former PWID in prison", "non-PWID in prison")
  
  subpop_N_plot <- list() 
  
  for(x in seq_along(names(subpop_N))){ 
    
    subpop_N_plot[[x]] <- indicatorPlot(POC_AU, subpop_N[[x]] , 
                                        ylabel = "Number",
                                        xlimits = c(POC_AU$startYear, 
                                                    POC_AU$startYear+30, 5),
                                        calibration_Y = POC_AU$cabY,
                                        rangeun = NULL, 
                                        groupPlot = NULL, 
                                        facetPlot = NULL,
                                        observationData = NULL, 
                                        simulateYear = NULL) +
      theme_bw() + ggtitle(paste0(Namelab[x])) 
    
  }
  
  subpop_N_plot[[1]] <- subpop_N_plot[[1]] + 
    scale_y_continuous(limits = c(0, 110000), breaks = seq(0, 110000, 5000)) + 
    geom_point(aes(
      y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[1], "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 5, size = 2) + 
    geom_segment(
      aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[1], "low"], 
          yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[1], "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  
  subpop_N_plot[[2]] <- subpop_N_plot[[2]] + 
    scale_y_continuous(limits = c(0, 700000), breaks = seq(0, 700000, 50000)) + 
    
    geom_point(aes(
      y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[2], "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 5, size = 2) + 
    geom_segment(
      aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[2], "low"], 
          yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[2], "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  
  subpop_N_plot[[3]] <- subpop_N_plot[[3]] + 
    scale_y_continuous(limits = c(0, 15000), breaks = seq(0, 15000, 1000)) + 
    geom_point(aes(
      y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[3], "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 5, size = 2) + 
    geom_segment(
      aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[3], "low"], 
          yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[3], "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  
  
  subpop_N_plot[[4]] <- subpop_N_plot[[4]] + 
    scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 20000, 2000)) + 
    geom_point(aes(
      y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[4], "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 5, size = 2) + 
    geom_segment(
      aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[4], "low"], 
          yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[4], "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  
  subpop_N_plot[[5]] <- subpop_N_plot[[5]] + 
    scale_y_continuous(limits = c(0, 25000), breaks = seq(0, 25000, 2500)) + 
    geom_point(aes(
      y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[5], "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 5, size = 2) + 
    geom_segment(
      aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[5], "low"], 
          yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[5], "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  
  
  
  #####plot: % of PWID/fPWID in community/prison  ####
  # P_PWID 
  frac_PPWID_plot <- indicatorPlot(POC_AU, prison_proP , 
                                   ylabel = "Percentage (%)",
                                   xlimits = c(POC_AU$startYear, 
                                               POC_AU$startYear+30, 5),
                                   calibration_Y = POC_AU$cabY,
                                   rangeun = NULL, 
                                   groupPlot = NULL, 
                                   facetPlot = NULL,
                                   observationData = NULL, 
                                   simulateYear = NULL) + 
    theme_bw() + 
    scale_y_continuous(limits = c(0,20), breaks = seq(0, 20,5)) + ggtitle("PWID in prison")
  frac_PPWID_plot <- frac_PPWID_plot + 
    geom_point(aes(
      y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_proP_fit", "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 18) +
    geom_segment(
      aes(y = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_proP_fit", "low"], 
          yend = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_proP_fit", "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  
  # P_fPWID
  frac_PfPWID_plot <- indicatorPlot(POC_AU, prison_profP , 
                                    ylabel = "Percentage (%)",
                                    xlimits = c(POC_AU$startYear, 
                                                POC_AU$startYear+30, 5),
                                    calibration_Y = POC_AU$cabY,
                                    rangeun = NULL, 
                                    groupPlot = NULL, 
                                    facetPlot = NULL,
                                    observationData = NULL, 
                                    simulateYear = NULL) + 
    theme_bw() + 
    scale_y_continuous(limits = c(0,50), breaks = seq(0, 50,5)) + ggtitle("Former PWID in prison") + 
    geom_point(aes(
      y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_profP_fit", "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black",shape = 18) + 
    geom_segment(
      aes(y = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_profP_fit", "low"], 
          yend = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_profP_fit", "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  
  # C_PWID 
  frac_CPWID_plot <- indicatorPlot(POC_AU, commu_proP , 
                                   ylabel = "Percentage (%)",
                                   xlimits = c(POC_AU$startYear, 
                                               POC_AU$startYear+30, 5),
                                   calibration_Y = POC_AU$cabY,
                                   rangeun = NULL, 
                                   groupPlot = NULL, 
                                   facetPlot = NULL,
                                   observationData = NULL, 
                                   simulateYear = NULL) + 
    theme_bw() + 
    scale_y_continuous(limits = c(0,20), breaks = seq(0, 20,5)) + ggtitle("PWID in community")
  
  frac_CPWID_plot <- frac_CPWID_plot + geom_point(aes(
    y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "commu_proP_fit", "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape= 5, size = 0.6) + 
    geom_segment(
      aes(y = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "commu_proP_fit", "low"], 
          yend = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "commu_proP_fit", "up"], 
          x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  
  
  #####plot: annual number entry to/release from prisons #### 
  New_incar_plot <-  indicatorPlot(POC_AU, incar_total,
                                   ylabel = "Number",
                                   xlimits = c(POC_AU$cabY,
                                               POC_AU$cabY+30, 5),
                                   calibration_Y = POC_AU$cabY,
                                   rangeu = NULL,
                                   groupPlot = NULL,
                                   facetPlot = NULL,
                                   observationData = NULL, 
                                   simulateYear = NULL) +
    ggtitle("Annual number entry prisons") + theme_bw() + 
    scale_y_continuous(limits = c(0, 65000), breaks = seq(0, 65000, 5000)) + 
    geom_point(aes(
      y=datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "incar", "realpop"], 
      x = POC_AU$simY ), colour = "black")
  
  New_release_plot <-  indicatorPlot(POC_AU, release_total,
                                     ylabel = "Number",
                                     xlimits = c(POC_AU$cabY,
                                                 POC_AU$cabY+30, 5),
                                     calibration_Y = POC_AU$cabY,
                                     rangeu = NULL,
                                     groupPlot = NULL,
                                     facetPlot = NULL,
                                     observationData = NULL, 
                                     simulateYear = NULL) +
    ggtitle("Annual number release from prisons") + theme_bw() + 
    scale_y_continuous(limits = c(0, 65000), breaks = seq(0, 65000, 5000)) + 
    geom_point(aes(
      y=datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "release", "realpop"], 
      x = POC_AU$simY ), colour = "black")
  
  #####plot: annual number stopping/relapse injection ####
  # stopping injection
  New_stopinj_plot <-  indicatorPlot(POC_AU, inj_stop_bind,
                                     ylabel = "Number",
                                     xlimits = c(POC_AU$cabY,
                                                 POC_AU$cabY+30, 5),
                                     calibration_Y = POC_AU$cabY,
                                     rangeu = NULL,
                                     groupPlot = NULL,
                                     facetPlot = population,
                                     observationData = NULL, 
                                     simulateYear = NULL) +
    ggtitle("Annual number stop injection") + theme_bw() + 
    scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 20000, 5000)) 
  # relapse injection
  New_relapinj_plot <-  indicatorPlot(POC_AU, inj_relap_bind,
                                      ylabel = "Number",
                                      xlimits = c(POC_AU$cabY,
                                                  POC_AU$cabY+30, 5),
                                      calibration_Y = POC_AU$cabY,
                                      rangeu = NULL,
                                      groupPlot = NULL,
                                      facetPlot = population,
                                      observationData = NULL, 
                                      simulateYear = NULL) +
    ggtitle("Annual number relapse injection") + theme_bw() + 
    scale_y_continuous(limits = c(0,2500), breaks = seq(0, 2500, 500)) 
  
  ##============================================================================##
  
  
  
  #####plot: incidence of entry to/release from prisons #### 
  Inc_incar_plot <- indicatorPlot(POC_AU,incar_inc_bind ,
                                  ylabel = "Incidence",
                                  xlimits = c(POC_AU$cabY,
                                              POC_AU$cabY+30, 5),
                                  calibration_Y = POC_AU$cabY,
                                  rangeu = NULL,
                                  groupPlot = NULL,
                                  facetPlot = population,
                                  observationData = NULL, 
                                  simulateYear = NULL) +
    ggtitle("Incidence of incarceration") + theme_bw() + 
    scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 1)) 
  
  Inc_release_plot <-  indicatorPlot(POC_AU, release_inc_bind,
                                     ylabel = "Incidence",
                                     xlimits = c(POC_AU$cabY,
                                                 POC_AU$cabY+30, 5),
                                     calibration_Y = POC_AU$cabY,
                                     rangeu = NULL,
                                     groupPlot = NULL,
                                     facetPlot = population,
                                     observationData = NULL, 
                                     simulateYear = NULL) +
    ggtitle("Incidence of released") + theme_bw() + 
    scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 10)) 
  
  #####plot: incidence of stopping/relapse injection ####
  # stopping injection
  Inc_stopinj_plot <- list()
  
  
  Inc_stopinj_plot <- lapply(inj_stop_inc, function(x){ 
    
    a <- indicatorPlot(POC_AU, x,
                       ylabel = "Incidence",
                       xlimits = c(POC_AU$startYear ,
                                   POC_AU$startYear+30, 5),
                       calibration_Y = POC_AU$cabY,
                       rangeu = NULL,
                       groupPlot = NULL,
                       facetPlot = NULL,
                       observationData = NULL, 
                       simulateYear = NULL) +
      ggtitle("Incidence of stop injection") + theme_bw() + 
      scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 5))
    
    
    
    return(a)
    
  })
  
  Inc_stopinj_plot[["community"]] <-  Inc_stopinj_plot[["community"]] + 
    geom_point(aes(
      y=datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "commu_stopinj", "realpop"], 
      x = POC_AU$simY - POC_AU$cabY + 1 ), 
      colour = "black") + 
    geom_segment(
      aes(y = datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "commu_stopinj", "low"], 
          yend = datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "commu_stopinj", "up"], 
          x = POC_AU$simY - POC_AU$cabY + 1, xend = POC_AU$simY - POC_AU$cabY + 1 )) + 
    scale_y_continuous(limits = c(0,20), breaks = seq(0, 20, 1))
  
  
  
  # relapse injection
  Inc_relapinj_plot <-  indicatorPlot(POC_AU, inj_relap_inc_bind,
                                      ylabel = "Incidence",
                                      xlimits = c(POC_AU$cabY,
                                                  POC_AU$cabY+30, 5),
                                      calibration_Y = POC_AU$cabY,
                                      rangeu = NULL,
                                      groupPlot = NULL,
                                      facetPlot = population,
                                      observationData = NULL, 
                                      simulateYear = NULL) +
    ggtitle("Incidence of relapse injection") + theme_bw() + 
    scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 1))
  
  endY_plot <- 2030- POC_AU$cabY
  
  #### Results: prevalence & incidence #### 
  # Structure 
  #####Prevalence#####
  # (A)antibody prevalence 
  #    (I) each subpop
  #   (II) community
  #   (III) Prison
  #   (IV) Prison_f/PWID (P_fPWID, P_PWID)
  # (B)RNA prevalence 
  #    (I) each subpop
  #   (II) community
  #   (III) Prison
  #   (IV) Prison_f/PWID (P_fPWID, P_PWID)
  
  ##### (A)antibody prevalence ##### 
  #    (I) each subpop
  # number of susceptible 
  pop_state <- popResults_MidYear(POC_AU, calibrateInit,
                                  Population = POC_AU$popNames,
                                  Disease_prog = POC_AU$progress_name, 
                                  Cascade = POC_AU$cascade_name, param = NULL, 
                                  endYear = endY)%>%ungroup()%>%
    as_tibble()
  
  tempNOTInfected_subpop <- pop_state%>%filter(disease_prog!= "a")%>%
    filter(state == "s")%>%group_by(year, population)%>%
    summarise(best = sum(best))%>%arrange(year, population)
  
  tempChronic_subpop <- pop_state%>%filter(disease_prog!= "a")%>%
    group_by(year, population)%>%
    summarise(best = sum(best))%>%arrange(year, population)
  
  # arrange order to align with other dts 
  pop_N <- pop_N%>%arrange(year, population)
  
  
  tempPrev_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                           population = POC_AU$popNames,
                           
                           as.data.frame(100*(tempChronic_subpop[, -c(1,2)] - 
                                                tempNOTInfected_subpop[ ,-c(1,2)])/ 
                                           tempChronic_subpop[ ,-c(1,2)]))%>%
    tibble::as_tibble()  
  
  
  #prevalence in setting  
  # tempPrev_setting[["setting]]
  
  # community 
  commu_N <- commu_N%>%arrange(year)
  
  tempNOTInfected_commu <- popResults_MidYear(POC_AU, calibrateInit,
                                              Population = c("C_PWID", "C_fPWID"),
                                              Disease_prog = NULL, 
                                              Cascade = c("s"), 
                                              param = NULL ,endYear = endY)%>%
    as_tibble()%>%group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  tempPrev_setting <- list()
  
  tempPrev_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                       as.data.frame(100*(commu_N[, -c(1)] - 
                                                            tempNOTInfected_commu[ ,-c(1)])/ 
                                                       commu_N[ ,-c(1)]))%>%tibble::as_tibble()
  
  # prison 
  ## PWID + former PWID + nonPWID 
  prison_N <- prison_N%>%arrange(year)
  
  tempNOTInfected_prison <- popResults_MidYear(POC_AU, calibrateInit,
                                               Population = c("P_PWID", "P_fPWID", 
                                                              "P_nPWID"),
                                               Disease_prog = NULL, 
                                               Cascade = c("s"), 
                                               param = NULL ,endYear = endY)%>%
    as_tibble()%>%group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  
  tempPrev_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                         as.data.frame(100*(prison_N[, -c(1)] - 
                                                              tempNOTInfected_prison[ ,-c(1)])/ 
                                                         prison_N[ ,-c(1)]))%>%tibble::as_tibble()
  
  
  # prison_PWID experienced  
  ## PWID + former PWID + nonPWID 
  prisonPWID_N <- prisonPWID_N%>%arrange(year)
  
  tempNOTInfected_prisonPWID <- popResults_MidYear(POC_AU, calibrateInit,
                                                   Population = c("P_PWID", "P_fPWID"),
                                                   Disease_prog = NULL, 
                                                   Cascade = c("s"), 
                                                   param = NULL ,endYear = endY)%>%
    as_tibble()%>%group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  
  tempPrev_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                             as.data.frame(100*(prisonPWID_N[, -c(1)] - 
                                                                  tempNOTInfected_prisonPWID[ ,-c(1)])/ 
                                                             prisonPWID_N[, -c(1)]))%>%tibble::as_tibble()
  
  
  
  ##### (B)RNA prevalence ##### 
  #    (I) each subpop
  #   (II) community
  #   (III) Prison
  #   (IV) Prison_f/PWID (P_fPWID, P_PWID)
  tempNOTInfectedRNA_subpop <- popResults_MidYear(POC_AU, calibrateInit,Population = POC_AU$popNames,
                                                  Disease_prog = NULL, 
                                                  Cascade = c("s", "cured"), 
                                                  param = NULL ,endYear = endY)%>%
    as_tibble()%>%group_by(year, population)%>%
    summarise(best = sum(best))%>%arrange(year, population)
  
  
  
  tempPrevRNA_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                              population = POC_AU$popNames,
                              
                              as.data.frame(100*(pop_N[, -c(1,2)] - 
                                                   tempNOTInfectedRNA_subpop[ ,-c(1,2)])/ 
                                              pop_N[ ,-c(1,2)]))%>%
    tibble::as_tibble()  
  
  
  #prevalence in setting  
  # tempPrev_setting[["setting]]
  
  
  tempNOTInfectedRNA_commu <- popResults_MidYear(POC_AU, calibrateInit,
                                                 Population = c("C_PWID", "C_fPWID"),
                                                 Disease_prog = NULL, 
                                                 Cascade = c("s", "cured"), 
                                                 param = NULL ,endYear = endY)%>%
    as_tibble()%>%group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  tempPrevRNA_setting <- list()
  
  tempPrevRNA_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                          as.data.frame(100*(commu_N[, -c(1)] - 
                                                               tempNOTInfectedRNA_commu[ ,-c(1)])/ 
                                                          commu_N[ ,-c(1)]))%>%tibble::as_tibble()
  
  # prison 
  ## PWID + former PWID + nonPWID 
  
  
  tempNOTInfectedRNA_prison <- popResults_MidYear(POC_AU, calibrateInit,
                                                  Population = c("P_PWID", "P_fPWID", 
                                                                 "P_nPWID"),
                                                  Disease_prog = NULL, 
                                                  Cascade = c("s", "cured"), 
                                                  param = NULL ,endYear = endY)%>%
    as_tibble()%>%group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  
  tempPrevRNA_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                            as.data.frame(100*(prison_N[, -c(1)] - 
                                                                 tempNOTInfectedRNA_prison[ ,-c(1)])/ 
                                                            prison_N[ ,-c(1)]))%>%tibble::as_tibble()
  
  
  # prison_PWID experienced  
  ## PWID + former PWID + nonPWID 
  
  
  tempNOTInfectedRNA_prisonPWID <- popResults_MidYear(POC_AU, calibrateInit,
                                                      Population = c("P_PWID", "P_fPWID"),
                                                      Disease_prog = NULL, 
                                                      Cascade = c("s", "cured"), 
                                                      param = NULL ,endYear = endY)%>%
    as_tibble()%>%group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  
  tempPrevRNA_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                                as.data.frame(100*(prisonPWID_N[, -c(1)] - 
                                                                     tempNOTInfectedRNA_prisonPWID[ ,-c(1)])/ 
                                                                prison_N[ ,-c(1)]))%>%tibble::as_tibble()
  
  
  
  ##### HCV incidence ##### 
  HCVInfect_subpop <- indicatorResults(POC_AU, calibrateInit, "newInfections", 
                                       pop=POC_AU$popNames,
                                       paramR = NULL, range = NULL,
                                       endY = endY)
  
  HCVInc_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                    each = POC_AU$npops),
                         population = POC_AU$popNames,
                         as.data.frame(HCVInfect_subpop[, -c(1,2)] / 
                                         pop_N[ ,-c(1,2)]*100))%>%
    tibble::as_tibble() 
  
  
  
  HCVInfect_setting <- list()
  
  HCVInfect_setting[["commu"]] <- indicatorResults(POC_AU, 
                                                   calibrateInit, 
                                                   "newInfections", 
                                                   pop= c("C_PWID", "C_fPWID"),
                                                   paramR = NULL, range = NULL,
                                                   endY = endY)%>%group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  HCVInfect_setting[["prisons"]] <- indicatorResults(POC_AU, 
                                                     calibrateInit, 
                                                     "newInfections", 
                                                     pop= c("P_PWID", "P_fPWID",
                                                            "P_nPWID"), 
                                                     paramR = NULL, 
                                                     range = NULL,
                                                     endY = endY)%>%group_by(year)%>%
    summarise(best = sum(best))%>%arrange(year)
  
  HCVInc_setting <- list()
  
  HCVInc_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(100*(HCVInfect_setting[["commu"]][ , -c(1)]/ 
                                                          commu_N[ ,-c(1)])))%>%
    tibble::as_tibble()
  HCVInc_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                       as.data.frame(100*(HCVInfect_setting[["prisons"]][ , -c(1)]/ 
                                                            prison_N[ ,-c(1)])))%>%
    tibble::as_tibble()
  
  
  ####Plot: incidence and prevalence  ####
  #####plot: prevalence #####
  pop_labname <- c("PWID in community",  "Former PWID in community", 
                   "PWID in prisons",  "Former PWID in prisons", 
                   "nonPWID in prisons")
  HCVPrev <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrev_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = HCV.seroprevalence*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = POC_AU$popNames, 
                                                 labels = pop_labname ))
  
  tempPrev_subpop <- tempPrev_subpop%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  
  popPrevPlot <- indicatorPlot(POC_AU, tempPrev_subpop, 
                               ylabel = "HCV prevalence (%)",
                               xlimits = c(POC_AU$startYear, 
                                           (POC_AU$startYear+endY_plot), 5),
                               calibration_Y = POC_AU$cabY,
                               rangeun = NULL, 
                               groupPlot = NULL, 
                               facetPlot = population,
                               observationData = HCVPrev, 
                               simulateYear = POC_AU$simY) +
    ggtitle("HCV prevalence by population") 
  
  popPrevPlot <- popPrevPlot + 
    facet_custom (~population,
                  scales = "free", ncol = 3,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 80))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 80))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, 80))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, 10)))
                    )) + theme_bw()
  
  popPrevPlot 
  # setting 
  tempPrev_setting_bind <- dplyr::bind_rows(tempPrev_setting, .id = 'population')%>%
    mutate(population = factor(population, 
                               levels = c("commu", "prisons", "prisonsPWID"), 
                               labels = c("Community", "Prisons", 
                                          "Current & former PWID in prisons")))%>%
    as.data.frame()
  
  HCVPrev_setting <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrev_setting_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = HCVPrev*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = c("commu", "prisons"), 
                                                 labels = c("Community", "Prisons")))
  
  
  settingPrevPlot <- indicatorPlot(POC_AU, tempPrev_setting_bind, 
                                   ylabel = "HCV prevalence (%)",
                                   xlimits = c(POC_AU$startYear, 
                                               (POC_AU$startYear+endY_plot - 1), 5),
                                   calibration_Y = POC_AU$cabY,
                                   rangeun = NULL, 
                                   groupPlot = NULL, 
                                   facetPlot = population,
                                   observationData = HCVPrev_setting, 
                                   simulateYear = POC_AU$simY) +
    ggtitle("HCV prevalence by setting") 
  
  settingPrevPlot <- settingPrevPlot + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 35))),
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 100)))
                    )) + theme_bw() 
  
  ##### B: RNA prevalence #####
  tempPrevRNA_subpop <- tempPrevRNA_subpop%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  HCVPrevRNA <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrevRNA_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = HCV.RNA.prevalence*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = POC_AU$popNames, 
                                                 labels = pop_labname ))
  popPrevRNAPlot <- indicatorPlot(POC_AU, tempPrevRNA_subpop, 
                                  ylabel = "HCV prevalence (%)",
                                  xlimits = c(POC_AU$startYear, 
                                              (POC_AU$startYear+endY_plot), 5),
                                  calibration_Y = POC_AU$cabY,
                                  rangeun = NULL, 
                                  groupPlot = NULL, 
                                  facetPlot = population,
                                  observationData = HCVPrevRNA , 
                                  simulateYear = POC_AU$simY) +
    ggtitle("HCV prevalence by population") 
  
  popPrevRNAPlot <- popPrevRNAPlot + 
    facet_custom (~population,
                  scales = "free", ncol = 3,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 50))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 50))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 50))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, 50))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, 2)))
                    )) + theme_bw()
  
  popPrevRNAPlot 
  # setting 
  tempPrevRNA_setting_bind <- dplyr::bind_rows(tempPrevRNA_setting, .id = 'population')%>%
    mutate(population = factor(population, 
                               levels = c("commu", "prisons", "prisonsPWID"), 
                               labels = c("Community", "Prisons", 
                                          "Current & former PWID in prisons")))%>%
    as.data.frame()
  
  HCVPrevRNAsetting <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrevRNA_setting_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = HCVPrev*100,
                             low = lower*100,
                             up = upper*100,
                             population = factor(population, 
                                                 levels = c("commu", "prisons"), 
                                                 labels = c("Community", "Prisons" )))
  
  
  settingPrevRNAPlot <- indicatorPlot(POC_AU, tempPrevRNA_setting_bind, 
                                      ylabel = "HCV prevalence (%)",
                                      xlimits = c(POC_AU$startYear, 
                                                  (POC_AU$startYear+endY_plot - 1), 5),
                                      calibration_Y = POC_AU$cabY,
                                      rangeun = NULL, 
                                      groupPlot = NULL, 
                                      facetPlot = population,
                                      observationData = HCVPrevRNAsetting, 
                                      simulateYear = POC_AU$simY) +
    ggtitle("HCV prevalence by setting") 
  
  settingPrevRNAPlot <- settingPrevRNAPlot + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 50))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 50))),
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 50)))
                    )) + theme_bw()
  
  
  settingPrevRNAPlot
  #####plot: incidence #####
  HCVInc_subpop <- HCVInc_subpop%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname )) 
  HCVInc <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVInc_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = HCVInc*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = POC_AU$popNames, 
                                                 labels = pop_labname ))
  
  popIncPlot <- indicatorPlot(POC_AU, HCVInc_subpop, 
                              ylabel = "HCV incidence",
                              xlimits = c(POC_AU$startYear, 
                                          (POC_AU$startYear+endY_plot), 5),
                              calibration_Y = POC_AU$cabY,
                              rangeun = NULL, 
                              groupPlot = NULL, 
                              facetPlot = population,
                              observationData = HCVInc, 
                              simulateYear = POC_AU$simY) +
    ggtitle("HCV incidence by population") 
  
  popIncPlot <- popIncPlot + 
    facet_custom (~population,
                  scales = "free", ncol = 3,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 20))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 20))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 20))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, 20))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, 10)))
                    )) + theme_bw()
  
  # setting 
  HCVInc_setting_bind <- dplyr::bind_rows(HCVInc_setting, .id = 'population')%>%
    mutate(population = factor(population, 
                               levels = c("commu", "prisons"), 
                               labels = c("Community", "Prisons")))%>%
    as.data.frame()
  
  HCVInc_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVInc_setting_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = HCVInc*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = c("commu", "prisons"), 
                                                 labels = c("Community", "Prisons")))
  
  settingIncPlot <- indicatorPlot(POC_AU, HCVInc_setting_bind, 
                                  ylabel = "HCV incidence",
                                  xlimits = c(POC_AU$startYear, 
                                              (POC_AU$startYear+endY_plot - 1), 5),
                                  calibration_Y = POC_AU$cabY,
                                  rangeun = NULL, 
                                  groupPlot = NULL, 
                                  facetPlot = population,
                                  observationData = HCVInc_setting_fit, 
                                  simulateYear = POC_AU$simY) +
    ggtitle("HCV incidence by setting") 
  
  settingIncPlot <- settingIncPlot + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 20))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 10))))) + theme_bw()
  
  ####Results: death ####
  #### HCV death ####
  
  HCVdeath <- indicatorResults(POC_AU, calibrateInit, "newHCVdeaths", 
                               pop=POC_AU$popNames,
                               paramR = NULL, range = NULL,
                               endY = endY)%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  
  death <- indicatorResults(POC_AU, calibrateInit, "newDeath", 
                            pop=POC_AU$popNames,
                            paramR = NULL, range = NULL,
                            endY = endY)%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  
  
  HCVdeath_plot <- indicatorPlot(POC_AU, HCVdeath, 
                                 ylabel = "HCV-related deaths",
                                 xlimits = c(POC_AU$startYear, 
                                             (POC_AU$startYear+endY_plot), 5),
                                 calibration_Y = POC_AU$cabY,
                                 rangeun = NULL, 
                                 groupPlot = NULL, 
                                 facetPlot = population,
                                 observationData = NULL, 
                                 simulateYear = POC_AU$simY) +
    ggtitle("HCV-related deaths") + theme_bw()
  
  
  HCVdeath_plot <- HCVdeath_plot + 
    facet_custom(~population,
                 scales = "free", ncol = 3,
                 scale_overrides = 
                   list(
                     scale_new(1,
                               scale_y_continuous(limits = 
                                                    c(0, 1000))),
                     scale_new(2,
                               scale_y_continuous(limits = 
                                                    c(0, 10000))),
                     scale_new(3,
                               scale_y_continuous(limits = 
                                                    c(0, 100))),
                     scale_new(4,
                               scale_y_continuous(limits = 
                                                    c(0, 1000))),
                     scale_new(5,
                               scale_y_continuous(limits = 
                                                    c(0, 5))))) 
  
  death_plot <- indicatorPlot(POC_AU, death, 
                              ylabel = "Deaths",
                              xlimits = c(POC_AU$startYear, 
                                          (POC_AU$startYear+endY_plot), 5),
                              calibration_Y = POC_AU$cabY,
                              rangeun = NULL, 
                              groupPlot = NULL, 
                              facetPlot = population,
                              observationData = NULL, 
                              simulateYear = POC_AU$simY) +
    ggtitle("Deaths") 
  
  
  death_plot <- death_plot + 
    facet_custom(~population,
                 scales = "free", ncol = 3,
                 scale_overrides = 
                   list(
                     scale_new(1,
                               scale_y_continuous(limits = 
                                                    c(0, 2000))),
                     scale_new(2,
                               scale_y_continuous(limits = 
                                                    c(0, 10000))),
                     scale_new(3,
                               scale_y_continuous(limits = 
                                                    c(0, 1000))),
                     scale_new(4,
                               scale_y_continuous(limits = 
                                                    c(0, 500))),
                     scale_new(5,
                               scale_y_continuous(limits = 
                                                    c(0, 5))))) 
  
  #### Results: all component ####
  pop_state <- popResults_MidYear(POC_AU, calibrateInit,
                                  Population = POC_AU$popNames,
                                  Disease_prog = POC_AU$progress_name, 
                                  Cascade = POC_AU$cascade_name, param = NULL, 
                                  endYear = endY)%>%ungroup()
  
  pop_stateAll <- pop_state%>%group_by(year, state)%>%
    summarise_at("best",sum)
  
  pop_stateSub <- pop_state%>%group_by(year, population, state)%>%
    summarise_at("best", sum)
  
  # plot 
  Ppop_stateAll <-  ggplot(data = pop_stateAll, aes(x = year, y = best)) + 
    geom_area() + 
    facet_wrap(.~ state, scale = "free") + 
    labs(y = "Number") + 
    scale_x_continuous( 
      limits = c(1, 31), 
      breaks = seq(1,31, 5),
      labels = seq(POC_AU$cabY,POC_AU$cabY + 30, 5)) 
  
  #### Results: disease progression ####
  
  pop_dispro <- popResults_MidYear(POC_AU, calibrateInit,
                                   Population = POC_AU$popNames,
                                   Disease_prog = POC_AU$progress_name, 
                                   Cascade = NULL, param = NULL, 
                                   endYear = endY)%>%ungroup()
  
  pop_disproAll <- pop_dispro%>%group_by(year, disease_prog)%>%
    summarise_at("best",sum)
  
  pop_disproSub <- pop_dispro%>%group_by(year, population, disease_prog)%>%
    summarise_at("best", sum)
  
  # tidy up the code for presenting plot 
  #  Area plot 
  ## disease_prog by order 
  
  Ppop_disproAll <- ggplot(data = pop_disproAll, aes(x = year, y = best,group)) + 
    geom_area(aes(fill = disease_prog)) + 
    scale_x_continuous( 
      limits = c(1, 31), 
      breaks = seq(1,31, 5),
      labels = seq(POC_AU$cabY,POC_AU$cabY + 30, 5))
  
  #### Results: Cascade ####
  pop_casca <- popResults_MidYear(POC_AU, calibrateInit,
                                  Population = POC_AU$popNames,
                                  Disease_prog = NULL, 
                                  Cascade = POC_AU$cascade_name, 
                                  param = NULL, 
                                  endYear = endY)%>%ungroup()
  
  pop_cascaAll <- pop_casca%>%group_by(year, cascade)%>%
    summarise_at("best",sum)
  
  pop_cascaSub <- pop_casca%>%group_by(year, population, cascade)%>%
    summarise_at("best", sum)
  
  Ppop_cascaAll <- ggplot(data = pop_cascaAll, aes(x = year, y = best,group)) + 
    geom_area(aes(fill = cascade)) + 
    scale_x_continuous( 
      limits = c(1, 31), 
      breaks = seq(1,31, 5),
      labels = seq(POC_AU$cabY,POC_AU$cabY + 30, 5))
  
  
  #### flow
  # select flow outcomes from the output list 
  calibrateFlow <- calibrateInit[!names(calibrateInit)%in%
                                   c("allPops", "newpop_tran", "HCVdeathState",
                                     "newDeathState", "death_hcv")]
  
  flow_sub <- list()
  flow_all <- list()
  
  flow_sub <- lapply(names(calibrateFlow), function(x){ 
    a <- indicatorResults(POC_AU, calibrateFlow, x, 
                          pop=POC_AU$popNames,
                          paramR = NULL, range = NULL,
                          endY = endY)
  })
  
  names(flow_sub) <- names(calibrateFlow)
  
  flow_all <- lapply(flow_sub, function(x){ 
    
    a <- x%>%group_by(year)%>%summarise_at("best", sum)
  }) 
  
  names(flow_all) <- names(calibrateFlow)
  
  # plot for flow 
  pflow_all <- list()
  pflow_all <- lapply(names(flow_all), function(x){ 
    
    a <- indicatorPlot(POC_AU, flow_all[[x]], 
                       ylabel = "Numbers",
                       xlimits = c(POC_AU$startYear, 
                                   (POC_AU$startYear+endY_plot), 5),
                       calibration_Y = POC_AU$cabY,
                       rangeun = NULL, 
                       groupPlot = NULL, 
                       facetPlot = NULL,
                       observationData = NULL, 
                       simulateYear = POC_AU$simY) +
      ggtitle(x) + theme_Publication()
    
    a <- a + expandy_man(flow_all[[x]], POC_AU$startYear+endY_plot)
    
    return(a)
  })
  
  names(pflow_all) <- names(flow_all)
  
  
  pflow_sub <- lapply(names(flow_sub), function(x){ 
    
    a <- indicatorPlot(POC_AU, flow_sub[[x]], 
                       ylabel = "Numbers",
                       xlimits = c(POC_AU$startYear, 
                                   (POC_AU$startYear+endY_plot), 5),
                       calibration_Y = POC_AU$cabY,
                       rangeun = NULL, 
                       groupPlot = NULL, 
                       facetPlot = population,
                       observationData = NULL, 
                       simulateYear = POC_AU$simY) +
      ggtitle(x)
    
  })
  
  names(pflow_sub) <- names(flow_sub)
  
  #cascade ab
  ## calibration datasets 
  HCVab <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVAbtest_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = abCoverage*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = POC_AU$popNames, 
                                                 labels = pop_labname ))
  
  HCVrna <- read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVRNAtest_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = RNACoverage*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = POC_AU$popNames, 
                                                 labels = pop_labname ))
  
  HCVrnaonly <- read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVRNAonlytest_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = RNAonlyCoverage*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = POC_AU$popNames, 
                                                 labels = pop_labname ))
  
  
  
  
  # newtestingab/chronic infected + cured from chronic stage 
  testab <- pop_stateSub%>%
    filter(state%in% c("f0_undiag", "f1_undiag", "f2_undiag", "f3_undiag",
                       "f4_undiag", "dc_undiag", "hcc_undiag", "lt_undiag",
                       "plt_undiag"))%>%
    group_by(year,population)%>%
    summarize(best = sum(best))
  
  
  
  
  Tab_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                               each = POC_AU$npops),
                    population = POC_AU$popNames,
                    as.data.frame(100*flow_sub$newTestingAb[, -c(1,2)] / 
                                    testab[ ,-c(1,2)]))%>%
    tibble::as_tibble()%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  
  
  Tpoct_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                 each = POC_AU$npops),
                      population = POC_AU$popNames,
                      as.data.frame(100*flow_sub$newTestingPOCT[, -c(1,2)] / 
                                      testab[ ,-c(1,2)]))%>%
    tibble::as_tibble()%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  
  testrna <- pop_stateSub%>%
    filter(state%in% c("f0_diag_ab", "f1_diag_ab", "f2_diag_ab", "f3_diag_ab",
                       "f4_diag_ab", "dc_diag_ab", "hcc_diag_ab", "lt_diag_ab",
                       "plt_diag_ab"))%>%
    group_by(year,population)%>%
    summarize(best = sum(best))
  
  Trna_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                each = POC_AU$npops),
                     population = POC_AU$popNames,
                     as.data.frame(100*flow_sub$newTestingAg[, -c(1,2)] / 
                                     testrna[ ,-c(1,2)]))%>%
    tibble::as_tibble()%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  
  treatint <- pop_stateSub%>%
    filter(state%in% c("a_diag_RNA", "f0_diag_RNA", "f1_diag_RNA", "f2_diag_RNA", "f3_diag_RNA",
                       "f4_diag_RNA", "dc_diag_RNA", "hcc_diag_RNA", "lt_diag_RNA",
                       "plt_diag_RNA"))%>%
    group_by(year,population)%>%
    summarize(best = sum(best))
  
  Tint_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                each = POC_AU$npops),
                     population = POC_AU$popNames,
                     as.data.frame(100*flow_sub$newTreatment[, -c(1,2)] / 
                                     treatint[ ,-c(1,2)]))%>%
    tibble::as_tibble()%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  
  Tab_frag_p <- indicatorPlot(POC_AU, Tab_frag, 
                              ylabel = "%",
                              xlimits = c(POC_AU$startYear, 
                                          (POC_AU$startYear+endY_plot), 5),
                              calibration_Y = POC_AU$cabY,
                              rangeun = NULL, 
                              groupPlot = NULL, 
                              facetPlot = population,
                              observationData = HCVab, 
                              simulateYear = POC_AU$simY) + 
    ggtitle("antibody testing coverage") 
  
  Tab_frag_p <- Tab_frag_p + 
    facet_custom (~population,
                  scales = "free", ncol = 3,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 40))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 40))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 80))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, 40))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, 20)))
                    )) + theme_bw()
  
  Tpoct_frag_p <- indicatorPlot(POC_AU, Tpoct_frag, 
                                ylabel = "%",
                                xlimits = c(POC_AU$startYear, 
                                            (POC_AU$startYear+endY_plot), 5),
                                calibration_Y = POC_AU$cabY,
                                rangeun = NULL, 
                                groupPlot = NULL, 
                                facetPlot = population,
                                observationData = HCVrnaonly, 
                                simulateYear = POC_AU$simY) + 
    ggtitle("point-of-care RNA testing coverage")
  
  Tpoct_frag_p <- Tpoct_frag_p + 
    facet_custom (~population,
                  scales = "free", ncol = 3,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 40))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 40))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 30))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, 40))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, 20)))
                    )) + theme_bw()
  
  
  Trna_frag_p <- indicatorPlot(POC_AU, Trna_frag, 
                               ylabel = "%",
                               xlimits = c(POC_AU$startYear, 
                                           (POC_AU$startYear+endY_plot), 5),
                               calibration_Y = POC_AU$cabY,
                               rangeun = NULL, 
                               groupPlot = NULL, 
                               facetPlot = population,
                               observationData = HCVrna, 
                               simulateYear = POC_AU$simY) + 
    ggtitle("RNA testing coverage") 
  
  Trna_frag_p <- Trna_frag_p + 
    facet_custom (~population,
                  scales = "free", ncol = 3,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, 100)))
                    )) + theme_bw()
  
  Tint_frag_p <- indicatorPlot(POC_AU, Tint_frag, 
                               ylabel = "%",
                               xlimits = c(POC_AU$startYear, 
                                           (POC_AU$startYear+endY_plot), 5),
                               calibration_Y = POC_AU$cabY,
                               rangeun = NULL, 
                               groupPlot = NULL, 
                               facetPlot = population,
                               observationData = NULL, 
                               simulateYear = POC_AU$simY) + 
    ggtitle("Percentage of treatment init") 
  
  Tint_frag_p <- Tint_frag_p + 
    facet_custom (~population,
                  scales = "free", ncol = 3,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, 100)))
                    )) + theme_bw()
  
  
  ### treatment init 
  
  # lifetime treatment init in C_PWID 
  
  
  # treatment init annual numbers in prisons from 2016 to 2022 
  
  pflow_sub$newTreatment + theme_bw()
  
  # by setting 
  flow_setting <- lapply(flow_sub, function(x){ 
    
    a <- x%>%
      mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), 
                              "commu", "prison"))
    
    a <- a%>%group_by(year, setting)%>%summarise_at("best", sum)%>%
      mutate(population = setting)%>%select(-setting)
  }) 
  N_treatment <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                  each = 2),
                       population = flow_setting$newTreatment[ ,3],
                       as.data.frame(flow_setting$newTreatment[, -c(1,3)] + 
                                       flow_setting$newRetreat[, -c(1,3)] ))%>%
    tibble::as_tibble()
  
  N_treatment_setting_p <- indicatorPlot(POC_AU, N_treatment, 
                                         ylabel = "N",
                                         xlimits = c(POC_AU$startYear, 
                                                     (POC_AU$startYear+endY_plot), 5),
                                         calibration_Y = POC_AU$cabY,
                                         rangeun = NULL, 
                                         groupPlot = NULL, 
                                         facetPlot = population,
                                         observationData = NULL, 
                                         simulateYear = POC_AU$simY) + 
    ggtitle("Number of treatment init") 
  
  N_treatment_setting_p <- N_treatment_setting_p + 
    facet_custom (~population,
                  scales = "free", ncol = 2,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 60000))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 5000)))))
  N_treatment_setting_p
  
  
  # lifetime treatment initated 
  HCVtreatever <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatever_POC_AU.csv")), header = TRUE)%>%
    as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                             realPop = realpop*100,
                             up = upper*100,
                             low = lower*100,
                             population = factor(population, 
                                                 levels = POC_AU$popNames, 
                                                 labels = pop_labname ))
  
  
  
  Tever <- pop_state%>%
    filter(cascade %in% c("treat", "treat_f", "cured") & disease_prog!= "a")%>%
    group_by(year, population)%>%summarise(best = sum(best))
  HCVever <- pop_state%>%
    filter(!disease_prog%in% c ("s", "a"))%>%
    group_by(year, population)%>%summarise(best = sum(best))
  
  treatInit_ever <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                     each = POC_AU$npops),
                          population = POC_AU$popNames,
                          as.data.frame(100*Tever[, -c(1,2)] / 
                                          HCVever[ ,-c(1,2)]))%>%
    tibble::as_tibble()%>%
    mutate(population = factor(population, 
                               levels = POC_AU$popNames, 
                               labels = pop_labname ))
  
  treatInit_ever_p <- indicatorPlot(POC_AU, treatInit_ever, 
                                    ylabel = "%",
                                    xlimits = c(POC_AU$startYear, 
                                                (POC_AU$startYear+endY_plot), 5),
                                    calibration_Y = POC_AU$cabY,
                                    rangeun = NULL, 
                                    groupPlot = NULL, 
                                    facetPlot = population,
                                    observationData = HCVtreatever, 
                                    simulateYear = POC_AU$simY) + 
    ggtitle("Lifetime treatment initiated") 
  
  treatInit_ever_p <- treatInit_ever_p + 
    facet_custom (~population,
                  scales = "free", ncol = 3,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, 100))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, 100)))
                    )) + theme_bw()
  treatInit_ever_p
  # check 
  g1 <- marrangeGrob(grobs=c(list(subpop_N_plot[[1]], subpop_N_plot[[2]], 
                                  subpop_N_plot[[3]], subpop_N_plot[[4]],
                                  subpop_N_plot[[5]]),
                             list(popPrevPlot), 
                             list(popPrevRNAPlot),
                             list(settingPrevPlot),
                             list(settingPrevRNAPlot),
                             list(settingIncPlot),
                             list(popIncPlot),
                             list(Tab_frag_p),
                             list(Tpoct_frag_p),
                             list(Trna_frag_p),
                             list(Tint_frag_p),
                             list(treatInit_ever_p),
                             list(N_treatment_setting_p)
                             
  ), 
  nrow=1, ncol=1)
  ggsave(file=paste0("calibrationPlot20231023_",i, ".pdf"), g1, width = 15 , height = 10.2) 
  
  }




x

# show plots 
# ggpubr::ggarrange(plotlist = subpop_N_plot) 
# popPrevPlot
# popPrevRNAPlot
# settingPrevPlot
# settingPrevRNAPlot
# settingIncPlot
# popIncPlot
# Tab_frag_p
# Tpoct_frag_p
# Trna_frag_p
# Tint_frag_p
# treatInit_ever_p
# N_treatment_setting_p
