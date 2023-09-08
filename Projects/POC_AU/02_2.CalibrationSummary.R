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

source(file.path(Rcode, "/Functions/plotFunctions.R")) 

# setting end Year 

endY <- 100
#### tidy up datapoint #### 

datapoint <- list()
datapoint[["N"]] <- 
  cbind.data.frame(year = rep(POC_AU$simY, length(POC_AU$popNames) +3), 
                   indicator = c(POC_AU$popNames, "Community", "Prison", "Total"), 
                   realpop = c(80000, 400000, NA, NA, NA, 480000, 40000, 520000),
                   low= c(60000, 30000, NA, NA, NA, 360000, NA, NA),
                   up = c(100000, 600000, NA, NA, NA, 700000, NA, NA))

datapoint[["frac"]] <- 
  cbind.data.frame(year = rep(POC_AU$simY, 3), 
                   indicator = c("commu_proP_fit", "prison_proP_fit", 
                                 "prison_profP_fit"), 
                   realpop = c(16.67, 16.12, 35.88),
                   low= c(NA, NA, NA),
                   up = c(NA, NA, NA))

datapoint[["flow"]] <- 
  cbind.data.frame(year = rep(POC_AU$simY, 3), 
                   indicator = c("incar", "release", 
                                 "commu_stopinj"), 
                   realpop = c(63753, 63113, 5.4),
                   low= c(NA, NA, 4.3),
                   up = c(NA, NA, 6.7))



#### Result:population ####
##### N: each subpop/ total  #####
# subpop as list
subpop_N <- lapply(POC_AU$popNames, function(x){ 
  
  a <- popResults_MidYear(POC_AU, calibrateInit,
                          Population = x,
                          Disease_prog = NULL, 
                          Cascade = NULL, param = NULL, 
                          endYear = endY)%>%ungroup()
})

names(subpop_N) <- POC_AU$popNames

ggplot(data = as.data.frame(subpop_N[[5]]), aes(x = year, y = best)) + 
  geom_line() + scale_y_continuous(limits = c(20800, 21000))

# all subpop in one list 
pop_N <- dplyr::bind_rows(subpop_N, .id = 'population')

# total number 

total_N <- pop_N%>%group_by(year)%>%summarise(best = sum(best))

##### N: community #####
commu_N <- popResults_MidYear(POC_AU, calibrateInit,
                              Population = c("C_PWID", "C_fPWID"),
                              Disease_prog = NULL, 
                              Cascade = NULL, param = NULL, 
                              endYear = endY)%>%ungroup()%>%
  dplyr::group_by(year)%>%summarise_at("best",sum)

##### N: prison #####
prison_N <- popResults_MidYear(POC_AU, calibrateInit,
                               Population = c("P_PWID", "P_fPWID", "P_nPWID"),
                               Disease_prog = NULL, 
                               Cascade = NULL, param = NULL, 
                               endYear = endY)%>%ungroup()%>%
  dplyr::group_by(year)%>%summarise_at("best",sum)

##### %: PWID in community/prison #####
commu_proP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                    as.data.frame(pop_N[pop_N$population =="C_PWID",3] / commu_N[ ,-1])*100)%>%
  tibble::as_tibble()

prison_proP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                     as.data.frame(pop_N[pop_N$population =="P_PWID",3] / prison_N[ ,-1])*100)%>%
  tibble::as_tibble()

##### % fPWID in prison ##### 
prison_profP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                      as.data.frame(pop_N[pop_N$population =="P_fPWID",3] / prison_N[ ,-1])*100)%>%
  tibble::as_tibble()

##### number of annual leaving in each subpop #####
# this is for getting annual number of released in non-PWID in prison 
rel <- indicatorResult_uno(POC_AU, calibrateInit, "newLeave",
                           populations = POC_AU$popNames, endYear= endY)%>%
  mutate(year = year + POC_AU$cabY - 1) 

##### number of annual entry model in each subpop #####
# to get the annual number of incarceration in non-PWID in prison 
entry <- indicatorResult_uno(POC_AU, calibrateInit, "newEntry",
                             populations = POC_AU$popNames, endYear= endY) %>%
  mutate(year = year + POC_AU$cabY - 1)

##### annual number of incarceration/release/injection relapse/ stopping injection ##### 
PopTransition <- as.data.frame.table(calibrateInit$newpop_tran)%>%
  mutate(timestep = c(rep(seq(POC_AU$startYear, endY-POC_AU$timestep,
                              POC_AU$timestep),each = POC_AU$npops*POC_AU$npops)),
         from = Var1,
         To = Var2)%>%select(-c(Var1, Var2, Var3))
# giving the average number to first time step 
impute <- PopTransition%>%filter(timestep >1 &timestep <2)%>%
  group_by(from, To)%>%
  mutate(total_pop = sum(Freq)/length(Freq))%>%
  select(total_pop)%>%ungroup()

impute <- impute[c(1:as.numeric(POC_AU$npops*POC_AU$npops)),]

PopTransition[c(1:as.numeric(POC_AU$npops*POC_AU$npops)), "Freq"] <- impute$total_pop 

PopTransition_all <- cbind.data.frame(timestep = PopTransition$timestep, 
                                      from = PopTransition$from,  
                                      to = PopTransition$To, 
                                      best = PopTransition$Freq)%>%
  as_tibble()%>%
  mutate(year = rep(rep(seq(1, 100-1,1),each = 1/POC_AU$timestep),
                    each = POC_AU$npops*POC_AU$npops))

PPTranTo <- PopTransition_all%>%
  group_by(year, from, to)%>%summarise_at(.vars = "best", sum)

# incarceration 
incarce <- list()

incarce[["PWID"]] <- PPTranTo%>%filter(from == "C_PWID" & to == "P_PWID")

incarce[["fPWID"]] <- PPTranTo%>%filter(from == "C_fPWID" & to == "P_fPWID")

incarce_bind <- dplyr::bind_rows(incarce, .id = 'population')%>%
  mutate(year = POC_AU$cabY + year - 1)%>%ungroup()%>%
  select(year, population, best)

entry_nonPWID <- entry%>%filter(population =="P_nPWID") 

# total incarceration = incarceration in C_PWID + incarceration in C_fPWID + incarceration in non-PWID (entry:non-PWID)
incar_total <- rbind(incarce_bind, entry_nonPWID)%>%group_by(year)%>%
  summarize(best = sum(best))

# release 
release <- list()

release[["PWID"]] <- PPTranTo%>%filter(from == "P_PWID" & to == "C_PWID")

release[["fPWID"]] <- PPTranTo%>%filter(from == "P_fPWID" & to == "C_fPWID")

release_bind <- dplyr::bind_rows(release, .id = 'population')%>%
  mutate(year = POC_AU$cabY + year - 1)%>%ungroup()%>%select(year, population, best)

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



#####  incidence of stopping injecting ##### 
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
#####  incidence of injection relapse ##### 
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

##### incidence of incarceration  #####
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

##### incidence of release  #####
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
              ylabel = "Number",
              xlimits = c(POC_AU$startYear, 
                          POC_AU$startYear+30, 5),
              calibration_Y = POC_AU$cabY,
              rangeun = NULL, 
              groupPlot = NULL, 
              facetPlot = NULL,
              observationData = NULL, 
              simulateYear = (POC_AU$simY - POC_AU$cabY + 1 )) + theme_bw() +
  scale_y_continuous(limits = c(515000,525000)) + 
  ggtitle("Number of total population") + 
  geom_point(aes(
    y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Total", "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), 
             colour = "black") 
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
    colour = "black") + 
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
    colour = "black") + 
  geom_segment(
    aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Prison", "low"], 
        yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == "Prison", "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1))) + 
  scale_y_continuous(limits = c(0,50000), breaks = seq(0, 50000, 5000))

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
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black") + 
  geom_segment(
    aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[1], "low"], 
        yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[1], "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))

subpop_N_plot[[2]] <- subpop_N_plot[[2]] + 
  scale_y_continuous(limits = c(0, 700000), breaks = seq(0, 700000, 50000)) + 
  
  geom_point(aes(
    y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[2], "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black") + 
  geom_segment(
    aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[2], "low"], 
        yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[2], "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))

subpop_N_plot[[3]] <- subpop_N_plot[[3]] + 
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 500)) + 
  geom_point(aes(
    y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[3], "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black") + 
  geom_segment(
    aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[3], "low"], 
        yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[3], "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))
  

subpop_N_plot[[4]] <- subpop_N_plot[[4]] + 
  scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 20000, 1000)) + 
  geom_point(aes(
    y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[4], "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black") + 
  geom_segment(
    aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[4], "low"], 
        yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[4], "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))

subpop_N_plot[[5]] <- subpop_N_plot[[5]] + 
  scale_y_continuous(limits = c(0, 25000), breaks = seq(0, 25000, 1000)) + 
  geom_point(aes(
    y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[5], "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black") + 
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
  scale_y_continuous(limits = c(0,20), breaks = seq(0, 20,1)) + ggtitle("PWID in prison") + 
  geom_point(aes(
    y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_proP_fit", "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black") 
  
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
  scale_y_continuous(limits = c(0,40), breaks = seq(0, 40,5)) + ggtitle("Former PWID in prison") + 
  geom_point(aes(
    y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_profP_fit", "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black") 

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
  scale_y_continuous(limits = c(0,20), breaks = seq(0, 20,1)) + ggtitle("PWID in community") + 
  geom_point(aes(
    y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "commu_proP_fit", "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black") 

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
  scale_y_continuous(limits = c(0, 50000), breaks = seq(0, 50000, 1000)) 

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

inj_stop_inc_bind[[1]]
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
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 1))

endY_plot <- 2050- POC_AU$cabY

#### Results: incidence and prevalence #### 
##### HCV prevalence #####
# prevalence in each group: tempPrev_subpop

# number of susceptible 
tempNOTInfected_subpop <- popResults_MidYear(POC_AU, calibrateInit,Population = POC_AU$popNames,
                                          Disease_prog = NULL, 
                                          Cascade = c("s", "cured"), 
                                          param = NULL ,endYear = endY)%>%
  as_tibble()%>%group_by(year, population)%>%
  summarise(best = sum(best))%>%arrange(year, population)

# arrange order to align with other dts 
pop_N <- pop_N%>%arrange(year, population)


tempPrev_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                         population = POC_AU$popNames,
                         
                         as.data.frame(100*(pop_N[, -c(1,2)] - 
                                              tempNOTInfected_subpop[ ,-c(1,2)])/ 
                                         pop_N[ ,-c(1,2)]))%>%
  tibble::as_tibble()  


#prevalence in setting  
# tempPrev_setting[["setting]]

# community 
commu_N <- commu_N%>%arrange(year)

tempNOTInfected_commu <- popResults_MidYear(POC_AU, calibrateInit,
                                            Population = c("C_PWID", "C_fPWID"),
                                             Disease_prog = NULL, 
                                             Cascade = c("s", "cured"), 
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
                                            Cascade = c("s", "cured"), 
                                            param = NULL ,endYear = endY)%>%
  as_tibble()%>%group_by(year)%>%
  summarise(best = sum(best))%>%arrange(year)


tempPrev_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(100*(prison_N[, -c(1)] - 
                                                          tempNOTInfected_prison[ ,-c(1)])/ 
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
              observationData = NULL, 
              simulateYear = POC_AU$simY) +
  ggtitle("HCV prevalence by population") 

popPrevPlot <- popPrevPlot + 
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
                                                   c(0, 5)))
                   )) + theme_bw()


# setting 
tempPrev_setting_bind <- dplyr::bind_rows(tempPrev_setting, .id = 'population')%>%
  mutate(population = factor(population, 
                             levels = c("commu", "prisons"), 
                             labels = c("Community", "Prisons")))%>%
  as.data.frame()

str(tempPrev_setting_bind)


settingPrevPlot <- indicatorPlot(POC_AU, tempPrev_setting_bind, 
                             ylabel = "HCV prevalence (%)",
                             xlimits = c(POC_AU$startYear, 
                                         (POC_AU$startYear+endY_plot - 1), 5),
                             calibration_Y = POC_AU$cabY,
                             rangeun = NULL, 
                             groupPlot = NULL, 
                             facetPlot = population,
                             observationData = NULL, 
                             simulateYear = POC_AU$simY) +
  ggtitle("HCV prevalence by setting") 

settingPrevPlot <- settingPrevPlot + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 50))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 50))))) + theme_bw()
#####plot: incidence #####
HCVInc_subpop <- HCVInc_subpop%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname )) 

popIncPlot <- indicatorPlot(POC_AU, HCVInc_subpop, 
                             ylabel = "HCV incidence ",
                             xlimits = c(POC_AU$startYear, 
                                         (POC_AU$startYear+endY_plot), 5),
                             calibration_Y = POC_AU$cabY,
                             rangeun = NULL, 
                             groupPlot = NULL, 
                             facetPlot = population,
                             observationData = NULL, 
                             simulateYear = POC_AU$simY) +
  ggtitle("HCV incidence by population") 

popIncPlot <- popIncPlot + 
  facet_custom (~population,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 5))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 5))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 5))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 5))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 1)))
                  )) + theme_bw()

# setting 
HCVInc_setting_bind <- dplyr::bind_rows(HCVInc_setting, .id = 'population')%>%
  mutate(population = factor(population, 
                             levels = c("commu", "prisons"), 
                             labels = c("Community", "Prisons")))%>%
  as.data.frame()

settingIncPlot <- indicatorPlot(POC_AU, HCVInc_setting_bind, 
                                 ylabel = "HCV incidence",
                                 xlimits = c(POC_AU$startYear, 
                                             (POC_AU$startYear+endY_plot - 1), 5),
                                 calibration_Y = POC_AU$cabY,
                                 rangeun = NULL, 
                                 groupPlot = NULL, 
                                 facetPlot = population,
                                 observationData = NULL, 
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
                                                     c(0, 20))))) + theme_bw()



#### Results: disease progression ####

#### Results: HCV test and treatment ####




#### debugging 


x_entry <- indicatorResults(POC_AU, calibrateInit, "newEntry", 
                 pop=POC_AU$popNames,
                 paramR = NULL, range = NULL,
                 endY = endY)

x_leave <- indicatorResults(POC_AU, calibrateInit, "newLeave", 
                            pop=POC_AU$popNames,
                            paramR = NULL, range = NULL,
                            endY = endY)

x_death <- indicatorResults(POC_AU, calibrateInit, "newDeath", 
                            pop=POC_AU$popNames,
                            paramR = NULL, range = NULL,
                            endY = endY)

x_death_hcv <- indicatorResults(POC_AU, calibrateInit, "newHCVdeaths", 
                           pop=POC_AU$popNames,
                           paramR = NULL, range = NULL,
                           endY = endY)

xt <- cbind(x_entry, x_leave[ ,3], x_death[ ,3], x_death_hcv[, 3])

colnames(xt) <- c("year", "population", "entry", "leave", "death", "death_hcv")

xt <- xt%>%mutate(toleave = leave + death + death_hcv)

ggplot(xt, aes(x =year, y = entry, group = population)) + 
  geom_line(aes(color = population))

xt_pnPWID <- xt%>%filter(population =="P_nPWID")%>%mutate(diff = entry - toleave)


xt_CPWID <- xt%>%filter(population != "P_nPWID")%>%group_by(year)%>%
  summarize(entry = sum(entry), 
            leave = sum(leave),
            death = sum(death),
            death_hcv = sum(death_hcv), 
            toleave = sum(toleave))
xt_CPWID <- xt_CPWID%>%mutate(diff = entry - toleave)


xtt <- indicatorResults(POC_AU, calibrateInit, "newS", 
                        pop=POC_AU$popNames,
                        paramR = NULL, range = NULL,
                        endY = endY)
ggplot(xtt, aes(x = year , y = best, group = population)) +
  geom_line(aes(color = population))

xtt_l <-
indicatorResults(POC_AU, calibrateInit, "newLeave", 
                 pop=POC_AU$popNames,
                 paramR = NULL, range = NULL,
                 endY = endY)
ggplot(xtt_l, aes(x = year , y = best, group = population)) +
  geom_line(aes(color = population))
