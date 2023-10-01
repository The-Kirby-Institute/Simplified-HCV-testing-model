# copy from POC_AU 
# indicators need to revise

rm(list = ls())

library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library("readxl")

project_name <- "TWPrisoners"

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

OutputFolder <- file.path(data_path, "02. Output/RDA")

load(file.path(OutputFolder, paste0(project_name, ".rda")))

load(file.path(OutputFolder, paste0(project_name, "cali_init", ".rda")))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 

# setting end Year 

endY <- 100
####import calibration data points #### 
CaliFolder <- DataFolder%>%dirname()
files <- list.files(path = paste0(CaliFolder , 
                                  "/Calibration_dt/", sep = ""),
                    pattern = '*.xlsx')

calib_dt <- lapply(files, function(f) {
  df <- read_xlsx(file.path(paste0(CaliFolder , 
                                  "/Calibration_dt/", f, sep = "")))
  
  df <- df%>%as_tibble()
  
})
names(calib_dt) <- c("N", "frac")

calib_dt[["N"]] <- calib_dt[["N"]]%>%
  mutate(indicators = factor(indicators, levels = c(unique(.$indicators)),
                              label = c("prison", 
                                        "rehab", 
                                        "incar_P", 
                                        "incar_DI", 
                                        "incar_D", 
                                        "rel_P", 
                                        "rel_D")))

calib_dt[[2]] <- calib_dt[[2]]%>%mutate(time = time - TWPrisoners$cabY +1 )
calib_dt[[1]] <- calib_dt[[1]]%>%mutate(time = time - TWPrisoners$cabY +1 )
endY <- 15
#### Result:population ####
##### N: each subpop/ total  #####
# subpop as list
subpop_N <- lapply(TWPrisoners$popNames, function(x){ 
  
  a <- popResults_MidYear(TWPrisoners, calibrateInit,
                          Population = x,
                          Disease_prog = NULL, 
                          Cascade = NULL, param = NULL, 
                          endYear = endY, YearCut = "END")%>%ungroup()
})

names(subpop_N) <- TWPrisoners$popNames
subpop_N$N_inca_NCID
ggplot(data = as.data.frame(subpop_N[[6]]), aes(x = year, y = best)) + 
  geom_line() 



# all subpop in one list 
pop_N <- dplyr::bind_rows(subpop_N, .id = 'population')

# total number 

total_N <- pop_N%>%group_by(year)%>%summarise(best = sum(best))

##### N: community #####
commu_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                              Population = c("N_inca_NCID", 
                                             "N_inca_CID",
                                             "E_inca_NCID",
                                             "E_inca_CID"),
                              Disease_prog = NULL, 
                              Cascade = NULL, param = NULL, 
                              endYear = endY, YearCut = "END")%>%ungroup()%>%
  dplyr::group_by(year)%>%summarise_at("best",sum)

##### N: prison #####
prison_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                               Population = c("D_inca", "P_inca"),
                               Disease_prog = NULL, 
                               Cascade = NULL, param = NULL, 
                               endYear = endY, YearCut = "END")%>%ungroup()%>%
  dplyr::group_by(year)%>%summarise_at("best",sum)

##### %: PWID in community/prison #####
commu_proP <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                    as.data.frame((pop_N[pop_N$population =="N_inca_CID",3] + 
                                    pop_N[pop_N$population =="E_inca_CID",3])/ commu_N[ ,-1])*100)%>%
  tibble::as_tibble()

prison_proD <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                     as.data.frame(pop_N[pop_N$population =="D_inca",3] / prison_N[ ,-1])*100)%>%
  tibble::as_tibble() 

commu_proPNinca <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                    as.data.frame((pop_N[pop_N$population =="N_inca_CID",3])/ (pop_N[pop_N$population =="N_inca_CID",3] + 
                                                                                   pop_N[pop_N$population =="E_inca_CID",3]))*100)%>%
  tibble::as_tibble()

##### % prison in total pops ##### 
prison_pro <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                      as.data.frame(prison_N[ ,-1] / total_N[ ,-1])*100)%>%
  tibble::as_tibble()


####========================================================================####

##### number of annual leaving in each subpop #####
# this is for getting annual number of released in non-PWID in prison 
rel <- indicatorResult_uno(TWPrisoners, calibrateInit, "newLeave",
                           populations = TWPrisoners$popNames, endYear= endY)%>%
  mutate(year = year + TWPrisoners$cabY - 1) 

reD <- indicatorResult_uno(TWPrisoners, calibrateInit, "newDeath",
                           populations = TWPrisoners$popNames, endYear= endY)%>%
  mutate(year = year + TWPrisoners$cabY - 1) 

reHCVD <- indicatorResult_uno(TWPrisoners, calibrateInit, "newHCVdeaths",
                           populations = TWPrisoners$popNames, endYear= endY)%>%
  mutate(year = year + TWPrisoners$cabY - 1) 

##### number of annual entry model in each subpop #####
# to get the annual number of incarceration in non-PWID in prison 
entry <- indicatorResult_uno(TWPrisoners, calibrateInit, "newEntry",
                             populations = TWPrisoners$popNames, endYear= endY) %>%
  mutate(year = year + TWPrisoners$cabY - 1)

##### annual number of incarceration/release/injection relapse/ stopping injection ##### 
PopTransition <- as.data.frame.table(calibrateInit$newpop_tran)%>%
  mutate(timestep = c(rep(seq(TWPrisoners$startYear, endY-TWPrisoners$timestep,
                              TWPrisoners$timestep),each = TWPrisoners$npops*TWPrisoners$npops)),
         from = Var1,
         To = Var2)%>%dplyr::select(-c(Var1, Var2, Var3))
# giving the average number to first time step 
impute <- PopTransition%>%filter(timestep >1 &timestep <2)%>%
  group_by(from, To)%>%
  mutate(total_pop = sum(Freq)/length(Freq))%>%
  dplyr::select(total_pop)%>%ungroup()

impute <- impute[c(1:as.numeric(TWPrisoners$npops*TWPrisoners$npops)),]

PopTransition[c(1:as.numeric(TWPrisoners$npops*TWPrisoners$npops)), "Freq"] <- impute$total_pop 

PopTransition_all <- cbind.data.frame(timestep = PopTransition$timestep, 
                                      from = PopTransition$from,  
                                      to = PopTransition$To, 
                                      best = PopTransition$Freq)%>%
  as_tibble()%>%
  mutate(year = rep(rep(seq(1, endY-1,1),each = 1/TWPrisoners$timestep),
                    each = TWPrisoners$npops*TWPrisoners$npops))

PPTranTo <- PopTransition_all%>%
  group_by(year, from, to)%>%summarise_at(.vars = "best", sum)

# incarceration 
incarce_D <- list()

incarce_D[["nonPWID"]] <- PPTranTo%>%
  filter(from %in% c("N_inca_NCID", "E_inca_NCID") & to == "D_inca")

incarce_D[["PWID"]] <- PPTranTo%>%
  filter(from %in% c("N_inca_CID", "E_inca_CID") & to == "D_inca")%>%
  group_by(year)%>%summarise(best = sum(best))

incarce_D_bind <- dplyr::bind_rows(incarce_D, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()%>%
  dplyr::select(year, population, best)%>%group_by(year)%>%
  summarise(best = sum(best))

incarce_P <- list()

incarce_P[["nonPWID"]] <- PPTranTo%>%
  filter(from %in% c("N_inca_NCID", "E_inca_NCID") & to == "P_inca")

incarce_P[["PWID"]] <- PPTranTo%>%
  filter(from %in% c("N_inca_CID", "E_inca_CID") & to == "P_inca")



incarce_P_bind <- dplyr::bind_rows(incarce_P, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()%>%
  dplyr::select(year, population, best)%>%group_by(year)%>%
  summarise(best = sum(best))


incarce_exp <- list()



incarce_exp[["N"]] <- PPTranTo%>%
  filter(from %in% c("N_inca_NCID", "N_inca_CID") & to %in% c("D_inca", "P_inca"))

incarce_exp[["E"]] <- PPTranTo%>%
  filter(from %in% c("E_inca_NCID", "E_inca_CID") & to %in% c("D_inca", "P_inca"))

incarce_exp_bind <- dplyr::bind_rows(incarce_exp, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()%>%
  dplyr::select(year, population, best)

incarce <- PPTranTo%>%
  filter(to %in% c("D_inca", "P_inca"))%>%
  filter(!from %in% c("D_inca", "P_inca"))

# release 
release <- list()

release[["D"]] <- PPTranTo%>%
  filter(from == "D_inca" & to %in% c("E_inca_NCID"))

release[["P"]] <- PPTranTo%>%
  filter(from == "P_inca" & to %in% c("E_inca_NCID"))

release_bind <- dplyr::bind_rows(release, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()%>%
  dplyr::select(year, population, best)

release_bind <- release_bind%>%group_by(year, population)%>%summarize(best = sum(best))

release_total <- rbind(release_bind)%>%group_by(year)%>%
  summarize(best = sum(best))

# injection 
inj_relap <- list()

inj_relap[["N"]] <- PPTranTo%>%filter(from == "N_inca_NCID" & to == "N_inca_CID")

inj_relap[["E"]] <- PPTranTo%>%filter(from == "E_inca_NCID" & to == "E_inca_CID")

inj_relap_bind <- dplyr::bind_rows(inj_relap, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()


# stopping injection 
inj_stop <- list()

inj_stop[["N"]] <- PPTranTo%>%filter(from == "N_inca_CID" & to == "N_inca_NCID")

inj_stop[["E"]] <- PPTranTo%>%filter(from == "E_inca_CID" & to == "E_inca_NCID")

inj_stop_bind <- dplyr::bind_rows(inj_stop, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()



#####  incidence of stopping injecting ##### 
inj_stop_inc <- list()

inj_stop_inc[["N"]] <- 
  cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
        as.data.frame(inj_stop[["N"]][ , "best"] / 
                        subpop_N[["N_inca_CID"]][ ,"best"])*100)%>%
  tibble::as_tibble()

inj_stop_inc[["E"]] <- 
  cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
        as.data.frame(inj_stop[["E"]][ , "best"] / 
                        subpop_N[["E_inca_CID"]][ ,"best"])*100)%>%
  tibble::as_tibble()

inj_stop_inc_bind <- dplyr::bind_rows(inj_stop_inc, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()
#####  incidence of injection relapse ##### 
inj_relap_inc <- list()

inj_relap_inc[["N"]] <- 
  cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
        as.data.frame(inj_relap[["N"]][ , "best"] / 
                        subpop_N[["N_inca_NCID"]][ ,"best"])*100)%>%
  tibble::as_tibble()

inj_relap_inc[["E"]] <- 
  cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
        as.data.frame(inj_relap[["E"]][ , "best"] / 
                        subpop_N[["E_inca_NCID"]][ ,"best"])*100)%>%
  tibble::as_tibble()

inj_relap_inc_bind <- dplyr::bind_rows(inj_relap_inc, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()

##### incidence of incarceration  #####
# combine subpop number into one tibble
incar_dominator <- subpop_N[c(unique(incarce$from))]

incar_dom <- rbind(incar_dominator$N_inca_NCID, incar_dominator$N_inca_NCID, 
              incar_dominator$N_inca_CID, incar_dominator$N_inca_CID, 
              incar_dominator$E_inca_NCID, incar_dominator$E_inca_NCID, 
              incar_dominator$E_inca_CID, incar_dominator$E_inca_CID)%>%
  arrange(year,population)%>%mutate(dominator = best)



incarce <- incarce%>%arrange(year, from)%>%ungroup()%>%
  mutate(dominator = incar_dom$dominator,
         incar_inc = best/dominator*100)

rel
##### incidence of release  #####
release_inc <- list()

release_inc[["D"]] <- 
  cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
        as.data.frame(release[["D"]][ , "best"] / 
                        subpop_N[["D_inca"]][ ,"best"])*100)%>%
  tibble::as_tibble()

release_inc[["P"]] <- 
  cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
        as.data.frame(release[["P"]][ , "best"] /
                        subpop_N[["P_inca"]][ ,"best"])*100)%>%
  tibble::as_tibble()

release_inc_bind <- dplyr::bind_rows(release_inc, .id = 'population')%>%
  mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()%>%
  dplyr::select(year, population, best)

#==============================================================================#

#==============================================================================#


####Plot: population ####
#####plot: number of total pop  #####

totalPop_plot <- indicatorPlot(TWPrisoners, total_N, 
                               ylabel = "Number",
                               xlimits = c(TWPrisoners$startYear, 
                                           TWPrisoners$startYear+30, 5),
                               calibration_Y = TWPrisoners$cabY,
                               rangeun = NULL, 
                               groupPlot = NULL, 
                               facetPlot = NULL,
                               observationData = calib_dt[[2]]%>%filter(indicator == "Total_pop"), 
                               simulateYear = (TWPrisoners$simY - TWPrisoners$cabY + 1 )) + theme_bw() +
  scale_y_continuous(limits = c(19500000,21000000)) + 
  ggtitle("Number of total population") 

totalPop_plot


#####plot: number of pop in community/prison   #####
commuPop_plot <- indicatorPlot(TWPrisoners, commu_N, 
                               ylabel = "Number",
                               xlimits = c(TWPrisoners$startYear, 
                                           TWPrisoners$startYear+30, 5),
                               calibration_Y = TWPrisoners$cabY,
                               rangeun = NULL, 
                               groupPlot = NULL, 
                               facetPlot = NULL,
                               observationData = NULL, 
                               simulateYear = NULL) + theme_bw() +
  ggtitle("Number in community") + 
  scale_y_continuous(limits = c(19000000, 21000000))

prisonPop_plot <- indicatorPlot(TWPrisoners, prison_pro, 
                                ylabel = "Number",
                                xlimits = c(TWPrisoners$startYear, 
                                            TWPrisoners$startYear+30, 5),
                                calibration_Y = TWPrisoners$cabY,
                                rangeun = NULL, 
                                groupPlot = NULL, 
                                facetPlot = NULL,
                                observationData = calib_dt[[2]]%>%
                                  filter(indicator == "prison_pro_fit"), 
                                simulateYear = NULL) + theme_bw() +
  ggtitle("Number of people in correctional settings") 
  


#####plot: number of subpop   ####
# sub pop full names 
Namelab <- c("Never experienced incarceration non current PWID in community",
             "Never experienced incarceration PWID in community", 
             "Ever incarcerated non current PWID in community",
             "Ever incarcerated PWID in community", 
             "People in detention", 
             "People in prisons")

subpop_N_plot <- list() 

for(x in seq_along(names(subpop_N))){ 
  
  subpop_N_plot[[x]] <- indicatorPlot(TWPrisoners, subpop_N[[x]] , 
                                      ylabel = "Number",
                                      xlimits = c(TWPrisoners$startYear, 
                                                  TWPrisoners$startYear+30, 5),
                                      calibration_Y = TWPrisoners$cabY,
                                      rangeun = NULL, 
                                      groupPlot = NULL, 
                                      facetPlot = NULL,
                                      observationData = NULL, 
                                      simulateYear = NULL) +
    theme_bw() + ggtitle(paste0(Namelab[x])) 
  
}

subpop_N_plot[[1]] <- subpop_N_plot[[1]] + 
  scale_y_continuous(limits = c(0, 21000000))
  
subpop_N_plot[[2]] <- subpop_N_plot[[2]] + 
  scale_y_continuous(limits = c(0, 3000)) 

subpop_N_plot[[3]] <- subpop_N_plot[[3]] + 
  scale_y_continuous(limits = c(0, 300000)) 

subpop_N_plot[[4]] <- subpop_N_plot[[4]] + 
  scale_y_continuous(limits = c(0, 50000)) 

subpop_N_plot[[5]] <- subpop_N_plot[[5]] + 
  scale_y_continuous(limits = c(0, 100000)) + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "rehab"),aes(
    y = realPop,
    x = time), colour = "black")
  
subpop_N_plot[[6]] <- subpop_N_plot[[6]] + 
  scale_y_continuous(limits = c(0,70000)) + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "prison"),aes(
    y = realPop,
    x = time), colour = "black")

#####plot: % ####
# 
frac_prison_plot <- indicatorPlot(TWPrisoners, prison_pro , 
                                 ylabel = "Percentage (%)",
                                 xlimits = c(TWPrisoners$startYear, 
                                             TWPrisoners$startYear+30, 5),
                                 calibration_Y = TWPrisoners$cabY,
                                 rangeun = NULL, 
                                 groupPlot = NULL, 
                                 facetPlot = NULL,
                                 observationData = calib_dt[[2]]%>%filter(indicator == "prison_pro_fit"), 
                                 simulateYear = NULL) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1,0.1)) + ggtitle("% of people in correctional settings") 

# % of people in detention among people in correctional settings 
frac_D_plot <- indicatorPlot(TWPrisoners, prison_proD , 
                                  ylabel = "Percentage (%)",
                                  xlimits = c(TWPrisoners$startYear, 
                                              TWPrisoners$startYear+30, 5),
                                  calibration_Y = TWPrisoners$cabY,
                                  rangeun = NULL, 
                                  groupPlot = NULL, 
                                  facetPlot = NULL,
                                  observationData = calib_dt[[2]]%>%filter(indicator == "prison_proD_fit"), 
                                  simulateYear = NULL) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100,10)) + 
  ggtitle("% of people in detention among people in correctional settings")


# % of PWID in the community 

frac_commuP_plot <- indicatorPlot(TWPrisoners, commu_proP , 
                                 ylabel = "Percentage (%)",
                                 xlimits = c(TWPrisoners$startYear, 
                                             TWPrisoners$startYear+30, 5),
                                 calibration_Y = TWPrisoners$cabY,
                                 rangeun = NULL, 
                                 groupPlot = NULL, 
                                 facetPlot = NULL,
                                 observationData = calib_dt[[2]]%>%filter(indicator == "commu_proP_fit"), 
                                 simulateYear = NULL) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0, 0.5,0.05)) + 
  ggtitle("PWID in community") 


#####plot: annual number entry to/release from detention #### 
incarce_D[["PWID"]] <- incarce_D[["PWID"]]%>%
  mutate(year = year + TWPrisoners$cabY - 1)

New_incarDI_plot <- indicatorPlot(TWPrisoners, incarce_D[["PWID"]],
                                  ylabel = "Number",
                                  xlimits = c(TWPrisoners$cabY,
                                              TWPrisoners$cabY+30, 5),
                                  calibration_Y = TWPrisoners$cabY,
                                  rangeu = NULL,
                                  groupPlot = NULL,
                                  facetPlot = NULL,
                                  observationData = NULL, 
                                  simulateYear = NULL) +
  ggtitle("Annual number entry detention") + theme_bw() + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "incar_DI"),aes(
    y = realPop,
    x = time + TWPrisoners$cabY - 1), colour = "black") + 
  scale_y_continuous(limits = c(0, 8000))


New_incarD_plot <-  indicatorPlot(TWPrisoners, incarce_D_bind,
                                 ylabel = "Number",
                                 xlimits = c(TWPrisoners$cabY,
                                             TWPrisoners$cabY+30, 5),
                                 calibration_Y = TWPrisoners$cabY,
                                 rangeu = NULL,
                                 groupPlot = NULL,
                                 facetPlot = NULL,
                                 observationData = NULL, 
                                 simulateYear = NULL) +
  ggtitle("Annual number entry detention") + theme_bw() + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "incar_D"),aes(
    y = realPop,
    x = time + TWPrisoners$cabY - 1), colour = "black")

New_incarP_plot <-  indicatorPlot(TWPrisoners, incarce_P_bind,
                                  ylabel = "Number",
                                  xlimits = c(TWPrisoners$cabY,
                                              TWPrisoners$cabY+30, 5),
                                  calibration_Y = TWPrisoners$cabY,
                                  rangeu = NULL,
                                  groupPlot = NULL,
                                  facetPlot = NULL,
                                  observationData = NULL, 
                                  simulateYear = NULL) +
  ggtitle("Annual number entry prisons") + theme_bw() + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "incar_P"),aes(
    y = realPop,
    x = time + TWPrisoners$cabY - 1), colour = "black")

New_releaseD_plot <-  indicatorPlot(TWPrisoners, release_bind%>%
                                      filter(population == "D"),
                                   ylabel = "Number",
                                   xlimits = c(TWPrisoners$cabY,
                                               TWPrisoners$cabY+30, 5),
                                   calibration_Y = TWPrisoners$cabY,
                                   rangeu = NULL,
                                   groupPlot = NULL,
                                   facetPlot = NULL,
                                   observationData = NULL, 
                                   simulateYear = NULL) +
  ggtitle("Annual number release from detention") + theme_bw() + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "rel_D"),aes(
    y = realPop,
    x = time + TWPrisoners$cabY - 1), colour = "black")

New_releaseP_plot <-  indicatorPlot(TWPrisoners, release_bind%>%
                                      filter(population == "P"),
                                    ylabel = "Number",
                                    xlimits = c(TWPrisoners$cabY,
                                                TWPrisoners$cabY+30, 5),
                                    calibration_Y = TWPrisoners$cabY,
                                    rangeu = NULL,
                                    groupPlot = NULL,
                                    facetPlot = NULL,
                                    observationData = NULL, 
                                    simulateYear = NULL) +
  ggtitle("Annual number release from detention") + theme_bw() + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "rel_P"),aes(
    y = realPop,
    x = time + TWPrisoners$cabY - 1), colour = "black")

#####plot: annual number stopping/relapse injection ####
# stopping injection
New_stopinj_plot <-  indicatorPlot(TWPrisoners, inj_stop_bind,
                                   ylabel = "Number",
                                   xlimits = c(TWPrisoners$cabY,
                                               TWPrisoners$cabY+30, 5),
                                   calibration_Y = TWPrisoners$cabY,
                                   rangeu = NULL,
                                   groupPlot = NULL,
                                   facetPlot = population,
                                   observationData = NULL, 
                                   simulateYear = NULL) +
  ggtitle("Annual number stop injection") + theme_bw() + 
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 500)) 
# relapse injection
New_relapinj_plot <-  indicatorPlot(TWPrisoners, inj_relap_bind,
                                    ylabel = "Number",
                                    xlimits = c(TWPrisoners$cabY,
                                                TWPrisoners$cabY+30, 5),
                                    calibration_Y = TWPrisoners$cabY,
                                    rangeu = NULL,
                                    groupPlot = NULL,
                                    facetPlot = population,
                                    observationData = NULL, 
                                    simulateYear = NULL) +
  ggtitle("Annual number relapse injection") + theme_bw() + 
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, 5000)) 

##============================================================================##


#####plot: incidence of entry to/release from prisons #### 
Inc_incar_plot <- indicatorPlot(TWPrisoners,incar_inc_bind,
                                ylabel = "Incidence",
                                xlimits = c(TWPrisoners$cabY,
                                            TWPrisoners$cabY+30, 5),
                                calibration_Y = TWPrisoners$cabY,
                                rangeu = NULL,
                                groupPlot = NULL,
                                facetPlot = population,
                                observationData = NULL, 
                                simulateYear = NULL) +
  ggtitle("Incidence of incarceration") + theme_bw() + 
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 1)) 

Inc_release_plot <-  indicatorPlot(TWPrisoners, release_inc_bind,
                                   ylabel = "Incidence",
                                   xlimits = c(TWPrisoners$cabY,
                                               TWPrisoners$cabY+30, 5),
                                   calibration_Y = TWPrisoners$cabY,
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
  
  a <- indicatorPlot(TWPrisoners, x,
                     ylabel = "Incidence",
                     xlimits = c(TWPrisoners$startYear ,
                                 TWPrisoners$startYear+30, 5),
                     calibration_Y = TWPrisoners$cabY,
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
    x = TWPrisoners$simY - TWPrisoners$cabY + 1 ), 
    colour = "black") + 
  geom_segment(
    aes(y = datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "commu_stopinj", "low"], 
        yend = datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "commu_stopinj", "up"], 
        x = TWPrisoners$simY - TWPrisoners$cabY + 1, xend = TWPrisoners$simY - TWPrisoners$cabY + 1 )) + 
  scale_y_continuous(limits = c(0,20), breaks = seq(0, 20, 1))



# relapse injection
Inc_relapinj_plot <-  indicatorPlot(TWPrisoners, inj_relap_inc_bind,
                                    ylabel = "Incidence",
                                    xlimits = c(TWPrisoners$cabY,
                                                TWPrisoners$cabY+30, 5),
                                    calibration_Y = TWPrisoners$cabY,
                                    rangeu = NULL,
                                    groupPlot = NULL,
                                    facetPlot = population,
                                    observationData = NULL, 
                                    simulateYear = NULL) +
  ggtitle("Incidence of relapse injection") + theme_bw() + 
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 1))

endY_plot <- 2050- TWPrisoners$cabY

#### Results: incidence and prevalence #### 
##### HCV prevalence #####
# prevalence in each group: tempPrev_subpop

# number of susceptible 
tempNOTInfected_subpop <- popResults_MidYear(TWPrisoners, calibrateInit,Population = TWPrisoners$popNames,
                                             Disease_prog = NULL, 
                                             Cascade = c("s", "cured"), 
                                             param = NULL ,endYear = endY)%>%
  as_tibble()%>%group_by(year, population)%>%
  summarise(best = sum(best))%>%arrange(year, population)

# arrange order to align with other dts 
pop_N <- pop_N%>%arrange(year, population)


tempPrev_subpop <- cbind(year = rep(seq(TWPrisoners$startYear , endY-1 ,1), each = TWPrisoners$npops),
                         population = TWPrisoners$popNames,
                         
                         as.data.frame(100*(pop_N[, -c(1,2)] - 
                                              tempNOTInfected_subpop[ ,-c(1,2)])/ 
                                         pop_N[ ,-c(1,2)]))%>%
  tibble::as_tibble()  


#prevalence in setting  
# tempPrev_setting[["setting]]

# community 
commu_N <- commu_N%>%arrange(year)

tempNOTInfected_commu <- popResults_MidYear(TWPrisoners, calibrateInit,
                                            Population = c("C_PWID", "C_fPWID"),
                                            Disease_prog = NULL, 
                                            Cascade = c("s", "cured"), 
                                            param = NULL ,endYear = endY)%>%
  as_tibble()%>%group_by(year)%>%
  summarise(best = sum(best))%>%arrange(year)

tempPrev_setting <- list()

tempPrev_setting[["commu"]] <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                                     as.data.frame(100*(commu_N[, -c(1)] - 
                                                          tempNOTInfected_commu[ ,-c(1)])/ 
                                                     commu_N[ ,-c(1)]))%>%tibble::as_tibble()

# prison 
## PWID + former PWID + nonPWID 
prison_N <- prison_N%>%arrange(year)

tempNOTInfected_prison <- popResults_MidYear(TWPrisoners, calibrateInit,
                                             Population = c("P_PWID", "P_fPWID", 
                                                            "P_nPWID"),
                                             Disease_prog = NULL, 
                                             Cascade = c("s", "cured"), 
                                             param = NULL ,endYear = endY)%>%
  as_tibble()%>%group_by(year)%>%
  summarise(best = sum(best))%>%arrange(year)


tempPrev_setting[["prisons"]] <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                                       as.data.frame(100*(prison_N[, -c(1)] - 
                                                            tempNOTInfected_prison[ ,-c(1)])/ 
                                                       prison_N[ ,-c(1)]))%>%tibble::as_tibble()



##### HCV incidence ##### 
HCVInfect_subpop <- indicatorResults(TWPrisoners, calibrateInit, "newInfections", 
                                     pop=TWPrisoners$popNames,
                                     paramR = NULL, range = NULL,
                                     endY = endY)

HCVInc_subpop <- cbind(year = rep(seq(TWPrisoners$startYear , endY-1 ,1), 
                                  each = TWPrisoners$npops),
                       population = TWPrisoners$popNames,
                       as.data.frame(HCVInfect_subpop[, -c(1,2)] / 
                                       pop_N[ ,-c(1,2)]*100))%>%
  tibble::as_tibble() 



HCVInfect_setting <- list()

HCVInfect_setting[["commu"]] <- indicatorResults(TWPrisoners, 
                                                 calibrateInit, 
                                                 "newInfections", 
                                                 pop= c("C_PWID", "C_fPWID"),
                                                 paramR = NULL, range = NULL,
                                                 endY = endY)%>%group_by(year)%>%
  summarise(best = sum(best))%>%arrange(year)

HCVInfect_setting[["prisons"]] <- indicatorResults(TWPrisoners, 
                                                   calibrateInit, 
                                                   "newInfections", 
                                                   pop= c("P_PWID", "P_fPWID",
                                                          "P_nPWID"), 
                                                   paramR = NULL, 
                                                   range = NULL,
                                                   endY = endY)%>%group_by(year)%>%
  summarise(best = sum(best))%>%arrange(year)

HCVInc_setting <- list()

HCVInc_setting[["commu"]] <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                                   as.data.frame(100*(HCVInfect_setting[["commu"]][ , -c(1)]/ 
                                                        commu_N[ ,-c(1)])))%>%
  tibble::as_tibble()
HCVInc_setting[["prisons"]] <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
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
                             levels = TWPrisoners$popNames, 
                             labels = pop_labname ))

popPrevPlot <- indicatorPlot(TWPrisoners, tempPrev_subpop, 
                             ylabel = "HCV prevalence (%)",
                             xlimits = c(TWPrisoners$startYear, 
                                         (TWPrisoners$startYear+endY_plot), 5),
                             calibration_Y = TWPrisoners$cabY,
                             rangeun = NULL, 
                             groupPlot = NULL, 
                             facetPlot = population,
                             observationData = NULL, 
                             simulateYear = TWPrisoners$simY) +
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


settingPrevPlot <- indicatorPlot(TWPrisoners, tempPrev_setting_bind, 
                                 ylabel = "HCV prevalence (%)",
                                 xlimits = c(TWPrisoners$startYear, 
                                             (TWPrisoners$startYear+endY_plot - 1), 5),
                                 calibration_Y = TWPrisoners$cabY,
                                 rangeun = NULL, 
                                 groupPlot = NULL, 
                                 facetPlot = population,
                                 observationData = NULL, 
                                 simulateYear = TWPrisoners$simY) +
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
                             levels = TWPrisoners$popNames, 
                             labels = pop_labname )) 

popIncPlot <- indicatorPlot(TWPrisoners, HCVInc_subpop, 
                            ylabel = "HCV incidence ",
                            xlimits = c(TWPrisoners$startYear, 
                                        (TWPrisoners$startYear+endY_plot), 5),
                            calibration_Y = TWPrisoners$cabY,
                            rangeun = NULL, 
                            groupPlot = NULL, 
                            facetPlot = population,
                            observationData = NULL, 
                            simulateYear = TWPrisoners$simY) +
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

settingIncPlot <- indicatorPlot(TWPrisoners, HCVInc_setting_bind, 
                                ylabel = "HCV incidence",
                                xlimits = c(TWPrisoners$startYear, 
                                            (TWPrisoners$startYear+endY_plot - 1), 5),
                                calibration_Y = TWPrisoners$cabY,
                                rangeun = NULL, 
                                groupPlot = NULL, 
                                facetPlot = population,
                                observationData = NULL, 
                                simulateYear = TWPrisoners$simY) +
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


x_entry <- indicatorResults(TWPrisoners, calibrateInit, "newEntry", 
                            pop=TWPrisoners$popNames,
                            paramR = NULL, range = NULL,
                            endY = endY)

x_leave <- indicatorResults(TWPrisoners, calibrateInit, "newLeave", 
                            pop=TWPrisoners$popNames,
                            paramR = NULL, range = NULL,
                            endY = endY)

x_death <- indicatorResults(TWPrisoners, calibrateInit, "newDeath", 
                            pop=TWPrisoners$popNames,
                            paramR = NULL, range = NULL,
                            endY = endY)

x_death_hcv <- indicatorResults(TWPrisoners, calibrateInit, "newHCVdeaths", 
                                pop=TWPrisoners$popNames,
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


xtt <- indicatorResults(TWPrisoners, calibrateInit, "newS", 
                        pop=TWPrisoners$popNames,
                        paramR = NULL, range = NULL,
                        endY = endY)
ggplot(xtt, aes(x = year , y = best, group = population)) +
  geom_line(aes(color = population))

xtt_l <-
  indicatorResults(TWPrisoners, calibrateInit, "newLeave", 
                   pop=TWPrisoners$popNames,
                   paramR = NULL, range = NULL,
                   endY = endY)
ggplot(xtt_l, aes(x = year , y = best, group = population)) +
  geom_line(aes(color = population))
