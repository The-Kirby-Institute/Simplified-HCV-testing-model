# mannually fit 
rm(list = ls())
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library("readxl")
library(ggpubr)

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



source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/check_steady.R"))

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

# get median number over 2019-2022 as data point for 2022 
calib_dt[["N"]] <- calib_dt[["N"]]%>%group_by(indicators)%>%
  mutate(low = quantile(realPop, 0.25),
         up = quantile(realPop, 0.75),
         realPop = median(realPop))%>%slice(n())%>%
  mutate(YearCut = ifelse(is.na(YearCut), "Mid", YearCut))

endY <- 15

#### Steps #### 
#### 1. adjust pop_array ####
rownames(pop_array) <- TWPrisoners$popNames
colnames(pop_array) <- TWPrisoners$popNames
pop_array[1,2 ,] <- 0.00002
pop_array[1,5 ,] <- 3.043484e-04
pop_array[1,6,] <- 2.5e-04

pop_array[2,1 ,] <- 0.058
pop_array[2,5 ,] <- 0.04
pop_array[2,6 ,] <- 0.3

pop_array[3,4 ,] <- 0.09
pop_array[3,5 ,] <- 7.543484e-05
pop_array[3,6,  ] <- 0.015

pop_array[4,3 ,] <- 0.01
pop_array[4,5 ,] <- 0.1
pop_array[4,6 ,] <- 0.7

pop_array[5,3,] <- 0.995
pop_array[6,3,] <- 0.48
#### 2. run steady and extract initial pops ####
tic <- proc.time()

steady <- HCVMSM(TWPrisoners, best_estimates, best_initial_pop,
                 disease_progress, pop_array, dfList, fib, end_Y = NULL, 
                 modelrun = "steady", proj = "TWPrisoners")

toc <- proc.time() - tic 
toc



check_steady(model_result = steady, endY = TWPrisoners$endYear,
             timestep = TWPrisoners$timestep, 
             Ncomp = TWPrisoners$ncomponent*TWPrisoners$npops, 
             Tequilibrium = 1500)


df_list <- lapply(steady, as.data.frame.table)

popPro_extract <- df_list$allPops%>%
  mutate(time = rep(seq(TWPrisoners$startYear, (TWPrisoners$endYear - TWPrisoners$timestep), TWPrisoners$timestep), 
                    each=TWPrisoners$ncomponent * TWPrisoners$npops))%>%
  filter(time==1500)%>%
  mutate(cascade_status = sub("^[^_]*_", "", Var2), 
         dis_prog = sub("\\_.*", "", Var2),
         SI = ifelse(cascade_status%in%c("s", "cured"), "S","I"),
         parameter =Var2)%>%group_by(Var1 ,SI)%>%
  mutate(total = sum(Freq),
         value = ifelse(Freq==0, 0, Freq/total))%>%
  ungroup()%>%group_by(Var1)%>%mutate(pop_prop = ifelse(
    Freq==0, 0, Freq/sum(Freq)))%>%
  ungroup()%>%dplyr::select(Var1,parameter, value, SI)



write.csv(popPro_extract, 
          file.path(DataFolder,"/Estimate_initial_pop.csv")) 

#### number of people in each population ####

estPops<- read.csv(file.path(DataFolder, "Estimate_initial_pop.csv"), 
                   header = TRUE)%>%dplyr::select(-"X")

init_pop <- filter(initialPops, parameter == "init_pop")$value

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4", 
                                                 "pop_prop5", "pop_prop6"))%>%
  dplyr::select(value)%>%unlist()%>%as.vector()

sum(pop_prop)
popProp <- as.numeric(init_pop)*pop_prop 

sum(popProp)
# prevalence at initial
init_prop_I <- c(constantsDf$HCVP1[1], constantsDf$HCVP2[1], 
                 constantsDf$HCVP3[1], constantsDf$HCVP4[1],
                 constantsDf$HCVP5[1], constantsDf$HCVP6[1])

init_prop_S <-c(1 - init_prop_I)

estPops <- estPops%>%mutate(
  pop_group = rep(c(popProp),dim(estPops)[1]/TWPrisoners$npops),
  SIprop = case_when(Var1 == "N_inca_NCID" & SI == "S" ~ init_prop_S[1],
                     Var1 == "N_inca_NCID" & SI == "I" ~ init_prop_I[1],
                     Var1 == "N_inca_CID" & SI == "S" ~ init_prop_S[2],
                     Var1 == "N_inca_CID" & SI == "I" ~ init_prop_I[2],
                     Var1 == "E_inca_NCID" & SI == "S" ~ init_prop_S[3],
                     Var1 == "E_inca_NCID" & SI == "I" ~ init_prop_I[3],
                     Var1 == "E_inca_CID" & SI == "S" ~ init_prop_S[4],
                     Var1 == "E_inca_CID" & SI == "I" ~ init_prop_I[4],
                     Var1 == "D_inca" & SI == "S" ~ init_prop_S[5],
                     Var1 == "D_inca" & SI == "I" ~ init_prop_I[5],
                     Var1 == "P_inca" & SI == "S" ~ init_prop_S[6],
                     Var1 == "P_inca" & SI == "I" ~ init_prop_I[6]),
  
  est_pop = value*pop_group*SIprop)

best_est_pop <- as.matrix(as.data.frame(matrix(estPops$est_pop, 
                                               ncol = TWPrisoners$ncomponent,  
                                               nrow = TWPrisoners$npops)))

colnames(best_est_pop) <- c(TWPrisoners$component_name)
sum(best_est_pop)



#### 3. simulate with new initial pops ####
tic <- proc.time()
endY <- 15
calibrateInit <- HCVMSM(TWPrisoners, best_estimates, best_est_pop,
                        disease_progress,  pop_array, dfList, fib,
                        modelrun="UN", proj = "TWPrisoners", end_Y = endY)


toc <- proc.time() - tic 
toc

#### 4. summarise the dt ####
##### Result:population #####
###### N: each subpop/ total  ######
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


# inflow ever incarcerated 
incarce_exp <- PPTranTo%>%
  filter(from %in% c("E_inca_NCID", "E_inca_CID") & to %in% c("P_inca"))%>%
  group_by(year)%>%summarise(best =sum(best))%>%
  mutate(year = TWPrisoners$cabY + year - 1)

incarce_all_bind<- dplyr::bind_rows(incarce_P_bind, incarce_D_bind, .id = 'population')%>%
  dplyr::select(year, population, best)%>%group_by(year)%>%
  summarise(best = sum(best))

incarce_exp_pro <- cbind(year = incarce_all_bind$year,
                              as_tibble(incarce_exp[, -1] / incarce_P_bind[ , -1])*100)%>%
  tibble::as_tibble()
#==============================================================================#

#==============================================================================#

# 5. plot with calib_dt 


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
  ggtitle("% of people in correctional settings")  + 
  scale_y_continuous(limits = c(0, 1))



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
  scale_y_continuous(limits = c(0, 600000)) 

subpop_N_plot[[4]] <- subpop_N_plot[[4]] + 
  scale_y_continuous(limits = c(0, 30000)) 

subpop_N_plot[[5]] <- subpop_N_plot[[5]] + 
  scale_y_continuous(limits = c(0, 5000)) + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "rehab"),aes(
    y = realPop,
    x = time), colour = "black") + 
  geom_segment(data = calib_dt[[1]]%>%filter(indicators == "rehab"),aes(
    y = low,yend = up,
    x= time, xend = time), colour = "black")

subpop_N_plot[[6]] <- subpop_N_plot[[6]] + 
  scale_y_continuous(limits = c(0,70000)) + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "prison"),aes(
    y = realPop,
    x = time), colour = "black") +
  geom_segment(data = calib_dt[[1]]%>%filter(indicators == "prison"),aes(
    y = low,yend = up,
    x= time, xend = time), colour = "black")

#####plot: % ####
# 
frac_prison_plot <- indicatorPlot(TWPrisoners, prison_pro, 
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
  scale_y_continuous(limits = c(0,15), breaks = seq(0, 15,5)) + 
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
  scale_y_continuous(limits = c(0,0.2), breaks = seq(0, 0.2,0.05)) + 
  ggtitle("PWID in community") 


#####plot: annual number entry to/release from detention #### 
incarce_D[["PWID"]] <- incarce_D[["PWID"]]%>%
  mutate( year = year + TWPrisoners$cabY - 1)
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
  geom_segment(data = calib_dt[[1]]%>%filter(indicators == "incar_DI"), aes(
    y = low,yend = up,
    x= time + TWPrisoners$cabY - 1, 
    xend = time + TWPrisoners$cabY - 1), colour = "black") +
  scale_y_continuous(limits = c(0, 5000))


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
    x = time + TWPrisoners$cabY - 1), colour = "black") + 
  geom_segment(data = calib_dt[[1]]%>%filter(indicators == "incar_D"),aes(
    y = low, yend = up,
    x= time + TWPrisoners$cabY - 1, 
    xend = time + TWPrisoners$cabY - 1), colour = "black")+ 
  scale_y_continuous(limits = c(0, 50000))

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
    x = time + TWPrisoners$cabY - 1), colour = "black") + 
  geom_segment(data = calib_dt[[1]]%>%filter(indicators == "incar_P"),aes(
    y = low, yend = up,
    x= time + TWPrisoners$cabY - 1, 
    xend = time + TWPrisoners$cabY - 1), colour = "black") + 
  scale_y_continuous(limits = c(20000, 50000))

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
    x = time + TWPrisoners$cabY - 1), colour = "black") + 
  geom_segment(data = calib_dt[[1]]%>%filter(indicators == "rel_D"),aes(
    y = low, yend = up,
    x= time + TWPrisoners$cabY - 1, 
    xend = time + TWPrisoners$cabY - 1), colour = "black")+ 
  scale_y_continuous(limits = c(0, 15000))

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
  ggtitle("Annual number release from Prisons") + theme_bw() + 
  geom_point(data = calib_dt[[1]]%>%filter(indicators == "rel_P"),aes(
    y = realPop,
    x = time + TWPrisoners$cabY - 1), colour = "black") + 
  geom_segment(data = calib_dt[[1]]%>%filter(indicators == "rel_P"),aes(
    y = low, yend = up,
    x= time + TWPrisoners$cabY - 1, 
    xend = time + TWPrisoners$cabY - 1), colour = "black")+ 
  scale_y_continuous(limits = c(0, 50000))


incarce_exp_pro_plot <- indicatorPlot(TWPrisoners, incarce_exp_pro,
                                      ylabel = "Number",
                                      xlimits = c(TWPrisoners$cabY,
                                                  TWPrisoners$cabY+30, 5),
                                      calibration_Y = TWPrisoners$cabY,
                                      rangeu = NULL,
                                      groupPlot = NULL,
                                      facetPlot = NULL,
                                      observationData = NULL, 
                                      simulateYear = NULL) +
  ggtitle("incarce_exp") + theme_bw() + 
  geom_point(data = calib_dt[[2]]%>%filter(indicator == "reincar_fit"),aes(
    y = realPop,
    x = time + TWPrisoners$cabY - 1), colour = "black") +
  scale_y_continuous(limits = c(0, 100, 10))

# 6. visually judgement

N_plot_incar  <- ggarrange(New_incarDI_plot, New_incarD_plot, New_incarP_plot,
                           New_releaseD_plot, New_releaseP_plot)
N_plot_subpop <- ggarrange(plotlist = subpop_N_plot)
frac_plot <- ggarrange(frac_commuP_plot, frac_D_plot, frac_prison_plot, 
                        incarce_exp_pro_plot)

N_plot_incar
N_plot_subpop
frac_plot 

New_releaseP_plot + scale_y_continuous(limits = c(0,50000))
subpop_N_plot[[6]] + scale_y_continuous(limits = c(0,120000))
subpop_N_plot[[3]] + scale_y_continuous(limits = c(200000,400000))
