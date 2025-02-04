gc()
rm(list = ls())
gc()
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
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev", ".rda")))
load(file.path(OutputFolder, paste0(project_name, "param_simulation", ".rda")))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 

# setting end Year 

endY <- 100
#### tidy up datapoint #### 

datapoint <- list()
datapoint[["N"]] <- 
  cbind.data.frame(year = rep(POC_AU$simY, length(POC_AU$popNames) +3), 
                   indicator = c(POC_AU$popNames, "Community", "Prison", "Total"), 
                   realpop = c(80000, 400000, 6993, 12696, 20311, 480000, 40000, 520000),
                   low= c(60000, 300000, 6163, 11819, 19332, 360000, 37314, 397314),
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



#### Result:population ####
#### N: each subpop/ total  #### 
# subpop as list
subpop_N <- lapply(POC_AU$popNames, function(x){ 
  
  a <- popResults_MidYear(POC_AU, calibrateInit,
                          Population = x,
                          Disease_prog = NULL, 
                          Cascade = NULL, param = param_sq, 
                          endYear = endY)%>%ungroup() 
  
  a <- popResults_range(POC_AU, a, Population = x,
                        Disease_prog = NULL , 
                        Cascade = NULL, end_Y = 100) 
})

names(subpop_N) <- POC_AU$popNames
# all subpop in one list 
pop_N <- dplyr::bind_rows(subpop_N, .id = 'population')

# colnames for the parameset and best estimation
name_parset <- colnames(pop_N)[c(3: 1003)]


# total number 

total_N <- pop_N%>%group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))

#### N: community ####
commu_N <- pop_N%>%filter(population %in% c("C_PWID", "C_fPWID"))%>%
  dplyr::group_by(year)%>%summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))

#### N: prison ####
prison_N <- pop_N%>%filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))%>%
  dplyr::group_by(year)%>%summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))

#### N: prisonPWID ####
prisonPWID_N <- pop_N%>%filter(population %in% c("P_PWID", "P_fPWID"))%>%
  dplyr::group_by(year)%>%summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))

#### %: PWID in community/prison ####
commu_proP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                    as.data.frame(pop_N[pop_N$population =="C_PWID",c(name_parset)] / commu_N[ ,c(name_parset)])*100)%>%
  tibble::as_tibble()

prison_proP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                     as.data.frame(pop_N[pop_N$population =="P_PWID",c(name_parset)] / prison_N[ ,c(name_parset)])*100)%>%
  tibble::as_tibble()

#### % fPWID in prison #### 
prison_profP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                      as.data.frame(pop_N[pop_N$population =="P_fPWID",c(name_parset)] / prison_N[ ,c(name_parset)])*100)%>%
  tibble::as_tibble()



#### number of annual leaving in each subpop ####
# this is for getting annual number of released in non-PWID in prison 
rel <- indicatorResults(POC_AU, calibrateInit, "newLeave", paramR = param_sq,
                        pop = POC_AU$popNames, endY= endY)%>%
  mutate(year = year + POC_AU$cabY - 1)

#### number of annual entry model in each subpop ####
# to get the annual number of incarceration in non-PWID in prison 
entry <- indicatorResults(POC_AU, calibrateInit, "newEntry", paramR = param_sq,
                          pop = POC_AU$popNames, endY= endY)%>%
  mutate(year = year + POC_AU$cabY - 1)

pop_tran <- list()
pop_tran[["best"]] <- calibrateInit$newpop_tran
parset_name <- c(paste0("set", seq(1,1000,1)))
for(i in 1: 1000){ 
  
  pop_tran[[parset_name[i]]] <- param_sq[[i]]$newpop_tran
  
  
}
PopTransition <- list()
PopTransition_all <- list()
for(name in names(pop_tran)){ 
  PopTransition[[name]] <- as.data.frame.table(pop_tran[[name]])%>%
    mutate(timestep = c(rep(seq(POC_AU$startYear, endY-POC_AU$timestep,
                                POC_AU$timestep),each = POC_AU$npops*POC_AU$npops)),
           from = Var1,
           To = Var2)%>%dplyr::select(-c(Var1, Var2, Var3))%>%
    filter(timestep != 1) 
  
  PopTransition_all[[name]] <- cbind.data.frame(timestep = PopTransition[[name]]$timestep, 
                                        from = PopTransition[[name]]$from,  
                                        to = PopTransition[[name]]$To, 
                                        best = PopTransition[[name]]$Freq)%>%
    as_tibble()%>%
    mutate(year = c(rep(seq(POC_AU$startYear, endY - 2, 1), 
                        each = (1/POC_AU$timestep)*POC_AU$npops*POC_AU$npops),
                    rep(endY - 1 , each = (1/POC_AU$timestep -1)*POC_AU$npops*POC_AU$npops)))
  

  }
name_flow_parset <- names(PopTransition_all)

val_col <- data.frame(year = PopTransition_all[[1]]$year,
                 timestep = PopTransition_all[[1]]$timestep,
                 from = PopTransition_all[[1]]$from,
                 to = PopTransition_all[[1]]$to,
                 lapply(PopTransition_all, function(x) x$best)%>%do.call(cbind, .))

PPTranTo <- val_col%>%
  group_by(year, from, to)%>%
  summarise(across(c(name_flow_parset ),~ sum(.x, na.rm = FALSE)))


incarce <- list()

incarce[["PWID"]] <- PPTranTo%>%filter(from == "C_PWID" & to == "P_PWID")

incarce[["fPWID"]] <- PPTranTo%>%filter(from == "C_fPWID" & to == "P_fPWID")

incarce_bind <- dplyr::bind_rows(incarce, .id = 'population')%>%
  mutate(year = POC_AU$cabY + year - 1)%>%ungroup()%>%
  dplyr::select(year, population, name_flow_parset)

entry_nonPWID <- entry%>%filter(population =="P_nPWID") 

# total incarceration = incarceration in C_PWID + incarceration in C_fPWID + incarceration in non-PWID (entry:non-PWID)
incar_total <- rbind(incarce_bind, entry_nonPWID)%>%group_by(year)%>%
  summarise(across(c(name_flow_parset ),~ sum(.x, na.rm = FALSE)))

release <- list()

release[["PWID"]] <- PPTranTo%>%filter(from == "P_PWID" & to == "C_PWID")

release[["fPWID"]] <- PPTranTo%>%filter(from == "P_fPWID" & to == "C_fPWID")

release_bind <- dplyr::bind_rows(release, .id = 'population')%>%
  mutate(year = POC_AU$cabY + year - 1)%>%ungroup()%>%
  dplyr::select(year, population, name_flow_parset)

release_nonPWID <- rel%>%filter(population =="P_nPWID") 

release_total <- rbind(release_bind, release_nonPWID)%>%group_by(year)%>%
  summarise(across(c(name_flow_parset ),~ sum(.x, na.rm = FALSE)))

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
        as.data.frame(inj_stop[["community"]][ , name_flow_parset ] / 
                        subpop_N[["C_PWID"]][ ,name_flow_parset ])*100)%>%
  tibble::as_tibble()


inj_stop_inc[["prison"]] <- 
  cbind(year = seq(POC_AU$startYear , endY-1 ,1),
        as.data.frame(inj_stop[["prison"]][ , name_flow_parset ] / 
                        subpop_N[["P_PWID"]][ ,name_flow_parset ])*100)%>%
  tibble::as_tibble()

inj_stop_inc_bind <- dplyr::bind_rows(inj_stop_inc, .id = 'population')%>%
  mutate(year = POC_AU$cabY + year - 1)%>%ungroup()
####  incidence of injection relapse #### 
inj_relap_inc <- list()

inj_relap_inc[["community"]] <- 
  cbind(year = seq(POC_AU$startYear , endY-1 ,1),
        as.data.frame(inj_relap[["community"]][ , name_flow_parset ] / 
                        subpop_N[["C_fPWID"]][ ,name_flow_parset ])*100)%>%
  tibble::as_tibble()

inj_relap_inc[["prison"]] <- 
  cbind(year = seq(POC_AU$startYear , endY-1 ,1),
        as.data.frame(inj_relap[["prison"]][ , name_flow_parset ] / 
                        subpop_N[["P_fPWID"]][ ,name_flow_parset ])*100)%>%
  tibble::as_tibble()

inj_relap_inc_bind <- dplyr::bind_rows(inj_relap_inc, .id = 'population')%>%
  mutate(year = POC_AU$cabY + year - 1)%>%ungroup()

#### incidence of incarceration  ####
incar_inc <- list()

incar_inc[["PWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                             as.data.frame(incarce[["PWID"]][ , name_flow_parset ] / 
                                             subpop_N[["C_PWID"]][ ,name_flow_parset ])*100)%>%
  tibble::as_tibble()

incar_inc[["fPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                              as.data.frame(incarce[["fPWID"]][ , name_flow_parset ] / 
                                              subpop_N[["C_fPWID"]][ ,name_flow_parset ])*100)%>%
  tibble::as_tibble()

incar_inc_bind <- dplyr::bind_rows(incar_inc, .id = 'population')%>%
  mutate(year = POC_AU$cabY + year - 1)%>%ungroup()

#### incidence of release  ####
release_inc <- list()

release_inc[["PWID"]] <- 
  cbind(year = seq(POC_AU$startYear , endY-1 ,1),
        as.data.frame(release[["PWID"]][ , name_flow_parset ] / 
                        subpop_N[["P_PWID"]][ ,name_flow_parset ])*100)%>%
  tibble::as_tibble()

release_inc[["fPWID"]] <- 
  cbind(year = seq(POC_AU$startYear , endY-1 ,1),
        as.data.frame(release[["fPWID"]][ , name_flow_parset ] /
                        subpop_N[["P_fPWID"]][ ,name_flow_parset ])*100)%>%
  tibble::as_tibble()

release_inc_bind <- dplyr::bind_rows(release_inc, .id = 'population')%>%
  mutate(year = POC_AU$cabY + year - 1)%>%ungroup()

# adding aggregate columns to each dtset for plots generation 


# plots for pops 


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
                                Cascade = POC_AU$cascade_name, param = param_sq, 
                                endYear = endY)%>%ungroup()%>%
  as_tibble()

tempNOTInfected_subpop <- pop_state%>%filter(disease_prog!= "a")%>%
  filter(state == "s")%>%group_by(year, population)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)

tempChronic_subpop <- pop_state%>%filter(disease_prog!= "a")%>%
  group_by(year, population)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)

# arrange order to align with other dts 
pop_N <- pop_N%>%arrange(year, population)


tempPrev_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                         population = POC_AU$popNames,
                         
                         as.data.frame(100*(pop_N[, name_flow_parset] - 
                                              tempNOTInfected_subpop[ ,name_flow_parset])/ 
                                         pop_N[ ,name_flow_parset]))%>%
  tibble::as_tibble()  


#prevalence in setting  
# tempPrev_setting[["setting]]

# community 
commu_N <- commu_N%>%arrange(year)

tempNOTInfected_commu <- pop_state%>%
  filter(state == "s" & population %in% c("C_PWID", "C_fPWID"))%>%group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  

tempPrev_setting <- list()

tempPrev_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(100*(commu_N[, name_flow_parset] - 
                                                          tempNOTInfected_commu[ ,name_flow_parset])/ 
                                                     commu_N[ ,name_flow_parset]))%>%tibble::as_tibble()

# prison 
## PWID + former PWID + nonPWID 
prison_N <- prison_N%>%arrange(year)

tempNOTInfected_prison <- pop_state%>%
  filter(state == "s" & population %in% c("P_PWID", "P_fPWID", 
                                          "P_nPWID"))%>%group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  


tempPrev_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                       as.data.frame(100*(prison_N[, name_flow_parset] - 
                                                            tempNOTInfected_prison[ ,name_flow_parset])/ 
                                                       prison_N[ ,name_flow_parset]))%>%tibble::as_tibble()


# prison_PWID experienced  
## PWID + former PWID + nonPWID 
prisonPWID_N <- prisonPWID_N%>%arrange(year)

tempNOTInfected_prisonPWID <- pop_state%>%
  filter(state == "s" & population %in% c("P_PWID", "P_fPWID"))%>%group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)


tempPrev_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                           as.data.frame(100*(prisonPWID_N[, name_flow_parset] - 
                                                                tempNOTInfected_prisonPWID[ ,name_flow_parset])/ 
                                                           prisonPWID_N[, name_flow_parset]))%>%tibble::as_tibble()



##### (B)RNA prevalence ##### 
#    (I) each subpop
#   (II) community
#   (III) Prison
#   (IV) Prison_f/PWID (P_fPWID, P_PWID)
tempNOTInfectedRNA_subpop <- pop_state%>%
  filter(cascade %in%  c("s", "cured" ) )%>%group_by(year, population)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)


tempPrevRNA_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                            population = POC_AU$popNames,
                            
                            as.data.frame(100*(pop_N[, name_flow_parset] - 
                                                 tempNOTInfectedRNA_subpop[ ,name_flow_parset])/ 
                                            pop_N[ ,name_flow_parset]))%>%
  tibble::as_tibble()  


#prevalence in setting  
# tempPrev_setting[["setting]]


tempNOTInfectedRNA_commu <- pop_state%>%
  filter(cascade %in%  c("s", "cured" ) & population %in% c("C_PWID", "C_fPWID"))%>%
  group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  

tempPrevRNA_setting <- list()

tempPrevRNA_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                        as.data.frame(100*(commu_N[, name_flow_parset] - 
                                                             tempNOTInfectedRNA_commu[ ,name_flow_parset])/ 
                                                        commu_N[ ,name_flow_parset]))%>%tibble::as_tibble()

# prison 
## PWID + former PWID + nonPWID 


tempNOTInfectedRNA_prison <- pop_state%>%
  filter(cascade %in%  c("s", "cured" ) & population %in% c("P_PWID", "P_fPWID", 
                                                          "P_nPWID"))%>%
  group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  


tempPrevRNA_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                          as.data.frame(100*(prison_N[, name_flow_parset] - 
                                                               tempNOTInfectedRNA_prison[ ,name_flow_parset])/ 
                                                          prison_N[ ,name_flow_parset]))%>%tibble::as_tibble()


# prison_PWID experienced  
## PWID + former PWID + nonPWID 


tempNOTInfectedRNA_prisonPWID <- pop_state%>%
  filter(cascade %in%  c("s", "cured" ) & population %in% c("P_PWID", "P_fPWID"))%>%
  group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)


tempPrevRNA_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                              as.data.frame(100*(prisonPWID_N[, name_flow_parset] - 
                                                                   tempNOTInfectedRNA_prisonPWID[ ,name_flow_parset])/ 
                                                              prisonPWID_N[ ,name_flow_parset]))%>%tibble::as_tibble()



##### HCV incidence ##### 
HCVInfect_subpop <- indicatorResults(POC_AU, calibrateInit, "newInfections", 
                                     pop=POC_AU$popNames,
                                     paramR = param_sq, range = NULL,
                                     endY = endY)

HCVInfect_subpop_P <- HCVInfect_subpop%>%
  filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))

HCVInfect_subpop_C <- HCVInfect_subpop%>%
  filter(population %in% c("C_PWID", "C_fPWID"))

pop_N_P <- pop_N%>%
  filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))

pop_N_C <- pop_N%>%
  filter(population %in% c("C_PWID", "C_fPWID"))


HCVInc_subpop_C <- cbind(year = HCVInfect_subpop_C$year,
                       population = HCVInfect_subpop_C$population,
                       as.data.frame(100*HCVInfect_subpop_C[, name_flow_parset]/
                                       pop_N_C[ ,name_flow_parset]))
HCVInc_subpop_P <- cbind(year = HCVInfect_subpop_P$year,
                         population = HCVInfect_subpop_P$population,
                         as.data.frame(100*HCVInfect_subpop_P[, name_flow_parset]/
                                         (2*pop_N_P[ ,name_flow_parset])))

HCVInc_subpop <- rbind(HCVInc_subpop_C, HCVInc_subpop_P)%>%arrange(year, population)%>%
  tibble::as_tibble()


HCVInfect_setting <- list()

HCVInfect_setting[["commu"]] <- HCVInfect_subpop_C%>%group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

HCVInfect_setting[["prisons"]] <- HCVInfect_subpop_P%>%group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

HCVInfect_setting[["prisonsPWID"]] <- HCVInfect_subpop%>%
  filter(population %in% c("P_PWID", "P_fPWID"))%>%group_by(year)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

HCVInc_setting <- list()

HCVInc_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                   as.data.frame(100*(HCVInfect_setting[["commu"]][ , name_flow_parset]/ 
                                                        commu_N[ ,name_flow_parset])))%>%
  tibble::as_tibble()
HCVInc_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(100*(HCVInfect_setting[["prisons"]][ , name_flow_parset]/ 
                                                          (2*prison_N[ ,name_flow_parset]))))%>%
  tibble::as_tibble()

HCVInc_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                         as.data.frame(100*(HCVInfect_setting[["prisonsPWID"]][ , name_flow_parset]/ 
                                                              (2*prisonPWID_N[ ,name_flow_parset]))))%>%
  tibble::as_tibble()

# reinfection
Cured <- pop_state%>%
  filter(cascade == "cured")%>%
  group_by(year, population)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year, population)

pop_cured_C <- Cured%>%filter(population %in% c("C_PWID", "C_fPWID"))
pop_cured_P <- Cured%>%filter(!population %in% c("C_PWID", "C_fPWID"))

HCVInfectRE_subpop <- indicatorResults(POC_AU, calibrateInit, "newreinfection", 
                                       pop=POC_AU$popNames,
                                       paramR = param_sq, range = NULL,
                                       endY = endY)
HCVInfectRE_subpop_C <- HCVInfectRE_subpop%>%filter(population %in% c("C_PWID", "C_fPWID"))
HCVInfectRE_subpop_P <- HCVInfectRE_subpop%>%filter(!population %in% c("C_PWID", "C_fPWID")) 

HCVIncre_subpop_C <- cbind(year = HCVInfectRE_subpop_C$year,
                         population = HCVInfectRE_subpop_C$population,
                         as.data.frame(100*HCVInfectRE_subpop_C[, name_flow_parset]/
                                         pop_cured_C[ ,name_flow_parset]))
HCVIncre_subpop_P <- cbind(year = HCVInfectRE_subpop_P$year,
                         population = HCVInfectRE_subpop_P$population,
                         as.data.frame(100*HCVInfectRE_subpop_P[, name_flow_parset]/
                                         (2*pop_cured_P[ ,name_flow_parset])))

HCVIncre_subpop <- rbind(HCVIncre_subpop_C, HCVIncre_subpop_P)%>%arrange(year, population)%>%
  tibble::as_tibble()

pop_N
primary_N  <- cbind(year = pop_N$year, 
                    population = pop_N$population,
                    as.data.frame(pop_N[ ,name_flow_parset] - Cured[, name_flow_parset]))
primary_N_C <- primary_N%>%filter(population %in% c("C_PWID", "C_fPWID"))
primary_N_P <- primary_N%>%filter(!population %in% c("C_PWID", "C_fPWID"))

HCVInfectp_subpop_C <- HCVInfect_subpop_C[, name_flow_parset] - HCVInfectRE_subpop_C[, name_flow_parset]
HCVInfectp_subpop_P <- HCVInfect_subpop_P[, name_flow_parset] - HCVInfectRE_subpop_P[, name_flow_parset]

HCVIncp_subpop_C <- cbind(year = primary_N_C$year,
                           population = primary_N_C$population,
                           as.data.frame(100*HCVInfectp_subpop_C[, name_flow_parset]/
                                           primary_N_C[ ,name_flow_parset]))
HCVIncp_subpop_P <- cbind(year = primary_N_P$year,
                           population = primary_N_P$population,
                           as.data.frame(100*HCVInfectp_subpop_P[, name_flow_parset]/
                                           (2*primary_N_P[ ,name_flow_parset])))

HCVIncp_subpop <- rbind(HCVIncp_subpop_C, HCVIncp_subpop_P)%>%arrange(year, population)%>%
  tibble::as_tibble()

#### plots #### 
pop_labname <- c("PWID in community",  "Former PWID in community", 
                 "PWID in prisons",  "Former PWID in prisons", 
                 "nonPWID in prisons")

####Results: death ####
#### HCV death ####

HCVdeath <- indicatorResults(POC_AU, calibrateInit, "newHCVdeaths", 
                             pop=POC_AU$popNames,
                             paramR = param_sq, range = NULL,
                             endY = endY)%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))

death <- indicatorResults(POC_AU, calibrateInit, "newDeath", 
                          pop=POC_AU$popNames,
                          paramR = param_sq, range = NULL,
                          endY = endY)%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))


#### Results: all component ####
pop_state <- popResults_MidYear(POC_AU, calibrateInit,
                                Population = POC_AU$popNames,
                                Disease_prog = POC_AU$progress_name, 
                                Cascade = POC_AU$cascade_name, param = param_sq, 
                                endYear = endY)%>%ungroup()

pop_stateAll <- pop_state%>%group_by(year, state)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

pop_stateSub <- pop_state%>%group_by(year, population, state)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)



pop_disproAll <- pop_state%>%group_by(year, disease_prog)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

pop_disproSub <- pop_state%>%group_by(year, population, disease_prog)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year, population) 


#### Results: Cascade ####

pop_cascaAll <- pop_state%>%group_by(year, cascade)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

pop_cascaSub <- pop_state%>%group_by(year, population, cascade)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year, population) 


#### flow
# select flow outcomes from the output list 
calibrateFlow <- calibrateInit[!names(calibrateInit)%in%
                                 c("allPops", "newpop_tran", "newpop_tranState", "HCVdeathState",
                                   "newDeathState", "death_hcv")]

flow_sub <- list()
flow_all <- list()

flow_sub <- lapply(names(calibrateFlow), function(x){ 
  a <- indicatorResults(POC_AU, calibrateFlow, x, 
                        pop=POC_AU$popNames,
                        paramR = param_sq, range = NULL,
                        endY = endY)
})

names(flow_sub) <- names(calibrateFlow)

flow_all <- lapply(flow_sub, function(x){ 
  
  a <- x%>%group_by(year)%>%
    summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%
    arrange(year)
}) 

names(flow_all) <- names(calibrateFlow)

testab <- pop_stateSub%>%
  filter(state%in% c("f0_undiag", "f1_undiag", "f2_undiag", "f3_undiag",
                     "f4_undiag", "dc_undiag", "hcc_undiag", "lt_undiag",
                     "plt_undiag"))%>%
  group_by(year,population)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year)




Tab_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                             each = POC_AU$npops),
                  population = POC_AU$popNames,
                  as.data.frame(100*flow_sub$newTestingAb[, name_flow_parset] / 
                                  testab[ , name_flow_parset]))%>%
  tibble::as_tibble()%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))


Tpoct_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                               each = POC_AU$npops),
                    population = POC_AU$popNames,
                    as.data.frame(100*flow_sub$newTestingPOCT[, name_flow_parset] / 
                                    testab[ ,name_flow_parset]))%>%
  tibble::as_tibble()%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))

testrna <- pop_stateSub%>%
  filter(state%in% c("f0_diag_ab", "f1_diag_ab", "f2_diag_ab", "f3_diag_ab",
                     "f4_diag_ab", "dc_diag_ab", "hcc_diag_ab", "lt_diag_ab",
                     "plt_diag_ab"))%>%
  group_by(year,population)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year)

Trna_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                              each = POC_AU$npops),
                   population = POC_AU$popNames,
                   as.data.frame(100*flow_sub$newTestingAg[, name_flow_parset] / 
                                   testrna[ , name_flow_parset]))%>%
  tibble::as_tibble()%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))

treatint <- pop_stateSub%>%
  filter(state%in% c("a_diag_RNA", "f0_diag_RNA", "f1_diag_RNA", "f2_diag_RNA", "f3_diag_RNA",
                     "f4_diag_RNA", "dc_diag_RNA", "hcc_diag_RNA", "lt_diag_RNA",
                     "plt_diag_RNA"))%>%
  group_by(year,population)%>%
  summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year)

Tint_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                              each = POC_AU$npops),
                   population = POC_AU$popNames,
                   as.data.frame(100*flow_sub$newTreatment[, name_flow_parset] / 
                                   treatint[ ,name_flow_parset]))%>%
  tibble::as_tibble()%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))


### treatment init 

# lifetime treatment init in C_PWID 


# treatment init annual numbers in prisons from 2016 to 2022 


# by setting 
flow_setting <- lapply(flow_sub, function(x){ 
  
  a <- x%>%
    mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), 
                            "commu", "prisons"))
  
  a <- a%>%group_by(year, setting)%>%
    summarise(across(c(name_flow_parset),~ sum(.x, na.rm = FALSE)))%>%
    mutate(population = setting)%>%select(-setting)
}) 

N_treatment <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                each = 2),
                     population = flow_setting$newTreatment[ , "population"],
                     as.data.frame(flow_setting$newTreatment[, name_flow_parset] + 
                                     flow_setting$newRetreat[, name_flow_parset] ))%>%
  tibble::as_tibble()%>%
  mutate(population = factor(population, 
                             levels = c("commu", "prisons"), 
                             labels = c("Community", "Prisons" )))


####Plot: population ####
#####plot: number of total pop  #####
total_N_range <- popResults_range(POC_AU, total_N, Population = NULL, 
                                      Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
totalPop_plot <- indicatorPlot(POC_AU, total_N_range, 
                               ylabel = "Number (thousands)",
                               xlimits = c(POC_AU$startYear, 
                                           POC_AU$startYear+30, 5),
                               calibration_Y = POC_AU$cabY,
                               rangeun = "y", 
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
        x = (POC_AU$simY - POC_AU$cabY + 1), 
        xend = (POC_AU$simY - POC_AU$cabY + 1)))
#####plot: number of pop in community/prison   #####
commu_N_range <- popResults_range(POC_AU, commu_N, Population = NULL, 
                                  Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)

commuPop_plot <- indicatorPlot(POC_AU, commu_N_range, 
                               ylabel = "Number",
                               xlimits = c(POC_AU$startYear, 
                                           POC_AU$startYear+30, 5),
                               calibration_Y = POC_AU$cabY,
                               rangeun = "y", 
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


prison_N_range <- popResults_range(POC_AU, prison_N, Population = NULL, 
                                  Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)


prisonPop_plot <- indicatorPlot(POC_AU, prison_N_range, 
                                ylabel = "Number",
                                xlimits = c(POC_AU$startYear, 
                                            POC_AU$startYear + 30, 5),
                                calibration_Y = POC_AU$cabY,
                                rangeun = "y", 
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

##### plot: number of subpop ####
# sub pop full names 
Namelab <- c("PWID in community", "former PWID in community",
             "PWID in prison", "former PWID in prison", "non-PWID in prison")

subpop_N_plot <- list() 

subpop_N_range <- list()
for(i in names(subpop_N)){ 
  subpop_N_range[[i]] <- popResults_range(POC_AU, subpop_N[[i]], Population = i,
                                          Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
  
  }

for(x in seq_along(names(subpop_N))){ 
  
  subpop_N_plot[[x]] <- indicatorPlot(POC_AU, subpop_N_range[[x]] , 
                                      ylabel = "Number",
                                      xlimits = c(POC_AU$startYear, 
                                                  POC_AU$startYear + 30, 5),
                                      calibration_Y = POC_AU$cabY,
                                      rangeun = "y", 
                                      groupPlot = NULL, 
                                      facetPlot = NULL,
                                      observationData = NULL, 
                                      simulateYear = NULL) +
    theme_bw() + ggtitle(paste0(Namelab[x])) 
  
}

subpop_N_plot[[1]] <- subpop_N_plot[[1]] + 
  scale_y_continuous(limits = c(0, 125000), breaks = seq(0, 125000, 5000)) + 
  geom_point(aes(
    y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[1], "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 5, size = 2) + 
  geom_segment(
    aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[1], "low"], 
        yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[1], "up"],         x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))

subpop_N_plot[[2]] <- subpop_N_plot[[2]] + 
  scale_y_continuous(limits = c(0, 600000), breaks = seq(0, 600000, 50000)) + 
  
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
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, 2500)) + 
  geom_point(aes(
    y=datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[5], "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 5, size = 2) + 
  geom_segment(
    aes(y = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[5], "low"], 
        yend = datapoint[["N"]][datapoint[["N"]][, "indicator"] == names(subpop_N)[5], "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))



#####plot: % of PWID/fPWID in community/prison  ####
# P_PWID 

prison_proP_range <- popResults_range(POC_AU, prison_proP, Population = NULL,
                                      Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
frac_PPWID_plot <- indicatorPlot(POC_AU, prison_proP_range , 
                                 ylabel = "Percentage (%)",
                                 xlimits = c(POC_AU$startYear, 
                                             POC_AU$startYear+30, 5),
                                 calibration_Y = POC_AU$cabY,
                                 rangeun = "y", 
                                 groupPlot = NULL, 
                                 facetPlot = NULL,
                                 observationData = NULL, 
                                 simulateYear = NULL) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,25), breaks = seq(0, 25,5)) + ggtitle("PWID in prison")
frac_PPWID_plot <- frac_PPWID_plot + 
  geom_point(aes(
    y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_proP_fit", "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape = 18) +
  geom_segment(
    aes(y = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_proP_fit", "low"], 
        yend = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_proP_fit", "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))

# P_fPWID
prison_profP_range <- popResults_range(POC_AU, prison_profP, Population = NULL,
                                      Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
frac_PfPWID_plot <- indicatorPlot(POC_AU, prison_profP_range, 
                                  ylabel = "Percentage (%)",
                                  xlimits = c(POC_AU$startYear, 
                                              POC_AU$startYear+30, 5),
                                  calibration_Y = POC_AU$cabY,
                                  rangeun = "y", 
                                  groupPlot = NULL, 
                                  facetPlot = NULL,
                                  observationData = NULL, 
                                  simulateYear = NULL) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,40), breaks = seq(0, 40,5)) + ggtitle("Former PWID in prison") + 
  geom_point(aes(
    y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_profP_fit", "realpop"], 
    x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black",shape = 18) + 
  geom_segment(
    aes(y = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_profP_fit", "low"], 
        yend = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "prison_profP_fit", "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))

# C_PWID 
commu_proP_range <- popResults_range(POC_AU, commu_proP, Population = NULL,
                                       Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)

frac_CPWID_plot <- indicatorPlot(POC_AU, commu_proP_range, 
                                 ylabel = "Percentage (%)",
                                 xlimits = c(POC_AU$startYear, 
                                             POC_AU$startYear+30, 5),
                                 calibration_Y = POC_AU$cabY,
                                 rangeun = "y", 
                                 groupPlot = NULL, 
                                 facetPlot = NULL,
                                 observationData = NULL, 
                                 simulateYear = NULL) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,25), breaks = seq(0, 25,5)) + ggtitle("PWID in community")

frac_CPWID_plot <- frac_CPWID_plot + geom_point(aes(
  y=datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "commu_proP_fit", "realpop"], 
  x = POC_AU$simY - POC_AU$cabY + 1 ), colour = "black", shape= 5, size = 0.6) + 
  geom_segment(
    aes(y = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "commu_proP_fit", "low"], 
        yend = datapoint[["frac"]][datapoint[["frac"]][, "indicator"] == "commu_proP_fit", "up"], 
        x = (POC_AU$simY - POC_AU$cabY + 1) , xend = (POC_AU$simY - POC_AU$cabY + 1)))


#####plot: annual number entry to/release from prisons #### 
incar_total_range <- popResults_range(POC_AU, incar_total, Population = NULL,
                                     Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)

New_incar_plot <-  indicatorPlot(POC_AU, incar_total_range,
                                 ylabel = "Number",
                                 xlimits = c(POC_AU$cabY,
                                             POC_AU$cabY+30, 5),
                                 calibration_Y = POC_AU$cabY,
                                 rangeu = "y",
                                 groupPlot = NULL,
                                 facetPlot = NULL,
                                 observationData = NULL, 
                                 simulateYear = NULL) +
  ggtitle("Annual number entry prisons") + theme_bw() + 
  scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 5000)) + 
  geom_point(aes(
    y=datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "incar", "realpop"], 
    x = POC_AU$simY ), colour = "black")

release_total_range <- popResults_range(POC_AU, release_total, Population = NULL,
                                      Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)


New_release_plot <-  indicatorPlot(POC_AU, release_total_range,
                                   ylabel = "Number",
                                   xlimits = c(POC_AU$cabY,
                                               POC_AU$cabY+30, 5),
                                   calibration_Y = POC_AU$cabY,
                                   rangeu = "y",
                                   groupPlot = NULL,
                                   facetPlot = NULL,
                                   observationData = NULL, 
                                   simulateYear = NULL) +
  ggtitle("Annual number release from prisons") + theme_bw() + 
  scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, 5000)) + 
  geom_point(aes(
    y=datapoint[["flow"]][datapoint[["flow"]][, "indicator"] == "release", "realpop"], 
    x = POC_AU$simY ), colour = "black")

#####plot: annual number stopping/relapse injection ####
# stopping injection
inj_stop_bind_range <- popResults_range(POC_AU, inj_stop_bind, Population = c("community", "prisons"),
                                        Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)

New_stopinj_plot <-  indicatorPlot(POC_AU, inj_stop_bind_range,
                                   ylabel = "Number",
                                   xlimits = c(POC_AU$cabY,
                                               POC_AU$cabY+30, 5),
                                   calibration_Y = POC_AU$cabY,
                                   rangeu = "y",
                                   groupPlot = NULL,
                                   facetPlot = population,
                                   observationData = NULL, 
                                   simulateYear = NULL) +
  ggtitle("Annual number stop injection") + theme_bw() + 
  scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 20000, 5000)) 
# relapse injection
inj_relap_bind_range <- popResults_range(POC_AU, inj_relap_bind, Population = c("community", "prisons"),
                                        Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)


New_relapinj_plot <-  indicatorPlot(POC_AU, inj_relap_bind_range,
                                    ylabel = "Number",
                                    xlimits = c(POC_AU$cabY,
                                                POC_AU$cabY+30, 5),
                                    calibration_Y = POC_AU$cabY,
                                    rangeu = "y",
                                    groupPlot = NULL,
                                    facetPlot = population,
                                    observationData = NULL, 
                                    simulateYear = NULL) +
  ggtitle("Annual number relapse injection") + theme_bw() + 
  scale_y_continuous(limits = c(0,2500), breaks = seq(0, 2500, 500)) 

##============================================================================##
#####plot: incidence of entry to/release from prisons #### 
incar_inc_bind_range <- popResults_range(POC_AU, incar_inc_bind, Population = c("PWID", "fPWID"),
                                         Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)

Inc_incar_plot <- indicatorPlot(POC_AU,incar_inc_bind_range ,
                                ylabel = "Incidence",
                                xlimits = c(POC_AU$cabY,
                                            POC_AU$cabY+30, 5),
                                calibration_Y = POC_AU$cabY,
                                rangeu = "y",
                                groupPlot = NULL,
                                facetPlot = population,
                                observationData = NULL, 
                                simulateYear = NULL) +
  ggtitle("Incidence of incarceration") + theme_bw() + 
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 1)) 

release_inc_bind_range <- popResults_range(POC_AU, release_inc_bind, Population = c("PWID", "fPWID"),
                                         Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)


Inc_release_plot <-  indicatorPlot(POC_AU, release_inc_bind_range,
                                   ylabel = "Incidence",
                                   xlimits = c(POC_AU$cabY,
                                               POC_AU$cabY+30, 5),
                                   calibration_Y = POC_AU$cabY,
                                   rangeu = "y",
                                   groupPlot = NULL,
                                   facetPlot = population,
                                   observationData = NULL, 
                                   simulateYear = NULL) +
  ggtitle("Incidence of released") + theme_bw() + 
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 10)) 

#####plot: incidence of stopping/relapse injection ####
# stopping injection
Inc_stopinj_plot <- list()


Inc_stopinj_plot <- lapply(inj_stop_inc, function(x){ 
  
  x_range <- popResults_range(POC_AU, x, Population = NULL,
                        Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
  
  
  a <- indicatorPlot(POC_AU, x_range,
                     ylabel = "Incidence",
                     xlimits = c(POC_AU$startYear ,
                                 POC_AU$startYear+30, 5),
                     calibration_Y = POC_AU$cabY,
                     rangeu = "y",
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
inj_relap_inc_bind_range <- popResults_range(POC_AU, inj_relap_inc_bind, 
                                             Population = c("community", "prisons"),
                                             Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)

Inc_relapinj_plot <-  indicatorPlot(POC_AU, inj_relap_inc_bind_range,
                                    ylabel = "Incidence",
                                    xlimits = c(POC_AU$cabY,
                                                POC_AU$cabY+30, 5),
                                    calibration_Y = POC_AU$cabY,
                                    rangeu = "y",
                                    groupPlot = NULL,
                                    facetPlot = population,
                                    observationData = NULL, 
                                    simulateYear = NULL) +
  ggtitle("Incidence of relapse injection") + theme_bw() + 
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 1))

endY_plot <- 2030- POC_AU$cabY

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

tempPrev_subpop_range <- popResults_range(POC_AU, tempPrev_subpop, 
                                             Population = POC_AU$popNames,
                                             Disease_prog = NULL, Cascade = NULL, 
                                          end_Y = endY-1)%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))
  


tempPrev_subpop_range%>%filter(population == "nonPWID in prisons")%>%select(q5, q95)
popPrevPlot <- indicatorPlot(POC_AU, tempPrev_subpop_range, 
                             ylabel = "HCV seroprevalence (%)",
                             xlimits = c(POC_AU$startYear, 
                                         (POC_AU$startYear+endY_plot), 5),
                             calibration_Y = POC_AU$cabY,
                             rangeun = "y", 
                             groupPlot = NULL, 
                             facetPlot = pop_labname,
                             observationData = HCVPrev, 
                             simulateYear = POC_AU$simY) +
  ggtitle("HCV seroprevalence by population") 

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
                                                   c(0, 20)))
                  )) + theme_bw()

popPrevPlot 
# setting 
tempPrev_setting_bind_range <- dplyr::bind_rows(tempPrev_setting, .id = 'population')%>%
  popResults_range(POC_AU, ., 
                   Population = c("commu", "prisons", "prisonsPWID"),
                   Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)



tempPrev_setting_bind_range <- tempPrev_setting_bind_range%>%
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


settingPrevPlot <- indicatorPlot(POC_AU, tempPrev_setting_bind_range, 
                                 ylabel = "HCV seroprevalence (%)",
                                 xlimits = c(POC_AU$startYear, 
                                             (POC_AU$startYear+endY_plot - 1), 5),
                                 calibration_Y = POC_AU$cabY,
                                 rangeun = "y", 
                                 groupPlot = NULL, 
                                 facetPlot = population,
                                 observationData = HCVPrev_setting, 
                                 simulateYear = POC_AU$simY) +
  ggtitle("HCV seroprevalence by setting") 

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
                                                   c(0, 100))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 100)))
                  )) + theme_bw() 

##### B: RNA prevalence #####
tempPrevRNA_subpop_range <- popResults_range(POC_AU, tempPrevRNA_subpop, 
                                          Population = POC_AU$popNames,
                                          Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)


tempPrevRNA_subpop_range <- tempPrevRNA_subpop_range%>%
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
popPrevRNAPlot <- indicatorPlot(POC_AU, tempPrevRNA_subpop_range, 
                                ylabel = "HCV prevalence (%)",
                                xlimits = c(POC_AU$startYear, 
                                            (POC_AU$startYear+endY_plot), 5),
                                calibration_Y = POC_AU$cabY,
                                rangeun = "y", 
                                groupPlot = NULL, 
                                facetPlot = population,
                                observationData = HCVPrevRNA , 
                                simulateYear = POC_AU$simY) +
  ggtitle("HCV RNA prevalence by population") 

popPrevRNAPlot <- popPrevRNAPlot + 
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
                                                   c(0, 20)))
                  )) + theme_bw()

popPrevRNAPlot 
# setting 
tempPrevRNA_setting_bind_range <- dplyr::bind_rows(tempPrevRNA_setting, .id = 'population')%>%
  popResults_range(POC_AU, ., 
                   Population = c("commu", "prisons", "prisonsPWID"),
                   Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
  
tempPrevRNA_setting_bind_range <- tempPrevRNA_setting_bind_range%>%
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


settingPrevRNAPlot <- indicatorPlot(POC_AU, tempPrevRNA_setting_bind_range, 
                                    ylabel = "HCV RNA prevalence (%)",
                                    xlimits = c(POC_AU$startYear, 
                                                (POC_AU$startYear+endY_plot - 1), 5),
                                    calibration_Y = POC_AU$cabY,
                                    rangeun = "y", 
                                    groupPlot = NULL, 
                                    facetPlot = population,
                                    observationData = HCVPrevRNAsetting, 
                                    simulateYear = POC_AU$simY) +
  ggtitle("HCV RNA prevalence by setting") 

settingPrevRNAPlot <- settingPrevRNAPlot + 
  facet_custom (~population,
                scales = "free", ncol = 1,
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
                                                   c(0, 100)))
                  )) + theme_bw()


settingPrevRNAPlot
#####plot: incidence #####
HCVInc_subpop_range <- popResults_range(POC_AU, HCVInc_subpop, 
                                        Population = POC_AU$popNames,
                                        Disease_prog = NULL, Cascade = NULL, 
                                        end_Y = endY-1)

HCVInc_subpop_range <- HCVInc_subpop_range%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname )) 
# primary infection
HCVIncp_subpop_range <- popResults_range(POC_AU, HCVIncp_subpop, 
                                        Population = POC_AU$popNames,
                                        Disease_prog = NULL, Cascade = NULL, 
                                        end_Y = endY-1)
HCVIncp_subpop_range <- HCVIncp_subpop_range%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))
# reinfection
HCVIncre_subpop_range <- popResults_range(POC_AU, HCVIncre_subpop, 
                                         Population = POC_AU$popNames,
                                         Disease_prog = NULL, Cascade = NULL, 
                                         end_Y = endY-1)
HCVIncre_subpop_range <- HCVIncre_subpop_range%>%
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

HCVIncp <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVIncP_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVInc,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))


HCVIncre <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVIncRe_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVInc,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))

popIncPlot <- indicatorPlot(POC_AU, HCVInc_subpop_range, 
                            ylabel = "HCV incidence",
                            xlimits = c(POC_AU$startYear, 
                                        (POC_AU$startYear+endY_plot), 5),
                            calibration_Y = POC_AU$cabY,
                            rangeun = "y", 
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
                                                   c(0, 25))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 25))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 40))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 25))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 10)))
                  )) + theme_bw()


popIncpPlot <- indicatorPlot(POC_AU, HCVIncp_subpop_range, 
                             ylabel = "HCV incidence",
                             xlimits = c(POC_AU$startYear, 
                                         (POC_AU$startYear+endY_plot), 5),
                             calibration_Y = POC_AU$cabY,
                             rangeun = "y", 
                             groupPlot = NULL, 
                             facetPlot = population,
                             observationData = HCVIncp, 
                             simulateYear = POC_AU$simY) +
  ggtitle("HCV primary incidence by population") 

popIncpPlot <- popIncpPlot + 
  facet_custom (~population,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 25))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 25))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 80))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 50))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 5)))
                  )) + theme_bw()

popIncrePlot <- indicatorPlot(POC_AU, HCVIncre_subpop_range, 
                              ylabel = "HCV incidence",
                              xlimits = c(POC_AU$startYear, 
                                          (POC_AU$startYear+endY_plot), 5),
                              calibration_Y = POC_AU$cabY,
                              rangeun = "y", 
                              groupPlot = NULL, 
                              facetPlot = population,
                              observationData = HCVIncre, 
                              simulateYear = POC_AU$simY) +
  ggtitle("HCV reinfection incidence by population") 

popIncrePlot <- popIncrePlot + 
  facet_custom (~population,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 20))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 5))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 80))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 20))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 80)))
                  )) + theme_bw()

# setting 
HCVInc_setting_bind_range <- dplyr::bind_rows(HCVInc_setting, .id = 'population')%>%
  popResults_range(POC_AU, ., 
                   Population = c("commu", "prisons", "prisonsPWID"),
                   Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)

HCVInc_setting_bind_range <-  HCVInc_setting_bind_range%>%mutate(population = factor(population, 
                             levels = c("commu", "prisons", "prisonsPWID"), 
                             labels = c("Community", "Prisons", 
                                        "Current & former PWID in prisons")))%>%
  as.data.frame()

HCVInc_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVInc_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVInc*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))

settingIncPlot <- indicatorPlot(POC_AU, HCVInc_setting_bind_range, 
                                ylabel = "HCV incidence",
                                xlimits = c(POC_AU$startYear, 
                                            (POC_AU$startYear+endY_plot - 1), 5),
                                calibration_Y = POC_AU$cabY,
                                rangeun = "y", 
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
                                                   c(0, 5))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 25))),
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 25))))) + theme_bw()

#### flow plots ####

pflow_sub <- lapply(names(flow_sub), function(x){ 
  x_range <- popResults_range(POC_AU, flow_sub[[x]], 
                              Population = POC_AU$popNames,
                              Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)
  a <- indicatorPlot(POC_AU, x_range, 
                     ylabel = "Numbers",
                     xlimits = c(POC_AU$startYear, 
                                 (POC_AU$startYear+endY_plot), 5),
                     calibration_Y = POC_AU$cabY,
                     rangeun = "y", 
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


HCVTinit <- read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVTinit_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = abCoverage*100,
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
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year, population)




Tab_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                             each = POC_AU$npops),
                  population = POC_AU$popNames,
                  as.data.frame(100*flow_sub$newTestingAb[, name_parset] / 
                                  testab[ ,name_parset]))%>%
  tibble::as_tibble()

Tab_frag_range <- popResults_range(POC_AU, Tab_frag, 
                                   Population = POC_AU$popNames,
                                   Disease_prog = NULL, Cascade = NULL, 
                                   end_Y = endY-1)%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))


Tpoct_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                               each = POC_AU$npops),
                    population = POC_AU$popNames,
                    as.data.frame(100*flow_sub$newTestingPOCT[, name_parset] / 
                                    testab[ ,name_parset]))%>%
  tibble::as_tibble()
  
Tpoct_frag_range <- popResults_range(POC_AU, Tpoct_frag, 
                                   Population = POC_AU$popNames,
                                   Disease_prog = NULL, Cascade = NULL, 
                                   end_Y = endY-1)%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))  


testrna <- pop_stateSub%>%
  filter(state%in% c("f0_diag_ab", "f1_diag_ab", "f2_diag_ab", "f3_diag_ab",
                     "f4_diag_ab", "dc_diag_ab", "hcc_diag_ab", "lt_diag_ab",
                     "plt_diag_ab"))%>%
  group_by(year,population)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year, population)

Trna_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                              each = POC_AU$npops),
                   population = POC_AU$popNames,
                   as.data.frame(100*flow_sub$newTestingAg[, name_parset] / 
                                   testrna[ ,name_parset]))%>%
  tibble::as_tibble()

Trna_frag_range <- popResults_range(POC_AU, Trna_frag, 
                                     Population = POC_AU$popNames,
                                     Disease_prog = NULL, Cascade = NULL, 
                                     end_Y = endY-1)%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))  

treatint <- pop_stateSub%>%
  filter(state%in% c("a_diag_RNA", "f0_diag_RNA", "f1_diag_RNA", "f2_diag_RNA", "f3_diag_RNA",
                     "f4_diag_RNA", "dc_diag_RNA", "hcc_diag_RNA", "lt_diag_RNA",
                     "plt_diag_RNA"))%>%
  group_by(year,population)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year, population)


Tint_frag <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                              each = POC_AU$npops),
                   population = POC_AU$popNames,
                   as.data.frame(100*flow_sub$newTreatment[, name_parset] / 
                                   treatint[ ,name_parset]))%>%
  tibble::as_tibble()

Tint_frag_range <- popResults_range(POC_AU, Tint_frag, 
                                    Population = POC_AU$popNames,
                                    Disease_prog = NULL, Cascade = NULL, 
                                    end_Y = endY-1)%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))  

# treatment 

Prevtreatinit_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                         each = POC_AU$npops),
                              population = POC_AU$popNames,
                              as.data.frame((flow_sub$newTreatment[, name_parset] + 
                                               flow_sub$newRetreat[, name_parset]) / 
                                              pop_N[ ,name_parset]*100))%>%
  tibble::as_tibble()

Prevtreatinit_subpop_range <- popResults_range(POC_AU, Prevtreatinit_subpop, 
                                               Population = POC_AU$popNames,
                                               Disease_prog = NULL, Cascade = NULL, 
                                               end_Y = endY-1)%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))


treatrecentPrev <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatrecentPrev_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = realpop*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))

treatrecentPrev_p <- indicatorPlot(POC_AU, Prevtreatinit_subpop_range, 
                                   ylabel = "%",
                                   xlimits = c(POC_AU$startYear, 
                                               (POC_AU$startYear+endY_plot), 5),
                                   calibration_Y = POC_AU$cabY,
                                   rangeun = "y", 
                                   groupPlot = NULL, 
                                   facetPlot = population,
                                   observationData = treatrecentPrev, 
                                   simulateYear = POC_AU$simY) + 
  ggtitle("treatment initiated in the past year among population") 

treatrecentPrev_p  <- treatrecentPrev_p  + 
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
                                                   c(0, 50))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 30))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 1)))
                  )) + theme_bw()
treatrecentPrev_p 



N_treatment_range <- popResults_range(POC_AU, N_treatment, Population = c("Community", "Prisons" ), 
                             Disease_prog = NULL, Cascade = NULL, end_Y = endY-1)


HCVtreatinitN_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))

endY_plot <- 30
N_treatment_setting_p <- indicatorPlot(POC_AU, N_treatment_range, 
                                       ylabel = "N",
                                       xlimits = c(POC_AU$startYear, 
                                                   (POC_AU$startYear+endY_plot), 5),
                                       calibration_Y = POC_AU$cabY,
                                       rangeun = "y", 
                                       groupPlot = NULL, 
                                       facetPlot = population,
                                       observationData = HCVtreatinitN_setting_fit, 
                                       simulateYear = POC_AU$simY) + 
  ggtitle("Number of treatment init") + theme_bw()

N_treatment_setting_p <- N_treatment_setting_p + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 50000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 5000)))))


subpop_N_plot[[1]] <- subpop_N_plot[[1]] + theme_Publication()
subpop_N_plot[[2]] <- subpop_N_plot[[2]] + theme_Publication()
subpop_N_plot[[3]] <- subpop_N_plot[[3]] + theme_Publication()
subpop_N_plot[[4]] <- subpop_N_plot[[4]] + theme_Publication()
subpop_N_plot[[5]] <- subpop_N_plot[[5]] + theme_Publication()

popPrevPlot <- popPrevPlot + theme_Publication_facet() 
popPrevRNAPlot <- popPrevRNAPlot + theme_Publication_facet() 
settingPrevPlot <- settingPrevPlot + theme_Publication_facet() 
settingPrevRNAPlot <- settingPrevRNAPlot + theme_Publication_facet()

popIncPlot <- popIncPlot + theme_Publication_facet()
settingIncPlot <- settingIncPlot + theme_Publication_facet() 
popIncpPlot <- popIncpPlot + theme_Publication_facet()
popIncrePlot <- popIncrePlot + theme_Publication_facet()
N_treatment_setting_p <- N_treatment_setting_p + theme_Publication_facet() + 
  ggtitle("Number of treatment initiation") + labs( y = "Numbers")
treatrecentPrev_p <- treatrecentPrev_p + theme_Publication_facet()

totalPop_plot <- totalPop_plot + theme_Publication()
commuPop_plot <- commuPop_plot + theme_Publication() + ggtitle("Community")
prisonPop_plot <- prisonPop_plot +theme_Publication() + ggtitle("Prisons")
frac_PPWID_plot <- frac_PPWID_plot + theme_Publication()
frac_PfPWID_plot <- frac_PfPWID_plot + theme_Publication()
frac_CPWID_plot <- frac_CPWID_plot + theme_Publication()
New_incar_plot <- New_incar_plot + theme_Publication()
New_release_plot <- New_release_plot + theme_Publication()
New_stopinj_plot <- New_stopinj_plot + theme_Publication_facet()
New_relapinj_plot <- New_relapinj_plot + theme_Publication_facet()
Inc_incar_plot <- Inc_incar_plot + theme_Publication_facet()
Inc_release_plot <- Inc_release_plot + theme_Publication_facet()

Inc_stopinj_plot_c <- Inc_stopinj_plot[["community"]] + theme_Publication() + 
  ggtitle("Incidence of stopping injection")
Inc_stopinj_plot_p <- Inc_stopinj_plot[["prison"]] + theme_Publication() + 
  ggtitle("Incidence of stopping injection")
Inc_relapinj_plot <- Inc_relapinj_plot + theme_Publication_facet()

Inc_stopinj_plot <- Inc_stopinj_plot + theme_Publication_facet()

save(subpop_N_plot, popPrevPlot, popPrevRNAPlot, settingPrevPlot, settingPrevRNAPlot, 
     popIncPlot, settingIncPlot, popIncpPlot, popIncrePlot, 
     N_treatment_setting_p, treatrecentPrev_p, 
     
     totalPop_plot, commuPop_plot, prisonPop_plot, frac_PPWID_plot, 
     frac_PfPWID_plot, frac_CPWID_plot, New_incar_plot, 
     New_release_plot, New_stopinj_plot, New_relapinj_plot, Inc_incar_plot, 
     Inc_release_plot, Inc_stopinj_plot_c, Inc_stopinj_plot_p, Inc_relapinj_plot,
     
     file = file.path(OutputFolder,
                      paste0(project_name, "calibra_plot_param.rda")))


