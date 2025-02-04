gc()
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
OutputFig <- file.path(paste0(OutputFolder, "/Figs/PrevInc"))
# project specific code path 
Proj_code <- file.path(codefun_path, paste0("projects/", project_name))



load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "param_simulation.rda")))
load(file.path(OutputFolder, paste0(project_name, "Simulations.rda")))

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Proj_code, "/model_timestep.R")) 
# simulation outcomes 
# Sce_sq: pre-national program scenario 
# Sce_np: national program scenarios 

#### epi outcomes ####  
# HCV prevalence, incidence 
endY <- 100
indicator_flow <- Sce_sq[!names(Sce_sq)%in% c("allPops", "newpop_tran", 
                                              "newpop_tranState", "HCVdeathState",
                                              "newDeathState", "death_hcv", 
                                              "costPops", "QALYPops")]

subpop_N <- lapply(POC_AU$popNames, function(x){ 
  
  a <- popResults_MidYear(POC_AU, Sce_sq,
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
name_parset <- c("best", paste0("set", seq(1, 1000, 1)))


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


pop_state <- popResults_MidYear(POC_AU, Sce_sq,
                                Population = POC_AU$popNames,
                                Disease_prog = POC_AU$progress_name, 
                                Cascade = POC_AU$cascade_name, param = param_sq, 
                                endYear = endY)%>%ungroup()%>%
  as_tibble()

tempNOTInfected_subpop <- pop_state%>%filter(disease_prog!= "a")%>%
  filter(state == "s")%>%group_by(year, population)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)

tempChronic_subpop <- pop_state%>%filter(disease_prog!= "a")%>%
  group_by(year, population)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)

# arrange order to align with other dts 
pop_N <- pop_N%>%arrange(year, population)


tempPrev_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                         population = POC_AU$popNames,
                         
                         as.data.frame(100*(pop_N[, name_parset] - 
                                              tempNOTInfected_subpop[ ,name_parset])/ 
                                         pop_N[ ,name_parset]))%>%
  tibble::as_tibble()  


#prevalence in setting  
# tempPrev_setting[["setting]]

# community 
commu_N <- commu_N%>%arrange(year)

tempNOTInfected_commu <- pop_state%>%
  filter(state == "s" & population %in% c("C_PWID", "C_fPWID"))%>%group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)



tempPrev_setting <- list()

tempPrev_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(100*(commu_N[, name_parset] - 
                                                          tempNOTInfected_commu[ ,name_parset])/ 
                                                     commu_N[ ,name_parset]))%>%tibble::as_tibble()

# prison 
## PWID + former PWID + nonPWID 
prison_N <- prison_N%>%arrange(year)

tempNOTInfected_prison <- pop_state%>%
  filter(state == "s" & population %in% c("P_PWID", "P_fPWID", 
                                          "P_nPWID"))%>%group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)




tempPrev_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                       as.data.frame(100*(prison_N[, name_parset] - 
                                                            tempNOTInfected_prison[ ,name_parset])/ 
                                                       prison_N[ ,name_parset]))%>%tibble::as_tibble()


# prison_PWID experienced  
## PWID + former PWID + nonPWID 
prisonPWID_N <- prisonPWID_N%>%arrange(year)

tempNOTInfected_prisonPWID <- pop_state%>%
  filter(state == "s" & population %in% c("P_PWID", "P_fPWID"))%>%group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)


tempPrev_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                           as.data.frame(100*(prisonPWID_N[, name_parset] - 
                                                                tempNOTInfected_prisonPWID[ ,name_parset])/ 
                                                           prisonPWID_N[, name_parset]))%>%tibble::as_tibble()



##### (B)RNA prevalence ##### 
#    (I) each subpop
#   (II) community
#   (III) Prison
#   (IV) Prison_f/PWID (P_fPWID, P_PWID)
tempNOTInfectedRNA_subpop <- pop_state%>%
  filter(cascade %in%  c("s", "cured" ) )%>%group_by(year, population)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)


tempPrevRNA_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                            population = POC_AU$popNames,
                            
                            as.data.frame(100*(pop_N[, name_parset] - 
                                                 tempNOTInfectedRNA_subpop[ ,name_parset])/ 
                                            pop_N[ ,name_parset]))%>%
  tibble::as_tibble()  
View(tempPrevRNA_subpop)
x <- popResults_range(POC_AU, tempPrevRNA_subpop, Population =  POC_AU$popNames,
                                                   Disease_prog = NULL , 
                                                   Cascade = NULL, end_Y = 100) 
View(x)
#prevalence in setting  
# tempPrev_setting[["setting]]


tempNOTInfectedRNA_commu <- pop_state%>%
  filter(cascade %in%  c("s", "cured" ) & population %in% c("C_PWID", "C_fPWID"))%>%
  group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)


tempPrevRNA_setting <- list()

tempPrevRNA_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                        as.data.frame(100*(commu_N[, name_parset] - 
                                                             tempNOTInfectedRNA_commu[ ,name_parset])/ 
                                                        commu_N[ ,name_parset]))%>%tibble::as_tibble()

# prison 
## PWID + former PWID + nonPWID 


tempNOTInfectedRNA_prison <- pop_state%>%
  filter(cascade %in%  c("s", "cured" ) & population %in% c("P_PWID", "P_fPWID", 
                                                            "P_nPWID"))%>%
  group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)




tempPrevRNA_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                          as.data.frame(100*(prison_N[, name_parset] - 
                                                               tempNOTInfectedRNA_prison[ ,name_parset])/ 
                                                          prison_N[ ,name_parset]))%>%tibble::as_tibble()


# prison_PWID experienced  
## PWID + former PWID + nonPWID 


tempNOTInfectedRNA_prisonPWID <- pop_state%>%
  filter(cascade %in%  c("s", "cured" ) & population %in% c("P_PWID", "P_fPWID"))%>%
  group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)


tempPrevRNA_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                              as.data.frame(100*(prisonPWID_N[, name_parset] - 
                                                                   tempNOTInfectedRNA_prisonPWID[ ,name_parset])/ 
                                                              prisonPWID_N[ ,name_parset]))%>%tibble::as_tibble()



##### HCV incidence ##### 
HCVInfect_subpop <- indicatorResults(POC_AU, Sce_sq, "newInfections", 
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
                         as.data.frame(100*HCVInfect_subpop_C[, name_parset]/
                                         pop_N_C[ ,name_parset]))
HCVInc_subpop_P <- cbind(year = HCVInfect_subpop_P$year,
                         population = HCVInfect_subpop_P$population,
                         as.data.frame(100*HCVInfect_subpop_P[, name_parset]/
                                         (2*pop_N_P[ ,name_parset])))

HCVInc_subpop <- rbind(HCVInc_subpop_C, HCVInc_subpop_P)%>%arrange(year, population)%>%
  tibble::as_tibble()


HCVInfect_setting <- list()

HCVInfect_setting[["commu"]] <- HCVInfect_subpop_C%>%group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

HCVInfect_setting[["prisons"]] <- HCVInfect_subpop_P%>%group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

HCVInfect_setting[["prisonsPWID"]] <- HCVInfect_subpop%>%
  filter(population %in% c("P_PWID", "P_fPWID"))%>%group_by(year)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)

HCVInc_setting <- list()

HCVInc_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                   as.data.frame(100*(HCVInfect_setting[["commu"]][ , name_parset]/ 
                                                        commu_N[ ,name_parset])))%>%
  tibble::as_tibble()
HCVInc_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(100*(HCVInfect_setting[["prisons"]][ , name_parset]/ 
                                                          (2*prison_N[ ,name_parset]))))%>%
  tibble::as_tibble()

HCVInc_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                         as.data.frame(100*(HCVInfect_setting[["prisonsPWID"]][ , name_parset]/ 
                                                              (2*prisonPWID_N[ ,name_parset]))))%>%
  tibble::as_tibble()

# reinfection
Cured <- pop_state%>%
  filter(cascade == "cured")%>%
  group_by(year, population)%>%
  summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%
  arrange(year, population)

pop_cured_C <- Cured%>%filter(population %in% c("C_PWID", "C_fPWID"))
pop_cured_P <- Cured%>%filter(!population %in% c("C_PWID", "C_fPWID"))

HCVInfectRE_subpop <- indicatorResults(POC_AU, Sce_sq, "newreinfection", 
                                       pop=POC_AU$popNames,
                                       paramR = param_sq, range = NULL,
                                       endY = endY)
HCVInfectRE_subpop_C <- HCVInfectRE_subpop%>%filter(population %in% c("C_PWID", "C_fPWID"))
HCVInfectRE_subpop_P <- HCVInfectRE_subpop%>%filter(!population %in% c("C_PWID", "C_fPWID")) 

HCVIncre_subpop_C <- cbind(year = HCVInfectRE_subpop_C$year,
                           population = HCVInfectRE_subpop_C$population,
                           as.data.frame(100*HCVInfectRE_subpop_C[, name_parset]/
                                           pop_cured_C[ ,name_parset]))
HCVIncre_subpop_P <- cbind(year = HCVInfectRE_subpop_P$year,
                           population = HCVInfectRE_subpop_P$population,
                           as.data.frame(100*HCVInfectRE_subpop_P[, name_parset]/
                                           (2*pop_cured_P[ ,name_parset])))

HCVIncre_subpop <- rbind(HCVIncre_subpop_C, HCVIncre_subpop_P)%>%arrange(year, population)%>%
  tibble::as_tibble()


primary_N  <- cbind(year = pop_N$year, 
                    population = pop_N$population,
                    as.data.frame(pop_N[ ,name_parset] - Cured[, name_parset]))
primary_N_C <- primary_N%>%filter(population %in% c("C_PWID", "C_fPWID"))
primary_N_P <- primary_N%>%filter(!population %in% c("C_PWID", "C_fPWID"))

HCVInfectp_subpop_C <- HCVInfect_subpop_C[, name_parset] - HCVInfectRE_subpop_C[, name_parset]
HCVInfectp_subpop_P <- HCVInfect_subpop_P[, name_parset] - HCVInfectRE_subpop_P[, name_parset]

HCVIncp_subpop_C <- cbind(year = primary_N_C$year,
                          population = primary_N_C$population,
                          as.data.frame(100*HCVInfectp_subpop_C[, name_parset]/
                                          primary_N_C[ ,name_parset]))
HCVIncp_subpop_P <- cbind(year = primary_N_P$year,
                          population = primary_N_P$population,
                          as.data.frame(100*HCVInfectp_subpop_P[, name_parset]/
                                          (2*primary_N_P[ ,name_parset])))

HCVIncp_subpop <- rbind(HCVIncp_subpop_C, HCVIncp_subpop_P)%>%arrange(year, population)%>%
  tibble::as_tibble()

save(tempPrev_subpop, tempPrev_setting,
     tempPrevRNA_subpop, tempPrevRNA_setting,
     HCVInc_subpop, HCVInc_setting, HCVInfect_subpop,
     HCVIncre_subpop, HCVIncp_subpop, HCVInfectRE_subpop, 
     file = file.path(OutputFolder,
                      paste0(project_name,"PrevInc_sq" ,".rda")))


rm(tempPrev_subpop, tempPrev_setting,
   tempPrevRNA_subpop, tempPrevRNA_setting,
   HCVInc_subpop, HCVInc_setting, 
   HCVIncre_subpop, HCVIncp_subpop, param_sq, HCVInfectRE_subpop,HCVInfect_subpop ) 
gc()



################################# scenarios ####################################


################################ scenarios #####################################
Sce_sq <- list()
indicator_flow <- list()
endY <- 100
subpop_N <- list()
pop_N <- list()
name_parset <- c()
total_N <- list()
commu_N <- list()
prison_N <- list()
prisonPWID_N <- list()
pop_state <- list()
tempNOTInfected_subpop <- list()
tempChronic_subpop <- list()
tempPrev_subpop <- list()
tempNOTInfected_commu <- list()
tempPrev_setting <- list()
tempNOTInfected_prison <- list()
tempNOTInfected_prisonPWID <- list()
tempNOTInfectedRNA_subpop <- list()
tempPrevRNA_subpop <- list()
tempNOTInfectedRNA_commu <- list() 
tempPrevRNA_setting <- list()
tempNOTInfectedRNA_commu <- list()
tempPrevRNA_setting <- list()
tempNOTInfectedRNA_prison <- list()
tempNOTInfectedRNA_prisonPWID <- list()
HCVInfect_subpop <- list()

tempPrev_subpop <- list() 
tempPrev_setting <- list()
tempPrevRNA_subpop  <- list()
tempPrevRNA_setting <- list()
HCVInc_subpop <- list() 
HCVInc_setting <- list() 
HCVIncre_subpop <- list() 
HCVIncp_subpop <- list()

for(i in names(Sce_np)){
  load(file.path(OutputFolder, paste0(project_name, "param_sc_",i ,".rda"))) 
  
  indicator_flow <- Sce_np[[i]][!names(Sce_np[[i]])%in% c("allPops", "newpop_tran", 
                                                "newpop_tranState", "HCVdeathState",
                                                "newDeathState", "death_hcv", 
                                                "costPops", "QALYPops")]
  
  subpop_N <- lapply(POC_AU$popNames, function(x){ 
    
    a <- popResults_MidYear(POC_AU, Sce_np[[i]],
                            Population = x,
                            Disease_prog = NULL, 
                            Cascade = NULL, param = param_scenario, 
                            endYear = endY)%>%ungroup() 
    
    a <- popResults_range(POC_AU, a, Population = x,
                          Disease_prog = NULL , 
                          Cascade = NULL, end_Y = 100) 
  })
  
  names(subpop_N) <- POC_AU$popNames
  
  pop_N <- dplyr::bind_rows(subpop_N, .id = 'population')
  name_parset <- colnames(pop_N)[c(3: 1003)]
  
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
  
  pop_state <- popResults_MidYear(POC_AU, Sce_np[[i]],
                                  Population = POC_AU$popNames,
                                  Disease_prog = POC_AU$progress_name, 
                                  Cascade = POC_AU$cascade_name, param = param_scenario, 
                                  endYear = endY)%>%ungroup()%>%
    as_tibble()
  
  tempNOTInfected_subpop <- pop_state%>%filter(disease_prog!= "a")%>%
    filter(state == "s")%>%group_by(year, population)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)
  
  tempChronic_subpop <- pop_state%>%filter(disease_prog!= "a")%>%
    group_by(year, population)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)
  
  # arrange order to align with other dts 
  pop_N <- pop_N%>%arrange(year, population)
  
  tempPrev_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                           population = POC_AU$popNames,
                           
                           as.data.frame(100*(pop_N[, name_parset] - 
                                                tempNOTInfected_subpop[ ,name_parset])/ 
                                           pop_N[ ,name_parset]))%>%
    tibble::as_tibble()  
  
  
  #prevalence in setting  
  # tempPrev_setting[["setting]]
  
  # community 
  commu_N <- commu_N%>%arrange(year)
  
  tempNOTInfected_commu <- pop_state%>%
    filter(state == "s" & population %in% c("C_PWID", "C_fPWID"))%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  tempPrev_setting <- list()
  
  tempPrev_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                       as.data.frame(100*(commu_N[, name_parset] - 
                                                            tempNOTInfected_commu[ ,name_parset])/ 
                                                       commu_N[ ,name_parset]))%>%tibble::as_tibble()
  
  # prison 
  ## PWID + former PWID + nonPWID 
  prison_N <- prison_N%>%arrange(year)
  
  tempNOTInfected_prison <- list()
  
  tempNOTInfected_prison <- pop_state%>%
    filter(state == "s" & population %in% c("P_PWID", "P_fPWID", 
                                            "P_nPWID"))%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  
  
  
  tempPrev_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                         as.data.frame(100*(prison_N[, name_parset] - 
                                                              tempNOTInfected_prison[ ,name_parset])/ 
                                                         prison_N[ ,name_parset]))%>%tibble::as_tibble()
  
  
  # prison_PWID experienced  
  ## PWID + former PWID + nonPWID 
  prisonPWID_N <- prisonPWID_N%>%arrange(year)
  
  tempNOTInfected_prisonPWID <- list()
  
  tempNOTInfected_prisonPWID <- pop_state%>%
    filter(state == "s" & population %in% c("P_PWID", "P_fPWID"))%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  
  tempPrev_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                             as.data.frame(100*(prisonPWID_N[, name_parset] - 
                                                                  tempNOTInfected_prisonPWID[ ,name_parset])/ 
                                                             prisonPWID_N[, name_parset]))%>%tibble::as_tibble()
  
  
  
  ##### (B)RNA prevalence ##### 
  #    (I) each subpop
  #   (II) community
  #   (III) Prison
  #   (IV) Prison_f/PWID (P_fPWID, P_PWID)
  
  tempNOTInfectedRNA_subpop <- list()
  
  tempNOTInfectedRNA_subpop <- pop_state%>%
    filter(cascade %in%  c("s", "cured" ) )%>%group_by(year, population)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year, population)
  
  tempPrevRNA_subpop <- list()
  tempPrevRNA_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                              population = POC_AU$popNames,
                              
                              as.data.frame(100*(pop_N[, name_parset] - 
                                                   tempNOTInfectedRNA_subpop[ ,name_parset])/ 
                                              pop_N[ ,name_parset]))%>%
    tibble::as_tibble()  
  
  
  #prevalence in setting  
  # tempPrev_setting[["setting]]
  
  tempNOTInfectedRNA_commu <- list()
  tempNOTInfectedRNA_commu <- pop_state%>%
    filter(cascade %in%  c("s", "cured" ) & population %in% c("C_PWID", "C_fPWID"))%>%
    group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  
  tempPrevRNA_setting <- list()
  
  tempPrevRNA_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                          as.data.frame(100*(commu_N[, name_parset] - 
                                                               tempNOTInfectedRNA_commu[ ,name_parset])/ 
                                                          commu_N[ ,name_parset]))%>%tibble::as_tibble()
  
  # prison 
  ## PWID + former PWID + nonPWID 
  
  tempNOTInfectedRNA_prison <- list()
  tempNOTInfectedRNA_prison <- pop_state%>%
    filter(cascade %in%  c("s", "cured" ) & population %in% c("P_PWID", "P_fPWID", 
                                                              "P_nPWID"))%>%
    group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  
  
  
  tempPrevRNA_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                            as.data.frame(100*(prison_N[, name_parset] - 
                                                                 tempNOTInfectedRNA_prison[ ,name_parset])/ 
                                                            prison_N[ ,name_parset]))%>%tibble::as_tibble()
  
  
  # prison_PWID experienced  
  ## PWID + former PWID + nonPWID 
  
  tempNOTInfectedRNA_prisonPWID <- list()
  
  tempNOTInfectedRNA_prisonPWID <- pop_state%>%
    filter(cascade %in%  c("s", "cured" ) & population %in% c("P_PWID", "P_fPWID"))%>%
    group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  
  tempPrevRNA_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                                as.data.frame(100*(prisonPWID_N[, name_parset] - 
                                                                     tempNOTInfectedRNA_prisonPWID[ ,name_parset])/ 
                                                                prisonPWID_N[ ,name_parset]))%>%tibble::as_tibble()
  
 
  
  ##### HCV incidence ##### 
  HCVInfect_subpop <- list()
  HCVInfect_subpop <- indicatorResults(POC_AU, Sce_np[[i]], "newInfections", 
                                       pop=POC_AU$popNames,
                                       paramR = param_scenario, range = NULL,
                                       endY = endY)
  HCVInfect_subpop_P <- list()
  HCVInfect_subpop_P <- HCVInfect_subpop%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))
  
  HCVInfect_subpop_C <- list()
  HCVInfect_subpop_C <- HCVInfect_subpop%>%
    filter(population %in% c("C_PWID", "C_fPWID"))
  
  pop_N_P <- list()
 
  pop_N_P <- pop_N%>%
    filter(population %in% c("P_PWID", "P_fPWID", "P_nPWID"))
  
  pop_N_C <- list()
  pop_N_C <- pop_N%>%
    filter(population %in% c("C_PWID", "C_fPWID"))
  
  HCVInc_subpop_C <- list()
  HCVInc_subpop_C <- cbind(year = HCVInfect_subpop_C$year,
                           population = HCVInfect_subpop_C$population,
                           as.data.frame(100*HCVInfect_subpop_C[, name_parset]/
                                           pop_N_C[ ,name_parset]))
  HCVInc_subpop_P <- list()
  HCVInc_subpop_P <- cbind(year = HCVInfect_subpop_P$year,
                           population = HCVInfect_subpop_P$population,
                           as.data.frame(100*HCVInfect_subpop_P[, name_parset]/
                                           (2*pop_N_P[ ,name_parset])))
  HCVInc_subpop <- list()
  HCVInc_subpop <- rbind(HCVInc_subpop_C, HCVInc_subpop_P)%>%arrange(year, population)%>%
    tibble::as_tibble()
  
  
  HCVInfect_setting <- list()
  
  HCVInfect_setting[["commu"]] <- HCVInfect_subpop_C%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  HCVInfect_setting[["prisons"]] <- HCVInfect_subpop_P%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  HCVInfect_setting[["prisonsPWID"]] <- HCVInfect_subpop%>%
    filter(population %in% c("P_PWID", "P_fPWID"))%>%group_by(year)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%arrange(year)
  
  HCVInc_setting <- list()
  
  HCVInc_setting[["commu"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                     as.data.frame(100*(HCVInfect_setting[["commu"]][ , name_parset]/ 
                                                          commu_N[ ,name_parset])))%>%
    tibble::as_tibble()
  HCVInc_setting[["prisons"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                       as.data.frame(100*(HCVInfect_setting[["prisons"]][ , name_parset]/ 
                                                            (2*prison_N[ ,name_parset]))))%>%
    tibble::as_tibble()
  
  HCVInc_setting[["prisonsPWID"]] <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                                           as.data.frame(100*(HCVInfect_setting[["prisonsPWID"]][ , name_parset]/ 
                                                                (2*prisonPWID_N[ ,name_parset]))))%>%
    tibble::as_tibble()
  
  # reinfection
  Cured <- list()
  Cured <- pop_state%>%
    filter(cascade == "cured")%>%
    group_by(year, population)%>%
    summarise(across(c(name_parset),~ sum(.x, na.rm = FALSE)))%>%
    arrange(year, population)
  pop_cured_C <- list()
  pop_cured_C <- Cured%>%filter(population %in% c("C_PWID", "C_fPWID"))
  pop_cured_P <- list()
  pop_cured_P <- Cured%>%filter(!population %in% c("C_PWID", "C_fPWID"))
  
  HCVInfectRE_subpop <- list()
  HCVInfectRE_subpop <- indicatorResults(POC_AU, Sce_np[[i]], "newreinfection", 
                                         pop=POC_AU$popNames,
                                         paramR = param_scenario, range = NULL,
                                         endY = endY)
  HCVInfectRE_subpop_C <- list()
  HCVInfectRE_subpop_C <- HCVInfectRE_subpop%>%filter(population %in% c("C_PWID", "C_fPWID"))
  HCVInfectRE_subpop_P <- list()
  HCVInfectRE_subpop_P <- HCVInfectRE_subpop%>%filter(!population %in% c("C_PWID", "C_fPWID")) 
  HCVIncre_subpop_C <- list()
  HCVIncre_subpop_C <- cbind(year = HCVInfectRE_subpop_C$year,
                             population = HCVInfectRE_subpop_C$population,
                             as.data.frame(100*HCVInfectRE_subpop_C[, name_parset]/
                                             pop_cured_C[ ,name_parset]))
  HCVIncre_subpop_P <- list()
  HCVIncre_subpop_P <- cbind(year = HCVInfectRE_subpop_P$year,
                             population = HCVInfectRE_subpop_P$population,
                             as.data.frame(100*HCVInfectRE_subpop_P[, name_parset]/
                                             (2*pop_cured_P[ ,name_parset])))
  HCVIncre_subpop <- list()
  HCVIncre_subpop <- rbind(HCVIncre_subpop_C, HCVIncre_subpop_P)%>%arrange(year, population)%>%
    tibble::as_tibble()
  
  
  primary_N  <- cbind(year = pop_N$year, 
                      population = pop_N$population,
                      as.data.frame(pop_N[ ,name_parset] - Cured[, name_parset]))
  primary_N_C <- primary_N%>%filter(population %in% c("C_PWID", "C_fPWID"))
  primary_N_P <- primary_N%>%filter(!population %in% c("C_PWID", "C_fPWID"))
  
  HCVInfectp_subpop_C <- HCVInfect_subpop_C[, name_parset] - HCVInfectRE_subpop_C[, name_parset]
  HCVInfectp_subpop_P <- HCVInfect_subpop_P[, name_parset] - HCVInfectRE_subpop_P[, name_parset]
  
  HCVIncp_subpop_C <- cbind(year = primary_N_C$year,
                            population = primary_N_C$population,
                            as.data.frame(100*HCVInfectp_subpop_C[, name_parset]/
                                            primary_N_C[ ,name_parset]))
  HCVIncp_subpop_P <- cbind(year = primary_N_P$year,
                            population = primary_N_P$population,
                            as.data.frame(100*HCVInfectp_subpop_P[, name_parset]/
                                            (2*primary_N_P[ ,name_parset])))
  
  HCVIncp_subpop <- rbind(HCVIncp_subpop_C, HCVIncp_subpop_P)%>%arrange(year, population)%>%
    tibble::as_tibble()
  
  save(tempPrev_subpop, tempPrev_setting,
       tempPrevRNA_subpop, tempPrevRNA_setting,
       HCVInc_subpop, HCVInc_setting, 
       HCVIncre_subpop, HCVIncp_subpop, HCVInfect_subpop, HCVInfectRE_subpop,
       file = file.path(OutputFolder,
                        paste0(project_name,"PrevInc_", i ,".rda")))
  
  
 gc() 
  
}

# load files in a list 
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

files <- list.files(OutputFolder, pattern = paste0(project_name,"PrevInc_"))

PrevInc <- Map(rda2list, file.path(OutputFolder, files))

name_file <- sub("POC_AUPrevInc_", "", files)

names(PrevInc) <- tools::file_path_sans_ext(name_file)

View(PrevInc$dfList_NP_2023$tempPrevRNA_subpop%>%mutate(year = year +2014)%>%filter(year %in% c(2022,2030)))
options(scipen = 999)

#dropping plot_dt file 
# names(PrevInc)[-7]
PrevInc_range <- list()

for(i in names(PrevInc)[-7]){
  for( n in names(PrevInc[[1]])){ 
    if(n%in% c("tempPrevRNA_setting" ,"tempPrev_setting", "HCVInc_setting")){ 
      PrevInc_range[[i]][[n]] <- 
        lapply(PrevInc[[i]][[n]], function(x){ 
          popResults_range(POC_AU, x, Population = NULL, Disease_prog = NULL, 
                           Cascade = NULL, end_Y = endY - 1)}
          )
      
      names(PrevInc_range[[i]][[n]]) <- names(PrevInc[[i]][[n]])
      
      PrevInc_range[[i]][[n]] <- bind_rows(PrevInc_range[[i]][[n]], .id = 'setting')
    }
    else{ 
      PrevInc_range[[i]][[n]] <- popResults_range(POC_AU, PrevInc[[i]][[n]], 
                                                  Population = POC_AU$popNames, 
                                                  Disease_prog = NULL, 
                                                  Cascade = NULL, end_Y = endY - 1)
      }
  }
}

# turn list inside out 
PrevInc_range_bind <- PrevInc_range%>%purrr::transpose()%>%
  lapply(., function(x) dplyr::bind_rows(x, .id = 'scenario'))
names(PrevInc_range_bind) <- names(PrevInc_range[[1]])

write.xlsx(x = PrevInc_range_bind, file.path(OutputFolder, "PrevInc_epi.xlsx"), append = TRUE)

pop_labname <- c("PWID in community",  "Former PWID in community", 
                 "PWID in prisons",  "Former PWID in prisons", 
                 "nonPWID in prisons")

PrevInc_trajectory <- list()
# bind
for(i in names(PrevInc_range_bind)){
  PrevInc_range_bind[[i]] <- PrevInc_range_bind[[i]]%>%
    mutate(scenario = factor(scenario, levels = c("sq", "dfList_NP_2023", "dfList_NP_2024", 
                                                  "dfList_NPexp_A", "dfList_NPexp_B", "dfList_NPexp_C",
                                                  "dfList_NPexp_D"), 
                             labels = c("Pre national program", "Achievement 2023", 
                                        "Achievement 2024", "NP expand 2024", 
                                        "NP expand 2025", "NP expand 2026", 
                                        "NP expand 2027")))
  if(i %in% c("tempPrevRNA_setting" ,"tempPrev_setting", "HCVInc_setting")){ 
    
    PrevInc_range_bind[[i]] <-  
      PrevInc_range_bind[[i]]%>%
      mutate(population = factor(setting, 
                                 levels = c("commu", "prisons", "prisonsPWID"), 
                                 labels = c("Community", "Prisons", 
                                            "Current & former PWID in prisons")))
  }
  else{ 
    PrevInc_range_bind[[i]] <-  
      PrevInc_range_bind[[i]]%>%mutate(
        population = factor(population, levels = POC_AU$popNames, 
                            labels = pop_labname))
    
    }
  PrevInc_trajectory[[i]] <-  PrevInc_range_bind[[i]]%>%
    filter(scenario %in% c("Pre national program", "Achievement 2023"))
  
}
PrevInc_trajectory$tempPrevRNA_setting$population

# epi data for calibration 

files <- list.files(DataFolder%>%dirname(), pattern = "\\.csv$")
obserdt_name <- files[grep('HCVPrev', files)] 
obserdt_name <- c(obserdt_name, c(files[grep("HCVInc", files)]))

observedt_lst <- list() 

observedt_lst[[names(PrevInc_trajectory)[2]]] <- 
  read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrevRNA_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVPrev*100,
                           low = lower*100,
                           up = upper*100,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons" )))

observedt_lst[[names(PrevInc_trajectory)[3]]] <- 
  read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVIncP_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVInc,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))


observedt_lst[[names(PrevInc_trajectory)[5]]] <- 
  read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVInc_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVInc*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))

observedt_lst[[names(PrevInc_trajectory)[6]]] <- 
  read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrevRNA_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCV.RNA.prevalence*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))
  
observedt_lst[[names(PrevInc_trajectory)[7]]] <- 
  read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrev_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVPrev*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))

observedt_lst[[names(PrevInc_trajectory)[8]]] <- 
  read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVInc_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVInc*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))

observedt_lst[[names(PrevInc_trajectory)[9]]] <- 
  read.csv(file.path(paste0(DataFolder%>%dirname() ,"/HCVPrev_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCV.seroprevalence*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))

observedt_lst[[names(PrevInc_trajectory)[10]]] <- 
  read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVIncRe_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCVInc,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))
endY_plot <- 50 


# hex color code for Hep C website: "#1a979d"
xlimits <- c(POC_AU$startYear, POC_AU$startYear + 15, 5)
# simple function for plot generate 


PrevInc_plot <- function(pj, dt, obdt =NULL, xlimits, UI = NULL){ 
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
      geom_ribbon(aes(ymin = q5, ymax = q95, fill = scenario), alpha = 0.2) +
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = col_pal ) + 
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
      geom_ribbon(aes(ymin = q5, ymax = q95, fill = scenario), alpha = 0.2) +
      facet_wrap(~ population, scale ="free", ncol = 2 ) + 
      scale_color_manual(name = "Scenarios", values = col_pal ) + 
      scale_fill_manual(name = "Scenarios", values = col_pal ) + 
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


# function for getting ceiling number for plots 
lim_ident <- function(dt, group_index, year_range, summar_col){ 
  if(summar_col == "max"){
    if(group_index == "pop"){ 
      lim <- dt%>%group_by(population)%>%
        filter(year %in% year_range)%>%
        summarise(x = max(max))%>%
        mutate(lim = case_when( 
          x <1 ~5, 
          x>=1 & x <10 ~ 10, 
          x>=10 & x< 30 ~ 30, 
          x>=30 & x<60 ~ 60, 
          x>=60 & x<80 ~ 80, 
          x>80 & x<=100 ~100,
          x>100 & x<=1000 ~ (x%/%100 + 1)*100,
          x >1000 & x<10000 ~ (x%/%1000 + 1)*1000
        ))
    }
    else if(group_index == "setting"){ 
      lim <- dt%>%group_by(setting)%>%
        filter(year %in% year_range)%>%
        summarise(x = max(max))%>%
        mutate(lim = case_when( 
          x <1 ~5, 
          x>=1 & x <10 ~ 10, 
          x>=10 & x< 30 ~ 30, 
          x>=30 & x<60 ~ 60, 
          x>=60 & x<80 ~ 80, 
          x>80 & x<=100 ~100,
          x>100 & x<=1000 ~ (x%/%100 + 1)*100,
          x >1000 & x<10000 ~ (x%/%1000 + 1)*1000
        ))
    }
  return(lim)  
  }
  else if(summar_col == "best"){ 
    if(group_index == "pop"){ 
      lim <- dt%>%group_by(population)%>%
        filter(year %in% year_range)%>%
        summarise(x = max(best))%>%
        mutate(lim = case_when( 
          x <1 ~5, 
          x>=1 & x <10 ~ 10, 
          x>=10 & x< 30 ~ 30, 
          x>=30 & x<60 ~ 60, 
          x>=60 & x<80 ~ 80, 
          x>80 & x<=100 ~100,
          x>100 & x<=1000 ~ (x%/%100 + 1)*100,
          x >1000 & x<10000 ~ (x%/%1000 + 1)*1000
        ))
    }
    else if(group_index == "setting"){ 
      lim <- dt%>%group_by(setting)%>%
        filter(year %in% year_range)%>%
        summarise(x = max(best))%>%
        mutate(lim = case_when( 
          x <1 ~5, 
          x>=1 & x <10 ~ 10, 
          x>=10 & x< 30 ~ 30, 
          x>=30 & x<60 ~ 60, 
          x>=60 & x<80 ~ 80, 
          x>80 & x<=100 ~100,
          x>100 & x<=1000 ~ (x%/%100 + 1)*100,
          x >1000 & x<10000 ~ (x%/%1000 + 1)*1000
        ))
    }
    
  }

  return(lim)
} 

limx <- list()
for(i in c("tempPrevRNA_setting" ,"tempPrev_setting", "HCVInc_setting")){ 
    limx[[i]] <- 
      lim_ident(dt = PrevInc_trajectory[[i]], group_index = "setting", 
                year_range = seq(POC_AU$startYear, POC_AU$startYear + 15,1), 
                summar_col = "max")
}

for(i in c("HCVIncp_subpop","HCVInc_subpop","tempPrevRNA_subpop",
           "tempPrev_subpop", "HCVIncre_subpop") ){ 
    limx[[i]] <- 
      lim_ident(dt = PrevInc_trajectory[[i]], group_index = "pop", 
                year_range = seq(POC_AU$startYear, POC_AU$startYear + 15,1), 
                summar_col = "max")
  }

PrevInc_p <- list()

ylab_PrevInc <- list("HCV RNA prevalence (%)",
                     "Incidence of primary HCV infeciton (100 PY)",
                     "Number of HCV reinfection",
                     "HCV incidence (100 PY)", 
                     "HCV RNA prevalecne (%)", 
                     "HCV seroprevalecne (%)", 
                     "HCV incidence (100 PY)", 
                     "HCV seroprevalence (%)", 
                     "HCV reinfection incidence (100 PY)") 

names(ylab_PrevInc) <- names(PrevInc_trajectory)[-c(1,10)]

for(i in names(PrevInc_trajectory)[-c(1,10)]){ 
  PrevInc_p[[i]] <- PrevInc_plot(pj = POC_AU, 
                                 dt = PrevInc_trajectory[[i]], 
                                 obdt = observedt_lst[[i]], 
                                 xlimits = c(1, 16, 5), 
                                 UI = "y") + 
    labs(x = "Year", y = ylab_PrevInc[[i]]) 

    
} 

for(i in c("tempPrevRNA_setting" ,"tempPrev_setting", "HCVInc_setting")){ 
  PrevInc_p[[i]] <- PrevInc_p[[i]] + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][1,"lim"])))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][2,"lim"])))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][3,"lim"]))))
                    )) 
    
    
    
}

for(i in c("HCVIncp_subpop","HCVInc_subpop","tempPrevRNA_subpop",
           "tempPrev_subpop", "HCVIncre_subpop") ){ 
  PrevInc_p[[i]] <- PrevInc_p[[i]] + 
    facet_custom (~population,
                  scales = "free", ncol = 2,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][1,"lim"])))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][2,"lim"])))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][3,"lim"])))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][4,"lim"])))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][5,"lim"]))))
                    )) 
}



for(i in names(PrevInc_p)){ 
  ggsave(file=file.path(OutputFig, paste0(i,".png")), 
         PrevInc_p[[i]], 
         width = 6, height = 8, bg = "white", dpi = 300)
}

#### reinfection number #### 

unique(PrevInc_trajectory$HCVInfect_subpop$scenario)


x <- rbind(PrevInc_trajectory$HCVInfect_subpop%>%mutate(indicator = "Total new infections"), 
      PrevInc_trajectory$HCVInfectRE_subpop%>%mutate(indicator = "Reinfections"))
ggplot(data = x%>%filter(scenario == "Achievement 2023")%>%
         mutate(year = year + POC_AU$cabY - 1), 
       aes(x = year, y = best, colour = indicator)) + 
  geom_line() + facet_wrap(.~population, scale = "free") + theme_bw() + 
  scale_x_continuous(limits = c(2015, 2051), breaks = seq(2015, 2050,5))


# other scenarios 
PrevInc_range_sce <- list()
PrevInc_sce_p <- list()
for(i in names(PrevInc_range_bind)){ 
  PrevInc_range_sce[[i]] <- PrevInc_range_bind[[i]]%>%
    filter(!scenario %in% c("NP expand 2025", "NP expand 2026"))
  
  PrevInc_sce_p[[i]] <- PrevInc_plot(pj = POC_AU, 
                                     dt = PrevInc_range_sce[[i]], 
                                     obdt = NULL, 
                                     xlimits = c(8, 28, 5), 
                                     UI = NULL) + 
    labs(x = "Year", y = ylab_PrevInc[[i]]) +
    guides(colour = guide_legend(override.aes = list(alpha = 5))) + 
    theme(legend.position = "right", legend.direction="vertical")
}

for(i in c("tempPrevRNA_setting" ,"tempPrev_setting", "HCVInc_setting")){ 
  PrevInc_sce_p[[i]] <- PrevInc_sce_p[[i]] + 
    facet_custom (~population,
                  scales = "free", ncol = 1,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][1,"lim"])))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][2,"lim"])))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][3,"lim"]))))
                    )) 
  
  
  
}

for(i in c("HCVIncp_subpop","HCVInc_subpop","tempPrevRNA_subpop",
           "tempPrev_subpop", "HCVIncre_subpop") ){ 
  PrevInc_sce_p[[i]] <- PrevInc_sce_p[[i]] + 
    facet_custom (~population,
                  scales = "free", ncol = 2,
                  scale_overrides = 
                    list(
                      scale_new(1,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][1,"lim"])))),
                      scale_new(2,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][2,"lim"])))),
                      
                      scale_new(3,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][3,"lim"])))),
                      scale_new(4,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][4,"lim"])))),
                      scale_new(5,
                                scale_y_continuous(limits = 
                                                     c(0, as.numeric(limx[[i]][5,"lim"]))))
                    )) 
}


PrevInc_sce_p$tempPrevRNA_setting <- PrevInc_sce_p$tempPrevRNA_setting + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 15))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 15))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 20)))
                  )) 
  
PrevInc_sce_p[[2]] <- PrevInc_sce_p[[2]] + 
  facet_custom (~population,
                scales = "free", ncol = 2,
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
                                                   c(0, 20)))))

PrevInc_sce_p[[3]] <- PrevInc_sce_p[[3]] + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 3))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 0.5))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 20))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 3))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 1))))) 
  
PrevInc_sce_p[[4]] <- PrevInc_sce_p[[4]] + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 500))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 500))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 500))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 500))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 10))))) 

PrevInc_sce_p[[5]] <- PrevInc_sce_p[[5]] + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 80))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 40))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 80)))
                  ))
 
PrevInc_sce_p[[6]] <- PrevInc_sce_p[[6]] + 
  facet_custom (~population,
                scales = "free", ncol = 1,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 1))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 5))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 10)))
                  ))

PrevInc_sce_p[[7]] <- PrevInc_sce_p[[7]] + 
  facet_custom (~population,
                scales = "free", ncol = 1,
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
                                                   c(0, 80)))))

PrevInc_sce_p[[8]] <- PrevInc_sce_p[[8]] + 
  facet_custom (~population,
                scales = "free", ncol = 2,
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
                                                   c(0, 10))))) 


for(i in names(PrevInc_sce_p)){ 
  ggsave(file=file.path(OutputFig, paste0(i,"_sce" ,".png")), 
         PrevInc_sce_p[[i]], 
         width = 9, height = 6, bg = "white", dpi = 300)
}

save(PrevInc_range_bind,
     file = file.path(OutputFolder,
                      paste0(project_name,"PrevInc_plot_dt" ,".rda"))) 


