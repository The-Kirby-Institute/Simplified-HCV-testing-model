# optimize calibration 
# this script is for finding the initial condition using optimizing methods 
# ISSUES remaine: slow process 

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


urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- FALSE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))


N <- 10 # number of accepted particles 

maxrun <- 5000 # maxmium run to avoid the infinte loop 

epsilon <- 50 # epison value 

n_par <- 11 # how many parameters will be estimated 


res <- list() # empty matrix to restore results



# loading original dataframe 
pop_transitions <- read.csv(file.path(DataFolder, "population_transitions.csv"), 
                            header = TRUE)

pop_transitions <- pop_transitions[ , -1] # drop the row number

# get first row only 
pop_transitions <- pop_transitions[1,]

xx <- array(unlist(pop_transitions), c(POC_AU$npops, POC_AU$npops, 1))

xx_array <- aperm(xx, c(2, 1, 3))

colnames(xx_array) <- POC_AU$popNames

rownames(xx_array) <- POC_AU$popNames

# calibration y 
endY <- 10

tic <- proc.time()

for(j in 1:maxrun){ 
  C_PWIDXC_fPWID <- runif(1,0, 0.1)
  C_PWIDXP_PWID <- runif(1, 0, 0.25)
  C_fPWIDXC_PWID <- runif(1, 0, 0.005)
  C_fPWIDXP_fPWID <- runif(1, 0, 0.07)
  P_PWIDXC_PWID <- runif(1, 0, 0.6)
  P_PWIDXP_fPWID <- runif(1, 0.2, 0.8)
  P_fPWIDXC_fPWID <- runif(1, 0, 0.8)
  P_fPWIDXP_PWID <- runif(1, 0, 0.5)
  
  xx_array["C_PWID","C_fPWID", 1] <- C_PWIDXC_fPWID 
  xx_array["C_PWID","P_PWID", 1] <- C_PWIDXP_PWID 
  xx_array["C_fPWID","C_PWID", 1] <- C_fPWIDXC_PWID
  xx_array["C_fPWID","P_fPWID", 1] <- C_fPWIDXP_fPWID
  xx_array["P_PWID","C_PWID", 1] <- P_PWIDXC_PWID
  xx_array["P_PWID","P_fPWID", 1] <- P_PWIDXP_fPWID
  xx_array["P_fPWID","C_fPWID", 1] <- P_fPWIDXC_fPWID
  xx_array["P_fPWID","P_PWID", 1] <- P_fPWIDXP_PWID
  
  # extend to same length 
  pop_array  <- array(matrix(xx_array[ , ,1]), 
                      c(POC_AU$npops, POC_AU$npops,POC_AU$npts))
  
  
  
  
  #### finding steady stage ####

  
  steady <- HCVMSM(POC_AU, best_estimates, best_initial_pop,
                   disease_progress, pop_array, dfList, fib,
                   modelrun = "steady", proj = "POC_AU")

  
  df_list <- lapply(steady, as.data.frame.table)
  
  allpop <- df_list$allPops%>%mutate(time = rep(seq(1,(POC_AU$endYear - POC_AU$timestep),POC_AU$timestep), 
                                                each=POC_AU$ncomponent*POC_AU$npops),
                                     Frequency=round(Freq, digits = 3))
  
  popPro_extract <- df_list$allPops%>%
    mutate(time = rep(seq(POC_AU$startYear, (POC_AU$endYear - POC_AU$timestep), POC_AU$timestep), 
                      each=POC_AU$ncomponent * POC_AU$npops),
           Frequency=round(Freq, digits = 3))%>%
    filter(time == 1000)%>%
    mutate(cascade_status = sub("^[^_]*_", "", Var2), 
           dis_prog = sub("\\_.*", "", Var2),
           SI = ifelse(cascade_status%in%c("s", "cured"), "S","I"),
           parameter =Var2)%>%group_by(Var1 ,SI)%>%
    mutate(total = sum(Frequency),
           value = ifelse(Frequency==0, 0, round(Frequency/total, digits = 4)))%>%
    ungroup()%>%group_by(Var1)%>%mutate(pop_prop = ifelse(
      Frequency==0, 0, round(Frequency/sum(Frequency), digits = 4)))%>%
    ungroup()%>%select(Var1,parameter, value, SI)
  
  estPops <- popPro_extract
  
  init_pop <- filter(initialPops, parameter == "init_pop")$value
  
  pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                   "pop_prop3", "pop_prop4", 
                                                   "pop_prop5"))%>%
    select(value)%>%unlist()%>%as.vector()
  
  popProp <- as.numeric(init_pop)*pop_prop 
  
  
  # prevalence at initial
  init_prop_I <- c(constantsDf$HCVP1[1], constantsDf$HCVP2[1], 
                   constantsDf$HCVP3[1], constantsDf$HCVP4[1], constantsDf$HCVP5[1])
  
  init_prop_S <-c(1 - init_prop_I)
  
  estPops <- estPops%>%mutate(
    pop_group = rep(c(popProp),dim(estPops)[1]/POC_AU$npops),
    SIprop = case_when(Var1 == "C_PWID" & SI == "S" ~ init_prop_S[1],
                       Var1 == "C_PWID" & SI == "I" ~ init_prop_I[1],
                       Var1 == "C_fPWID" & SI == "S" ~ init_prop_S[2],
                       Var1 == "C_fPWID" & SI == "I" ~ init_prop_I[2],
                       Var1 == "P_PWID" & SI == "S" ~ init_prop_S[3],
                       Var1 == "P_PWID" & SI == "I" ~ init_prop_I[3],
                       Var1 == "P_fPWID" & SI == "S" ~ init_prop_S[4],
                       Var1 == "P_fPWID" & SI == "I" ~ init_prop_I[4], 
                       Var1 == "P_nPWID" & SI == "S" ~ init_prop_S[5],
                       Var1 == "P_nPWID" & SI == "I" ~ init_prop_I[5]
    ),
    est_pop = value*pop_group*SIprop)
  
  best_est_pop <- as.matrix(as.data.frame(matrix(estPops$est_pop, 
                                                 ncol = POC_AU$ncomponent,  
                                                 nrow = POC_AU$npops)))
  
  colnames(best_est_pop) <- c(POC_AU$component_name)
  
  
  # calibration 
  
  
  calibrateInit <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                          disease_progress,  pop_array, dfList, fib,
                          modelrun="UN", proj = "POC_AU", end_Y = endY)

  
  # outcomes calculation 
  # tidy model outputs 
  # total pops 
  # subpops 
  pop_N <- popResults_MidYear(POC_AU, calibrateInit,
                              Population = POC_AU$popNames,
                              Disease_prog = NULL, 
                              Cascade = NULL, param = NULL, 
                              endYear = endY)%>%ungroup()%>%
    dplyr::group_by(year, population)%>%summarise_at("best",sum)
  subpop_N <- lapply(POC_AU$popNames, function(x){ 
    
    a <- popResults_MidYear(POC_AU, calibrateInit,
                            Population = x,
                            Disease_prog = NULL, 
                            Cascade = NULL, param = NULL, 
                            endYear = endY)%>%ungroup()
  })
  
  names(subpop_N) <- POC_AU$popNames
  commu_N <- popResults_MidYear(POC_AU, calibrateInit,
                                Population = c("C_PWID", "C_fPWID"),
                                Disease_prog = NULL, 
                                Cascade = NULL, param = NULL, 
                                endYear = endY)%>%ungroup()%>%
    dplyr::group_by(year)%>%
    summarise_at("best",sum)
  
  prison_N <- popResults_MidYear(POC_AU, calibrateInit,
                                 Population = c("P_PWID", "P_fPWID", "P_nPWID"),
                                 Disease_prog = NULL, 
                                 Cascade = NULL, param = NULL, 
                                 endYear = endY)%>%ungroup()%>%
    dplyr::group_by(year)%>%
    summarise_at("best",sum)
  
  
  commu_proP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                      as_tibble(pop_N[pop_N$population =="C_PWID",3] / commu_N[ ,-1])*100)%>%
    tibble::as_tibble()
  
  prison_proP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                       as_tibble(pop_N[pop_N$population =="P_PWID",3] / prison_N[ ,-1])*100)%>%
    tibble::as_tibble()
  
  prison_profP <- cbind(year = seq(POC_AU$startYear , endY-1 ,1),
                        as_tibble(pop_N[pop_N$population =="P_fPWID",3] / prison_N[ ,-1])*100)%>%
    tibble::as_tibble()
  
  entry <- indicatorResult_uno(POC_AU, calibrateInit, "newEntry",
                               populations = POC_AU$popNames, endYear= endY) %>%
    mutate(year = year + POC_AU$cabY - 1)
  
  ##### annual number of incarceration/release/injection relapse/ stopping injection ##### 
  PopTransition <- as.data.frame.table(calibrateInit$newpop_tran)%>%
    mutate(timestep = c(rep(seq(POC_AU$startYear, endY-POC_AU$timestep,
                                POC_AU$timestep),each = POC_AU$npops*POC_AU$npops)),
           from = Var1,
           To = Var2)%>%select(-c(Var1, Var2, Var3))%>%ungroup()
  # giving the average number to first time step 
  impute <- PopTransition%>%filter(timestep >0 &timestep <2)%>%
    group_by(from, To)%>%
    mutate(total_pop = mean(Freq))%>%
    select(total_pop)%>%ungroup()
  
  impute <- impute[c(1:as.numeric(POC_AU$npops*POC_AU$npops)),]
  
  PopTransition[c(1:as.numeric(POC_AU$npops*POC_AU$npops)), "Freq"] <- impute$total_pop 
  
  PopTransition_all <- cbind.data.frame(timestep = PopTransition$timestep, 
                                        from = PopTransition$from,  
                                        to = PopTransition$To, 
                                        best = PopTransition$Freq)%>%
    as_tibble()%>%
    mutate(year = rep(rep(seq(1, endY-1,1),each = 1/POC_AU$timestep),
                      each = POC_AU$npops*POC_AU$npops))
  
  PPTranTo <- PopTransition_all%>%
    group_by(year, from, to)%>%summarise_at(.vars = "best", sum)
  PPTranTo%>%filter(from == "C_fPWID" & to == "C_PWID")
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
  
  release <- list()
  
  release[["PWID"]] <- PPTranTo%>%filter(from == "P_PWID" & to == "C_PWID")
  
  release[["fPWID"]] <- PPTranTo%>%filter(from == "P_fPWID" & to == "C_fPWID")
  
  release_bind <- dplyr::bind_rows(release, .id = 'population')%>%
    mutate(year = POC_AU$cabY + year - 1)%>%ungroup()%>%select(year, population, best)
  
  release_nonPWID <- rel%>%filter(population =="P_nPWID") 
  
  release_total <- rbind(release_bind, release_nonPWID)%>%group_by(year)%>%
    summarize(best = sum(best))
  
  inj_stop_inc <- list()
  
  inj_stop_inc[["community"]] <- 
    cbind(year = seq(POC_AU$startYear , endY-1 ,1),
          as.data.frame(inj_stop[["community"]][ , "best"] / 
                          subpop_N[["C_PWID"]][ ,"best"])*100)%>%
    tibble::as_tibble()
  
  n_row <- POC_AU$simY - POC_AU$cabY +1 

  
  res[[j]] <- c(C_PWIDXC_fPWID, C_PWIDXP_PWID, C_fPWIDXC_PWID, 
                C_fPWIDXP_fPWID, P_PWIDXC_PWID, P_PWIDXP_fPWID, 
                P_fPWIDXC_fPWID, P_fPWIDXP_PWID, 
                commu_proP[n_row ,2], 
                prison_proP[n_row ,2], 
                prison_profP[n_row ,2],
                subpop_N[[1]][n_row, 3],
                subpop_N[[2]][n_row, 3],
                subpop_N[[3]][n_row, 3],
                subpop_N[[4]][n_row, 3],
                subpop_N[[5]][n_row, 3],
                incar_total[n_row, 2],
                release_total[n_row, 2],
                inj_stop_inc[["community"]][n_row, 2]) # store result
  
  
  }
toc <- proc.time() - tic
toc

save(res, 
     file = file.path(OutputFolder,
                      paste0("poptransit", ".rda")))

# sampling and refill in 
#calculating distance 
resx <- list()
resx <- lapply(res, function(x){ 
  
  a <- cbind(year = seq(POC_AU$cabY, (POC_AU$cabY + endY -1), 1),
             commu_proP = x[[9]], prison_proP = x[[10]], 
             prison_profP = x[[11]], commu_proP_fit = 16.67,
             prison_proP_fit = 17.48, prison_profP_fit = 31.74, 
             commu_stopinj = 5.4)
  return(a)
})


# function to calculate sum of squared errors 
calc_distance <- function( D, D_star){
  # Define sum of squared errors
  # Vectorised function
  dist <- sqrt(sum((D - D_star ) ^ 2 ))
  return(dist)
}

distance <- lapply(resx, function(x){ 
  
  a <- calc_distance(x[, "commu_proP"], x[, "commu_proP_fit"]) + 
    calc_distance(x[, "prison_proP"], x[, "prison_proP_fit"]) + 
    calc_distance(x[, "prison_profP"], x[, "prison_profP_fit"]) + 
    calc_distance(x[, "commu_stopinj"], x[, "commu_stopinj_fit"])
  
  })

distance_bind <- do.call("rbind", distance)%>%as_tibble%>%
  mutate(sample = row_number())%>%arrange(V1)

op_pop <- lapply(res, function(x){ 
  
  a <- do.call("cbind",x)
  
  b <- matrix(a[1, c(1:8)], nrow = 1)
  colnames(b) <- c("C_PWIDXC_fPWID", "C_PWIDXP_PWID", "C_fPWIDXC_PWID", 
                   "C_fPWIDXP_fPWID", "P_PWIDXC_PWID", "P_PWIDXP_fPWID", 
                   "P_fPWIDXC_fPWID", "P_fPWIDXP_PWID")
  return(b)
  }) 

save(op_pop, 
     file = file.path(OutputFolder,
                      paste0("pop_test", ".rda")))






#####################finding here ##############################################


head(distance_bind,1)
# summary the data 
Hmisc::describe(reDmin)
# summary data with hist 
pop_postdist <- skim(resx[[370]][, c(1:4)])

save(pop_postdist, 
     file = file.path(OutputFolder,
                      paste0("pop_postdist", ".rda")))



#### mannual change the pop_array 


# test pop array # 

xx_array <- pop_array[ , ,1]
xx_array["C_PWID","C_fPWID"] <- op_pop[[370]][ , "C_PWIDXC_fPWID"]
xx_array["C_PWID","P_PWID"] <- op_pop[[370]][ , "C_PWIDXP_PWID"]
xx_array["C_fPWID","C_PWID"] <- op_pop[[370]][ , "C_fPWIDXC_PWID"]
xx_array["C_fPWID","P_fPWID"] <- op_pop[[370]][ , "C_fPWIDXP_fPWID"]
xx_array["P_PWID","C_PWID"] <- op_pop[[370]][ , "P_PWIDXC_PWID"]
xx_array["P_PWID","P_fPWID"] <- op_pop[[370]][ , "P_PWIDXP_fPWID"]
xx_array["P_fPWID","C_fPWID"] <- op_pop[[370]][ , "P_fPWIDXC_fPWID"]
xx_array["P_fPWID","P_PWID"] <- op_pop[[370]][ , "P_fPWIDXP_PWID"]

xx_array["C_PWID","C_fPWID"] <- 0.135
xx_array["C_fPWID","C_PWID"] <- 0.00025
xx_array["C_PWID","P_PWID"] <- 0.035
xx_array["P_fPWID","C_fPWID"] <- 0.3
xx_array["P_PWID","P_fPWID"] <- 0.75

pop_array  <- array(matrix(xx_array), 
                    c(POC_AU$npops, POC_AU$npops,POC_AU$npts))

#####