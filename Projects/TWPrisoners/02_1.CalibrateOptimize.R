# optimize calibration 
# this script is for finding the initial condition using optimizing methods 
# ISSUES remaine: slow process 

rm(list = ls())

library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)

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
load(file.path(OutputFolder, paste0(project_name, "cali", ".rda")))

urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- FALSE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))


N <- 1000 # number of accepted particles 

maxrun <- 2500 # maxmium run to avoid the infinte loop 

epsilon <- 20 # epison value 

n_par <- 11 # how many parameters will be estimated 


res <- list() # empty matrix to restore results



# loading original dataframe 
pop_transitions <- read.csv(file.path(DataFolder, "population_transitions.csv"), 
                            header = TRUE)

pop_transitions <- pop_transitions[ , -1] # drop the row number

# get first row only 
pop_transitions <- pop_transitions[1,]

xx <- array(unlist(pop_transitions), c(TWPrisoners$npops, TWPrisoners$npops, 1))

xx_array <- aperm(xx, c(2, 1, 3))

colnames(xx_array) <- TWPrisoners$popNames

rownames(xx_array) <- TWPrisoners$popNames

# calibration y 
endY <- 10

tic <- proc.time()
maxrun <- 10
for(j in 1:maxrun){ 
  N_inca_NCIDXN_inca_CID <- runif(1,0, 0.00001)
  N_inca_NCIDXP_inca <- runif(1, 0, 0.0001)
  N_inca_CIDXN_inca_NCID <- runif(1, 0, 0.1)
  N_inca_CIDXD_inca <- runif(1, 0, 0.1)
  N_inca_CIDXP_inca <- runif(1, 0.1, 0.5)
  E_inca_NCIDXE_inca_CID <- runif(1, 0, 0.12)
  E_inca_NCIDXP_inca <- runif(1, 0.1, 0.8)
  E_inca_CIDXE_inca_NCID <- runif(1, 0, 0.1)
  
  E_inca_CIDXP_inca <- runif(1, 0.2, 0.8)
  
  
  D_incaXE_inca_NCID <- runif(1, 0.4, 0.5)
  P_incaXE_inca_NCID <- runif(1, 0, 0.2)
  
  xx_array["N_inca_NCID", "N_inca_CID", 1] <- N_inca_NCIDXN_inca_CID
  xx_array["N_inca_NCID", "P_inca", 1] <- N_inca_NCIDXP_inca
  xx_array["N_inca_CID", "N_inca_NCID", 1] <- N_inca_CIDXN_inca_NCID
  xx_array["N_inca_CID", "D_inca", 1] <- N_inca_CIDXD_inca
  xx_array["N_inca_CID", "P_inca", 1] <- N_inca_CIDXP_inca
  xx_array["E_inca_NCID", "E_inca_CID", 1] <- E_inca_NCIDXE_inca_CID
  
  
  xx_array["E_inca_NCID", "P_inca", 1] <- E_inca_NCIDXP_inca
  
  xx_array["E_inca_CID", "E_inca_NCID", 1] <- E_inca_CIDXE_inca_NCID
  xx_array["E_inca_CID", "P_inca", 1] <- E_inca_CIDXP_inca
  
  
  xx_array["D_inca", "E_inca_NCID", 1] <- D_incaXE_inca_NCID
  xx_array["P_inca", "E_inca_NCID", 1] <- P_incaXE_inca_NCID
  
  
  # extend to same length 
  pop_array  <- array(matrix(xx_array[ , ,1]), 
                      c(TWPrisoners$npops, TWPrisoners$npops,TWPrisoners$npts))
  
  
  
  
  #### finding steady stage ####
  # calibration 
  
  
  calibrateInit <- HCVMSM(TWPrisoners, best_estimates, best_est_pop,
                          disease_progress,  pop_array, dfList, fib,
                          modelrun="UN", proj = "TWPrisoners", end_Y = endY)
  

  # outcomes calculation 
  # tidy model outputs 
  # total pops 
  # subpops 
  pop_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                              Population = TWPrisoners$popNames,
                              Disease_prog = NULL, 
                              Cascade = NULL, param = NULL, 
                              endYear = endY)%>%ungroup()%>%
    dplyr::group_by(year, population)%>%summarise_at("best",sum)%>%tibble::as.tibble()
  
  commu_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                                Population = c("N_inca_NCID", "N_inca_CID",
                                               "E_inca_NCID", "E_inca_CID"),
                                Disease_prog = NULL, 
                                Cascade = NULL, param = NULL, 
                                endYear = endY)%>%ungroup()%>%
    dplyr::group_by(year)%>%
    summarise_at("best",sum)%>%tibble::as.tibble()
  
  prison_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                                 Population = c("P_inca", "D_inca"),
                                 Disease_prog = NULL, 
                                 Cascade = NULL, param = NULL, 
                                 endYear = endY)%>%ungroup()%>%
    dplyr::group_by(year)%>%
    summarise_at("best",sum)%>%tibble::as.tibble()
  
  
  commu_proP <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                      as_tibble(
                        sum(pop_N[pop_N$population =="N_inca_CID",3],
                            pop_N[pop_N$population =="E_inca_CID",3])/ commu_N[ ,-1])*100)%>%
    tibble::as_tibble()
  
  prison_pro <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                      as_tibble((prison_N[ , -1]) /sum(commu_N[ ,-1],prison_N[ ,-1])*100))%>%
    tibble::as_tibble()
  
  prison_proD <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                       as_tibble(pop_N[pop_N$population =="D_inca",3] / prison_N[ ,-1])*100)%>%
    tibble::as_tibble()
  
  PWID_profNINCA <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                        as_tibble(pop_N[pop_N$population =="N_inca_CID",3] /
                                    sum(pop_N[pop_N$population =="N_inca_CID",3], 
                                        pop_N[pop_N$population =="E_inca_CID",3]))*100)%>%
    tibble::as_tibble()
  
  output <- c(commu_proP[1,2],  prison_pro[1,2], prison_proD[1,2], PWID_profNINCA[1,2])
  
  res[[j]] <- c(N_inca_NCIDXN_inca_CID, N_inca_NCIDXP_inca, 
                N_inca_CIDXN_inca_NCID, N_inca_CIDXD_inca,N_inca_CIDXP_inca,
                E_inca_NCIDXE_inca_CID, E_inca_NCIDXP_inca,
                E_inca_CIDXE_inca_NCID,  E_inca_CIDXP_inca,
                D_incaXE_inca_NCID,  P_incaXE_inca_NCID, 
                commu_proP[ ,2],
                prison_pro[ ,2],
                prison_proD[ ,2], 
                PWID_profNINCA[ ,2]) # store result
  
  }
toc <- proc.time() - tic
toc

save(res, 
     file = file.path(OutputFolder,
                      paste0("poptransit", ".rda")))
resx[[2]]

# sampling and refill in 
#calculating distance 
resx <- list()

resx <- lapply(res, function(x){ 
  
  a <- cbind(year = seq(TWPrisoners$cabY, (TWPrisoners$cabY + endY -1), 1),
             commu_proP = x[[12]], prison_pro = x[[13]],
             prison_proD = x[[14]], PWID_profNINCA = x[[15]],
             commu_proP_fit = 0.09,
             prison_pro_fit = 0.27,prison_proD_fit = 3.59, PWID_profNINCA_fit = 13.3)
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
    calc_distance(x[, "prison_pro"], x[, "prison_pro_fit"]) + 
    calc_distance(x[, "prison_proD"], x[, "prison_proD_fit"]) 
  
})

distance_bind <- do.call("rbind", distance)%>%as_tibble%>%
  mutate(sample = row_number())%>%arrange(V1)

op_pop <- lapply(res, function(x){ 
  
  a <- do.call("cbind",x)
  
  b <- matrix(a[1, c(1:11)], nrow = 1)
  colnames(b) <- c("N_inca_NCIDXN_inca_CID", "N_inca_NCIDXP_inca", 
                   "N_inca_CIDXN_inca_NCID", "N_inca_CIDXD_inca","N_inca_CIDXP_inca",
                   "E_inca_NCIDXE_inca_CID", "E_inca_NCIDXP_inca",
                   "E_inca_CIDXE_inca_NCID",  "E_inca_CIDXP_inca",
                   "D_incaXE_inca_NCID",  "P_incaXE_inca_NCID")
  return(b)
}) 

save(op_pop, 
     file = file.path(OutputFolder,
                      paste0("pop_test", ".rda")))
