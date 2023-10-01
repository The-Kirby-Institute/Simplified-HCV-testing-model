# optimize calibration 
# this script is for finding the initial condition using optimizing methods 
# ISSUES remaine: slow process 

rm(list = ls())

library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library(readxl)
library(tmvtnorm)

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
library("readxl")
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
# calibration indicators 
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


N <- 1000 # number of accepted particles 

maxrun <- 50000 # maxmium run to avoid the infinte loop 

# Epsilon values for numbers
epsilon_N <- c(20000, 17500, 15000, 12500, 11000, 10000)  

# Epsilon values for percentage 
epsilon_P <- epsilon_N/1000  

G <- length(epsilon_N) # Number of generations

# Number of simulations for each parameter set
n <- 1

#  Lower and upper boundaries for priors
lm.low <-c(0,0,0,
          0, 0.1, 0.1, 
          0, 0, 0,
          0, 0.1, 0.3,
          0.2, 0.2)
lm.upp <-c(0.00001,0.00001,0.00001,
          0.1, 0.6, 0.6, 
          0.15, 0.00001, 0.5, 
          0.1, 0.8, 0.8,
          0.8, 0.8)


n_par <- 14 # how many parameters will be estimated 

# Empty matrices to store results 
res.old<-matrix(ncol=n_par,nrow=N)
res.new<-matrix(ncol=n_par,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)


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
endY <- 15

tic <- proc.time()

for(g in 1:G){  
  
  #Initiate counter
  i<-1    
  while(i <= N){ # While the number of accepted particles is less than N_particles
    if(g==1){
      # Sample from prior distributions 
      N_inca_NCIDXN_inca_CID <- runif(1,lm.low[1], lm.upp[1])
      N_inca_NCIDXD_inca <- runif(1, lm.low[2], lm.upp[2])
      N_inca_NCIDXP_inca <- runif(1, lm.low[3], lm.upp[3])
      
      N_inca_CIDXN_inca_NCID <- runif(1, lm.low[4], lm.upp[4])
      N_inca_CIDXD_inca <- runif(1, lm.low[5], lm.upp[5])
      N_inca_CIDXP_inca <- runif(1, lm.low[6], lm.upp[6])
      
      E_inca_NCIDXE_inca_CID <- runif(1, lm.low[7], lm.upp[7])
      E_inca_NCIDXD_inca <- runif(1, lm.low[8], lm.upp[8])
      E_inca_NCIDXP_inca <- runif(1, lm.low[9], lm.upp[9])
      
      E_inca_CIDXE_inca_NCID <- runif(1, lm.low[10], lm.upp[10])
      E_inca_CIDXD_inca <- runif(1, lm.low[11], lm.upp[11])
      E_inca_CIDXP_inca <- runif(1, lm.low[12], lm.upp[12])
      
      D_incaXE_inca_NCID <- runif(1, lm.low[13], lm.upp[13])
      
      P_incaXE_inca_NCID <- runif(1, lm.low[14], lm.upp[14])
    }
    else {
      #  Select particle from previous generation
      p <-sample(seq(1,N),1,prob=w.old)
      par <- rK(res.old[p,],sigma)
      N_inca_NCIDXN_inca_CID <- par[1]
      N_inca_NCIDXD_inca <- par[2]
      N_inca_NCIDXP_inca <- par[3]
      
      N_inca_CIDXN_inca_NCID <- par[4]
      N_inca_CIDXD_inca <- par[5]
      N_inca_CIDXP_inca <- par[6]
      
      E_inca_NCIDXE_inca_CID <- par[7]
      E_inca_NCIDXD_inca <- par[8]
      E_inca_NCIDXP_inca <- par[9]
      
      E_inca_CIDXE_inca_NCID <- par[10]
      E_inca_CIDXD_inca <- par[11]
      E_inca_CIDXP_inca <- par[12]
      
      D_incaXE_inca_NCID <- par[13]
      
      P_incaXE_inca_NCID <- par[14]
    }
    xx_array["N_inca_NCID", "N_inca_CID", 1] <- N_inca_NCIDXN_inca_CID
    xx_array["N_inca_NCID", "D_inca", 1] <- N_inca_NCIDXD_inca
    xx_array["N_inca_NCID", "P_inca", 1] <- N_inca_NCIDXP_inca
    
    xx_array["N_inca_CID", "N_inca_NCID", 1] <- N_inca_CIDXN_inca_NCID
    xx_array["N_inca_CID", "D_inca", 1] <- N_inca_CIDXD_inca
    xx_array["N_inca_CID", "P_inca", 1] <- N_inca_CIDXP_inca
    
    xx_array["E_inca_NCID", "E_inca_CID", 1] <- E_inca_NCIDXE_inca_CID
    xx_array["E_inca_NCID", "D_inca", 1] <- E_inca_NCIDXD_inca
    xx_array["E_inca_NCID", "P_inca", 1] <- E_inca_NCIDXP_inca
    
    xx_array["E_inca_CID", "E_inca_NCID", 1] <- E_inca_CIDXE_inca_NCID
    xx_array["E_inca_CID", "D_inca", 1] <- E_inca_CIDXD_inca
    xx_array["E_inca_CID", "P_inca", 1] <- E_inca_CIDXP_inca
    
    xx_array["D_inca", "E_inca_NCID", 1] <- D_incaXE_inca_NCID
    
    xx_array["P_inca", "E_inca_NCID", 1] <- P_incaXE_inca_NCID
    
    
    # extend to same length 
    pop_array  <- array(matrix(xx_array[ , ,1]), 
                        c(TWPrisoners$npops, TWPrisoners$npops,TWPrisoners$npts))
    #  Test if prior non zero
    if(prior.non.zero(c(N_inca_NCIDXN_inca_CID,
                        N_inca_NCIDXD_inca,
                        N_inca_NCIDXP_inca,
                    
                        N_inca_CIDXN_inca_NCID,
                        N_inca_CIDXD_inca,
                        N_inca_CIDXP_inca,
                        
                        E_inca_NCIDXE_inca_CID,
                        E_inca_NCIDXD_inca,
                        E_inca_NCIDXP_inca,
                        
                        E_inca_CIDXE_inca_NCID,
                        E_inca_CIDXD_inca,
                        E_inca_CIDXP_inca,
                        
                        D_incaXE_inca_NCID,
                        
                        P_incaXE_inca_NCID))) 
      # Set number of accepted simulations to zero
      m <- 0
      distance <-matrix(ncol=2,nrow=n)
    for(j in 1:n){
      calibrateInit <- HCVMSM(TWPrisoners, best_estimates, best_est_pop,
                              disease_progress,  pop_array, dfList, fib,
                              modelrun="UN", proj = "TWPrisoners", end_Y = endY)
      
      pop_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                                  Population = TWPrisoners$popNames,
                                  Disease_prog = NULL, 
                                  Cascade = NULL, param = NULL, 
                                  endYear = endY)%>%ungroup()%>%
        dplyr::group_by(year, population)%>%summarise_at("best",sum)%>%
        mutate(year = TWPrisoners$cabY + year - 1)%>%
        tibble::as.tibble()
      
      commu_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                                    Population = c("N_inca_NCID", "N_inca_CID",
                                                   "E_inca_NCID", "E_inca_CID"),
                                    Disease_prog = NULL, 
                                    Cascade = NULL, param = NULL, 
                                    endYear = endY)%>%ungroup()%>%
        dplyr::group_by(year)%>%
        summarise_at("best",sum)%>%mutate(year = TWPrisoners$cabY + year - 1)%>%
        tibble::as.tibble()
      
      prison_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                                     Population = c("P_inca"),
                                     Disease_prog = NULL, 
                                     Cascade = NULL, param = NULL, 
                                     endYear = endY)%>%ungroup()%>%
        dplyr::group_by(year)%>%
        summarise_at("best",sum)%>%
        mutate(year = TWPrisoners$cabY + year - 1)%>%
        tibble::as.tibble()%>%filter(year %in% seq(2019, 2022,1))
      
      rehab_N <- popResults_MidYear(TWPrisoners, calibrateInit,
                                    Population = c("D_inca"),
                                    Disease_prog = NULL, 
                                    Cascade = NULL, param = NULL, 
                                    endYear = endY)%>%ungroup()%>%
        dplyr::group_by(year)%>%
        summarise_at("best",sum)%>%
        mutate(year = TWPrisoners$cabY + year - 1)%>%
        tibble::as.tibble()%>%filter(year %in% seq(2019, 2022,1))
      
      
      commu_proP <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                          as_tibble(
                            sum(pop_N[pop_N$population =="N_inca_CID",3],
                                pop_N[pop_N$population =="E_inca_CID",3])/ commu_N[ ,-1])*100)%>%
        tibble::as_tibble()%>%
        mutate(year = TWPrisoners$cabY + year - 1)%>%
        filter(year == 2022)
      
      prison_pro <- cbind(year = seq(TWPrisoners$startYear , endY-1 ,1),
                          best = as_tibble(sum(prison_N[ , -1], rehab_N[, -1]) /sum(commu_N[ ,-1],prison_N[ ,-1], rehab_N[, -1])*100))%>%
        tibble::as_tibble()%>%mutate(year = TWPrisoners$cabY + year - 1)%>%
        rename(best = value)%>%
        filter(year == 2022)
      
      prison_proD <- cbind(year = rehab_N$year,
                           as_tibble(rehab_N[, -1] /(prison_N[ ,-1]+ rehab_N[, -1]))*100)%>%
        tibble::as_tibble()%>%
        filter(year == 2022)
      
      
      
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
      
      # inflow prison 
      
      incarce_P <- list()
      
      incarce_P[["nonPWID"]] <- PPTranTo%>%
        filter(from %in% c("N_inca_NCID", "E_inca_NCID") & to == "P_inca")
      
      incarce_P[["PWID"]] <- PPTranTo%>%
        filter(from %in% c("N_inca_CID", "E_inca_CID") & to == "P_inca")
      
      
      
      incarce_P_bind <- dplyr::bind_rows(incarce_P, .id = 'population')%>%
        mutate(year = TWPrisoners$cabY + year - 1)%>%ungroup()%>%
        dplyr::select(year, population, best)%>%group_by(year)%>%
        summarise(best = sum(best))
      
      
      
      # inflow rehab /rehab PWID 
      incarce_D <- list()
      
      incarce_D[["nonPWID"]] <- PPTranTo%>%
        filter(from %in% c("N_inca_NCID", "E_inca_NCID") & to == "D_inca")%>%
        group_by(year)%>%summarise(best = sum(best))%>%
        mutate(year = TWPrisoners$cabY + year - 1)
      
      incarce_D[["PWID"]] <- PPTranTo%>%
        filter(from %in% c("N_inca_CID", "E_inca_CID") & to == "D_inca")%>%
        group_by(year)%>%summarise(best = sum(best))%>%
        mutate(year = TWPrisoners$cabY + year - 1)
      
      incarce_D_bind <- dplyr::bind_rows(incarce_D, .id = 'population')%>%
        dplyr::select(year, population, best)%>%group_by(year)%>%
        summarise(best = sum(best))
      
      # outflow prison 
      release <- list()
      release[["P"]] <- PPTranTo%>%
        filter(from == "P_inca" & to %in% c("E_inca_NCID"))%>%
        mutate(year = TWPrisoners$cabY + year - 1)%>%
        filter(year%in% seq(2019, 2022,1))%>%ungroup()%>%dplyr::select(year, best)
      # outflow rehab
      release[["D"]] <- PPTranTo%>%
        filter(from == "D_inca" & to %in% c("E_inca_NCID"))%>%
        mutate(year = TWPrisoners$cabY + year - 1)%>%
        filter(year%in% seq(2019, 2022,1))%>%ungroup()%>%dplyr::select(year, best)
      
      # inflow ever incarcerated 
      incarce_exp <- PPTranTo%>%
        filter(from %in% c("E_inca_NCID", "E_inca_CID") & to %in% c("D_inca", "P_inca"))%>%
        group_by(year)%>%summarise(best =sum(best))%>%
        mutate(year = TWPrisoners$cabY + year - 1)
      
      incarce_all_bind <- dplyr::bind_rows(incarce_P_bind, incarce_D_bind, .id = 'population')%>%
        dplyr::select(year, population, best)%>%group_by(year)%>%
        summarise(best = sum(best))
      
      incarce_exp_pro <- cbind(year = incarce_all_bind$year,
            as_tibble(incarce_exp[, -1] / incarce_all_bind[ , -1])*100)%>%
        tibble::as_tibble()%>%
        filter(year == 2017)
      
    # filter years 
      incarce_P_bind <- incarce_P_bind%>%filter(year%in% seq(2019, 2022,1))
      incarce_D$PWID <- incarce_D$PWID%>%filter(year%in% seq(2019, 2022,1))
      incarce_D_bind <- incarce_D_bind%>%filter(year%in% seq(2019, 2022,1))
      
        
      output <- list(N = list(prison = prison_N,
                              rehab = rehab_N, 
                              incar_P = incarce_P_bind,
                              incar_DI = incarce_D$PWID,
                              incar_D = incarce_D_bind,
                              rel_P = release[["P"]],
                              rel_D = release[["D"]]),
                     frac = c(commu_proP[, 2] ,prison_pro[, 2], 
                              prison_proD[ ,2],incarce_exp_pro[, 2])
                     )
      
      # Calculate distances 
      calc.dist <- calc_distance(output, calib_dt)
      distance[j,] <-calc.dist    
      if((calc.dist[1] <= epsilon_N[g]) & (calc.dist[2] <= epsilon_P[g])){ # If both distances are less than their tolerances
        m<-m+1
      }
    }   
    if (m>0){
      # Store results
      res.new[i,]<-c(N_inca_NCIDXN_inca_CID,
                     N_inca_NCIDXD_inca,
                     N_inca_NCIDXP_inca,
                     
                     N_inca_CIDXN_inca_NCID,
                     N_inca_CIDXD_inca,
                     N_inca_CIDXP_inca,
                     
                     E_inca_NCIDXE_inca_CID,
                     E_inca_NCIDXD_inca,
                     E_inca_NCIDXP_inca,
                     
                     E_inca_CIDXE_inca_NCID,
                     E_inca_CIDXD_inca,
                     E_inca_CIDXP_inca,
                     
                     D_incaXE_inca_NCID,
                     
                     P_incaXE_inca_NCID)  
      # Calculate weights
      w1<-prod(sapply(1:14, function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])), na.rm = TRUE)
      if(g==1){
        w2<-1
      } else {
        w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
      }
      w.new[i] <- (m/n)*w1/w2
      # Update counter
      i <- i+1
      print(paste0('Generation: ', g, ", particle: ", i))
    }
  } 
  
  sigma <- cov(res.new) 
  
  
  res.old<-res.new
  
  w.old<-w.new/sum(w.new)
  
  write.csv(res.new, file = paste("results_case_2_ABC_SMC_gen_",g,".csv",sep=""), row.names=FALSE)


}



