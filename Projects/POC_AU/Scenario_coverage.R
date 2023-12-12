# this script is using the coverage to estimate the costs around test and treat in the scenarios 
# run PAF_coverage_cal.R before this script 
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

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
load(file.path(OutputFolder, paste0(project_name, "Sce_np.rda")))
load(file.path(OutputFolder, paste0(project_name, "S_NP_test.rda")))
load(file.path(OutputFolder, paste0(project_name, "S_NPscale_test.rda")))
load(file.path(OutputFolder, paste0(project_name, "Sce_npscale.rda")))

source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))

# coverage of cascade estimation 

S_Yint <- 2022
S_Yfir <- 2024
S_Ymid <- 2027
S_Yend <-2050
SYpoint_int <- (S_Yint - POC_AU$cabY)/POC_AU$timestep + 1
SYpoint_fir <- (S_Yfir - POC_AU$cabY)/POC_AU$timestep + 1
SYpoint_mid <- (S_Ymid - POC_AU$cabY)/POC_AU$timestep + 1
SYpoint_end <- (S_Yend - POC_AU$cabY)/POC_AU$timestep + 1
SY_leng <- SYpoint_end - SYpoint_int
SY_lengmid <- SYpoint_mid - SYpoint_int 
SY_lengfir <- SYpoint_fir - SYpoint_int

Ccal_C <- 0.01208

Ccal_P <- 0.05424902

Ccal_scale <- list("C" = Ccal_C*2, "P" = Ccal_P*2)  

# coverage change
# status qup 
coverage_sq <- matrix(0, ncol = POC_AU$npts, nrow = POC_AU$npops)

save(coverage_sq,
     file = file.path(OutputFolder ,
                      paste0(project_name,"coverage_sq" ,".rda")))
coverage_np <- coverage_sq
# nationa program 
coverage_np[1, SYpoint_int:POC_AU$npts] <- 
  c(seq(0, Ccal$C, length = (SY_lengfir)), rep(Ccal$C, (POC_AU$npts - SYpoint_fir +1)))
coverage_np[2, SYpoint_int:POC_AU$npts] <- 
  c(seq(0, Ccal$C, length = (SY_lengfir)), rep(Ccal$C, (POC_AU$npts - SYpoint_fir +1)))
coverage_np[3, SYpoint_int:POC_AU$npts] <- 
  c(seq(0, Ccal$P, length = (SY_lengfir)), rep(Ccal$P, (POC_AU$npts - SYpoint_fir +1)))
coverage_np[4, SYpoint_int:POC_AU$npts] <- 
  c(seq(0, Ccal$P, length = (SY_lengfir)), rep(Ccal$P, (POC_AU$npts - SYpoint_fir +1)))


coverage_np[5, ] <- 0 # assuming non_PWID in prison not the targeted population

save(coverage_np,
     file = file.path(OutputFolder ,
                      paste0(project_name,"coverage_np" ,".rda")))

S_Yint = 2022 
S_Yfir = 2024
S_Ymid = 2027
S_Yend <-2050
SYpoint_int <- (S_Yint - POC_AU$cabY)/POC_AU$timestep + 1
SYpoint_fir <- (S_Yfir - POC_AU$cabY)/POC_AU$timestep + 1
SYpoint_mid <- (S_Ymid - POC_AU$cabY)/POC_AU$timestep + 1
SYpoint_end <- (S_Yend - POC_AU$cabY)/POC_AU$timestep + 1
SY_leng <- SYpoint_end - SYpoint_int
SY_lengmid <- SYpoint_mid - SYpoint_int
SY_lengfirtomid <- SYpoint_mid - SYpoint_fir


coverage_npscale <- coverage_np
intVal <- coverage_np[, SYpoint_fir]
coverage_npscale[1, SYpoint_fir:POC_AU$npts] <- 
  c(seq(intVal[1], Ccal_scale$C, length = (SY_lengfirtomid)), rep(Ccal_scale$C, (POC_AU$npts - SYpoint_mid +1)))
coverage_npscale[2, SYpoint_fir:POC_AU$npts] <- 
  c(seq(intVal[2], Ccal_scale$C, length = (SY_lengfirtomid)), rep(Ccal_scale$C, (POC_AU$npts - SYpoint_mid +1)))
coverage_npscale[3, SYpoint_fir:POC_AU$npts] <- 
  c(seq(intVal[3], Ccal_scale$P, length = (SY_lengfirtomid)), rep(Ccal_scale$P, (POC_AU$npts - SYpoint_mid +1)))
coverage_npscale[4, SYpoint_fir:POC_AU$npts] <- 
  c(seq(intVal[4], Ccal_scale$P, length = (SY_lengfirtomid)), rep(Ccal_scale$P, (POC_AU$npts - SYpoint_mid +1)))

coverage_npscale[5, ] <- 0 # assuming non_PWID in prison not the targeted population

save(coverage_npscale,
     file = file.path(OutputFolder ,
                      paste0(project_name,"coverage_npscale" ,".rda")))

coverage_npscale[, 110]
coverage_np[, 110]
